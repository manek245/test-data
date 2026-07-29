# ===============================================================
# Oncology synthetic: ADSL + ADTR + ADRS + ADTTE
# - Variable visit calendars + jitter
# - Two-visit missed blocks → censor at last assessment before block
# - PFS = min(PD, Death), else censored (including missed-block rule)
# - No subject has tumors in both liver and kidney
# ===============================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(stringr)
})

set.seed(20251101)

# --------------------------- PARAMS ----------------------------
STUDYID    <- "ONCOLOGY-001"
N_SUBJ     <- 80
SITE_IDS   <- sprintf("%03d", 1:6)
RAND_START <- as.Date("2024-01-15")
RAND_END   <- as.Date("2024-06-30")
ABS_PD_CM  <- 0.5                 # 5 mm absolute rise
VIS_INT    <- 42                  # nominal interval (days)
VIS_JIT    <- 10                  # +/- jitter in days
MIN_VIS    <- 3; MAX_VIS <- 8     # number of visits per subject
P_DROPOUT  <- 0.28                # may truncate schedule
P_NEWLES   <- 0.25                # true new lesion
P_DEATH    <- 0.20                # death probability
P_MISSBLK  <- 0.25                # create a 2-visit missed block
# --------------------------- PARAMS (add) ---------------------------
P_RESP         <- 0.45  # % of subjects who respond (PR or CR)
P_CR_IN_RESP   <- 0.30  # among responders, % that go to CR (rest PR)


# ------------------------ HELPERS ------------------------------
sample_date <- function(n, start, end) start + sample.int(as.integer(end - start + 1L), n, TRUE) - 1L
merge_parts <- function(x) strsplit(x, "/", fixed = TRUE)

# Per-subject schedule (with jitter and optional dropout)
make_schedule <- function() {
  n_vis <- sample(MIN_VIS:MAX_VIS, 1)
  deltas <- pmax(7, round(rnorm(n_vis - 1, VIS_INT, VIS_JIT)))
  ADY <- c(0, cumsum(deltas))
  if (runif(1) < P_DROPOUT && n_vis >= 4) {
    keep_to <- sample(2:(n_vis-1), 1)
    ADY <- ADY[seq_len(keep_to)]
  }
  tibble(AVISITN = seq_along(ADY),
         AVISIT  = paste0("Visit ", seq_along(ADY)),
         ADY     = ADY,
         ASSESS  = 1L)  # we will zero-out missed visits later
}

# Create a 2-visit missed block (when feasible), return censor cutoff date
apply_missed_block <- function(sched, randdt) {
  s <- sched
  censor_adt <- as.Date(NA)
  # choose start >= 3 and ensure start+1 exists
  cand <- which(s$AVISITN >= 3 & s$AVISITN <= max(s$AVISITN) - 1)
  if (length(cand) >= 1 && runif(1) < P_MISSBLK) {
    i0 <- sample(cand, 1)
    s$ASSESS[s$AVISITN %in% c(i0, i0 + 1)] <- 0L
    # censor at last assessed visit before the block
    prev <- max(s$AVISITN[s$AVISITN < i0 & s$ASSESS == 1L], na.rm = TRUE)
    if (is.finite(prev)) censor_adt <- as.Date(randdt) + s$ADY[s$AVISITN == prev]
  }
  list(sched = s, censor_adt = censor_adt)
}

# Put near your helpers
derived_from_baseline <- function(id, baseline_ids) {
  if (id %in% baseline_ids) return(TRUE)
  if (grepl("\\.", id)) return(sub("\\..*$", "", id) %in% baseline_ids)  # split child -> parent in baseline
  if (grepl("/", id, fixed = TRUE)) {
    parts <- strsplit(id, "/", fixed = TRUE)[[1]]
    return(all(parts %in% baseline_ids))                                  # merge of baseline parts
  }
  FALSE
}


# ----------------------------- ADSL ----------------------------
arm_df <- tibble::tibble(
  ARM   = c("Placebo", "Drug A 200 mg QD", "Drug B 100 mg BID"),
  ARMCD = c("PBO", "A200QD", "B100BID")
)

adsl <- tibble::tibble(
  SUBJSEQ = 1:N_SUBJ,
  SITEID  = sample(SITE_IDS, N_SUBJ, TRUE)
) %>%
  mutate(
    SUBJID   = sprintf("%03d", SUBJSEQ),
    USUBJID  = paste(SITEID, SUBJID, sep = "-"),
    SEX      = sample(c("M","F"), N_SUBJ, TRUE, c(0.55,0.45)),
    AGE      = pmin(pmax(round(rnorm(N_SUBJ, 60, 9)), 22), 88),
    RACE     = sample(c("WHITE","ASIAN","BLACK OR AFRICAN AMERICAN","OTHER"),
                      N_SUBJ, TRUE, c(0.65,0.20,0.10,0.05)),
    HEIGHT   = if_else(SEX=="M", round(rnorm(N_SUBJ, 176, 7)), round(rnorm(N_SUBJ, 164, 6))),
    WEIGHTBL = if_else(SEX=="M", round(rnorm(N_SUBJ, 82, 12),1), round(rnorm(N_SUBJ, 68, 10),1)),
    RANDDT   = sample_date(N_SUBJ, RAND_START, RAND_END),
    ARMCD    = sample(arm_df$ARMCD, N_SUBJ, TRUE)
  ) %>%
  left_join(arm_df, by = "ARMCD") %>%
  mutate(
    ACTARM   = ARM,  ACTARMCD = ARMCD,
    TRT01P   = ARM,  TRT01A   = ACTARM,
    TRTSDT   = RANDDT,
    TRTENDT  = TRTSDT + sample(120:360, N_SUBJ, TRUE),
    DIED     = rbinom(N_SUBJ, 1, P_DEATH),
    DTHDT    = ifelse(DIED==1, RANDDT + sample(35:420, N_SUBJ, TRUE), as.Date(NA)),
    DTHFL    = ifelse(DIED==1, "Y", NA_character_)
  ) %>%
  mutate(RANDDT = as.Date(RANDDT), DTHDT = as.Date(DTHDT)) %>%
  transmute(
    STUDYID = STUDYID,
    USUBJID, SUBJID, SITEID, AGE, SEX, RACE, HEIGHT, WEIGHTBL,
    ARM, ARMCD, ACTARM, ACTARMCD, TRT01P, TRT01A, RANDDT, TRTSDT, TRTENDT,
    DTHDT, DTHFL
  )

trt_vars <- adsl %>%
  dplyr::select(
    USUBJID,
    ARM, ARMCD,
    ACTARM, ACTARMCD,
    TRT01P, TRT01A,
    RANDDT, TRTSDT, TRTENDT
  )

# ---------------------- Response profile per subject ----------------
resp_profile <- adsl %>%
  dplyr::transmute(
    USUBJID,
    RESP = sample(
      c("NR","PR","CR"),
      size = dplyr::n(),
      replace = TRUE,
      prob = c(1 - P_RESP, P_RESP * (1 - P_CR_IN_RESP), P_RESP * P_CR_IN_RESP)
    )
  )

# handy named lookup for later
resp_map <- setNames(resp_profile$RESP, resp_profile$USUBJID)


# ---------------------- Base tumor evolutions ------------------
make_evolution <- function(...) setNames(list(...), paste0("v", seq_along(list(...))))
evo_liver <- make_evolution(
  list(t=c("T01","T02","T03"), d=c(3.5,2.8,2.1), o="liver"),
  list(t=c("T01","T02","T03","T04"), d=c(3.8,3.0,2.3,1.8), o="liver"),
  list(t=c("T01","T02","T03.1","T03.2","T04"), d=c(4.0,3.2,1.3,1.2,1.9), o="liver"),
  list(t=c("T01","T02","T03.1","T03.2","T04"), d=c(4.2,3.5,1.5,1.3,2.0), o="liver"),
  list(t=c("T01/T02","T03.1","T03.2","T04"), d=c(6.0,1.7,1.5,2.1), o="liver"),
  list(t=c("T01/T02","T03.1","T04"), d=c(6.5,1.9,2.3), o="liver"),
  list(t=c("T01/T02","T03.1"), d=c(10,2.1), o="liver")
)
evo_lung <- make_evolution(
  list(t=c("T06","T07"), d=c(2.5,1.8), o="lung"),
  list(t=c("T06","T07","T08"), d=c(2.7,1.9,1.2), o="lung"),
  list(t=c("T06","T07","T08"), d=c(2.9,2.0,1.3), o="lung"),
  list(t=c("T06/T07","T08","T09"), d=c(4.0,1.4,0.9), o="lung"),
  list(t=c("T06/T07","T08","T09"), d=c(4.3,1.5,1.0), o="lung"),
  list(t=c("T06/T07","T08"), d=c(4.6,1.6), o="lung"),
  list(t=c("T06/T07","T08"), d=c(4.9,1.7), o="lung")
)
evo_brain <- make_evolution(
  list(t=c("T10"), d=c(1.2), o="brain"),
  list(t=c("T10"), d=c(1.3), o="brain"),
  list(t=c("T10.1","T10.2","T11"), d=c(0.7,0.6,0.6), o="brain"),
  list(t=c("T10.1","T10.2","T11","T12"), d=c(0.8,0.7,0.7,0.5), o="brain"),
  list(t=c("T10.1","T10.2","T11","T12"), d=c(0.9,0.8,0.8,0.6), o="brain"),
  list(t=c("T10.1","T10.2","T11"), d=c(1.0,0.9,0.9), o="brain"),
  list(t=c("T10.1","T10.2"), d=c(1.1,1.0), o="brain")
)
evo_kidney <- make_evolution(
  list(t=c("T20","T21"), d=c(2.2,1.6), o="kidney"),
  list(t=c("T20","T21","T22"), d=c(2.4,1.7,1.1), o="kidney"),
  list(t=c("T20","T21","T22"), d=c(2.6,1.8,1.2), o="kidney"),
  list(t=c("T20/T21","T22","T23"), d=c(3.6,1.3,0.9), o="kidney"),
  list(t=c("T20/T21","T22","T23"), d=c(3.8,1.4,1.0), o="kidney"),
  list(t=c("T20/T21","T22"), d=c(4.0,1.5), o="kidney"),
  list(t=c("T20/T21","T22"), d=c(4.2,1.6), o="kidney")
)

make_rows_sched <- function(evo, usubjid, sched) {
  out <- vector("list", 0L); i <- 0L
  K <- min(length(evo), nrow(sched))
  for (k in seq_len(K)) {
    if (sched$ASSESS[k] == 0L) next  # skip missed visit
    visit <- k; t <- evo[[k]]$t; d <- evo[[k]]$d; o <- evo[[k]]$o
    for (j in seq_along(t)) {
      i <- i + 1L
      out[[i]] <- data.frame(
        USUBJID = usubjid, STUDYID = STUDYID,
        AVISIT  = sched$AVISIT[visit], AVISITN = sched$AVISITN[visit], ADY = sched$ADY[visit],
        TRLOC   = o, TRGRPID = t[j], TRLNKID = t[j],
        TRTESTCD="LDIAM", TRTEST="Longest Diameter",
        TRORRES = as.character(d[j]), TRORRESU="cm",
        TRSTRESN = d[j], TRSTRESU="cm", stringsAsFactors = FALSE
      )
    }
  }
  dplyr::bind_rows(out)
}

# organ assignment w/ constraint: never both liver & kidney
assign_organs <- function(n) {
  base <- c("liver","lung","brain","kidney")
  p <- c(liver=0.40, lung=0.30, brain=0.18, kidney=0.30)
  size_probs <- c(`1`=0.45, `2`=0.45, `3`=0.10)
  res <- vector("list", n)
  for (i in seq_len(n)) {
    k <- sample(1:3, 1, prob = size_probs)
    combos <- combn(base, k, simplify = FALSE)
    combos <- Filter(function(s) !("liver" %in% s && "kidney" %in% s), combos)
    w <- sapply(combos, function(s) prod(p[s])); w <- w/sum(w)
    res[[i]] <- combos[[sample.int(length(combos), 1, prob = w)]]
  }
  res
}

# subject variation + optional new lesion + drift
# subject variation + optional new lesion + drift + response shaping
apply_subject_variation <- function(df_subj, profile = "NR") {
  # --- lock responders to baseline-derived lesions only (prevents NEWL) ---
  if (profile %in% c("PR","CR")) {
    base_ids <- df_subj %>%
      dplyr::filter(AVISITN == 1) %>%
      dplyr::distinct(TRLNKID) %>% dplyr::pull(TRLNKID)
    keep_row <- vapply(df_subj$TRLNKID, derived_from_baseline, logical(1), baseline_ids = base_ids)
    df_subj  <- df_subj %>% dplyr::filter(AVISITN == 1 | keep_row)
  }

  # --- base scale + drift ---
  scale_s <- rlnorm(1, 0, 0.20)
  v       <- sort(unique(df_subj$AVISITN))
  if (length(v) == 0) return(df_subj)

  base_drift <- cumprod(1 + rnorm(length(v), mean = rnorm(1, 0.00, 0.03), sd = 0.05))
  base_drift <- base_drift / base_drift[1]; names(base_drift) <- v

  # --- response shaping ---
  resp_curve <- rep(1, length(v)); names(resp_curve) <- v
  v_start <- NA_integer_; v_nadir <- NA_integer_; target_frac <- 1

  if (profile %in% c("PR","CR") && length(v) >= 3) {
    v_start <- sample(v[v >= 2 & v <= max(v) - 1], 1)
    v_nadir <- sample(v[v > v_start], 1)
    # PR: ensure ≤ 70% of baseline; give it some room (25–60%) to see clear PRs
    target_frac <- if (profile == "PR") runif(1, 0.25, 0.60) else 0.0

    for (i in seq_along(v)) {
      if (v[i] < v_start) {
        resp_curve[i] <- 1
      } else if (v[i] <= v_nadir) {
        alpha <- (v[i] - v_start) / (v_nadir - v_start)
        resp_curve[i] <- 1 + alpha * (target_frac - 1)
      } else {
        resp_curve[i] <- target_frac
      }
    }
  }

  mult      <- base_drift * resp_curve
  floor_low <- if (identical(profile, "CR")) 0.0 else 0.1  # allow true 0 only for CR

  df_subj <- df_subj %>%
    dplyr::mutate(
      mult      = mult[as.character(AVISITN)],
      TRSTRESN  = pmax(floor_low, TRSTRESN * scale_s * mult),
      TRORRES   = sprintf("%.1f", TRSTRESN)
    )

  # --- no injected new lesions for responders ---
  allow_newles <- !(profile %in% c("PR","CR"))
  if (allow_newles && runif(1) < P_NEWLES && max(df_subj$AVISITN) >= 3) {
    v_new <- sample(unique(df_subj$AVISITN)[unique(df_subj$AVISITN) >= 3], 1)
    org   <- sample(unique(df_subj$TRLOC), 1)
    size0 <- runif(1, 0.8, 1.4); grow <- runif(1, 1.05, 1.20)
    sch   <- df_subj %>% dplyr::distinct(AVISITN, AVISIT, ADY) %>% dplyr::arrange(AVISITN)
    add <- purrr::map_dfr(seq_len(nrow(sch)), function(i) {
      vv <- sch$AVISITN[i]
      if (vv < v_new) return(NULL)
      xx <- size0 * (grow^(vv - v_new))
      data.frame(
        USUBJID  = df_subj$USUBJID[1], STUDYID = STUDYID,
        AVISIT   = sch$AVISIT[i], AVISITN = vv, ADY = sch$ADY[i],
        TRLOC    = org, TRGRPID = "TN1", TRLNKID = "TN1",
        TRTESTCD = "LDIAM", TRTEST = "Longest Diameter",
        TRORRES  = sprintf("%.1f", xx), TRORRESU = "cm",
        TRSTRESN = xx, TRSTRESU = "cm",
        stringsAsFactors = FALSE
      )
    })
    df_subj <- dplyr::bind_rows(df_subj, add)
  }

  # --- enforce CR zeros at/after nadir ---
  if (identical(profile, "CR") && !is.na(v_nadir)) {
    df_subj <- df_subj %>%
      dplyr::mutate(
        TRSTRESN = dplyr::if_else(AVISITN >= v_nadir, 0, TRSTRESN),
        TRORRES  = dplyr::if_else(AVISITN >= v_nadir, "0.0", TRORRES)
      )
  }

  df_subj %>% dplyr::select(-mult) %>% dplyr::arrange(AVISITN, TRLOC, TRLNKID)
}

# --------------------------- Build ADTR ------------------------
organ_list <- assign_organs(nrow(adsl))

# schedules + missed-block censor map
sched_list <- vector("list", nrow(adsl))
censor_map <- tibble(USUBJID = adsl$USUBJID, CENSOR_ADT = as.Date(NA))
for (i in seq_len(nrow(adsl))) {
  sc0 <- make_schedule()
  blk <- apply_missed_block(sc0, adsl$RANDDT[i])
  sched_list[[i]] <- blk$sched
  censor_map$CENSOR_ADT[i] <- blk$censor_adt
}

adtr_base <- map2_dfr(seq_len(nrow(adsl)), organ_list, function(i, organs) {
  usubjid <- adsl$USUBJID[i]; sched <- sched_list[[i]]
  rbind(
    if ("liver"  %in% organs) make_rows_sched(evo_liver,  usubjid, sched) else NULL,
    if ("lung"   %in% organs) make_rows_sched(evo_lung,   usubjid, sched) else NULL,
    if ("brain"  %in% organs) make_rows_sched(evo_brain,  usubjid, sched) else NULL,
    if ("kidney" %in% organs) make_rows_sched(evo_kidney, usubjid, sched) else NULL
  )
})

adtr <- adtr_base %>%
  dplyr::group_split(USUBJID, .keep = TRUE) %>%
  purrr::map_dfr(function(df_subj) {
    prof <- resp_map[df_subj$USUBJID[1]]
    apply_subject_variation(df_subj, profile = unname(prof))
  })

adtr <- adtr %>%
  dplyr::left_join(trt_vars, by = "USUBJID") %>%
  dplyr::relocate(ARM, ARMCD, ACTARM, ACTARMCD, TRT01P, TRT01A, RANDDT, TRTSDT, TRTENDT,
                  .after = USUBJID)


# --------------------------- ADRS ------------------------------
is_visit_new_lesion <- function(ids_at_visit, baseline_ids) {
  sapply(ids_at_visit, function(id) {
    if (id %in% baseline_ids) return(FALSE)
    if (grepl("/", id, fixed = TRUE)) {
      parts <- unlist(merge_parts(id))
      return(!all(parts %in% baseline_ids))
    }
    if (grepl("\\.", id)) {
      parent <- sub("\\..*$", "", id)
      return(!(parent %in% baseline_ids))
    }
    TRUE
  })
}
get_baseline_targets <- function(df_subj) df_subj %>% filter(AVISITN == 1) %>% distinct(TRLNKID) %>% pull(TRLNKID)

compute_sld_one_visit <- function(dv, baseline_ids) {
  ids <- dv$TRLNKID; vals <- dv$TRSTRESN
  contrib <- vapply(baseline_ids, function(b) {
    v <- sum(vals[ids == b], na.rm = TRUE)
    child_idx <- which(startsWith(ids, paste0(b, ".")))
    if (length(child_idx)) v <- v + sum(vals[child_idx], na.rm = TRUE)
    merge_idx <- which(grepl("/", ids, fixed = TRUE) & grepl(paste0("(^|/)", b, "(/|$)"), ids))
    if (length(merge_idx)) {
      add <- mapply(function(id, x) x/length(unlist(merge_parts(id))), id = ids[merge_idx], x = vals[merge_idx])
      v <- v + sum(add, na.rm = TRUE)
    }
    v
  }, numeric(1))
  sum(contrib, na.rm = TRUE)
}

adtr_l <- adtr %>% filter(TRTESTCD == "LDIAM") %>% mutate(TRSTRESN = as.numeric(TRSTRESN))

sld_by_subj <- adtr_l %>%
  arrange(USUBJID, AVISITN) %>%
  group_split(USUBJID, .keep = TRUE) %>%
  map_dfr(function(df_subj) {
    usubjid <- df_subj$USUBJID[1]
    baseline_ids <- get_baseline_targets(df_subj)
    visits <- df_subj %>% distinct(AVISITN, AVISIT, ADY) %>% arrange(AVISITN)
    map_dfr(seq_len(nrow(visits)), function(i) {
      v <- visits$AVISITN[i]
      dv <- df_subj %>% filter(AVISITN == v)
      sld <- compute_sld_one_visit(dv, baseline_ids)
      ids_here <- dv$TRLNKID
      new_flag <- any(is_visit_new_lesion(ids_here, baseline_ids))
      tibble(STUDYID=STUDYID, USUBJID=usubjid, AVISITN=v, AVISIT=visits$AVISIT[i],
             ADY=visits$ADY[i], SLD=sld, NEWL=as.integer(new_flag))
    })
  }) %>%
  group_by(USUBJID) %>%
  arrange(AVISITN, .by_group = TRUE) %>%
  mutate(BSLD = first(SLD), NADIR = cummin(SLD),
         CHG = SLD - BSLD, PCHG = ifelse(BSLD>0, 100*CHG/BSLD, NA_real_)) %>%
  ungroup()

recist_call <- function(SLD, NADIR, NEWL, PCHG) {
  if (isTRUE(SLD == 0)) return("CR")
  if (isTRUE(NEWL == 1)) return("PD")
  if (!is.na(NADIR) && !is.na(SLD) && (SLD >= 1.2*NADIR) && (SLD - NADIR >= ABS_PD_CM)) return("PD")
  if (!is.na(PCHG) && (PCHG <= -30)) return("PR")
  "SD"
}

ADRS <- sld_by_subj %>%
  rowwise() %>%
  mutate(PARAMCD="OVR", PARAM="Overall Response (RECIST 1.1)",
         AVAL=SLD, AVALC=recist_call(SLD, NADIR, NEWL, PCHG)) %>%
  ungroup() %>%
  left_join(adsl %>% select(USUBJID, RANDDT), by="USUBJID") %>%
  mutate(ADT = as.Date(RANDDT) + ADY) %>%
  select(STUDYID, USUBJID, PARAMCD, PARAM, AVISITN, AVISIT, ADY, ADT,
         BSLD, NADIR, SLD=AVAL, CHG, PCHG, NEWL, AVALC)

ADRS <- ADRS %>%
  dplyr::left_join(trt_vars, by = "USUBJID") %>%
  dplyr::relocate(ARM, ARMCD, ACTARM, ACTARMCD, TRT01P, TRT01A, RANDDT, TRTSDT, TRTENDT,
                  .after = USUBJID)


# ---------------------------- ADTTE ----------------------------
# PD table
adtte_pd <- ADRS %>%
  arrange(USUBJID, ADT, AVISITN) %>%
  group_by(USUBJID) %>%
  summarise(PD_ADT = first(ADT[AVALC=="PD"]),
            PD_ADY = first(ADY[AVALC=="PD"]),
            LAST_ADT = last(ADT), LAST_ADY = last(ADY), .groups="drop")

# Death table
adtte_death <- adsl %>%
  select(USUBJID, RANDDT, DTHDT) %>%
  mutate(RANDDT = as.Date(RANDDT), DTHDT = as.Date(DTHDT),
         DTH_ADY = ifelse(is.na(DTHDT), NA_real_, as.numeric(DTHDT - RANDDT)))

# Missed-block censor cutoff (from schedule step)
censor_map <- censor_map %>% mutate(CENSOR_ADT = as.Date(CENSOR_ADT))

# Build ADTTE: PFS = earliest(PD, Death), unless a missed-block censor date exists earlier.
ADTTE <- adtte_pd %>%
  left_join(adtte_death, by = "USUBJID") %>%
  left_join(censor_map,  by = "USUBJID") %>%
  # ensure all date-like columns are Date
  mutate(
    PD_ADT     = as.Date(PD_ADT),
    DTHDT      = as.Date(DTHDT),
    LAST_ADT   = as.Date(LAST_ADT),
    CENSOR_ADT = as.Date(CENSOR_ADT),
    RANDDT     = as.Date(RANDDT)
  ) %>%
  # ---- SAFE earliest-of(PD, Death) as Date ----
mutate(
  .early_num     = pmin(as.numeric(PD_ADT), as.numeric(DTHDT), na.rm = TRUE),
  .early_num     = ifelse(is.infinite(.early_num), NA_real_, .early_num),
  EARLY_EVT_ADT  = as.Date(.early_num, origin = "1970-01-01")
) %>%
  select(-.early_num) %>%
  # ---- Apply missed-block censoring if it is earlier than the earliest event ----
mutate(
  USE_CENSOR = !is.na(CENSOR_ADT) & (is.na(EARLY_EVT_ADT) | CENSOR_ADT < EARLY_EVT_ADT),
  CNSR       = ifelse(USE_CENSOR | is.na(EARLY_EVT_ADT), 1L, 0L),
  # Final ADT (Date) and AVAL (days)
  ADT  = dplyr::case_when(
    USE_CENSOR ~ CENSOR_ADT,
    CNSR == 0L ~ EARLY_EVT_ADT,
    TRUE       ~ LAST_ADT
  ),
  AVAL = dplyr::case_when(
    USE_CENSOR ~ as.numeric(CENSOR_ADT - RANDDT),
    CNSR == 0L & !is.na(PD_ADT)  & (is.na(DTHDT) | PD_ADT <= DTHDT) ~ PD_ADY,
    CNSR == 0L & !is.na(DTHDT)   & (is.na(PD_ADT) | DTHDT <  PD_ADT) ~ as.numeric(DTHDT - RANDDT),
    TRUE ~ LAST_ADY
  ),
  AVALU = "days",
  EVNTDESC = dplyr::case_when(
    USE_CENSOR ~ "Censored at last assessment before missed-visit block",
    CNSR==1L   ~ "Censored at last tumor assessment",
    !is.na(PD_ADT) & (is.na(DTHDT) | PD_ADT <= DTHDT) ~ "First progression (PD)",
    TRUE ~ "Death (no prior PD)"
  ),
  STUDYID = STUDYID,
  PARAMCD = "PFS",
  PARAM   = "Progression-Free Survival (RECIST 1.1)"
) %>%
  select(STUDYID, USUBJID, PARAMCD, PARAM, CNSR, AVAL, AVALU, ADT, EVNTDESC)

ADTTE <- ADTTE %>%
  dplyr::left_join(trt_vars, by = "USUBJID") %>%
  dplyr::relocate(ARM, ARMCD, ACTARM, ACTARMCD, TRT01P, TRT01A, RANDDT, TRTSDT, TRTENDT,
                  .after = USUBJID)

# ======================= ADLB + ADAE ==========================
# Assumes you already have: STUDYID, adsl, trt_vars, sched_list, RAND_START/END

# ---------- ADLB (ALT/AST/TBIL with Hy's Law) ----------------
# ULN constants (tweak if you prefer per-subject ULNs)
ULN <- list(
  ALT  = 40,    # U/L
  AST  = 40,    # U/L
  TBIL = 1.2    # mg/dL
)

# Lab parameter metadata
lb_meta <- tibble::tribble(
  ~PARAMCD, ~PARAM,                          ~UNIT,   ~ULN,
  "ALT",    "Alanine Aminotransferase",      "U/L",   ULN$ALT,
  "AST",    "Aspartate Aminotransferase",    "U/L",   ULN$AST,
  "TBIL",   "Total Bilirubin",               "mg/dL", ULN$TBIL
)

# Helper to generate subject's lab trajectory per PARAM
gen_param_vals <- function(n_vis, paramcd) {
  # Baseline around 0.8–1.2 x ULN (ALT/AST) or 0.6–1.0 x ULN (TBIL)
  if (paramcd %in% c("ALT","AST")) {
    base_mult <- runif(1, 0.8, 1.2)
    drift <- cumprod(1 + rnorm(n_vis, 0.00, 0.05))
    vals <- pmax(0.05, base_mult * lb_meta$ULN[lb_meta$PARAMCD==paramcd] * drift)
  } else { # TBIL
    base_mult <- runif(1, 0.6, 1.0)
    drift <- cumprod(1 + rnorm(n_vis, 0.00, 0.03))
    vals <- pmax(0.05, base_mult * lb_meta$ULN[lb_meta$PARAMCD==paramcd] * drift)
  }
  vals
}

# We will spike a subset of subjects to meet Hy's Law at one assessment
# Hy's Law (simplified): ALT or AST >= 3x ULN *AND* TBIL >= 2x ULN at the same visit
P_HY_SUBJ <- 0.10  # ~10% of subjects will have at least one Hy's Law visit

ADLB <- purrr::map2_dfr(seq_len(nrow(adsl)), adsl$USUBJID, function(i, usubjid) {
  sched <- sched_list[[i]] %>% dplyr::filter(ASSESS == 1L) %>% dplyr::arrange(AVISITN)
  if (nrow(sched) == 0) return(NULL)

  # Seed baseline/trajectory
  vals_ALT  <- gen_param_vals(nrow(sched), "ALT")
  vals_AST  <- gen_param_vals(nrow(sched), "AST")
  vals_TBIL <- gen_param_vals(nrow(sched), "TBIL")

  # Possibly inject a Hy's Law spike
  if (runif(1) < P_HY_SUBJ && nrow(sched) >= 3) {
    hy_idx <- sample(2:nrow(sched), 1)
    # Force ALT spike on that visit (or AST with 50%)
    if (runif(1) < 0.5) {
      vals_ALT[hy_idx] <- max(vals_ALT[hy_idx], 3.2 * ULN$ALT * runif(1, 1.0, 1.6))
    } else {
      vals_AST[hy_idx] <- max(vals_AST[hy_idx], 3.2 * ULN$AST * runif(1, 1.0, 1.6))
    }
    vals_TBIL[hy_idx] <- max(vals_TBIL[hy_idx], 2.1 * ULN$TBIL * runif(1, 1.0, 1.4))
  }

  # Build long rows
  dplyr::bind_rows(
    tibble::tibble(USUBJID=usubjid, PARAMCD="ALT",  AVAL=vals_ALT,  AVISITN=sched$AVISITN, AVISIT=sched$AVISIT, ADY=sched$ADY),
    tibble::tibble(USUBJID=usubjid, PARAMCD="AST",  AVAL=vals_AST,  AVISITN=sched$AVISITN, AVISIT=sched$AVISIT, ADY=sched$ADY),
    tibble::tibble(USUBJID=usubjid, PARAMCD="TBIL", AVAL=vals_TBIL, AVISITN=sched$AVISITN, AVISIT=sched$AVISIT, ADY=sched$ADY)
  )
}) %>%
  dplyr::left_join(adsl %>% dplyr::select(USUBJID, RANDDT), by="USUBJID") %>%
  dplyr::left_join(lb_meta, by="PARAMCD") %>%
  dplyr::mutate(
    STUDYID = STUDYID,
    ADT     = as.Date(RANDDT) + ADY,
    AVALU   = UNIT,
    A1LO    = NA_real_,                     # (optional LLOQ/ULOQ placeholders)
    A1HI    = NA_real_,
    ANRLO   = 0,                            # lower ref (for illustration)
    ANRHI   = ULN,
    ABLFL   = dplyr::if_else(AVISITN == 1, "Y", NA_character_),
    AVISITN = as.integer(AVISITN)
  ) %>%
  dplyr::select(STUDYID, USUBJID, PARAMCD, PARAM, AVISITN, AVISIT, ADY, ADT,
                AVAL, AVALU, ANRLO, ANRHI, ABLFL) %>%
  dplyr::left_join(trt_vars, by="USUBJID") %>%
  dplyr::relocate(ARM, ARMCD, ACTARM, ACTARMCD, TRT01P, TRT01A, RANDDT, TRTSDT, TRTENDT,
                  .after = USUBJID)

# Hy's Law visit-level flag (simple)
hy_visit <- ADLB %>%
  dplyr::select(USUBJID, AVISITN, ADT, PARAMCD, AVAL) %>%
  tidyr::pivot_wider(names_from = PARAMCD, values_from = AVAL) %>%
  dplyr::mutate(
    ALT_ULN  = ULN$ALT,
    AST_ULN  = ULN$AST,
    TBIL_ULN = ULN$TBIL,
    HYSLAW   = ((ALT >= 3*ALT_ULN | AST >= 3*AST_ULN) & (TBIL >= 2*TBIL_ULN))
  ) %>%
  dplyr::select(USUBJID, AVISITN, ADT, HYSLAW)

ADLB <- ADLB %>%
  dplyr::left_join(hy_visit, by=c("USUBJID","AVISITN","ADT")) %>%
  dplyr::mutate(HYSLAWFL = dplyr::if_else(HYSLAW, "Y", NA_character_)) %>%
  dplyr::select(-HYSLAW)

# Subject-level earliest Hy's Law date (optional helper table)
HYSL <- ADLB %>%
  dplyr::filter(HYSLAWFL == "Y") %>%
  dplyr::group_by(USUBJID) %>%
  dplyr::summarise(HYADT = min(ADT), .groups="drop")

# --------------- ADAE (with some SAEs) ------------------------
# Some common AE terms/codes
ae_terms <- tibble::tribble(
  ~AEDECOD,                 ~AETERM,
  "NAUSEA",                 "Nausea",
  "DIARRHOEA",              "Diarrhea",
  "FATIGUE",                "Fatigue",
  "AST INCREASED",          "Aspartate aminotransferase increased",
  "ALT INCREASED",          "Alanine aminotransferase increased",
  "BILIRUBIN INCREASED",    "Blood bilirubin increased",
  "HEADACHE",               "Headache",
  "RASH",                   "Rash",
  "HYPERTRANSAMINASAEMIA",  "Hypertransaminasaemia",
  "INFUSION RELATED REACTION", "Infusion-related reaction"
)

# Probabilities
P_HAS_AE   <- 0.70        # subject has at least one AE
P_SAE      <- 0.15        # each AE has chance to be serious
P_DRUGREL  <- 0.55
SEV_LVL    <- c("MILD","MODERATE","SEVERE")
SEV_PROB   <- c(0.45, 0.40, 0.15)

ADAE <- purrr::map_dfr(seq_len(nrow(adsl)), function(i) {
  usubjid <- adsl$USUBJID[i]
  if (runif(1) > P_HAS_AE) return(NULL)

  # 1–4 AEs with a skew toward fewer AEs
  n_ae <- sample(1:4, 1, prob = c(0.5, 0.3, 0.15, 0.05))

  randdt <- as.Date(adsl$RANDDT[i]); trtend <- as.Date(adsl$TRTENDT[i])
  if (is.na(trtend) || trtend < randdt) trtend <- randdt + 120

  purrr::map_dfr(seq_len(n_ae), function(j) {
    sel   <- ae_terms[sample.int(nrow(ae_terms), 1), ]
    dthdt <- as.Date(adsl$DTHDT[i])  # define BEFORE using it

    # sample a start within treatment + 30 days (guard upper bound ≥ 1)
    span_days <- max(1L, as.integer(trtend - randdt) + 30L)
    start     <- randdt + sample.int(span_days, 1)

    # duration 1–14 days; provisional end
    dur  <- sample.int(14L, 1)
    end  <- start + dur
    aeser <- ifelse(runif(1) < P_SAE, "Y", "N")

    # If died, sometimes make a fatal AE ending exactly on DTHDT
    if (!is.na(dthdt) && runif(1) < 0.35) {
      aeser <- "Y"
      end   <- dthdt
    }

    # Safety: never let AE end before it starts
    if (!is.na(end) && end < start) end <- start

    tibble::tibble(
      STUDYID = STUDYID,
      USUBJID = usubjid,
      AESEQ   = j,
      AETERM  = sel$AETERM,
      AEDECOD = sel$AEDECOD,
      AESTDT  = start,
      AEENDT  = end,
      AESTDY  = as.integer(AESTDT - randdt),
      AEENDY  = as.integer(AEENDT - randdt),
      AEREL   = ifelse(runif(1) < P_DRUGREL, "RELATED", "NOT RELATED"),
      AESEV   = sample(SEV_LVL, 1, prob = SEV_PROB),
      AESER   = aeser,
      AEOUT   = dplyr::case_when(
        !is.na(dthdt) && end == dthdt ~ "FATAL",
        runif(1) < 0.10              ~ "ONGOING",
        runif(1) < 0.10              ~ "NOT RECOVERED/NOT RESOLVED",
        TRUE                         ~ "RECOVERED/RESOLVED"
      )
    )
  })

}) %>%
  dplyr::left_join(trt_vars, by="USUBJID") %>%
  dplyr::mutate(
    AEBODSYS = NA_character_,   # (placeholder)
    TRTEMFL  = dplyr::if_else(AESTDT >= TRTSDT & AESTDT <= TRTENDT, "Y", NA_character_),
    AEENDT = dplyr::case_when(
      AEOUT %in% c("ONGOING", "NOT RECOVERED/NOT RESOLVED") ~ as.Date(NA),
      TRUE ~ as.Date(AEENDT)
    ),
    AEENDY = dplyr::if_else(is.na(AEENDT), NA_integer_, as.integer(AEENDT - as.Date(RANDDT)))
  ) %>%
  dplyr::relocate(ARM, ARMCD, ACTARM, ACTARMCD, TRT01P, TRT01A, RANDDT, TRTSDT, TRTENDT,
                  .after = USUBJID)

# Ensure at least some SAEs exist
if (!any(ADAE$AESER == "Y")) {
  if (nrow(ADAE) > 0) ADAE$AESER[sample.int(nrow(ADAE), 1)] <- "Y"
}

# ------------------------- SAVE & QC --------------------------
save(adsl,  file="ADSL.rda")
save(adtr,  file="ADTR.rda")
save(ADRS,  file="ADRS.rda")
save(ADTTE, file="ADTTE.rda")
save(ADLB,  file="ADLB.rda")
save(ADAE,  file="ADAE.rda")
cat("Saved ADSL.rda, ADTR.rda, ADRS.rda, ADTTE.rda, ADLB.rda, ADAE.rda\n")

# Quick QCs:
# --- ADLB / Hy's Law ---
# with(ADLB, table(PARAMCD))
# subset(ADLB, HYSLAWFL == "Y") |> dplyr::arrange(USUBJID, ADT) |> head(10)
# HYSL |> head()
#
# --- ADAE / SAEs ---
# with(ADAE, table(AESER))
# ADAE |> dplyr::arrange(USUBJID, AESTDT) |> head(10)
