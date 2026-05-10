library(survival)
library(Hmisc)
library(splines)
library(here)
library(dplyr)

prostate <- read.csv("prostate.csv")

## data cleaning
prostate$allCause <- prostate$status != "alive"
prostate$eventType <- as.factor(prostate$status)
levels(prostate$eventType) <-
  list(alive = "alive",
       pdeath = "dead - prostatic ca",
       odeath = c(levels(prostate$eventType)[c(2:5, 7:10)]))

prostate$rx <- as.factor(prostate$rx)
prost_red <- prostate[prostate$rx %in% levels(prostate$rx)[3:4], ]
prost_red$temprx <- as.integer(prost_red$rx) - 3
prost_red$rx <- abs(1 - prost_red$temprx) #DES is A=1 and placebo A=0
prost_red$eventType <- as.integer(prost_red$eventType) - 1 #0 is censoring, 1 is pdeath, 2 is odeath
prost_red$hgBinary <- prost_red$hg < 12
prost_red$ageCat <- cut2(prost_red$age, c(0, 60, 70, 80, 100))
prost_red$normalAct <- prost_red$pf == "normal activity"
prost_red$eventCens <- prost_red$eventType == 0
prost_red$eventProst <- prost_red$eventType == 1

prost_red$Tstart <- -0.01
cut_times <- c(0:50)


## build long-format dataset from a (wide) prost_red-style data frame
build_long <- function(data_wide) {
  data_wide$Tstart <- -0.01
  long <- survSplit(data = data_wide, cut = cut_times, start = "Tstart", end = "dtime", event = "allCause")
  long_cens <- survSplit(data = data_wide, cut = cut_times, start = "Tstart", end = "dtime", event = "eventCens")
  long$eventCens <- long_cens$eventCens
  long$prostateDeath <- long$allCause == 1 & long$eventType == 1
  long$otherDeath <- long$allCause == 1 & long$eventType == 2
  long$prostateDeath[long$eventCens == 1] <- NA
  long$otherDeath[long$eventCens == 1] <- NA
  long$prostateDeath[long$otherDeath == 1] <- NA
  long <- long[long$dtime < length(cut_times), ]
  #both estimators need a copy of rx under another name
  long$Orx <- long$rx
  long$Mrx <- long$rx
  return(long)
}

long_prost_red <- build_long(prost_red)
baseline <- long_prost_red[long_prost_red$dtime == 0, ]
n <- length(unique(long_prost_red$patno))


## cumulative incidence (g-formula)
calculateCumInc <- function(inputData, timepts = cut_times, competing = FALSE) {
  cumulativeIncidence <- matrix(NA, ncol = length(unique(inputData$patno)), nrow = length(cut_times))
  #insert event probabilities at the first time interval
  cumulativeIncidence[1, ] <- inputData[inputData$dtime == 0, ]$hazardP * inputData[inputData$dtime == 0, ]$hazardO
  #survival probabilities at each time
  survivalProb <- t(aggregate(s ~ patno, data = inputData, FUN = cumprod)$s)
  for (i in 2:length(cut_times)) {
    subInputDataP <- inputData[inputData$dtime == (i - 1), ]$hazardP
    subInputDataO <- (1 - inputData[inputData$dtime == (i - 1), ]$hazardO)
    if (!competing) {
      cumulativeIncidence[i, ] <- subInputDataO * subInputDataP * survivalProb[(i - 1), ]
    } else {
      cumulativeIncidence[i, ] <- (1 - subInputDataO) * survivalProb[(i - 1), ]
    }
  }
  meanCumulativeIncidence <- rowMeans(apply(cumulativeIncidence, MARGIN = 2, cumsum))
  return(meanCumulativeIncidence)
}

## nonparametric cumulative hazard (IPW)
nonParametricCumHaz <- function(weightVector, inputdata, grp, outcomeProstate = TRUE) {
  outputHazards <- rep(NA, length.out = length(cut_times))
  counter <- 1
  for (i in cut_times) {
    if (outcomeProstate) {
      indices <- inputdata$dtime == i & inputdata$rx == grp & inputdata$eventCens == 0 & inputdata$otherDeath == 0
      eventIndicator <- indices & inputdata$prostateDeath == 1
    } else {
      indices <- inputdata$dtime == i & inputdata$rx == grp & inputdata$eventCens == 0
      eventIndicator <- indices & inputdata$otherDeath == 1
    }
    outputHazards[counter] <- sum(weightVector[eventIndicator]) / sum(weightVector[indices])
    counter <- counter + 1
  }
  return(outputHazards)
}

## nonparametric cumulative incidence (IPW)
nonParametricCumInc <- function(hazard1, hazard2, competing = FALSE) {
  inc <- rep(NA, length.out = length(cut_times))
  cumulativeSurvival <- c(1, cumprod((1 - hazard1) * (1 - hazard2)))
  for (i in 1:length(cut_times)) {
    if (!competing) {
      inc[i] <- hazard1[i] * (1 - hazard2[i]) * cumulativeSurvival[i]
    } else {
      inc[i] <- hazard1[i] * cumulativeSurvival[i]
    }
  }
  cumInc <- cumsum(inc)
  return(cumInc)
}

## discrete cumulative incidence (IPW)
discrete_cuminc_prost <- function(weight_vector, inputdata, grp = 0, outcome_y = TRUE, follow_up = 1:50) {
  event_vec <- rep(NA, length.out = length(follow_up))
  counter <- 1
  n_grp <- sum(inputdata$dtime == 0 & inputdata$rx == grp)
  for (i in follow_up) {
    if (outcome_y) {
      indices <- inputdata$dtime == i & inputdata$rx == grp & inputdata$eventCens == 0 & inputdata$otherDeath == 0
      eventIndicator <- indices & inputdata$prostateDeath == 1
    } else {
      indices <- inputdata$dtime == i & inputdata$rx == grp & inputdata$eventCens == 0
      eventIndicator <- indices & inputdata$otherDeath == 1
    }
    event_vec[counter] <- sum(weight_vector[eventIndicator]) / n_grp
    counter <- counter + 1
  }
  output_cuminc <- cumsum(event_vec)
  return(output_cuminc)
}


## g-estimator: point estimates from a single long dataset
compute_g_estimates <- function(long_data, base_data, n_pat) {
  
  #pooled logistic regression for prostate cancer death
  plrFitP <- glm(prostateDeath ~ rx * (dtime + I(dtime^2) + I(dtime^3)) + normalAct + ageCat + hx + hgBinary,
                 data = long_data, family = binomial())
  
  #pooled logistic regression for death due to other causes
  plrFitO <- glm(otherDeath ~ Orx * (dtime + I(dtime^2) + I(dtime^3)) + normalAct + ageCat + hx + hgBinary,
                 data = long_data, family = binomial())
  
  #counterfactual datasets
  #Ay=1, Ad=1
  treated <- base_data[rep(1:n_pat, each = length(cut_times)), ]
  treated$dtime <- rep(cut_times, n_pat); treated$rx <- 1; treated$Orx <- 1
  #Ay=0, Ad=0
  placebo <- base_data[rep(1:n_pat, each = length(cut_times)), ]
  placebo$dtime <- rep(cut_times, n_pat); placebo$rx <- 0; placebo$Orx <- 0
  #Ay=1, Ad=0
  treatAy <- base_data[rep(1:n_pat, each = length(cut_times)), ]
  treatAy$dtime <- rep(cut_times, n_pat); treatAy$rx <- 1; treatAy$Orx <- 0
  #Ay=0, Ad=1
  treatAd <- base_data[rep(1:n_pat, each = length(cut_times)), ]
  treatAd$dtime <- rep(cut_times, n_pat); treatAd$rx <- 0; treatAd$Orx <- 1
  
  #predict hazards and survival
  treated$hazardP <- predict(plrFitP, newdata = treated, type = "response")
  treated$hazardO <- predict(plrFitO, newdata = treated, type = "response")
  treated$s <- (1 - treated$hazardP) * (1 - treated$hazardO)
  
  placebo$hazardP <- predict(plrFitP, newdata = placebo, type = "response")
  placebo$hazardO <- predict(plrFitO, newdata = placebo, type = "response")
  placebo$s <- (1 - placebo$hazardP) * (1 - placebo$hazardO)
  
  treatAy$hazardP <- predict(plrFitP, newdata = treatAy, type = "response")
  treatAy$hazardO <- predict(plrFitO, newdata = treatAy, type = "response")
  treatAy$s <- (1 - treatAy$hazardP) * (1 - treatAy$hazardO)
  
  treatAd$hazardP <- predict(plrFitP, newdata = treatAd, type = "response")
  treatAd$hazardO <- predict(plrFitO, newdata = treatAd, type = "response")
  treatAd$s <- (1 - treatAd$hazardP) * (1 - treatAd$hazardO)
  
  cumIncTreated <- calculateCumInc(treated)
  cumIncPlacebo <- calculateCumInc(placebo)
  cumIncTreatAy <- calculateCumInc(treatAy)
  cumIncTreatAd <- calculateCumInc(treatAd)
  
  #CDE(0): hazardO set to 0
  treated_de <- base_data[rep(1:n_pat, each = length(cut_times)), ]
  treated_de$dtime <- rep(cut_times, n_pat); treated_de$rx <- 1
  placebo_de <- base_data[rep(1:n_pat, each = length(cut_times)), ]
  placebo_de$dtime <- rep(cut_times, n_pat); placebo_de$rx <- 0
  
  treated_de$hazardP <- predict(plrFitP, newdata = treated_de, type = "response")
  treated_de$hazardO <- 0
  treated_de$s <- (1 - treated_de$hazardP) * (1 - treated_de$hazardO)
  
  placebo_de$hazardP <- predict(plrFitP, newdata = placebo_de, type = "response")
  placebo_de$hazardO <- 0
  placebo_de$s <- (1 - placebo_de$hazardP) * (1 - placebo_de$hazardO)
  
  cumIncTreated_de <- calculateCumInc(treated_de)
  cumIncPlacebo_de <- calculateCumInc(placebo_de)
  
  #extract at month 50
  r_treated <- cumIncTreated[length(cut_times)]
  r_placebo <- cumIncPlacebo[length(cut_times)]
  r_treatAy <- cumIncTreatAy[length(cut_times)]
  r_treatAd <- cumIncTreatAd[length(cut_times)]
  r_treated_de <- cumIncTreated_de[length(cut_times)]
  r_placebo_de <- cumIncPlacebo_de[length(cut_times)]
  
  #effect estimates (risk difference)
  CDE_rd <- r_treated_de - r_placebo_de
  SDE_0_rd <- r_treatAy - r_placebo
  SIE_0_rd <- r_treatAd - r_placebo
  SDE_1_rd <- r_treated - r_treatAd
  SIE_1_rd <- r_treated - r_treatAy
  
  return(c(CDE_rd = CDE_rd, SDE_0_rd = SDE_0_rd, SIE_0_rd = SIE_0_rd,
           SDE_1_rd = SDE_1_rd, SIE_1_rd = SIE_1_rd))
}


## IPW estimator: point estimates from a single long dataset
compute_ipw_estimates <- function(long_data, base_data, n_pat) {
  
  #CDE (M set to 0)
  plrFitM_de <- glm(otherDeath ~ dtime + I(dtime^2) + normalAct + ageCat + hx + hgBinary + rx,
                    data = long_data, family = binomial())
  
  predM <- 1 - predict(plrFitM_de, newdata = long_data, type = "response")
  cumPredM <- unlist(aggregate(predM ~ long_data$patno, FUN = cumprod)$predM, use.names = FALSE)
  weights_cde <- 1 / cumPredM
  
  Treated_cde_hazP <- nonParametricCumHaz(weights_cde, long_data, grp = 1, outcomeProstate = TRUE)
  Placebo_cde_hazP <- nonParametricCumHaz(weights_cde, long_data, grp = 0, outcomeProstate = TRUE)
  hazM_zero <- rep(0, length(cut_times))
  
  cumInc_treated_cde <- nonParametricCumInc(Treated_cde_hazP, hazM_zero)
  cumInc_placebo_cde <- nonParametricCumInc(Placebo_cde_hazP, hazM_zero)
  CDE_rd <- cumInc_treated_cde[length(cut_times)] - cumInc_placebo_cde[length(cut_times)]
  
  #separable effects
  plrFitM_sep <- glm(otherDeath ~ Mrx * (dtime + I(dtime^2) + I(dtime^3)) + normalAct + ageCat + hx + hgBinary,
                     data = long_data, family = binomial())
  
  #counterfactual datasets
  treated <- base_data[rep(1:n_pat, each = length(cut_times)), ]
  treated$dtime <- rep(cut_times, n_pat); treated$rx <- 1; treated$Mrx <- 1
  
  placebo <- base_data[rep(1:n_pat, each = length(cut_times)), ]
  placebo$dtime <- rep(cut_times, n_pat); placebo$rx <- 0; placebo$Mrx <- 0
  
  treated$hazardM <- predict(plrFitM_sep, newdata = treated, type = "response")
  placebo$hazardM <- predict(plrFitM_sep, newdata = placebo, type = "response")
  
  #weight frame
  weight_frame <- data.frame(
    dtime = treated$dtime,
    patno = treated$patno,
    predMx1 = (1 - treated$hazardM),
    predMx0 = (1 - placebo$hazardM)
  )
  
  comp_weights_M <- merge(long_data, weight_frame, by = c("patno", "dtime"), all.x = TRUE, sort = FALSE)
  cum_pred_M_0 <- unlist(aggregate(comp_weights_M$predMx0 ~ comp_weights_M$patno, FUN = cumprod)[[2]], use.names = FALSE)
  cum_pred_M_1 <- unlist(aggregate(comp_weights_M$predMx1 ~ comp_weights_M$patno, FUN = cumprod)[[2]], use.names = FALSE)
  
  #weights for shifted-M counterfactuals
  ipw_sep_0 <- cum_pred_M_1 / cum_pred_M_0
  ipw_sep_1 <- cum_pred_M_0 / cum_pred_M_1
  ipw_1 <- rep(1, nrow(long_data))
  
  #cumulative incidence estimates at month 50
  risk_x1 <- discrete_cuminc_prost(ipw_1, long_data, grp = 1, follow_up = cut_times)[length(cut_times)]
  risk_x0 <- discrete_cuminc_prost(ipw_1, long_data, grp = 0, follow_up = cut_times)[length(cut_times)]
  risk_xY1_xM0 <- discrete_cuminc_prost(ipw_sep_0, long_data, grp = 1, follow_up = cut_times)[length(cut_times)]
  risk_xY0_xM1 <- discrete_cuminc_prost(ipw_sep_1, long_data, grp = 0, follow_up = cut_times)[length(cut_times)]
  
  SDE_0_rd <- risk_xY1_xM0 - risk_x0
  SIE_0_rd <- risk_xY0_xM1 - risk_x0
  SDE_1_rd <- risk_x1 - risk_xY0_xM1
  SIE_1_rd <- risk_x1 - risk_xY1_xM0
  
  return(c(CDE_rd = CDE_rd, SDE_0_rd = SDE_0_rd, SIE_0_rd = SIE_0_rd,
           SDE_1_rd = SDE_1_rd, SIE_1_rd = SIE_1_rd,
           risk_x1 = risk_x1, risk_x0 = risk_x0))
}


## point estimates: g-estimator
g_point <- compute_g_estimates(long_prost_red, baseline, n)

cat("=== Point Estimates at 50 months (Risk Difference) - g-estimator ===\n\n")
cat("CDE(0):", round(g_point["CDE_rd"], 4), "\n")
cat("NDE(0)/SDE(0):", round(g_point["SDE_0_rd"], 4), "\n")
cat("NIE(0)/SIE(0):", round(g_point["SIE_0_rd"], 4), "\n")
cat("NDE(1)/SDE(1):", round(g_point["SDE_1_rd"], 4), "\n")
cat("NIE(1)/SIE(1):", round(g_point["SIE_1_rd"], 4), "\n\n")

cat("Decomposition check: TE = NDE(0)/SDE(0) + NIE(1)/SIE(1) =", round(g_point["SDE_0_rd"] + g_point["SIE_1_rd"], 4), "\n")
cat("Decomposition check: TE = NIE(0)/SIE(0) + NDE(1)/SDE(1) =", round(g_point["SIE_0_rd"] + g_point["SDE_1_rd"], 4), "\n\n")


## point estimates: IPW estimator
ipw_point <- compute_ipw_estimates(long_prost_red, baseline, n)

cat("=== Point Estimates at 50 months (Risk Difference) - IPW estimator ===\n\n")
cat("CDE(0):", round(ipw_point["CDE_rd"], 4), "\n")
cat("NDE(0)/SDE(0):", round(ipw_point["SDE_0_rd"], 4), "\n")
cat("NIE(0)/SIE(0):", round(ipw_point["SIE_0_rd"], 4), "\n")
cat("NDE(1)/SDE(1):", round(ipw_point["SDE_1_rd"], 4), "\n")
cat("NIE(1)/SIE(1):", round(ipw_point["SIE_1_rd"], 4), "\n\n")

cat("Decomposition check: TE = NDE(0)/SDE(0) + NIE(1)/SIE(1) =", round(ipw_point["SDE_0_rd"] + ipw_point["SIE_1_rd"], 4), "\n")
cat("Decomposition check: TE = NIE(0)/SIE(0) + NDE(1)/SDE(1) =", round(ipw_point["SIE_0_rd"] + ipw_point["SDE_1_rd"], 4), "\n")
cat("Total Effect:", round(ipw_point["risk_x1"] - ipw_point["risk_x0"], 4), "\n\n")


## bootstrap confidence intervals (shared loop for both estimators)
set.seed(12345)
nboot <- 1000

#g-estimator
boot_g_CDE_rd <- rep(NA, nboot)
boot_g_SDE_0_rd <- rep(NA, nboot)
boot_g_SIE_0_rd <- rep(NA, nboot)
boot_g_SDE_1_rd <- rep(NA, nboot)
boot_g_SIE_1_rd <- rep(NA, nboot)

#IPW estimator
boot_ipw_CDE_rd <- rep(NA, nboot)
boot_ipw_SDE_0_rd <- rep(NA, nboot)
boot_ipw_SIE_0_rd <- rep(NA, nboot)
boot_ipw_SDE_1_rd <- rep(NA, nboot)
boot_ipw_SIE_1_rd <- rep(NA, nboot)

unique_ids <- unique(prost_red$patno)
n_patients <- length(unique_ids)

for (b in 1:nboot) {
  
  if (b %% 10 == 0) cat("Bootstrap iteration:", b, "\n")
  
  boot_ids <- sample(unique_ids, size = n_patients, replace = TRUE)
  
  boot_data_list <- lapply(seq_along(boot_ids), \(i) {
    temp <- prost_red[prost_red$patno == boot_ids[i], ]
    temp$patno <- i
    temp
  })
  boot_prost_red <- do.call(rbind, boot_data_list)
  
  boot_long <- tryCatch(build_long(boot_prost_red), error = \(e) NULL)
  if (is.null(boot_long)) next
  
  boot_baseline <- boot_long[boot_long$dtime == 0, ]
  boot_n <- length(unique(boot_long$patno))
  
  #g-estimator
  boot_g <- tryCatch(compute_g_estimates(boot_long, boot_baseline, boot_n), error = \(e) NULL)
  if (!is.null(boot_g)) {
    boot_g_CDE_rd[b] <- boot_g["CDE_rd"]
    boot_g_SDE_0_rd[b] <- boot_g["SDE_0_rd"]
    boot_g_SIE_0_rd[b] <- boot_g["SIE_0_rd"]
    boot_g_SDE_1_rd[b] <- boot_g["SDE_1_rd"]
    boot_g_SIE_1_rd[b] <- boot_g["SIE_1_rd"]
  }
  
  #IPW estimator
  boot_ipw <- tryCatch(compute_ipw_estimates(boot_long, boot_baseline, boot_n), error = \(e) NULL)
  if (!is.null(boot_ipw)) {
    boot_ipw_CDE_rd[b] <- boot_ipw["CDE_rd"]
    boot_ipw_SDE_0_rd[b] <- boot_ipw["SDE_0_rd"]
    boot_ipw_SIE_0_rd[b] <- boot_ipw["SIE_0_rd"]
    boot_ipw_SDE_1_rd[b] <- boot_ipw["SDE_1_rd"]
    boot_ipw_SIE_1_rd[b] <- boot_ipw["SIE_1_rd"]
  }
}


## results
alpha <- 0.05
lower_q <- alpha / 2
upper_q <- 1 - alpha / 2

#g-estimator
cat("\nSuccessful bootstrap iterations (g-estimator):", sum(!is.na(boot_g_CDE_rd)), "out of", nboot, "\n\n")

ci_g_CDE_rd <- quantile(boot_g_CDE_rd, probs = c(lower_q, upper_q), na.rm = TRUE)
ci_g_SDE_0_rd <- quantile(boot_g_SDE_0_rd, probs = c(lower_q, upper_q), na.rm = TRUE)
ci_g_SIE_0_rd <- quantile(boot_g_SIE_0_rd, probs = c(lower_q, upper_q), na.rm = TRUE)
ci_g_SDE_1_rd <- quantile(boot_g_SDE_1_rd, probs = c(lower_q, upper_q), na.rm = TRUE)
ci_g_SIE_1_rd <- quantile(boot_g_SIE_1_rd, probs = c(lower_q, upper_q), na.rm = TRUE)

results_g <- data.frame(
  Estimand = c("CDE(0)", "NDE(0)/SDE(0)", "NIE(0)/SIE(0)", "NDE(1)/SDE(1)", "NIE(1)/SIE(1)"),
  Estimate = round(c(g_point["CDE_rd"], g_point["SDE_0_rd"], g_point["SIE_0_rd"], g_point["SDE_1_rd"], g_point["SIE_1_rd"]), 4),
  CI_lower = round(c(ci_g_CDE_rd[1], ci_g_SDE_0_rd[1], ci_g_SIE_0_rd[1], ci_g_SDE_1_rd[1], ci_g_SIE_1_rd[1]), 4),
  CI_upper = round(c(ci_g_CDE_rd[2], ci_g_SDE_0_rd[2], ci_g_SIE_0_rd[2], ci_g_SDE_1_rd[2], ci_g_SIE_1_rd[2]), 4)
)

cat("=== Results at 50 months with 95% Bootstrap CIs (Risk Difference) - g-estimator ===\n\n")
print(knitr::kable(results_g, col.names = c("Estimand", "Estimate", "95% CI Lower", "95% CI Upper")))

#IPW estimator
cat("\n\nSuccessful bootstrap iterations (IPW estimator):", sum(!is.na(boot_ipw_CDE_rd)), "out of", nboot, "\n\n")

ci_ipw_CDE_rd <- quantile(boot_ipw_CDE_rd, probs = c(lower_q, upper_q), na.rm = TRUE)
ci_ipw_SDE_0_rd <- quantile(boot_ipw_SDE_0_rd, probs = c(lower_q, upper_q), na.rm = TRUE)
ci_ipw_SIE_0_rd <- quantile(boot_ipw_SIE_0_rd, probs = c(lower_q, upper_q), na.rm = TRUE)
ci_ipw_SDE_1_rd <- quantile(boot_ipw_SDE_1_rd, probs = c(lower_q, upper_q), na.rm = TRUE)
ci_ipw_SIE_1_rd <- quantile(boot_ipw_SIE_1_rd, probs = c(lower_q, upper_q), na.rm = TRUE)

results_ipw <- data.frame(
  Estimand = c("CDE(0)", "NDE(0)/SDE(0)", "NIE(0)/SIE(0)", "NDE(1)/SDE(1)", "NIE(1)/SIE(1)"),
  Estimate = round(c(ipw_point["CDE_rd"], ipw_point["SDE_0_rd"], ipw_point["SIE_0_rd"], ipw_point["SDE_1_rd"], ipw_point["SIE_1_rd"]), 4),
  CI_lower = round(c(ci_ipw_CDE_rd[1], ci_ipw_SDE_0_rd[1], ci_ipw_SIE_0_rd[1], ci_ipw_SDE_1_rd[1], ci_ipw_SIE_1_rd[1]), 4),
  CI_upper = round(c(ci_ipw_CDE_rd[2], ci_ipw_SDE_0_rd[2], ci_ipw_SIE_0_rd[2], ci_ipw_SDE_1_rd[2], ci_ipw_SIE_1_rd[2]), 4)
)

cat("=== Results at 50 months with 95% Bootstrap CIs (Risk Difference) - IPW estimator ===\n\n")
print(knitr::kable(results_ipw, col.names = c("Estimand", "Estimate", "95% CI Lower", "95% CI Upper")))


#bounds

#for all patients who are even-free at month 50, set dtime to 50 and eventType to 0 (censored)
prostRed$eventType[prostRed$dtime > 50] <- 0
prostRed$dtime[prostRed$dtime > 50] <- 50

table(prostRed$eventType[prostRed$rx==0])/length(prostRed$eventType[prostRed$rx==0])
table(prostRed$eventType[prostRed$rx==1])/length(prostRed$eventType[prostRed$rx==1])

prostRed_treated <- prostRed[prostRed$rx == 1,]
prostRed_placebo <- prostRed[prostRed$rx == 0,]

data_treated_eventAtTime_50 <- prostRed_treated %>%
  group_by(patno) %>%  
  filter(dtime <= 50) %>%
  slice_max(dtime, n = 1, with_ties = FALSE) %>%
  ungroup()

data_placebo_eventAtTime_50 <- prostRed_placebo %>%
  group_by(patno) %>%  
  filter(dtime <= 50) %>%
  slice_max(dtime, n = 1, with_ties = FALSE) %>%
  ungroup()


prostate_p00.0 <- nrow(data_placebo_eventAtTime_50[
  data_placebo_eventAtTime_50$eventType == 0,
]) / nrow(data_placebo_eventAtTime_50)

prostate_p00.1 <- nrow(data_treated_eventAtTime_50[
  data_treated_eventAtTime_50$eventType == 0,
]) / nrow(data_treated_eventAtTime_50)

prostate_p01.0 <- nrow(data_placebo_eventAtTime_50[
  data_placebo_eventAtTime_50$eventType == 1,
]) / nrow(data_placebo_eventAtTime_50)

prostate_p01.1 <- nrow(data_treated_eventAtTime_50[
  data_treated_eventAtTime_50$eventType == 1,
]) / nrow(data_treated_eventAtTime_50)

prostate_p10.0 <- nrow(data_placebo_eventAtTime_50[
  data_placebo_eventAtTime_50$eventType == 2,
]) / nrow(data_placebo_eventAtTime_50)

prostate_p10.1 <- nrow(data_treated_eventAtTime_50[
  data_treated_eventAtTime_50$eventType == 2,
]) / nrow(data_treated_eventAtTime_50)


prostate_p00.0 + prostate_p01.0 + prostate_p10.0
prostate_p00.1 + prostate_p01.1 + prostate_p10.1

#CDE(0)
prostate_p00.0 + prostate_p01.1 - 1
1 - prostate_p00.1 - prostate_p01.0 

#NDE(0)
max((-1 + prostate_p00.0 + prostate_p01.1), (-prostate_p01.0))
min(prostate_p00.0, (1 - prostate_p00.1 - prostate_p01.0))

#NDE(1)
max((-prostate_p00.1), (-1 + prostate_p00.0 + prostate_p01.1))
min(prostate_p01.1, (1 - prostate_p00.1 - prostate_p01.0))

#NIE(0)
max((-1 + prostate_p00.1 + prostate_p01.1), (-prostate_p01.0))
min((1 - prostate_p00.0 - prostate_p01.0), (prostate_p00.1 + prostate_p01.1 - prostate_p01.0))

#NIE(1)
max((-prostate_p00.0 - prostate_p01.0 + prostate_p01.1), (-1 + prostate_p00.1 + prostate_p01.1))
min((1 - prostate_p00.0 - prostate_p01.0), prostate_p01.1)

#Frechet bounds for the NDE
max(0, nrow(prostRed[prostRed$rx == 1 & prostRed$eventType == 1,]) / nrow(prostRed[prostRed$rx == 1 & prostRed$eventType != 2,]) + nrow(prostRed[prostRed$rx == 0 & prostRed$eventType != 2,]) / nrow(prostRed[prostRed$rx == 0,]) -1) - nrow(prostRed[prostRed$rx == 0 & prostRed$eventType == 1,]) / nrow(prostRed[prostRed$rx == 0,]) 
min(nrow(prostRed[prostRed$rx == 1 & prostRed$eventType == 1,]) / nrow(prostRed[prostRed$rx == 1 & prostRed$eventType != 2,]), nrow(prostRed[prostRed$rx == 0 & prostRed$eventType != 2,]) / nrow(prostRed[prostRed$rx == 0,]) - nrow(prostRed[prostRed$rx == 0 & prostRed$eventType == 1,]) / nrow(prostRed[prostRed$rx == 0,])) 



