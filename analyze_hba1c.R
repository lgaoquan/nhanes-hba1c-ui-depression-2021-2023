suppressPackageStartupMessages({library(haven); library(survey); library(broom)})
options(survey.lonely.psu = "adjust", survey.adjust.domain.lonely = TRUE)
args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) normalizePath(args[[1]], winslash = "/") else normalizePath(getwd(), winslash = "/")
out_dir <- file.path(root, "results"); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
stems <- c("DEMO_L", "KIQ_U_L", "DPQ_L", "GHB_L", "BMX_L", "DIQ_L", "HIQ_L")
dat <- lapply(stems, function(x) read_xpt(file.path(root, paste0(x, ".xpt")))); names(dat) <- stems
stopifnot(all(vapply(dat, function(x) !anyDuplicated(x$SEQN), logical(1))))
d <- Reduce(function(a,b) merge(a,b,by="SEQN",all.x=TRUE,sort=FALSE), dat)
stopifnot(nrow(d) == nrow(dat$DEMO_L), !anyDuplicated(d$SEQN))
phq_items <- sprintf("DPQ%03d", seq(10,90,10))
phq <- as.data.frame(lapply(d[phq_items], function(x) {x <- as.numeric(x); ifelse(x %in% 0:3,x,NA_real_)}))
d$phq_complete <- rowSums(!is.na(phq)) == length(phq_items)
d$phq9 <- ifelse(d$phq_complete,rowSums(phq),NA_real_)
d$depressive_symptoms <- ifelse(d$phq_complete,as.integer(d$phq9 >= 10),NA_integer_)
d$ui_current <- d$KIQ005 %in% 2:5
d$hba1c_valid <- is.finite(d$LBXGH) & d$LBXGH >= 3 & d$LBXGH <= 20
d$design_ok <- d$RIAGENDR == 2 & d$RIDAGEYR >= 20 & is.finite(d$WTPH2YR) & d$WTPH2YR > 0 & !is.na(d$SDMVSTRA) & !is.na(d$SDMVPSU)
d$diabetes_status <- factor(ifelse(d$DIQ010==1,"Diagnosed diabetes",ifelse(d$DIQ010==2,"No diagnosed diabetes",ifelse(d$DIQ010==3,"Borderline diabetes",NA_character_))), levels=c("No diagnosed diabetes","Borderline diabetes","Diagnosed diabetes"))
d$race_group <- factor(ifelse(d$RIDRETH3==3,"Non-Hispanic White",ifelse(d$RIDRETH3==4,"Non-Hispanic Black",ifelse(d$RIDRETH3 %in% c(1,2),"Hispanic","Other/multiracial"))), levels=c("Non-Hispanic White","Hispanic","Non-Hispanic Black","Other/multiracial"))
d$education_group <- factor(ifelse(d$DMDEDUC2 %in% 1:3,"High school or less",ifelse(d$DMDEDUC2==4,"Some college/AA",ifelse(d$DMDEDUC2==5,"College or above",NA_character_))), levels=c("College or above","Some college/AA","High school or less"))
d$marital_group <- factor(ifelse(d$DMDMARTZ==1,"Married/partnered",ifelse(d$DMDMARTZ %in% 2:3,"Not married/partnered",NA_character_)), levels=c("Married/partnered","Not married/partnered"))
d$hba1c_group <- cut(d$LBXGH, breaks=c(-Inf,5.7,6.5,Inf), right=FALSE, labels=c("<5.7%","5.7-<6.5%",">=6.5%"))
d$bmi <- as.numeric(d$BMXBMI)
analysis <- d$design_ok & d$ui_current & d$phq_complete & d$hba1c_valid
complete_adj <- analysis & complete.cases(d[,c("RIDAGEYR","race_group","bmi","diabetes_status")])
d$analysis <- analysis; d$complete_adj <- complete_adj
sample_flow <- data.frame(step=c("NHANES records","Women aged >=20","Women aged >=20 with valid exam design and phlebotomy weight","Current urinary leakage","Complete PHQ-9 and valid HbA1c","Complete prespecified adjustment covariates"), n=c(nrow(d),sum(d$RIAGENDR==2 & d$RIDAGEYR>=20,na.rm=TRUE),sum(d$design_ok,na.rm=TRUE),sum(d$design_ok & d$ui_current,na.rm=TRUE),sum(analysis,na.rm=TRUE),sum(complete_adj,na.rm=TRUE)))
write.csv(sample_flow,file.path(out_dir,"sample_flow.csv"),row.names=FALSE,na="")
design_full <- svydesign(ids=~SDMVPSU,strata=~SDMVSTRA,weights=~WTPH2YR,nest=TRUE,data=d[d$design_ok,])
design_main <- subset(design_full,analysis); design_adj <- subset(design_full,complete_adj)
stopifnot(nrow(design_main$variables) > 0)
model_cont <- svyglm(depressive_symptoms ~ LBXGH,design=design_main,family=quasipoisson(link="log"))
model_cont_adj <- svyglm(depressive_symptoms ~ LBXGH + RIDAGEYR + race_group + bmi + diabetes_status,design=design_adj,family=quasipoisson(link="log"))
model_cat <- svyglm(depressive_symptoms ~ hba1c_group,design=design_main,family=quasipoisson(link="log"))
model_cat_adj <- svyglm(depressive_symptoms ~ hba1c_group + RIDAGEYR + race_group + bmi + diabetes_status,design=design_adj,family=quasipoisson(link="log"))
model_phq_linear <- svyglm(phq9 ~ LBXGH,design=design_main,family=gaussian())
model_phq_linear_adj <- svyglm(phq9 ~ LBXGH + RIDAGEYR + race_group + bmi + diabetes_status,design=design_adj,family=gaussian())
tidy_model <- function(m,name,exponentiate=TRUE) {out <- broom::tidy(m,conf.int=TRUE,exponentiate=exponentiate); out$model <- name; out$design_df <- degf(m$survey.design); out}
model_results <- rbind(tidy_model(model_cont,"Modified Poisson, continuous HbA1c, unadjusted"),tidy_model(model_cont_adj,"Modified Poisson, continuous HbA1c, adjusted"),tidy_model(model_cat,"Modified Poisson, categorical HbA1c, unadjusted"),tidy_model(model_cat_adj,"Modified Poisson, categorical HbA1c, adjusted"))
write.csv(model_results,file.path(out_dir,"model_results_prevalence_ratios.csv"),row.names=FALSE,na="")
linear_results <- rbind(tidy_model(model_phq_linear,"PHQ-9 linear model, unadjusted",FALSE),tidy_model(model_phq_linear_adj,"PHQ-9 linear model, adjusted",FALSE))
write.csv(linear_results,file.path(out_dir,"model_results_phq9_linear.csv"),row.names=FALSE,na="")
weighted_mean <- svymean(~LBXGH + phq9,design_main,na.rm=TRUE)
weighted_props <- svymean(~depressive_symptoms,design_main,na.rm=TRUE)
weighted_prop_ci <- svyciprop(~I(depressive_symptoms == 1), design_main, method="logit")
summary_results <- data.frame(estimand=c("Unweighted analytic sample","Unweighted elevated PHQ-9 count","Weighted mean HbA1c (%)","Weighted mean PHQ-9","Weighted prevalence PHQ-9 >=10"), estimate=c(nrow(design_main$variables),sum(design_main$variables$depressive_symptoms==1,na.rm=TRUE),coef(weighted_mean)["LBXGH"],coef(weighted_mean)["phq9"],coef(weighted_props)["depressive_symptoms"]), se=c(NA,NA,SE(weighted_mean)["LBXGH"],SE(weighted_mean)["phq9"],SE(weighted_props)["depressive_symptoms"]), design_df=c(degf(design_main),NA,degf(design_main),degf(design_main),degf(design_main)))
summary_results$ci_low <- NA_real_; summary_results$ci_high <- NA_real_
summary_results$ci_low[5] <- confint(weighted_prop_ci)[1]; summary_results$ci_high[5] <- confint(weighted_prop_ci)[2]
write.csv(summary_results,file.path(out_dir,"summary_results.csv"),row.names=FALSE,na="")
table1_rows <- list()
for (g in c(0,1)) {
  dg <- subset(design_main, depressive_symptoms == g)
  gl <- if (g == 1) "PHQ-9 >=10" else "PHQ-9 <10"
  add_mean <- function(v, label) {
    z <- svymean(as.formula(paste0("~",v)), dg, na.rm=TRUE)
    data.frame(characteristic=label, level="Mean (SE)", group=gl, estimate=as.numeric(coef(z)[1]), se=as.numeric(SE(z)[1]))
  }
  table1_rows[[length(table1_rows)+1]] <- add_mean("RIDAGEYR", "Age, years")
  table1_rows[[length(table1_rows)+1]] <- add_mean("bmi", "BMI, kg/m2")
  table1_rows[[length(table1_rows)+1]] <- add_mean("LBXGH", "HbA1c, %")
  table1_rows[[length(table1_rows)+1]] <- add_mean("phq9", "PHQ-9 total score")
  for (v in c("race_group","diabetes_status","hba1c_group")) {
    z <- svymean(as.formula(paste0("~",v)), dg, na.rm=TRUE)
    for (nm in names(coef(z))) table1_rows[[length(table1_rows)+1]] <- data.frame(characteristic=v, level=nm, group=gl, estimate=as.numeric(coef(z)[nm]), se=as.numeric(SE(z)[nm]))
  }
}
write.csv(do.call(rbind, table1_rows), file.path(out_dir,"table1_weighted_characteristics.csv"), row.names=FALSE, na="")
models <- list(model_cont,model_cont_adj,model_cat,model_cat_adj,model_phq_linear,model_phq_linear_adj)
diagnostics <- data.frame(model=c("cont_unadjusted","cont_adjusted","cat_unadjusted","cat_adjusted","phq_linear_unadjusted","phq_linear_adjusted"), finite_coefficients=vapply(models,function(m) all(is.finite(coef(m))),logical(1)), finite_se=vapply(models,function(m) all(is.finite(sqrt(diag(vcov(m))))),logical(1)), residual_df=vapply(models,function(m) m$df.residual,numeric(1)), design_df=vapply(models,function(m) degf(m$survey.design),numeric(1)))
write.csv(diagnostics,file.path(out_dir,"model_diagnostics.csv"),row.names=FALSE,na="")
writeLines(capture.output(sessionInfo()),file.path(out_dir,"session_info.txt"))
writeLines(capture.output(str(list(data_files=stems,exposure="LBXGH with WTPH2YR",outcome="PHQ-9 total and >=10",design="SDMVPSU, SDMVSTRA, WTPH2YR",analysis_n=nrow(design_main$variables),adjustment_n=nrow(design_adj$variables)))),file.path(out_dir,"analysis_provenance.txt"))
message("Completed independent HbA1c analysis. Main n=",nrow(design_main$variables),"; adjusted n=",nrow(design_adj$variables),"; design df=",degf(design_main))
