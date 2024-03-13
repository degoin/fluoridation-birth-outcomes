library(tidyverse)
library(readxl)
library(splines)
library(mgcViz)

# fluoride and birth outcomes analysis file 
ls <- list()

for (yr in 2000:2018) {
  print(yr)
ls[[yr-1999]] <- read.csv(paste0(".../data/analytic_file_birth_fluoride_",yr,".csv"))

 }

df <- data.frame(do.call(rbind, ls))

df$macrosomia <- ifelse(df$birthweight_g>=4000,1,0)

# complete cases 
df0 <- df[complete.cases(df),]

# keep only singletons 
df0 <- df0 %>% filter(mult_births==1)


# overall analyses 

o_result <- function(outcome) {
  
  df0$y <- df0[[outcome]]

m2 <- glm(y ~ fluoride10 + mat_age + black_nh + asian_nh + amerind_nh + 
            hipacis_nh + other_nh + multi + hispanic_single + 
            lths + hs_grad + some_col + col_grad + govt_ins_pay + other_pay +
            Income_Ratio + unemp_rate + avg_temp + Total_pop +   
            factor(urban_code) + factor(con_month) + factor(con_year) + 
          avg.finding_Arsenic,
          data =df0)


vars <- vcovHC(m2)

result <- data.frame(est = coef(m2)["fluoride10"], 
                     lb = summary(m2)$coefficients["fluoride10","Estimate"] - 1.96*sqrt(vars["fluoride10","fluoride10"]), 
                     ub = summary(m2)$coefficients["fluoride10","Estimate"] + 1.96*sqrt(vars["fluoride10","fluoride10"])) 
                     

return(result) 
}


outcomes <- c("birthweight_g", "est_gest_wks","gzscore", "ptb", "sga", "lga", "macrosomia")

out_ls <- list()
for (i in 1:length(outcomes)) {
  print(outcomes[i])
 out_ls[[i]] <- o_result(outcomes[i]) 
}

df_o <- do.call(rbind, out_ls)

rownames(df_o) <- outcomes

write.csv(df_o, file=".../results/fluoride_overall_associations.csv")

# compare to sibling analysis 

# fluoride and birth outcomes analysis file 
ls <- list()

# only have sibling data from 2001-2020
for (yr in 2001:2018) {
  print(yr)
  ls[[yr-1999]] <- read.csv(paste0(".../data/analytic_file_birth_fluoride_",yr,".csv"))
  
}

df <- data.frame(do.call(rbind, ls))
df <- df %>% mutate(SFN_t = paste0("'", as.character(SFN_1)))

load(".../data/sibling-data/sibling-ids-01-20")

df <- left_join(df, df_s)

df_sibs <- df %>% filter(nsibs>1 & nsibs<1000)
dim(df_sibs)
table(df_sibs$nsibs)


df_sibs$macrosomia <- ifelse(df_sibs$birthweight_g>=4000,1,0)

# complete cases 
df0 <- df_sibs[complete.cases(df_sibs),]
dim(df0)
dim(df_sibs)

# keep only singletons 
df0 <- df0 %>% filter(mult_births==1)

dim(df0)

# verify you still have at least 2 births per mother 

df0 <- df0 %>% group_by(MotherID) %>% arrange(con_year, con_month) %>% mutate(nsibs = n(), sib_order = row_number()) %>% mutate(flag = as.numeric(nsibs==1))
table(df0$nsibs)

# because there are a lot missing fluoride and covariates, need to re-restrict to multiple births with exposure info
df0 <- df0 %>% filter(nsibs>1)

dim(df0)

# create data frame that includes maternal fixed effects, ie. subtract the mean across siblings for each variable

df_s1 <- df0 %>% group_by(MotherID) %>% 
  mutate(fluoride10_diff1 = fluoride10 -mean(fluoride10), 
         mfluoride10 = mean(fluoride10),
         mat_age_diff1 = mat_age - mean(mat_age), 
         lths_diff1 = lths - mean(lths), 
         hs_grad_diff1 = hs_grad - mean(hs_grad), 
         some_col_diff1 = some_col - mean(some_col), 
         col_grad_diff1 = col_grad - mean(col_grad), 
         grad_grad_diff1 = grad_grad - mean(grad_grad),
         govt_ins_pay_diff1 = govt_ins_pay - mean(govt_ins_pay), 
         private_ins_pay_diff1 = private_ins_pay - mean(private_ins_pay),
         other_pay_diff1 = other_pay - mean(other_pay),
         Income_Ratio_diff1 = Income_Ratio - mean(Income_Ratio),
         unemp_rate_diff1 = unemp_rate - mean(unemp_rate),
         avg_temp_diff1 = avg_temp - mean(avg_temp),
         Total_pop_diff1 = Total_pop - mean(Total_pop),
         avg.finding_Arsenic_diff1 = avg.finding_Arsenic - mean(avg.finding_Arsenic), 
         bw_diff1 = birthweight_g - mean(birthweight_g), 
         gz_diff1 = gzscore - mean(gzscore),
         est_gest_wks_diff1 = est_gest_wks - mean(est_gest_wks),
         ptb_diff1 = ptb - mean(ptb),
         sga_diff1 = sga - mean(sga),
         lga_diff1 = lga - mean(lga), 
         mac_diff1 = macrosomia - mean(macrosomia)) 

# sibling analyses 

s_result <- function(outcome) {
  
  df_s1$y <- df_s1[[outcome]]
  
  m2 <- glm(y ~ fluoride10_diff1 + mat_age_diff1 + 
              lths_diff1 + hs_grad_diff1 + some_col_diff1 + col_grad_diff1 + govt_ins_pay_diff1 + other_pay_diff1 +
              Income_Ratio_diff1 + unemp_rate_diff1 + avg_temp_diff1 + Total_pop_diff1 +   
              factor(urban_code) + factor(con_month) + factor(con_year) + 
              avg.finding_Arsenic_diff1,
            data =df_s1)
  
  
  vars <- vcovHC(m2)
  
  result <- data.frame(est = coef(m2)["fluoride10_diff1"], 
                       lb = summary(m2)$coefficients["fluoride10_diff1","Estimate"] - 1.96*sqrt(vars["fluoride10_diff1","fluoride10_diff1"]), 
                       ub = summary(m2)$coefficients["fluoride10_diff1","Estimate"] + 1.96*sqrt(vars["fluoride10_diff1","fluoride10_diff1"])) 
  
  
  return(result) 
}


outcomes <- c("bw_diff1", "est_gest_wks_diff1","gz_diff1", "ptb_diff1", "sga_diff1", "lga_diff1", "mac_diff1")

out_ls <- list()
for (i in 1:length(outcomes)) {
  print(outcomes[i])
  out_ls[[i]] <- s_result(outcomes[i]) 
}

df_os <- do.call(rbind, out_ls)


rownames(df_os) <- outcomes

write.csv(df_os, file=".../results/fluoride_sibling_associations.csv")


# unmatched sibling results 

s2_result <- function(outcome) {
  
  df_s1$y <- df_s1[[outcome]]
  
  m2 <- glm(y ~ fluoride10 + mat_age + black_nh + asian_nh + amerind_nh + 
              hipacis_nh + other_nh + multi + hispanic_single + 
              lths + hs_grad + some_col + col_grad + govt_ins_pay + other_pay +
              Income_Ratio + unemp_rate + avg_temp + Total_pop +   
              factor(urban_code) + factor(con_month) + factor(con_year) + 
              avg.finding_Arsenic,
            data =df_s1)
  
  
  vars <- vcovHC(m2)
  
  result <- data.frame(est = coef(m2)["fluoride10"], 
                       lb = summary(m2)$coefficients["fluoride10","Estimate"] - 1.96*sqrt(vars["fluoride10","fluoride10"]), 
                       ub = summary(m2)$coefficients["fluoride10","Estimate"] + 1.96*sqrt(vars["fluoride10","fluoride10"])) 
  

  return(result) 
}


outcomes <- c("birthweight_g", "est_gest_wks","gzscore", "ptb", "sga", "lga", "macrosomia")

out_ls <- list()
for (i in 1:length(outcomes)) {
  print(outcomes[i])
  out_ls[[i]] <- s2_result(outcomes[i]) 
}

df_os2 <- do.call(rbind, out_ls)


rownames(df_os2) <- outcomes

write.csv(df_os2, file=".../results/fluoride_sibling_associations_unmatched.csv")



