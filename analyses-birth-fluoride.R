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

# keep complete cases
df0 <- df[complete.cases(df),]

# keep only singletons 
df0 <- df0 %>% filter(mult_births==1)

# differences in fluoride levels
summary(df0$fluoride_ppm)

# by race
summary(df0$fluoride_ppm[df0$white_nh==1])
summary(df0$fluoride_ppm[df0$black_nh==1])

summary(df0$fluoride_ppm[df0$asian_nh==1])
summary(df0$fluoride_ppm[df0$amerind_nh==1])
summary(df0$fluoride_ppm[df0$hipacis_nh==1])
summary(df0$fluoride_ppm[df0$other_nh==1])
summary(df0$fluoride_ppm[df0$multi==1])
summary(df0$fluoride_ppm[df0$hispanic_single==1])

table(df0$white_nh)
table(df0$black_nh)
table(df0$asian_nh)
table(df0$amerind_nh)
table(df0$hipacis_nh)
table(df0$other_nh)
table(df0$multi)
table(df0$hispanic_single)

table(df0$lths)
summary(df0$fluoride_ppm[df0$lths==1])

table(df0$hs_grad)
summary(df0$fluoride_ppm[df0$hs_grad==1])

table(df0$some_col)
summary(df0$fluoride_ppm[df0$some_col==1])

table(df0$col_grad)
summary(df0$fluoride_ppm[df0$col_grad==1])

table(df0$grad_grad)
summary(df0$fluoride_ppm[df0$grad_grad==1])


table(df0$govt_ins_pay)
summary(df0$fluoride_ppm[df0$govt_ins_pay==1])

table(df0$private_ins_pay)
summary(df0$fluoride_ppm[df0$private_ins_pay==1])

table(df0$other_pay)
summary(df0$fluoride_ppm[df0$other_pay==1])

summary(df0$avg.finding_Arsenic)
quantile(df0$avg.finding_Arsenic, 0.5)

df0$arsenic_mcl <- ifelse(df0$avg.finding_Arsenic>=10,1,0)
df0$arsenic_75q <- ifelse(df0$avg.finding_Arsenic>=quantile(df0$avg.finding_Arsenic,0.75),1,0)

df0$macrosomia <- ifelse(df0$birthweight_g>=4000,1,0)


# check distribution of birthweight at earlier gestational ages 
hist(df0$birthweight_g[df0$est_gest_wks==24])

# identify those with implausible combinations of birthweight and gestational age 
df0$implausible_bw_ga <- ifelse(df0$est_gest_wks>36 & df0$gzscore>5, 1, 
                                ifelse(df0$est_gest_wks>36 & df0$gzscore < -5, 1, 
                                ifelse(df0$est_gest_wks<37 & df0$gzscore < -4, 1, 
                                       ifelse(df0$est_gest_wks<37 & df0$gzscore > 3, 1,0 ))))

table(df0$implausible_bw_ga)
prop.table(table(df0$implausible_bw_ga))

hist(df0$birthweight_g[df0$est_gest_wks==24 & df0$implausible_bw_ga==0])

df0b <- df0 %>% filter(implausible_bw_ga==0)
dim(df0b)
dim(df0) - dim(df0b)

# non-linear analyses 
# ----- birthweight 
gam1 <- gam(birthweight_g ~ s(fluoride_ppm) + mat_age + black_nh + asian_nh + amerind_nh + 
             hipacis_nh + other_nh + multi + hispanic_single + 
             lths + hs_grad + some_col + col_grad + govt_ins_pay + other_pay +
             Income_Ratio + unemp_rate + avg_temp + Total_pop +   
             factor(urban_code) + factor(con_month) + factor(con_year) + 
             avg.finding_Arsenic,
           data =df0, method="REML")

b1 <- getViz(gam1)

o <- plot( sm(b1, 1) )
o + l_fitLine(colour = "red") +
  l_ciLine(mul = 5, colour = "blue", linetype = 2) + labs(x="fluoride", y="birthweight") 
  theme_classic()

pdf(".../results/gam_plot_bw.pdf")
o + l_fitLine(colour = "red") +
  l_ciLine(mul = 5, colour = "blue", linetype = 2) +  labs(x="fluoride (ppm)", y="change in birthweight (g)") 

dev.off()


# ----- gestational age 
gam2 <- gam(est_gest_wks ~ s(fluoride_ppm) + mat_age + black_nh + asian_nh + amerind_nh + 
              hipacis_nh + other_nh + multi + hispanic_single + 
              lths + hs_grad + some_col + col_grad + govt_ins_pay + other_pay +
              Income_Ratio + unemp_rate + avg_temp + Total_pop +   
              factor(urban_code) + factor(con_month) + factor(con_year) + 
              avg.finding_Arsenic,
            data =df0, method="REML")

b2 <- getViz(gam2)

o2 <- plot( sm(b2, 1) )


pdf(".../results/gam_plot_ga.pdf")
o2 + l_fitLine(colour = "red") +
  l_ciLine(mul = 5, colour = "blue", linetype = 2) +  labs(x="fluoride (ppm)", y="change in gestational age (wks)") 

dev.off()


# ----- birthweight for gestational age z score
gam3 <- gam(gzscore ~ s(fluoride_ppm) + mat_age + black_nh + asian_nh + amerind_nh + 
              hipacis_nh + other_nh + multi + hispanic_single + 
              lths + hs_grad + some_col + col_grad + govt_ins_pay + other_pay +
              Income_Ratio + unemp_rate + avg_temp + Total_pop +   
              factor(urban_code) + factor(con_month) + factor(con_year) + 
              avg.finding_Arsenic,
            data =df0, method="REML")

b3 <- getViz(gam3)

o3 <- plot( sm(b3, 1) )


pdf(".../results/gam_plot_gz.pdf")
o3 + l_fitLine(colour = "red") +
  l_ciLine(mul = 5, colour = "blue", linetype = 2) +  labs(x="fluoride (ppm)", y="change in z-score") 

dev.off()





# main intervention analyses 

o_result <- function(outcome, quantile, os_df) {
  
  df0$y <- df0[[outcome]]

m2 <- glm(y ~ ns(fluoride_ppm, df=os_df) + mat_age + black_nh + asian_nh + amerind_nh + 
            hipacis_nh + other_nh + multi + hispanic_single + 
            lths + hs_grad + some_col + col_grad + govt_ins_pay + other_pay +
            Income_Ratio + unemp_rate + avg_temp + Total_pop +   
            factor(urban_code) + factor(con_month) + factor(con_year) + 
          avg.finding_Arsenic,
          data =df0)


df1 <- df0
df1$fluoride_ppm <- ifelse(df1$fluoride_ppm> quantile, quantile, df1$fluoride_ppm)

p0 <- predict(m2)
p1 <- predict(m2, newdata = df1)

p_result <- mean(p1-p0)

b_ls <- lapply(1:200, function(x) f_boot(x, outcome = outcome, quantile = quantile, os_df = os_df))
b_df <- do.call(rbind, b_ls)

result <- data.frame(est = p_result, lb = p_result - 1.96*sqrt(var(b_df)), ub = p_result + 1.96*sqrt(var(b_df)))
colnames(result) <- c("est", "lb", "ub")

return(result) 
}



# bootstrap function


f_boot <- function(x, outcome, quantile, os_df) {
  print(x)
  
  water_systems <- unique(df0$System.ID)
  
  resample_t <- sample(water_systems, length(water_systems), replace=T)
  
  boot_ls <- list()
  for (i in 1:length(resample_t)) {
    boot_ls[[i]] <- df0 %>% filter(System.ID==resample_t[i])
  }
  boot_df <- data.frame(do.call(rbind, boot_ls)) 
  
  boot_df$y <- boot_df[[outcome]]
  
  
  b_fit <- glm(y ~ ns(fluoride_ppm, df=os_df) + mat_age + black_nh + asian_nh + amerind_nh + 
                 hipacis_nh + other_nh + multi + hispanic_single + 
                 lths + hs_grad + some_col + col_grad + govt_ins_pay + other_pay +
                 Income_Ratio + unemp_rate + avg_temp + Total_pop +
                 factor(urban_code) + factor(con_month) + factor(con_year) +
                 avg.finding_Arsenic,
               data =boot_df)
  
  boot_df1 <- boot_df
  boot_df1$fluoride_ppm <- ifelse(boot_df1$fluoride_ppm> quantile, quantile, boot_df1$fluoride_ppm)
  
  p0 <- predict(b_fit)
  p1 <- predict(b_fit, newdata = boot_df1)
  
  p_bw <- p1-p0
  
  b_result <- mean(p_bw)
  return(b_result)
}

set.seed(620573585)


ma_07 <- o_result(outcome = "macrosomia", quantile = 0.7, os_df=5)
ma_07
write.csv(ma_07, file=".../results/ma_07.csv", row.names = F)

ma_05 <- o_result(outcome = "macrosomia", quantile = 0.5, os_df=5)
ma_05
write.csv(ma_05, file=".../results/ma_05.csv", row.names = F)


bw_07 <- o_result(outcome = "birthweight_g", quantile = 0.7, os_df=5)
bw_07
write.csv(bw_07, file=".../results/bw_07.csv", row.names = F)

bw_05 <- o_result(outcome = "birthweight_g", quantile = 0.5, os_df=5)
bw_05
write.csv(bw_05, file=".../results/bw_05.csv", row.names = F)

ga_07 <- o_result(outcome = "est_gest_wks", quantile = 0.7, os_df=10)
ga_07
write.csv(ga_07, file=".../results/ga_07.csv", row.names = F)

ga_05 <- o_result(outcome = "est_gest_wks", quantile = 0.5, os_df=10)
ga_05
write.csv(ga_05, file=".../results/ga_05.csv", row.names = F)


bz_07 <- o_result(outcome = "gzscore", quantile = 0.7, os_df=7)
bz_07
write.csv(bz_07, file=".../results/bz_07.csv", row.names = F)

bz_05 <- o_result(outcome = "gzscore", quantile = 0.5, os_df=7)
bz_05
write.csv(bz_05, file=".../results/bz_05.csv", row.names = F)


ptb_07 <- o_result(outcome = "ptb", quantile = 0.7, os_df=10)
ptb_07
write.csv(ptb_07, file=".../results/ptb_07.csv", row.names = F)

ptb_05 <- o_result(outcome = "ptb", quantile = 0.5, os_df=10)
ptb_05
write.csv(ptb_05, file=".../results/ptb_05.csv", row.names = F)


sga_07 <- o_result(outcome = "sga", quantile = 0.7, os_df=7)
sga_07
write.csv(sga_07, file=".../results/sga_07.csv", row.names = F)

sga_05 <- o_result(outcome = "sga", quantile = 0.5, os_df=7)
sga_05
write.csv(sga_05, file=".../results/sga_05.csv", row.names = F)


lga_07 <- o_result(outcome = "lga", quantile = 0.7, os_df=7)
lga_07
write.csv(lga_07, file=".../results/lga_07.csv", row.names = F)

lga_05 <- o_result(outcome = "lga", quantile = 0.5, os_df=7)
lga_05
write.csv(lga_05, file=".../results/lga_05.csv", row.names = F)


# intervention analyses by race/ethnicity

df0$race_eth <- ifelse(df0$white_nh==1, "white", 
                       ifelse(df0$black_nh==1, "black", 
                              ifelse(df0$asian_nh==1, "asian", 
                                     ifelse(df0$amerind_nh==1, "american_indian", 
                                            ifelse(df0$hipacis_nh==1, "hawaiian_pacific_islander", 
                                                   ifelse(df0$other_nh==1, "other", 
                                                          ifelse(df0$multi==1, "multirace", ifelse(df0$hispanic_single==1, "hispanic_latinx", NA))))))))

r_result <- function(outcome, quantile, os_df, race_cat) {
  
  df0$y <- df0[[outcome]]
  
  df_r <- df0 %>% filter(race_eth==race_cat)
  
  m2 <- glm(y ~ ns(fluoride_ppm, df=os_df) + mat_age + 
              lths + hs_grad + some_col + col_grad + govt_ins_pay + other_pay +
              Income_Ratio + unemp_rate + avg_temp + Total_pop +   
              factor(urban_code) + factor(con_month) + factor(con_year) + 
              avg.finding_Arsenic,
            data =df_r)
  
  
  
  df1 <- df_r
  df1$fluoride_ppm <- ifelse(df1$fluoride_ppm> quantile, quantile, df1$fluoride_ppm)
  
  p0 <- predict(m2)
  p1 <- predict(m2, newdata = df1)
  
  p_result <- mean(p1-p0)
  
  b_ls <- lapply(1:200, function(x) r_boot(x, outcome = outcome, quantile = quantile, os_df = os_df, race_cat=race_cat))
  b_df <- do.call(rbind, b_ls)
  
  result <- data.frame(pop = race_cat, est = p_result, lb = p_result - 1.96*sqrt(var(b_df)), ub = p_result + 1.96*sqrt(var(b_df)))
  colnames(result) <- c("pop", "est", "lb", "ub")
  
  return(result) 
}



# bootstrap function by race


r_boot <- function(x, outcome, quantile, os_df, race_cat) {
  print(x)
  
  df_r <- df0 %>% filter(race_eth==race_cat)
  water_systems <- unique(df_r$System.ID)
  
  resample_t <- sample(water_systems, length(water_systems), replace=T)
  
  boot_ls <- list()
  for (i in 1:length(resample_t)) {
    boot_ls[[i]] <- df_r %>% filter(System.ID==resample_t[i])
  }
  boot_df <- data.frame(do.call(rbind, boot_ls)) 
  
  boot_df$y <- boot_df[[outcome]]
  
  
  b_fit <- glm(y ~ ns(fluoride_ppm, df=os_df) + mat_age + 
                 lths + hs_grad + some_col + col_grad + govt_ins_pay + other_pay +
                 Income_Ratio + unemp_rate + avg_temp + Total_pop +
                 factor(urban_code) + factor(con_month) + factor(con_year) +
                 avg.finding_Arsenic,
               data =boot_df)
  
  boot_df1 <- boot_df
  boot_df1$fluoride_ppm <- ifelse(boot_df1$fluoride_ppm> quantile, quantile, boot_df1$fluoride_ppm)
  
  p0 <- predict(b_fit)
  p1 <- predict(b_fit, newdata = boot_df1)
  
  p_bw <- p1-p0
  
  b_result <- mean(p_bw)
  return(b_result)
}

set.seed(620573585)

# race_eth categories
# 1: white 
# 2: black
# 3: asian
# 4: american indian
# 5: hawaiian or pacific islander
# 6: other 
# 7: muli-race
# 8: hispanic/latino single race 

race_groups <- c("white", "black",  "asian",  "american_indian", "hawaiian_pacific_islander", "other",  "multirace", "hispanic_latinx")


ma_r <- list()
for (r in 1:length(race_groups)) {
  ma_r[[r]] <- r_result(outcome = "macrosomia", quantile = 0.7, os_df=5, race_cat=race_groups[r])
}

ma_df_r_07 <- do.call(rbind, ma_r)

write.csv(ma_df_r_07, file=".../results/ma_07_by_race.csv", row.names = F)



ma_r <- list()
for (r in 1:length(race_groups)) {
  ma_r[[r]] <- r_result(outcome = "macrosomia", quantile = 0.5, os_df=5, race_cat=race_groups[r])
}

ma_df_r_05 <- do.call(rbind, ma_r)

write.csv(ma_df_r_05, file=".../results/ma_05_by_race.csv", row.names = F)


                 
bw_r <- list()
for (r in 1:length(race_groups)) {
bw_r[[r]] <- r_result(outcome = "birthweight_g", quantile = 0.7, os_df=5, race_cat=race_groups[r])
}

bw_df_r_07 <- do.call(rbind, bw_r)

write.csv(bw_df_r_07, file=".../results/bw_07_by_race.csv", row.names = F)


bw_r <- list()
for (r in 1:length(race_groups)) {
  bw_r[[r]] <- r_result(outcome = "birthweight_g", quantile = 0.5, os_df=5, race_cat=race_groups[r])
}

bw_df_r_05 <- do.call(rbind, bw_r)

write.csv(bw_df_r_05, file=".../results/bw_05_by_race.csv", row.names = F)



ga_r <- list()
for (r in 1:length(race_groups)) {
  ga_r[[r]] <- r_result(outcome = "est_gest_wks", quantile = 0.7, os_df=10, race_cat=race_groups[r])
}

ga_df_r_07 <- do.call(rbind, ga_r)

write.csv(ga_df_r_07, file=".../results/ga_07_by_race.csv", row.names = F)



ga_r <- list()
for (r in 1:length(race_groups)) {
  ga_r[[r]] <- r_result(outcome = "est_gest_wks", quantile = 0.5, os_df=10, race_cat=race_groups[r])
}

ga_df_r_05 <- do.call(rbind, ga_r)

write.csv(ga_df_r_05, file=".../results/ga_05_by_race.csv", row.names = F)



bz_r <- list()
for (r in 1:length(race_groups)) {
  bz_r[[r]] <- r_result(outcome = "gzscore", quantile = 0.7, os_df=7, race_cat=race_groups[r])
}

bz_df_r_07 <- do.call(rbind, bz_r)

write.csv(bz_df_r_07, file=".../results/bz_07_by_race.csv", row.names = F)



bz_r <- list()
for (r in 1:length(race_groups)) {
  bz_r[[r]] <- r_result(outcome = "gzscore", quantile = 0.5, os_df=7, race_cat=race_groups[r])
}

bz_df_r_05 <- do.call(rbind, bz_r)

write.csv(bz_df_r_05, file=".../results/bz_05_by_race.csv", row.names = F)




ptb_r <- list()
for (r in 1:length(race_groups)) {
  ptb_r[[r]] <- r_result(outcome = "ptb", quantile = 0.7, os_df=10, race_cat=race_groups[r])
}

ptb_df_r_07 <- do.call(rbind, ptb_r)

write.csv(ptb_df_r_07, file=".../results/ptb_07_by_race.csv", row.names = F)



ptb_r <- list()
for (r in 1:length(race_groups)) {
  ptb_r[[r]] <- r_result(outcome = "ptb", quantile = 0.5, os_df=10, race_cat=race_groups[r])
}

ptb_df_r_05 <- do.call(rbind, ptb_r)

write.csv(ptb_df_r_05, file=".../results/ptb_05_by_race.csv", row.names = F)




sga_r <- list()
for (r in 1:length(race_groups)) {
  sga_r[[r]] <- r_result(outcome = "sga", quantile = 0.7, os_df=7, race_cat=race_groups[r])
}

sga_df_r_07 <- do.call(rbind, sga_r)

write.csv(sga_df_r_07, file=".../results/sga_07_by_race.csv", row.names = F)



sga_r <- list()
for (r in 1:length(race_groups)) {
  sga_r[[r]] <- r_result(outcome = "sga", quantile = 0.5, os_df=7, race_cat=race_groups[r])
}

sga_df_r_05 <- do.call(rbind, sga_r)

write.csv(sga_df_r_05, file=".../results/sga_05_by_race.csv", row.names = F)




lga_r <- list()
for (r in 1:length(race_groups)) {
  lga_r[[r]] <- r_result(outcome = "lga", quantile = 0.7, os_df=7, race_cat=race_groups[r])
}

lga_df_r_07 <- do.call(rbind, lga_r)

write.csv(lga_df_r_07, file=".../results/lga_07_by_race.csv", row.names = F)



lga_r <- list()
for (r in 1:length(race_groups)) {
  lga_r[[r]] <- r_result(outcome = "lga", quantile = 0.5, os_df=7, race_cat=race_groups[r])
}

lga_df_r_05 <- do.call(rbind, lga_r)

write.csv(lga_df_r_05, file=".../results/lga_05_by_race.csv", row.names = F)




# intervention analyses by insurance status

df0$insur <- ifelse(df0$private_ins_pay==1, "private", 
                       ifelse(df0$govt_ins_pay==1, "public", 
                              ifelse(df0$other_pay==1, "other", NA)))

# 6,183 are still unclassified -- recode to other
df0$insur <- ifelse(is.na(df0$insur), "other", df0$insur)

i_result <- function(outcome, quantile, os_df, ins_cat) {
  
  df0$y <- df0[[outcome]]
  
  df_i <- df0 %>% filter(insur==ins_cat)
  
  m2 <- glm(y ~ ns(fluoride_ppm, df=os_df) + mat_age + black_nh + asian_nh + amerind_nh + 
              hipacis_nh + other_nh + multi + hispanic_single + 
              lths + hs_grad + some_col + col_grad +
              Income_Ratio + unemp_rate + avg_temp + Total_pop +   
              factor(urban_code) + factor(con_month) + factor(con_year) + 
              avg.finding_Arsenic,
            data =df_i)
  
  
  
  df1 <- df_i
  df1$fluoride_ppm <- ifelse(df1$fluoride_ppm> quantile, quantile, df1$fluoride_ppm)
  
  p0 <- predict(m2)
  p1 <- predict(m2, newdata = df1)
  
  p_result <- mean(p1-p0)
  
  b_ls <- lapply(1:200, function(x) i_boot(x, outcome = outcome, quantile = quantile, os_df = os_df, ins_cat=ins_cat))
  b_df <- do.call(rbind, b_ls)
  
  result <- data.frame(pop = ins_cat, est = p_result, lb = p_result - 1.96*sqrt(var(b_df)), ub = p_result + 1.96*sqrt(var(b_df)))
  colnames(result) <- c("pop", "est", "lb", "ub")
  
  return(result) 
}



# bootstrap function by insurance status


i_boot <- function(x, outcome, quantile, os_df, ins_cat) {
  print(x)
  
  df_i <- df0 %>% filter(insur==ins_cat)
  water_systems <- unique(df_i$System.ID)
  
  resample_t <- sample(water_systems, length(water_systems), replace=T)
  
  boot_ls <- list()
  for (i in 1:length(resample_t)) {
    boot_ls[[i]] <- df_i %>% filter(System.ID==resample_t[i])
  }
  boot_df <- data.frame(do.call(rbind, boot_ls)) 
  
  boot_df$y <- boot_df[[outcome]]
  
  
  b_fit <- glm(y ~ ns(fluoride_ppm, df=os_df) + mat_age + black_nh + asian_nh + amerind_nh + 
                 hipacis_nh + other_nh + multi + hispanic_single +
                 lths + hs_grad + some_col + col_grad +
                 Income_Ratio + unemp_rate + avg_temp + Total_pop +
                 factor(urban_code) + factor(con_month) + factor(con_year) +
                 avg.finding_Arsenic,
               data =boot_df)
  
  boot_df1 <- boot_df
  boot_df1$fluoride_ppm <- ifelse(boot_df1$fluoride_ppm> quantile, quantile, boot_df1$fluoride_ppm)
  
  p0 <- predict(b_fit)
  p1 <- predict(b_fit, newdata = boot_df1)
  
  p_bw <- p1-p0
  
  b_result <- mean(p_bw)
  return(b_result)
}

set.seed(620573585)


insur_groups <- c("private", "public",  "other")




ma_i <- list()
for (i in 1:length(insur_groups)) {
  ma_i[[i]] <- i_result(outcome = "macrosomia", quantile = 0.7, os_df=5, ins_cat=insur_groups[i])
}

ma_df_i_07 <- do.call(rbind, ma_i)

write.csv(ma_df_i_07, file=".../results/ma_07_by_insur.csv", row.names = F)

ma_i <- list()
for (i in 1:length(insur_groups)) {
  ma_i[[i]] <- i_result(outcome = "macrosomia", quantile = 0.5, os_df=5, ins_cat=insur_groups[i])
}

ma_df_i_05 <- do.call(rbind, ma_i)

write.csv(ma_df_i_05, file=".../results/ma_05_by_insur.csv", row.names = F)


lga_i <- list()
for (i in 1:length(insur_groups)) {
  lga_i[[i]] <- i_result(outcome = "lga", quantile = 0.7, os_df=7, ins_cat=insur_groups[i])
}

lga_df_i_07 <- do.call(rbind, lga_i)

write.csv(lga_df_i_07, file=".../results/lga_07_by_insur.csv", row.names = F)


lga_i <- list()
for (i in 1:length(insur_groups)) {
  lga_i[[i]] <- i_result(outcome = "lga", quantile = 0.5, os_df=7, ins_cat=insur_groups[i])
}

lga_df_i_05 <- do.call(rbind, lga_i)

write.csv(lga_df_i_05, file=".../results/lga_05_by_insur.csv", row.names = F)


sga_i <- list()
for (i in 1:length(insur_groups)) {
  sga_i[[i]] <- i_result(outcome = "sga", quantile = 0.7, os_df=7, ins_cat=insur_groups[i])
}

sga_df_i_07 <- do.call(rbind, sga_i)

write.csv(sga_df_i_07, file=".../results/sga_07_by_insur.csv", row.names = F)


sga_i <- list()
for (i in 1:length(insur_groups)) {
  sga_i[[i]] <- i_result(outcome = "sga", quantile = 0.5, os_df=7, ins_cat=insur_groups[i])
}

sga_df_i_05 <- do.call(rbind, sga_i)

write.csv(sga_df_i_05, file=".../results/sga_05_by_insur.csv", row.names = F)



ptb_i <- list()
for (i in 1:length(insur_groups)) {
  ptb_i[[i]] <- i_result(outcome = "ptb", quantile = 0.7, os_df=10, ins_cat=insur_groups[i])
}

ptb_df_i_07 <- do.call(rbind, ptb_i)

write.csv(ptb_df_i_07, file=".../results/ptb_07_by_insur.csv", row.names = F)




ptb_i <- list()
for (i in 1:length(insur_groups)) {
  ptb_i[[i]] <- i_result(outcome = "ptb", quantile = 0.5, os_df=10, ins_cat=insur_groups[i])
}

ptb_df_i_05 <- do.call(rbind, ptb_i)

write.csv(ptb_df_i_05, file=".../results/ptb_05_by_insur.csv", row.names = F)


# continuous measures

bw_i <- list()
for (i in 1:length(insur_groups)) {
  bw_i[[i]] <- i_result(outcome = "birthweight_g", quantile = 0.7, os_df=5, ins_cat=insur_groups[i])
}

bw_df_i_07 <- do.call(rbind, bw_i)

write.csv(bw_df_i_07, file=".../results/bw_07_by_insur.csv", row.names = F)



bw_i <- list()
for (i in 1:length(insur_groups)) {
  bw_i[[i]] <- i_result(outcome = "birthweight_g", quantile = 0.5, os_df=5, ins_cat=insur_groups[i])
}

bw_df_i_05 <- do.call(rbind, bw_i)

write.csv(bw_df_i_05, file=".../results/bw_05_by_insur.csv", row.names = F)





ga_i <- list()
for (i in 1:length(insur_groups)) {
  ga_i[[i]] <- i_result(outcome = "est_gest_wks", quantile = 0.7, os_df=10, ins_cat=insur_groups[i])
}

ga_df_i_07 <- do.call(rbind, ga_i)

write.csv(ga_df_i_07, file=".../results/ga_07_by_insur.csv", row.names = F)



ga_i <- list()
for (i in 1:length(insur_groups)) {
  ga_i[[i]] <- i_result(outcome = "est_gest_wks", quantile = 0.5, os_df=10, ins_cat=insur_groups[i])
}

ga_df_i_05 <- do.call(rbind, ga_i)

write.csv(ga_df_i_05, file=".../results/ga_05_by_insur.csv", row.names = F)



bz_i <- list()
for (i in 1:length(insur_groups)) {
  bz_i[[i]] <- i_result(outcome = "gzscore", quantile = 0.7, os_df=7, ins_cat=insur_groups[i])
}

bz_df_i_07 <- do.call(rbind, bz_i)

write.csv(bz_df_i_07, file=".../results/bz_07_by_insur.csv", row.names = F)



bz_i <- list()
for (i in 1:length(insur_groups)) {
  bz_i[[i]] <- i_result(outcome = "gzscore", quantile = 0.5, os_df=7, ins_cat=insur_groups[i])
}

bz_df_i_05 <- do.call(rbind, bz_i)

write.csv(bz_df_i_05, file=".../results/bz_05_by_insur.csv", row.names = F)


# intervention analyses by fetal sex


f_result <- function(outcome, quantile, os_df, fs_cat) {
  
  df0$y <- df0[[outcome]]
  
  df_f <- df0 %>% filter(child_sex==fs_cat)
  
  m2 <- glm(y ~ ns(fluoride_ppm, df=os_df) + mat_age + black_nh + asian_nh + amerind_nh + 
              hipacis_nh + other_nh + multi + hispanic_single + 
              lths + hs_grad + some_col + col_grad +
              Income_Ratio + unemp_rate + avg_temp + Total_pop +   
              factor(urban_code) + factor(con_month) + factor(con_year) + 
              avg.finding_Arsenic,
            data =df_f)
  
  
  
  df1 <- df_f
  df1$fluoride_ppm <- ifelse(df1$fluoride_ppm> quantile, quantile, df1$fluoride_ppm)
  
  p0 <- predict(m2)
  p1 <- predict(m2, newdata = df1)
  
  p_result <- mean(p1-p0)
  
  b_ls <- lapply(1:200, function(x) f_boot(x, outcome = outcome, quantile = quantile, os_df = os_df, fs_cat=fs_cat))
  b_df <- do.call(rbind, b_ls)
  
  result <- data.frame(pop = fs_cat, est = p_result, lb = p_result - 1.96*sqrt(var(b_df)), ub = p_result + 1.96*sqrt(var(b_df)))
  colnames(result) <- c("pop", "est", "lb", "ub")
  
  return(result) 
}



# bootstrap function by fetal sex 


f_boot <- function(x, outcome, quantile, os_df, fs_cat) {
  print(x)
  
  df_f <- df0 %>% filter(child_sex==fs_cat)
  water_systems <- unique(df_f$System.ID)
  
  resample_t <- sample(water_systems, length(water_systems), replace=T)
  
  boot_ls <- list()
  for (i in 1:length(resample_t)) {
    boot_ls[[i]] <- df_f %>% filter(System.ID==resample_t[i])
  }
  boot_df <- data.frame(do.call(rbind, boot_ls)) 
  
  boot_df$y <- boot_df[[outcome]]
  
  
  b_fit <- glm(y ~ ns(fluoride_ppm, df=os_df) + mat_age + black_nh + asian_nh + amerind_nh + 
                 hipacis_nh + other_nh + multi + hispanic_single +
                 lths + hs_grad + some_col + col_grad +
                 Income_Ratio + unemp_rate + avg_temp + Total_pop +
                 factor(urban_code) + factor(con_month) + factor(con_year) +
                 avg.finding_Arsenic,
               data =boot_df)
  
  boot_df1 <- boot_df
  boot_df1$fluoride_ppm <- ifelse(boot_df1$fluoride_ppm> quantile, quantile, boot_df1$fluoride_ppm)
  
  p0 <- predict(b_fit)
  p1 <- predict(b_fit, newdata = boot_df1)
  
  p_bw <- p1-p0
  
  b_result <- mean(p_bw)
  return(b_result)
}

set.seed(620573585)


fs_groups <- c("M", "F")



ma_f <- list()
for (i in 1:length(fs_groups)) {
  ma_f[[i]] <- f_result(outcome = "macrosomia", quantile = 0.7, os_df=5, fs_cat=fs_groups[i])
}

ma_df_f_07 <- do.call(rbind, ma_f)

write.csv(ma_df_f_07, file=".../results/ma_07_by_sex.csv", row.names = F)

ma_f <- list()
for (i in 1:length(fs_groups)) {
  ma_f[[i]] <- f_result(outcome = "macrosomia", quantile = 0.5, os_df=5, fs_cat=fs_groups[i])
}

ma_df_f_05 <- do.call(rbind, ma_f)

write.csv(ma_df_f_05, file=".../results/ma_05_by_sex.csv", row.names = F)


lga_f <- list()
for (i in 1:length(fs_groups)) {
  lga_f[[i]] <- f_result(outcome = "lga", quantile = 0.7, os_df=7, fs_cat=fs_groups[i])
}

lga_df_f_07 <- do.call(rbind, lga_f)

write.csv(lga_df_f_07, file=".../results/lga_07_by_sex.csv", row.names = F)


lga_f <- list()
for (i in 1:length(fs_groups)) {
  lga_f[[i]] <- f_result(outcome = "lga", quantile = 0.5, os_df=7, fs_cat=fs_groups[i])
}

lga_df_f_05 <- do.call(rbind, lga_f)

write.csv(lga_df_f_05, file=".../results/lga_05_by_sex.csv", row.names = F)


sga_f <- list()
for (i in 1:length(fs_groups)) {
  sga_f[[i]] <- f_result(outcome = "sga", quantile = 0.7, os_df=7, fs_cat=fs_groups[i])
}

sga_df_f_07 <- do.call(rbind, sga_f)

write.csv(sga_df_f_07, file=".../results/sga_07_by_sex.csv", row.names = F)


sga_f <- list()
for (i in 1:length(fs_groups)) {
  sga_f[[i]] <- f_result(outcome = "sga", quantile = 0.5, os_df=7, fs_cat=fs_groups[i])
}

sga_df_f_05 <- do.call(rbind, sga_f)

write.csv(sga_df_f_05, file=".../results/sga_05_by_sex.csv", row.names = F)



ptb_f <- list()
for (i in 1:length(fs_groups)) {
  ptb_f[[i]] <- f_result(outcome = "ptb", quantile = 0.7, os_df=10, fs_cat=fs_groups[i])
}

ptb_df_f_07 <- do.call(rbind, ptb_f)

write.csv(ptb_df_f_07, file=".../results/ptb_07_by_sex.csv", row.names = F)




ptb_f <- list()
for (i in 1:length(fs_groups)) {
  ptb_f[[i]] <- f_result(outcome = "ptb", quantile = 0.5, os_df=10, fs_cat=fs_groups[i])
}

ptb_df_f_05 <- do.call(rbind, ptb_f)

write.csv(ptb_df_f_05, file=".../results/ptb_05_by_sex.csv", row.names = F)


# continuous measures

bw_f <- list()
for (i in 1:length(fs_groups)) {
  bw_f[[i]] <- f_result(outcome = "birthweight_g", quantile = 0.7, os_df=5, fs_cat=fs_groups[i])
}

bw_df_f_07 <- do.call(rbind, bw_f)

write.csv(bw_df_f_07, file=".../results/bw_07_by_sex.csv", row.names = F)



bw_f <- list()
for (i in 1:length(fs_groups)) {
  bw_f[[i]] <- f_result(outcome = "birthweight_g", quantile = 0.5, os_df=5, fs_cat=fs_groups[i])
}

bw_df_f_05 <- do.call(rbind, bw_f)

write.csv(bw_df_f_05, file=".../results/bw_05_by_sex.csv", row.names = F)



ga_f <- list()
for (i in 1:length(fs_groups)) {
  ga_f[[i]] <- f_result(outcome = "est_gest_wks", quantile = 0.7, os_df=10, fs_cat=fs_groups[i])
}

ga_df_f_07 <- do.call(rbind, ga_f)

write.csv(ga_df_f_07, file=".../results/ga_07_by_sex.csv", row.names = F)



ga_f <- list()
for (i in 1:length(fs_groups)) {
  ga_f[[i]] <- f_result(outcome = "est_gest_wks", quantile = 0.5, os_df=10, fs_cat=fs_groups[i])
}

ga_df_f_05 <- do.call(rbind, ga_f)

write.csv(ga_df_f_05, file=".../results/ga_05_by_sex.csv", row.names = F)



bz_f <- list()
for (i in 1:length(fs_groups)) {
  bz_f[[i]] <- f_result(outcome = "gzscore", quantile = 0.7, os_df=7, fs_cat=fs_groups[i])
}

bz_df_f_07 <- do.call(rbind, bz_f)

write.csv(bz_df_f_07, file=".../results/bz_07_by_sex.csv", row.names = F)



bz_f <- list()
for (i in 1:length(fs_groups)) {
  bz_f[[i]] <- f_result(outcome = "gzscore", quantile = 0.5, os_df=7, fs_cat=fs_groups[i])
}

bz_df_f_05 <- do.call(rbind, bz_f)

write.csv(bz_df_f_05, file=".../results/bz_05_by_sex.csv", row.names = F)



# stratify by arsenic levels 


as_result <- function(outcome, quantile, os_df, as_cat) {
  
  df0$y <- df0[[outcome]]
  
  df_f <- df0 %>% filter(arsenic_mcl==as_cat)
  
  m2 <- glm(y ~ ns(fluoride_ppm, df=os_df) + mat_age + black_nh + asian_nh + amerind_nh + 
              hipacis_nh + other_nh + multi + hispanic_single + 
              lths + hs_grad + some_col + col_grad +
              Income_Ratio + unemp_rate + avg_temp + Total_pop +   
              factor(urban_code) + factor(con_month) + factor(con_year) + 
              avg.finding_Arsenic,
            data =df_f)
  
  
  
  df1 <- df_f
  df1$fluoride_ppm <- ifelse(df1$fluoride_ppm> quantile, quantile, df1$fluoride_ppm)
  
  p0 <- predict(m2)
  p1 <- predict(m2, newdata = df1)
  
  p_result <- mean(p1-p0)
  
  b_ls <- lapply(1:200, function(x) as_boot(x, outcome = outcome, quantile = quantile, os_df = os_df, as_cat=as_cat))
  b_df <- do.call(rbind, b_ls)
  
  result <- data.frame(pop = as_cat, est = p_result, lb = p_result - 1.96*sqrt(var(b_df)), ub = p_result + 1.96*sqrt(var(b_df)))
  colnames(result) <- c("pop", "est", "lb", "ub")
  
  return(result) 
}



# bootstrap function by arsenic levels


as_boot <- function(x, outcome, quantile, os_df, as_cat) {
  print(x)
  
  df_f <- df0 %>% filter(arsenic_mcl==as_cat)
  water_systems <- unique(df_f$System.ID)
  
  resample_t <- sample(water_systems, length(water_systems), replace=T)
  
  boot_ls <- list()
  for (i in 1:length(resample_t)) {
    boot_ls[[i]] <- df_f %>% filter(System.ID==resample_t[i])
  }
  boot_df <- data.frame(do.call(rbind, boot_ls)) 
  
  boot_df$y <- boot_df[[outcome]]
  
  
  b_fit <- glm(y ~ ns(fluoride_ppm, df=os_df) + mat_age + black_nh + asian_nh + amerind_nh + 
                 hipacis_nh + other_nh + multi + hispanic_single +
                 lths + hs_grad + some_col + col_grad +
                 Income_Ratio + unemp_rate + avg_temp + Total_pop +
                 factor(urban_code) + factor(con_month) + factor(con_year) +
                 avg.finding_Arsenic,
               data =boot_df)
  
  boot_df1 <- boot_df
  boot_df1$fluoride_ppm <- ifelse(boot_df1$fluoride_ppm> quantile, quantile, boot_df1$fluoride_ppm)
  
  p0 <- predict(b_fit)
  p1 <- predict(b_fit, newdata = boot_df1)
  
  p_bw <- p1-p0
  
  b_result <- mean(p_bw)
  return(b_result)
}

set.seed(620573585)


as_groups <- c(0,1)


ma_as <- list()
for (i in 1:length(as_groups)) {
  ma_as[[i]] <- as_result(outcome = "macrosomia", quantile = 0.7, os_df=5, as_cat=as_groups[i])
}

ma_df_as_07 <- do.call(rbind, ma_as)

write.csv(ma_df_as_07, file=".../results/ma_07_by_as_mcl.csv", row.names = F)

ma_as <- list()
for (i in 1:length(as_groups)) {
  ma_as[[i]] <- as_result(outcome = "macrosomia", quantile = 0.5, os_df=5, as_cat=as_groups[i])
}

ma_df_as_05 <- do.call(rbind, ma_as)

write.csv(ma_df_as_05, file=".../results/ma_05_by_as_mcl.csv", row.names = F)

lga_as <- list()
for (i in 1:length(as_groups)) {
  lga_as[[i]] <- as_result(outcome = "lga", quantile = 0.7, os_df=7, as_cat=as_groups[i])
}

lga_df_as_07 <- do.call(rbind, lga_as)

write.csv(lga_df_as_07, file=".../results/lga_07_by_as_mcl.csv", row.names = F)


lga_as <- list()
for (i in 1:length(as_groups)) {
  lga_as[[i]] <- as_result(outcome = "lga", quantile = 0.5, os_df=7, as_cat=as_groups[i])
}

lga_df_as_05 <- do.call(rbind, lga_as)

write.csv(lga_df_as_05, file=".../results/lga_05_by_as_mcl.csv", row.names = F)


sga_as <- list()
for (i in 1:length(as_groups)) {
  sga_as[[i]] <- as_result(outcome = "sga", quantile = 0.7, os_df=7, as_cat=as_groups[i])
}

sga_df_as_07 <- do.call(rbind, sga_as)

write.csv(sga_df_as_07, file=".../results/sga_07_by_as_mcl.csv", row.names = F)


sga_as <- list()
for (i in 1:length(as_groups)) {
  sga_as[[i]] <- as_result(outcome = "sga", quantile = 0.5, os_df=7, as_cat=as_groups[i])
}

sga_df_as_05 <- do.call(rbind, sga_as)

write.csv(sga_df_as_05, file=".../results/sga_05_by_as_mcl.csv", row.names = F)



ptb_as <- list()
for (i in 1:length(as_groups)) {
  ptb_as[[i]] <- as_result(outcome = "ptb", quantile = 0.7, os_df=10, as_cat=as_groups[i])
}

ptb_df_as_07 <- do.call(rbind, ptb_as)

write.csv(ptb_df_as_07, file=".../results/ptb_07_by_as_mcl.csv", row.names = F)




ptb_as <- list()
for (i in 1:length(as_groups)) {
  ptb_as[[i]] <- as_result(outcome = "ptb", quantile = 0.5, os_df=10, as_cat=as_groups[i])
}

ptb_df_as_05 <- do.call(rbind, ptb_as)

write.csv(ptb_df_as_05, file=".../results/ptb_05_by_as_mcl.csv", row.names = F)


# continuous measures

bw_as <- list()
for (i in 1:length(as_groups)) {
  bw_as[[i]] <- as_result(outcome = "birthweight_g", quantile = 0.7, os_df=5, as_cat=as_groups[i])
}

bw_df_as_07 <- do.call(rbind, bw_as)

write.csv(bw_df_as_07, file=".../results/bw_07_by_as_mcl.csv", row.names = F)



bw_as <- list()
for (i in 1:length(as_groups)) {
  bw_as[[i]] <- as_result(outcome = "birthweight_g", quantile = 0.5, os_df=5, as_cat=as_groups[i])
}

bw_df_as_05 <- do.call(rbind, bw_as)

write.csv(bw_df_as_05, file=".../results/bw_05_by_as_mcl.csv", row.names = F)


ga_as <- list()
for (i in 1:length(as_groups)) {
  ga_as[[i]] <- as_result(outcome = "est_gest_wks", quantile = 0.7, os_df=10, as_cat=as_groups[i])
}

ga_df_as_07 <- do.call(rbind, ga_as)

write.csv(ga_df_as_07, file=".../results/ga_07_by_as_mcl.csv", row.names = F)



ga_as <- list()
for (i in 1:length(as_groups)) {
  ga_as[[i]] <- as_result(outcome = "est_gest_wks", quantile = 0.5, os_df=10, as_cat=as_groups[i])
}

ga_df_as_05 <- do.call(rbind, ga_as)

write.csv(ga_df_as_05, file=".../results/ga_05_by_as_mcl.csv", row.names = F)



bz_as <- list()
for (i in 1:length(as_groups)) {
  bz_as[[i]] <- as_result(outcome = "gzscore", quantile = 0.7, os_df=7, as_cat=as_groups[i])
}

bz_df_as_07 <- do.call(rbind, bz_as)

write.csv(bz_df_as_07, file=".../results/bz_07_by_as_mcl.csv", row.names = F)



bz_as <- list()
for (i in 1:length(as_groups)) {
  bz_as[[i]] <- as_result(outcome = "gzscore", quantile = 0.5, os_df=7, as_cat=as_groups[i])
}

bz_df_as_05 <- do.call(rbind, bz_as)

write.csv(bz_df_as_05, file=".../results/bz_05_by_as_mcl.csv", row.names = F)

