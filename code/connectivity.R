library(tidyverse)
library(ggpubr)
library(rstatix)
library(datarium)
library(plotly)
library(lme4) # for the analysis
library(haven) # to load the SPSS .sav file
library(RColorBrewer) # needed for some extra colours in one of the graphs
library(lmerTest)# to get p-value estimations that are not part of the standard lme4 packages
library(broom)
library(ggpmisc)
library(ggExtra)
library(psych)
library(lavaan)
library(readxl)
library(semPlot)


setwd("~/Library/CloudStorage/GoogleDrive-yuanygao@stanford.edu/My Drive/Reiss/Projects/Archives/CAMACOSconnectivity/code")
DAP_table = read.csv(file.path('..','data', 'DAP_conn_behav.csv'))

cbPalette <- c("#373f92", "#d54237")

g <- ggplot(data = DAP_table, 
            aes(y = lmem_9_12, x = ldap_tot_c_auc, col = sex, group = sex)) +
  geom_point(show.legend = FALSE, size = 1, alpha = 0.7)+
  geom_smooth(method = "lm",se = TRUE,show.legend = FALSE)+
  theme_bw()+
  theme(legend.position = "none")+
  ylab('Connectivity')+
  xlab('Childhood DAP')+
  scale_colour_manual(values=cbPalette)+
  theme_bw()+
  theme(text=element_text(size=14,  family="Arial"))

g <- ggMarginal(g, groupColour = TRUE,groupFill = TRUE)
g
save_path = file.path('..','figure',"scatter_lmem_9_12_ldap_tot_c_auc.png")
ggsave(save_path, g, width = 3.5, height = 3, dpi = 300)

DAP_table$educcat_mom = factor(DAP_table$educcat_mom)
DAP_table$sex = factor(DAP_table$sex)

DAP_table_boys = DAP_table %>% filter(sex=='Boy')
DAP_table_girls = DAP_table %>% filter(sex=='Girl')

# association between ldap_tot_i_sg2_pg/ldap_tot_c_auc and match_perseverativeErrors/acc/rt/meaning_acc
dataset = DAP_table
# preg
model = lm(formula = 'match_perseverativeErrors ~ 1 + ldap_tot_i_sg2_pg  + sex + age_fnirs + momdl_age2 + educcat_mom + z2home_6 + usyr3 + povcat2_m18y_i + sppstd_i', data = dataset)
summary(model)

model = lm(formula = 'acc ~ 1 + ldap_tot_i_sg2_pg + sex + age_fnirs + momdl_age2 + educcat_mom + z2home_6 + usyr3 + povcat2_m18y_i + sppstd_i', data = dataset)
summary(model)

model = lm(formula = 'rt ~ 1 + ldap_tot_i_sg2_pg + sex + age_fnirs + momdl_age2 + educcat_mom + z2home_6 + usyr3 + povcat2_m18y_i + sppstd_i', data = dataset)
summary(model)

model = lm(formula = 'meaning_size_acc ~ 1 + ldap_tot_i_sg2_pg + sex + age_fnirs + momdl_age2 + educcat_mom + z2home_6 + usyr3 + povcat2_m18y_i + sppstd_i', data = dataset)
summary(model)

# child
model = lm(formula = 'match_perseverativeErrors ~ 1  + ldap_tot_c_auc + sex + age_fnirs + momdl_age2 + educcat_mom + z2home_6 + usyr3 + povcat2_m18y_i + sppstd_i', data = dataset)
summary(model)

model = lm(formula = 'acc ~ 1  + ldap_tot_c_auc + sex + age_fnirs + momdl_age2 + educcat_mom + z2home_6 + usyr3 + povcat2_m18y_i + sppstd_i', data = dataset)
summary(model)

model = lm(formula = 'rt ~ 1  + ldap_tot_c_auc + sex + age_fnirs + momdl_age2 + educcat_mom + z2home_6 + usyr3 + povcat2_m18y_i + sppstd_i', data = dataset)
summary(model)

model = lm(formula = 'meaning_size_acc ~ 1  + ldap_tot_c_auc + sex + age_fnirs + momdl_age2 + educcat_mom + z2home_6 + usyr3 + povcat2_m18y_i + sppstd_i', data = dataset)
summary(model)

# stratified by sex: boys
dataset = DAP_table_boys
# preg
model = lm(formula = 'match_perseverativeErrors ~ 1 + ldap_tot_i_sg2_pg  + age_fnirs + momdl_age2 + educcat_mom + z2home_6 + usyr3 + povcat2_m18y_i + sppstd_i', data = dataset)
summary(model)

model = lm(formula = 'acc ~ 1 + ldap_tot_i_sg2_pg  + age_fnirs + momdl_age2 + educcat_mom + z2home_6 + usyr3 + povcat2_m18y_i + sppstd_i', data = dataset)
summary(model)

model = lm(formula = 'rt ~ 1 + ldap_tot_i_sg2_pg  + age_fnirs + momdl_age2 + educcat_mom + z2home_6 + usyr3 + povcat2_m18y_i + sppstd_i', data = dataset)
summary(model)

model = lm(formula = 'meaning_size_acc ~ 1 + ldap_tot_i_sg2_pg  + age_fnirs + momdl_age2 + educcat_mom + z2home_6 + usyr3 + povcat2_m18y_i + sppstd_i', data = dataset)
summary(model)

# child
model = lm(formula = 'match_perseverativeErrors ~ 1  + ldap_tot_c_auc  + age_fnirs + momdl_age2 + educcat_mom + z2home_6 + usyr3 + povcat2_m18y_i + sppstd_i', data = dataset)
summary(model)

model = lm(formula = 'acc ~ 1  + ldap_tot_c_auc  + age_fnirs + momdl_age2 + educcat_mom + z2home_6 + usyr3 + povcat2_m18y_i + sppstd_i', data = dataset)
summary(model)

model = lm(formula = 'rt ~ 1  + ldap_tot_c_auc  + age_fnirs + momdl_age2 + educcat_mom + z2home_6 + usyr3 + povcat2_m18y_i + sppstd_i', data = dataset)
summary(model)

model = lm(formula = 'meaning_size_acc ~ 1  + ldap_tot_c_auc  + age_fnirs + momdl_age2 + educcat_mom + z2home_6 + usyr3 + povcat2_m18y_i + sppstd_i', data = dataset)
summary(model)



# mediation analysis
DAP_table = read.csv(file.path('processed_data', 'DAP_conn_behav.csv'))
DAP_table_boys = DAP_table %>% filter(sex=='Boy')
mediation_model <- '
  # Direct effects
  lmem_9_12 ~ a * ldap_tot_c_auc + age_fnirs + momdl_age2 + educcat_mom + z2home_6 + usyr3 + povcat2_m18y_i + sppstd_i
  acc ~ c * ldap_tot_c_auc + b * lmem_9_12 + age_fnirs + momdl_age2 + educcat_mom + z2home_6 + usyr3 + povcat2_m18y_i + sppstd_i

  # Indirect effect (a * b)
  indirect := a * b

  # Total effect (c + indirect)
  total := c + indirect
'
mediation_results <- sem(mediation_model, data = DAP_table_boys)
summary(mediation_results, standardized = TRUE, fit.measures = TRUE)

mediation_model <- '
  # Direct effects
  lmem_9_12 ~ a * ldap_tot_c_auc + age_fnirs + momdl_age2 + educcat_mom + z2home_6 + usyr3 + povcat2_m18y_i + sppstd_i
  rt ~ c * ldap_tot_c_auc + b * lmem_9_12 + age_fnirs + momdl_age2 + educcat_mom + z2home_6 + usyr3 + povcat2_m18y_i + sppstd_i

  # Indirect effect (a * b)
  indirect := a * b

  # Total effect (c + indirect)
  total := c + indirect
'
mediation_results <- sem(mediation_model, data = DAP_table_boys)
summary(mediation_results, standardized = TRUE, fit.measures = TRUE)



## boys and girls scores
# Boxplot
ggplot(DAP_table, aes(x = sex, y = acc, color = sex)) +
  geom_boxplot() +
  labs(x = "Sex", y = "WM Accuracy") +
  scale_colour_manual(values=cbPalette)+
  theme_bw()

ggplot(DAP_table, aes(x = sex, y = rt, color = sex)) +
  geom_boxplot() +
  labs(x = "Sex", y = "WM reaction time") +
  scale_colour_manual(values=cbPalette)+
  theme_bw()


ggplot(DAP_table, aes(x = sex, y = match_perseverativeErrors, color = sex)) +
  geom_boxplot() +
  labs(x = "Sex", y = "Perseverative Errors") +
  scale_colour_manual(values=cbPalette)+
  theme_bw()


ggplot(DAP_table, aes(x = sex, y = meaning_size_acc, color = sex)) +
  geom_boxplot() +
  labs(x = "Sex", y = "Accuracy (meaning-size)") +
  scale_colour_manual(values=cbPalette)+
  theme_bw()

# t-test
make_demo <- function(x, y) {
  cat(sprintf("%.2f  %.2f  %.2f  %.2f\n %.2f  %.2f  %.2f  %.2f\n", mean(x, na.rm=TRUE), sd(x, na.rm=TRUE),min(x, na.rm=TRUE),max(x, na.rm=TRUE), mean(y, na.rm=TRUE), sd(y, na.rm=TRUE),min(y, na.rm=TRUE),max(y, na.rm=TRUE)))
  t.test(x, y,
         alternative = c("two.sided", "less", "greater"),
         mu = 0, paired = FALSE, var.equal = FALSE,
         conf.level = 0.95)
}
make_demo(DAP_table_boys$match_perseverativeErrors, DAP_table_girls$match_perseverativeErrors)
make_demo(DAP_table_boys$acc, DAP_table_girls$acc)
make_demo(DAP_table_boys$rt, DAP_table_girls$rt)
make_demo(DAP_table_boys$meaning_size_acc, DAP_table_girls$meaning_size_acc)
