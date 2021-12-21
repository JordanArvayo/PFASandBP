```{r}
##### LOAD PACKAGES
library("nlme")
library("lme4")
library("lmerTest")
library("tidyverse")
library("ggplot2")
library("gtools")
library("sitar")
library("gamlss")
library("haven")
library("foreign")
library("fdapace")
library("dplyr")
library("tidyr")
library("ggpubr")
library("growthstandards")
library("ggpubr")
library("hbgd")
library("data.table")
library("table1")

##### LOAD DATASET
PFAS <-read.csv("C:/Users/jarva/Desktop/James-Todd Lab/PFASandBP/Data/PFAS_L.csv")

```

```{r}
##### COUNT NUMBER OF DUPLICATED IDs
setDT(PFAS)[, count_w := uniqueN(ga_w), by = aid]

##### DESCRIPTIV STATISTICS
hist(PFAS$ga_days, main="Number of BP Measures per Day",
     xlab="Day of Gestation", 
     col="darkturquoise", freq=FALSE)
summary(PFAS$ga_days)
summary(PFAS$count)

hist(PFAS$ga_w, main="Number of BP Measures per Week",
     xlab="Week of Gestation", 
     col="darkturquoise", freq=FALSE)
summary(PFAS$ga_w)
summary(PFAS$count_w)

hist(PFAS$ga_m, main="Number of BP Measures per Month",
     xlab="Month of Gestation", 
     col="darkturquoise", freq=FALSE)
summary(PFAS$ga_m)
summary(PFAS$count_m)

```


```{r}
PFAS$exer <- PFAS$mod_pre_d + PFAS$vig_pre_d
PFAS$race <- PFAS$race2_mom_epi_epia_d
PFAS$edu <- PFAS$coll_grad
PFAS$bmi <- PFAS$bmi_mom_prepreg_d
PFAS$smoke <- PFAS$smokpreg_final_d
PFAS$income <- PFAS$income_hh_epq_epqa_d
PFAS$parity <- PFAS$parity_d

##### RE-LABEL 
PFAS$race <- ifelse(PFAS$race == "white", 0,
                    ifelse(PFAS$race == "black", 1,
                           ifelse(PFAS$race == "hispa", 2,
                                  ifelse(PFAS$race == "asian", 3,
                                         ifelse(PFAS$race == "other", 4,
                                                ifelse(PFAS$race == "more than 1 race", 5, 5))))))
PFAS$edu <- ifelse(PFAS$edu == 0, 1,
                   ifelse(PFAS$edu == 1, 0, 1))
PFAS$bmi <- ifelse(PFAS$bmi <25, 0, 1)
PFAS$smoke <- ifelse(PFAS$smoke == "xnever", 0,
                     ifelse(PFAS$smoke == "former", 1,
                            ifelse(PFAS$smoke == "smoke preg", 2, 2)))
PFAS$income <- ifelse(PFAS$income == 1, 2, 
                      ifelse(PFAS$income == 2, 2,
                             ifelse(PFAS$income == 3, 2,
                                    ifelse(PFAS$income == 4, 2,
                                           ifelse(PFAS$income == 5, 1,
                                                  ifelse(PFAS$income == 6, 0,2))))))
PFAS$parity <- ifelse(PFAS$parity == 0, 0,
                      ifelse(PFAS$parity > 0, 1, 0))

##### REFACTOR
PFAS$race <- factor(PFAS$race,
                    levels = c(0, 1, 2, 3, 4, 5),
                    labels = c("White", "Black", "Hispanic", "Asian", "Other", "More than 1 race"))

PFAS$edu <- factor(PFAS$edu,
                   levels = c(0, 1),
                   labels = c(">= college degree","< college degree"))
PFAS$bmi <- factor(PFAS$bmi,
                   levels = c(0, 1),
                   labels = c("< 25", ">= 25"))
PFAS$smoke <- factor(PFAS$smoke,
                     levels = c(0, 1, 2),
                     labels = c("Never smokes", "Smoked during pregnancy", "Former smoker"))
PFAS$income <- factor(PFAS$income,
                      levels = c(0, 1, 2),
                      labels = c("> $70,000", "$40,001-$70,000", "<= $40,000"))
PFAS$parity <- factor(PFAS$parity,
                      levels = c(0, 1),
                      labels = c("Nulliparous", "Parous"))

table1(~ mom_agey_bp + ga_days + exer + race + bmi + edu + smoke + income + parity, data = PFAS, overall="Total")

```

```{r}
##### CREATE QUARTILES OF PFAS
PFAS$PFOAQ <- quantcut(PFAS$PFOA, q = 4, na.rm = T)
PFAS$PFOSQ <- quantcut(PFAS$PFOS, q = 4, na.rm = T)
PFAS$PFNAQ <- quantcut(PFAS$PFNA2, q = 4, na.rm = T)
PFAS$PFHxSQ <- quantcut(PFAS$PFHxS, q = 4, na.rm = T)
PFAS$ETPFQ <- quantcut(PFAS$Et_PFOSA_AcOH, q = 4, na.rm = T)
PFAS$MEPFQ <- quantcut(PFAS$Me_PFOSA_AcOH2, q = 4, na.rm = T)

##### SUBSET DATA BY PFOA QUARTILE
summary(PFAS$PFOAQ)
PFAS$PFOAH <- ifelse(PFAS$PFOA < mean(PFAS$PFOA), 0, 1 )
#PFAS_PFOA1 <- subset(PFAS, PFAS$PFOA <= 4.3)
#PFAS_PFOA2 <- subset(PFAS, PFAS$PFOA > 4.3 & PFAS$PFOA <= 6.0)
#PFAS_PFOA3 <- subset(PFAS, PFAS$PFOA > 6.0 & PFAS$PFOA <= 8.0)
#PFAS_PFOA4 <- subset(PFAS, PFAS$PFOA > 8.0)

##### SUBSET DATA BY PFOS QUARTILE
summary(PFAS$PFOSQ)
#PFAS_PFOS1 <- subset(PFAS, PFAS$PFOS <= 19.2)
#PFAS_PFOS2 <- subset(PFAS, PFAS$PFOS > 19.2 & PFAS$PFOS <= 26.1)
#PFAS_PFOS3 <- subset(PFAS, PFAS$PFOS > 26.1 & PFAS$PFOS <= 35.6)
#PFAS_PFOS4 <- subset(PFAS, PFAS$PFOS > 35.6)

##### SUBSET DATA BY PFNA2 QUARTILE
summary(PFAS$PFNAQ)
#PFAS_PFNA1 <- subset(PFAS, PFAS$PFNA2 <= 0.5)
#PFAS_PFNA2 <- subset(PFAS, PFAS$PFNA2 > 0.5 & PFAS$PFNA2 <= 0.7)
#PFAS_PFNA3 <- subset(PFAS, PFAS$PFNA2 > 0.7 & PFAS$PFNA2 <= 0.9)
#PFAS_PFNA4 <- subset(PFAS, PFAS$PFNA2 > 0.9)

##### SUBSET DATA BY PFHxS QUARTILE
summary(PFAS$PFHxSQ)
#PFAS_PFHxS1 <- subset(PFAS, PFAS$PFHxS <= 1.7)
#PFAS_PFHxS2 <- subset(PFAS, PFAS$PFHxS > 1.7 & PFAS$PFHxS <= 2.5)
#PFAS_PFHxS3 <- subset(PFAS, PFAS$PFHxS > 2.5 & PFAS$PFHxS <= 3.8)
#PFAS_PFHxS4 <- subset(PFAS, PFAS$PFHxS > 3.8)

##### SUBSET DATA BY Et_PFOSA_AcOH QUARTILE
summary(PFAS$ETPFQ)
#PFAS_ETPF1 <- subset(PFAS, PFAS$Et_PFOSA_AcOH <= 0.8)
#PFAS_ETPF2 <- subset(PFAS, PFAS$Et_PFOSA_AcOH > 0.8 & PFAS$Et_PFOSA_AcOH <= 1.2)
#PFAS_ETPF3 <- subset(PFAS, PFAS$Et_PFOSA_AcOH > 1.2 & PFAS$Et_PFOSA_AcOH <= 1.9)
#PFAS_ETPF4 <- subset(PFAS, PFAS$Et_PFOSA_AcOH > 1.9)

##### SUBSET DATA BY Me_PFOSA_AcOH2 QUARTILE
summary(PFAS$MEPFQ)
#PFAS_MEPF1 <- subset(PFAS, PFAS$Me_PFOSA_AcOH2 <= 1.3)
#PFAS_MEPF2 <- subset(PFAS, PFAS$Me_PFOSA_AcOH2 > 1.3 & PFAS$Me_PFOSA_AcOH2 <= 2)
#PFAS_MEPF3 <- subset(PFAS, PFAS$Me_PFOSA_AcOH2 > 2 & PFAS$Me_PFOSA_AcOH2 <= 3.2)
#PFAS_MEPF4 <- subset(PFAS, PFAS$Me_PFOSA_AcOH2 > 3.2)

```


```{r}
##### IDENTIFY OUTLIERS FOR SYSTOLIC BP
outliers_sys <- velout(x=ga_w_sc, y=ga_w_sys, id=aid, data=PFAS, limit=3)

## SET 205 SYS OUTLIERS MISSING
PFAS_SYS <- zapvelout(outliers_sys, icode=c(4,6))
PFAS_SYS <- subset(PFAS_SYS, !is.na(PFAS_SYS$ga_w_sys))


##### IDENTIFY OUTLIERS FOR DIASTOLIC BP
outliers_dias <- velout(x=ga_w_sc, y=ga_w_dias, id=aid, data=PFAS, limit=3)

## SET 246 DIAS OUTLIERS MISSING
PFAS_DIAS <- zapvelout(outliers_dias, icode=c(4,6))
PFAS_DIAS <- subset(PFAS_DIAS, !is.na(PFAS_DIAS$ga_w_dias))

```

```{r}
##### CHECK DEGREES OF FREEDOM
dfset(ga_w_sc, ga_w_sys, PFAS_SYS, FUN=BIC) # df = 2
dfset(ga_w_sc, ga_w_dias, PFAS_DIAS, FUN=BIC) # df = 2
```

```{r}
##### FITS SITAR 

# SYSTOLIC BP
m11_w <- sitar(x=ga_w_sc, y=ga_w_sys, id=aid, df=2, data=PFAS_SYS, fixed='a', random='a+b+c') #Total pop.


# DIASTOLIC BP
m21_w <- sitar(x=ga_w_sc, y=ga_w_dias, id=aid, df=2, data=PFAS_DIAS, fixed='a', random='a+b+c') #Total pop.

# SUMMARIES
#print(m11_w)
#summary(m11_w)

#print(m21_w)
#summary(m21_w)


```

```{r}
# RANDOM EFFECTS FOR SYSTOLIC AND DIASTOLIC BP
#ranef(m11_w) # SITAR random effects for systolic bp
#ranef(m21_w) # SITAR random effects for diastolic bp

# Attaching the random effects a and c to the systolic data
random_coef <- as.data.frame(m11_w$coefficients$random$id)
random_coef$aid <- as.numeric(rownames(random_coef))
PFAS_SYS <- merge(PFAS_SYS, random_coef, by = "aid")
fit_sys<-m11_w$fitted

# Renaming the columns in m11_w$fitted from "fixed" and "id" to "Systolic fixed effect" and "Systolic predicted"
colnames(fit_sys)<-c("Sys_Fixed_Effect","Sys_Predicted")
PFAS_SYS_R <- cbind(PFAS_SYS,fit_sys)
PFAS_SYS_R$a_sys <- PFAS_SYS_R$a
PFAS_SYS_R$b_sys <- PFAS_SYS_R$b
PFAS_SYS_R$c_sys <- PFAS_SYS_R$c


# Attaching the random effects a and c to the diastolic data
random_coef2 <- as.data.frame(m21_w$coefficients$random$id)
random_coef2$aid <- as.numeric(rownames(random_coef2))
PFAS_DIAS <- merge(PFAS_DIAS, random_coef2, by = "aid")
fit_dias<-m21_w$fitted

# Renaming the columns in m1$fitted from "fixed" and "id" to "Male fixed effect" and "Male predicted height"
colnames(fit_dias)<-c("Dias_Fixed_Effect","Dias_Predicted")
PFAS_DIAS_R <- cbind(PFAS_DIAS,fit_dias)
PFAS_DIAS_R$a_dias <- PFAS_DIAS_R$a
PFAS_DIAS_R$b_dias <- PFAS_DIAS_R$b
PFAS_DIAS_R$c_dias <- PFAS_DIAS_R$c
```


```{r}
##### EFFECT OF PFAS ON SITAR SIZE PARAMTER ON SYSTOLIC BP

## RACE STRATA ##
pfoa_a_sys <- lm(a_sys ~ PFOAQ + mom_agey_bp + race + edu + bmi + smoke + income + parity, data = PFAS_SYS_R)
summary(pfoa_a_sys)
pfoa_a_sys_w <- lm(a_sys ~ PFOAQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_SYS_R, race == "White"))
summary(pfoa_a_sys_w)
pfoa_a_sys_b <- lm(a_sys ~ PFOAQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_SYS_R, race == "Black"))
summary(pfoa_a_sys_b)
pfoa_a_sys_h <- lm(a_sys ~ PFOAQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_SYS_R, race == "Hispanic"))
summary(pfoa_a_sys_h)
pfoa_a_sys_a <- lm(a_sys ~ PFOAQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_SYS_R, race == "Asian"))
summary(pfoa_a_sys_a)
pfoa_a_sys_om <- lm(a_sys ~ PFOAQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_SYS_R, race == "Other" | race == "More than 1 race"))
summary(pfoa_a_sys_om)


pfos_a_sys <- lm(a_sys ~ PFOSQ + mom_agey_bp + race + edu + bmi + smoke + income + parity, data = PFAS_SYS_R)
summary(pfos_a_sys)
pfos_a_sys_w <- lm(a_sys ~ PFOSQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_SYS_R, race == "White"))
summary(pfos_a_sys_w)
pfos_a_sys_b <- lm(a_sys ~ PFOSQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_SYS_R, race == "Black"))
summary(pfos_a_sys_b)
pfos_a_sys_h <- lm(a_sys ~ PFOSQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_SYS_R, race == "Hispanic"))
summary(pfos_a_sys_h)
pfos_a_sys_a <- lm(a_sys ~ PFOSQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_SYS_R, race == "Asian"))
summary(pfos_a_sys_a)
pfos_a_sys_om <- lm(a_sys ~ PFOSQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_SYS_R, race == "Other" | race == "More than 1 race"))
summary(pfos_a_sys_om)


pfna2_a_sys <- lm(a_sys ~ PFNAQ + mom_agey_bp + race + edu + bmi + smoke + income + parity, data = PFAS_SYS_R)
summary(pfna2_a_sys)
pfna2_a_sys_w <- lm(a_sys ~ PFNAQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_SYS_R, race == "White"))
summary(pfna2_a_sys_w)
pfna2_a_sys_b <- lm(a_sys ~ PFNAQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_SYS_R, race == "Black"))
summary(pfna2_a_sys_b)
pfna2_a_sys_h <- lm(a_sys ~ PFNAQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_SYS_R, race == "Hispanic"))
summary(pfna2_a_sys_h)
pfna2_a_sys_a <- lm(a_sys ~ PFNAQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_SYS_R, race == "Asian"))
summary(pfna2_a_sys_a)
pfna2_a_sys_om <- lm(a_sys ~ PFNAQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_SYS_R, race == "Other" | race == "More than 1 race"))
summary(pfna2_a_sys_om)


pfhxs_a_sys <- lm(a_sys ~ PFHxSQ + mom_agey_bp + race + edu + bmi + smoke + income + parity, data = PFAS_SYS_R)
summary(pfhxs_a_sys)
pfhxs_a_sys_w <- lm(a_sys ~ PFHxSQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_SYS_R, race == "White"))
summary(pfhxs_a_sys_w)
pfhxs_a_sys_b <- lm(a_sys ~ PFHxSQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_SYS_R, race == "Black"))
summary(pfhxs_a_sys_b)
pfhxs_a_sys_h <- lm(a_sys ~ PFHxSQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_SYS_R, race == "Hispanic"))
summary(pfhxs_a_sys_h)
pfhxs_a_sys_a <- lm(a_sys ~ PFHxSQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_SYS_R, race == "Asian"))
summary(pfhxs_a_sys_a)
pfhxs_a_sys_om <- lm(a_sys ~ PFHxSQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_SYS_R, race == "Other" | race == "More than 1 race"))
summary(pfhxs_a_sys_om)


etpfq_a_sys <- lm(a_sys ~ ETPFQ + mom_agey_bp + race + edu + bmi + smoke + income + parity, data = PFAS_SYS_R)
summary(etpfq_a_sys)
etpfq_a_sys_w <- lm(a_sys ~ ETPFQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_SYS_R, race == "White"))
summary(etpfq_a_sys_w)
etpfq_a_sys_b <- lm(a_sys ~ ETPFQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_SYS_R, race == "Black"))
summary(etpfq_a_sys_b)
etpfq_a_sys_h <- lm(a_sys ~ ETPFQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_SYS_R, race == "Hispanic"))
summary(etpfq_a_sys_h)
etpfq_a_sys_a <- lm(a_sys ~ ETPFQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_SYS_R, race == "Asian"))
summary(etpfq_a_sys_a)
etpfq_a_sys_om <- lm(a_sys ~ ETPFQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_SYS_R, race == "Other" | race == "More than 1 race"))
summary(etpfq_a_sys_om)


mepfq_a_sys <- lm(a_sys ~ MEPFQ + mom_agey_bp + race + edu + bmi + smoke + income + parity, data = PFAS_SYS_R)
summary(mepfq_a_sys)
mepfq_a_sys_w <- lm(a_sys ~ MEPFQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_SYS_R, race == "White"))
summary(mepfq_a_sys_w)
mepfq_a_sys_b <- lm(a_sys ~ MEPFQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_SYS_R, race == "Black"))
summary(mepfq_a_sys_b)
mepfq_a_sys_h <- lm(a_sys ~ MEPFQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_SYS_R, race == "Hispanic"))
summary(mepfq_a_sys_h)
mepfq_a_sys_a <- lm(a_sys ~ MEPFQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_SYS_R, race == "Asian"))
summary(mepfq_a_sys_a)
mepfq_a_sys_om <- lm(a_sys ~ MEPFQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_SYS_R, race == "Other" | race == "More than 1 race"))
summary(mepfq_a_sys_om)


## AGE STRATA ##

pfoa_a_sys <- lm(a_sys ~ PFOAQ + mom_agey_bp + race + edu + bmi + smoke + income + parity, data = PFAS_SYS_R)
summary(pfoa_a_sys)
pfoa_a_sys_y <- lm(a_sys ~ PFOAQ + race + edu + bmi + smoke + income + parity, data = subset(PFAS_SYS_R, mom_agey_bp < 35.0))
summary(pfoa_a_sys_y)
pfoa_a_sys_o <- lm(a_sys ~ PFOAQ + race + edu + bmi + smoke + income + parity, data = subset(PFAS_SYS_R, mom_agey_bp >= 35.0))
summary(pfoa_a_sys_o)


pfos_a_sys <- lm(a_sys ~ PFOSQ + mom_agey_bp + race + edu + bmi + smoke + income + parity, data = PFAS_SYS_R)
summary(pfos_a_sys)
pfos_a_sys_y <- lm(a_sys ~ PFOSQ + race + edu + bmi + smoke + income + parity, data = subset(PFAS_SYS_R, mom_agey_bp < 35.0))
summary(pfos_a_sys_y)
pfos_a_sys_o <- lm(a_sys ~ PFOSQ + race + edu + bmi + smoke + income + parity, data = subset(PFAS_SYS_R, mom_agey_bp >= 35.0))
summary(pfos_a_sys_o)


pfna2_a_sys <- lm(a_sys ~ PFNAQ + mom_agey_bp + race + edu + bmi + smoke + income + parity, data = PFAS_SYS_R)
summary(pfna2_a_sys)
pfna2_a_sys_y <- lm(a_sys ~ PFNAQ + race + edu + bmi + smoke + income + parity, data = subset(PFAS_SYS_R, mom_agey_bp < 35.0))
summary(pfna2_a_sys_y)
pfna2_a_sys_o <- lm(a_sys ~ PFNAQ + race + edu + bmi + smoke + income + parity, data = subset(PFAS_SYS_R, mom_agey_bp >= 35.0))
summary(pfna2_a_sys_o)


pfhxs_a_sys <- lm(a_sys ~ PFHxSQ + mom_agey_bp + race + edu + bmi + smoke + income + parity, data = PFAS_SYS_R)
summary(pfhxs_a_sys)
pfhxs_a_sys_y <- lm(a_sys ~ PFHxSQ + race + edu + bmi + smoke + income + parity, data = subset(PFAS_SYS_R, mom_agey_bp < 35.0))
summary(pfhxs_a_sys_y)
pfhxs_a_sys_o <- lm(a_sys ~ PFHxSQ + race + edu + bmi + smoke + income + parity, data = subset(PFAS_SYS_R, mom_agey_bp >= 35.0))
summary(pfhxs_a_sys_o)


etpfq_a_sys <- lm(a_sys ~ ETPFQ + mom_agey_bp + race + edu + bmi + smoke + income + parity, data = PFAS_SYS_R)
summary(etpfq_a_sys)
etpfq_a_sys_y <- lm(a_sys ~ ETPFQ + race + edu + bmi + smoke + income + parity, data = subset(PFAS_SYS_R, mom_agey_bp < 35.0))
summary(etpfq_a_sys_y)
etpfq_a_sys_o <- lm(a_sys ~ ETPFQ + race + edu + bmi + smoke + income + parity, data = subset(PFAS_SYS_R, mom_agey_bp >= 35.0))
summary(etpfq_a_sys_o)


mepfq_a_sys <- lm(a_sys ~ MEPFQ + mom_agey_bp + race + edu + bmi + smoke + income + parity, data = PFAS_SYS_R)
summary(mepfq_a_sys)
mepfq_a_sys_y <- lm(a_sys ~ MEPFQ + race + edu + bmi + smoke + income + parity, data = subset(PFAS_SYS_R, mom_agey_bp < 35.0))
summary(mepfq_a_sys_y)
mepfq_a_sys_o <- lm(a_sys ~ MEPFQ + race + edu + bmi + smoke + income + parity, data = subset(PFAS_SYS_R, mom_agey_bp >= 35.0))
summary(mepfq_a_sys_o)


## BMI STRATA ##

pfoa_a_sys <- lm(a_sys ~ PFOAQ + mom_agey_bp + race + edu + bmi + smoke + income + parity, data = PFAS_SYS_R)
summary(pfoa_a_sys)
pfoa_a_sys_l <- lm(a_sys ~ PFOAQ + mom_agey_bp + race + edu + smoke + income + parity, data = subset(PFAS_SYS_R, bmi_mom_prepreg_d < 25.0))
summary(pfoa_a_sys_l)
pfoa_a_sys_h <- lm(a_sys ~ PFOAQ + mom_agey_bp + race + edu + smoke + income + parity, data = subset(PFAS_SYS_R, bmi_mom_prepreg_d >= 25.0))
summary(pfoa_a_sys_h)


pfos_a_sys <- lm(a_sys ~ PFOSQ + mom_agey_bp + race + edu + bmi + smoke + income + parity, data = PFAS_SYS_R)
summary(pfos_a_sys)
pfos_a_sys_l <- lm(a_sys ~ PFOSQ + mom_agey_bp + race + edu + smoke + income + parity, data = subset(PFAS_SYS_R, bmi_mom_prepreg_d < 25.0))
summary(pfos_a_sys_l)
pfos_a_sys_h <- lm(a_sys ~ PFOSQ + mom_agey_bp + race + edu + smoke + income + parity, data = subset(PFAS_SYS_R, bmi_mom_prepreg_d >= 25.0))
summary(pfos_a_sys_h)


pfna2_a_sys <- lm(a_sys ~ PFNAQ + mom_agey_bp + race + edu + bmi + smoke + income + parity, data = PFAS_SYS_R)
summary(pfna2_a_sys)
pfna2_a_sys_l <- lm(a_sys ~ PFNAQ + mom_agey_bp + race + edu + smoke + income + parity, data = subset(PFAS_SYS_R, bmi_mom_prepreg_d < 25.0))
summary(pfna2_a_sys_l)
pfna2_a_sys_h <- lm(a_sys ~ PFNAQ + mom_agey_bp + race + edu + smoke + income + parity, data = subset(PFAS_SYS_R, bmi_mom_prepreg_d >= 25.0))
summary(pfna2_a_sys_h)


pfhxs_a_sys <- lm(a_sys ~ PFHxSQ + mom_agey_bp + race + edu + bmi + smoke + income + parity, data = PFAS_SYS_R)
summary(pfhxs_a_sys)
pfhxs_a_sys_l <- lm(a_sys ~ PFHxSQ + mom_agey_bp + race + edu + smoke + income + parity, data = subset(PFAS_SYS_R, bmi_mom_prepreg_d < 25.0))
summary(pfhxs_a_sys_l)
pfhxs_a_sys_h <- lm(a_sys ~ PFHxSQ + mom_agey_bp + race + edu + smoke + income + parity, data = subset(PFAS_SYS_R, bmi_mom_prepreg_d >= 25.0))
summary(pfhxs_a_sys_h)


etpfq_a_sys <- lm(a_sys ~ ETPFQ + mom_agey_bp + race + edu + bmi + smoke + income + parity, data = PFAS_SYS_R)
summary(etpfq_a_sys)
etpfq_a_sys_l <- lm(a_sys ~ ETPFQ + mom_agey_bp + race + edu + smoke + income + parity, data = subset(PFAS_SYS_R, bmi_mom_prepreg_d < 25.0))
summary(etpfq_a_sys_l)
etpfq_a_sys_h <- lm(a_sys ~ ETPFQ + mom_agey_bp + race + edu + smoke + income + parity, data = subset(PFAS_SYS_R, bmi_mom_prepreg_d >= 25.0))
summary(etpfq_a_sys_h)


mepfq_a_sys <- lm(a_sys ~ MEPFQ + mom_agey_bp + race + edu + bmi + smoke + income + parity, data = PFAS_SYS_R)
summary(mepfq_a_sys)
mepfq_a_sys_l <- lm(a_sys ~ MEPFQ + mom_agey_bp + race + edu + smoke + income + parity, data = subset(PFAS_SYS_R, bmi_mom_prepreg_d < 25.0))
summary(mepfq_a_sys_l)
mepfq_a_sys_h <- lm(a_sys ~ MEPFQ + mom_agey_bp + race + edu + smoke + income + parity, data = subset(PFAS_SYS_R, bmi_mom_prepreg_d >= 25.0))
summary(mepfq_a_sys_h)

```

```{r}
## EFFECT OF PFAS ON SITAR SIZE PARAMTER ON DIASTOLIC BP

## RACE STRATA ##
pfoa_a_dias <- lm(a_dias ~ PFOAQ + mom_agey_bp + race + edu + bmi + smoke + income + parity, data = PFAS_DIAS_R)
summary(pfoa_a_dias)
pfoa_a_dias_w <- lm(a_dias ~ PFOAQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_DIAS_R, race == "White"))
summary(pfoa_a_dias_w)
pfoa_a_diass_b <- lm(a_dias ~ PFOAQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_DIAS_R, race == "Black"))
summary(pfoa_a_diass_b)
pfoa_a_diass_h <- lm(a_dias ~ PFOAQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_DIAS_R, race == "Hispanic"))
summary(pfoa_a_diass_h)
pfoa_a_diass_a <- lm(a_dias ~ PFOAQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_DIAS_R, race == "Asian"))
summary(pfoa_a_diass_a)
pfoa_a_diass_om <- lm(a_dias ~ PFOAQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_DIAS_R, race == "Other" | race == "More than 1 race"))
summary(pfoa_a_diass_om)


pfos_a_dias <- lm(a_dias ~ PFOSQ + mom_agey_bp + race + edu + bmi + smoke + income + parity, data = PFAS_DIAS_R)
summary(pfos_a_dias)
pfos_a_dias_w <- lm(a_dias ~ PFOSQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_DIAS_R, race == "White"))
summary(pfos_a_dias_w)
pfos_a_diass_b <- lm(a_dias ~ PFOSQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_DIAS_R, race == "Black"))
summary(pfos_a_diass_b)
pfos_a_diass_h <- lm(a_dias ~ PFOSQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_DIAS_R, race == "Hispanic"))
summary(pfos_a_diass_h)
pfos_a_diass_a <- lm(a_dias ~ PFOSQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_DIAS_R, race == "Asian"))
summary(pfos_a_diass_a)
pfos_a_diass_om <- lm(a_dias ~ PFOSQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_DIAS_R, race == "Other" | race == "More than 1 race"))
summary(pfos_a_diass_om)


pfna2_a_dias <- lm(a_dias ~ PFNAQ + mom_agey_bp + race + edu + bmi + smoke + income + parity, data = PFAS_DIAS_R)
summary(pfna2_a_dias)
pfna2_a_dias_w <- lm(a_dias ~ PFNAQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_DIAS_R, race == "White"))
summary(pfna2_a_dias_w)
pfna2_a_diass_b <- lm(a_dias ~ PFNAQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_DIAS_R, race == "Black"))
summary(pfna2_a_diass_b)
pfna2_a_diass_h <- lm(a_dias ~ PFNAQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_DIAS_R, race == "Hispanic"))
summary(pfna2_a_diass_h)
pfna2_a_diass_a <- lm(a_dias ~ PFNAQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_DIAS_R, race == "Asian"))
summary(pfna2_a_diass_a)
pfna2_a_diass_om <- lm(a_dias ~ PFNAQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_DIAS_R, race == "Other" | race == "More than 1 race"))
summary(pfna2_a_diass_om)


pfhxs_a_dias <- lm(a_dias ~ PFHxSQ + mom_agey_bp + race + edu + bmi + smoke + income + parity, data = PFAS_DIAS_R)
summary(pfhxs_a_dias)
pfhxs_a_dias_w <- lm(a_dias ~ PFHxSQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_DIAS_R, race == "White"))
summary(pfhxs_a_dias_w)
pfhxs_a_diass_b <- lm(a_dias ~ PFHxSQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_DIAS_R, race == "Black"))
summary(pfhxs_a_diass_b)
pfhxs_a_diass_h <- lm(a_dias ~ PFHxSQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_DIAS_R, race == "Hispanic"))
summary(pfhxs_a_diass_h)
pfhxs_a_diass_a <- lm(a_dias ~ PFHxSQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_DIAS_R, race == "Asian"))
summary(pfhxs_a_diass_a)
pfhxs_a_diass_om <- lm(a_dias ~ PFHxSQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_DIAS_R, race == "Other" | race == "More than 1 race"))
summary(pfhxs_a_diass_om)


etpfq_a_dias <- lm(a_dias ~ ETPFQ + mom_agey_bp + race + edu + bmi + smoke + income + parity, data = PFAS_DIAS_R)
summary(etpfq_a_dias)
etpfq_a_dias_w <- lm(a_dias ~ ETPFQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_DIAS_R, race == "White"))
summary(etpfq_a_dias_w)
etpfq_a_diass_b <- lm(a_dias ~ ETPFQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_DIAS_R, race == "Black"))
summary(etpfq_a_diass_b)
etpfq_a_diass_h <- lm(a_dias ~ ETPFQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_DIAS_R, race == "Hispanic"))
summary(etpfq_a_diass_h)
etpfq_a_diass_a <- lm(a_dias ~ ETPFQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_DIAS_R, race == "Asian"))
summary(etpfq_a_diass_a)
etpfq_a_diass_om <- lm(a_dias ~ ETPFQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_DIAS_R, race == "Other" | race == "More than 1 race"))
summary(etpfq_a_diass_om)


mepfq_a_dias <- lm(a_dias ~ MEPFQ + mom_agey_bp + race + edu + bmi + smoke + income + parity, data = PFAS_DIAS_R)
summary(mepfq_a_dias)
mepfq_a_dias_w <- lm(a_dias ~ MEPFQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_DIAS_R, race == "White"))
summary(mepfq_a_dias_w)
mepfq_a_diass_b <- lm(a_dias ~ MEPFQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_DIAS_R, race == "Black"))
summary(mepfq_a_diass_b)
mepfq_a_diass_h <- lm(a_dias ~ MEPFQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_DIAS_R, race == "Hispanic"))
summary(mepfq_a_diass_h)
mepfq_a_diass_a <- lm(a_dias ~ MEPFQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_DIAS_R, race == "Asian"))
summary(mepfq_a_diass_a)
mepfq_a_diass_om <- lm(a_dias ~ MEPFQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_DIAS_R, race == "Other" | race == "More than 1 race"))
summary(mepfq_a_diass_om)


## AGE STRATA ##

pfoa_a_dias <- lm(a_dias ~ PFOAQ + mom_agey_bp + race + edu + bmi + smoke + income + parity, data = PFAS_DIAS_R)
summary(pfoa_a_dias)
pfoa_a_dias_y <- lm(a_dias ~ PFOAQ + race + edu + bmi + smoke + income + parity, data = subset(PFAS_DIAS_R, mom_agey_bp < 35.0))
summary(pfoa_a_dias_y)
pfoa_a_dias_o <- lm(a_dias ~ PFOAQ + race + edu + bmi + smoke + income + parity, data = subset(PFAS_DIAS_R, mom_agey_bp >= 35.0))
summary(pfoa_a_dias_o)


pfos_a_dias <- lm(a_dias ~ PFOSQ + mom_agey_bp + race + edu + bmi + smoke + income + parity, data = PFAS_DIAS_R)
summary(pfos_a_dias)
pfos_a_dias_y <- lm(a_dias ~ PFOSQ + race + edu + bmi + smoke + income + parity, data = subset(PFAS_DIAS_R, mom_agey_bp < 35.0))
summary(pfos_a_dias_y)
pfos_a_dias_o <- lm(a_dias ~ PFOSQ + race + edu + bmi + smoke + income + parity, data = subset(PFAS_DIAS_R, mom_agey_bp >= 35.0))
summary(pfos_a_dias_o)


pfna2_a_dias <- lm(a_dias ~ PFNAQ + mom_agey_bp + race + edu + bmi + smoke + income + parity, data = PFAS_DIAS_R)
summary(pfna2_a_dias)
pfna2_a_dias_y <- lm(a_dias ~ PFNAQ + race + edu + bmi + smoke + income + parity, data = subset(PFAS_DIAS_R, mom_agey_bp < 35.0))
summary(pfna2_a_dias_y)
pfna2_a_dias_o <- lm(a_dias ~ PFNAQ + race + edu + bmi + smoke + income + parity, data = subset(PFAS_DIAS_R, mom_agey_bp >= 35.0))
summary(pfna2_a_dias_o)


pfhxs_a_dias <- lm(a_dias ~ PFHxSQ + mom_agey_bp + race + edu + bmi + smoke + income + parity, data = PFAS_DIAS_R)
summary(pfhxs_a_dias)
pfhxs_a_dias_y <- lm(a_dias ~ PFHxSQ + race + edu + bmi + smoke + income + parity, data = subset(PFAS_DIAS_R, mom_agey_bp < 35.0))
summary(pfhxs_a_dias_y)
pfhxs_a_dias_o <- lm(a_dias ~ PFHxSQ + race + edu + bmi + smoke + income + parity, data = subset(PFAS_DIAS_R, mom_agey_bp >= 35.0))
summary(pfhxs_a_dias_o)


etpfq_a_dias <- lm(a_dias ~ ETPFQ + mom_agey_bp + race + edu + bmi + smoke + income + parity, data = PFAS_DIAS_R)
summary(etpfq_a_dias)
etpfq_a_dias_y <- lm(a_dias ~ ETPFQ + race + edu + bmi + smoke + income + parity, data = subset(PFAS_DIAS_R, mom_agey_bp < 35.0))
summary(etpfq_a_dias_y)
etpfq_a_dias_o <- lm(a_dias ~ ETPFQ + race + edu + bmi + smoke + income + parity, data = subset(PFAS_DIAS_R, mom_agey_bp >= 35.0))
summary(etpfq_a_dias_o)


mepfq_a_dias <- lm(a_dias ~ MEPFQ + mom_agey_bp + race + edu + bmi + smoke + income + parity, data = PFAS_DIAS_R)
summary(mepfq_a_dias)
mepfq_a_dias_y <- lm(a_dias ~ MEPFQ + race + edu + bmi + smoke + income + parity, data = subset(PFAS_DIAS_R, mom_agey_bp < 35.0))
summary(mepfq_a_dias_y)
mepfq_a_dias_o <- lm(a_dias ~ MEPFQ + race + edu + bmi + smoke + income + parity, data = subset(PFAS_DIAS_R, mom_agey_bp >= 35.0))
summary(mepfq_a_dias_o)


## BMI STRATA ##

pfoa_a_dias <- lm(a_dias ~ PFOAQ + mom_agey_bp + race + edu + bmi + smoke + income + parity, data = PFAS_DIAS_R)
summary(pfoa_a_dias)
pfoa_a_dias_l <- lm(a_dias ~ PFOAQ + mom_agey_bp + race + edu + smoke + income + parity, data = subset(PFAS_DIAS_R, bmi_mom_prepreg_d < 25.0))
summary(pfoa_a_dias_l)
pfoa_a_dias_h <- lm(a_dias ~ PFOAQ + mom_agey_bp + race + edu + smoke + income + parity, data = subset(PFAS_DIAS_R, bmi_mom_prepreg_d >= 25.0))
summary(pfoa_a_dias_h)


pfos_a_dias <- lm(a_dias ~ PFOSQ + mom_agey_bp + race + edu + bmi + smoke + income + parity, data = PFAS_DIAS_R)
summary(pfos_a_dias)
pfos_a_dias_l <- lm(a_dias ~ PFOSQ + mom_agey_bp + race + edu + smoke + income + parity, data = subset(PFAS_DIAS_R, bmi_mom_prepreg_d < 25.0))
summary(pfos_a_dias_l)
pfos_a_dias_h <- lm(a_dias ~ PFOSQ + mom_agey_bp + race + edu + smoke + income + parity, data = subset(PFAS_DIAS_R, bmi_mom_prepreg_d >= 25.0))
summary(pfos_a_dias_h)


pfna2_a_dias <- lm(a_dias ~ PFNAQ + mom_agey_bp + race + edu + bmi + smoke + income + parity, data = PFAS_DIAS_R)
summary(pfna2_a_dias)
pfna2_a_dias_l <- lm(a_dias ~ PFNAQ + mom_agey_bp + race + edu + smoke + income + parity, data = subset(PFAS_DIAS_R, bmi_mom_prepreg_d < 25.0))
summary(pfna2_a_dias_l)
pfna2_a_dias_h <- lm(a_dias ~ PFNAQ + mom_agey_bp + race + edu + smoke + income + parity, data = subset(PFAS_DIAS_R, bmi_mom_prepreg_d >= 25.0))
summary(pfna2_a_dias_h)


pfhxs_a_dias <- lm(a_dias ~ PFHxSQ + mom_agey_bp + race + edu + bmi + smoke + income + parity, data = PFAS_DIAS_R)
summary(pfhxs_a_dias)
pfhxs_a_dias_l <- lm(a_dias ~ PFHxSQ + mom_agey_bp + race + edu + smoke + income + parity, data = subset(PFAS_DIAS_R, bmi_mom_prepreg_d < 25.0))
summary(pfhxs_a_dias_l)
pfhxs_a_dias_h <- lm(a_dias ~ PFHxSQ + mom_agey_bp + race + edu + smoke + income + parity, data = subset(PFAS_DIAS_R, bmi_mom_prepreg_d >= 25.0))
summary(pfhxs_a_dias_h)


etpfq_a_dias <- lm(a_dias ~ ETPFQ + mom_agey_bp + race + edu + bmi + smoke + income + parity, data = PFAS_DIAS_R)
summary(etpfq_a_dias)
etpfq_a_dias_l <- lm(a_dias ~ ETPFQ + mom_agey_bp + race + edu + smoke + income + parity, data = subset(PFAS_DIAS_R, bmi_mom_prepreg_d < 25.0))
summary(etpfq_a_dias_l)
etpfq_a_dias_h <- lm(a_dias ~ ETPFQ + mom_agey_bp + race + edu + smoke + income + parity, data = subset(PFAS_DIAS_R, bmi_mom_prepreg_d >= 25.0))
summary(etpfq_a_dias_h)


mepfq_a_dias <- lm(a_dias ~ MEPFQ + mom_agey_bp + race + edu + bmi + smoke + income + parity, data = PFAS_DIAS_R)
summary(mepfq_a_dias)
mepfq_a_dias_l <- lm(a_dias ~ MEPFQ + mom_agey_bp + race + edu + smoke + income + parity, data = subset(PFAS_DIAS_R, bmi_mom_prepreg_d < 25.0))
summary(mepfq_a_dias_l)
mepfq_a_dias_h <- lm(a_dias ~ MEPFQ + mom_agey_bp + race + edu + smoke + income + parity, data = subset(PFAS_DIAS_R, bmi_mom_prepreg_d >= 25.0))
summary(mepfq_a_dias_h)


```

```{r}
##### EFFECT OF PFAS ON SITAR VELOCITY PARAMTER ON SYSTOLIC BP

## RACE STRATA ##
pfoa_c_sys <- lm(c_sys ~ PFOAQ + mom_agey_bp + race + edu + bmi + smoke + income + parity, data = PFAS_SYS_R)
summary(pfoa_c_sys)
pfoa_c_sys_w <- lm(c_sys ~ PFOAQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_SYS_R, race == "White"))
summary(pfoa_c_sys_w)
pfoa_c_sys_b <- lm(c_sys ~ PFOAQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_SYS_R, race == "Black"))
summary(pfoa_c_sys_b)
pfoa_c_sys_h <- lm(c_sys ~ PFOAQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_SYS_R, race == "Hispanic"))
summary(pfoa_c_sys_h)
pfoa_c_sys_a <- lm(c_sys ~ PFOAQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_SYS_R, race == "Asian"))
summary(pfoa_c_sys_a)
pfoa_c_sys_om <- lm(c_sys ~ PFOAQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_SYS_R, race == "Other" | race == "More than 1 race"))
summary(pfoa_c_sys_om)


pfos_c_sys <- lm(c_sys ~ PFOSQ + mom_agey_bp + race + edu + bmi + smoke + income + parity, data = PFAS_SYS_R)
summary(pfos_c_sys)
pfos_c_sys_w <- lm(c_sys ~ PFOSQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_SYS_R, race == "White"))
summary(pfos_c_sys_w)
pfos_c_sys_b <- lm(c_sys ~ PFOSQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_SYS_R, race == "Black"))
summary(pfos_c_sys_b)
pfos_c_sys_h <- lm(c_sys ~ PFOSQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_SYS_R, race == "Hispanic"))
summary(pfos_c_sys_h)
pfos_c_sys_a <- lm(c_sys ~ PFOSQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_SYS_R, race == "Asian"))
summary(pfos_c_sys_a)
pfos_c_sys_om <- lm(c_sys ~ PFOSQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_SYS_R, race == "Other" | race == "More than 1 race"))
summary(pfos_c_sys_om)


pfna2_c_sys <- lm(c_sys ~ PFNAQ + mom_agey_bp + race + edu + bmi + smoke + income + parity, data = PFAS_SYS_R)
summary(pfna2_c_sys)
pfna2_c_sys_w <- lm(c_sys ~ PFNAQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_SYS_R, race == "White"))
summary(pfna2_c_sys_w)
pfna2_c_sys_b <- lm(c_sys ~ PFNAQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_SYS_R, race == "Black"))
summary(pfna2_c_sys_b)
pfna2_c_sys_h <- lm(c_sys ~ PFNAQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_SYS_R, race == "Hispanic"))
summary(pfna2_c_sys_h)
pfna2_c_sys_a <- lm(c_sys ~ PFNAQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_SYS_R, race == "Asian"))
summary(pfna2_c_sys_a)
pfna2_c_sys_om <- lm(c_sys ~ PFNAQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_SYS_R, race == "Other" | race == "More than 1 race"))
summary(pfna2_c_sys_om)


pfhxs_c_sys <- lm(c_sys ~ PFHxSQ + mom_agey_bp + race + edu + bmi + smoke + income + parity, data = PFAS_SYS_R)
summary(pfhxs_c_sys)
pfhxs_c_sys_w <- lm(c_sys ~ PFHxSQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_SYS_R, race == "White"))
summary(pfhxs_c_sys_w)
pfhxs_c_sys_b <- lm(c_sys ~ PFHxSQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_SYS_R, race == "Black"))
summary(pfhxs_c_sys_b)
pfhxs_c_sys_h <- lm(c_sys ~ PFHxSQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_SYS_R, race == "Hispanic"))
summary(pfhxs_c_sys_h)
pfhxs_c_sys_a <- lm(c_sys ~ PFHxSQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_SYS_R, race == "Asian"))
summary(pfhxs_c_sys_a)
pfhxs_c_sys_om <- lm(c_sys ~ PFHxSQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_SYS_R, race == "Other" | race == "More than 1 race"))
summary(pfhxs_c_sys_om)


etpfq_c_sys <- lm(c_sys ~ ETPFQ + mom_agey_bp + race + edu + bmi + smoke + income + parity, data = PFAS_SYS_R)
summary(etpfq_c_sys)
etpfq_c_sys_w <- lm(c_sys ~ ETPFQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_SYS_R, race == "White"))
summary(etpfq_c_sys_w)
etpfq_c_sys_b <- lm(c_sys ~ ETPFQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_SYS_R, race == "Black"))
summary(etpfq_c_sys_b)
etpfq_c_sys_h <- lm(c_sys ~ ETPFQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_SYS_R, race == "Hispanic"))
summary(etpfq_c_sys_h)
etpfq_c_sys_a <- lm(c_sys ~ ETPFQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_SYS_R, race == "Asian"))
summary(etpfq_c_sys_a)
etpfq_c_sys_om <- lm(c_sys ~ ETPFQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_SYS_R, race == "Other" | race == "More than 1 race"))
summary(etpfq_c_sys_om)


mepfq_c_sys <- lm(c_sys ~ MEPFQ + mom_agey_bp + race + edu + bmi + smoke + income + parity, data = PFAS_SYS_R)
summary(mepfq_c_sys)
mepfq_c_sys_w <- lm(c_sys ~ MEPFQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_SYS_R, race == "White"))
summary(mepfq_c_sys_w)
mepfq_c_sys_b <- lm(c_sys ~ MEPFQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_SYS_R, race == "Black"))
summary(mepfq_c_sys_b)
mepfq_c_sys_h <- lm(c_sys ~ MEPFQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_SYS_R, race == "Hispanic"))
summary(mepfq_c_sys_h)
mepfq_c_sys_a <- lm(c_sys ~ MEPFQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_SYS_R, race == "Asian"))
summary(mepfq_c_sys_a)
mepfq_c_sys_om <- lm(c_sys ~ MEPFQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_SYS_R, race == "Other" | race == "More than 1 race"))
summary(mepfq_c_sys_om)


## AGE STRATA ##

pfoa_c_sys <- lm(c_sys ~ PFOAQ + mom_agey_bp + race + edu + bmi + smoke + income + parity, data = PFAS_SYS_R)
summary(pfoa_c_sys)
pfoa_c_sys_y <- lm(c_sys ~ PFOAQ + race + edu + bmi + smoke + income + parity, data = subset(PFAS_SYS_R, mom_agey_bp < 35.0))
summary(pfoa_c_sys_y)
pfoa_c_sys_o <- lm(c_sys ~ PFOAQ + race + edu + bmi + smoke + income + parity, data = subset(PFAS_SYS_R, mom_agey_bp >= 35.0))
summary(pfoa_c_sys_o)


pfos_c_sys <- lm(c_sys ~ PFOSQ + mom_agey_bp + race + edu + bmi + smoke + income + parity, data = PFAS_SYS_R)
summary(pfos_c_sys)
pfos_c_sys_y <- lm(c_sys ~ PFOSQ + race + edu + bmi + smoke + income + parity, data = subset(PFAS_SYS_R, mom_agey_bp < 35.0))
summary(pfos_c_sys_y)
pfos_c_sys_o <- lm(c_sys ~ PFOSQ + race + edu + bmi + smoke + income + parity, data = subset(PFAS_SYS_R, mom_agey_bp >= 35.0))
summary(pfos_c_sys_o)


pfna2_c_sys <- lm(c_sys ~ PFNAQ + mom_agey_bp + race + edu + bmi + smoke + income + parity, data = PFAS_SYS_R)
summary(pfna2_c_sys)
pfna2_c_sys_y <- lm(c_sys ~ PFNAQ + race + edu + bmi + smoke + income + parity, data = subset(PFAS_SYS_R, mom_agey_bp < 35.0))
summary(pfna2_c_sys_y)
pfna2_c_sys_o <- lm(c_sys ~ PFNAQ + race + edu + bmi + smoke + income + parity, data = subset(PFAS_SYS_R, mom_agey_bp >= 35.0))
summary(pfna2_c_sys_o)


pfhxs_c_sys <- lm(c_sys ~ PFHxSQ + mom_agey_bp + race + edu + bmi + smoke + income + parity, data = PFAS_SYS_R)
summary(pfhxs_c_sys)
pfhxs_c_sys_y <- lm(c_sys ~ PFHxSQ + race + edu + bmi + smoke + income + parity, data = subset(PFAS_SYS_R, mom_agey_bp < 35.0))
summary(pfhxs_c_sys_y)
pfhxs_c_sys_o <- lm(c_sys ~ PFHxSQ + race + edu + bmi + smoke + income + parity, data = subset(PFAS_SYS_R, mom_agey_bp >= 35.0))
summary(pfhxs_c_sys_o)


etpfq_c_sys <- lm(c_sys ~ ETPFQ + mom_agey_bp + race + edu + bmi + smoke + income + parity, data = PFAS_SYS_R)
summary(etpfq_c_sys)
etpfq_c_sys_y <- lm(c_sys ~ ETPFQ + race + edu + bmi + smoke + income + parity, data = subset(PFAS_SYS_R, mom_agey_bp < 35.0))
summary(etpfq_c_sys_y)
etpfq_c_sys_o <- lm(c_sys ~ ETPFQ + race + edu + bmi + smoke + income + parity, data = subset(PFAS_SYS_R, mom_agey_bp >= 35.0))
summary(etpfq_c_sys_o)


mepfq_c_sys <- lm(c_sys ~ MEPFQ + mom_agey_bp + race + edu + bmi + smoke + income + parity, data = PFAS_SYS_R)
summary(mepfq_c_sys)
mepfq_c_sys_y <- lm(c_sys ~ MEPFQ + race + edu + bmi + smoke + income + parity, data = subset(PFAS_SYS_R, mom_agey_bp < 35.0))
summary(mepfq_c_sys_y)
mepfq_c_sys_o <- lm(c_sys ~ MEPFQ + race + edu + bmi + smoke + income + parity, data = subset(PFAS_SYS_R, mom_agey_bp >= 35.0))
summary(mepfq_c_sys_o)


## BMI STRATA ##

pfoa_c_sys <- lm(c_sys ~ PFOAQ + mom_agey_bp + race + edu + bmi + smoke + income + parity, data = PFAS_SYS_R)
summary(pfoa_c_sys)
pfoa_c_sys_l <- lm(c_sys ~ PFOAQ + mom_agey_bp + race + edu + smoke + income + parity, data = subset(PFAS_SYS_R, bmi_mom_prepreg_d < 25.0))
summary(pfoa_c_sys_l)
pfoa_c_sys_h <- lm(c_sys ~ PFOAQ + mom_agey_bp + race + edu + smoke + income + parity, data = subset(PFAS_SYS_R, bmi_mom_prepreg_d >= 25.0))
summary(pfoa_c_sys_h)


pfos_c_sys <- lm(c_sys ~ PFOSQ + mom_agey_bp + race + edu + bmi + smoke + income + parity, data = PFAS_SYS_R)
summary(pfos_c_sys)
pfos_c_sys_l <- lm(c_sys ~ PFOSQ + mom_agey_bp + race + edu + smoke + income + parity, data = subset(PFAS_SYS_R, bmi_mom_prepreg_d < 25.0))
summary(pfos_c_sys_l)
pfos_c_sys_h <- lm(c_sys ~ PFOSQ + mom_agey_bp + race + edu + smoke + income + parity, data = subset(PFAS_SYS_R, bmi_mom_prepreg_d >= 25.0))
summary(pfos_c_sys_h)


pfna2_c_sys <- lm(c_sys ~ PFNAQ + mom_agey_bp + race + edu + bmi + smoke + income + parity, data = PFAS_SYS_R)
summary(pfna2_c_sys)
pfna2_c_sys_l <- lm(c_sys ~ PFNAQ + mom_agey_bp + race + edu + smoke + income + parity, data = subset(PFAS_SYS_R, bmi_mom_prepreg_d < 25.0))
summary(pfna2_c_sys_l)
pfna2_c_sys_h <- lm(c_sys ~ PFNAQ + mom_agey_bp + race + edu + smoke + income + parity, data = subset(PFAS_SYS_R, bmi_mom_prepreg_d >= 25.0))
summary(pfna2_c_sys_h)


pfhxs_c_sys <- lm(c_sys ~ PFHxSQ + mom_agey_bp + race + edu + bmi + smoke + income + parity, data = PFAS_SYS_R)
summary(pfhxs_c_sys)
pfhxs_c_sys_l <- lm(c_sys ~ PFHxSQ + mom_agey_bp + race + edu + smoke + income + parity, data = subset(PFAS_SYS_R, bmi_mom_prepreg_d < 25.0))
summary(pfhxs_c_sys_l)
pfhxs_c_sys_h <- lm(c_sys ~ PFHxSQ + mom_agey_bp + race + edu + smoke + income + parity, data = subset(PFAS_SYS_R, bmi_mom_prepreg_d >= 25.0))
summary(pfhxs_c_sys_h)


etpfq_c_sys <- lm(c_sys ~ ETPFQ + mom_agey_bp + race + edu + bmi + smoke + income + parity, data = PFAS_SYS_R)
summary(etpfq_c_sys)
etpfq_c_sys_l <- lm(c_sys ~ ETPFQ + mom_agey_bp + race + edu + smoke + income + parity, data = subset(PFAS_SYS_R, bmi_mom_prepreg_d < 25.0))
summary(etpfq_c_sys_l)
etpfq_c_sys_h <- lm(c_sys ~ ETPFQ + mom_agey_bp + race + edu + smoke + income + parity, data = subset(PFAS_SYS_R, bmi_mom_prepreg_d >= 25.0))
summary(etpfq_c_sys_h)


mepfq_c_sys <- lm(c_sys ~ MEPFQ + mom_agey_bp + race + edu + bmi + smoke + income + parity, data = PFAS_SYS_R)
summary(mepfq_c_sys)
mepfq_c_sys_l <- lm(c_sys ~ MEPFQ + mom_agey_bp + race + edu + smoke + income + parity, data = subset(PFAS_SYS_R, bmi_mom_prepreg_d < 25.0))
summary(mepfq_c_sys_l)
mepfq_c_sys_h <- lm(c_sys ~ MEPFQ + mom_agey_bp + race + edu + smoke + income + parity, data = subset(PFAS_SYS_R, bmi_mom_prepreg_d >= 25.0))
summary(mepfq_a_dias_h)

```


```{r}
## EFFECT OF PFAS ON SITAR VELOCITY PARAMTER ON DIASTOLIC BP

## RACE STRATA ##
pfoa_c_dias <- lm(c_dias ~ PFOAQ + mom_agey_bp + race + edu + bmi + smoke + income + parity, data = PFAS_DIAS_R)
summary(pfoa_c_dias)
pfoa_c_dias_w <- lm(c_dias ~ PFOAQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_DIAS_R, race == "White"))
summary(pfoa_c_dias_w)
pfoa_c_diass_b <- lm(c_dias ~ PFOAQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_DIAS_R, race == "Black"))
summary(pfoa_c_diass_b)
pfoa_c_diass_h <- lm(c_dias ~ PFOAQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_DIAS_R, race == "Hispanic"))
summary(pfoa_c_diass_h)
pfoa_c_diass_a <- lm(c_dias ~ PFOAQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_DIAS_R, race == "Asian"))
summary(pfoa_c_diass_a)
pfoa_c_diass_om <- lm(c_dias ~ PFOAQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_DIAS_R, race == "Other" | race == "More than 1 race"))
summary(pfoa_c_diass_om)


pfos_c_dias <- lm(c_dias ~ PFOSQ + mom_agey_bp + race + edu + bmi + smoke + income + parity, data = PFAS_DIAS_R)
summary(pfos_c_dias)
pfos_c_dias_w <- lm(c_dias ~ PFOSQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_DIAS_R, race == "White"))
summary(pfos_c_dias_w)
pfos_c_diass_b <- lm(c_dias ~ PFOSQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_DIAS_R, race == "Black"))
summary(pfos_c_diass_b)
pfos_c_diass_h <- lm(c_dias ~ PFOSQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_DIAS_R, race == "Hispanic"))
summary(pfos_c_diass_h)
pfos_c_diass_a <- lm(c_dias ~ PFOSQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_DIAS_R, race == "Asian"))
summary(pfos_c_diass_a)
pfos_c_diass_om <- lm(c_dias ~ PFOSQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_DIAS_R, race == "Other" | race == "More than 1 race"))
summary(pfos_c_diass_om)


pfna2_c_dias <- lm(c_dias ~ PFNAQ + mom_agey_bp + race + edu + bmi + smoke + income + parity, data = PFAS_DIAS_R)
summary(pfna2_c_dias)
pfna2_c_dias_w <- lm(c_dias ~ PFNAQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_DIAS_R, race == "White"))
summary(pfna2_c_dias_w)
pfna2_c_dias_b <- lm(c_dias ~ PFNAQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_DIAS_R, race == "Black"))
summary(pfna2_c_dias_b)
pfna2_c_dias_h <- lm(c_dias ~ PFNAQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_DIAS_R, race == "Hispanic"))
summary(pfna2_c_dias_h)
pfna2_c_dias_a <- lm(c_dias ~ PFNAQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_DIAS_R, race == "Asian"))
summary(pfna2_c_dias_a)
pfna2_c_dias_om <- lm(c_dias ~ PFNAQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_DIAS_R, race == "Other" | race == "More than 1 race"))
summary(pfna2_c_diass_om)


pfhxs_c_dias <- lm(c_dias ~ PFHxSQ + mom_agey_bp + race + edu + bmi + smoke + income + parity, data = PFAS_DIAS_R)
summary(pfhxs_c_dias)
pfhxs_c_dias_w <- lm(c_dias ~ PFHxSQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_DIAS_R, race == "White"))
summary(pfhxs_c_dias_w)
pfhxs_c_diass_b <- lm(c_dias ~ PFHxSQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_DIAS_R, race == "Black"))
summary(pfhxs_c_diass_b)
pfhxs_c_diass_h <- lm(c_dias ~ PFHxSQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_DIAS_R, race == "Hispanic"))
summary(pfhxs_c_diass_h)
pfhxs_c_diass_a <- lm(c_dias ~ PFHxSQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_DIAS_R, race == "Asian"))
summary(pfhxs_c_diass_a)
pfhxs_c_diass_om <- lm(c_dias ~ PFHxSQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_DIAS_R, race == "Other" | race == "More than 1 race"))
summary(pfhxs_c_diass_om)


etpfq_c_dias <- lm(c_dias ~ ETPFQ + mom_agey_bp + race + edu + bmi + smoke + income + parity, data = PFAS_DIAS_R)
summary(etpfq_c_dias)
etpfq_c_dias_w <- lm(c_dias ~ ETPFQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_DIAS_R, race == "White"))
summary(etpfq_c_dias_w)
etpfq_c_diass_b <- lm(c_dias ~ ETPFQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_DIAS_R, race == "Black"))
summary(etpfq_c_diass_b)
etpfq_c_diass_h <- lm(c_dias ~ ETPFQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_DIAS_R, race == "Hispanic"))
summary(etpfq_c_diass_h)
etpfq_c_diass_a <- lm(c_dias ~ ETPFQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_DIAS_R, race == "Asian"))
summary(etpfq_c_diass_a)
etpfq_c_diass_om <- lm(c_dias ~ ETPFQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_DIAS_R, race == "Other" | race == "More than 1 race"))
summary(etpfq_c_diass_om)


mepfq_c_dias <- lm(c_dias ~ MEPFQ + mom_agey_bp + race + edu + bmi + smoke + income + parity, data = PFAS_DIAS_R)
summary(mepfq_c_dias)
mepfq_c_dias_w <- lm(c_dias ~ MEPFQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_DIAS_R, race == "White"))
summary(mepfq_c_dias_w)
mepfq_c_diass_b <- lm(c_dias ~ MEPFQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_DIAS_R, race == "Black"))
summary(mepfq_c_diass_b)
mepfq_c_diass_h <- lm(c_dias ~ MEPFQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_DIAS_R, race == "Hispanic"))
summary(mepfq_c_diass_h)
mepfq_c_diass_a <- lm(c_dias ~ MEPFQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_DIAS_R, race == "Asian"))
summary(mepfq_c_diass_a)
mepfq_c_diass_om <- lm(c_dias ~ MEPFQ + mom_agey_bp + edu + bmi + smoke + income + parity, data = subset(PFAS_DIAS_R, race == "Other" | race == "More than 1 race"))
summary(mepfq_c_diass_om)


## AGE STRATA ##

pfoa_c_dias <- lm(c_dias ~ PFOAQ + mom_agey_bp + race + edu + bmi + smoke + income + parity, data = PFAS_DIAS_R)
summary(pfoa_c_dias)
pfoa_c_dias_y <- lm(c_dias ~ PFOAQ + race + edu + bmi + smoke + income + parity, data = subset(PFAS_DIAS_R, mom_agey_bp < 35.0))
summary(pfoa_c_dias_y)
pfoa_c_dias_o <- lm(c_dias ~ PFOAQ + race + edu + bmi + smoke + income + parity, data = subset(PFAS_DIAS_R, mom_agey_bp >= 35.0))
summary(pfoa_c_dias_o)


pfos_c_dias <- lm(c_dias ~ PFOSQ + mom_agey_bp + race + edu + bmi + smoke + income + parity, data = PFAS_DIAS_R)
summary(pfos_c_dias)
pfos_c_dias_y <- lm(c_dias ~ PFOSQ + race + edu + bmi + smoke + income + parity, data = subset(PFAS_DIAS_R, mom_agey_bp < 35.0))
summary(pfos_c_dias_y)
pfos_c_dias_o <- lm(c_dias ~ PFOSQ + race + edu + bmi + smoke + income + parity, data = subset(PFAS_DIAS_R, mom_agey_bp >= 35.0))
summary(pfos_c_dias_o)


pfna2_c_dias <- lm(c_dias ~ PFNAQ + mom_agey_bp + race + edu + bmi + smoke + income + parity, data = PFAS_DIAS_R)
summary(pfna2_c_dias)
pfna2_c_dias_y <- lm(c_dias ~ PFNAQ + race + edu + bmi + smoke + income + parity, data = subset(PFAS_DIAS_R, mom_agey_bp < 35.0))
summary(pfna2_c_dias_y)
pfna2_c_dias_o <- lm(c_dias ~ PFNAQ + race + edu + bmi + smoke + income + parity, data = subset(PFAS_DIAS_R, mom_agey_bp >= 35.0))
summary(pfna2_c_dias_o)


pfhxs_c_dias <- lm(c_dias ~ PFHxSQ + mom_agey_bp + race + edu + bmi + smoke + income + parity, data = PFAS_DIAS_R)
summary(pfhxs_c_dias)
pfhxs_c_dias_y <- lm(c_dias ~ PFHxSQ + race + edu + bmi + smoke + income + parity, data = subset(PFAS_DIAS_R, mom_agey_bp < 35.0))
summary(pfhxs_c_dias_y)
pfhxs_c_dias_o <- lm(c_dias ~ PFHxSQ + race + edu + bmi + smoke + income + parity, data = subset(PFAS_DIAS_R, mom_agey_bp >= 35.0))
summary(pfhxs_c_dias_o)


etpfq_c_dias <- lm(c_dias ~ ETPFQ + mom_agey_bp + race + edu + bmi + smoke + income + parity, data = PFAS_DIAS_R)
summary(etpfq_c_dias)
etpfq_c_dias_y <- lm(c_dias ~ ETPFQ + race + edu + bmi + smoke + income + parity, data = subset(PFAS_DIAS_R, mom_agey_bp < 35.0))
summary(etpfq_c_dias_y)
etpfq_c_dias_o <- lm(c_dias ~ ETPFQ + race + edu + bmi + smoke + income + parity, data = subset(PFAS_DIAS_R, mom_agey_bp >= 35.0))
summary(etpfq_c_dias_o)


mepfq_c_dias <- lm(c_dias ~ MEPFQ + mom_agey_bp + race + edu + bmi + smoke + income + parity, data = PFAS_DIAS_R)
summary(mepfq_c_dias)
mepfq_c_dias_y <- lm(c_dias ~ MEPFQ + race + edu + bmi + smoke + income + parity, data = subset(PFAS_DIAS_R, mom_agey_bp < 35.0))
summary(mepfq_c_dias_y)
mepfq_c_dias_o <- lm(c_dias ~ MEPFQ + race + edu + bmi + smoke + income + parity, data = subset(PFAS_DIAS_R, mom_agey_bp >= 35.0))
summary(mepfq_c_dias_o)


## BMI STRATA ##

pfoa_c_dias <- lm(c_dias ~ PFOAQ + mom_agey_bp + race + edu + bmi + smoke + income + parity, data = PFAS_DIAS_R)
summary(pfoa_c_dias)
pfoa_c_dias_l <- lm(c_dias ~ PFOAQ + mom_agey_bp + race + edu + smoke + income + parity, data = subset(PFAS_DIAS_R, bmi_mom_prepreg_d < 25.0))
summary(pfoa_c_dias_l)
pfoa_c_dias_h <- lm(c_dias ~ PFOAQ + mom_agey_bp + race + edu + smoke + income + parity, data = subset(PFAS_DIAS_R, bmi_mom_prepreg_d >= 25.0))
summary(pfoa_c_dias_h)


pfos_c_dias <- lm(c_dias ~ PFOSQ + mom_agey_bp + race + edu + bmi + smoke + income + parity, data = PFAS_DIAS_R)
summary(pfos_c_dias)
pfos_c_dias_l <- lm(c_dias ~ PFOSQ + mom_agey_bp + race + edu + smoke + income + parity, data = subset(PFAS_DIAS_R, bmi_mom_prepreg_d < 25.0))
summary(pfos_c_dias_l)
pfos_c_dias_h <- lm(c_dias ~ PFOSQ + mom_agey_bp + race + edu + smoke + income + parity, data = subset(PFAS_DIAS_R, bmi_mom_prepreg_d >= 25.0))
summary(pfos_c_dias_h)


pfna2_c_dias <- lm(c_dias ~ PFNAQ + mom_agey_bp + race + edu + bmi + smoke + income + parity, data = PFAS_DIAS_R)
summary(pfna2_c_dias)
pfna2_c_dias_l <- lm(c_dias ~ PFNAQ + mom_agey_bp + race + edu + smoke + income + parity, data = subset(PFAS_DIAS_R, bmi_mom_prepreg_d < 25.0))
summary(pfna2_c_dias_l)
pfna2_c_dias_h <- lm(c_dias ~ PFNAQ + mom_agey_bp + race + edu + smoke + income + parity, data = subset(PFAS_DIAS_R, bmi_mom_prepreg_d >= 25.0))
summary(pfna2_c_dias_h)


pfhxs_c_dias <- lm(c_dias ~ PFHxSQ + mom_agey_bp + race + edu + bmi + smoke + income + parity, data = PFAS_DIAS_R)
summary(pfhxs_c_dias)
pfhxs_c_dias_l <- lm(c_dias ~ PFHxSQ + mom_agey_bp + race + edu + smoke + income + parity, data = subset(PFAS_DIAS_R, bmi_mom_prepreg_d < 25.0))
summary(pfhxs_c_dias_l)
pfhxs_c_dias_h <- lm(c_dias ~ PFHxSQ + mom_agey_bp + race + edu + smoke + income + parity, data = subset(PFAS_DIAS_R, bmi_mom_prepreg_d >= 25.0))
summary(pfhxs_c_dias_h)


etpfq_c_dias <- lm(c_dias ~ ETPFQ + mom_agey_bp + race + edu + bmi + smoke + income + parity, data = PFAS_DIAS_R)
summary(etpfq_c_dias)
etpfq_c_dias_l <- lm(c_dias ~ ETPFQ + mom_agey_bp + race + edu + smoke + income + parity, data = subset(PFAS_DIAS_R, bmi_mom_prepreg_d < 25.0))
summary(etpfq_c_dias_l)
etpfq_c_dias_h <- lm(c_dias ~ ETPFQ + mom_agey_bp + race + edu + smoke + income + parity, data = subset(PFAS_DIAS_R, bmi_mom_prepreg_d >= 25.0))
summary(etpfq_c_dias_h)


mepfq_c_dias <- lm(c_dias ~ MEPFQ + mom_agey_bp + race + edu + bmi + smoke + income + parity, data = PFAS_DIAS_R)
summary(mepfq_c_dias)
mepfq_c_dias_l <- lm(c_dias ~ MEPFQ + mom_agey_bp + race + edu + smoke + income + parity, data = subset(PFAS_DIAS_R, bmi_mom_prepreg_d < 25.0))
summary(mepfq_c_dias_l)
mepfq_c_dias_h <- lm(c_dias ~ MEPFQ + mom_agey_bp + race + edu + smoke + income + parity, data = subset(PFAS_DIAS_R, bmi_mom_prepreg_d >= 25.0))
summary(mepfq_c_dias_h)


```


```{r}
##### ADJUSTED PFAS PLOTS FOR SYSTOLIC BP

#par(mfrow=c(3,2))

mplot(x = ga_w_sc, y = ga_w_sys, id = aid, col = aid, data = PFAS_SYS,
      ylim=c(70, 200), xlim=c(-20, 20), las = 1,
      ylab="Systolic Blood Pressure (mmHg)", xlab="Gestational age (weeks)",
      main="Undjusted SITAR Curve (PFAS), Systolic BP")

plot(m11_w, opt = 'a', las = 1, col = aid, ylim = c(70, 200), xlim=c(-20, 20),
     ylab="Systolic Blood Pressure (mmHg)", xlab="Gestational age (weeks)",
     main="Adjusted SITAR Curve (PFAS), Systolic BP")

```

```{r}
##### ADJUSTED PFAS PLOTS FOR DIASTOLIC BP

#par(mfrow=c(3,2))

mplot(x = ga_w_sc, y = ga_w_dias, id = aid, col = aid, data = PFAS_DIAS,
      las = 1, ylim = c(30, 110), xlim=c(-20, 20),
      ylab="Diastolic Blood Pressure (mmHg)", xlab="Gestational age (weeks)",
      main="Undjusted SITAR Curve (PFAS), Diastolic BP")

plot(m21_w, opt = 'a', las = 1, col = aid, ylim = c(30, 110), xlim=c(-20, 20),
     ylab="Diastolic Blood Pressure (mmHg)", xlab="Gestational age (weeks)",
     main="Adjusted SITAR Curve (PFAS), Diastolic BP")

```

```{r}
##### PFAS MEAN VELOCITY PLOTS FOR DIASTOLIC BP

#par(mfrow=c(2,2))

#plot(m11_w, opt = 'd', las = 2, col = "deepskyblue3", lwd = 0.5, 
#     ylab="Mean systolic blood pressure (mmHg)", xlab="Gestational age (weeks)",
#     main="Mean Growth Curve, Systolic BP", ylim=c(110,120))
#legend('topleft', c('Size: a = 6.39', 'Intensity: c = 8.17^-14'), cex = 0.8, inset=0.04)

#plot(m11_w, opt = 'v', las = 2, col = "deepskyblue3", lwd = 0.5, 
#     ylab="Mean systolic blood pressure velocity (mmHg/day)", xlab="Gestational age (weeks)",
#     main="Mean Velocity Curve,Systolic BP")

plot(m11_w, opt = 'dv', las = 2, col = "deepskyblue3", lwd = 2, ylim = c(105, 122), xlim=c(-20, 20),
     ylab="Mean systolic blood pressure (mmHg)", xlab="Gestational age (weeks)",
     main="Mean Growth/Velocity Curve, Systolic BP")
legend('bottomright', c('Size: a = 6.39', 'Intensity: c = 8.17^-14'), cex = 0.8, inset=0.04)



#plot(m21_w, opt = 'd', las = 2, col = "red1", lwd = 0.5, 
#     ylab="Mean diastolic blood pressure (mmHg/day)", xlab="Gestational age (weeks)",
#     main="Mean Growth Curve, Diastolic BP")

#plot(m21_w, opt = 'v', las = 2, col = "red1", lwd = 0.5, 
#     ylab="Mean diastolic blood pressure velocity (mmHg/day)", xlab="Gestational age (weeks)",
#     main="Mean Velocity Curve, Diastolic BP")

plot(m21_w, opt = 'dv', las = 2, col = "red1", lwd = 2, ylim = c(65, 77), xlim=c(-20, 20),
     ylab="Mean diastolic blood pressure (mmHg)", xlab="Gestational age (weeks)",
     main="Mean Growth/Velocity Curve, Diastolic BP")
legend('bottomright', c('Size: a = 4.69', 'Intensity: c = 8.30^-5'), cex = 0.8, inset=0.04)



```

```{r}
plot(m11_w, opt = 'd', las = 1, apv = F, legend = NULL, vlab="", ylab="", xlab="", ylim=c(100,130),
     y2par=list(lwd=2), lwd=2, main="Mean Growth Curve (95% CI), Sytolic BP")
lines(m11_w, opt = 'd', y2par=list(col='lightblue', lwd=2), 
      apv = F, lwd=2, lty=2, col='deeppink', abc=(sqrt(diag(getVarCov(m11_w)))*1.96))
lines(m11_w, opt = 'd', y2par=list(col='turquoise3', lwd=2), 
      apv = F, lwd=1.5, col='gray', lty=2, abc=(sqrt(diag(getVarCov(m11_w)))*1.96))
mtext("Gestational age (weeks)", side = 1, line = 3)
mtext("Systolic BP (mmHg)", side = 2, line = 3.5)
mtext("Systolic BP velocity (mmHg/day)", side = 4, line = 3)
legend("topleft", c( "Mean growth curve", "Upper 95% CI", "Lower 95% CI"),
       col = c("black", "deeppink", "turquoise3"), lty = c(1, 2, 2))


plot(m21_w, opt = 'd', las = 1, apv = F, legend = NULL, vlab="", ylab="", xlab="", ylim=c(60,85),
     y2par=list(lwd=2), lwd=2, main="Mean Growth Curve (95% CI), Diastolic BP")
lines(m21_w, opt = 'd', y2par=list(col='lightblue', lwd=2), 
      apv = F, lwd=2, lty=2, col='deeppink', abc=(sqrt(diag(getVarCov(m21_w)))*1.96))
lines(m21_w, opt = 'd', y2par=list(col='turquoise3', lwd=2), 
      apv = F, lwd=1.5, col='gray', lty=2, abc=(sqrt(diag(getVarCov(m21_w)))*1.96))
mtext("Gestational age (weeks)", side = 1, line = 3)
mtext("Diastolic BP (mmHg)", side = 2, line = 3.5)
mtext("Diastolic BP velocity (mmHg/day)", side = 4, line = 3)
legend("topleft", c( "Mean growth curve", "Upper 95% CI", "Lower 95% CI"),
       col = c("black", "deeppink", "turquoise3"), lty = c(1, 2, 2))

```


```{r}
##### GOODNESS OF FIT PLOT FOR SITAR MODELLING

plot(sys ~ ga_m_sc, data = PFAS_SYS_R, las=1, col="lightblue4", ylim=c(105,120),
     ylab="Systolic BP (mmHg)", xlab="Gestational age (weeks)",
     main="Goodness of Fit Plot for SITAR Modelling, Systolic BP")
lines(m11_w, opt = 'D', lty = 1, col="darkturquoise")
lines(m11_w, opt = 'd', lwd=2, lty = 2, col="red")
legend('bottomright', c('Observed', 'Fitted trajectory', 'SITAR average-mean curve'),
       lty = c(NA, 1, 2), pch=c(1,NA,NA), col = c("lightblue4", "turquoise4", "red"), cex = 0.8, inset=0.04) 


plot(dias ~ ga_m_sc, data = PFAS_DIAS_R, las=1, col="lightblue4", ylim=c(65,75),
     ylab="Diastolic BP (mmHg)", xlab="Gestational age (weeks)",
     main="Goodness of Fit Plot for SITAR Modelling, Diastolic BP")
lines(m21_w, opt = 'D', lty = 1, col="darkturquoise")
lines(m21_w, opt = 'd', lwd=2, lty = 2, col="red")
legend('bottomright', c('Observed', 'Fitted trajectory', 'SITAR average-mean curve'),
       lty = c(NA, 1, 2), pch=c(1,NA,NA), col = c("lightblue4", "turquoise4", "red"), cex = 0.8, inset=0.04) 


```

```{r}
##### SAVE DATASET
PFAS_SYS_BKMR <- subset(PFAS_SYS_R, dup_num == 1)
PFAS_SYS_BKMR <- subset(PFAS_SYS_BKMR, select = c(aid, a_sys, b_sys, c_sys, ga_days, ga_days_sc, ga_w, ga_w_sc, ga_w_sys,
                                                  mod_pre_d, vig_pre_d, exer, race2_mom_epi_epia_d, coll_grad, 
                                                  bmi_mom_prepreg_d, smokpreg_final_d, income_hh_epq_epqa_d, parity_d, 
                                                  mom_agey_bp, dup_num, count, LNPFOS, LNPFOA, LNPFNA2, LNPFHxS, LNME_PFOSA, LNET_PFOSA))

PFAS_DIAS_BKMR <- subset(PFAS_DIAS_R, dup_num == 1)
PFAS_DIAS_BKMR <- subset(PFAS_DIAS_BKMR, select = c(aid, a_dias, b_dias, c_dias, ga_days, ga_days_sc, ga_w, ga_w_sc,
                                                  ga_w_dias, mod_pre_d, vig_pre_d, exer, race2_mom_epi_epia_d, coll_grad, 
                                                  bmi_mom_prepreg_d, smokpreg_final_d, income_hh_epq_epqa_d, parity_d, 
                                                  mom_agey_bp, dup_num, count, LNPFOS, LNPFOA, LNPFNA2, LNPFHxS, LNME_PFOSA, LNET_PFOSA))


write.csv(PFAS_DIAS_BKMR, file = "C:/Users/jarva/Desktop/James-Todd Lab/PFASandBP/Data/PFAS_DIAS_BKMR.csv", row.names = T)
write.csv(PFAS_SYS_BKMR, file = "C:/Users/jarva/Desktop/James-Todd Lab/PFASandBP/Data/PFAS_SYS_BKMR.csv", row.names = T)

```


