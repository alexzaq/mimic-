
library(dplyr)
library(mice)
base <- read.csv("基线数据及随访时间处理完.csv", header = TRUE, sep = ",")
a1 <- read.csv("neutrophils.csv", header = TRUE, sep = ",")
a2 <- read.csv("cerebral_seq.csv", header = TRUE, sep = ",")
a3 <- read.csv("albumin.csv", header = TRUE, sep = ",")
a4 <- read.csv("creatinine.csv", header = TRUE, sep = ",")
a5 <- read.csv("hemoglobin.csv", header = TRUE, sep = ",")
a6 <- read.csv("white_blood_cells.csv", header = TRUE, sep = ",")
a7 <- read.csv("platelet_count.csv", header = TRUE, sep = ",")
a8 <- read.csv("glucose.csv",header = TRUE, sep = ",")
a9<-read.csv("hypertension.csv",header = TRUE, sep = ",")
a10<-read.csv("PT.csv",header = TRUE, sep = ",")
b1<-read.csv("ldh.csv",header = TRUE, sep = ",")
b2<-read.csv("ptt.csv",header = TRUE, sep = ",")
b3<-read.csv("potassium.csv",header = TRUE, sep = ",")
b4<-read.csv("fe.csv",header = TRUE, sep = ",")
b6<-read.csv("cholesterol.csv",header = TRUE, sep = ",")
b7<-read.csv("calcium.csv",header = TRUE, sep = ",")
b8<-read.csv("alt.csv",header = TRUE, sep = ",")
b9<-read.csv("ast.csv",header = TRUE, sep = ",")
a9 <- a9 %>%
  distinct(subject_id1, .keep_all = TRUE)
total <- base %>%
  left_join(a1, by = c("hadm_id1" = "hadm_id1")) %>%
  left_join(a2, by = c("hadm_id1" = "hadm_id1")) %>%
  left_join(a3, by = c("hadm_id1" = "hadm_id1")) %>%
  left_join(a4, by = c("hadm_id1" = "hadm_id1")) %>%
  left_join(a5, by = c("hadm_id1" = "hadm_id1")) %>%
  left_join(a6, by = c("hadm_id1" = "hadm_id1")) %>%
  left_join(a7, by = c("hadm_id1" = "hadm_id1")) %>%
  left_join(a8, by = c("hadm_id1" = "hadm_id1")) %>%
  left_join(a10, by = c("hadm_id1" = "hadm_id1")) %>%
  left_join(b1, by = c("hadm_id1" = "hadm_id1")) %>%
  left_join(b2, by = c("hadm_id1" = "hadm_id1")) %>%
  left_join(b3, by = c("hadm_id1" = "hadm_id1")) %>%
  left_join(b4, by = c("hadm_id1" = "hadm_id1")) %>%
  left_join(b6, by = c("hadm_id1" = "hadm_id1")) %>%
  left_join(b7, by = c("hadm_id1" = "hadm_id1")) %>%
  left_join(b8, by = c("hadm_id1" = "hadm_id1")) %>%
  left_join(b9, by = c("hadm_id1" = "hadm_id1")) %>%
  left_join(a9, by = c("subject_id1" = "subject_id1"))%>%
  dplyr::select(
    subject_id1, hadm_id1 = hadm_id1.x, icu_intime, icu_outtime, 
    age, gender, dod,time_icuD,time_30d,time_90d,time_365d,
    neutrophils, albumin,                
    sbp_mean, dbp_mean, gcs_min,         
    sapsii,stau_icu,
    stau_30d, stau_90d, stau_365d,       
    icd_code, charlson_comorbidity_index,
    PLT,WBC,Hb,Cr,weight,height,
    LDH,PTT,potassium,Fe,Chol,CA,ALT,AST
    ,severe_liver_disease,
    malignant_cancer,aids,diabetes,
    hypertension,PT
  )
total <- total %>%
 mutate(hypertension = ifelse(!is.na(hypertension), 1, 0))
write.csv(total,file = "total.csv")
rm(list = ls()) #清空工作盘
total <- read.csv("total.csv")
total <- total %>%
  # 按入院时间排序（假设保留最后一次入院记录，若需保留第一次入院记录则用升序）
  arrange(desc(icu_intime)) %>% 
  distinct(subject_id1, .keep_all = TRUE)
#筛选icu时间小于1天
cleant2 <- total %>%
  mutate(
    icu_intime = as.POSIXct(icu_intime, format = "%Y/%m/%d %H:%M:%S"),
    icu_outtime = as.POSIXct(icu_outtime, format = "%Y/%m/%d %H:%M:%S"),
    icu_duration = as.numeric(difftime(icu_outtime, icu_intime, units = "hours"))
  ) %>%
  filter(!is.na(icd_code) & icu_duration >= 24) %>%
  dplyr::select(-icu_duration, -subject_id1,-hadm_id1)
cleant2_reduced <- cleant2 %>%
  mutate(
    stau_30d = factor(stau_30d, levels = c(0,1), labels = c("no","yes")),
    stau_90d = factor(stau_90d, levels = c(0,1), labels = c("no","yes")),
    stau_365d = factor(stau_365d, levels = c(0,1), labels = c("no","yes")),
    gender = factor(gender)
  )
sum(is.na(cleant2_reduced$dod))
cleant2_reduced <- cleant2_reduced %>%
  mutate(
    icu_duration_days = time_icuD,
    survival_time = time_365d,
    status = stau_365d,
    survival_time = pmax(pmin(survival_time, 365), 0.04)
  )
cleant2_reduced <- cleant2_reduced %>%
  mutate(
    icu_mortality = stau_icu
    )

names(cleant2_reduced)
cleant2_reduced <- cleant2_reduced%>%
  dplyr::select( -icu_duration_days,-icd_code,
                 -icu_outtime,-icu_intime,-dod,-X,)
cleant2_reduced <- cleant2_reduced %>%
  mutate(
    neutrophils = as.numeric(gsub("[^0-9.]", "", as.character(neutrophils))),
    albumin = as.numeric(gsub("[^0-9.]", "", as.character(albumin)))
  )
cleant2_reduced <- cleant2_reduced %>%
  filter((severe_liver_disease!=1))%>%
  dplyr::select(-severe_liver_disease)
cleant2_reduced <- cleant2_reduced %>%
  filter((malignant_cancer!=1))%>%
  dplyr::select(-malignant_cancer)
cleant2_reduced <- cleant2_reduced %>%
  filter((aids!=1))%>%
  dplyr::select(-aids)
cleant2_reduced <- cleant2_reduced %>%
  mutate(
    npar_ratio = neutrophils / albumin,
    npar_ratio = ifelse(albumin == 0, NA, npar_ratio),
    # 按分位数分组
    npar= cut(npar_ratio,
              breaks = quantile(npar_ratio, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE),
              labels = c("T1", "T2", "T3"),
              include.lowest = TRUE)
  )
table(cleant2_reduced$npar)
summary(cleant2_reduced$npar_ratio)
cleant2_reduced <- cleant2_reduced %>%
  filter(!is.na(npar_ratio)) %>%
  dplyr::select(-status)
table(is.na(cleant2_reduced$npar_ratio))
str(cleant2_reduced)
cleant2_reduced <- cleant2_reduced %>%
  mutate(
    across(c(gender, stau_icu, stau_30d, stau_90d, stau_365d, 
             icu_mortality, diabetes, hypertension, npar), 
           factor)
  )
write.csv(cleant2_reduced,file = "cleant2_reduced.csv")
names(cleant2_reduced)
library(epiDisplay)
cleant2_reduced <- cleant2_reduced %>%
  mutate(
    PLT = as.numeric(gsub("[^0-9.]", "", PLT)),
    WBC = as.numeric(gsub("[^0-9.]", "", WBC)),
    Hb = as.numeric(gsub("[^0-9.]", "", Hb)),
    Cr = as.numeric(gsub("[^0-9.]", "", Cr)),
    AST = as.numeric(gsub("[^0-9.]", "", AST)),
    ALT = as.numeric(gsub("[^0-9.]", "", ALT)),
    PT = as.numeric(gsub("[^0-9.]", "", PT)),
    weight= as.numeric(gsub("[^0-9.]", "",weight)),
    neutrophils= as.numeric(gsub("[^0-9.]", "",neutrophils)),
    albumin= as.numeric(gsub("[^0-9.]", "",albumin)),
    height= as.numeric(gsub("[^0-9.]", "",height)),
    LDH= as.numeric(gsub("[^0-9.]", "",LDH)),
    PTT= as.numeric(gsub("[^0-9.]", "",PTT)),
    potassium= as.numeric(gsub("[^0-9.]", "",potassium)),
    Fe= as.numeric(gsub("[^0-9.]", "",Fe)),
    Chol= as.numeric(gsub("[^0-9.]", "",Chol)),
    CA= as.numeric(gsub("[^0-9.]", "",CA)),
    ALT= as.numeric(gsub("[^0-9.]", "",ALT)),
    AST= as.numeric(gsub("[^0-9.]", "",AST))
  )
table_mod1<-tableStack(vars = age:npar_ratio,by=npar,dataFrame =cleant2_reduced) #基本统计
write.csv(table_mod1,file = "table1基本统计描述npar分1.csv")
cleant2_reduced <- cleant2_reduced %>%
  dplyr::select(-albumin,-neutrophils)
missing_summary <- cleant2_reduced %>% 
  summarise(across(everything(), ~ sum(is.na(.x)) / n() * 100)) %>%
  t() %>% 
  as.data.frame() %>% 
  setNames("Missing_Percentage")
cleant2_reduced <- cleant2_reduced %>%
  dplyr::select(-Chol,-Fe,-LDH,-height,-PTT,-potassium)
imp_method <- make.method(cleant2_reduced)
tempData <- mice(cleant2_reduced,m=5,method = "pmm",maxit=50,seed=2025)
summary(tempData)
stripplot(tempData, col=c("grey",mdc(2)),pch=c(1,20))
densityplot(tempData)
names(cleant2_reduced)
fit <- with(tempData, 
            glm(stau_90d ~ npar_ratio+ age + gender
                + sapsii+hypertension
                + gcs_min+diabetes
                + charlson_comorbidity_index+
                  icu_mortality+  Cr+
                  PLT + WBC + Hb +PT+
                  sbp_mean + dbp_mean + CA+ALT+AST, 
                family = binomial))
summary(fit)
summary(fit,type="glance")
xs<-pool(fit)
summary(xs) 
tidy(xs)
completedData_mod <- complete(tempData,action =5)
write.csv(completedData_mod ,file = "completedData_mod.csv")
rm(list = ls()) #清空工作盘
library(survival)
library(survminer)
library(readxl)

completedData_mod <- read.csv("completedData_mod.csv",
                              header = TRUE,sep = ",")

#因子化hypertension diabetes
completedData_mod <- completedData_mod %>%
  mutate(diabetes =factor(diabetes,levels = c(0,1),
                     labels = c("no","yes")),
    hypertension =factor(hypertension,levels = c(0,1),
                         labels = c("no","yes")),
  )

str(completedData_mod)
splots <- list()
fit_Surgery <- survfit(Surv(survival_time,icu_mortality) ~ npar,
                       data =completedData_mod) 
splots[[1]]<- ggsurvplot(fit_Surgery, 
                         data =completedData_mod,
                         mark.time=FALSE,
                         conf.int =FALSE,
                         pval.method=T,
                         pval.method.coord=c(0.3,0.25),pval.method.size =3.3,
                         pval = TRUE,
                         pval.coord=c(0.2,0.2),pval.size =3,
                         risk.table = TRUE,
                         title="ICU survival",
                         xlab = "Follow up time(Days)",
                         xlim=c(0, 60),
                         legend.title="NPAR",
                         legend.labs= c("T1","T2","T3"),
                         censor=T, censor.shape=126, censor.size=2,
                         break.x.by = 5)
ggtheme = theme_bw()#边框+底纹
completedData_mod <- completedData_mod %>%
  mutate(stau_30d = as.integer(stau_30d == "yes")) %>%
  mutate(stau_90d = as.integer(stau_90d == "yes")) %>%
  mutate(stau_365d = as.integer(stau_365d == "yes"))
#30天死亡
fit_Surgery3 <- survfit(Surv(survival_time,stau_30d) ~ npar,
                        data =completedData_mod) 
splots[[2]]<- ggsurvplot(fit_Surgery3, 
                         data =completedData_mod,
                         mark.time=FALSE,
                         conf.int =FALSE,
                         pval.method=T,
                         pval.method.coord=c(0.3,0.25),pval.method.size =3.3,
                         pval = TRUE,
                         pval.coord=c(0.2,0.2),pval.size =3,
                         risk.table = TRUE,
                         title="30d survival",
                         xlab = "Follow up time(Days)",
                         xlim=c(0, 30),
                         legend.title="NPAR",
                         legend.labs= c("T1","T2","T3"),
                         censor=T, censor.shape=126, censor.size=2,
                         break.x.by =3)

#90天死亡
fit_Surgery4 <- survfit(Surv(survival_time,stau_90d) ~ npar,
                        data =completedData_mod) 
splots[[3]]<- ggsurvplot(fit_Surgery4, 
                         data =completedData_mod,
                         mark.time=FALSE,
                         conf.int =FALSE,
                         pval.method=T,
                         pval.method.coord=c(0.3,0.25),pval.method.size =3.3,
                         pval = TRUE,
                         pval.coord=c(0.2,0.2),pval.size =3,
                         risk.table = TRUE,
                         title="90d survival",
                         xlab = "Follow up time(Days)",
                         xlim=c(0, 90),
                         legend.title="NPAR",
                         legend.labs= c("T1","T2","T3"),
                         censor=T, censor.shape=126, censor.size=2,
                         break.x.by = 10)
fit_Surgery5 <- survfit(Surv(survival_time,stau_365d) ~ npar,
                        data =completedData_mod) 
splots[[4]]<- ggsurvplot(fit_Surgery5, 
                         data =completedData_mod,
                         mark.time=FALSE,
                         conf.int =FALSE,
                         pval.method=T,
                         pval.method.coord=c(0.3,0.25),pval.method.size =3.3,
                         pval = TRUE,
                         pval.coord=c(0.2,0.2),pval.size =3,
                         risk.table = TRUE,
                         title="365d survival",
                         xlab = "Follow up time(Days)",
                         xlim=c(0, 365),
                         legend.title="NPAR",
                         legend.labs= c("T1","T2","T3"),
                         censor=T, censor.shape=126, censor.size=2,
                         break.x.by = 30)
res <- arrange_ggsurvplots(splots, print = TRUE,
                           ncol = 2, nrow = 2,
                           risk.table.height = 0.25)
ggsave("K-M_mod.pdf", dpi = 600,width =14, height = 12,res)
ggsave("K-M_mod.png", dpi = 600,width =14, height = 12,res)

library(dplyr)
library(tidyverse)
library(survey)
library(tableone)
library(reshape2)
rm(list = ls()) #清空工作盘
bc<- read.csv("completedData_mod.csv", header=TRUE, sep=",")
numeric_cols_bc <- sapply(bc, is.numeric)
numeric_vars <- names(numeric_cols_bc)[numeric_cols_bc]
shapiro_test_results <- list()
for (col in numeric_vars) {
  shapiro_test_result <- shapiro.test(bc[[col]])
  shapiro_test_results[[col]] <- shapiro_test_result
}
print(shapiro_test_results)
str(bc)
df <-bc%>%
  mutate(gender=factor(gender, labels = c("F" ,"M")), 
         npar=factor(npar,labels=c("T1","T2","T3")),
         stau_icu=factor(stau_icu,levels = c(0, 1),
                       labels = c("survive","dead")),
         stau_30d=factor(stau_30d,levels = c("no", "yes"),
                       labels = c("survive","dead")),
         stau_90d=factor(stau_90d,levels = c("no", "yes"),
                       labels = c("survive","dead")),
         stau_365d=factor(stau_365d,levels = c("no", "yes"),
                         labels = c("survive","dead")),
         icu_mortality=factor(icu_mortality,levels = c(0,1),
                                labels = c("survive","dead")),
         diabetes=factor(diabetes,levels = c(0,1),
                              labels = c("no","yes")),
         hypertension=factor(hypertension,levels = c(0,1),
                              labels = c("no","yes")),
         )
y1 <- df
str(df)
glimpse(df)
names(df)
df_stats <- df %>%
  dplyr::select(npar, AST, ALT, PT, Cr) %>%    
  pivot_longer(cols = -npar,              
               names_to = "variable",
               values_to = "value") %>%
  group_by(npar, variable) %>%           
  summarise(
    median = round(median(value, na.rm = TRUE), 1), 
    Q1 = round(quantile(value, 0.25, na.rm = TRUE), 1), 
    Q3 = round(quantile(value, 0.75, na.rm = TRUE), 1), 
    .groups = "drop"
  ) %>%
  mutate(
    median_iqr = sprintf("%.2f (%.2f-%.2f)", median, Q1, Q3) 
  ) %>%
  dplyr::select(-median, -Q1, -Q3) %>% 
  pivot_wider(names_from = npar, 
              values_from = median_iqr)
print(df_stats)
nonnormal_vars <- c("AST", "ALT", "Cr", "PT")
library(broom)
p_values <- df %>%
  dplyr::select(npar, all_of(nonnormal_vars)) %>%
  pivot_longer(cols = -npar, names_to = "variable", values_to = "value") %>%
  nest(data = -variable) %>%
  mutate(
    test = map(data, ~ {
      if(length(unique(.x$npar)) > 2) {
        tidy(kruskal.test(value ~ npar, data = .x))
      } else {
        tidy(wilcox.test(value ~ npar, data = .x))
      }
    })
  ) %>%
  unnest(test) %>%
  dplyr::select(variable, p.value) %>%
  mutate(
    p.value = format.pval(p.value, 
                          digits = 2, 
                          eps = 0.001,
                          na.form = "NA")
  )
final_table <- df_stats %>%
  left_join(p_values, by = "variable") %>%
  relocate(p.value, .after = last_col())
print(final_table)
covariates <- c("age", "gender", 
  "sapsii", "charlson_comorbidity_index","diabetes",
  "sbp_mean", "dbp_mean", "gcs_min",
  "PLT", "WBC", "weight",
  "Hb", "Cr", "PT",
  "CA", "ALT","hypertension","npar_ratio","AST"
)
factorVars <- c("gender", "diabetes", "hypertension", 
                "stau_30d", "stau_90d", "stau_365d", "icu_mortality")
tab_Unmatched <- CreateTableOne(vars = covariates, 
                                strata = "npar", 
                                data = df, 
                                factorVars = factorVars,
                                test =T)
print(tab_Unmatched,showAllLevels=TRUE,smd=TRUE)
library(ipw)
w1 <- ipwpoint(
  exposure = npar,
  family = "multinomial",
  numerator = ~ 1,
  denominator = ~age+gender
  +sbp_mean+dbp_mean+gcs_min+PT+weight
  +sapsii+charlson_comorbidity_index
  +PLT+WBC+Hb+Cr+CA+AST+ALT,
  data = df)

#将计算的权重带入数据
df$w1<-w1$ipw.weights
str(df)

#提取IPTW后的数据
dataIPTW<-svydesign(ids = ~1,strata = ~npar,weights = ~w1,
                    nest = TRUE,data=df)
summary(dataIPTW)
factorVars <- c("gender","diabetes","hypertension")
tab_IPTW<-svyCreateTableOne(vars=covariates, 
                            strata="npar",
                            factorVars=factorVars,
                            data=dataIPTW
                            )
print(tab_IPTW,showAllLevels=TRUE,smd=TRUE)

#4\标准化平均差SMD可视化
#提取作图数据
dataPlot <- data.frame(variable=rownames(ExtractSmd(tab_Unmatched)),
                       Unmatched=as.numeric(ExtractSmd(tab_Unmatched)),
                       IPTW=as.numeric(ExtractSmd(tab_IPTW))
)
#指定将要出现在图中的变量

library(reshape2)
dataPlotMelt<-melt(data= dataPlot,
                   id.vars=c("variable"),
                   variable.name= "Method",
                   value.name= "SMD")

#
varNames <- as.character(dataPlot$variable)[order(dataPlot$Unmatched)]
varNames <- varNames[!duplicated(varNames)]  # 去重但保留顺序
dataPlotMelt$variable <- factor(dataPlotMelt$variable, levels = varNames)
#画图
ggplot(data = dataPlotMelt,
       mapping = aes(x = variable, y = SMD, 
                     group = Method, 
                     color = Method,
                     shape = Method )) +
  #geom_line() +
  geom_point(size=4) +
  geom_hline(yintercept = 0.1, 
             color = "red",
             lty=2,
             size = 0.1) +
  coord_flip() +
  theme_bw(base_size = 18)

#5\两个Table 1合并输出
#1.提取两个结果
table1<- cbind(print(tab_Unmatched,printToggle =F,showAllLevels=T,),
               print(tab_IPTW,printToggle =F,showAllLevels=T,)
)

# 插入一行分组
table1<- rbind(Group=rep(c("Level","T1","T2","T3","P","test"),2),table1)

#更改列名
colnames(table1) <- c("Level","Unmatched",NA,NA,NA,NA,
                      "Level","IPTW",NA,NA,NA,NA)

#打印或导出Excel
print(table1, quote = FALSE)
write.csv(table1, file = "tablemod1IPTW.csv")

library(survival)
library(survminer)
library(readxl)
library(broom)
library(dplyr)
library(ggplot2)
library(dplyr)
y1 <- y1 %>%
  mutate(
    across(c(stau_icu,stau_30d, stau_90d, stau_365d,icu_mortality), 
           ~ case_when(.x == "dead" ~ 1, .x == "survive" ~ 0))
  )

outcomes_ordered <- c("icu_mortality", 
                      "stau_30d", "stau_90d", "stau_365d")
predictors <- c("npar_ratio", "npar")

results_df <- data.frame()

for (outcome in outcomes_ordered) { 
  time_var <- if(outcome == "icu_mortality") {
    "time_icuD"
  } else {
    gsub("stau", "time", outcome)  # stau_30d → time_30d
  }
  for (pred in predictors) {
    surv.obj <- Surv(time = y1[[time_var]],
                     event = y1[[outcome]])
    COX1 <- coxph(as.formula(paste("surv.obj ~", pred)), 
                  data = y1)
    
    tidy_res <- tidy(COX1, conf.int = TRUE) %>%
      mutate(
        Outcome = factor(outcome, levels = outcomes_ordered), 
        Predictor = pred,
        HR_CI = sprintf("%.2f (%.2f-%.2f)", 
                        exp(estimate), exp(conf.low), exp(conf.high)) 
      )
    
    results_df <- rbind(results_df, tidy_res)
  }
}
print(results_df)



#有多个变量时可进行p值矫正
# 提取原始p值（排除NA值）
raw_p <- summary(COX1)$coefficients[, "Pr(>|z|)"]
raw_p <- na.omit(raw_p)  # 删除缺失值
# Bonferroni校正（控制FWER）
bonferroni_p <- p.adjust(raw_p, method = "bonferroni") 
# FDR校正（Benjamini-Hochberg）
fdr_p <- p.adjust(raw_p, method = "fdr")
# 构建结果表
results1_pjust <- data.frame(
  Variable = names(raw_p),
  Raw_p = round(raw_p, 3),
  Bonferroni_p = round(bonferroni_p, 4),
  FDR_p = round(fdr_p, 3),
  Significance = ifelse(fdr_p < 0.05, "*", "")
)
print(results1_pjust)
event_vars <- c("icu_mortality","stau_30d", 
                "stau_90d", "stau_365d")
predictors <- c("npar_ratio", "npar")
analy <- y1
for (event_var in event_vars) {
  time_var <- if(event_var == "icu_mortality") {
    "time_icuD"
  } else {gsub("stau", "time", event_var)
  }
  surv.obj <- Surv(
    time = analy[[time_var]], 
    event = analy[[event_var]]
  )
  for (pred in predictors) {
    formula_str <- paste("surv.obj ~", pred, "+ age + gender")
    cox_model <- coxph(as.formula(formula_str), data = analy)
    cat("\n", strrep("-", 40), 
        "\n", toupper(event_var), "with predictor:", pred,
        "\n", strrep("-", 40), "\n")
    print( cbind(
      summary(cox_model)$coefficients[, c("coef", "exp(coef)", "Pr(>|z|)")], 
      "95% CI Lower" = summary(cox_model)$conf.int[, "lower .95"],
      "95% CI Upper" = summary(cox_model)$conf.int[, "upper .95"]
    ) )
  }
}
#Model3(全部调整)
event_vars <- c("icu_mortality", 
                "stau_30d", "stau_90d", "stau_365d")
predictors <- c("npar_ratio", "npar") 
covariates <- c("age", "gender", "sapsii", "hypertension", "weight" ,
                "sbp_mean", "dbp_mean", "gcs_min", "diabetes",
                "charlson_comorbidity_index", "PLT", "WBC", 
                "Hb", "CA", "Cr","PT","ALT","AST")
analysis<- y1
for (event_var in event_vars) { 
  time_var <- if(event_var == "icu_mortality") {
  "time_icuD"
      } else {gsub("stau", "time", event_var)
}
  surv.obj <- Surv(
    time = analysis[[time_var]], 
    event = analysis[[event_var]]
  )
  for (pred in predictors) {
    formula_str <- paste("surv.obj ~", pred, "+", paste(covariates, collapse = " + "))
    cox_model <- coxph(as.formula(formula_str), data = analysis)
    res <- summary(cox_model)
    coef_table <- res$coefficients
    confint_table <- res$conf.int
    cat("\n", strrep("=", 60), 
        "\n 模型: ", toupper(event_var), 
        "\n 主变量: ", pred,
        "\n", strrep("-", 60),
        "\n 所有变量核心指标:\n")
    for(var in rownames(coef_table)){
      hr <- exp(coef_table[var, "coef"])
      ci_low <- confint_table[var, "lower .95"]
      ci_high <- confint_table[var, "upper .95"]
      p_value <- format.pval(coef_table[var, "Pr(>|z|)"], eps = 0.001)
      
      cat(sprintf(" - %-20s HR=%.2f (95%% CI %.2f-%.2f), p=%s\n",
                  var, hr, ci_low, ci_high, p_value))
    }
    cat(strrep("=", 60), "\n")
  }
}

library(jstable)
library(survival)
library(epiDisplay)
library(tidyverse)
npar_breaks <- quantile(y1$npar_ratio, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
print(npar_breaks)
df <-y1%>%
  mutate(gender=factor(gender,labels=c("F","M")),
         age=ifelse(age >74,">74","<=74"),
         age=factor(age, levels=c(">74","<=74")),
         sapsii=ifelse(sapsii >37,">37","<=37"),
         sapsii=factor(sapsii, levels=c(">37","<=37")),
         gcs_min=ifelse(gcs_min >14,">14","<=14"),
         gcs_min=factor(gcs_min, levels=c(">14","<=14")),
         charlson_comorbidity_index= ifelse(charlson_comorbidity_index >6,">6","<=6"),
         charlson_comorbidity_index=factor(charlson_comorbidity_index, levels=c(">6","<=6")))
str(df)
library(dplyr)
interaction_results <- list()
# age分析ICU死亡 --------------------------------------------------------------
cox_icu1 <- coxph(Surv(time_icuD,icu_mortality) 
                  ~ npar + age + sbp_mean + dbp_mean + gender+
                    hypertension+diabetes+charlson_comorbidity_index+
                    gcs_min+sapsii+weight+
                    PT + PLT + WBC + Hb + Cr +CA+ AST + ALT,
                  data = df)
print(summary(cox_icu1))
cox_icu2 <- coxph(Surv(time_icuD, icu_mortality) 
                  ~ npar + age + npar*age + sbp_mean + dbp_mean + weight+
                    PT + PLT + WBC + Hb + Cr + AST + ALT+CA,
                  data = df)
print(summary(cox_icu2))
print(anova(cox_icu1, cox_icu2, test = "Chisq"))
interaction_results[["age_icuD"]] <- anova(cox_icu1, cox_icu2, test = "Chisq")$P[2]
## 分析30天死亡 -------------------------------------------------------------
cox_30d1 <- coxph(Surv(time_30d, stau_30d) 
                  ~ npar + age + sbp_mean + dbp_mean + gender+
                    hypertension+diabetes+charlson_comorbidity_index+
                    gcs_min+sapsii+weight+
                    PT + PLT + WBC + Hb + Cr + AST + ALT+CA,
                  data = df)
print(summary(cox_30d1))
cox_30d2 <- coxph(Surv(time_30d, stau_30d) 
                  ~ npar + age + npar*age + sbp_mean + dbp_mean + gender+
                    hypertension+diabetes+charlson_comorbidity_index+
                    gcs_min+sapsii+weight+
                    PT + PLT + WBC + Hb + Cr + AST + ALT+CA,
                  data = df)
print(summary(cox_30d2))
print(anova(cox_30d1, cox_30d2, test = "Chisq"))
interaction_results[["age_30d"]] <- anova(cox_30d1, cox_30d2, test = "Chisq")$P[2]
## 分析90天死亡 -------------------------------------------------------------
cox_90d1 <- coxph(Surv(time_90d, stau_90d) 
                  ~ npar + age + sbp_mean + dbp_mean + gender+
                    hypertension+diabetes+charlson_comorbidity_index+
                    gcs_min+sapsii+weight+
                    PT + PLT + WBC + Hb + Cr + AST + ALT+CA,
                  data = df)
print(summary(cox_90d1))
cox_90d2 <- coxph(Surv(time_90d, stau_90d) 
                  ~ npar + age + npar*age + sbp_mean + dbp_mean + gender+
                    hypertension+diabetes+charlson_comorbidity_index+
                    gcs_min+sapsii+weight+
                    PT + PLT + WBC + Hb + Cr + AST + ALT+CA,
                  data = df)
print(summary(cox_90d2))
print(anova(cox_90d1, cox_90d2, test = "Chisq"))
interaction_results[["age_90d"]] <- anova(cox_90d1, cox_90d2, test = "Chisq")$P[2]
## 分析365天死亡 ------------------------------------------------------------
cox_365d1 <- coxph(Surv(time_365d, stau_365d) 
                   ~ npar + age + sbp_mean + dbp_mean + gender+
                     hypertension+diabetes+charlson_comorbidity_index+
                     gcs_min+sapsii+weight+
                     PT + PLT + WBC + Hb + Cr + AST + ALT+CA,
                   data = df)
print(summary(cox_365d1))
cox_365d2 <- coxph(Surv(time_365d, stau_365d) 
                   ~ npar + age + npar*age + sbp_mean + dbp_mean + gender+
                     hypertension+diabetes+charlson_comorbidity_index+
                     gcs_min+sapsii+weight+
                     PT + PLT + WBC + Hb + Cr + AST + ALT+CA,
                   data = df)
print(summary(cox_365d2))
print(anova(cox_365d1, cox_365d2, test = "Chisq"))
interaction_results[["age_365d"]] <- anova(cox_365d1, cox_365d2, test = "Chisq")$P[2]
# gender分析ICU死亡 ------------------------------------------------------------
cox_icu1 <- coxph(Surv(time_icuD, icu_mortality) 
                  ~ npar + gender + sbp_mean + dbp_mean + age+
                    hypertension+diabetes+charlson_comorbidity_index+
                    gcs_min+sapsii+weight+
                    PT + PLT + WBC + Hb + Cr + AST + ALT+CA,
                  data = df)
print(summary(cox_icu1))
cox_icu2 <- coxph(Surv(time_icuD, icu_mortality) 
                  ~ npar + gender + npar*gender + sbp_mean + dbp_mean + age+
                    hypertension+diabetes+charlson_comorbidity_index+
                    gcs_min+sapsii+weight+
                    PT + PLT + WBC + Hb + Cr + AST + ALT+CA,
                  data = df)
print(summary(cox_icu2))
print(anova(cox_icu1, cox_icu2, test = "Chisq"))
interaction_results[["gender_icuD"]] <- anova(cox_icu1, cox_icu2, test = "Chisq")$P[2]
## 分析30天死亡 -------------------------------------------------------------
cox_30d1 <- coxph(Surv(time_30d, stau_30d) 
                  ~ npar + gender + sbp_mean + dbp_mean + age+
                    hypertension+diabetes+charlson_comorbidity_index+
                    gcs_min+sapsii+weight+
                    PT + PLT + WBC + Hb + Cr + AST + ALT+CA,
                  data = df)
print(summary(cox_30d1))
cox_30d2 <- coxph(Surv(time_30d, stau_30d) 
                  ~ npar + gender + npar*gender + sbp_mean + dbp_mean + age+
                    hypertension+diabetes+charlson_comorbidity_index+
                    gcs_min+sapsii+weight+
                    PT + PLT + WBC + Hb + Cr + AST + ALT+CA,
                  data = df)
print(summary(cox_30d2))
print(anova(cox_30d1, cox_30d2, test = "Chisq"))
interaction_results[["gender_30d"]] <- anova(cox_30d1, cox_30d2, test = "Chisq")$P[2]
## 分析90天死亡 -------------------------------------------------------------
cox_90d1 <- coxph(Surv(time_90d, stau_90d) 
                  ~ npar + gender + sbp_mean + dbp_mean + age+
                    hypertension+diabetes+charlson_comorbidity_index+
                    gcs_min+sapsii+weight+
                    PT + PLT + WBC + Hb + Cr + AST + ALT+CA,
                  data = df)
print(summary(cox_90d1))
cox_90d2 <- coxph(Surv(time_90d, stau_90d) 
                  ~ npar + gender + npar*gender + sbp_mean + dbp_mean +age+
                    hypertension+diabetes+charlson_comorbidity_index+
                    gcs_min+sapsii+ weight+
                    PT + PLT + WBC + Hb + Cr + AST + ALT+CA,
                  data = df)
print(summary(cox_90d2))
print(anova(cox_90d1, cox_90d2, test = "Chisq"))
interaction_results[["gender_90d"]] <- anova(cox_90d1, cox_90d2, test = "Chisq")$P[2]
## 分析365天死亡 ------------------------------------------------------------
cox_365d1 <- coxph(Surv(time_365d, stau_365d) 
                   ~ npar + gender + sbp_mean + dbp_mean + age+
                     hypertension+diabetes+charlson_comorbidity_index+
                     gcs_min+sapsii+weight+
                     PT + PLT + WBC + Hb + Cr + AST + ALT+CA,
                   data = df)
print(summary(cox_365d1))
cox_365d2 <- coxph(Surv(time_365d, stau_365d) 
                   ~ npar + gender + npar*gender + sbp_mean + dbp_mean + age+
                     hypertension+diabetes+charlson_comorbidity_index+
                     gcs_min+sapsii+weight+
                     PT + PLT + WBC + Hb + Cr + AST + ALT+CA,
                   data = df)
print(summary(cox_365d2))
print(anova(cox_365d1, cox_365d2, test = "Chisq"))
interaction_results[["gender_365d"]] <- anova(cox_365d1, cox_365d2, test = "Chisq")$P[2]
#npar*hypertension
# hypertension分析ICU死亡 --------------------------------------------------------------
cox_icu1 <- coxph(Surv(time_icuD, icu_mortality) 
                  ~ npar + hypertension + sbp_mean + dbp_mean + age+
                    gender+diabetes+charlson_comorbidity_index+
                    gcs_min+sapsii+weight+
                    PT + PLT + WBC + Hb + Cr + AST + ALT+CA,
                  data = df)
print(summary(cox_icu1))
cox_icu2 <- coxph(Surv(time_icuD, icu_mortality) 
                  ~ npar + hypertension + npar*hypertension + sbp_mean + dbp_mean + 
                    age+gender+diabetes+charlson_comorbidity_index+weight+
                    gcs_min+sapsii+ PT + PLT + WBC + Hb + Cr + AST + ALT+CA,
                  data = df)
print(summary(cox_icu2))
print(anova(cox_icu1, cox_icu2, test = "Chisq"))
interaction_results[["hypertension_icuD"]] <- anova(cox_icu1, cox_icu2, test = "Chisq")$P[2]
## 分析30天死亡 -------------------------------------------------------------
cox_30d1 <- coxph(Surv(time_30d, stau_30d) 
                  ~ npar + hypertension + sbp_mean + dbp_mean + weight+
                    age+gender+diabetes+charlson_comorbidity_index+
                    gcs_min+sapsii+PT + PLT + WBC + Hb + Cr + AST + ALT+CA,
                  data = df)
print(summary(cox_30d1))
cox_30d2 <- coxph(Surv(time_30d, stau_30d) 
                  ~ npar + hypertension + npar*hypertension + sbp_mean + dbp_mean + 
                    age+weight+
                    gender+diabetes+charlson_comorbidity_index+
                    gcs_min+sapsii+PT + PLT + WBC + Hb + Cr + AST + ALT+CA,
                  data = df)
print(summary(cox_30d2))
print(anova(cox_30d1, cox_30d2, test = "Chisq"))
interaction_results[["hypertension_30d"]] <- anova(cox_30d1, cox_30d2, test = "Chisq")$P[2]
## 分析90天死亡 -------------------------------------------------------------
cox_90d1 <- coxph(Surv(time_90d, stau_90d) 
                  ~ npar + hypertension + sbp_mean + dbp_mean + 
                    age+weight+
                    gender+diabetes+charlson_comorbidity_index+
                    gcs_min+sapsii+PT + PLT + WBC + Hb + Cr + AST + ALT+CA,
                  data = df)
print(summary(cox_90d1))
cox_90d2 <- coxph(Surv(time_90d, stau_90d) 
                  ~ npar + hypertension + npar*hypertension + sbp_mean + dbp_mean + 
                    age+weight+
                    gender+diabetes+charlson_comorbidity_index+
                    gcs_min+sapsii+PT + PLT + WBC + Hb + Cr + AST + ALT+CA,
                  data = df)
print(summary(cox_90d2))
print(anova(cox_90d1, cox_90d2, test = "Chisq"))
interaction_results[["hypertension_90d"]] <- anova(cox_90d1, cox_90d2, test = "Chisq")$P[2]
## 分析365天死亡 ------------------------------------------------------------
cox_365d1 <- coxph(Surv(time_365d, stau_365d) 
                   ~ npar + hypertension + sbp_mean + dbp_mean + 
                     age+weight+
                     gender+diabetes+charlson_comorbidity_index+
                     gcs_min+sapsii+PT + PLT + WBC + Hb + Cr + AST + ALT+CA,
                   data = df)
print(summary(cox_365d1))
cox_365d2 <- coxph(Surv(time_365d, stau_365d) 
                   ~ npar + hypertension + npar*hypertension + sbp_mean + dbp_mean + 
                     age+weight+
                     gender+diabetes+charlson_comorbidity_index+
                     gcs_min+sapsii+PT + PLT + WBC + Hb + Cr + AST + ALT+CA,
                   data = df)
print(summary(cox_365d2))
print(anova(cox_365d1, cox_365d2, test = "Chisq"))
interaction_results[["hypertension_365d"]] <- anova(cox_365d1, cox_365d2, test = "Chisq")$P[2]
# diabetes分析ICU死亡
cox_icu1 <- coxph(Surv(time_icuD, icu_mortality) 
                  ~ npar + diabetes + sbp_mean + dbp_mean + age+
                    gender+hypertension+charlson_comorbidity_index+
                    gcs_min+sapsii+weight+
                    PT + PLT + WBC + Hb + Cr + AST + ALT+CA,
                  data = df)
print(summary(cox_icu1))
cox_icu2 <- coxph(Surv(time_icuD, icu_mortality) 
                  ~ npar + diabetes + npar*diabetes + sbp_mean + dbp_mean + 
                    age+weight+
                    gender+hypertension+charlson_comorbidity_index+
                    gcs_min+sapsii+PT + PLT + WBC + Hb + Cr + AST + ALT+CA,
                  data = df)
print(summary(cox_icu2))
print(anova(cox_icu1, cox_icu2, test = "Chisq"))
interaction_results[["diabetes_icuD"]] <- anova(cox_icu1, cox_icu2, test = "Chisq")$P[2]
## 分析30天死亡
cox_30d1 <- coxph(Surv(time_30d, stau_30d) 
                  ~ npar + diabetes + sbp_mean + dbp_mean + 
                    age+weight+
                    gender+hypertension+charlson_comorbidity_index+
                    gcs_min+sapsii+ PT + PLT + WBC + Hb + Cr + AST + ALT+CA,
                  data = df)
print(summary(cox_30d1))
cox_30d2 <- coxph(Surv(time_30d, stau_30d) 
                  ~ npar + diabetes + npar*diabetes + sbp_mean + dbp_mean + 
                    age+weight+
                    gender+hypertension+charlson_comorbidity_index+
                    gcs_min+sapsii+PT + PLT + WBC + Hb + Cr + AST + ALT+CA,
                  data = df)
print(summary(cox_30d2))
print(anova(cox_30d1, cox_30d2, test = "Chisq"))
interaction_results[["diabetes_30d"]] <- anova(cox_30d1, cox_30d2, test = "Chisq")$P[2]
## 分析90天死亡 -------------------------------------------------------------
cox_90d1 <- coxph(Surv(time_90d, stau_90d) 
                  ~ npar + diabetes + sbp_mean + dbp_mean + 
                    age+weight+
                    gender+hypertension+charlson_comorbidity_index+
                    gcs_min+sapsii+PT + PLT + WBC + Hb + Cr + AST + ALT+CA,
                  data = df)
print(summary(cox_90d1))
cox_90d2 <- coxph(Surv(time_90d, stau_90d) 
                  ~ npar + diabetes + npar*diabetes + sbp_mean + dbp_mean + 
                    age+weight+
                    gender+hypertension+charlson_comorbidity_index+
                    gcs_min+sapsii+PT + PLT + WBC + Hb + Cr + AST + ALT+CA,
                  data = df)
print(summary(cox_90d2))
print(anova(cox_90d1, cox_90d2, test = "Chisq"))
interaction_results[["diabetes_90d"]] <- anova(cox_90d1, cox_90d2, test = "Chisq")$P[2]
## 分析365天死亡 ------------------------------------------------------------
cox_365d1 <- coxph(Surv(time_365d, stau_365d) 
                   ~ npar + diabetes + sbp_mean + dbp_mean + 
                     age+weight+
                     gender+hypertension+charlson_comorbidity_index+
                     gcs_min+sapsii+ PT + PLT + WBC + Hb + Cr + AST + ALT+CA,
                   data = df)
print(summary(cox_365d1))
cox_365d2 <- coxph(Surv(time_365d, stau_365d) 
                   ~ npar + diabetes + npar*diabetes + sbp_mean + dbp_mean + 
                     age+weight+
                     gender+hypertension+charlson_comorbidity_index+
                     gcs_min+sapsii+PT + PLT + WBC + Hb + Cr + AST + ALT+CA,
                   data = df)
print(summary(cox_365d2))
print(anova(cox_365d1, cox_365d2, test = "Chisq"))
interaction_results[["diabetes_365d"]] <- anova(cox_365d1, cox_365d2, test = "Chisq")$P[2]
#  npar*sapsii
cox_icu1 <- coxph(Surv(time_icuD,icu_mortality) 
                  ~ npar + sapsii + sbp_mean + dbp_mean + 
                    age+weight+
                    gender+hypertension+charlson_comorbidity_index+
                    gcs_min+diabetes+ PT + PLT + WBC + Hb + Cr + AST + ALT+CA,
                  data = df)
print(summary(cox_icu1))
cox_icu2 <- coxph(Surv(time_icuD, icu_mortality) 
                  ~ npar + sapsii + npar*sapsii + sbp_mean + dbp_mean + 
                    age+weight+
                    gender+hypertension+charlson_comorbidity_index+
                    gcs_min+diabetes+PT + PLT + WBC + Hb + Cr + AST + ALT+CA,
                  data = df)
print(summary(cox_icu2))
print(anova(cox_icu1, cox_icu2, test = "Chisq"))
interaction_results[["sapsii_icuD"]] <- anova(cox_icu1, cox_icu2, test = "Chisq")$P[2]
# sapsii分析30天死亡 -------------------------------------------------------------
cox_30d1 <- coxph(Surv(time_30d, stau_30d) 
                  ~ npar + sapsii + sbp_mean + dbp_mean + 
                    age+weight+
                    gender+hypertension+charlson_comorbidity_index+
                    gcs_min+diabetes+PT + PLT + WBC + Hb + Cr + AST + ALT+CA,
                  data = df)
print(summary(cox_30d1))
cox_30d2 <- coxph(Surv(time_30d, stau_30d) 
                  ~ npar + sapsii + npar*sapsii + sbp_mean + dbp_mean + 
                    age+weight+
                    gender+hypertension+charlson_comorbidity_index+
                    gcs_min+diabetes+PT + PLT + WBC + Hb + Cr + AST + ALT+CA,
                  data = df)
print(summary(cox_30d2))
print(anova(cox_30d1, cox_30d2, test = "Chisq"))
interaction_results[["sapsii_30d"]] <- anova(cox_30d1, cox_30d2, test = "Chisq")$P[2]
## 分析90天死亡 -------------------------------------------------------------
cox_90d1 <- coxph(Surv(time_90d, stau_90d) 
                  ~ npar + sapsii + sbp_mean + dbp_mean + 
                    age+weight+
                    gender+hypertension+charlson_comorbidity_index+
                    gcs_min+diabetes+PT + PLT + WBC + Hb + Cr + AST + ALT+CA,
                  data = df)
print(summary(cox_90d1))
cox_90d2 <- coxph(Surv(time_90d, stau_90d) 
                  ~ npar + sapsii + npar*sapsii + sbp_mean + dbp_mean + 
                    age+weight+
                    gender+hypertension+charlson_comorbidity_index+
                    gcs_min+diabetes+PT + PLT + WBC + Hb + Cr + AST + ALT+CA,
                  data = df)
print(summary(cox_90d2))
print(anova(cox_90d1, cox_90d2, test = "Chisq"))
interaction_results[["sapsii_90d"]] <- anova(cox_90d1, cox_90d2, test = "Chisq")$P[2]
## 分析365天死亡 ------------------------------------------------------------
cox_365d1 <- coxph(Surv(time_365d, stau_365d) 
                   ~ npar + sapsii + sbp_mean + dbp_mean + 
                     age+weight+
                     gender+hypertension+charlson_comorbidity_index+
                     gcs_min+diabetes+PT + PLT + WBC + Hb + Cr + AST + ALT+CA,
                   data = df)
print(summary(cox_365d1))
cox_365d2 <- coxph(Surv(time_365d, stau_365d) 
                   ~ npar + sapsii + npar*sapsii + sbp_mean + dbp_mean + 
                     age+weight+
                     gender+hypertension+charlson_comorbidity_index+
                     gcs_min+diabetes+ PT + PLT + WBC + Hb + Cr + AST + ALT+CA,
                   data = df)
print(summary(cox_365d2))
print(anova(cox_365d1, cox_365d2, test = "Chisq"))
interaction_results[["sapsii_365d"]] <- anova(cox_365d1, cox_365d2, test = "Chisq")$P[2]
#npar*gcs_min
cox_icu1 <- coxph(Surv(time_icuD, icu_mortality) 
                  ~ npar + gcs_min + sbp_mean + dbp_mean + age+
                    gender+hypertension+charlson_comorbidity_index+
                    sapsii+diabetes+weight+
                    PT + PLT + WBC + Hb + Cr + AST + ALT+CA,
                  data = df)
print(summary(cox_icu1))
cox_icu2 <- coxph(Surv(time_icuD, icu_mortality) 
                  ~ npar + gcs_min + npar*gcs_min + sbp_mean + dbp_mean + age+
                    gender+hypertension+charlson_comorbidity_index+
                    sapsii+diabetes+weight+
                    PT + PLT + WBC + Hb + Cr + AST + ALT+CA,
                  data = df)
print(summary(cox_icu2))
print(anova(cox_icu1, cox_icu2, test = "Chisq"))
interaction_results[["gcs_min_icuD"]] <- anova(cox_icu1, cox_icu2, test = "Chisq")$P[2]
#gcs_min 分析30天死亡 -------------------------------------------------------------
cox_30d1 <- coxph(Surv(time_30d, stau_30d) 
                  ~ npar + gcs_min + sbp_mean + dbp_mean + age+
                    gender+hypertension+charlson_comorbidity_index+
                    sapsii+diabetes+weight+
                    PT + PLT + WBC + Hb + Cr + AST + ALT+CA,
                  data = df)
print(summary(cox_30d1))
cox_30d2 <- coxph(Surv(time_30d, stau_30d) 
                  ~ npar + gcs_min + npar*gcs_min + sbp_mean + dbp_mean + 
                    age+weight+
                    gender+hypertension+charlson_comorbidity_index+
                    sapsii+diabetes+PT + PLT + WBC + Hb + Cr + AST + ALT+CA,
                  data = df)
print(summary(cox_30d2))
print(anova(cox_30d1, cox_30d2, test = "Chisq"))
interaction_results[["gcs_min_30d"]] <- anova(cox_30d1, cox_30d2, test = "Chisq")$P[2]
## 分析90天死亡 -------------------------------------------------------------
cox_90d1 <- coxph(Surv(time_90d, stau_90d) 
                  ~ npar + gcs_min + sbp_mean + dbp_mean + 
                    age+weight+
                    gender+hypertension+charlson_comorbidity_index+
                    sapsii+diabetes+PT + PLT + WBC + Hb + Cr + AST + ALT+CA,
                  data = df)
print(summary(cox_90d1))
cox_90d2 <- coxph(Surv(time_90d, stau_90d) 
                  ~ npar + gcs_min + npar*gcs_min + sbp_mean + dbp_mean + 
                    age+weight+
                    gender+hypertension+charlson_comorbidity_index+
                    sapsii+diabetes+PT + PLT + WBC + Hb + Cr + AST + ALT+CA,
                  data = df)
print(summary(cox_90d2))
print(anova(cox_90d1, cox_90d2, test = "Chisq"))
interaction_results[["gcs_min_90d"]] <- anova(cox_90d1, cox_90d2, test = "Chisq")$P[2]
## 分析365天死亡 ------------------------------------------------------------
cox_365d1 <- coxph(Surv(time_365d, stau_365d) 
                   ~ npar + gcs_min + sbp_mean + dbp_mean + 
                     age+weight+
                     gender+hypertension+charlson_comorbidity_index+
                     sapsii+diabetes+ PT + PLT + WBC + Hb + Cr + AST + ALT+CA,
                   data = df)
print(summary(cox_365d1))
cox_365d2 <- coxph(Surv(time_365d, stau_365d) 
                   ~ npar + gcs_min + npar*gcs_min + sbp_mean + dbp_mean + 
                     age+weight+
                     gender+hypertension+charlson_comorbidity_index+
                     sapsii+diabetes+PT + PLT + WBC + Hb + Cr + AST + ALT+CA,
                   data = df)
print(summary(cox_365d2))
print(anova(cox_365d1, cox_365d2, test = "Chisq"))
interaction_results[["gcs_min_365d"]] <- anova(cox_365d1, cox_365d2, test = "Chisq")$P[2]
#npar*charlson_comorbidity_index
cox_icu1 <- coxph(Surv(time_icuD, icu_mortality) 
                  ~ npar + charlson_comorbidity_index + sbp_mean + dbp_mean + 
                    age+weight+
                    gender+hypertension+gcs_min+
                    sapsii+PT + PLT + WBC + Hb + Cr + AST + ALT+CA,
                  data = df)
print(summary(cox_icu2))
print(anova(cox_icu1, cox_icu2, test = "Chisq"))
interaction_results[["charlson_comorbidity_index_icuD"]] <- anova(cox_icu1, cox_icu2, test = "Chisq")$P[2]

# charlson_comorbidity_index分析30天死亡 -------------------------------------------------------------
cox_30d1 <- coxph(Surv(time_30d, stau_30d) 
                  ~ npar + charlson_comorbidity_index + sbp_mean + dbp_mean + 
                    age+weight+
                    gender+hypertension+gcs_min+
                    sapsii+PT + PLT + WBC + Hb + Cr + AST + ALT+CA,
                  data = df)
print(summary(cox_30d1))
cox_30d2 <- coxph(Surv(time_30d, stau_30d) 
                  ~ npar + charlson_comorbidity_index + npar*charlson_comorbidity_index + sbp_mean + dbp_mean + 
                    age+weight+
                    gender+hypertension+gcs_min+
                    sapsii+PT + PLT + WBC + Hb + Cr + AST + ALT+CA,
                  data = df)
print(summary(cox_30d2))
print(anova(cox_30d1, cox_30d2, test = "Chisq"))
interaction_results[["charlson_comorbidity_index_30d"]] <- anova(cox_30d1, cox_30d2, test = "Chisq")$P[2]
## 分析90天死亡 -------------------------------------------------------------
cox_90d1 <- coxph(Surv(time_90d, stau_90d) 
                  ~ npar + charlson_comorbidity_index + sbp_mean + dbp_mean + 
                    age+weight+
                    gender+hypertension+gcs_min+
                    sapsii+PT + PLT + WBC + Hb + Cr + AST + ALT+CA,
                  data = df)
print(summary(cox_90d1))
cox_90d2 <- coxph(Surv(time_90d, stau_90d) 
                  ~ npar + charlson_comorbidity_index + npar*charlson_comorbidity_index + sbp_mean + dbp_mean + 
                    age+weight+
                    gender+hypertension+gcs_min+
                    sapsii+PT + PLT + WBC + Hb + Cr + AST + ALT+CA,
                  data = df)
print(summary(cox_90d2))
print(anova(cox_90d1, cox_90d2, test = "Chisq"))
interaction_results[["charlson_comorbidity_index_90d"]] <- anova(cox_90d1, cox_90d2, test = "Chisq")$P[2]
## 分析365天死亡 ------------------------------------------------------------
cox_365d1 <- coxph(Surv(time_365d, stau_365d) 
                   ~ npar + charlson_comorbidity_index + sbp_mean + dbp_mean + 
                     age+weight+
                     gender+hypertension+gcs_min+
                     sapsii+PT + PLT + WBC + Hb + Cr + AST + ALT+CA,
                   data = df)
print(summary(cox_365d1))
cox_365d2 <- coxph(Surv(time_365d, stau_365d) 
                   ~ npar + charlson_comorbidity_index + npar*charlson_comorbidity_index + sbp_mean + dbp_mean + 
                     age+weight+
                     gender+hypertension+gcs_min+
                     sapsii+ PT + PLT + WBC + Hb + Cr + AST + ALT+CA,
                   data = df)
print(summary(cox_365d2))
print(anova(cox_365d1, cox_365d2, test = "Chisq"))
interaction_results[["charlson_comorbidity_index_365d"]] <- anova(cox_365d1, cox_365d2, test = "Chisq")$P[2]
interaction_df <- data.frame(
  Term = names(interaction_results),
  P_Value = unlist(interaction_results)
)
group_vars <- c("age", "gender", "sapsii", "gcs_min", 
                "charlson_comorbidity_index", "hypertension", "diabetes")

# 自动生成分组统计
generate_group_counts <- function(var) {
  if(var %in% c("hypertension", "diabetes")) { # 二分类变量处理
    counts <- table(df[[var]]) %>% paste(names(.), ., sep = ": ")
  } else { # 其他分组变量处理
    counts <- table(df[[var]]) %>% 
      as.data.frame() %>% 
      mutate(Label = paste0(Var1, ": ", Freq)) %>% 
      pull(Label)
  }
  return(paste(counts, collapse = ", "))
}

# 生成统计字符串
group_counts <- sapply(group_vars, generate_group_counts)

# 创建统计结果数据框
counts_df <- data.frame(
  Term = paste0(group_vars, "_group"),
  Group_Counts = group_counts
)

# 合并结果 -----------------------------------------------------------------
# 创建对应关系（根据你实际检验的变量名）
term_mapping <- c(
  "age_icuD" = "age_group",
  "gender_icuD" = "gender_group",
  "sapsii_icuD" = "sapsii_group",
  "gcs_min_icuD" = "gcs_min_group",
  "charlson_comorbidity_index_icuD" = "charlson_comorbidity_index_group",
  "hypertension_icuD" = "hypertension_group",
  "diabetes_icuD" = "diabetes_group"
)

# 添加分组统计到结果
final_df <- interaction_df %>% 
  mutate(Group_Var = term_mapping[Term]) %>% 
  left_join(counts_df, by = c("Group_Var" = "Term")) %>% 
  dplyr::select(-Group_Var)
write.csv(final_df, file ="interaction_results.csv")

#3\亚组分析计算可信区间
y2 <- y1%>%
  dplyr::select(-X)
df <-y2%>%
  mutate(age=ifelse(age >74,">74","<=74"),
         age=factor(age, levels=c(">74","<=74")),
         sapsii=ifelse(sapsii >37,">37","<=37"),
         sapsii=factor(sapsii, levels=c(">37","<=37")),
         gcs_min=ifelse(gcs_min >14,">14","<=14"),
         gcs_min=factor(gcs_min, levels=c(">14","<=14")),
         charlson_comorbidity_index= ifelse(charlson_comorbidity_index >6,">6","<=6"),
         charlson_comorbidity_index=factor(charlson_comorbidity_index, levels=c(">6","<=6")),
         npar_ratio=cut(npar_ratio,
                        breaks =c(3.870968,20.769673,
                                  26.359122 ,100.625000 ),
                        levels = c(1,2,3)))
library(jstable)  
time_points <- c("icu","30d", "90d", "365d")
analyses <- list(
  list(name = "exclude3", exclude_value = 3),
  list(name = "exclude2", exclude_value = 2)
)
if(exists("res")) rm(res)
for (time_suffix in time_points) {
  time_num <- if (time_suffix == "icu") "0" else gsub("\\D", "", time_suffix)
  time_var <- if(time_suffix == "icu") {
    "time_icuD"
  } else {
    paste0("time_", time_suffix)
  }
  status_var <- paste0("stau_", time_suffix)
  for (analysis in analyses) {
    obj_name <- paste0("res", time_num, "_", analysis$name)
    assign(obj_name, 
           local({
             filtered_df <- df %>%
               mutate(npar_ratio = as.numeric(npar_ratio)) %>%
               filter(npar_ratio != analysis$exclude_value)
             formula_str <- paste0(
               "Surv(", time_var, ", ", status_var, ") ~ npar_ratio"
             )
             TableSubgroupMultiCox(
               as.formula(formula_str),
               var_subgroups = c("age","gender","hypertension","diabetes",
                                 "sapsii","gcs_min","charlson_comorbidity_index"
                                 ),
               var_cov = c("sbp_mean","dbp_mean","WBC","Hb","PLT","Cr","gender",
                           "hypertension","diabetes","sapsii","gcs_min","age","weight",
                           "charlson_comorbidity_index","ALT","AST","PT","CA"),
               data = filtered_df
             )
           }))
    cat("\n===== Results for", obj_name, "======\n")
    print(get(obj_name))  # 直接打印结果对象
    write.csv(get(obj_name), 
              file = paste0(time_suffix, "_", analysis$name, "_subgroup.csv"))
  }
}
ls(pattern = "^res\\d+")

rm(list = ls())
#7\限制性立方样条(RCS)曲线图------------
#导入包和数据
library(dplyr)
library(rms)
library(foreign)
library(survival)
library(survminer)#曲线
library(ggplot2) #作图
library(ggsci) #调色板
library(haven) #读取数据集
completedData_mod<- read.csv("completedData_mod.csv",header = TRUE,sep = ",")
df <-completedData_mod%>%
  mutate(gender=factor(gender,labels=c("F","M")),
         hypertension=factor(hypertension,levels=c(0,1),labels=c("No","Yes")),
         diabetes=factor(diabetes,levels=c(0,1),labels=c("No","Yes")),
         npar=factor(npar,labels = c("T1","T2","T3")))
df
dd<-datadist(df) 
options(datadist='dd')
head(df)
str(df)
splots <- list()

#ICU_mortality
fit1<- cph(Surv(time_icuD,icu_mortality)~rcs(npar_ratio,3),
           data=df)
anova(fit1) 
fit1

HR1<-Predict(fit1,npar_ratio,fun=exp,ref.zero=TRUE)
HR1
find_hr1_cross <- function(hr_data) {
  # 找到第一个跨越HR=1的区间
  below <- hr_data$yhat < 1
  cross_idx <- which(diff(below) != 0)[1]
  
  if(is.na(cross_idx)) return(NA)  # 无交叉时返回NA
  
  # 线性插值
  x1 <- hr_data$npar_ratio[cross_idx]
  x2 <- hr_data$npar_ratio[cross_idx+1]
  y1 <- hr_data$yhat[cross_idx]
  y2 <- hr_data$yhat[cross_idx+1]
  
  x_cross <- x1 + (1 - y1)*(x2 - x1)/(y2 - y1)
  round(x_cross, 1)  # 保留1位小数
}
cross1 <- find_hr1_cross(HR1)
splots[[1]]<-ggplot()+
  geom_line(data=HR1,aes(npar_ratio,yhat),linetype="solid", 
            size=1,alpha=0.7,colour="red")+
  geom_ribbon(data=HR1,
              aes(npar_ratio,ymin=lower,ymax=upper),
              alpha=0.1,fill="red")+
  theme_classic()+
  geom_hline(yintercept=1,linetype=2,size=0.75)+
  geom_segment(aes(x = cross1, y =0, xend = cross1, yend =1),  
               colour = "black", linetype=2, size=0.75) +
  annotate("text", x =1.2, y=3, 
           label = paste0("Threshold = ", cross1, "\nP for non-linearity=0.109"),  # 添加数值
           colour="black", size=4)+
  scale_x_continuous(breaks=c(seq(10,60,10),23.3),expand = expansion(add = c(2, 5))
                     ,guide = guide_axis(angle = 30))+
  scale_y_continuous(breaks=seq(0,4,1),expand = c(0,0),limits = c(0, 4) )+
  labs(title="ICU mortality",x="npar",
       y="HR (95%CI)")

splots[[1]]  
#stau_30d
df <- df %>%
  mutate(across(c(stau_30d, stau_90d, stau_365d), ~ ifelse(.x == "yes", 1, 0)))
fit3<- cph(Surv(time_30d,stau_30d)~rcs(npar_ratio,4),
           data=df)
anova(fit3) #p<0.001
fit3

HR3<-Predict(fit3,npar_ratio,fun=exp,ref.zero=TRUE)
HR3
cross3 <- find_hr1_cross(HR3)
splots[[3]]<-ggplot()+
  geom_line(data=HR3,aes(npar_ratio,yhat),linetype="solid", 
            size=1,alpha=0.7,colour="red")+
  geom_ribbon(data=HR3,
              aes(npar_ratio,ymin=lower,ymax=upper),
              alpha=0.1,fill="red")+
  theme_classic()+
  geom_hline(yintercept=1,linetype=2,size=0.75)+
  geom_segment(aes(x = cross3, y =0, xend = cross3, yend =1),  # 同理替换
               colour = "black", linetype=2, size=0.75) +
  annotate("text", x =1.2, y=3,
           label = paste0("Threshold = ", cross3, "\nP <0.001"))+
  scale_x_continuous(breaks=c(seq(10,60,10),23.3),expand = expansion(add = c(2, 5))
                     ,guide = guide_axis(angle = 30))+
  scale_y_continuous(breaks=seq(0,4,1),expand = c(0,0),limits = c(0, 4) )+
  labs(title="30d mortality",x="npar",
       y="HR (95%CI)")
splots[[3]] 

#stau——90d
fit4<- cph(Surv(time_90d,stau_90d)~rcs(npar_ratio,3),
           data=df)
anova(fit4) #p<0.001
fit4

HR4<-Predict(fit4,npar_ratio,fun=exp,ref.zero=TRUE)
HR4
cross4 <- find_hr1_cross(HR4)
splots[[4]]<-ggplot()+
  geom_line(data=HR4,aes(npar_ratio,yhat),linetype="solid", 
            size=1,alpha=0.7,colour="red")+
  geom_ribbon(data=HR4,
              aes(npar_ratio,ymin=lower,ymax=upper),
              alpha=0.1,fill="red")+
  theme_classic()+
  geom_hline(yintercept=1,linetype=2,size=0.75)+
  geom_segment(aes(x = cross4, y =0, xend = cross4, yend =1),  # 同理替换
               colour = "black", linetype=2, size=0.75) +
  annotate("text", x =1.2, y=3,
           label = paste0("Threshold = ", cross4, "\nP <0.001"))+
  scale_x_continuous(breaks=c(seq(10,60,10),23.3),expand = expansion(add = c(2, 5))
                     ,guide = guide_axis(angle = 30))+
  scale_y_continuous(breaks=seq(0,4,1),expand = c(0,0),limits = c(0, 4) )+
  labs(title="90d mortality",x="npar",
       y="HR (95%CI)")
splots[[4]]
fit5<- cph(Surv(time_365d,stau_365d)~rcs(npar_ratio,3),
           data=df)
anova(fit5) #p<0.001
fit5

HR5<-Predict(fit5,npar_ratio,fun=exp,ref.zero=TRUE)
HR5
cross5 <- find_hr1_cross(HR5)
splots[[5]]<-ggplot()+
  geom_line(data=HR5,aes(npar_ratio,yhat),linetype="solid", 
            size=1,alpha=0.7,colour="red")+
  geom_ribbon(data=HR5,
              aes(npar_ratio,ymin=lower,ymax=upper),
              alpha=0.1,fill="red")+
  theme_classic()+
  geom_hline(yintercept=1,linetype=2,size=0.75)+
  geom_segment(aes(x = cross5, y =0, xend = cross5, yend =1),  # 同理替换
               colour = "black", linetype=2, size=0.75) +
  annotate("text", x =1.2, y=3,
           label = paste0("Threshold = ", cross5, "\nP <0.001"))+
  scale_x_continuous(breaks=c(seq(10,60,10),23.3),expand = expansion(add = c(2, 5))
                     ,guide = guide_axis(angle = 30))+
  scale_x_continuous(breaks=c(seq(10,60,10),23.3),expand = expansion(add = c(2, 5))
                     ,guide = guide_axis(angle = 30))+
  scale_y_continuous(breaks=seq(0,4,1),expand = c(0,0),limits = c(0, 4) )+
  labs(title="365d mortality",x="npar",
       y="HR (95%CI)")
splots[[5]]

#组合图片
rcs <- ggarrange(splots[[1]],splots[[3]],splots[[4]],splots[[5]])
rcs
ggsave("RCS_mod.pdf", dpi = 600,width =14, height = 12,rcs)
ggsave("RCS_mod.png", dpi = 600,width =14, height = 12,rcs)

