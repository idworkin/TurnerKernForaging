
## ----libraries-----------------------------------------------------------
require(Rcpp)
require(lme4)
require(lattice)


## ----get_data_in---------------------------------------------------------
tk_path_lengths <- read.csv("../TurnerKern_data/TurnerKernForaging.csv", h=T)

# Checking data input
summary(tk_path_lengths)
str(tk_path_lengths)
head(tk_path_lengths)
tail(tk_path_lengths)

# Checking to see if there is any missing data
colSums(is.na(tk_path_lengths))


## ----tables_date---------------------------------------------------------
table(tk_path_lengths$Date) 


## ----tables_line---------------------------------------------------------
table(tk_path_lengths$Line) 


## ----tables_line_date----------------------------------------------------
table(tk_path_lengths[,1:2]) 


## ----boxplot_line--------------------------------------------------------
plot(Distance ~ Line, data=tk_path_lengths, las=2)


## ----boxplot_date--------------------------------------------------------
plot(Distance ~ as.factor(Date), data=tk_path_lengths, las=2)


## ----bwplot_line_date----------------------------------------------------
bwplot(Distance ~ as.factor(Date)|Line, data=tk_path_lengths, las=2)


## ----mixed_mod1----------------------------------------------------------
tk_path_lengths$Temp_c <- scale(tk_path_lengths$Temp, center=TRUE, scale=FALSE)

mod1 <- lmer(Distance ~ Temp_c + (1|Date) + (1|Line), data=tk_path_lengths, REML=TRUE)

# Let's take a look at the summary from this model.

summary(mod1)


## ----blup_plot_line------------------------------------------------------
dotplot(ranef(mod1,condVar=TRUE),
              scales = list(x = list(relation = 'free')), strip=FALSE)[["Line"]]


## ----blup_plot_date------------------------------------------------------
dotplot(ranef(mod1,condVar=TRUE),
              scales = list(x = list(relation = 'free')), strip=FALSE)[["Date"]]


## ----model_checks--------------------------------------------------------
mod1_noDate <- lmer(Distance ~ Temp_c + (1|Line), data=tk_path_lengths, REML=TRUE)
mod1_noLine <- lmer(Distance ~ Temp_c + (1|Date), data=tk_path_lengths, REML=TRUE)
mod1_noTemp <- lmer(Distance ~ 1 + (1|Date) + (1|Line), data=tk_path_lengths, REML=TRUE)


## ----model_check2--------------------------------------------------------
mod2 <- lmer(Distance ~ Temp_c + (1|Line) + (1|Date) + (1|Date:Line), data=tk_path_lengths, REML=TRUE)


## ----date_line_plot------------------------------------------------------
dotplot(ranef(mod2,condVar=TRUE),
              scales = list(x = list(relation = 'free')), strip=FALSE)[["Date:Line"]]


## ----plots_from_interaction_model----------------------------------------
dotplot(ranef(mod2,condVar=TRUE),
              scales = list(x = list(relation = 'free')), strip=FALSE)[["Date"]]
dotplot(ranef(mod2,condVar=TRUE),
              scales = list(x = list(relation = 'free')), strip=FALSE)[["Line"]]


## ----model_examination---------------------------------------------------
AIC(mod1, mod1_noDate, mod1_noLine, mod1_noTemp, mod2)
BIC(mod1, mod1_noDate, mod1_noLine, mod1_noTemp, mod2)
logLik(mod1)
logLik(mod1_noDate)
logLik(mod1_noLine)


## ----confint-------------------------------------------------------------
CI_mod1 <- confint(mod1)
CI_mod1


## ----internal_checks-----------------------------------------------------
coef_mod1 <- coef(mod1)
lines_mod1 <- as.data.frame(coef_mod1$Line[,1])
rownames(lines_mod1) <- rownames(coef_mod1$Line)

### confirming that the above coefs are the BLUPs + intercept. Confirmed...
ranef_mod1 <- ranef(mod1)
lines_ranef_mod1 <- unlist(ranef_mod1$Line)
cor(lines_mod1, lines_ranef_mod1)
plot(lines_mod1[,1], lines_ranef_mod1, pch=20)


## ----blups_line_means----------------------------------------------------
# This computes raw line means, ignoring sample sizes or dates.
line_means <- tapply(tk_path_lengths$Distance, INDEX=tk_path_lengths$Line, FUN=mean)

cor(line_means,lines_mod1)
plot(line_means,lines_mod1[,], pch=20)


## ----multicolin----------------------------------------------------------
# estimate in lm
mod3_lm <- lm(Distance ~ 1 + Temp_c + Line + Date, data=tk_path_lengths)
mod3_lm_noDate <- lm(Distance ~ 1 + Temp_c + Line, data=tk_path_lengths)
mod3_lm_noLine <- lm(Distance ~ 1 + Temp_c + Date, data=tk_path_lengths)


## ----condNum-------------------------------------------------------------
ConditionNumber <- function(model.x=model.object){
	x <- model.matrix(model.x)
    eigen.x <- eigen(t(x)%*%x)
    sqrt(max(eigen.x$values)/min(eigen.x$values))
}


## ----condNum_checks------------------------------------------------------
ConditionNumber(mod3_lm)
ConditionNumber(mod3_lm_noDate)
ConditionNumber(mod3_lm_noLine)


## ----vif-----------------------------------------------------------------
car::vif(mod3_lm) 
car::vif(mod3_lm_noDate)
car::vif(mod3_lm_noLine) 


## ----read_check_ms-------------------------------------------------------
ms_path_lengths <- read.csv("../TurnerKern_data/MS_DGRP_foraging.csv", h=T)

# Checking data input
summary(ms_path_lengths)
str(ms_path_lengths)
head(ms_path_lengths)
tail(ms_path_lengths)

# Change BB to 'r' for rover
levels(ms_path_lengths$Line)[34] <- "r"

#Change ee to 's' for sitter
levels(ms_path_lengths$Line)[35] <- "s"

# Checking to see if there is any missing data
colSums(is.na(ms_path_lengths))


## ----tables_data_ms------------------------------------------------------
table(ms_path_lengths$Day) 


## ----tables_lines_ms-----------------------------------------------------
table(ms_path_lengths$Line) 


## ----tables_lines_days_ms------------------------------------------------
table(ms_path_lengths[,c(1,4)])


## ----boxplot_line_ms-----------------------------------------------------
plot(Pathlengthx2 ~ Line, data=ms_path_lengths, las=2, 
    main="boxplot of path lengths for each line, unadjusted")


## ----boxplot_Day_ms------------------------------------------------------
plot(Pathlengthx2 ~ as.factor(Day), data=ms_path_lengths, las=2,
    main="boxplot of path lengths for each Day, unadjusted")


## ----bwplot_line_day_ms--------------------------------------------------
bwplot(Pathlengthx2 ~ as.factor(Day)|Line, data=ms_path_lengths, las=2)


## ----mod1_ms-------------------------------------------------------------
ms_mod1 <- lmer(Pathlengthx2 ~ 1 + (1|Day) + (1|Line), data=ms_path_lengths, REML=TRUE)

# Let's look at some summaries.
summary(ms_mod1)


## ----plots_blups1_ms-----------------------------------------------------
dotplot(ranef(ms_mod1,condVar=TRUE), strip=FALSE,
              scales = list(x = list(relation = 'free')))[["Line"]]


## ----plots_blups2_ms-----------------------------------------------------
dotplot(ranef(ms_mod1,condVar=TRUE), strip=FALSE,
              scales = list(x = list(relation = 'free')))[["Day"]]


## ----colinearity_ms------------------------------------------------------
ms_mod1_lm <- lm(Pathlengthx2 ~ 1 + as.factor(Day) + Line, data=ms_path_lengths)

ConditionNumber(ms_mod1_lm)
car::vif(ms_mod1_lm) 


## ----compare_between_sets------------------------------------------------
# Some set up (not really necessary)
ms_mod1_coef <- coef(ms_mod1)
ms_line <- as.data.frame(ms_mod1_coef$Line)
colnames(ms_line) <- "blup"

line_names <- as.character(rownames(ms_line))

ms_line <- as.data.frame(cbind(ms_line, line_names))
tt_line <- as.data.frame(cbind(lines_mod1, names =rownames(lines_mod1)))

merged_data <- merge(ms_line,lines_mod1,by="row.names")


## ----compare_plot1-------------------------------------------------------
plot(merged_data[,c(2,4)], xlab="Sokolowski", ylab="Turner_Kern", pch=20,
    main="Correlation of BLUPs for DGRP lines as measured across two labs")

# Consistent with plot the data is only moderately correlated.
cor.test(x=merged_data[,2], y=merged_data[,4])


## ----compare_line_means_across-------------------------------------------
tt_line_means <- tapply(tk_path_lengths$Distance, INDEX=tk_path_lengths$Line, FUN=mean)
ms_line_means <- tapply(ms_path_lengths$Pathlengthx2, INDEX=ms_path_lengths$Line, FUN=mean)

merged_line_means <- merge(tt_line_means, ms_line_means, by="row.names")
cor.test(merged_line_means[,2], merged_line_means[,3])

plot(merged_line_means[,c(3,2)], pch=20,
    main="Correlation of line means (unadjusted) for DGRP lines as measured across two labs",
    xlab="Sokolowski", ylab="Turner_Kern")


