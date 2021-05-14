# Srikar Kavirayuni
#
# clears and set working directory
rm(list=ls())
setwd("C:/Users/babys/Desktop/R")

#
#  Code goes here

#Loading the library after installing the QTL package
library(qtl)
library(dplyr)

################# START PROJECT

########## Initial Data checking
data <- read.cross("csv", 
                    file = "B6xA.csv", 
                    genotypes = c("A","H","B"), 
                    na.strings = "-", 
                    alleles = c("A","B"))
jittermap(data)
class(data)
names(data)
names(data)
names(data$pheno)
names(data$geno)

#checking if data is clean
#recombinant frequency plot
data <- est.rf(data)
plotRF(data)
#genetic map
plot.map(data)
#missing data map
plotMissing(data)

#phenotype distributions
TG <- data$pheno$TG
TG

## normalize TG column
TG_normal = log10(TG+1)
hist(log10(TG+1))

#data quality checks
qqnorm(log10(TG+1))
qqline(log10(TG+1))

hist(log10(data$pheno$TG))
hist(data$pheno$TG)

########### single qtl analysis -- interval mapping

data <- calc.genoprob(data, step=1)
data <- calc.genoprob(data, step=1,error.prob =0.0001)
out <- scanone(data,pheno.col = 5)
plot(out)
iplotScanone(out)
iplot_data <- iplotScanone(out, data)

## em
out.em <- scanone(data, pheno.col = 5, method ="em")
plot(out,out.em)
plot(out,out.em, chr=c(11))
plot(out.em,ylab="LOD SCore")

## ehk
out.ehk <- scanone(data, pheno.col = 5, method ="ehk")
plot(out,out.ehk)
plot(out,out.ehk, chr=c(11))
plot(out.ehk,ylab="LOD SCore")

## main scan
operm.em <- scanone(data, method="em",pheno.col = 5, n.perm = 100)
operm.em
plot(operm.em)
thresh = summary(operm.em,alpha = c(0.63,0.10,0.05))
thresh
plot(out.em)

abline(h=thresh[1], col = "blue")
abline(h=thresh[2], col = "red")
abline(h=thresh[3], col = "green")

## main scan using the extended haley knot
operm.ehk <- scanone(data, method="ehk",pheno.col = 5, n.perm = 100)
operm.ehk
plot(operm.ehk)
summary(operm.ehk,alpha = c(0.63,0.10,0.05))
plot(out.ehk)

########## comparing various methods

### em
data.em=calc.genoprob(data, step=1,error.prob = .00001)
out.em=scanone(data, pheno.col=5,method="em")
plot(out.em, ylab="LOD Score")

### hk
data.hk=calc.genoprob(data, step=1,error.prob = .00001)
out.hk=scanone(data, pheno.col=5,method="hk")
plot(out.hk, ylab="LOD Score")

### comparing em and hk
plot(out.em, out.hk, col=c("blue","red"),ylab="LOD score")
plot(out.hk - out.em, ylim=c(-0.5, 1.0), ylab=expression(LOD[HK] - LOD[EM]))
abline(h=0, lty=3)

### ehk
data.ehk=calc.genoprob(data, step=1,error.prob = .00001)
out.ehk=scanone(data, pheno.col=5,method="ehk")
plot(out.ehk, ylab="LOD Score")

### comparing em,hk, and ehk
plot(out.em, out.hk, out.ehk, ylab="LOD score",lty=c(1,2,3))

plot(out.ehk - out.em,ylim=c(-0.5, 1.0), ylab=expression(LOD[EHK] - LOD[EM]))
abline(h=0, lty=3)

plot(out.ehk - out.hk, ylim=c(-0.5, 1.0), ylab=expression(LOD[EHK] - LOD[HK]))
abline(h=0, lty=3)

### individual chromosomes -- big spike on chromosome 11
data=calc.genoprob(data,step=1,error.prob = 0.00001)
out11.em=scanone(data,chr=11,method="em")
out11.hk=scanone(data,chr=11,method="hk")
out11.ehk=scanone(data,chr=11,method="ehk")
plot(out11.em,out11.hk,out11.ehk,ylab="LOD Score",lty=c(1,2,3))



### Even though there is a spike in chromosome 11, it's not statistically significant
### This can be confirmed with the p-value and LOD scores of the mainscan.
summary(out.em,perms=operm.em,alpha=0.05)
# There were no LOD peaks above the threshold for .05 alpha
summary(out.em,perms=operm.em,alpha=0.10)
# There were no LOD peaks above the threshold for .10 alpha
summary(out.em,perms=operm.em,alpha=0.63, pvalues = TRUE)
# There is one LOD peak with significance of .63 meaning there is genetic loci SOMEWHERE in the chromosome.
# The p-value is .31 which isn't significant. 


###############################################333#######
## Single QTL ANALYSIS INTERVAL MAPPING WITH HALEY_KNOT REGRESSION
### hk -- mainscan using haley-knot regression
data <- calc.genoprob(data, step = 2.0, off.end = 0.0, error.prob = 1.0e-4, map.function = "haldane", stepwidth = "fixed")
data <- sim.geno(data, step = 2.0, off.end = 0.0, error.prob = 1.0e-4, map.function = "haldane", stepwidth = "fixed")
#performing the mainscan      
out.hk <- scanone(data, pheno.col = 5, model = "normal", method = "hk")
operm.hk <- scanone(data, pheno.col = 5, model = "normal", method = "hk", n.perm = 1000)
#plotting the mainscan
plot(out.hk, main = "Mainscan plot of Triglycerides in B6 x A/J")
#threshold lines
thresh <- summary(operm.hk, alpha = c(0.63,0.10,0.05))
thresh
#LOD thresholds (1000 permutations)
#      lod
# 63% 2.14
# 10% 3.15
# 5%  3.41
abline(h=thresh[1], col = "blue")
abline(h=thresh[2], col = "red")
abline(h=thresh[3], col = "green")
# summary of mainscan output
summary(out.hk,perms=operm.hk,alpha=.05,pvalues=TRUE)
## There were no LOD peaks above the threshold.


##### QTL EFFECTS, ADDITIVE COVARIATES AND INTERACTIVE COVARIATES


### TEST FOR COVARIATE PRESENCE OF SEX
SEX_PHENOTYPE = data$pheno$sex
boxplot(TG ~ SEX_PHENOTYPE, data=data$pheno,
          horizontal=TRUE, xlab="TG Levels",
          col=c("red","blue", 
          main = "Boxplot of TG levels grouped by sex"))
### Males have more triglyceride levels than females on average.
anova(aov(TG ~ SEX_PHENOTYPE, data=data$pheno))
### VERY SIGNIFICANT!!! SEX IS DEFINITELY HAVING AN EFFECT AND 
# CAN BE CONSIDERED A COVARIATE IN THIS QTL MAPPING STUDY!!!

#### Analysis of Variance Table

# Response: TG
#                Df  Sum Sq Mean Sq F value    Pr(>F)    
# SEX_PHENOTYPE   1  931302  931302  395.32 < 2.2e-16 ***
# Residuals      612 1441765    2356                      
# ---------------------------------------------------------------

### ADDED SEX COVARIATE
data <- calc.genoprob(data, step = 2.0, off.end = 0.0, error.prob = 1.0e-4, map.function = "haldane", stepwidth = "fixed")
data <- sim.geno(data, step = 2.0, off.end = 0.0, error.prob = 1.0e-4, map.function = "haldane", stepwidth = "fixed")
#performing the mainscan      
out.addsex <- scanone(data, pheno.col = 5, model = "normal", method = "hk", addcovar=data$pheno$sex)
operm.addsex <- scanone(data, pheno.col = 5, model = "normal", method = "hk", n.perm = 1000, addcovar=data$pheno$sex)
#plotting the mainscan
plot(out.addsex, main = "Mainscan plot of Triglycerides in B6 x A/J with Additive Sex Covariate Model")

#threshold lines
thresh <- summary(operm.addsex, alpha = c(0.63,0.10,0.05))
thresh
#LOD thresholds (1000 permutations)
#     lod
#63% 2.11
#10% 3.20
#5%  3.38
abline(h=thresh[1], col = "blue")
abline(h=thresh[2], col = "red")
abline(h=thresh[3], col = "green")

#summary of mainscan (text output)
summary(out.addsex,perm=operm.addsex,lodcolumn = 1, alpha=0.05)

#          chr  pos  lod
#c10.loc56  10 58.1 3.73


## There was 1 LOD peak above the threshold!

#finding markers - do an effect plot
mname1 <- find.marker(data,chr=10,pos=58.1)
effectplot(data,pheno.col=5, mname1 = mname1)

#finding markers - confidence intervals
CIchr10 <- bayesint(out.addsex, chr=10, prob=.05)
plot(out.addsex, chr=10, lodcolumn = 1, main = "Confidence interval for c10.loc56")
lines(x=CIchr10[c(1,3),2], y=c(0,0), type="l", col="green", lwd=4)
CIchr10[c(1,3),2]
### So the position of the genetics loci is 58.061 cM in chromosme 10.

### INTERACTIVE SEX COVARIATE
data <- calc.genoprob(data, step = 2.0, off.end = 0.0, error.prob = 1.0e-4, map.function = "haldane", stepwidth = "fixed")
data <- sim.geno(data, step = 2.0, off.end = 0.0, error.prob = 1.0e-4, map.function = "haldane", stepwidth = "fixed")
#performing the mainscan      
out.INTsex <- scanone(data, pheno.col = 5, model = "normal", method = "hk", intcovar=data$pheno$sex)
operm.INTsex <- scanone(data, pheno.col = 5, model = "normal", method = "hk", n.perm = 1000, intcovar=data$pheno$sex)
#plotting the mainscan
plot(out.INTsex, main = "Mainscan plot of Triglycerides in B6 x A/J with INTERACTIVE Sex Covariate Model")

#threshold lines
thresh <- summary(operm.INTsex, alpha = c(0.63,0.10,0.05))
thresh

#LOD thresholds (1000 permutations)
#     lod
#63% 3.33
#10% 4.80
#5%  5.34

abline(h=thresh[1], col = "blue")
abline(h=thresh[2], col = "red")
abline(h=thresh[3], col = "green")

#summary of mainscan (text output)
summary(out.INTsex,perm=operm.INTsex,lodcolumn = 1, alpha=0.05)
#           chr  pos   lod
# c1.loc72   1  73.9  5.92 

## THERE WAS ONE SIGNIFICANT QTL! 

## Now plot the standard interval mapping (hk), additive and interactive sex covariate mainscans together:
plot(out.hk,out.addsex,out.INTsex,main="Plot of Standard Single QTL Analysis with Additive and Interactive Sex Covariates")

#finding markers - do an effect plot
mname1 <- find.marker(data,chr=1,pos=73.9)
effectplot(data,pheno.col=5, mname1 = mname1)

#finding markers - confidence intervals
CIchr1 <- bayesint(out.INTsex, chr=1, prob=.05)
plot(out.INTsex, chr=1, lodcolumn = 1, main = "Confidence interval for c1.loc72")
lines(x=CIchr1[c(1,3),2], y=c(0,0), type="l", col="green", lwd=4)
CIchr1[c(1,3),2]
### So the position of the genetics loci is 73.89 cM in chromosme 1.

###### Two-dimesnsional, Two-QTL SCANS
data <- calc.genoprob(data, step = 2.0, off.end = 0.0, error.prob = 1.0e-4, map.function = "haldane", stepwidth = "fixed")
data <- sim.geno(data, step = 2.0, off.end = 0.0, error.prob = 1.0e-4, map.function = "haldane", stepwidth = "fixed")
#performing the mainscan      
out.two <- scantwo(data, pheno.col = 5, model = "normal", method = "hk")
operm.two <- scantwo(data, pheno.col = 5, model = "normal", method = "hk", n.perm = 100) # Note only 100 permutations were run for computational time
#plotting the mainscan
plot(out.two, main = "Mainscan plot of Triglycerides in B6 x A/J with two-dimensional QTL model")
plot(out.two,chr=c(1,10,11)) # the main chromosomes in the 1D scan
plot(out.two, chr=c(1,10,11),lower="fv1")
plot(out.two, chr=c(1,10,11),lower="fv1", upper="av1")
## We see a strong link for pair of QTL in chromosomes 1 and 10
## Let's check its significance:
summary(operm.two) # gives us thresholds for significance
#
# TG (100 permutations)
# full  fv1  int  add  av1  one
# 5%  9.27 7.78 6.76 5.64 3.26 3.43
# 10% 8.79 6.68 6.04 5.34 2.91 3.23
#
#
summary(out.two, thresholds = c(9.27,7.78,6.76,5.64,3.26))
summary(out.two)
## Don't observe any QTL Interactions or pairs of linked QTL
#
## HENCE, WILL NOW MOVE ON TO COVARIATE AND ADDITIVE VARIATES
## Additive Covariate Study of SEX
out.two <- scantwo(data, pheno.col = 5, model = "normal", method = "hk",addcovar = data$pheno$sex)
# plotting
plot(out.two, main = "Mainscan plot of Triglycerides in B6 x A/J with two-dimensional QTL model")
plot(out.two,chr=c(1,10,11)) # the main chromosomes in the 1D scan
## we can observe the same results as the one in the additive one-dimesional scan -- however, we also can see a possibility for a linked QTL between chromosomes 1 and 10.
## Try to figure out the significance thresholds first -- can visualize the results after we know significance
## RUN PERMUTATION TEST FOR SIGNIFCANCE THRESHOLDS
operm.two <- scantwo(data, pheno.col = 5, model = "normal", method = "hk", n.perm = 100, addcovar = data$pheno$sex)
summary(operm.two)
# TG (100 permutations)
# full  fv1  int  add  av1  one
# 5%  10.77 9.03 7.31 6.01 3.40 3.58
# 10%  9.62 7.71 6.86 5.50 3.14 3.18
summary(out.two)
summary(out.two, thresholds = c(10.77,9.03,7.31,6.01,3.40))
#        pos1f pos2f lod.full lod.fv1 lod.int     pos1a pos2a lod.add
# c1:c10  67.9  54.1     8.73    4.99    1.49      67.9  58.1    7.23
#       lod.av1
# c1:c10  3.5
## We can observe linked QTL in chromosmes 1 and 10.
# LET'S TRY TO VISUALIZE THIS
plot(out.two, chr=c(1,10),lower="fv1", upper="av1")
# A clear relationship is observed!!!
# FIND EFFECTS
effectplot(data,pheno.col = 5, mname1 = "1@73.9",mname2 = "10@58.1", main="Effect Plot for Sex Additive Covariate Effect")
## There is an ADDITIVE effect amongst both QTL when "addcovar = sex"
# EFFECTPLOT: the QTL on chromosomes 1 and 10 are seen to act
# approximately ADDITIVELY: the effect of the chromosome 10 locus is 
# the SAME for each of the THREE genotypes at the chromosome 1 locus, 
# and VICE-VERSA. Also, both loci have effects of the SAME SIGN:
# BB > AB > AA -- in regards to individuals' TG.
#####################################################################
#
#
#
# 
###### MULTIPLE-QTL MODEL FITTING
# Make QTL
data <- calc.genoprob(data, step = 2.0, off.end = 0.0, error.prob = 1.0e-4, map.function = "haldane", stepwidth = "fixed")
qtl <- makeqtl(data, chr=c(1, 10), pos=c(73.9,58.1),what="prob")
qtl
### QTL object containing imputed genotypes, with 16 imputations. 
#      name chr     pos n.gen
# Q1  1@73.9   1 73.890     3
# Q2 10@58.1  10 58.061     3
plot(qtl)
# Fit QTL to create a model
out.fq <- fitqtl(data, qtl=qtl, formula = y~Q1+Q2, method = "hk")
summary(out.fq)
## getting estimated effects
out.fq.effects <- fitqtl(data, qtl=qtl, formula = y~Q1+Q2,method = "hk",dropone = TRUE,get.ests = TRUE)
summary(out.fq.effects)
## addcovar = sex
Sex <- data.frame(Sex=pull.pheno(data, "sex"))
out.fq.addsex <- fitqtl(data, pheno.col=5, qtl, formula=y~Q1+Q2+Sex,
                  cov=Sex, method="hk")
summary(out.fq.addsex)
## intcovar = sex
out.fq.intsex <- fitqtl(data, pheno.col=5, qtl, formula=y~(Q1*Sex)+(Q2*Sex),
                     cov=Sex, method="hk")
summary(out.fq.intsex)
## intcovar = sex BUT just with "1@73.9"
out.fq.intsex <- fitqtl(data, pheno.col=5, qtl, formula=y~(Q1*Sex)+Q2,
                        cov=Sex, method="hk")
summary(out.fq.intsex)
## intcovar = sex AND intQTL = TRUE
out.fq.intsex <- fitqtl(data, pheno.col=5, qtl, formula=y~Q1*Q2*Sex,
                        cov=Sex, method="hk")
summary(out.fq.intsex)
# 
# The best fitted "initial" model was found and now the positions need to be confirmed in order to finish model analysis.
# Hence, will "refine" the QTL positions and run all models again.

################################### 

##REFINE QTL POSITIONS
rqtl <- refineqtl(data, qtl=qtl, formula=y~Q1+Q2, verbose=FALSE)
rqtl
#       name chr    pos n.gen
# Q1  1@25.9   1 25.890     3
# Q2 10@48.7  10 48.674     3
rqtl_addsex <- refineqtl(data, qtl=qtl, formula=y~Q1+Q2+Sex, verbose=FALSE,cov=Sex)
rqtl_addsex
#       name chr    pos n.gen
# Q1  1@50.7   1 50.703     3
# Q2 10@68.1  10 68.061     3
rqtl_intsex <- refineqtl(data, qtl=qtl, formula=y~(Q1+Q2)*Sex, verbose=FALSE,cov=Sex)
rqtl_intsex
#       name chr    pos n.gen
# Q1  1@50.7   1 50.703     3
# Q2 10@68.1  10 68.061     3
rqtl_intsex_q1 <- refineqtl(data, qtl=qtl, formula=y~(Q1*Sex)+Q2, verbose=FALSE,cov=Sex)
rqtl_intsex_q1
#       name chr    pos n.gen
# Q1  1@59.9   1 59.890     3
# Q2 10@40.7  10 40.665     3

## NOW WE NEED TO OBSERVE ALL THE QTL FITTED MODELS WITH EACH MODEL 
# AND SEE WHICH INCREASES THE LOD SCORE MOST SIGNIFICANTLY ON BOTH 
# THE INITIAL AND ADDITIVE SEX COVARIATE MODELS.

############# RQTL
out.rfq.effects <- fitqtl(data, qtl=rqtl, formula = y~Q1+Q2,method = "hk")
summary(out.rfq.effects)
## addcovar = sex
Sex <- data.frame(Sex=pull.pheno(data, "sex"))
out.rfq.addsex <- fitqtl(data, pheno.col=5, qtl=rqtl, formula=y~Q1+Q2+Sex,
                        cov=Sex, method="hk")
summary(out.rfq.addsex)

############# RQTL_ADDSEX
out.rfq.effects.ADD <- fitqtl(data, qtl=rqtl_addsex, formula = y~Q1+Q2,method = "hk")
summary(out.rfq.effects.ADD)
## addcovar = sex
Sex <- data.frame(Sex=pull.pheno(data, "sex"))
out.rfq.addsex.ADD <- fitqtl(data, pheno.col=5, qtl=rqtl_addsex, formula=y~Q1+Q2+Sex,
                             cov=Sex, method="hk")
summary(out.rfq.addsex.ADD)

############# RQTL_INTSEX
out.rfq.effects.INT <- fitqtl(data, qtl=rqtl_intsex, formula = y~Q1+Q2,method = "hk")
summary(out.rfq.effects.INT)
## addcovar = sex
Sex <- data.frame(Sex=pull.pheno(data, "sex"))
out.rfq.addsex.INT <- fitqtl(data, pheno.col=5, qtl=rqtl_intsex, formula=y~Q1+Q2+Sex,
                               cov=Sex, method="hk")
summary(out.rfq.addsex.INT)

############# RQTL_INTSEX_Q1
out.rfq.effects.INTQ1 <- fitqtl(data, qtl=rqtl_intsex_q1, formula = y~Q1+Q2,method = "hk")
summary(out.rfq.effects.INTQ1)
## addcovar = sex
Sex <- data.frame(Sex=pull.pheno(data, "sex"))
out.rfq.addsex.INTQ1 <- fitqtl(data, pheno.col=5, qtl=rqtl_intsex_q1, formula=y~Q1+Q2+Sex,
                         cov=Sex, method="hk")
summary(out.rfq.addsex.INTQ1)

########################################

## It seems that although, refine QTL will improve estimates of the QTL,
# the model fit of the initial QTL position observations are better than 
# that of the RefineQTL models. But to statistically speaking, its only 
# marginally worse for the ADD_sex, but very good for the initial model.

## BUT THERE IS ONE STRIKING EFFECT THAT IS SEEN. tHE ORIGINAL MODEL
# FROM THE QTL=QTL THAT WAS FOUND TO BE THE BEST, INT_Q1, WAS THE ONLY
# MODEL WITH ALL STATISTISCALLY SIGNIFICANT PARAMETERS INDICATING 
# REFINE QTL FOR THE THAT MODEL WILL MAYBE BE THE BEST. 

## LET'S TEST IT!!!

##########################################
out.rfq.effects <- fitqtl(data, qtl=rqtl, formula = y~Q1+Q2,method = "hk")
summary(out.rfq.effects)
# ---- lod SCORE IS 1.47 BUT NOT REALLY SIGNIFICANT

out.rfq.effects.INTQ1 <- fitqtl(data, qtl=rqtl_intsex_q1, formula = y~Q1+Q2,method = "hk")
summary(out.rfq.effects.INTQ1)
## addcovar = sex
Sex <- data.frame(Sex=pull.pheno(data, "sex"))
out.rfq.addsex.INTQ1 <- fitqtl(data, pheno.col=5, qtl=rqtl_intsex_q1, formula=y~Q1+Q2+Sex,
                               cov=Sex, method="hk")
summary(out.rfq.addsex.INTQ1)
## INTcovar = sex 
Sex <- data.frame(Sex=pull.pheno(data, "sex"))
out.rfq.intsex.INTQ1 <- fitqtl(data, pheno.col=5, qtl=rqtl_intsex_q1, formula=y~(Q1+Q2)*Sex,
                               cov=Sex, method="hk")
summary(out.rfq.intsex.INTQ1)
## INTcovar = sex BUT JUST FOR QTL1
Sex <- data.frame(Sex=pull.pheno(data, "sex"))
out.rfq.intsexQ1.INTQ1 <- fitqtl(data, pheno.col=5, qtl=rqtl_intsex_q1, formula=y~(Q1*Sex)+Q2,
                               cov=Sex, method="hk")
summary(out.rfq.intsexQ1.INTQ1)

#############################################

## IT LOOKS LIKE THE ORIGINAL MODEL IS THE BEST MODEL WITH REGARDS TO 
# LOD SCORE AND P-VALUE SIGNIFICANCE, BUT THIS NEW MODEL WITH RQTL
# HAS A POSITION 1@59.9*SEX AS NON-SIGNIFICANT, WHICH JUST SUGGESTS
# THE ADDITIVE MODEL FOR SEX ---- THAT'S HOW THE EFFECT WAS AS WELL.

#####################################################

# DOWN TO TWO MODEL POSSIBILITIES -- ADDSEX(RQTL) OR INTSEX_Q1(QTL)

## ADDSEX(RQTL)
Sex <- data.frame(Sex=pull.pheno(data, "sex"))
out.rfq.addsex.INTQ1 <- fitqtl(data, pheno.col=5, qtl=rqtl_intsex_q1, formula=y~Q1+Q2+Sex,
                               cov=Sex, method="hk")
summary(out.rfq.addsex.INTQ1)

## INTSEX_Q1(QTL) -- THIS IS THE BEST MODEL THAT EXPLAINS QTL
out.fq.intsex <- fitqtl(data, pheno.col=5, qtl=qtl, formula=y~(Q1*Sex)+Q2,
                        cov=Sex, method="hk")
summary(out.fq.intsex)

## But the one that makes sense is the intsex_q1(qtl) model. Why?
# Because, 
# 1. In the other model there is no significant interactive variate
# of sex which is offputting as the single-qtl scan with the intcovar
# = sex showed the interactive co-variate bring rise to QTL(chr=1).
# 2. Even though the effect plot of the two-qtl scan showed a strictly
# additive effect for the two linked QTLs in the presence of a sex
# additive covariate model, the interacting QTL model makes more sense
# because of the interaction significance shown both in the 1-qtl
# (mainscan) and multiple qtl model (anova variance table matrix).
# Also, the interactive Q1 model accounts for the additive effect but 
# also places an emphasis on Sex inteactinng with QTL 1, which was 
# statistically significant as shown in the anova above.
# 3. Plus the positions of the other one just don't make sense!


### AddInt Function -- Test for all QTL * QTL Interaction
out.ai = addint(data,pheno.col=5,qtl=qtl,formula = y~ Q1+Q2)
summary(out.ai)

## intcovar = sex
out.ai.intsex = addint(data, pheno.col=5, qtl=qtl, formula=y~(Q1*Sex)+Q2,
       cov=Sex, method="hk")
summary(out.ai.intsex)

###### We see in both cases -- there is no interaction of the QTL
## sex with mname1 is also not majorly significant. So maybe the 
# correct model is just all the additive cases -- which we will check. 

## addcovar = sex
Sex <- data.frame(Sex=pull.pheno(data, "sex"))
out.ai.addsex <- addint(data, pheno.col=5, qtl=qtl, formula=y~Q1+Q2+Sex,
                        cov=Sex, method="hk")
summary(out.ai.addsex)
### We see the significant interaction between chromosome 1 QTL
# and the sex covariate -- meaning the model holds true -- 
# interaction between QTL 1 and sex covariate...

## But which is the Q1 and Q2 -- Q1 = chr 1 and Q2 = chr 10. 


############# ADD QTL APPROACH
out.aq = addqtl(data,pheno.col=5,qtl=qtl,formula = y~ Q1+Q2)
summary(out.aq)

## intcovar = sex
out.aq.intsex = addqtl(data, pheno.col=5, qtl=qtl, formula=y~(Q1*Sex)+Q2,
                       cov=Sex, method="hk")
summary(out.aq.intsex)

### No additional QTL, one on the X-chromosome is tested, but...
# NOT significant. 
plot(out.aq, ylab = "LOD Score")
plot(out.aq.intsex, ylab = "LOD Score")

## It looks like there are two peaks on the chr 1 loci -- let's test!
############## ADD Pair -- Two dimensional approach
out.ap <- addpair(data, qtl=qtl, chr=1, formula=y~Q2, verbose=FALSE)
summary(out.ap)
## Little to no evidence of an additional QTL pair.

