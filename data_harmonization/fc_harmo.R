# install.packages("devtools")
# devtools::install_github("andy1764/CovBat_Harmonization/R")
# devtools::install_github("andy1764/FCHarmony")

library(FCHarmony)
library(R.matlab)


#-------------------------- Please Load Your Data, including FC, site label, age, and sex
FC = readMat("Your_FC_data.mat") # FC is a d * d * S matrix, where d is the number of ROIs and S is the number of subjects.
Age =  readMat("Your_age_data.mat") # Age is a S * 1 vector.
Sex = readMat("Your_sex_data.mat") # Sex is a S * 1 vector.
Site = readMat("Your_site_data.mat") # Site is a S * 1 vector.

#-------------------------- Run the following code.
batch=factor(Site)

X = data.frame(Age=Age, Age2 = Age^2, Sex=Sex, Sex_Age = Sex * Age, Sex_Age2 = Sex * Age^2)
XX = model.matrix(~.,X)

res=fcComBat(FC,batch, mod=XX, to.corr=FALSE, fisher=FALSE) # if your FC is not Fisher-transformed, please set fisher = TRUE.

FC_aligned = res$dat.out






