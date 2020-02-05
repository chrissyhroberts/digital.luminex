library(tidyr)
library(psych)
library(ranger)

#apply the lxb.parse function to all lxb files in directory
lapply(lxbs,lxb.parse)

#aggregate the data
a<-lxb.aggregate.data(path = "data_files")
a$col<-substr(a$well,start = 1,stop = 1)
a$row<-substr(a$well,start = 2,stop = 3)

# get rid of standards
a<-subset(a, subset = row %in% c(3:12))

# trim data to beads we actually used
a<-subset(a, subset = Bead.ID %in% c(12:15,18,19,20,21,22,25:30))

#give beads better names
names(a)<-gsub(names(a),pattern = " ",replacement = ".")

bead.ag<-read.csv(file = "bead.ids.csv")

a<-merge(a,bead.ag)


a<- a %>% group_by(well,Bead.Ag) %>% summarise(mean.RPI=mean(RP1L.Peak,na.rm = T))

#spread data to get all on one line
dat <- a %>% spread(key = Bead.Ag, value = mean.RPI,fill = NA, convert = FALSE)
dat$row<-substr(dat$well,start = 2,stop = 3)


# start ranger analysis
modelling<-dat[,-1]
modelling$PHENOTYPE<-"CONTROL"
modelling$PHENOTYPE[which(modelling$row>7)]<-"CASE"
modelling$row<-NULL

names(modelling)<-gsub(x = names(modelling),pattern = "-",replacement = "")
ranger_results<-ranger(formula = PHENOTYPE ~ .,data=modelling,write.forest=T,importance="impurity",classification = T,num.trees = 1000000)


plot(sort(ranger_results$variable.importance,decreasing = T),pch=16,ylab="Variable Importance",xlab="SNPs")
xx<-sort(ranger_results$variable.importance,decreasing = T)
xx<-as.data.frame(xx)
xx


ranger_results$confusion.matrix
cohen.kappa(ranger_results$confusion.matrix)


ranger_results1<-ranger(formula = PHENOTYPE ~ IL6+IP10+CRP+PCT,data=modelling,write.forest=T,importance="impurity")



ranger_results1$confusion.matrix
cohen.kappa(ranger_results1$confusion.matrix)


ranger_results2$classification.table
cohen.kappa(ranger_results2$classification.table)
ranger_results3$classification.table
cohen.kappa(ranger_results3$classification.table)
ranger_results4$classification.table
cohen.kappa(ranger_results4$classification.table)
ranger_results5$classification.table
cohen.kappa(ranger_results5$classification.table)




#load the main data set from Gambia

markers_for_model<-names(sort(ranger_results$variable.importance,decreasing = T))[1:3]
write.table(markers_for_model,file = "markers_for_model.txt",quote = F,sep = "\t",row.names = F,col.names = F)
#open markers_for_model file and remove allele codes at end
system(command = "plink --bfile chr12_imputed_final --extract markers_for_model.txt --recodeA --out gambian_12_for_ranger_predict")

modelling2<-fread("gambian_12_for_ranger_predict.raw",header = T,stringsAsFactors = F, showProgress = T)#pull in the data
dim(modelling2)
modelling2<-as.data.frame(modelling2)
modelling2<-modelling2[-(which(is.na(modelling2$rs2734561_1))),]
modelling2<-modelling2[-(which(is.na(modelling2$rs2246809_1))),]
modelling2<-modelling2[-(which(is.na(modelling2$rs12318583_2))),]



modelling2.prediction<-predict(object = ranger_results3,data = modelling2)
modelling2$NKG2C.prediction<-modelling2.prediction$predictions

summary(modelling2$NKG2C.prediction)
detach("package:ranger")
detach("package:Rcpp")
require(HardyWeinberg)

HWAlltests(x = c(928,1138,359))

#> HWAlltests(x = c(928,1138,359))
#Statistic   p-value
#Chi-square test:                            0.11081432 0.7392189
#Chi-square test with continuity correction: 0.09015181 0.7639843
#Likelihood-ratio test:                      0.11075740 0.7392835
#Exact test with selome p-value:                     NA 0.7311665
#Exact test with dost p-value:                       NA 0.7635623
#Exact test with mid p-value:                        NA 0.7150154
#Permutation test:                                   NA        NA

modelling2$PHENOTYPE<-factor(modelling2$PHENOTYPE)

summary(modelling2$PHENOTYPE)

ages<-read.table("EGAS00001001516_AGE.txt",header = T)
names(ages)[1]="IID"
modelling2<-merge(modelling2,ages,by = "IID")
modelling2$AGE<-as.numeric(modelling2$AGE)
fit<-glm(modelling2$PHENOTYPE~modelling2$NKG2C.prediction+modelling2$SEX+modelling2$AGE,family = "binomial")
summary(fit)
fit
exp(coef(fit))







#### test on Keneba specimens
system("head markers_for_model.txt")
#in keneba data, the markers are
#rs2246809_1 = 2:10557044:A:G
#rs2734561_1 = 12:10554985:T:C
#rs12318583_2 = 12:10579084:A:G
#create a markers file
system(command = "plink --bfile keneba_12_impute --extract markers_for_model_keneba.txt --recodeA --out keneba_data")

keneba_raw<-fread("keneba_data.raw",header = T,stringsAsFactors = F, showProgress = T)#pull in the data
dim(keneba_raw)

keneba_raw<-as.data.frame(keneba_raw)
head(keneba_raw)
names(keneba_raw)[7:9]<-c("rs2246809_1","rs2734561_1","rs12318583_2")
head(keneba_raw)
keneba_raw<-keneba_raw[-(which(is.na(keneba_raw$rs2734561_1))),]
keneba_raw<-keneba_raw[-(which(is.na(keneba_raw$rs2246809_1))),]
keneba_raw<-keneba_raw[-(which(is.na(keneba_raw$rs12318583_2))),]
keneba.prediction<-predict(object = ranger_results3,data = keneba_raw)
keneba_raw$NKG2C.prediction<-keneba.prediction$predictions


keneba_nkg2c<-read.table("001_Keneba_NKG2C_Genotypes.txt",header = T,fill = T)
keneba_nkg2c<-merge(keneba_nkg2c,keneba_raw,by.x = "cWKNO2",by.y="IID")
keneba_nkg2c<-subset(keneba_nkg2c,keneba_nkg2c$NKG2C!="BLANK")
keneba_nkg2c<-subset(keneba_nkg2c,keneba_nkg2c$NKG2C!="")
keneba_nkg2c$NKG2C<-factor(keneba_nkg2c$NKG2C)
summary(keneba_nkg2c$NKG2C)

table(keneba_nkg2c$NKG2C,keneba_nkg2c$NKG2C.prediction)


summary(modelling2$NKG2C.prediction)
detach("package:ranger")
detach("package:Rcpp")
require(HardyWeinberg)

HWAlltests(x = c(928,1138,359))
keneba_classification<-table(keneba_nkg2c$NKG2C,keneba_nkg2c$NKG2C.prediction)

cohen.kappa(keneba_classification)



