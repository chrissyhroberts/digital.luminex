library(lxb)
library(data.table)
library(mclust)
library(ggplot2)
library(mixtools)
library(dplyr)



lxbs<-list.files(path = "data_files",pattern = "*.lxb",full.names = T)


#create function to import an lxb file, add the well name, convert to a data frame and export a CSV with similar name to original file

lxb.parse<-function(lxbfile)
{
#get the file
outputfile<-gsub(x = lxbfile,pattern = ".lxb",replacement = ".csv")
well<-unlist(strsplit(lxbfile,split = "_"))
well<-well[length(well)]
well<-gsub(well,pattern = ".lxb",replacement = "")
lxbdata <- lxb::readLxb(lxbfile)  
lxbdata <- as.data.frame(lxbdata) 
names(lxbdata)<-gsub(names(lxbdata),pattern = " ",replacement = ".")
lxbdata$well<-well
#remove invalid beads
lxbdata<-lxbdata[which(lxbdata$Bead.Valid==1),]
write.csv(x = lxbdata,file = outputfile,row.names = F)
}



#define a function to aggregate the various csv files
lxb.aggregate.data<-function(path)
{
  #list only the csv files
  lxbs<-list.files(path = path,pattern = "*.csv",full.names = T)
  tables <- lapply(lxbs, fread)
  data.table::rbindlist(tables)
}


lxb.digital.thresholds<-function(df,pdf.out=FALSE,sd.threshold=3)
{
  df2<-df %>% group_by(Bead.ID) %>% summarise(n=n())
  #get bead values from DF
  df2$EM.threshold<-NA
  
  for(i in 1:length(df2$Bead.ID))
  {
    fit<-normalmixEM(df$RP1L.Peak[which(df$Bead.ID==df2$Bead.ID[i])])
    #plot(fit,which=2,breaks=40)
    
    df2$EM.threshold[i]<-fit$mu[which(fit$mu==min(fit$mu))]+(sd.threshold*fit$sigma[which(fit$mu==min(fit$mu))])
    if(pdf.out==TRUE){
                      pdf(file = paste(df2$Bead.ID[i],".pdf",sep=""))
                      par(mfrow=c(1,2))
                      plot(sort(df$RP1L.Peak[which(df$Bead.ID==df2$Bead.ID[i])]))
                      abline(h=df2$EM.threshold[i],lty=2,col="red")
                      plot(fit,whichplots = 2)
                      abline(v=df2$EM.threshold[i],lty=2,col="red")
                                        dev.off()
                      }
  }
  df<-merge(df,df2,by="Bead.ID")
  print(df2)
  df$classification<-as.numeric(df$RP1L.Peak>df$EM.threshold)  
  return(df)
}

lxb.digital.thresholds.locator<-function(df,pdf.out=FALSE,sd.threshold=3)
{
  df2<-df %>% group_by(Bead.ID) %>% summarise(n=n())
  #get bead values from DF
  df2$EM.threshold<-NA
  df2$EM.threshold<-as.numeric(df2$EM.threshold)
  
  for(i in 1:length(df2$Bead.ID))
  {
    fit<-normalmixEM(df$RP1L.Peak[which(df$Bead.ID==df2$Bead.ID[i])])
    #plot(fit,which=2,breaks=40)
    
    
      plot(fit,whichplots = 2)
      EM.threshold<-locator(n = 1)
      df2$EM.threshold[i]<-EM.threshold$x
      abline(v=df2$EM.threshold[i],lty=2,col="red")
  }  
    
  
  df<-merge(df,df2,by="Bead.ID")
  print(df2)
  df$classification<-as.numeric(df$RP1L.Peak>df$EM.threshold)  
  return(df)
}



############################################################################################################################################################
################################################################################################################################################
############################################################################################################################################################



#apply the lxb.parse function to all lxb files in directory
lapply(lxbs,lxb.parse)

#aggregate the data
a<-lxb.aggregate.data(path = "data_files")
a$col<-substr(a$well,start = 1,stop = 1)
a$row<-substr(a$well,start = 2,stop = 3)

# trim data to beads we actually used

a<-subset(a, subset = Bead.ID %in% c(12:15,18,19,20,21,22,25:30))


#plot data for all beads
ggplot(data = a,aes(x=CL1.Peak,y=CL2.Peak,color=factor(Bead.ID)))+geom_point()+
  scale_x_continuous(trans = "log")

#show table of bead counts
table(factor(a$Bead.ID),a$Bead.Valid)

#threshold classify everything based on EM algorithm
a<-lxb.digital.thresholds(a,pdf.out = T,sd.threshold = 1)

#summarise based on bead counts
aa<-a %>% group_by(Bead.ID,well) %>% summarise(count=n(),bead.pos=sum(classification==1),bead.neg=sum(classification==0),bead.fraction=(bead.pos=sum(classification==1)/n()))
aa$col<-substr(aa$well,start = 1,stop = 1)
aa$row<-substr(aa$well,start = 2,stop = 3)

#1 - standards
#2 - controls
#3 - cases
aa$welltype<-NA
aa$welltype[which(aa$row==1)]<-1
aa$welltype[which(aa$row==2)]<-1
aa$welltype[which(aa$row==3)]<-2
aa$welltype[which(aa$row==4)]<-2
aa$welltype[which(aa$row==5)]<-2
aa$welltype[which(aa$row==6)]<-2
aa$welltype[which(aa$row==7)]<-2
aa$welltype[which(aa$row==8)]<-3
aa$welltype[which(aa$row==9)]<-3
aa$welltype[which(aa$row==10)]<-3
aa$welltype[which(aa$row==11)]<-3
aa$welltype[which(aa$row==12)]<-3


bead.ag<-read.csv(file = "bead.ids.csv")

aa<-merge(aa,bead.ag)

aa<-aa[-(which(aa$welltype==1)),]
ggplot(aa,aes(factor(welltype),bead.fraction))+geom_boxplot()+facet_wrap(.~Bead.Ag,scales = "free")




##################


