install.packages("lxb")
library(lxb)
library(data.table)
library(mclust)
library(ggplot2)





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


############################################################################################################################################################
################################################################################################################################################
############################################################################################################################################################



#apply the lxb.parse function to all lxb files in directory
lapply(lxbs,lxb.parse)

#aggregate the data
a<-lxb.aggregate.data(path = "data_files")


# trim data to beads we actually used

a<-subset(a, subset = Bead.ID %in% c(12:15,18,19,20,21,22,25:30))

#plot data for all beads
ggplot(data = a,aes(x=CL1.Peak,y=CL2.Peak,color=factor(Bead.ID)))+geom_point()+
  scale_x_continuous(trans = "log")

#show table of bead counts

table(factor(a$Bead.ID),a$Bead.Valid)
