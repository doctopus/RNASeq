args <- commandArgs(trailingOnly = TRUE)
wrkg<- args[1]
input_dir<-"/mnt/beegfs/training/CITIWorkshops/RNASeq/data/"
fc_dir<-paste0(input_dir,"subread/featureCounts/")
#setwd(paste0(wrkg,"/subread/featureCounts"))
filelst<-paste0(input_dir,"Bamfiles.txt")
p<- scan(file = filelst, what = "character")
samples=gsub(".Aligned.sortedByCoord.out.bam","",basename(p))
lst<-as.numeric(c(1:length(samples)))

fc<-c("featureCounts_0","featureCounts_s1","featureCounts_s2")
assigned<-data.frame(matrix(nrow=3,ncol=2))
colnames(assigned)<-c("file","assigned")
assigned$assigned<-0
for(f_file in 1:length(fc)){
  assigned$file[f_file]<-fc[f_file]
  for(i in 1:length(lst)){
    f<-paste(fc[f_file],".",i,".summary",sep='')
    x<-read.table(file = paste0(fc_dir,f),header = T, sep="\t",fill = T, stringsAsFactors = FALSE)
    
    assigned$assigned[f_file]<-assigned$assigned[f_file] + x[,2][which(x[,1]=="Assigned")]
  }
}

f_file<-which.max(assigned$assigned)

#### create matric for the chosen counts file
f<-paste(fc[f_file],".",lst[1],sep='')
Mat<-read.table(file = paste0(fc_dir,f),header = T, sep="\t",fill = T, stringsAsFactors = FALSE)
Idx<-c(1:dim(Mat)[1])
Mat<-data.frame(Idx,Mat)
colnames(Mat)[8]<-samples[1]
for(i in 2:length(lst)){
  f<-paste(fc[f_file],".",i,sep='')
  tmp<-read.table(file = paste0(fc_dir,f),header = T, sep="\t",fill = T, stringsAsFactors = FALSE)
  tmp<-tmp[,c(1,7)]
  colnames(tmp)[2]<-samples[i]
  
  Mat<-merge(Mat,tmp,by.x = "Geneid",by.y = "Geneid",all.x = TRUE,all.y = TRUE)
}
Mat<-Mat[order(Mat$Idx),]
write.table(Mat,file = paste0(wrkg,"/subread/featureCounts/",fc[f_file]),append = FALSE,quote = FALSE,sep="\t",row.names=FALSE,col.names = TRUE)
