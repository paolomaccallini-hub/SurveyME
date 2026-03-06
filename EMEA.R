# file: EMEA.R
#
source("EMEA_func.R")
#
#-------------------------------------------------------------------------------
# Read, edit, save edited data
#-------------------------------------------------------------------------------
#
EMEA.edit()
#
#-------------------------------------------------------------------------------
# Perform a few plots
#-------------------------------------------------------------------------------
#
EMEA.plot()
source("EMEA_Bimodal.R")
#
#-------------------------------------------------------------------------------
# Assign LoF to subjects
#-------------------------------------------------------------------------------
#
EMEA.LOF()
#
#-------------------------------------------------------------------------------
# Plot lof and PCA for custom cases
#-------------------------------------------------------------------------------
#
mycases<-c(nrow(mydata)-2,nrow(mydata)-1,nrow(mydata))
for (i in 1:length(mycases)) {
  LOF_ME(mycases[i])
}
#
#-------------------------------------------------------------------------------
# Correlations
#-------------------------------------------------------------------------------
#
mydata<-read.table(file="EMEA_edited.txt",sep="\t",header=T)
mydata<-subset.data.frame(mydata,select=-trigger_n)
mydata<-subset.data.frame(mydata,select=-course_n)
Correlations(mydata,data.select="EMEA",n.cex=10/ncol(mydata),alpha=0.05)
#
#-------------------------------------------------------------------------------
# Correlations
#-------------------------------------------------------------------------------
#
mydata<-read.table(file="EMEA_edited.txt",sep="\t",header=T)
mydata<-Expand(mydata)
Correlations(mydata,data.select="EMEA_course_triggers",n.cex=13/ncol(mydata),alpha=0.05)
#
#-------------------------------------------------------------------------------
# Comparison outliers vs rest of patients 
#-------------------------------------------------------------------------------
#
mydata<-read.table(file="EMEA_edited.txt",sep="\t",header=T)
mygroup.A<-subset.data.frame(mydata,mydata$lof<quantile(mydata$lof,probs=0.95))
mygroup.B<-subset.data.frame(mydata,mydata$lof>=quantile(mydata$lof,probs=0.95))
data.select<-c("common","outliers")
EMEA.comparison(mygroup.A,mygroup.B,data.select)
#
#-------------------------------------------------------------------------------
# Comparison females vs males 
#-------------------------------------------------------------------------------
#
mydata<-read.table(file="EMEA_edited.txt",sep="\t",header=T)
mygroup.A<-subset.data.frame(mydata,mydata$gender_n==2)
mygroup.B<-subset.data.frame(mydata,mydata$gender_n==1)
data.select<-c("females","males")
EMEA.comparison(mygroup.A,mygroup.B,data.select)
#
#-------------------------------------------------------------------------------
# Comparison affected relatives vs no affected relatives 
#-------------------------------------------------------------------------------
#
mydata<-read.table(file="EMEA_edited.txt",sep="\t",header=T)
mygroup.A<-subset.data.frame(mydata,mydata$relatives_n>0)
mygroup.B<-subset.data.frame(mydata,mydata$relatives_n==0)
data.select<-c("aff_relatives","no_aff_relatives")
EMEA.comparison(mygroup.A,mygroup.B,data.select)
#
#-------------------------------------------------------------------------------
# Comparison between two components of age at onset
#-------------------------------------------------------------------------------
#
mydata<-read.table(file="EMEA_edited.txt",sep="\t",header=T)
fit<-mixfit(mydata$age_ill+1,family="gamma",ncomp=2)
mydata$comp.prob.1<-fit$comp.prob[,1]
mydata$comp.prob.2<-fit$comp.prob[,2]
mygroup.A<-subset.data.frame(mydata,mydata$comp.prob.1>0.5&mydata$comp.prob.2<0.5)
mygroup.B<-subset.data.frame(mydata,mydata$comp.prob.1<0.5&mydata$comp.prob.2>0.5)
data.select<-c("comp_1","comp_2")
EMEA.comparison(mygroup.A,mygroup.B,data.select)
#
#-------------------------------------------------------------------------------
# Comparison between two components of age at onset for males only
#-------------------------------------------------------------------------------
#
mydata<-read.table(file="EMEA_edited.txt",sep="\t",header=T)
mydata<-subset.data.frame(mydata,gender_n==1)
fit<-mixfit(mydata$age_ill+1,family="gamma",ncomp=2)
mydata$comp.prob.1<-fit$comp.prob[,1]
mydata$comp.prob.2<-fit$comp.prob[,2]
mygroup.A<-subset.data.frame(mydata,mydata$comp.prob.1>mydata$comp.prob.2)
mygroup.B<-subset.data.frame(mydata,mydata$comp.prob.1<=mydata$comp.prob.2)
data.select<-c("M_comp_1","M_comp_2")
EMEA.comparison(mygroup.A,mygroup.B,data.select)
#
#-------------------------------------------------------------------------------
# Comparison between moderate-very severe and mild-recovered
#-------------------------------------------------------------------------------
#
mydata<-read.table(file="EMEA_edited.txt",sep="\t",header=T)
mygroup.A<-subset.data.frame(mydata,mydata$severity_n>=4)
mygroup.B<-subset.data.frame(mydata,mydata$severity_n<4)
data.select<-c("Housebound","notHousebound")
EMEA.comparison(mygroup.A,mygroup.B,data.select)
#
#-------------------------------------------------------------------------------
# Supervised classification (sensitivity, specificity) for severity and lof
#-------------------------------------------------------------------------------
#
# create a folder if absent
#
current_dir<-getwd()
folder_path<-file.path(current_dir,"Supervised")  
if(!dir.exists(folder_path)) {
  dir.create(folder_path) 
}
#
lof_cut_off<-quantile(mydata$lof,0.95)[[1]]
exit<-c("lof"=lof_cut_off,"severity_n"=3) # outcome variables and cut-off
data.select<-"Supervised"
model.type<-c("log","NB") # appropriate model among log, QDA, LDA, NB
#
# supervised classification
#
for (k in 1:length(exit)) {
  #
  mydata<-read.table(file="EMEA_edited.txt",sep="\t",header=T)
  mydata<-Expand(mydata) 
  #
  # classify 
  #
  temporary<-mydata[,names(exit)[k]]
  cut.off<-exit[k]
  for (i in 1:nrow(mydata)) {
    if (is.na(temporary[i])) {
      mydata$class[i]<-NA
    } else if (temporary[i]>cut.off) {
      mydata$class[i]<-1 # high score
    } else if (temporary[i]<=cut.off) {
      mydata$class[i]<-0 # low score
    }
  }
  #
  # remove unused variable
  #
  j<-which(colnames(mydata)==names(exit)[k])
  mydata<-subset.data.frame(mydata,select=-j)
  #
  # perform supervised classification
  #
  for (mod in model.type) {
    file.name<-paste0("Supervised/",data.select,"_",mod,"_Resampling_",names(exit)[k],".jpeg")
    jpeg(file.name,quality=100,res=100,width=1000,height=600)
    curves(mydata,model.type=mod,res=40)
    dev.off()
  }
}
#
#-------------------------------------------------------------------------------
# Supervised classification (sensitivity, specificity) for course: deteriorating
# course_n = 5 OR 7
#-------------------------------------------------------------------------------
#
current_dir<-getwd()
folder_path<-file.path(current_dir,"Supervised")  
if(!dir.exists(folder_path)) {
  dir.create(folder_path) 
}
#
exit<-c("course_n") # outcome variables and cut-off
data.select<-"Supervised"
model.type<-c("log","NB") # appropriate model among log, QDA, LDA, NB
#
# supervised classification
#
mydata<-read.table(file="EMEA_edited.txt",sep="\t",header=T)
#
# classify
#
temporary<-mydata[,exit]
cut.off<-c(5,7)
for (i in 1:nrow(mydata)) {
  if (is.na(temporary[i])) {
    mydata$class[i]<-NA
  } else if (temporary[i]%in%cut.off) {
    mydata$class[i]<-1 # high score
  } else {
    mydata$class[i]<-0 # low score
  }
}
#
# remove unused variable
#
j<-which(colnames(mydata)==exit)
mydata<-subset.data.frame(mydata,select=-j)
#
# perform supervised classification
#
for (mod in model.type) {
  file.name<-paste0("Supervised/",data.select,"_",mod,"_Resampling_",exit,".jpeg")
  jpeg(file.name,quality=100,res=100,width=1000,height=600)
  curves(mydata,model.type=mod,res=40)
  dev.off()
}
#
#-------------------------------------------------------------------------------
# Assign predicted probability to severity_n by logistic model to custom cases
#-------------------------------------------------------------------------------
#
exit<-c("severity_n"=3) # outcome variables and cut-off
data.select<-"Supervised"
model.type<-"log" # appropriate model among log, QDA, LDA, NB
#
mydata<-read.table(file="EMEA_edited.txt",sep="\t",header=T)
mydata<-Expand(mydata) 
#
# classify by severity
#
temporary<-mydata[,names(exit)]
cut.off<-exit
for (i in 1:nrow(mydata)) {
  if (is.na(temporary[i])) {
    mydata$class[i]<-NA
  } else if (temporary[i]>cut.off) {
    mydata$class[i]<-1 # high score
  } else if (temporary[i]<=cut.off) {
    mydata$class[i]<-0 # low score
  }
}
#
# remove unused variable
#
j<-which(colnames(mydata)==names(exit))
mydata<-subset.data.frame(mydata,select=-j)
mydata<-subset.data.frame(mydata,select=-lof)
#
# build model
#
model<-glm(class~.,data=mydata,family="binomial")
summary(model)
#
# make predictions    
#
pred.prob<-predict(model,newdata=mydata,type="response")
mydata$pred.prob<-pred.prob
#
# probability of being housebound for specific cases    
#
mycases<-c(nrow(mydata)-2,nrow(mydata)-1,nrow(mydata))
for (i in 1:length(mycases)) {
  print(paste(mycases[i],"is housebound with a probability of",
              round(mydata$pred.prob[mycases[i]],2)))
}
#
# ROC
#
file.name<-paste0("Supervised/ROC_",model.type,"_",names(exit),".jpeg")
jpeg(file.name,quality=100,res=100,width=1000,height=1000)
roc_obj<-pROC::roc(mydata$class,pred.prob)
plot(roc_obj,main=paste("AUC =",round(pROC::auc(roc_obj),2)))
dev.off()
#
#-------------------------------------------------------------------------------
# Assign predicted probability of being outlier by NB
#-------------------------------------------------------------------------------
#
# create a folder if absent
#
current_dir<-getwd()
folder_path<-file.path(current_dir,"Supervised")  
if(!dir.exists(folder_path)) {
  dir.create(folder_path) 
}
#
lof_cut_off<-quantile(mydata$lof,0.95)[[1]]
exit<-c("lof"=lof_cut_off) # outcome variables and cut-off
data.select<-"Supervised"
model.type<-c("NB","log") # appropriate model among log, QDA, LDA, NB
#
for (mod in model.type) {
  #
  # supervised classification
  #
  k<-1
  mydata<-read.table(file="EMEA_edited.txt",sep="\t",header=T)
  mydata<-Expand(mydata) 
  #
  # classify 
  #
  temporary<-mydata[,names(exit)[k]]
  cut.off<-exit[k]
  for (i in 1:nrow(mydata)) {
    if (is.na(temporary[i])) {
      mydata$class[i]<-NA
    } else if (temporary[i]>cut.off) {
      mydata$class[i]<-1 # high score
    } else if (temporary[i]<=cut.off) {
      mydata$class[i]<-0 # low score
    }
  }
  #
  # remove unused variable
  #
  j<-which(colnames(mydata)==names(exit)[k])
  mydata<-subset.data.frame(mydata,select=-j)
  #
  # build model
  #
  if (mod=="log") { # case of logistic regression 
    model<-glm(class~.,data=mydata,family="binomial")
  } else if (mod=="NB") { # case of Naive Bayes
    model<-e1071::naiveBayes(class~.,data=mydata)
  }
  summary(model)
  #
  # make predictions    
  #
  if (mod=="log") { # case of logistic regression 
    pred.prob<-predict(model,newdata=mydata,type="response")
  } else if (mod=="NB") { # case of naive Bayes
    pred.prob<-predict(model,mydata,type="raw")[,2]
  }
  mydata$pred.prob<-pred.prob
  #
  # probability of being outlier for specific cases    
  #
  mycases<-c(nrow(mydata)-2,nrow(mydata)-1,nrow(mydata))
  for (i in 1:length(mycases)) {
    print(paste(mycases[i],"is outlier with a probability of",
                round(mydata$pred.prob[mycases[i]],2)))
  }
  #
  # ROC
  #
  file.name<-paste0("Supervised/ROC_",mod,"_",names(exit),".jpeg")
  jpeg(file.name,quality=100,res=100,width=1000,height=1000)
  roc_obj<-pROC::roc(mydata$class,pred.prob)
  plot(roc_obj,main=paste("AUC =",round(pROC::auc(roc_obj),2)))
  dev.off()
}
#
#-------------------------------------------------------------------------------
# Probability of being deteriorating by logistic regression 
#-------------------------------------------------------------------------------
#
current_dir<-getwd()
folder_path<-file.path(current_dir,"Supervised")  
if(!dir.exists(folder_path)) {
  dir.create(folder_path) 
}
#
exit<-c("course_n") # outcome variables and cut-off
data.select<-"Supervised"
model.type<-c("log") # appropriate model among log, QDA, LDA, NB
#
# supervised classification
#
mydata<-read.table(file="EMEA_edited.txt",sep="\t",header=T)
#
# classify
#
temporary<-mydata[,exit]
cut.off<-c(5,7)
for (i in 1:nrow(mydata)) {
  if (is.na(temporary[i])) {
    mydata$class[i]<-NA
  } else if (temporary[i]%in%cut.off) {
    mydata$class[i]<-1 # high score
  } else {
    mydata$class[i]<-0 # low score
  }
}
#
# remove unused variable
#
j<-which(colnames(mydata)==exit)
mydata<-subset.data.frame(mydata,select=-j)
#
# build model
#
model<-glm(class~.,data=mydata,family="binomial")
summary(model)
#
# make predictions    
#
pred.prob<-predict(model,newdata=mydata,type="response")
mydata$pred.prob<-pred.prob
#
# probability of being deteriorating for specific cases    
#
mycases<-c(nrow(mydata)-2,nrow(mydata)-1,nrow(mydata))
for (i in 1:length(mycases)) {
  print(paste(mycases[i],"is outlier with a probability of",
              round(mydata$pred.prob[mycases[i]],2)))
}
#
# ROC
#
file.name<-paste0("Supervised/ROC_",model.type,"_",names(exit),".jpeg")
jpeg(file.name,quality=100,res=100,width=1000,height=1000)
roc_obj<-pROC::roc(mydata$class,pred.prob)
plot(roc_obj,main=paste("AUC =",round(pROC::auc(roc_obj),2)))
dev.off()
#
#-------------------------------------------------------------------------------
# Cluster analysis 
#-------------------------------------------------------------------------------
#
# create a folder if absent
#
current_dir<-getwd()
folder_path<-file.path(current_dir,"Unsupervised")  
if(!dir.exists(folder_path)) {
  dir.create(folder_path) 
}
#
# Read data
#
mydata<-fread("EMEA_edited.txt")
mydata<-as.data.frame(mydata)
#
# Clustering
#
for (c in 2:7) {
  mydata.CA<-mydata
  mydata.CA<-scale(mydata.CA)
  mydata.CA<-pam(mydata.CA,c,metric="euclidean")
  # 
  # Clusplot
  #
  file.name<-paste0("Unsupervised/PAM_",c,".jpeg")
  jpeg(file.name,quality=100,res=100,width=1300,height=500)
  par(mfrow=c(1,2))
  plot(mydata.CA,which.plot=1,main="")  
  plot(mydata.CA,which.plot=2,main="")  
  dev.off()
}
