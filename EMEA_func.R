# file: EMEA_func.R
#
library(readxl)
library(cluster)
library(data.table)
library(corrplot) # for corrplot
library(ellipse) # for generating ellipse
library(mixR) # study of bimodality
library(zoo)
library(pracma)
library(e1071) # for naiveBayes()
library(pROC) # for ROC curve
library(MASS) # for discriminant analysis
#
#-------------------------------------------------------------------------------
# Read, edit, save edited data
#-------------------------------------------------------------------------------
#
EMEA.edit<-function() {
  #
  # Read data
  #
  mydata<-read_excel("EMEA.xlsx",sheet=1)
  colnames(mydata)
  #
  # Select columns
  #
  # 2) Do you have ME? Binary (0-1)
  #
  variables<-("have_me_n")
  #
  # 4) Gender? 
  #
  # 1: male
  # 2: female
  # 3: prefer not to say
  #
  variables<-c(variables,"gender_n")
  #
  # 5) Year of birth? Numeric
  #
  variables<-c(variables,"year_born")
  #
  # 6) Year of onset? Numeric
  #
  variables<-c(variables,"year_ill")
  #
  # 7) Trigger? 1: infection, 2: accident/injury/surgery, 3: Traumatic life event
  #             4: vaccine, 5: pregnancy/birth, 6: none, 7:other
  #
  variables<-c(variables,"trigger_n")
  #
  # 9) Relatives with ME? Binary
  #
  variables<-c(variables,"relatives_n")
  #
  # 10) How much of these symptoms restrict your ability to participate in life? 
  #     1: not at all
  #     2: a little
  #     3: moderately
  #     4: a lot
  #     5: very much
  #     6: not relevant
  #
  variables<-c(variables,"fatigue_n","muscleorjointpain_n","dizzness_n",
               "sleepproblems_n","sensitivity_n","pem_n")
  #
  # 11) Course of the illness since onset 
  #
  #     1: mostly stable, with only small fluctuations
  #     2: major fluctuations throughout 
  #     3: flucuating initially, then mostly stable
  #     4: flucuating initially, then mostly improvement
  #     5: flucuating initially, then mostly deterioration
  #     6: mainly improvement
  #     7: mainly deterioration
  #
  variables<-c(variables,"course_n")
  #
  # 13) Comorbidities 
  #
  variables<-c(variables,"fibromyalgia_n","Hashimotosthyroiditis_n",
               "interstitialcystitis_n","siccasyndrome_n","irritablebowelsyndrome_n",
               "migraine_n","allergies_n","endometriosis_n","asthma_n",
               "multiplechemicalsensitivities_n")
  #
  # 14) Severity of the illness 
  #
  #     1: I am completely recovered
  #     2: I am better than "mild", but not completely recovered
  #     3: Mild: at least an approximate 50% reduction in pre illness activity level
  #     4: Moderate: mostly housebound
  #     5: Severe: mostly bedridden
  #     6: Very severe: totally bedridden and need help with basic functions
  #
  variables<-c(variables,"severity_n")
  #
  # Clean data
  #
  mydata<-subset.data.frame(mydata,select=variables)
  print(nrow(mydata))
  mydata<-subset.data.frame(mydata,have_me_n==1) # only subjects with ME diagnosis
  print(nrow(mydata))
  mydata<-subset.data.frame(mydata,select=-have_me_n) # not necessary anymore
  print(nrow(mydata))
  mydata<-subset.data.frame(mydata,gender_n!=3) # only subjects with specified sex
  print(nrow(mydata))
  mydata<-subset.data.frame(mydata,year_ill>=year_born) # patients must get sick after the have been born
  print(nrow(mydata))
  #
  # Remove male patients with trigger = pregnancy/birth with year_ill>year_born
  #
  mydata<-subset.data.frame(mydata,!(mydata$trigger_n==5&mydata$gender_n==1&(mydata$year_ill>mydata$year_born))) # patients must get sick after the have been born
  print(nrow(mydata))
  #
  # Edit data on comorbidities
  #
  variables<-c("fibromyalgia_n","Hashimotosthyroiditis_n",
               "interstitialcystitis_n","siccasyndrome_n","irritablebowelsyndrome_n",
               "migraine_n","allergies_n","endometriosis_n","asthma_n",
               "multiplechemicalsensitivities_n")
  #
  for (column in variables) {
    j<-which(colnames(mydata)==column)
    for (i in 1:nrow(mydata)) {
      if (!is.na(mydata[i,j])) {
        mydata[i,j]<-1
      } else {
        mydata[i,j]<-0
      }
    }
  }
  #
  # Edit data on symptoms
  #
  variables<-c("fatigue_n","muscleorjointpain_n","dizzness_n","sleepproblems_n",
               "sensitivity_n","pem_n")
  #
  for (column in variables) {
    j<-which(colnames(mydata)==column)
    for (i in 1:nrow(mydata)) {
      if (!is.na(mydata[i,j])) {
        if (mydata[i,j]==6) mydata[i,j]<-NA
      }
    }
  }
  #
  # Remove patients with missing data
  #
  mydata<-na.omit(mydata)
  print(nrow(mydata))
  #
  # Correct males with endometriosis
  #
  index<-which(mydata$endometriosis_n==1&mydata$gender_n==1)
  mydata$endometriosis_n[index]<-0
  #
  # Add age ill
  #
  mydata$age_ill<-mydata$year_ill-mydata$year_born
  mydata<-subset.data.frame(mydata,age_ill>-2)
  nrow(mydata)
  #
  # Add age 
  #
  mydata$age<-2021-mydata$year_born
  #
  # Add disease duration 
  #
  mydata$duration<-2021-mydata$year_ill
  #
  # Remove unused columns
  #
  mydata<-subset.data.frame(mydata,select=-year_born)
  mydata<-subset.data.frame(mydata,select=-year_ill)
  #
  # Add available custom cases
  #
  test<-0
  i<-1
  while(test==0) {
    file.name<-file.path("CustomCases",paste0("PT",i,".txt")) 
    if(file.exists(file.name)) {
      df<-read.table(file=file.name,sep="\t",header=T)
      rownames(df)<-paste0("PT",i)
      mydata<-rbind(mydata,df)
      i<-i+1
    } else {
      test<-1
    }
  }
  print(nrow(mydata))
  #
  # Save edited data
  #
  write.table(mydata,file="EMEA_edited.txt",sep="\t",col.names=T,row.names=T)
}
#
#-------------------------------------------------------------------------------
# Perform a few plots
#-------------------------------------------------------------------------------
#
EMEA.plot<-function() {
  #
  # create a folder if absent
  #
  current_dir<-getwd()
  folder_path<-file.path(current_dir,"Plot")  
  if(!dir.exists(folder_path)) {
    dir.create(folder_path) 
  }
  #
  # Read edited data
  #
  mydata<-read.table(file="EMEA_edited.txt",sep="\t",header=T)
  #
  # Plot distribution of age at first symptoms
  #
  mydata.F<-subset.data.frame(mydata,gender_n==2)
  mydata.M<-subset.data.frame(mydata,gender_n==1)
  #
  # Females and males
  #
  file.name<-paste0("Plot/Age_at_onset_Hist.jpeg")
  jpeg(file.name,quality=100,res=300,width=2000,height=2000)
  lay.mat<-matrix(c(1,3,2,3),2,2)
  layout(mat=lay.mat)
  hist(mydata.F$age_ill,freq=T,xlab="aget at first symptoms (F)",main="",col="pink")
  abline(v=mean(mydata.F$age_ill),col="red",lwd=2)
  hist(mydata.M$age_ill,freq=T,xlab="aget at first symptoms (M)",main="",col="blue")
  abline(v=mean(mydata.M$age_ill),col="red",lwd=2)
  hist(mydata$age_ill,freq=T,xlab="age at first symptoms (M+F)",main="")
  abline(v=mean(mydata$age_ill),col="red",lwd=2)
  dev.off()
}
#
#-------------------------------------------------------------------------------
# Assign LoF to subjects
#-------------------------------------------------------------------------------
#
EMEA.LOF<-function() {
  #
  mydata<-read.table(file="EMEA_edited.txt",sep="\t",header=T)
  #
  # Eliminate columns that won't be used in PCA
  #
  if ("lof"%in%colnames(mydata)) mydata<-subset.data.frame(mydata,select=-lof)
  #
  # Local Outlier Factor
  #
  if(!require("dbscan")) {
    install.packages("dbscan")
    library(dbscan)
  }
  minPts.LB<-10 # lower bound for minPts (minimum size of a cluster of patients)
  minPts.UB<-round(nrow(mydata)*0.05) # upper bound of minPts (maximum size of a cluster of outliers)
  lof<-matrix(nrow=(nrow(mydata)),ncol=(minPts.UB-minPts.LB+1))
  for (k in minPts.LB:minPts.UB) {
    lof[,(k-minPts.LB+1)]<-lof(as.matrix(mydata),minPts=k) # calculate lof for minPts=k
  }
  lof.result<-c()
  for (i in 1:(nrow(mydata))) {
    lof.result[i]<-max(lof[i,]) # algorithm by Breunig, Kriegel, Ng, Sander 2000
  }
  #
  mydata$lof<-lof.result
  #
  # Save edited data
  #
  write.table(mydata,file="EMEA_edited.txt",sep="\t",col.names=T,row.names=F)
}
#
#-------------------------------------------------------------------------------
# Expand columns trigger_n and course_n
#-------------------------------------------------------------------------------
#
Expand<-function(mydata) {
  #
  # Edit trigger_n
  #
  for (i in 1:nrow(mydata)) {
    if (mydata$trigger_n[i]==1) {
      mydata$infection_n[i]<-1
    } else {
      mydata$infection_n[i]<-0
    }
    if (mydata$trigger_n[i]==2) {
      mydata$AIS_n[i]<-1
    } else {
      mydata$AIS_n[i]<-0
    }
    if (mydata$trigger_n[i]==3) {
      mydata$trauma_n[i]<-1
    } else {
      mydata$trauma_n[i]<-0
    }
    if (mydata$trigger_n[i]==4) {
      mydata$vaccine_n[i]<-1
    } else {
      mydata$vaccine_n[i]<-0
    }
    if (mydata$trigger_n[i]==5) {
      mydata$pregnancy_n[i]<-1
    } else {
      mydata$pregnancy_n[i]<-0
    }
    if (mydata$trigger_n[i]==6) {
      mydata$no_trigger_n[i]<-1
    } else {
      mydata$no_trigger_n[i]<-0
    }
    if (mydata$trigger_n[i]==6) {
      mydata$other_n[i]<-1
    } else {
      mydata$other_n[i]<-0
    }
  }
  mydata<-subset.data.frame(mydata,select=-trigger_n)
  #
  # Edit course_n
  #
  for (i in 1:nrow(mydata)) {
    if (mydata$course_n[i]%in%c(1,3)) {
      mydata$stable_n[i]<-1
    } else {
      mydata$stable_n[i]<-0
    }
    if (mydata$course_n[i]%in%c(2)) {
      mydata$fluctuating_n[i]<-1
    } else {
      mydata$fluctuating_n[i]<-0
    }
    if (mydata$course_n[i]%in%c(4,6)) {
      mydata$improving_n[i]<-1
    } else {
      mydata$improving_n[i]<-0
    }
    if (mydata$course_n[i]%in%c(5,7)) {
      mydata$deteriorating_n[i]<-1
    } else {
      mydata$deteriorating_n[i]<-0
    }
  }
  mydata<-subset.data.frame(mydata,select=-course_n)
  #
  return(mydata)
}
#
#-------------------------------------------------------------------------------
# Comparisons of two subgroups
#-------------------------------------------------------------------------------
#
EMEA.comparison<-function(mygroup.A,mygroup.B,data.select) {
  #
  # create a folder if absent
  #
  current_dir<-getwd()
  folder_path<-file.path(current_dir,"Comparison")  
  if(!dir.exists(folder_path)) {
    dir.create(folder_path) 
  }
  #
  # Read and edit trigger_n and course_n
  #
  mygroup.A<-Expand(mygroup.A)
  mygroup.B<-Expand(mygroup.B)
  mydata<-Expand(mydata)
  #
  # Edit sex as 0-1 value
  #
  mygroup.A$gender_n<-mygroup.A$gender_n-1
  mygroup.B$gender_n<-mygroup.B$gender_n-1
  #
  # Build a matrix for metadata and statistical testing
  #
  mymetadata<-matrix(data=NA,nrow=ncol(mygroup.B),ncol=5)
  colnames(mymetadata)<-c("variables",
                          paste0(data.select[1]," (",nrow(mygroup.A),")"),
                          paste0(data.select[2]," (",nrow(mygroup.B),")"),
                          "p.val","p.adjust")
  for (j in 1:nrow(mymetadata)) {
    #
    mymetadata[j,1]<-colnames(mydata)[j]
    #
    if (length(unique(mydata[,j]))==2) { # binary data
      #
      A0<-length(which(mygroup.A[,j]==0))
      A1<-length(which(mygroup.A[,j]==1))
      mymetadata[j,2]<-paste0(round(100*A1/(A1+A0),2)," %")
      B0<-length(which(mygroup.B[,j]==0))
      B1<-length(which(mygroup.B[,j]==1))
      mymetadata[j,3]<-paste0(round(100*B1/(B1+B0),2)," %")
      #
      test<-fisher.test(rbind(c(A0,A1),c(B0,B1)))
      mymetadata[j,4]<-formatC(test$p.value,format="e",digits=2)
      #
    } else { #non-binary data
      #
      A<-mygroup.A[,j]
      B<-mygroup.B[,j]
      #
      m<-round(mean(A,na.rm=T),2)
      sd<-round(sd(A,na.rm=T),2)
      mymetadata[j,2]<-paste0(m," (",sd,")")
      #
      m<-round(mean(B,na.rm=T),2)
      sd<-round(sd(B,na.rm=T),2)
      mymetadata[j,3]<-paste0(m," (",sd,")")
      #
      test<-t.test(A,B,alternative="two.sided")
      mymetadata[j,4]<-formatC(test$p.value,format="e",digits=2)
      #
    }
  }
  mymetadata[,5]<-formatC(p.adjust(mymetadata[,4],method="BH"),format="e",digits=2)
  mymetadata<-mymetadata[order(as.numeric(mymetadata[,5])),]
  #
  # Save metadata
  #
  file.name<-paste0("Comparison/",data.select[1],"_vs_",data.select[2],"_comparison.tsv")
  write.table(mymetadata,file=file.name,sep="\t",col.names=T,
              row.names=F,quote=F)
  #
  return(mymetadata)
}
#
#-------------------------------------------------------------------------------
# It calculates and plot the correlation table 
#-------------------------------------------------------------------------------
#
Correlations<-function(mydata,data.select,n.cex,alpha) {
  #
  # create a folder if absent
  #
  current_dir<-getwd()
  folder_path<-file.path(current_dir,"Correlations")  
  if(!dir.exists(folder_path)) {
    dir.create(folder_path) 
  }
  #
  NC<-ncol(mydata)
  file.name<-paste0("Correlations/",data.select,"_corr.jpeg")
  jpeg(file.name,quality=100,res=100,width=1000,height=1000)
  p_mat<-matrix(nrow=NC,ncol=NC) # build a matrix with p-values for correlations
  for (i in 1:NC) {
    for (j in 1:NC) {
      temp<-cor.test(mydata[,i],mydata[,j],alternative="two.sided",method="spearman")
      p_mat[i,j]<-temp$p.value
    }
  }
  rownames(p_mat)<-colnames(mydata[,c(1:NC)])
  colnames(p_mat)<-colnames(mydata[,c(1:NC)])
  corrplot(cor(mydata[,1:NC],use="pairwise.complete.obs",method="pearson"),is.corr=F, # plot the correlation table
           method="color",addCoef.col="black",number.digits=2,number.cex=n.cex,
           mar=c(0,0,4,0),p.mat=p_mat,sig.level=alpha,insig="blank")
  dev.off()    
  #
}
#
#-------------------------------------------------------------------------------
# Plot lof and PCA for a new case
#-------------------------------------------------------------------------------
#
LOF_ME<-function(mycase) {
  #
  # create a folder if absent
  #
  current_dir<-getwd()
  folder_path<-file.path(current_dir,"LOF")  
  if(!dir.exists(folder_path)) {
    dir.create(folder_path) 
  }
  #
  mydata<-fread(file="EMEA_edited.txt",sep="\t",header=T)
  #
  # scale and PCA
  #
  NC<-ncol(mydata)
  mydata.pca<-scale(mydata) # scale data 
  mydata.pca<-princomp(mydata.pca) 
  Var1<-mydata.pca$sdev[1]^2 
  Var2<-mydata.pca$sdev[2]^2 
  Var3<-mydata.pca$sdev[3]^2
  VarTot<-sum(mydata.pca$sdev^2) 
  Var1str<-as.character(round(Var1*100/VarTot,1)) 
  Var2str<-as.character(round(Var2*100/VarTot,1))
  Var3str<-as.character(round(Var3*100/VarTot,1))
  #
  # save loadings
  # 
  loadings<-round(mydata.pca$loadings,4)
  colnames(loadings)<-seq(1,ncol(mydata))
  rownames(loadings)<-colnames(mydata)
  write.csv(loadings,file="LOF/loadings.csv",quote=F,row.names=T)
  #
  # PCA plot
  #
  N.R<-nrow(mydata)
  myrow<-which(rownames(mydata)==mycase)
  if (T) {
    random<-data.frame(PC1=rep(NA,N.R),PC2=rep(NA,N.R),PC3=rep(NA,N.R))
    col.vec<-c()
    pch.vec<-c()
    h<-0
    for (i in 1:N.R) {
      if (rownames(mydata)[i]!=mycase) {
        col.vec[i]<-"blue" 
        pch.vec[i]<-1
      } else {
        col.vec[i]<-"red" 
        pch.vec[i]<-16
      }
      random$PC1[i]<-mydata.pca$scores[i,1]
      random$PC2[i]<-mydata.pca$scores[i,2]
      random$PC3[i]<-mydata.pca$scores[i,3]
    }
    file.name<-paste0("LOF/PCA_plots_",mycase,".jpeg")
    jpeg(file.name,quality=100,res=300,width=2000,height=2000)
    par(mfrow=c(2,2))
    #
    # plot plane PC1-PC2
    #
    multy<-5 # multiplier for axes
    minX1<-2*min(mydata.pca$scores[,1])
    maxX1<-2*max(mydata.pca$scores[,1])
    minX2<-2*min(mydata.pca$scores[,2])
    maxX2<-2*max(mydata.pca$scores[,2])
    plot(mydata.pca$scores[,1],mydata.pca$scores[,2],pch=pch.vec,
         xlim=c(minX1,maxX1),ylim=c(minX2,maxX2),col=col.vec,cex=1.5,
         xlab=paste("PC1 - Var =",Var1str,"%"),ylab=paste("PC2 - Var =",Var2str,"%"),asp=1)
    ellipse_random<-ellipse(cov(random[,c(1,2)]),centre=colMeans(random[,c(1,2)]),level=0.95)
    lines(ellipse_random[,1],ellipse_random[,2],col="blue",lwd=3)
    #
    # plot plane PC1-PC3
    #
    minX3<-2*min(mydata.pca$scores[,3])
    maxX3<-2*max(mydata.pca$scores[,3])
    plot(mydata.pca$scores[,1],mydata.pca$scores[,3],pch=pch.vec,
         xlim=c(minX1,maxX1),ylim=c(minX3,maxX3),col=col.vec,cex=1.5,
         xlab=paste("PC1 - Var =",Var1str,"%"),ylab=paste("PC3 - Var =",Var3str,"%"),asp=1)
    ellipse_random<-ellipse(cov(random[,c(1,3)]),centre=colMeans(random[,c(1,3)]),level=0.95)
    lines(ellipse_random[,1],ellipse_random[,2],col="blue",lwd=3)
    #
    # plot plane PC2-PC3
    #
    plot(mydata.pca$scores[,2],mydata.pca$scores[,3],pch=pch.vec,
         xlim=c(minX2,maxX2),ylim=c(minX3,maxX3),col=col.vec,cex=1.5,
         xlab=paste("PC2 - Var =",Var2str,"%"),ylab=paste("PC3 - Var =",Var3str,"%"),asp=1)
    ellipse_random<-ellipse(cov(random[,c(2,3)]),centre=colMeans(random[,c(2,3)]),level=0.95)
    lines(ellipse_random[,1],ellipse_random[,2],col="blue",lwd=3)
    #
    # 3D plot
    #
    library(scatterplot3d)
    scatterplot3d(mydata.pca$scores[,c(1:3)],pch=pch.vec,color=col.vec,type="p",angle=30)
    dev.off()
  }
  #
  # alternative PCA plot
  #
  if (T) {
    file.name<-paste0("LOF/PCA_biplots_",mycase,".jpeg")
    jpeg(file.name,quality=100,res=300,width=2000,height=2000)
    par(mfrow=c(2,2))
    #
    biplot(mydata.pca,cex=1,col=c("blue"),choices=c(1,2),xlabs=rep("",N.R),
           xlab=paste("PC1 - Var =",Var1str,"%"), 
           ylab=paste("PC2 - Var =",Var2str,"%")) 
    points(mydata.pca$scores[1:N.R,1],
           mydata.pca$scores[1:N.R,2],pch=21,col="black",cex=1)
    points(mydata.pca$scores[myrow,1],
           mydata.pca$scores[myrow,2],pch=20,col="red",cex=2)
    ellipse_random<-ellipse(cov(random[,c(1,2)]),centre=colMeans(random[,c(1,2)]),level=0.95)
    lines(ellipse_random[,1],ellipse_random[,2],col="black",lwd=1.5)
    #
    biplot(mydata.pca,col=c("blue"),choices=c(1,3),xlabs=rep("",N.R),
           xlab=paste("PC1 - Var =",Var1str,"%"), 
           ylab=paste("PC3 - Var =",Var3str,"%"))  
    points(mydata.pca$scores[1:N.R,1],
           mydata.pca$scores[1:N.R,3],pch=21,col="black",cex=1)
    points(mydata.pca$scores[myrow,1],
           mydata.pca$scores[myrow,3],pch=20,col="red",cex=2)#
    ellipse_random<-ellipse(cov(random[,c(1,3)]),centre=colMeans(random[,c(1,3)]),level=0.95)
    lines(ellipse_random[,1],ellipse_random[,2],col="black",lwd=1.5)
    #
    biplot(mydata.pca,col=c("blue"),choices=c(2,3),xlabs=rep("",N.R),
           xlab=paste("PC2 - Var =",Var2str,"%"), 
           ylab=paste("PC3 - Var =",Var3str,"%"))
    points(mydata.pca$scores[1:N.R,2],
           mydata.pca$scores[1:N.R,3],pch=21,col="black",cex=1)
    points(mydata.pca$scores[myrow,2],
           mydata.pca$scores[myrow,3],pch=20,col="red",cex=2)
    ellipse_random<-ellipse(cov(random[,c(2,3)]),centre=colMeans(random[,c(2,3)]),level=0.95)
    lines(ellipse_random[,1],ellipse_random[,2],col="black",lwd=1.5)
    #
    plot(cumsum(mydata.pca$sdev^2)/sum(mydata.pca$sdev^2),type="b",ylim=c(0,1),
         xlab=("number of components"),
         ylab=("cumulative variance/total variance")) 
    segments(1,0.8,NC,0.8,col="red")
    dev.off()
  }
  #
  file.name<-paste0("LOF/LOF_histogram_",mycase,".jpeg")
  jpeg(file.name,quality=100,res=300,width=2000,height=1000)
  lof.result<-mydata$lof
  cutoff<-quantile(lof.result,0.95)  # Top 5% most anomalous points
  outliers<-which(lof.result>=lof.result[myrow])
  pval<-length(outliers)/nrow(mydata)
  pval<-format(pval,digit=2)
  par(mfrow=c(1,1))
  hist(lof.result,breaks=40,xlim=c(min(lof.result),3),
       main=paste0("Local Outlier Factor = ",round(lof.result[myrow],2)," (pval = ",pval,")")) 
  abline(v=lof.result[myrow],col="red",lwd=2) # lof of ME/CFS
  dev.off()
}
#
#-------------------------------------------------------------------------------
# Function for sensitivity, specificity, and curve FPR-TPR
#-------------------------------------------------------------------------------
#
curves<-function(mydata,model.type,res) {
  #
  # curves of specificity, sensitivity, accuracy, and error rates (type I and II)
  #
  PP<-seq(from=0.01,to=0.99,by=0.05)
  sensitivity<-c()
  specificity<-c()
  accuracy<-c()
  for (k in 1:length(PP)) {
    print(paste("working on PP",PP[k]))
    results<-mean.sens.spec(mydata,model.type,PP[k],times=res) 
    sensitivity[k]<-results[1]
    specificity[k]<-results[2]
    accuracy[k]<-results[3]
  }
  #
  indices<-seq(0,1,length.out=length(PP))
  sensitivity<-approx(indices,sensitivity,xout=seq(0,1,by=0.001))$y
  specificity<-approx(indices,specificity,xout=seq(0,1,by=0.001))$y
  accuracy<-approx(indices,accuracy,xout=seq(0,1,by = 0.001))$y
  PP<-approx(indices,PP,xout=seq(0,1,by=0.001))$y
  ERI<-1-specificity # Error Rate Type I
  ERII<-1-sensitivity # Error Rate Type II
  par(mfrow=c(1,2))
  title.str<-paste("Model type =", model.type)
  plot(PP,ERI,type="l",lty=1,col="red",xlab="predicted probability",ylab="",
       xlim=c(0.02,1.3),ylim=c(0,1),lwd=2)
  par(new=T)
  plot(PP,sensitivity,type="l",lty=2,col="red",xlab="",ylab="",
       xlim=c(0.02,1.3),ylim=c(0,1),lwd=2)
  par(new=T)
  plot(PP,ERII,type="l",lty=1,col="blue",xlab="",ylab="",
       xlim=c(0.02,1.3),ylim=c(0,1),lwd=2)
  par(new=T)
  plot(PP,specificity,type="l",lty=2,col="blue",xlab="",ylab="",
       xlim=c(0.02,1.3),ylim=c(0,1),lwd=2)
  par(new=T)
  plot(PP,accuracy,type="l",lty=2,col="green",xlab="",ylab="",
       xlim=c(0.02,1.3),ylim=c(0,1),lwd=2)
  abline(h=seq(0,1,0.1),v=seq(0,1,0.1),lty=2)
  legend("topright",legend=c("ERI","sens","ERII","spec","ACC"),col=c("red","red","blue"
                                                                     ,"blue","green"),lty=c(1,2,1,2,2),lwd=c(2,2,2,2,2))
  title(main=title.str)
  #
  # curve
  #
  ERI[1]<-1
  ERI[length(ERII)]<-0
  sensitivity[1]<-1
  sensitivity[length(sensitivity)]<-0
  ERI<-na.approx(ERI) # interpolate NAs
  AUC<-trapz(ERI,sensitivity)
  plot(ERI,sensitivity,type="l",lwd=2,col=1,xlab="False positive rate (1 - Specificity)",
       ylab="True positive rate (Sensitivity)",xlim=c(0,1),ylim=c(0,1))
  abline(h=seq(0,1,0.1),v=seq(0,1,0.1),lty=2)
  par(new=T)
  plot(c(0,1),c(0,1),col="black",type="l",lty=2,xlim=c(0,1),ylim=c(0,1),xlab="",ylab="",lwd=2)
  title(main=paste("AUC=",round(abs(AUC),2)))
}
#
#-------------------------------------------------------------------------------
# Function for mean sensitivity and specificity using a subset for training and 
# remaining data for predictions, repeating the algorithm several times
#-------------------------------------------------------------------------------
#
mean.sens.spec<-function(mydata,model.type,PP,times) {
  #
  mydata<-na.omit(mydata)
  K<-round(nrow(mydata)/2,0)-1
  sensitivity<-rep(NA,times)
  specificity<-rep(NA,times)
  accuracy<-rep(NA,times)
  #
  for (j in 1:times) {
    # print(paste("Attempt",j))
    pred.rows<-sample(nrow(mydata),K) # randomly select K rows 
    sub.mydata<-mydata[-pred.rows,] # subset of data (remove K rows)
    #
    # Regression on subset of data
    #
    if (model.type=="log") { # case of logistic regression 
      sub.model<-glm(class~.,data=sub.mydata,family="binomial")
    } else if (model.type=="LDA") { # case of LDA
      sub.model<-lda(class~.,data=sub.mydata)
    } else if (model.type=="QDA") { # case of QDA
      sub.model<-qda(class~.,data=sub.mydata)
    } else if (model.type=="NB") { # case of Naive Bayes
      sub.model<-e1071::naiveBayes(class~.,data=sub.mydata)
    }
    #
    # Make predictions on removed data based on training data-set   
    #
    pred.class<-c()
    if (model.type=="log") { # case of logistic regression 
      title.str<-paste("Log Reg,","cut off =",PP)
      pred.class<-predict(sub.model,newdata=mydata[pred.rows,],type="response")
      index<-which(is.na(pred.class))
      if(length(index)>0) pred.class[index]<-3
      index<-which(pred.class>=PP)
      if(length(index)>0) pred.class[index]<-1
      index<-which(pred.class<PP)
      if(length(index)>0) pred.class[index]<-0
    } else if (model.type=="LDA"|model.type=="QDA") { # case of DA
      title.str<-paste("QDA, PP=",PP)
      if (model.type=="LDA") title.str<-paste("LDA, PP=",PP)
      pred<-predict(sub.model,mydata[pred.rows,])
      index<-which(is.na(pred$posterior[,2]))
      if(length(index)>0) pred.class[index]<-3
      index<-which(pred$posterior[,2]>=PP)
      if(length(index)>0) pred.class[index]<-1
      index<-which(pred$posterior[,2]<PP)
      if(length(index)>0) pred.class[index]<-0
    } else if (model.type=="NB") { # case of naive Bayes
      title.str<-paste("Naive Bayes, PP=",PP)
      pred<-predict(sub.model,mydata[pred.rows,],type="raw")
      index<-which(is.na(pred[,2]))
      if(length(index)>0) pred.class[index]<-3
      index<-which(pred[,2]>=PP)
      if(length(index)>0) pred.class[index]<-1
      index<-which(pred[,2]<PP)
      if(length(index)>0) pred.class[index]<-0
    }
    #
    df<-data.frame(class=mydata$class[pred.rows],pred.class)
    CT<-table(df)
    if (ncol(CT)==2&nrow(CT)==2) {
      sensitivity[j]<-CT[2,2]/(CT[2,2]+CT[2,1])
      specificity[j]<-CT[1,1]/(CT[1,1]+CT[1,2])  
    }
    nc<-sum(pred.class==mydata[pred.rows,]$class) # number of correct predictions
    accuracy[j]<-nc/K
  }
  #
  result<-c(mean(sensitivity,na.rm=T),mean(specificity,na.rm=T),mean(accuracy,na.rm=T))
  if (F) {
    par(mfrow=c(1,3))
    hist(sensitivity,breaks=10,freq=F,xlab="sensitivity (%)",main="")
    abline(v=mean(sensitivity,na.rm=T),col="red",lwd=2)
    title(main=title.str)
    hist(specificity,breaks=10,freq=F,xlab="specificity (%)",main="")
    abline(v=mean(specificity,na.rm=T),col="red",lwd=2)
    hist(accuracy,breaks=10,freq=F,xlab="accuracy (%)",main="")
    abline(v=mean(accuracy,na.rm=T),col="red",lwd=2)  
  }
  return(result)
  #
}