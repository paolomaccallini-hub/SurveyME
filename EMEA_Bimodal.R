# file: EMEA_Bimodal.R
#
#-------------------------------------------------------------------------------
# Detection in bimodality of distributions of age at onset, stratified by sex
#-------------------------------------------------------------------------------
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
mydata.F<-subset.data.frame(mydata,gender_n==2)
mydata.M<-subset.data.frame(mydata,gender_n==1)
#
# Distributions 
#
Dist_F<-mydata.F$age_ill
summary(Dist_F)
Dist_M<-mydata.M$age_ill
summary(Dist_M)
#
Dist_F<-Dist_F+1 # avoid zeros
Dist_M<-Dist_M+1 # avoid zeros
Dist_MF<-c(Dist_F,Dist_M)
summary(Dist_MF)
#
Dist_list<-list(Dist_F,Dist_M,Dist_MF)
sexes<-c("Females"=1,"Males"=2,"Both sexes"=3)
age<-c("Younger","Older")
colv<-c("pink","cyan","gray")
xvals<-seq(0,100,1)
#
# Bimodal analysis
#
family<-c("gamma", "lnorm")
fit_list<-list()
for (fam in family) {
  print(fam)
  for (j in c(1,2,3)) {
    print(names(sexes)[j])
    fit_list[[j]]<-mixfit(Dist_list[[j]],family=fam,ncomp=2)
  }
  #
  # Plotting
  #
  file.name<-paste0("Plot/Age_at_onset_",fam,"_Bimodality.jpeg")
  jpeg(file.name,quality=100,res=300,width=2000,height=2000)
  lay.mat<-matrix(c(1,3,2,3),2,2)
  layout(mat=lay.mat)
  #
  myresults<-data.frame(group=rep(NA,6))
  counter<-0
  for (sex in c(1,2,3)) {
    #
    f_list<-list()
    if (fam=="lnorm") {
      for (gr in c(1,2)) {
        counter<-counter+1
        print(counter)
        prop<-fit_list[[sex]]$pi[gr]
        m<-fit_list[[sex]]$mulog[gr]
        s<-fit_list[[sex]]$sdlog[gr]
        myresults$group[counter]<-paste(age[gr],names(sexes[sex]))
        myresults$proportion[counter]<-round(prop,2)
        myresults$mean[counter]<-round(exp(m+0.5*s^2),2)  
        myresults$SD[counter]<-round(sqrt(exp(2*m+2*s^2)*(exp(s^2)-1)),2)
        f_list[[gr]]<-dlnorm(xvals,meanlog=m,sdlog=s)
      }
      m<-fit_list[[sex]]$pi[1]
      CDF_func<-function(q) {
        m*plnorm(q,meanlog=fit_list[[sex]]$mulog[1],sdlog=fit_list[[sex]]$sdlog[1]) +
          (1-m)*plnorm(q,meanlog=fit_list[[sex]]$mulog[2],sdlog=fit_list[[sex]]$sdlog[2])
      }
    } else if (fam=="gamma") {
      for (gr in c(1,2)) {
        counter<-counter+1
        print(counter)
        prop<-fit_list[[sex]]$pi[gr]
        a<-fit_list[[sex]]$alpha[gr]
        l<-fit_list[[sex]]$lambda[gr]
        myresults$group[counter]<-paste(age[gr],names(sexes[sex]))
        myresults$proportion[counter]<-round(prop,2)
        myresults$mean[counter]<-round(a/l,2)  
        myresults$SD[counter]<-round(sqrt(a/l^2),2)
        f_list[[gr]]<-dgamma(xvals,shape=a,rate=l)
      }
      m<-fit_list[[sex]]$pi[1]
      CDF_func<-function(q) {
        m*pgamma(q,shape=fit_list[[sex]]$alpha[1],rate=fit_list[[sex]]$lambda[1]) +
          (1-m)*pgamma(q,shape=fit_list[[sex]]$alpha[2],rate=fit_list[[sex]]$lambda[2])
      }
    }
    #
    # Kolmogorov test
    #
    a<-ks.test(Dist_list[[sex]],CDF_func) # KS test
    pval<-format(a$p.value,digits=3)
    #
    # plots
    #
    hist(Dist_list[[sex]],freq=F,xlim=c(5,90),ylim=c(0,0.04),xlab="age at onset",
         ylab="density of probability",main=paste0(names(sexes[sex])," (KS pval=",pval,")"),
         col=colv[sex])
    #
    f1<-f_list[[1]]
    f2<-f_list[[2]]
    m<-fit_list[[sex]]$pi[1]
    #
    par(new=T)
    plot(xvals,m*f1,type="l",col="black",lty=2,lwd=2,xlim=c(5,90),ylim=c(0,0.04),xlab="",ylab="")
    par(new=T)
    plot(xvals,(1-m)*f2,type="l",col="black",lty=2,lwd=2,xlim=c(5,90),ylim=c(0,0.04),xlab="",ylab="")
    par(new=T)
    plot(xvals,m*f1+(1-m)*f2,type="l",col="red",lty=1,lwd=2,xlim=c(5,90),ylim=c(0,0.04),xlab="",ylab="")
  }
  dev.off()
  #
  # Save metadata
  #
  file.name<-paste0("Plot/",fam,"_bimodal_fit.txt")
  write.table(myresults,file=file.name,sep="\t",col.names=T,row.names=F,quote=F)
}
