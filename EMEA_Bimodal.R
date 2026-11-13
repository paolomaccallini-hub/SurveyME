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
criteria_results<-data.frame() # goodness of fit tests
#
# Bimodal analysis
#
family<-c("gamma", "lnorm", "norm", "weibull")
fit_list<-list()
for (fam in family) {
  print(fam)
  for (j in c(1,2,3)) {
    print(names(sexes)[j])
    fit_list[[j]]<-mixfit(Dist_list[[j]],family=fam,ncomp=2)
  }
  #
  # Plotting empirical distributions and models
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
      density_func<-function(q) {
        m*dlnorm(q,meanlog=fit_list[[sex]]$mulog[1],sdlog=fit_list[[sex]]$sdlog[1]) +
          (1-m)*dlnorm(q,meanlog=fit_list[[sex]]$mulog[2],sdlog=fit_list[[sex]]$sdlog[2])
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
      density_func<-function(q) {
        m*dgamma(q,shape=fit_list[[sex]]$alpha[1],rate=fit_list[[sex]]$lambda[1]) +
          (1-m)*dgamma(q,shape=fit_list[[sex]]$alpha[2],rate=fit_list[[sex]]$lambda[2])
      }
    } else if (fam=="norm") {
      for (gr in c(1,2)) {
        counter<-counter+1
        print(counter)
        prop<-fit_list[[sex]]$pi[gr]
        mu<-fit_list[[sex]]$mu[gr]
        sd<-fit_list[[sex]]$sd[gr]
        myresults$group[counter]<-paste(age[gr],names(sexes[sex]))
        myresults$proportion[counter]<-round(prop,2)
        myresults$mean[counter]<-mu 
        myresults$SD[counter]<-sd
        f_list[[gr]]<-dnorm(xvals,mean=mu,sd=sd) 
      }
      m<-fit_list[[sex]]$pi[1]
      CDF_func<-function(q) {
        m*pnorm(q,mean=fit_list[[sex]]$mu[1],sd=fit_list[[sex]]$sd[1]) +
          (1-m)*pnorm(q,mean=fit_list[[sex]]$mu[2],sd=fit_list[[sex]]$sd[2])
      }
      density_func<-function(q) {
        m*dnorm(q,mean=fit_list[[sex]]$mu[1],sd=fit_list[[sex]]$sd[1],) +
          (1-m)*dnorm(q,mean=fit_list[[sex]]$mu[2],sd=fit_list[[sex]]$sd[2])
      }
    } else if (fam=="weibull") {
      for (gr in c(1,2)) {
        counter<-counter+1
        print(counter)
        prop<-fit_list[[sex]]$pi[gr]
        k<-fit_list[[sex]]$k[gr]
        lambda<-fit_list[[sex]]$lambda[gr]
        myresults$group[counter]<-paste(age[gr],names(sexes[sex]))
        myresults$proportion[counter]<-round(prop,2)
        myresults$mean[counter]<-round(lambda*gamma(1+1/k),2)
        myresults$SD[counter]<-round(lambda*sqrt(gamma(1+2/k)-gamma(1+1/k)^2),2)
        f_list[[gr]]<-dweibull(xvals,shape=k,scale=lambda)
      }
      m<-fit_list[[sex]]$pi[1]
      CDF_func<-function(q) {
        m*pweibull(q,shape=fit_list[[sex]]$k[1],scale=fit_list[[sex]]$lambda[1]) +
          (1-m)*pweibull(q,shape=fit_list[[sex]]$k[2],scale=fit_list[[sex]]$lambda[2])
      }
      density_func<-function(q) {
        m*dweibull(q,shape=fit_list[[sex]]$k[1],scale=fit_list[[sex]]$lambda[1]) +
          (1-m)*dweibull(q,shape=fit_list[[sex]]$k[2],scale=fit_list[[sex]]$lambda[2])
      }
    }
    #
    # Log-likelihood, AIC and BIC
    #
    fit<-fit_list[[sex]]
    #
    loglik<-sum(log(pmax(density_func(Dist_list[[sex]]),.Machine$double.xmin)))
    n<-length(Dist_list[[sex]])
    npar<-5
    AIC_value<-2*npar-2*loglik
    BIC_value<-npar*log(n)-2*loglik
    #
    # Kolmogorov test
    #
    a<-ks.test(Dist_list[[sex]],CDF_func) # KS test
    pval<-format(a$p.value,digits=3)
    #
    # Store results for goodness of fit
    #
    df<-data.frame(family=fam,group=names(sexes)[sex],n=n,logLik=loglik,parameters=npar,
      AIC=AIC_value,BIC=BIC_value,KS_D=unname(a$statistic),KS_p=a$p.value)
    criteria_results<-rbind(criteria_results,df)
    #
    # plots
    #
    hist(Dist_list[[sex]],breaks=seq(0,100,2),freq=F,xlim=c(5,90),ylim=c(0,0.04),xlab="age at onset",
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
  # Plotting QQplots: empirical distributions against models
  #
  file.name<-paste0("Plot/Age_at_onset_",fam,"_QQ.jpeg")
  jpeg(file.name,quality=100,res=300,width=3000,height=1000)
  lay.mat<-matrix(c(1,2,3),1,3)
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
        prop<-fit_list[[sex]]$pi[gr]
        m<-fit_list[[sex]]$mulog[gr]
        s<-fit_list[[sex]]$sdlog[gr]
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
        prop<-fit_list[[sex]]$pi[gr]
        a<-fit_list[[sex]]$alpha[gr]
        l<-fit_list[[sex]]$lambda[gr]
        f_list[[gr]]<-dgamma(xvals,shape=a,rate=l)
      }
      m<-fit_list[[sex]]$pi[1]
      CDF_func<-function(q) {
        m*pgamma(q,shape=fit_list[[sex]]$alpha[1],rate=fit_list[[sex]]$lambda[1]) +
          (1-m)*pgamma(q,shape=fit_list[[sex]]$alpha[2],rate=fit_list[[sex]]$lambda[2])
      }
    } else if (fam=="norm") {
      for (gr in c(1,2)) {
        counter<-counter+1
        prop<-fit_list[[sex]]$pi[gr]
        mu<-fit_list[[sex]]$mu[gr]
        sd<-fit_list[[sex]]$sd[gr]
        f_list[[gr]]<-dnorm(xvals,mean=mu,sd=sd) 
      }
      m<-fit_list[[sex]]$pi[1]
      CDF_func<-function(q) {
        m*pnorm(q,mean=fit_list[[sex]]$mu[1],sd=fit_list[[sex]]$sd[1]) +
          (1-m)*pnorm(q,mean=fit_list[[sex]]$mu[2],sd=fit_list[[sex]]$sd[2])
      }
    } else if (fam=="weibull") {
      for (gr in c(1,2)) {
        counter<-counter+1
        prop<-fit_list[[sex]]$pi[gr]
        k<-fit_list[[sex]]$k[gr]
        lambda<-fit_list[[sex]]$lambda[gr]
        f_list[[gr]]<-dweibull(xvals,shape=k,scale=lambda)
      }
      m<-fit_list[[sex]]$pi[1]
      CDF_func<-function(q) {
        m*pweibull(q,shape=fit_list[[sex]]$k[1],scale=fit_list[[sex]]$lambda[1]) +
          (1-m)*pweibull(q,shape=fit_list[[sex]]$k[2],scale=fit_list[[sex]]$lambda[2])
      }
    }
    #
    # plots
    #
    y<-Dist_list[[sex]]   # empirical distribution
    p<-ppoints(length(y)) # probabilities for the quantiles 
    x<-rep(NA,length(p))  # theoretical quantiles go here
    age<-seq(0,150,0.1)
    temp<-CDF_func(age)
    for (pi in c(1:length(y))) {
      x[pi]<-age[which.min(abs(temp-p[pi]))] # find the age corresponding to CDF=p[pi]
    }
    #
    # retrieve BIC
    #
    df<-subset.data.frame(criteria_results,family==fam)
    df<-subset.data.frame(df,group==names(sexes[sex]))
    BIC<-round(df$BIC,2)
    qqplot(x,y,xlab=paste0("theoretical quantiles (",fam,")"),
           ylab=paste0("empirical quantiles"),
           main=paste0(names(sexes[sex])," (BIC=",BIC,")"))
    abline(0,1,col="red",lty=3)
  }
  dev.off()
}
#
# Add -LOG(KS_p)
#
criteria_results$LOGP_KS<-(-log10(criteria_results$KS_p))
#
# Plot BIC values ranked from highest (worst) to lowest (best)
#
df<-criteria_results # copy criteria
Massimo<-max(df$LOGP_KS)
for (mycol in c("AIC","BIC","KS_D")) { # scale in LOGP_KS scale
  maximus<-max(df[[mycol]])
  df[[mycol]]<-Massimo*df[[mycol]]/maximus
}
df<-df[order(df$BIC,decreasing=T),] # rank by BIC (from worst to best)
models<-c()
for (r in 1:nrow(df)) models[r]<-paste(df$family[r],"-",df$group[r])
#
file.name<-paste0("Plot/Goodness_of_fit.jpeg")
jpeg(file.name,quality=100,res=300,width=2000,height=2000)
par(mar=c(11,4,4,2))
barplot(df$LOGP_KS~seq(1,length(models)),width=1,space=0.2,ylim=c(0,Massimo),
        xlab="",ylab="",axes=F,xaxt="n",col="gray")
abline(h=-log10(0.05),col="red",lwd=2,lty=2)
par(new=T)
plot(seq(1,length(models)),df$BIC,pch=20,type="b",col="black",lty=1,lwd=2,xlab="",
     ylab="",xlim=c(0.5,length(models)+0.5),ylim=c(0,Massimo),axes=FALSE)
par(new=T)
plot(seq(1,length(models)),df$KS_D,pch=17,type="b",col="blue",lty=1,lwd=2,xlab="",
     ylab="",xlim=c(0.5,length(models)+0.5),ylim=c(0,Massimo),axes=FALSE)
text(seq(1,length(models)),-0.05*Massimo,labels=models,srt=45,adj=1,xpd=TRUE,cex=0.8)
axis(2)
#
# Legend
#
legend("topright",legend=c("BIC","KS_D","LOGP_KS"),col=c("black","blue","gray"),pch=c(20,17,15),
       lty=c(1,2,NA),lwd=2,bty="n")
dev.off()
#
# Write results (KS, BIC, AIC)
# 
write.table(criteria_results,file="Plot/Bimodal_fit_criteria.csv",sep=",",
            col.names=T,row.names=F,quote=F)
