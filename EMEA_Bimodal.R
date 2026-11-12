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
family<-c("gamma", "lnorm", "norm", "weibull")
fit_list<-list()
criteria_results<-data.frame() # goodness of fit tests
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
        m*pnorm(q,mean=fit_list[[sex]]$mu[1],sd=fit_list[[sex]]$sd[1],) +
          (1-m)*pnorm(q,mean=fit_list[[sex]]$mu[2],sd=fit_list[[sex]]$sd[2])
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
    }
    #
    # Log-likelihood, AIC and BIC
    #
    fit<-fit_list[[sex]]
    #
    density_func<-switch(fam,
                         
                         gamma = function(x) {
                           fit$pi[1]*dgamma(x,shape=fit$alpha[1],rate=fit$lambda[1]) +
                             fit$pi[2]*dgamma(x,shape=fit$alpha[2],rate=fit$lambda[2])
                         },
                         
                         lnorm = function(x) {
                           fit$pi[1]*dlnorm(x,meanlog=fit$mulog[1],sdlog=fit$sdlog[1]) +
                             fit$pi[2]*dlnorm(x,meanlog=fit$mulog[2],sdlog=fit$sdlog[2])
                         },
                         
                         norm = function(x) {
                           fit$pi[1]*dnorm(x,mean=fit$mu[1],sd=fit$sd[1]) +
                             fit$pi[2]*dnorm(x,mean=fit$mu[2],sd=fit$sd[2])
                         },
                         
                         weibull = function(x) {
                           fit$pi[1]*dweibull(x,shape=fit$k[1],scale=fit$lambda[1]) +
                             fit$pi[2]*dweibull(x,shape=fit$k[2],scale=fit$lambda[2])
                         }
    )
    
    loglik<-sum(log(pmax(density_func(Dist_list[[sex]]),
                         .Machine$double.xmin)))
    
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
    criteria_results<-rbind(criteria_results,
                            data.frame(
                              family=fam,
                              group=names(sexes)[sex],
                              n=n,
                              logLik=loglik,
                              parameters=npar,
                              AIC=AIC_value,
                              BIC=BIC_value,
                              KS_D=unname(a$statistic),
                              KS_p=a$p.value
                            )
    )
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
  # Save metadata
  #
  file.name<-paste0("Plot/",fam,"_bimodal_fit.csv")
  write.table(myresults,file=file.name,sep=",",col.names=T,row.names=F,quote=F)
}
#
# Write results (KS, BIC, AIC)
# 
write.table(criteria_results,
            file="Plot/Bimodal_fit_criteria.csv",
            sep="t",
            col.names=T,
            row.names=F,
            quote=F)