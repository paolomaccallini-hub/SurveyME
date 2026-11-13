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
age_labels<-list(c("Single"),c("Younger","Older"),c("Younger","Middle","Older"))
colv<-c("pink","cyan","gray")
xvals<-seq(0,100,1)
criteria_results<-data.frame() # goodness of fit tests
#
# Bimodal analysis
#
family<-c("gamma","lnorm","norm","weibull")
fit_list<-list()
components<-c(1,2,3)
#
for (comp in components) {
  age<-age_labels[[comp]]
  for (fam in family) {
    print(fam)
    for (j in c(1,2,3)) {
      print(names(sexes)[j])
      fit_list[[j]]<-mixfit(Dist_list[[j]],family=fam,ncomp=comp)
    }
    #
    # Plotting empirical distributions and models
    #
    file.name<-paste0("Plot/Age_at_onset_",fam,"_",comp,".jpeg")
    jpeg(file.name,quality=100,res=300,width=2000,height=2000)
    lay.mat<-matrix(c(1,3,2,3),2,2)
    layout(mat=lay.mat)
    #
    myresults<-data.frame(group=rep(NA,3*comp))
    counter<-0
    qq_list<-list()
    for (sex in c(1,2,3)) {
      #
      f_list<-list()
      if (fam=="lnorm") {
        for (gr in seq(1,comp)) {
          counter<-counter+1
          print(counter)
          prop<-fit_list[[sex]]$pi[gr]
          m<-fit_list[[sex]]$mulog[gr]
          s<-fit_list[[sex]]$sdlog[gr]
          myresults$group[counter]<-paste(age[gr],names(sexes[sex]))
          myresults$proportion[counter]<-round(prop,2)
          myresults$mean[counter]<-round(exp(m+0.5*s^2),2)  
          myresults$SD[counter]<-round(sqrt(exp(2*m+s^2)*(exp(s^2)-1)),2)
          f_list[[gr]]<-dlnorm(xvals,meanlog=m,sdlog=s)
        }
        prop<-fit_list[[sex]]$pi
        for (gr in seq(1,comp)) {
          if (comp==1) {
            CDF_func<-function(q) {
              prop[1]*plnorm(q,meanlog=fit_list[[sex]]$mulog[1],sdlog=fit_list[[sex]]$sdlog[1]) 
            }
            density_func<-function(q) {
              prop[1]*dlnorm(q,meanlog=fit_list[[sex]]$mulog[1],sdlog=fit_list[[sex]]$sdlog[1]) 
            }  
          } else if (comp==2) {
            CDF_func<-function(q) {
              prop[1]*plnorm(q,meanlog=fit_list[[sex]]$mulog[1],sdlog=fit_list[[sex]]$sdlog[1]) +
                prop[2]*plnorm(q,meanlog=fit_list[[sex]]$mulog[2],sdlog=fit_list[[sex]]$sdlog[2])
            }
            density_func<-function(q) {
              prop[1]*dlnorm(q,meanlog=fit_list[[sex]]$mulog[1],sdlog=fit_list[[sex]]$sdlog[1]) +
                prop[2]*dlnorm(q,meanlog=fit_list[[sex]]$mulog[2],sdlog=fit_list[[sex]]$sdlog[2])
            }  
          } else if (comp==3) {
            CDF_func<-function(q) {
              prop[1]*plnorm(q,meanlog=fit_list[[sex]]$mulog[1],sdlog=fit_list[[sex]]$sdlog[1]) +
                prop[2]*plnorm(q,meanlog=fit_list[[sex]]$mulog[2],sdlog=fit_list[[sex]]$sdlog[2]) +
                prop[3]*plnorm(q,meanlog=fit_list[[sex]]$mulog[3],sdlog=fit_list[[sex]]$sdlog[3])
            }
            density_func<-function(q) {
              prop[1]*dlnorm(q,meanlog=fit_list[[sex]]$mulog[1],sdlog=fit_list[[sex]]$sdlog[1]) +
                prop[2]*dlnorm(q,meanlog=fit_list[[sex]]$mulog[2],sdlog=fit_list[[sex]]$sdlog[2]) +
                prop[3]*dlnorm(q,meanlog=fit_list[[sex]]$mulog[3],sdlog=fit_list[[sex]]$sdlog[3])
            }  
          }
        }
      } else if (fam=="gamma") {
        for (gr in seq(1,comp)) {
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
        prop<-fit_list[[sex]]$pi
        for (gr in seq(1,comp)) {
          if (comp==1) {
            CDF_func<-function(q) {
              prop[1]*pgamma(q,shape=fit_list[[sex]]$alpha[1],rate=fit_list[[sex]]$lambda[1]) 
            }
            density_func<-function(q) {
              prop[1]*dgamma(q,shape=fit_list[[sex]]$alpha[1],rate=fit_list[[sex]]$lambda[1]) 
            }  
          } else if (comp==2) {
            CDF_func<-function(q) {
              prop[1]*pgamma(q,shape=fit_list[[sex]]$alpha[1],rate=fit_list[[sex]]$lambda[1]) +
              prop[2]*pgamma(q,shape=fit_list[[sex]]$alpha[2],rate=fit_list[[sex]]$lambda[2])
            }
            density_func<-function(q) {
              prop[1]*dgamma(q,shape=fit_list[[sex]]$alpha[1],rate=fit_list[[sex]]$lambda[1]) +
              prop[2]*dgamma(q,shape=fit_list[[sex]]$alpha[2],rate=fit_list[[sex]]$lambda[2])
            }  
          } else if (comp==3) {
            CDF_func<-function(q) {
              prop[1]*pgamma(q,shape=fit_list[[sex]]$alpha[1],rate=fit_list[[sex]]$lambda[1]) +
              prop[2]*pgamma(q,shape=fit_list[[sex]]$alpha[2],rate=fit_list[[sex]]$lambda[2]) +
              prop[3]*pgamma(q,shape=fit_list[[sex]]$alpha[3],rate=fit_list[[sex]]$lambda[3])
            }
            density_func<-function(q) {
              prop[1]*dgamma(q,shape=fit_list[[sex]]$alpha[1],rate=fit_list[[sex]]$lambda[1]) +
              prop[2]*dgamma(q,shape=fit_list[[sex]]$alpha[2],rate=fit_list[[sex]]$lambda[2]) +
              prop[3]*dgamma(q,shape=fit_list[[sex]]$alpha[3],rate=fit_list[[sex]]$lambda[3])
            }  
          }
        }
      } else if (fam=="norm") {
        for (gr in seq(1,comp)) {
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
        prop<-fit_list[[sex]]$pi
        for (gr in seq(1,comp)) {
          if (comp==1) {
            CDF_func<-function(q) {
              prop[1]*pnorm(q,mean=fit_list[[sex]]$mu[1],sd=fit_list[[sex]]$sd[1]) 
            }
            density_func<-function(q) {
              prop[1]*dnorm(q,mean=fit_list[[sex]]$mu[1],sd=fit_list[[sex]]$sd[1]) 
            }  
          } else if (comp==2) {
            CDF_func<-function(q) {
              prop[1]*pnorm(q,mean=fit_list[[sex]]$mu[1],sd=fit_list[[sex]]$sd[1]) +
                prop[2]*pnorm(q,mean=fit_list[[sex]]$mu[2],sd=fit_list[[sex]]$sd[2])
            }
            density_func<-function(q) {
              prop[1]*dnorm(q,mean=fit_list[[sex]]$mu[1],sd=fit_list[[sex]]$sd[1]) +
                prop[2]*dnorm(q,mean=fit_list[[sex]]$mu[2],sd=fit_list[[sex]]$sd[2])
            }  
          } else if (comp==3) {
            CDF_func<-function(q) {
              prop[1]*pnorm(q,mean=fit_list[[sex]]$mu[1],sd=fit_list[[sex]]$sd[1]) +
                prop[2]*pnorm(q,mean=fit_list[[sex]]$mu[2],sd=fit_list[[sex]]$sd[2]) +
                prop[3]*pnorm(q,mean=fit_list[[sex]]$mu[3],sd=fit_list[[sex]]$sd[3])
            }
            density_func<-function(q) {
              prop[1]*dnorm(q,mean=fit_list[[sex]]$mu[1],sd=fit_list[[sex]]$sd[1]) +
                prop[2]*dnorm(q,mean=fit_list[[sex]]$mu[2],sd=fit_list[[sex]]$sd[2]) +
                prop[3]*dnorm(q,mean=fit_list[[sex]]$mu[3],sd=fit_list[[sex]]$sd[3])
            }  
          }
        }
      } else if (fam=="weibull") {
        for (gr in seq(1,comp)) {
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
        prop<-fit_list[[sex]]$pi
        for (gr in seq(1,comp)) {
          if (comp==1) {
            CDF_func<-function(q) {
              prop[1]*pweibull(q,shape=fit_list[[sex]]$k[1],scale=fit_list[[sex]]$lambda[1]) 
            }
            density_func<-function(q) {
              prop[1]*dweibull(q,shape=fit_list[[sex]]$k[1],scale=fit_list[[sex]]$lambda[1]) 
            }  
          } else if (comp==2) {
            CDF_func<-function(q) {
              prop[1]*pweibull(q,shape=fit_list[[sex]]$k[1],scale=fit_list[[sex]]$lambda[1]) +
                prop[2]*pweibull(q,shape=fit_list[[sex]]$k[2],scale=fit_list[[sex]]$lambda[2])
            }
            density_func<-function(q) {
              prop[1]*dweibull(q,shape=fit_list[[sex]]$k[1],scale=fit_list[[sex]]$lambda[1]) +
                prop[2]*dweibull(q,shape=fit_list[[sex]]$k[2],scale=fit_list[[sex]]$lambda[2])
            }  
          } else if (comp==3) {
            CDF_func<-function(q) {
              prop[1]*pweibull(q,shape=fit_list[[sex]]$k[1],scale=fit_list[[sex]]$lambda[1]) +
                prop[2]*pweibull(q,shape=fit_list[[sex]]$k[2],scale=fit_list[[sex]]$lambda[2]) +
                prop[3]*pweibull(q,shape=fit_list[[sex]]$k[3],scale=fit_list[[sex]]$lambda[3])
            }
            density_func<-function(q) {
              prop[1]*dweibull(q,shape=fit_list[[sex]]$k[1],scale=fit_list[[sex]]$lambda[1]) +
                prop[2]*dweibull(q,shape=fit_list[[sex]]$k[2],scale=fit_list[[sex]]$lambda[2]) +
                prop[3]*dweibull(q,shape=fit_list[[sex]]$k[3],scale=fit_list[[sex]]$lambda[3])
            }  
          }
        }
      }
      #
      # Log-likelihood, AIC and BIC
      #
      fit<-fit_list[[sex]]
      #
      loglik<-sum(log(pmax(density_func(Dist_list[[sex]]),.Machine$double.xmin)))
      n<-length(Dist_list[[sex]])
      npar<-2*comp+comp-1
      AIC_value<-2*npar-2*loglik
      BIC_value<-npar*log(n)-2*loglik
      #
      # Kolmogorov test
      #
      a<-ks.test(Dist_list[[sex]],CDF_func) # KS test
      pval<-format(a$p.value,digits=3)
      #
      # Store quantiles for QQ plots
      #
      p<-ppoints(length(Dist_list[[sex]]))
      age_grid<-seq(0,150,0.1)
      temp<-CDF_func(age_grid)
      qq_list[[sex]]<-sapply(p,function(pp) age_grid[which.min(abs(temp-pp))])
      #
      # Store results for goodness of fit
      #
      df<-data.frame(family=fam,group=names(sexes)[sex],n=n,logLik=loglik,components=comp,
                     parameters=npar,AIC=AIC_value,BIC=BIC_value,KS_D=unname(a$statistic),
                     KS_p=a$p.value)
      criteria_results<-rbind(criteria_results,df)
      #
      # plots
      #
      hist(Dist_list[[sex]],breaks=seq(0,100,2),freq=F,xlim=c(5,90),ylim=c(0,0.04),xlab="age at onset",
           ylab="density of probability",main=paste0(names(sexes[sex])," (KS pval=",pval,")"),
           col=colv[sex])
      #
      for (gr in seq(1,comp)) {
        par(new=T)
        plot(xvals,prop[gr]*f_list[[gr]],type="l",col="black",lty=2,lwd=2,
             xlim=c(5,90),ylim=c(0,0.04),xlab="",ylab="")
      }
      par(new=T)
      plot(xvals,density_func(xvals),type="l",col="red",lty=1,lwd=2,
           xlim=c(5,90),ylim=c(0,0.04),xlab="",ylab="")
    }
    dev.off()
    #
    # Write component parameters (proportion, mean, SD)
    #
    file.name<-paste0("Plot/Age_at_onset_",fam,"_",comp,"_components.tsv")
    write.table(myresults,file=file.name,sep="\t",col.names=T,row.names=F,quote=F)
    #
    # Plotting QQplots: empirical distributions against models
    #
    file.name<-paste0("Plot/Age_at_onset_",fam,"_",comp,"_QQ.jpeg")
    jpeg(file.name,quality=100,res=300,width=3000,height=1000)
    lay.mat<-matrix(c(1,2,3),1,3)
    layout(mat=lay.mat)
    #
    counter<-0
    for (sex in c(1,2,3)) {
      #
      # plots
      #
      y<-Dist_list[[sex]] # empirical distribution
      x<-qq_list[[sex]] # theoretical quantiles
      #
      # retrieve BIC
      #
      df<-subset.data.frame(criteria_results,family==fam)
      df<-subset.data.frame(df,group==names(sexes[sex]))
      df<-subset.data.frame(df,components==comp)
      BIC<-round(df$BIC,2)
      qqplot(x,y,xlab=paste0("theoretical quantiles (",fam,")"),
             ylab=paste0("empirical quantiles"),
             main=paste0(names(sexes[sex])," (BIC=",BIC,")"))
      abline(0,1,col="red",lty=3)
    }
    dev.off()
  }
}
#
# Add -LOG(KS_p)
#
criteria_results$LOGP_KS<-(-log10(criteria_results$KS_p))
#
# Plot BIC, KS_D and LOGP_10, one plot per dataset
#
for (comp in components) {
  for (sex in c(1,2,3)) {
    df<-criteria_results # copy criteria
    df<-subset.data.frame(df,group==names(sexes[sex]))
    df<-subset.data.frame(df,components==comp)
    Massimo<-max(df$LOGP_KS)
    for (mycol in c("AIC","BIC","KS_D")) { # scale in LOGP_KS scale
      maximus<-max(df[[mycol]])
      minimum<-min(df[[mycol]])
      df[[mycol]]<-Massimo*(df[[mycol]]-minimum)/(maximus-minimum)
    }
    models<-c()
    for (r in 1:nrow(df)) models[r]<-paste(df$family[r],"-",df$group[r])
    #
    file.name<-paste0("Plot/Goodness_of_fit_",names(sexes[sex]),"_",comp,".jpeg")
    jpeg(file.name,quality=100,res=300,width=1500,height=2000)
    par(mar=c(11,4,4,2))
    barplot(df$LOGP_KS~seq(1,length(models)),width=1,space=0.2,ylim=c(0,1.3*Massimo),
            xlab="",ylab="",axes=F,xaxt="n",col="gray")
    abline(h=-log10(0.05),col="red",lwd=2,lty=2)
    par(new=T)
    plot(seq(1,length(models)),df$BIC,pch=20,type="b",col="black",lty=1,lwd=2,xlab="",
         ylab="",xlim=c(0.5,length(models)+0.5),ylim=c(0,1.3*Massimo),axes=FALSE)
    par(new=T)
    plot(seq(1,length(models)),df$KS_D,pch=17,type="b",col="blue",lty=1,lwd=2,xlab="",
         ylab="",xlim=c(0.5,length(models)+0.5),ylim=c(0,1.3*Massimo),axes=FALSE)
    text(seq(1,length(models)),-0.1*Massimo,labels=models,srt=45,adj=1,xpd=TRUE,cex=0.8)
    axis(2)
    #
    # Legend
    #
    legend("topright",legend=c("BIC","KS_D","LOGP_KS"),col=c("black","blue","gray"),pch=c(20,17,15),
           lty=c(1,2,NA),lwd=2,bty="n")
    dev.off()  
  }  
}
#
# Write results (KS, BIC, AIC)
# 
write.table(criteria_results,file="Plot/Bimodal_fit_criteria.csv",sep=",",
            col.names=T,row.names=F,quote=F)
