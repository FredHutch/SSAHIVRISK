
#params.voice.orig=read.csv("calibratedparams_VOICE_Jberg.csv")

invert.theta=function(beta,p,phi){
  f=function(x){
    (sum(expit(x+beta)*p)-phi)^2
  }
  fp=function(x){
    sum(1/(+exp(-x-beta))^2*p)
  }
  
  lo=log(phi/sum(exp(beta)*p))
  
  hi=logit(phi)-min(beta)
  optimize(f, c(lo,hi))$minimum
}

invert.artprob=function(params,TARGET,VLS,...){
  relartinit=10
  f=function(x){
    p=params
    p['artinitasymptomatic']=x
    p['artinitlate']=relartinit*x
    if(VLS){
      return(getVLSprop(p)-TARGET)
    }
    else{
      return(getARTprop(p)-TARGET)
    }
  }
  
  #May have to adjust other parameters to match very high viral suppression
  while(f(.1)<0){
    params$artsuppressionrate = params$artsuppressionrate*1.1
    params$artquitrate = params$artquitrate*0.9
    params$artunsuppressionrate = params$artunsuppressionrate*0.9
  }
  x = uniroot(f, c(0,.1),tol=1e-7,...)$root
  params['artinitasymptomatic']=x
  params['artinitlate']=relartinit*x 
  params
}

convert.params=function(params,TARGET.ART,TARGET.incidence,TARGET.prevalence.main,TARGET.prevalence.casual,VLS){
  params=invert.artprob(params,TARGET.ART,VLS)
  relartinit=10
  #params['artinitasymptomatic']=x
  #params['artinitlate']=relartinit*x
  FS=FactorSex(params)
  params['prevalence']=logitcombine(params['prevalence'],invert.theta(logit(FS$prevalence),FS$populationfraction,TARGET.prevalence.main))
  p.casual=FS$casual*FS$casualsexrate*FS$populationfraction
  params['prevalence.casual']=logitcombine(params['prevalence.casual'],invert.theta(logit(FS$prevalence.casual),p.casual/sum(p.casual),TARGET.prevalence.casual))
  params['incidence']=logitcombine(params['incidence'],invert.theta(logit(FS$incidence),FS$populationfraction,TARGET.incidence))
  params
}


site.adjust.params=function(params.df,TARGET.ART.CI,TARGET.incidence.CI,TARGET.prevalence.main.CI,riskfrac,mp.effect.se=multipartnereffect, VLS=FALSE){
  Np=nrow(params.df)
  TARGET.ART=samplelogistic(TARGET.ART.CI,Np)
  TARGET.incidence=samplelogistic(TARGET.incidence.CI,Np)
  TARGET.prevalence.main=samplelogistic(TARGET.prevalence.main.CI,Np)
  
  mp.effect=rnorm(Np,mp.effect.se[1],mp.effect.se[2])
  TARGET.prevalence.casual=logitcombine(TARGET.prevalence.main,mp.effect)
  
  numinfections=numeric(Np)
  numfollowups=numeric(Np)
  for(i in seq(Np)){
    newparams=convert.params(params.df[i,],TARGET.ART[i], TARGET.incidence[i]/365, TARGET.prevalence.main[i], TARGET.prevalence.casual[i],VLS)
    Risktable=FactorRisk(newparams, riskfrac=riskfrac)
    Risktableyoung = subset(Risktable, AG==1)
    Risktableold = subset(Risktable, AG==0)
    numinfections[i]=sum(Risktable$numinfections*Risktable$populationfraction)
    numfollowups[i]=sum(Risktable$numfollowups*Risktable$populationfraction)
    #numinfectionsold[i]=sum(Risktableold$numinfections*Risktableold$populationfraction)
    #numfollowupsold[i]=sum(Risktableold$numfollowups*Risktableold$populationfraction)
    #numinfectionsyoung[i]=sum(Risktableyoung$numinfections*Risktableyoung$populationfraction)
    #numfollowupsyoung[i]=sum(Risktableyoung$numfollowups*Risktableyoung$populationfraction)
  }
  
  incidence=365*numinfections/numfollowups
  
  data.frame(numinfections=numinfections,
             numfollowups=numfollowups,
             incidence=incidence)
}

riskfrac.ASPIRE=c(AG=0.39, MP=0.58, DR=0.12, FN=0.46, SE=0.57, ST=0.21)
