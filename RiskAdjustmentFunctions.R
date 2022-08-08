#Setup tables for factorial analysis

riskfrac.VOICE=c(AG=0.51,MP=0.68,DR=0.26,FN=0.17,SE=0.75,ST=0.2)
getFactorTable=function(riskfrac){
  MP=rep(c(0,1),32)
  FN=rep(c(0,1),each=2,16)
  SE=rep(c(0,1),each=4,8)
  ST=rep(c(0,1),each=8,4)
  DR=rep(c(0,1),each=16,2)
  AG=rep(c(0,1),each=32)

  FactorTable=data.frame(MP=MP,FN=FN,SE=SE,ST=ST,DR=DR,AG=AG)
  components=names(FactorTable)
  

  
  
  
  multinomprob=function(x,y){
    prod(dbinom(as.numeric(x),size=1,prob=as.numeric(y)))
  }  
  FactorTable$populationfraction=apply(FactorTable,1,multinomprob,riskfrac[components])
  FactorTable$riskscore=2*MP+FN+2*SE+2*ST+DR+2*AG
  FactorTable
  
}


logitcombine=function(x,y){
  z=log(x/(1-x))+y
  ifelse(z>700,1,exp(z)/(1+exp(z)))
}

logcombine=function(x,y){
  z=log(x)+y
  ifelse(z>700,exp(700),exp(z))
}

parameterriskadjustment=function(params,risk){
  with(as.list(risk),{
    params$sexRiskV=logitcombine(params$sexRiskV,ST*params$sexriskST)
    params$sexRiskA=logitcombine(params$sexRiskA,ST*params$sexriskST)
    params$prevalence=logitcombine(params$prevalence,-(1-AG)*params$prevalenceAG)
    params$prevalence=logitcombine(params$prevalence,SE*params$prevalenceSE)
    params$prevalence=logitcombine(params$prevalence,DR*params$prevalenceDR)
    params$prevalence=logitcombine(params$prevalence,ST*params$prevalenceST)
    params$prevalence=logitcombine(params$prevalence,MP*params$prevalenceMP)
    params$prevalence=logitcombine(params$prevalence,FN*params$prevalenceFN)
    params$incidence=logitcombine(params$incidence,MP*params$prevalenceMP+FN*params$prevalenceFN+SE*params$prevalenceSE+ST*params$prevalenceST+DR*params$prevalenceDR)
    params$prevalence.casual=logitcombine(params$prevalence.casual,-(1-AG)*params$prevalencecasualAG)
    params$prevalence.casual=logitcombine(params$prevalence.casual,SE*params$prevalencecasualSE)
    params$prevalence.casual=logitcombine(params$prevalence.casual,DR*params$prevalencecasualDR)
    params$prevalence.casual=logitcombine(params$prevalence.casual,ST*params$prevalencecasualST)
    params$prevalence.casual=logitcombine(params$prevalence.casual,MP*params$prevalencecasualMP)
    params$prevalence.casual=logitcombine(params$prevalence.casual,FN*params$prevalencecasualFN)
    params$p0.low=logitcombine(params$p0.low,-(1-AG)*params$p0.lowAG)
    params$p0.low=logitcombine(params$p0.low,DR*params$p0.lowDR)
    params$p0.low=logitcombine(params$p0.low,SE*params$p0.lowSE)
    params$p0.low=logitcombine(params$p0.low,ST*params$p0.lowST)
    params$p0.low=logitcombine(params$p0.low,MP*params$p0.lowMP)
    params$p0.low=logitcombine(params$p0.low,FN*params$p0.lowFN)
    params$condomfrequency=logitcombine(params$condomfrequency,-(1-AG)*params$condomfrequencyAG)
    params$condomfrequency=logitcombine(params$condomfrequency,DR*params$condomfrequencyDR)
    params$condomfrequency=logitcombine(params$condomfrequency,SE*params$condomfrequencySE)
    params$condomfrequency=logitcombine(params$condomfrequency,ST*params$condomfrequencyST)
    params$condomfrequency=logitcombine(params$condomfrequency,MP*params$condomfrequencyMP)
    params$condomfrequency=logitcombine(params$condomfrequency,FN*params$condomfrequencyFN)
    params$analfrac=logitcombine(params$analfrac,-(1-AG)*params$analAG)
    params$analfrac=logitcombine(params$analfrac,DR*params$analDR)
    params$analfrac=logitcombine(params$analfrac,SE*params$analSE)
    params$analfrac=logitcombine(params$analfrac,ST*params$analST)
    params$analfrac=logitcombine(params$analfrac,MP*params$analMP)
    params$analfrac=logitcombine(params$analfrac,FN*params$analFN)
    params$analvaginal=logitcombine(params$analvaginal,-(1-AG)*params$analAG)
    params$analvaginal=logitcombine(params$analvaginal,DR*params$analDR)
    params$analvaginal=logitcombine(params$analvaginal,SE*params$analSE)
    params$analvaginal=logitcombine(params$analvaginal,ST*params$analST)
    params$analvaginal=logitcombine(params$analvaginal,MP*params$analMP)
    params$analvaginal=logitcombine(params$analvaginal,FN*params$analFN)
    params$anal.casual=logitcombine(params$analfrac,params$relanalcasual)
    params$sexrate=logcombine(params$sexrate,-(1-AG)*params$sexrateAG)
    params$sexrate=logcombine(params$sexrate,DR*params$sexrateDR)
    params$sexrate=logcombine(params$sexrate,SE*params$sexrateSE)
    params$sexrate=logcombine(params$sexrate,ST*params$sexrateST)
    params$sexrate=logcombine(params$sexrate,MP*params$sexrateMP)
    params$sexrate=logcombine(params$sexrate,FN*params$sexrateFN)
    params$sexcasual=logcombine(params$sexcasual,-(1-AG)*params$sexcasualAG)
    params$sexcasual=logcombine(params$sexcasual,DR*params$sexcasualDR)
    params$sexcasual=logcombine(params$sexcasual,SE*params$sexcasualSE)
    params$sexcasual=logcombine(params$sexcasual,ST*params$sexcasualST)
    params$sexcasual=logcombine(params$sexcasual,MP*params$sexcasualMP)
    params$sexcasual=logcombine(params$sexcasual,FN*params$sexcasualFN)
    params$vaginalonly=1-params$analvaginal
    params$cohabit=MP!=1
    params
  })
}

FactorRisk=function(params,full=TRUE,riskfrac=riskfrac.VOICE){
  FactorTable=getFactorTable(riskfrac)
  as.data.frame(cbind(FactorTable,t(apply(FactorTable,1,IncidenceRisk,params,full=full))))
}

FactorRisk.Fast=function(params,full=TRUE, riskfrac=riskfrac.VOICE){
  FactorTable=getFactorTable(riskfrac)
  as.data.frame(cbind(FactorTable,t(apply(FactorTable,1,IncidenceRisk.Fast,params,full=full))))
}

FactorSex=function(params, riskfrac=riskfrac.VOICE){
  FactorTable=getFactorTable(riskfrac)
  as.data.frame(cbind(FactorTable,t(apply(FactorTable,1,SexRisk,params))))
}

getIncidence.Risk=function(Risktable){
  365*sum(Risktable$numinfections*Risktable$populationfraction)/sum(Risktable$populationfraction*Risktable$numfollowups)
}

IncidenceRisk=function(risk,params,full=TRUE){
  params=parameterriskadjustment(params,risk)
  incidencecalc(params,full=full)
}


SexRisk=function(risk,params){
  params=parameterriskadjustment(params,risk)
  sexcalc(params)
}

incidencecalc=function(params,full=TRUE,eff=0,Ntimes=365){
  FT.all=fullTransition(params)
  fT.L=FT.all[[1]]
  fT.C=FT.all[[2]]
  p.low=InfectionPredict4(Ntimes,fT.L,eff,params)
  p.casual=InfectionPredict4(Ntimes,fT.C,eff,params)
  getIncidence(p.low,p.casual,params)
}

getIncidence=function(inc.low,inc.casual,params){
  with(as.list(params),{
    numinfections=inc.low[[1]]*p0.low+inc.casual[[1]]*(1-p0.low)
    numfollowups=inc.low[[2]]*p0.low+inc.casual[[2]]*(1-p0.low)
    c(incidence=365*numinfections/numfollowups,numinfections=numinfections,numfollowups=numfollowups)
  })
}

sexcalc=function(params){
  with(as.list(params),{
    lowa=sexrate*analvaginal*analfrac
    lowv=sexrate-lowa
    casa=lowa+anal.casual*sexcasual
    casv=lowv+(1-anal.casual)*sexcasual
    c(sexrateA.low=lowa*30,sexrateA.cas=casa*30,sexrateV.low=lowv*30,sexrateV.cas=casv*30,casualsexrate=sexcasual*90,casual=1-p0.low,condom=condomfrequency, prevalence=prevalence, incidence=incidence, prevalence.casual=prevalence.casual)
  })
}