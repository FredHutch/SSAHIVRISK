library(plyr)

dat=readRDS("math_model.rds")

myscore=2*(dat$VRSmp==2)+(dat$VRSfn==2)*1+(dat$VRSse!="No")*2+2*(dat$VRSst!=2)+(dat$VRSdr>0)*1+(dat$VRSag=="Yes")*2

#Fractions with each risk factor in VOICE cohort
riskfrac.VOICE=c(AG=0.51,MP=0.68,DR=0.26,FN=0.17,SE=0.75,ST=0.2)

#Standardize voice risk components
dat$VRSscore=myscore
dat$MP=(dat$VRSmp==2)*1
dat$FN=(dat$VRSfn==2)*1
dat$SE=(dat$VRSse!="No")*1
dat$ST=(dat$VRSst!=2)*1
dat$DR=(dat$VRSdr>0)*1
dat$AG=(dat$VRSag=="Yes")*1

#Setup tables for factorial analysis
MP=rep(c(0,1),32)
FN=rep(c(0,1),each=2,16)
SE=rep(c(0,1),each=4,8)
ST=rep(c(0,1),each=8,4)
DR=rep(c(0,1),each=16,2)
AG=rep(c(0,1),each=32)

FactorTable=data.frame(MP=MP,FN=FN,SE=SE,ST=ST,DR=DR,AG=AG)
components=names(FactorTable)
FactorList=split(FactorTable,seq(64))

multinomprob=function(x,y){
  prod(dbinom(as.numeric(x),size=1,prob=as.numeric(y)))
}  
populationfraction.VOICE=apply(FactorTable,1,multinomprob,riskfrac.VOICE[components])

FactorTable$populationfraction=populationfraction.VOICE

FactorTable$riskscore=2*MP+FN+2*SE+2*ST+DR+2*AG

test=rep(NA,length(dat$MPHIVSTAT))
test[dat$MPHIVSTAT=='HIV positive']=1
test[dat$MPHIVSTAT=='HIV negative']=0
dat$test=test



agecat=rep(NA,length(dat$MPHIVSTAT))
agecat[dat$MPAGE_TEXT<20&dat$MPAGE_TEXT>=15]=1
agecat[dat$MPAGE_TEXT<25&dat$MPAGE_TEXT>=20]=2
agecat[dat$MPAGE_TEXT<30&dat$MPAGE_TEXT>=25]=3
agecat[dat$MPAGE_TEXT<35&dat$MPAGE_TEXT>=30]=4
agecat[dat$MPAGE_TEXT<40&dat$MPAGE_TEXT>=35]=5
agecat[dat$MPAGE_TEXT<45&dat$MPAGE_TEXT>=40]=6
dat$agecat=agecat

N.age.082=sapply(seq(6),function(x){sum(agecat==x,na.rm=TRUE)})

condomusages=c("Never","Rarely","Sometimes","Often","Always")
tmp=length(condomusages)
condomfreq=numeric(tmp)
condomfreqA=numeric(tmp)
condomfreqV=numeric(tmp)
dat$CONDOMA=NA
dat$CONDOMV=NA
dat$CONDOM=NA
dat$CONDOMLASTA=NA
dat$CONDOMLASTV=NA
dat$CONDOMLASTA[dat$ASXCNDLAST=="Yes"]=1
dat$CONDOMLASTA[dat$ASXCNDLAST=="Part of the last time we had sex"]=0
dat$CONDOMLASTA[dat$ASXCNDLAST=="No"]=0
dat$CONDOMLASTV[dat$VSXCNDLAST=="Yes"]=1
dat$CONDOMLASTV[dat$VSXCNDLAST=="Part of the last time we had sex"]=0
dat$CONDOMLASTV[dat$VSXCNDLAST=="No"]=0

for(i in seq(tmp)){
  condomlastA=subset(dat,ASXCNDMO==condomusages[i])$ASXCNDLAST
  condomlastV=subset(dat,VSXCNDMO==condomusages[i])$VSXCNDLAST
  condomlast=c(condomlastA,condomlastV)
  condomfreq[i]=mean(condomlast==5)*1+mean(condomlast==3)*.5
  condomfreqA[i]=mean(condomlastA=="Yes")*1+mean(condomlastA=="Part of the last time we had sex")*.5
  condomfreqV[i]=mean(condomlastV=="Yes")*1+mean(condomlastV=="Part of the last time we had sex")*.5
  if(!condomusages[i]%in%c("","Prefer not to answer")){
    dat$CONDOMA[dat$ASXCNDMO==condomusages[i]]=condomfreqA[i]
    dat$CONDOMV[dat$VSXCNDMO==condomusages[i]]=condomfreqV[i]
  }
}

dat$CONDOM=(dat$CONDOMA*dat$ASXMONTH_TEXT+dat$CONDOMV*dat$VSXMONTH_TEXT)/(dat$ASXMONTH_TEXT+dat$VSXMONTH_TEXT)

riskfrac.082=sapply(components,function(x){mean(dat[,x])})


FitRisk1=function(x,...){
  lm(x~MP+FN+SE+ST+DR, dat,...)
}

FitRisk2=function(x,...){
  lm(x~(dat$VRSmp==2)*(dat$VRSfn==2)*(dat$VRSse!="No")*(dat$VRSst!=2)*(dat$VRSdr>0),...)
}

FitRisk3=function(x,...){
  glm(x~MP+FN+SE+ST+DR, dat, family=binomial,...)
}

FitRisk4=function(x,...){
  glm(x~MP+FN+SE+ST+DR+(SPNUM_TEXT>1), dat,family=binomial,...)
}


FitRisk5=function(x,...){
  glm(x~(MP+FN+SE+ST+DR)*(SPNUM_TEXT>1), dat,family=binomial,...)
}

FitRisk6=function(x,...){
  glm(x~(MP+FN+SE+ST+DR)*(SPNUM_TEXT>1), dat,family=poisson,...)
}

FitRisk7=function(x,...){
  glm(x~MP+FN+SE+ST+DR+(SPNUM_TEXT>1), dat,family=poisson,...)
}

ExtractEstimates=function(x){
  slm=summary(x)
  values=slm$coefficients[,1]
  stderrs=slm$coefficients[,2]
  data.frame(values=values,stderrs=stderrs)
}

#Extract regression parameters to use as starting guess for parameters

# MAINSEXA=exp(ExtractEstimates(FitRisk6(dat$ASXMONTH_TEXT))[1,1])/30 #Estimate anal sex rate in zero risk score/low risk
# MAINSEXmed=ExtractEstimates(FitRisk6(dat$ASXMONTH_TEXT/30+dat$VSXMONTH_TEXT/30,subset=dat$SPNUM_TEXT==1))[1,1] 
# MAINSEXse=ExtractEstimates(FitRisk6(dat$ASXMONTH_TEXT/30+dat$VSXMONTH_TEXT/30,subset=dat$SPNUM_TEXT==1))[1,2] 
# 
# SEXDIFFA=ExtractEstimates(FitRisk6(dat$ASXMONTH_TEXT))[2:6,1]/30 #How anal sex rate is modified by risk components
# CASUALSEXA=exp(ExtractEstimates(FitRisk6(dat$ASXMONTH_TEXT))[7,1])/30 #How anal sex rate increases with more than one partner
# CASUALSEXARISKS=ExtractEstimates(FitRisk6(dat$ASXMONTH_TEXT))[8:12,1]/30 #How anal sex rate changes with risk components
# 
# 
# #Repeat for vaginal sex
# MAINSEXV=exp(ExtractEstimates(FitRisk6(dat$VSXMONTH_TEXT))[1,1])/30
# MAINSEXVmed=ExtractEstimates(FitRisk6(dat$VSXMONTH_TEXT/30))[1,1]
# MAINSEXVse=ExtractEstimates(FitRisk6(dat$VSXMONTH_TEXT/30))[1,2]
# 
# SEXDIFFV=ExtractEstimates(FitRisk6(dat$VSXMONTH_TEXT))[2:6,1]/30
# CASUALSEXV=mean(dat$VSXMONTH_TEXT[dat$SPNUM_TEXT>1],na.rm=TRUE)-mean(dat$VSXMONTH_TEXT[dat$SPNUM_TEXT<2],na.rm=TRUE)
# CASUALSEXVRISKS=ExtractEstimates(FitRisk6(dat$VSXMONTH_TEXT))[8:12,1]/30 #How anal sex rate changes with risk components
# 
# MAINSEXRATE=MAINSEXA+MAINSEXV #Overall sex rate with main partners
# CASUALSEXRATE=CASUALSEXA+CASUALSEXV #Overall sex rate with casual partners
# 
# CASUALSEXmed=ExtractEstimates(FitRisk6(dat$SPNUM_TEXT,subset=dat$SPNUM_TEXT>1))[1,1] 
# CASUALSEXse=ExtractEstimates(FitRisk6(dat$SPNUM_TEXT,subset=dat$SPNUM_TEXT>1))[1,2]
# 
# #Estimate how risk components determines whether an individual practices anal sex
# LOGITPRACTICEANALMAIN=ExtractEstimates(FitRisk3(dat$ASXMONTH_TEXT>0,subset=dat$SPNUM_TEXT<2))[,1]
# LOGITPRACTICEANALMAINmed=ExtractEstimates(FitRisk3(dat$ASXMONTH_TEXT>0,subset=dat$SPNUM_TEXT<2))[1,1]
# LOGITPRACTICEANALMAINse=ExtractEstimates(FitRisk3(dat$ASXMONTH_TEXT>0,subset=dat$SPNUM_TEXT<2))[1,2]
# 
# #Estimate how risk components determines how frequent anal sex is (among those who practice)
# LOGITANALFRAC=ExtractEstimates(FitRisk3(dat$ASXMONTH_TEXT/(dat$VSXMONTH_TEXT+dat$ASXMONTH_TEXT),subset=dat$SPNUM_TEXT<2&dat$ASXMONTH_TEXT>0))[,1]
# LOGITANALFRACmed=ExtractEstimates(FitRisk3(dat$ASXMONTH_TEXT/(dat$VSXMONTH_TEXT+dat$ASXMONTH_TEXT),subset=dat$SPNUM_TEXT<2&dat$ASXMONTH_TEXT>0))[1,1]
# LOGITANALFRACse=ExtractEstimates(FitRisk3(dat$ASXMONTH_TEXT/(dat$VSXMONTH_TEXT+dat$ASXMONTH_TEXT),subset=dat$SPNUM_TEXT<2&dat$ASXMONTH_TEXT>0))[2,1]
# 
# 
# #Estimate how risk components determine probability of casual sex
# LOGITP0LOW=ExtractEstimates(FitRisk3(dat$SPNUM_TEXT>1))[,1]
# LOGITP0LOWmed=ExtractEstimates(FitRisk3(dat$SPNUM_TEXT>1))[1,1]
# LOGITP0LOWse=ExtractEstimates(FitRisk3(dat$SPNUM_TEXT>1))[1,2]
# 
# #Estimate how risk components determine probability of condom usage
# LOGITCONDOM=ExtractEstimates(FitRisk3(dat$CONDOM))[,1]
# LOGITCONDOMmed=ExtractEstimates(FitRisk3(dat$CONDOM))[1,1]
# LOGITCONDOMse=ExtractEstimates(FitRisk3(dat$CONDOM))[1,2]
# 
# ANALFRAC=as.numeric(exp(LOGITANALFRAC[1])/(1+exp(LOGITANALFRAC[1])))
# PRACTICEANALMAIN=as.numeric(exp(LOGITPRACTICEANALMAIN[1])/(1+exp(LOGITPRACTICEANALMAIN[1])))
# P0LOW=as.numeric(exp(LOGITP0LOW[1])/(1+exp(LOGITP0LOW[1])))
# CONDOM=exp(LOGITCONDOM[1])/(1+exp(LOGITCONDOM[1]))

#Assume for now that casual sex is 20% anal
ANALCAS=0.2


test=rep(NA,length(dat$MPHIVSTAT))
test[dat$MPHIVSTAT=='HIV positive']=1
test[dat$MPHIVSTAT=='HIV negative']=0
dat$test=test

agecat=rep(NA,length(dat$MPHIVSTAT))
agecat[dat$MPAGE_TEXT<20&dat$MPAGE_TEXT>=15]=1
agecat[dat$MPAGE_TEXT<25&dat$MPAGE_TEXT>=20]=2
agecat[dat$MPAGE_TEXT<30&dat$MPAGE_TEXT>=25]=3
agecat[dat$MPAGE_TEXT<35&dat$MPAGE_TEXT>=30]=4
agecat[dat$MPAGE_TEXT<40&dat$MPAGE_TEXT>=35]=5
agecat[dat$MPAGE_TEXT<45&dat$MPAGE_TEXT>=40]=6
dat$agecat=agecat

N.age.082=sapply(seq(6),function(x){sum(agecat==x,na.rm=TRUE)})

populationfraction.082=numeric(64)
X.082=1+dat$MP+2*dat$FN+4*dat$SE+8*dat$ST+16*dat$DR+32*dat$AG
populationfraction.082[count(X.082)[,1]]=count(X.082)[,2]/length(X.082)