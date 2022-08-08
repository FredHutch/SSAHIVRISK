VOICE.makeup=c(Durban=3110,Johannesburg=354+350,Klerksdorp=263,Kampala=322,Harare=205+425,Capetown=0)
HPTN082.makeup=c(Harare=133,Johannesburg=133,Capetown=133)
ASPIRE.makeup=c(Durban=244+180+150+117+103+103+150,Johannesburg=213,Kampala=253,Harare=230+224+224, 
                   Capetown=166, Blantyre=130, Lilongwe=142)
ECHO.makeup = c(Buffalo = 615, Capetown = 560, Durban = 861, Manzini = 502, Kisumu = 901, Tshwane = 407+810, Johannesburg = 697, Klerksdorp = 555, Kwazulu = 653+611, Lusaka = 658)

HPTN035.makeup=c(Blantyre=441, Lilongwe=596, Hlabisa=346,Durban=702,Lusaka=319,Harare=223+260)
FEMPREP.makeup=c(Bondo=364+356,Tshwane=378+372,Manguang=261+268)

HPTN084.makeup=c(Johannesburg = 201 + 175, Gaborone = 91, Harare = 162 + 166 + 153 + 160 + 138, Capetown = 148 + 223, Blantyre = 113, 
                 Lilongwe = 111, Kisumu = 66, Kampala = 210 + 182 + 204, Durban = 225 + 155 + 149, Siteki = 159)
#########Useful Functions

logit=function(x){log(x/(1-x))}
expit=function(x){exp(x)/(1+exp(x))}

CItoSE=function(x){
  with(as.list(x),{
    hi=log(hi/(1-hi))
    lo=log(lo/(1-lo))
    med=(hi+lo)/2
    c(med=med,se=(hi-med)/qnorm(.975))})
}

CombineSE2=function(x,y,ADD=FALSE){
  x=as.numeric(x)
  y=as.numeric(y)
  CombineSE(x[1],y[1],x[2],y[2],ADD=ADD)
}

CombineSE=function(med1,med2,se1,se2,ADD=FALSE){
  if(!ADD){
    return(c(med1-med2,sqrt(se1^2+se2^2)))
  }
  c(med=med1+med2,se=sqrt(se1^2+se2^2))
}

CombineCI=function(x1,x2,ADD=FALSE){
  y1=CItoSE(x1)
  y2=CItoSE(x2)
  CombineSE(y1[1],y2[1],y1[2],y2[2],ADD=ADD)
}

MultCI=function(x1,x2,ADD=FALSE){
  y1=CItoSE(x1)
  y2=CItoSE(x2)
  MultSE(y1[1],y2[1],y1[2],y2[2])
}

AddEffect=function(base,effect){
  y1=as.numeric(CItoSE(base))
  SEtoCI(CombineSE(y1[1],as.numeric(effect['med']),y1[2],as.numeric(effect['se']),ADD=TRUE))
}

SEtoCI=function(x){
  with(as.list(x),{
    hi=exp(med+qnorm(.975)*se)/(1+exp(med+qnorm(.975)*se))
    lo=exp(med-qnorm(.975)*se)/(1+exp(med-qnorm(.975)*se))
    c(lo=lo,hi=hi)
  })
}

SEtoCI2=function(x){
  with(as.list(x),{
    hi=exp(med+qnorm(.975)*se)
    lo=exp(med-qnorm(.975)*se)
    c(lo=lo,hi=hi)
  })
}



AgeEffectTable=function(med,lb,base){
  ages=c('15-19','20-24','25-29','30-34','35-39','40-44')
  base.SE=CItoSE(base)
  lo=log(lb/(1-lb))
  logitmed=log(med/(1-med))-base.SE[1]
  logitse=(log(med/(1-med))-lo)/qnorm(.975)
  data.frame(ages=ages,med=logitmed,se=sqrt(logitse^2+base.SE[2]^2))
}


AgeTable=function(med,lb,effect){
  ages=c('15-19','20-24','25-29','30-34','35-39','40-44')
  lo=log(lb/(1-lb))
  logitmed=log(med/(1-med))+effect[1]
  logitse=(log(med/(1-med))-lo)/qnorm(.975)
  newse=sqrt(logitse^2+effect[2]^2)
  logithi=logitmed+newse*qnorm(.975)
  logitlo=logitmed-newse*qnorm(.975)
  data.frame(ages=ages,lo=expit(logitlo),hi=expit(logithi))
}

AddAgeEffect=function(CITable,AgeEffectTable){
  f=function(x){
    out=adply(CITable,.margins=1,.fun=AddEffect,x)
  }
  adply(AgeEffectTable,.margins=1,.fun=f)
}


samplelogistic=function(x,N=100){
  y=CItoSE(x)
  expit(rnorm(N,mean=y[1],sd=y[2]))
}

combinelogistic.by.year=function(p1,p2,t1,t2,tout){
  s1=logit(samplelogistic(p1))
  s2=logit(samplelogistic(p2))
  samples=expit(s1+(tout-t1)*(s2-s1)/(t2-t1))
  out=quantile(samples,probs=c(.025,.975))
  names(out)=c('lo','hi')
  out
}

###############Botswana##############################
Male15p.pop = 208824 + 434258 + 59399 + 53708

Prevalence.Botswana.2019.15pm = c(lo = 140000, hi = 170000)/Male15p.pop
Incidence.Botswana.2019.15pm = c(lo = 3500, hi = 5000)/Male15p.pop

######################South Africa#######################

#Estimates of 2012 prevalence in cities from 2012 behavioral survey
Prevalence.2012.Capetown=c(lo=3.4, hi=7.8)/100
Prevalence.2012.Johannesburg=c(lo=8.3, hi=14.6)/100
Prevalence.2012.Durban=c(lo=11.2, hi=18.6)/100
Prevalence.2012.Klerksdorp=c(lo=9, hi=19)/100
Prevalence.2012.Tshwane = c(lo = 8.1, hi = 16.6)/100
Prevalence.2012.Buffalo = c(lo = 10.6, hi = 17.3)/100
Prevalence.2012.Manguang = c(lo = 5.3, hi = 11.6)/100
Prevalence.2012.Hlabisa = c(lo=13, hi=15)/100

Prevalence.2012.cities=as.data.frame(rbind(Prevalence.2012.Capetown,Prevalence.2012.Johannesburg,Prevalence.2012.Durban,Prevalence.2012.Klerksdorp))


Incidence.2012.SA.1549m=c(lo=.97,hi=1.45)/100
Incidence.2012.SA.OA=c(lo=.87,hi=1.27)/100

Prevalence.2008.SA.1549 = c(lo=15.5, hi=18.4)/100 #2012 SA survey
Incidence.2008.SA.1549 = c(lo=16.1, hi=17.4)/1000 #UNAIDS

Prevalence.2005.1549.wc.oa=c(lo=1.9, hi=5.3)/100
Prevalence.2005.1549.kz.oa=c(lo=18.3, hi=25.9)/100
Prevalence.2005.1549.nw.oa=c(lo=13.7, hi=23.2)/100
Prevalence.2005.1549.gt.oa=c(lo=13.0, hi=19.1)/100
Prevalence.2005.1549.ec.oa=c(lo=12.1, hi=19.8)/100
Prevalence.2005.1549.fs.oa=c(lo=13.3, hi=26.9)/100

Prevalence.2008.1549.wc.oa=c(lo=3.7, hi=7.5)/100
Prevalence.2008.1549.kz.oa=c(lo=22.1, hi=29.8)/100
Prevalence.2008.1549.nw.oa=c(lo=13.9, hi=22.3)/100
Prevalence.2008.1549.gt.oa=c(lo=12.1, hi=19.0)/100
Prevalence.2008.1549.ec.oa=c(lo=11.9, hi=19.1)/100
Prevalence.2008.1549.fs.oa=c(lo=15.2, hi=22.4)/100

Prevalence.2008.kz.oa = c(lo = 13.4, hi=18.6)/100

Prevalence.2012.1549.wc.oa=c(lo=5.5, hi=10.9)/100
Prevalence.2012.1549.kz.oa=c(lo=25.2, hi=30.8)/100
Prevalence.2012.1549.nw.oa=c(lo=17.5, hi=23.4)/100
Prevalence.2012.1549.gt.oa=c(lo=14.6, hi=21.6)/100
Prevalence.2012.1549.ec.oa=c(lo=17.1, hi=23.0)/100
Prevalence.2012.1549.fs.oa=c(lo=15.4, hi=26.5)/100

Prevalence.2012.wc.oa=c(lo=3.6, hi=7.2)/100
Prevalence.2012.gt.oa=c(lo=10.5, hi=15.5)/100
Prevalence.2012.kz.oa=c(lo=15.8, hi=19.2)/100
Prevalence.2012.ec.oa=c(lo=10.5, hi=14.1)/100
Prevalence.2012.nw.oa=c(lo=12.0, hi=16.1)/100

#South African HIV survey 2017
Prevalence.2017.1549.kz=c(lo=23.9, hi=30.4)/100
Prevalence.2017.1549.wc=c(lo=9.7, hi=16.1)/100
Prevalence.2017.1549.gt=c(lo=14.8, hi=20.7)/100
Prevalence.2017.1549.ec=c(lo=19.8, hi=31.5)/100
Prevalence.2017.1549.nw=c(lo=19.6, hi=26.2)/100



Prevalence.2017.SA.1549m=c(lo=13.3, hi=16.5)/100#South African HIV survey 2017
Prevalence.2017.SA.1549=c(lo=19.2, hi=22.0)/100#South African HIV survey 2017

Incidence.2017.SA.1549m=c(lo=.60,hi=.76)/100

#Comparison of overall prevalence with adult (15-49 year old) men (2012 behavioral survey)
Prevalence.2012.SA.oa=c(lo=11.4,hi=13.1)/100
Prevalence.2012.SA.1549m=c(lo=12.8,hi=16.3)/100
Prevalence.2012.SA.1549=c(lo=17.5,hi=20.3)/100

Incidence.2012.SA.1549=c(lo=1.38, hi=2.06)/100
Incidence.2012.SA.1549m=c(lo=0.97, hi=1.45)/100

#Comparison of overall prevalence with adult (15-49 year old) men (2017 behavioral survey)
Prevalence.2017.SA.oa=c(lo=13.1,hi=15.0)/100
Prevalence.2017.SA.1549m=c(lo=13.3,hi=16.9)/100
Prevalence.2017.SA.1549=c(lo=19.2,hi=22.0)/100

Incidence.2017.SA.1549=c(lo=0.67, hi=0.91)/100
Incidence.2017.SA.1549m=c(lo=0.60, hi=0.76)/100


#UNAIDS
Prevalence.2017.SA.am=c(lo=10.9,hi=15.8)/100
Incidence.2017.SA.1549=c(lo=.876,hi=1.046)/100

Prevalence.2019.SA.1549m=c(lo=9.3,hi=14.3)/100
Incidence.2019.SA.1549=c(lo=.637,hi=0.742)/100


incidence.onepartner=c(lo=1.33,hi=2.01)/100
incidence.multiplepartners=c(lo=1.95,hi=2.91)/100


#SA HIV survey 2017
prev.by.age.2017.SA=c(0.047,0.048,0.124,0.184,0.237,0.224)
prev.by.age.2017.SA.lb=c(0.03,0.03,0.1,0.15,0.2,0.18)


#SA HIV survey 2012
prev.by.age.2012.SA=c(0.007,0.051,0.173,0.256,0.288,0.158)
prev.by.age.2012.SA.lb=c(0.004,0.037,0.138,0.198,0.227,0.118)

inc.by.age.2012.SA=c(.55,.55,1.29,1.29,1.29,1.29)/100
inc.by.age.2012.SA.lb=c(.45,.45,0.91,0.91,0.91,0.91)/100



Capetown.affect.2012=CombineCI(Prevalence.2012.Capetown,Prevalence.2012.SA.oa)
Johannesburg.affect.2012=CombineCI(Prevalence.2012.Johannesburg,Prevalence.2012.SA.oa)
Klerksdorp.affect.2012=CombineCI(Prevalence.2012.Klerksdorp,Prevalence.2012.SA.oa)
Durban.affect.2012=CombineCI(Prevalence.2012.Durban,Prevalence.2012.SA.oa)
Tshwane.affect.2012=CombineCI(Prevalence.2012.Tshwane,Prevalence.2012.SA.oa)
Buffalo.affect.2012=CombineCI(Prevalence.2012.Buffalo,Prevalence.2012.SA.oa)
Kwazulu.affect.2012=CombineCI(Prevalence.2012.kz.oa, Prevalence.2012.SA.oa)
Manguang.affect.2012=CombineCI(Prevalence.2012.Manguang, Prevalence.2012.SA.oa)

westcape.affect.2017.1549=CombineCI(Prevalence.2017.1549.wc,Prevalence.2017.SA.1549)
eastcape.affect.2017.1549=CombineCI(Prevalence.2017.1549.ec,Prevalence.2017.SA.1549)
northwest.affect.2017.1549=CombineCI(Prevalence.2017.1549.nw,Prevalence.2017.SA.1549)
guateng.affect.2017.1549=CombineCI(Prevalence.2017.1549.gt,Prevalence.2017.SA.1549)
Kwazulu.affect.2017.1549=CombineCI(Prevalence.2017.1549.kz,Prevalence.2017.SA.1549)

westcape.affect.2012.oa=CombineCI(Prevalence.2012.wc.oa,Prevalence.2012.SA.oa)
eastcape.affect.2012.oa=CombineCI(Prevalence.2012.ec.oa,Prevalence.2012.SA.oa)
northwest.affect.2012.oa=CombineCI(Prevalence.2012.nw.oa,Prevalence.2012.SA.oa)
guateng.affect.2012.oa=CombineCI(Prevalence.2012.gt.oa,Prevalence.2012.SA.oa)
Kwazulu.affect.2012.oa=CombineCI(Prevalence.2012.kz.oa,Prevalence.2012.SA.oa)

Kwazulu.affect.2008.1549=CombineCI(Prevalence.2008.1549.kz.oa,Prevalence.2008.SA.1549)

Capetown.vs.westcape.2012.oa=CombineCI(Prevalence.2012.Capetown,Prevalence.2012.wc.oa)
Johannesburg.vs.guateng.2012.oa=CombineCI(Prevalence.2012.Johannesburg,Prevalence.2012.gt.oa)
Durban.vs.Kwazulu.2012.oa=CombineCI(Prevalence.2012.Durban,Prevalence.2012.kz.oa)
Hlabisa.vs.Kwazulu.2012.oa=CombineCI(Prevalence.2012.Hlabisa, Prevalence.2012.kz.oa)
Tshwane.vs.guateng.2012.oa=CombineCI(Prevalence.2012.Tshwane,Prevalence.2012.gt.oa)
Buffalo.vs.eastcape.2012.oa=CombineCI(Prevalence.2012.Buffalo,Prevalence.2012.ec.oa)
Klerksdorp.vs.northwest.oa=CombineCI(Prevalence.2012.Klerksdorp,Prevalence.2012.nw.oa)

adultmale.affect.2012.SA=CombineCI(Prevalence.2012.SA.1549m,Prevalence.2012.SA.oa)

multipartnereffect=CombineCI(incidence.multiplepartners,incidence.onepartner)

#Convert prevalence and incidence in 15-49 years to city specific prevalence and incidence in south africa
Prevalence.Capetown.2012.1549m=AddEffect(Prevalence.2012.SA.1549m,Capetown.affect.2012)
Prevalence.Klerksdorp.2012.1549m=AddEffect(Prevalence.2012.SA.1549m,Klerksdorp.affect.2012)
Prevalence.Johannesburg.2012.1549m=AddEffect(Prevalence.2012.SA.1549m,Johannesburg.affect.2012)
Prevalence.Durban.2012.1549m=AddEffect(Prevalence.2012.SA.1549m,Durban.affect.2012)
Prevalence.Tshwane.2012.1549m=AddEffect(Prevalence.2012.SA.1549m,Tshwane.affect.2012)
Prevalence.Manguang.2012.1549m=AddEffect(Prevalence.2012.SA.1549m,Manguang.affect.2012)
Prevalence.Buffalo.2012.1549m=AddEffect(Prevalence.2012.SA.1549m,Buffalo.affect.2012)
Prevalence.Kwazulu.2012.1549m=AddEffect(Prevalence.2012.SA.1549m,Kwazulu.affect.2012)



Incidence.Capetown.2012.1549m=AddEffect(Incidence.2012.SA.1549m,Capetown.affect.2012)
Incidence.Klerksdorp.2012.1549m=AddEffect(Incidence.2012.SA.1549m,Klerksdorp.affect.2012)
Incidence.Johannesburg.2012.1549m=AddEffect(Incidence.2012.SA.1549m,Johannesburg.affect.2012)
Incidence.Durban.2012.1549m=AddEffect(Incidence.2012.SA.1549m,Durban.affect.2012)
Incidence.Tshwane.2012.1549m=AddEffect(Incidence.2012.SA.1549m,Tshwane.affect.2012)
Incidence.Manguang.2012.1549m=AddEffect(Incidence.2012.SA.1549m,Manguang.affect.2012)
Incidence.Buffalo.2012.1549m=AddEffect(Incidence.2012.SA.1549m,Buffalo.affect.2012)
Incidence.Kwazulu.2012.1549m=AddEffect(Incidence.2012.SA.1549m,Kwazulu.affect.2012)

#For 2017 don't have city affects, use province levels and then assume the province to city affect holds true between 2012 and 2017
Prevalence.Capetown.2017.1549m=AddEffect(Prevalence.2017.SA.1549m, CombineSE2(westcape.affect.2017.1549,Capetown.vs.westcape.2012.oa,ADD=TRUE))
Prevalence.Johannesburg.2017.1549m=AddEffect(Prevalence.2017.SA.1549m,CombineSE2(guateng.affect.2017.1549,Johannesburg.vs.guateng.2012.oa,ADD=TRUE))
Prevalence.Durban.2017.1549m=AddEffect(Prevalence.2017.SA.1549m,CombineSE2(Kwazulu.affect.2017.1549,Durban.vs.Kwazulu.2012.oa,ADD=TRUE))
Prevalence.Tshwane.2017.1549m=AddEffect(Prevalence.2017.SA.1549m,CombineSE2(guateng.affect.2017.1549,Tshwane.vs.guateng.2012.oa,ADD=TRUE))
Prevalence.Buffalo.2017.1549m=AddEffect(Prevalence.2017.SA.1549m, CombineSE2(eastcape.affect.2017.1549,Buffalo.vs.eastcape.2012.oa,ADD=TRUE))
Prevalence.Klerksdorp.2017.1549m=AddEffect(Prevalence.2017.SA.1549m,CombineSE2(northwest.affect.2017.1549,Klerksdorp.vs.northwest.oa,ADD=TRUE))
Prevalence.Kwazulu.2017.1549m=AddEffect(Prevalence.2017.SA.1549m,Kwazulu.affect.2017.1549)


Incidence.Capetown.2017.1549m=AddEffect(Incidence.2017.SA.1549m,CombineSE2(westcape.affect.2017.1549,Capetown.vs.westcape.2012.oa,ADD=TRUE))
Incidence.Johannesburg.2017.1549m=AddEffect(Incidence.2017.SA.1549m,CombineSE2(guateng.affect.2017.1549,Johannesburg.vs.guateng.2012.oa,ADD=TRUE))
Incidence.Durban.2017.1549m=AddEffect(Incidence.2017.SA.1549m,CombineSE2(Kwazulu.affect.2017.1549,Durban.vs.Kwazulu.2012.oa,ADD=TRUE))
Incidence.Tshwane.2017.1549m=AddEffect(Incidence.2017.SA.1549m,CombineSE2(guateng.affect.2017.1549,Tshwane.vs.guateng.2012.oa,ADD=TRUE))
Incidence.Buffalo.2017.1549m=AddEffect(Incidence.2017.SA.1549m, CombineSE2(eastcape.affect.2017.1549,Buffalo.vs.eastcape.2012.oa,ADD=TRUE))
Incidence.Klerksdorp.2017.1549m=AddEffect(Incidence.2017.SA.1549m,CombineSE2(northwest.affect.2017.1549,Klerksdorp.vs.northwest.oa,ADD=TRUE))
Incidence.Kwazulu.2017.1549m=AddEffect(Incidence.2017.SA.1549m,Kwazulu.affect.2017.1549)

#Same with 2008


Prevalence.2008.SA.1549m = AddEffect(Prevalence.2008.SA.1549,CombineCI(Prevalence.2012.SA.1549m,Prevalence.2012.SA.1549))
Incidence.2008.SA.1549m = AddEffect(Incidence.2008.SA.1549,CombineCI(Incidence.2012.SA.1549m,Incidence.2012.SA.1549))
Prevalence.Durban.2008.1549m=AddEffect(Prevalence.2008.SA.1549m,CombineSE2(Kwazulu.affect.2008.1549,Durban.vs.Kwazulu.2012.oa,ADD=TRUE))
Prevalence.Hlabisa.2008.1549m=AddEffect(Prevalence.2008.SA.1549m,CombineSE2(Kwazulu.affect.2008.1549,Hlabisa.vs.Kwazulu.2012.oa,ADD=TRUE))
Incidence.Durban.2008.1549m=AddEffect(Incidence.2008.SA.1549m,CombineSE2(Kwazulu.affect.2008.1549,Durban.vs.Kwazulu.2012.oa,ADD=TRUE))
Incidence.Hlabisa.2008.1549m=AddEffect(Incidence.2008.SA.1549m,CombineSE2(Kwazulu.affect.2008.1549,Hlabisa.vs.Kwazulu.2012.oa,ADD=TRUE))



#For years other than 2017 and 2012 we just interpolate
Prevalence.Capetown.2014.1549m=combinelogistic.by.year(Prevalence.Capetown.2012.1549m,Prevalence.Capetown.2017.1549m,2012.5,2017.5,2015)
Prevalence.Johannesburg.2014.1549m=combinelogistic.by.year(Prevalence.Johannesburg.2012.1549m,Prevalence.Johannesburg.2017.1549m,2012.5,2017.5,2015)
Prevalence.Durban.2014.1549m=combinelogistic.by.year(Prevalence.Durban.2012.1549m,Prevalence.Durban.2017.1549m,2012.5,2017.5,2015)

Incidence.Capetown.2014.1549m=combinelogistic.by.year(Incidence.Capetown.2012.1549m,Incidence.Capetown.2017.1549m,2012.5,2017.5,2015)
Incidence.Johannesburg.2014.1549m=combinelogistic.by.year(Incidence.Johannesburg.2012.1549m,Incidence.Johannesburg.2017.1549m,2012.5,2017.5,2015)
Incidence.Durban.2014.1549m=combinelogistic.by.year(Incidence.Durban.2012.1549m,Incidence.Durban.2017.1549m,2012.5,2017.5,2015)

Prevalence.Capetown.2019.1549m=combinelogistic.by.year(Prevalence.Capetown.2012.1549m,Prevalence.Capetown.2017.1549m,2012.5,2017.5,2019.5)
Prevalence.Johannesburg.2019.1549m=combinelogistic.by.year(Prevalence.Johannesburg.2012.1549m,Prevalence.Johannesburg.2017.1549m,2012.5,2017.5,2019.5)
Prevalence.Durban.2019.1549m=combinelogistic.by.year(Prevalence.Durban.2012.1549m,Prevalence.Durban.2017.1549m,2012.5,2017.5,2019.5)

Incidence.Capetown.2019.1549m=combinelogistic.by.year(Incidence.Capetown.2012.1549m,Incidence.Capetown.2017.1549m,2012.5,2017.5,2015)
Incidence.Johannesburg.2019.1549m=combinelogistic.by.year(Incidence.Johannesburg.2012.1549m,Incidence.Johannesburg.2017.1549m,2012.5,2017.5,2019.5)
Incidence.Durban.2019.1549m=combinelogistic.by.year(Incidence.Durban.2012.1549m,Incidence.Durban.2017.1549m,2012.5,2017.5,2019.5)

Prevalence.Sites.SA.2012=as.data.frame(rbind(Prevalence.Capetown.2012.1549m,Prevalence.Klerksdorp.2012.1549m,Prevalence.Johannesburg.2012.1549m,Prevalence.Durban.2012.1549m))
Incidence.Sites.SA.2012=as.data.frame(rbind(Incidence.Capetown.2012.1549m,Incidence.Klerksdorp.2012.1549m,Incidence.Johannesburg.2012.1549m,Incidence.Durban.2012.1549m))
Prevalence.Sites.SA.2017=as.data.frame(rbind(Prevalence.Capetown.2017.1549m,Prevalence.Johannesburg.2017.1549m))
Incidence.Sites.SA.2017=as.data.frame(rbind(Incidence.Capetown.2017.1549m,Incidence.Johannesburg.2017.1549m))

Prevalence.Sites.SA.2012$site=c('Capetown','Klerksdorp','Johannesburg','Durban')
Prevalence.Sites.SA.2017$site=c('Capetown','Johannesburg')
Incidence.Sites.SA.2012$site=c('Capetown','Klerksdorp','Johannesburg','Durban')
Incidence.Sites.SA.2017$site=c('Capetown','Johannesburg')

AgeEffect.SA.Prev.2017=AgeEffectTable(prev.by.age.2017.SA,prev.by.age.2017.SA.lb,Prevalence.2017.SA.1549m)
AgeEffect.SA.Prev.2012=AgeEffectTable(prev.by.age.2012.SA,prev.by.age.2012.SA.lb,Prevalence.2012.SA.1549m)

AgeEffect.SA.Inc.2012=AgeEffectTable(inc.by.age.2012.SA,inc.by.age.2012.SA.lb,Incidence.2012.SA.1549m)

AgeSite.SA.2012.Prev=AddAgeEffect(Prevalence.Sites.SA.2012,AgeEffect.SA.Prev.2012)
AgeSite.SA.2017.Prev=AddAgeEffect(Prevalence.Sites.SA.2017,AgeEffect.SA.Prev.2017)

AgeSite.SA.2017.Inc=AddAgeEffect(Incidence.Sites.SA.2017,AgeEffect.SA.Inc.2012) #Don't have age effects on incidence for 2017, use 2012

AgeSite.SA.2012.Inc=AddAgeEffect(Incidence.Sites.SA.2012,AgeEffect.SA.Inc.2012)


##############Uganda########################
#.oa = adults
#Uganda Phia
Prevalence.2017.Kampala.oa=c(lo=5.6, hi=8.1)/100
Prevalence.2017.Ug.oa=c(lo=5.8, hi=6.7)/100
Prevalence.2017.Ug.1549m=c(lo=3.9, hi=4.7)/100
Incidence.2017.Ug.1549m=c(lo=0.12, hi=0.5)/100

#Uganda progress report 2013
prev.by.age.2011.Ug=c(1.7,2.8,4.0,9.1,11.0,11.3)/100
prev.by.age.2011.Ug.lb=c(1.7,2.8,4.0,9.1,11.0,11.3)/100


inc.by.age.2011.Ug=c(1.7,2.8,4.0,9.1,11.0,11.3)/100
inc.by.age.2011.Ug.lb=c(1.7,2.8,4.0,9.1,11.0,11.3)/100

Prevalence.2011.Kampala.1549m=c(lo=3.6,hi=4.6)/100 #Add some uncertainty as no confidence interval given
Prevalence.2011.Ug.1549m=c(lo=6.1,hi=6.1)/100


#UNAids
Incidence.2011.Ug.oa=c(lo=5.42,hi=6.98)/1000
Incidence.2019.Ug.oa=c(lo=2.08,hi=3.50)/1000

Prevalence.2019.Ug.oa=c(lo=5.4, hi=6.2)/100


Kampala.affect.2017=CombineCI(Prevalence.2017.Kampala.oa,Prevalence.2017.Ug.oa)
affect.1549m.2017=CombineCI(Prevalence.2017.Ug.1549m,Prevalence.2017.Ug.oa)

Prevalence.2017.Kampala.1549m=AddEffect(Prevalence.2017.Ug.1549m,Kampala.affect.2017)

Prevalence.2019.Kampala.1549m=AddEffect(AddEffect(Prevalence.2019.Ug.oa,Kampala.affect.2017),affect.1549m.2017)

Prevalence.2014.Kampala.1549m=combinelogistic.by.year(Prevalence.2017.Kampala.1549m,Prevalence.2011.Kampala.1549m,2017,2011,2015)

#Convert Overall incidence to site/age/gender specific assuming relative risks hold steady between 2011 and 2017
Incidence.2011.Kampala.1549m=AddEffect(Incidence.2011.Ug.oa,CombineSE2(Kampala.affect.2017,affect.1549m.2017,ADD=TRUE))

Incidence.2017.Kampala.1549m=AddEffect(Incidence.2017.Ug.1549m,Kampala.affect.2017)

Incidence.2019.Kampala.1549m=AddEffect(AddEffect(Incidence.2019.Ug.oa,Kampala.affect.2017),affect.1549m.2017)

Incidence.2014.Kampala.1549m=combinelogistic.by.year(Incidence.2017.Kampala.1549m,Incidence.2011.Kampala.1549m,2017,2011,2015)

#Get Age effects
AgeEffects.Ug.Prev=AgeEffectTable(prev.by.age.2011.Ug,prev.by.age.2011.Ug.lb,Prevalence.2011.Ug.1549m)
AgeEffects.Ug.Inc=AgeEffectTable(prev.by.age.2011.Ug,prev.by.age.2011.Ug.lb,Prevalence.2011.Ug.1549m)


Uganda.Sites.Prevalence=as.data.frame(as.list(Prevalence.2011.Kampala.1549m))
Uganda.Sites.Incidence=as.data.frame(as.list(Incidence.2011.Kampala.1549m))

Uganda.Sites.Prevalence$site='Kampala'
Uganda.Sites.Incidence$site='Kampala'

AgeSite.Ug.Prev=AddAgeEffect(Uganda.Sites.Prevalence,AgeEffects.Ug.Prev)
AgeSite.Ug.Inc=AddAgeEffect(Uganda.Sites.Incidence,AgeEffects.Ug.Inc)


############Malawi###############
#Malawi Phia 2016
Incidence.2016.Mal.1549m=c(lo=0.03, hi=0.46)/100
Incidence.2016.Mal.1549=c(lo=0.17, hi=0.49)/100
Incidence.2016.Mal.oa=c(lo=0.20, hi=0.53)/100

Prevalence.2016.Mal.1549m=c(lo=6.8, hi=8.3)/100
Prevalence.2016.Mal.1549=c(lo=9.3, hi=10.7)/100
Prevalence.2016.Mal.oa=c(lo=9.9, hi=11.2)/100

Prevalence.2016.Lilongwe.oa=c(lo=10.4,hi=13.1)/100
Prevalence.2016.Blantyre.oa=c(lo=16.4,hi=19.9)/100


effect.Lilongwe=CombineCI(Prevalence.2016.Lilongwe.oa,Prevalence.2016.Mal.oa)
effect.Blantyre=CombineCI(Prevalence.2016.Blantyre.oa,Prevalence.2016.Mal.oa)

effect.1549m.2016.prev = CombineCI(Prevalence.2016.Mal.1549m,Prevalence.2016.Mal.1549)
effect.1549m.2016.inc = CombineCI(Incidence.2016.Mal.1549m,Incidence.2016.Mal.1549)

#UNAIDS
Prevalence.2007.Mal.1549 = c(lo=9.7, hi=12.5)/100
Incidence.2007.Mal.1549 = c(lo=7.76, hi=8.83)/1000

Prevalence.2019.Mal.1549 = c(lo=7.6, hi=9.6)/100
Incidence.2019.Mal.1549 = c(lo=3.13, hi=4.23)/1000

Prevalence.2007.Lilongwe.1549m=AddEffect(Prevalence.2007.Mal.1549,CombineSE2(effect.Lilongwe,effect.1549m.2016.prev,ADD=TRUE))
Incidence.2007.Lilongwe.1549m=AddEffect(Incidence.2007.Mal.1549,CombineSE2(effect.Lilongwe,effect.1549m.2016.inc,ADD=TRUE))

Prevalence.2007.Blantyre.1549m=AddEffect(Prevalence.2007.Mal.1549,CombineSE2(effect.Blantyre,effect.1549m.2016.prev,ADD=TRUE))
Incidence.2007.Blantyre.1549m=AddEffect(Incidence.2007.Mal.1549,CombineSE2(effect.Blantyre,effect.1549m.2016.inc,ADD=TRUE))

Prevalence.2019.Lilongwe.1549m=AddEffect(Prevalence.2019.Mal.1549,CombineSE2(effect.Lilongwe,effect.1549m.2016.prev,ADD=TRUE))
Incidence.2019.Lilongwe.1549m=AddEffect(Incidence.2019.Mal.1549,CombineSE2(effect.Lilongwe,effect.1549m.2016.inc,ADD=TRUE))

Prevalence.2019.Blantyre.1549m=AddEffect(Prevalence.2019.Mal.1549,CombineSE2(effect.Blantyre,effect.1549m.2016.prev,ADD=TRUE))
Incidence.2019.Blantyre.1549m=AddEffect(Incidence.2019.Mal.1549,CombineSE2(effect.Blantyre,effect.1549m.2016.inc,ADD=TRUE))


Prevalence.2016.Lilongwe.1549m=AddEffect(Prevalence.2016.Mal.1549m,effect.Lilongwe)
Prevalence.2016.Blantyre.1549m=AddEffect(Prevalence.2016.Mal.1549m,effect.Blantyre)

Incidence.2016.Lilongwe.1549m=AddEffect(Incidence.2016.Mal.1549m,effect.Lilongwe)
Incidence.2016.Blantyre.1549m=AddEffect(Incidence.2016.Mal.1549m,effect.Blantyre)


############Zimbabwe###################
#Zimphia 2017
inc.by.age.2016.Zim=c(.14,.14,.48,.48,.38,.38)/100
inc.by.age.2016.Zim.ub=c(.37,.37,1.05,1.05,.91,.91)/100

prev.by.age.2016.Zim=c(.032,.027,.066,.122,.194,.254)
prev.by.age.2016.Zim.N=c(1950,1220,979,942,843,754)
prev.by.age.2016.Zim.lb=prev.by.age.2016.Zim-qnorm(.975)*sqrt(prev.by.age.2016.Zim*(1-prev.by.age.2016.Zim)/prev.by.age.2016.Zim.N)

Prevalence.2016.Zim.1549m=c(lo=10.0,hi=11.5)/100
Prevalence.2016.Zim.1549=c(lo=12.7,hi=13.1)/100

Incidence.2016.Zim.1549m=c(lo=.07,hi=.53)/100
Prevalence.2016.Harare.1549m=c(lo=.08,hi=.125)

Harare.Effect=CombineCI(Prevalence.2016.Harare.1549m,Prevalence.2016.Zim.1549m)
Harare.Effect2=CombineCI(Prevalence.2016.Harare.1549m,Prevalence.2016.Zim.1549)

#UNAIDS
Prevalence.2017.Zim.am=c(lo=8.7,hi=11.9)/100
Incidence.2017.Zim.1549=c(lo=.380,hi=0.692)/100 #Doesn't appear to be huge differences between men and women (UNAIDS)

Prevalence.2007.Zim.1549=c(lo=14.6,hi=19.4)/100
Prevalence.2010.Zim.1549=c(lo=13.1,hi=17.5)/100
Prevalence.2010.Zim.1549m=c(lo=11.5,hi=13)/100 #According to ZimPhia, Harare seems to match average for overall
Prevalence.2019.Zim.1549=c(lo=10.9,hi=14.7)/100

Incidence.2007.Zim.1549=c(lo=7.54, hi=13.97)/1000
Incidence.2010.Zim.oa=c(lo=52,hi=93)/13000
Incidence.2017.Zim.oa=c(lo=.218,hi=.396)/100
Incidence.2019.Zim.1549=c(lo=3.38, hi=6.72)/1000

Prevalence.2007.Harare.1549m=AddEffect(Prevalence.2007.Zim.1549,Harare.Effect2)
Prevalence.2010.Harare.1549m=AddEffect(Prevalence.2010.Zim.1549m,Harare.Effect)
Prevalence.2017.Harare.1549m=AddEffect(Prevalence.2017.Zim.am,Harare.Effect)
Prevalence.2019.Harare.1549m=AddEffect(Prevalence.2019.Zim.1549,Harare.Effect2)

Incidence.2007.Harare.1549m=AddEffect(Incidence.2007.Zim.1549, Harare.Effect2)
Incidence.2010.Harare.1549m=AddEffect(AddEffect(Incidence.2010.Zim.oa,CombineCI(Incidence.2017.Zim.1549,Incidence.2017.Zim.oa)), Harare.Effect) #No adjustment for gender
Incidence.2017.Harare.1549m=AddEffect(Incidence.2017.Zim.1549,Harare.Effect) #Appears to match
Incidence.2019.Harare.1549m=AddEffect(Incidence.2019.Zim.1549,Harare.Effect2) 

Prevalence.2014.Harare.1549m=combinelogistic.by.year(Prevalence.2010.Harare.1549m,Prevalence.2017.Harare.1549m,2010,2017,2015)
Incidence.2014.Harare.1549m=combinelogistic.by.year(Incidence.2010.Harare.1549m,Incidence.2017.Harare.1549m,2010,2017,2015)

AgeEffects.Zim.Prev.2016=AgeEffectTable(prev.by.age.2016.Zim,prev.by.age.2016.Zim.lb,Prevalence.2016.Zim.1549m)
AgeEffects.Zim.Inc.2016=AgeEffectTable(inc.by.age.2016.Zim,inc.by.age.2016.Zim.ub,Incidence.2016.Zim.1549m)

AgeEffects.Zim.Prev.2010=AgeEffect.SA.Prev.2012 #Assume age effects on prevalence change in time, are consistent in region
AgeEffects.Zim.Inc.2010=AgeEffects.Zim.Inc.2016 #Assume age effects on incidence are constant in time

Zimbabwe.Sites.Prevalence.2010=as.data.frame(as.list(Prevalence.2010.Harare.1549m))
Zimbabwe.Sites.Incidence.2010=as.data.frame(as.list(Incidence.2010.Harare.1549m))

Zimbabwe.Sites.Prevalence.2010$site='Harare'
Zimbabwe.Sites.Incidence.2010$site='Harare'

AgeSite.Zim.Prev.2010=AddAgeEffect(Zimbabwe.Sites.Prevalence.2010,AgeEffects.Zim.Prev.2010)
AgeSite.Zim.Inc.2010=AddAgeEffect(Zimbabwe.Sites.Incidence.2010,AgeEffects.Zim.Inc.2010)

Zimbabwe.Sites.Prevalence.2017=as.data.frame(as.list(Prevalence.2017.Harare.1549m))
Zimbabwe.Sites.Incidence.2017=as.data.frame(as.list(Incidence.2017.Harare.1549m))

Zimbabwe.Sites.Prevalence.2017$site='Harare'
Zimbabwe.Sites.Incidence.2017$site='Harare'

AgeSite.Zim.Prev.2017=AddAgeEffect(Zimbabwe.Sites.Prevalence.2017,AgeEffects.Zim.Prev.2016)
AgeSite.Zim.Inc.2017=AddAgeEffect(Zimbabwe.Sites.Incidence.2017,AgeEffects.Zim.Inc.2016)


################Eswatini#######################
#UNAIDS
Incidence.2017.Esw.1549 = c(lo = 12.66, hi = 17.88)/1000
Incidence.2019.Esw.1549 = c(lo = 7.79, hi = 12.44)/1000

Prevalence.2017.Esw.1549 = c(lo = 26.10, hi = 30.0)/100
Prevalence.2019.Esw.1549 = c(lo = 24.6, hi = 28.7)/100

#Eswatini Phia 2018 report covering 2016-17
Prevalence.2017.Esw.15p = c(lo = 25.7, hi=28.3)/100

Prevalence.2017.Esw.1549m = c(lo = 17.3, hi=20.4)/100
Incidence.2017.Esw.1549m = c(lo = 0.21, hi=1.49)/100

Prevalence.2017.Manzini.15p = c(lo = 24.3, hi=30.3)/100
Prevalence.2017.Lubombo.15p = c(lo = 26.7, hi=32.1)/100

Prevalence.2017.Manzini.1549m = AddEffect(Prevalence.2017.Esw.1549m, CombineCI(Prevalence.2017.Manzini.15p, Prevalence.2017.Esw.15p))
Incidence.2017.Manzini.1549m = AddEffect(Incidence.2017.Esw.1549m, CombineCI(Prevalence.2017.Manzini.15p, Prevalence.2017.Esw.15p))

Lubombo.Effect.2017.15p = CombineCI(Prevalence.2017.Lubombo.15p, Prevalence.2017.Esw.15p)
Male.Effect.2017.Esw = CombineCI(Prevalence.2017.Esw.1549m, Prevalence.2017.Esw.15p)

Siteki.male.effect = CombineSE2(Lubombo.Effect.2017.15p, Male.Effect.2017.Esw)
Prevalence.2019.Siteki.1549m = AddEffect(Prevalence.2019.Esw.1549, Siteki.male.effect)
Incidence.2019.Siteki.1549m = AddEffect(Incidence.2019.Esw.1549, Siteki.male.effect)

###############Zambia############################

# Zamphia Data 2019 report covering 2016

Prevalence.2016.Zam.1559 = c(lo = 11.3, hi=12.7)/100

Prevalence.2016.Zam.1549m = c(lo = 7.6, hi=9.0)/100
Prevalence.2016.Zam.1549 = c(lo = 10.7, hi=12.0)/100

Incidence.2016.Zam.1549m = c(lo = 0.07, hi=0.49)/100
Incidence.2016.Zam.1549 = c(lo = 0.42, hi=0.86)/100

Prevalence.2016.Lusaka.1559 = c(lo = 14.1, hi=17.3)/100

Lusaka.effect = CombineCI(Prevalence.2016.Lusaka.1559, Prevalence.2016.Zam.1559)
gender.effect.zambia.prev = CombineCI(Prevalence.2016.Zam.1549m,Prevalence.2016.Zam.1549)
gender.effect.zambia.inc = CombineCI(Incidence.2016.Zam.1549m,Incidence.2016.Zam.1549)
Prevalence.2016.Lusaka.1549m = AddEffect(Prevalence.2016.Zam.1549m, Lusaka.effect)
Incidence.2016.Lusaka.1549m = AddEffect(Incidence.2016.Zam.1549m, Lusaka.effect)

# UNAIDS
Prevalence.2007.Zam.1549 = c(lo = 11.9, hi=15.1)/100
Incidence.2007.Zam.1549 = c(lo = 6.92, hi=14.22)/1000

Prevalence.2007.Lusaka.1549m = AddEffect(Prevalence.2007.Zam.1549, CombineSE2(Lusaka.effect,gender.effect.zambia.prev, ADD = T))
Incidence.2007.Lusaka.1549m = AddEffect(Incidence.2007.Zam.1549, CombineSE2(Lusaka.effect,gender.effect.zambia.inc, ADD = T))

###############Kenya#####################

#No error range given so we will call this the raw value. Error will be interpreted from the national scale.
Prevalence.2017.Kisumu.1549m.raw = c(lo=15.0, hi=15.0)/100
Prevalence.2017.Siaya.1549m.raw = c(lo=19.4, hi=19.4)/100



Prevalence.2017.Kenya.1549 = c(lo=4.02, hi=5.80)/100


Incidence.2017.Kenya.1549 = c(lo=0.13, hi=0.29)/100


#Kenya report 2014
Prevalence.2012.Siaya.1549.raw = c(lo=17.8, hi=17.8)/100
Prevalence.2012.Kenya.1549m = c(lo=4.1, hi=4.1)/100

#UNAID
Prevalence.2010.Kenya.1549 = c(lo=5.1, hi=7.1)/100
Prevalence.2012.Kenya.1549 = c(lo=4.8, hi=6.7)/100
Prevalence.2019.Kenya.1549 = c(lo=4.0, hi=5.2)/100

Incidence.2010.Kenya.1549 = c(lo=1.54, hi=4.64)/1000
Incidence.2012.Kenya.1549 = c(lo=1.34, hi=4.04)/1000
Incidence.2019.Kenya.1549 = c(lo=0.94, hi=2.28)/1000



#Combined gender and regional effect
effect.2017.Kisumu.male.1549 = CombineCI(Prevalence.2017.Kisumu.1549m.raw, Prevalence.2017.Kenya.1549)
effect.2017.siaya.male.1549 = CombineCI(Prevalence.2017.Siaya.1549m.raw, Prevalence.2017.Kenya.1549)

#gender effect
effect.kenya.male.1549 = CombineCI(Prevalence.2012.Kenya.1549m,Prevalence.2012.Kenya.1549)

effect.siaya.1549.2012 = CombineCI(Prevalence.2012.Siaya.1549.raw,Prevalence.2012.Kenya.1549)

effect.2012.siaya.male.1549 = CombineSE2(effect.kenya.male.1549, effect.siaya.1549.2012, ADD = T)

Prevalence.2017.Kisumu.1549m=AddEffect(Prevalence.2017.Kenya.1549, effect.2017.Kisumu.male.1549) #This introduces error bars to the raw value
Incidence.2017.Kisumu.1549m=AddEffect(Incidence.2017.Kenya.1549, effect.2017.Kisumu.male.1549)

Prevalence.2019.Kisumu.1549m=AddEffect(Prevalence.2019.Kenya.1549, effect.2017.Kisumu.male.1549) #This introduces error bars to the raw value
Incidence.2019.Kisumu.1549m=AddEffect(Incidence.2019.Kenya.1549, effect.2017.Kisumu.male.1549)

Prevalence.2017.Siaya.1549m=AddEffect(Prevalence.2017.Kenya.1549, effect.2017.siaya.male.1549)
Incidence.2017.Siaya.1549m=AddEffect(Incidence.2017.Kenya.1549, effect.2017.siaya.male.1549)

Prevalence.2012.Siaya.1549m=AddEffect(Prevalence.2012.Kenya.1549, effect.2012.siaya.male.1549)
Incidence.2012.Siaya.1549m=AddEffect(Incidence.2012.Kenya.1549, effect.2012.siaya.male.1549)

Prevalence.2010.Siaya.1549m=combinelogistic.by.year(Prevalence.2012.Siaya.1549m, Prevalence.2017.Siaya.1549m, 12, 17, 10)
Incidence.2010.Siaya.1549m=combinelogistic.by.year(Incidence.2012.Siaya.1549m, Incidence.2017.Siaya.1549m, 12, 17, 10)

#####USA#########################

#HIV in Philly 2017
Nadultmen.philly = 0.5*(450/29.5-40/11.6)*100000 #Hundreds of thousands
Prevalence.2016.philly = c(lo=19199,hi=19199)
Incidence.2016.philly = c(lo=240,hi=660)

Prevalence.2016.philly.am = Prevalence.2016.philly*0.72*0.48/Nadultmen.philly
Incidence.2016.philly.am = Incidence.2016.philly*0.75*0.9/Nadultmen.philly


source('ArtCoverage.R')
#####Combine Prevalences###########

Prevalence.Voice.by.age=rbind(AgeSite.SA.2012.Prev,AgeSite.Ug.Prev,AgeSite.Zim.Prev.2010)
Incidence.Voice.by.age=rbind(AgeSite.SA.2012.Inc,AgeSite.Ug.Inc,AgeSite.Zim.Inc.2010)

Prevalence.082.by.age=rbind(AgeSite.SA.2017.Prev,AgeSite.Zim.Prev.2017)
Incidence.082.by.age=rbind(AgeSite.SA.2017.Inc,AgeSite.Zim.Inc.2017)

Prevalence.Voice=rbind(Prevalence.Sites.SA.2012,Uganda.Sites.Prevalence,Zimbabwe.Sites.Prevalence.2010)
Incidence.Voice=rbind(Incidence.Sites.SA.2012,Uganda.Sites.Incidence,Zimbabwe.Sites.Incidence.2010)

EpidemicStats.Voice = merge(Prevalence.Voice,Incidence.Voice,by='site',suffixes = c(".prev",".inc"))

EpidemicStats.Voice = merge(EpidemicStats.Voice,ART.VOICE,by='site')

EpidemicStats.Voice$VLS = F
EpidemicStats.Voice$study="VOICE"
EpidemicStats.Voice$country=c("South Africa","South Africa","South Africa","South Africa","Uganda","Zimbabwe")[order(Prevalence.Voice$site)]
EpidemicStats.Voice$N = VOICE.makeup[Prevalence.Voice$site]


Prevalence.082=rbind(Prevalence.Sites.SA.2017,Zimbabwe.Sites.Prevalence.2017)
Incidence.082=rbind(Incidence.Sites.SA.2017,Zimbabwe.Sites.Incidence.2017)




EpidemicStats.082 = merge(Prevalence.082,Incidence.082,by='site',suffixes = c(".prev",".inc"))
EpidemicStats.082$study="HPTN 082"
EpidemicStats.082$country=c("South Africa","South Africa","Zimbabwe")[order(Prevalence.082$site)]
EpidemicStats.082$N = HPTN082.makeup[Prevalence.082$site]

EpidemicStats.082.1 = merge(EpidemicStats.082,ART.HPTN082,by='site')
EpidemicStats.082.2 = merge(EpidemicStats.082,VLS.HPTN082,by='site')

EpidemicStats.082.1$VLS = F
EpidemicStats.082.2$VLS = T

EpidemicStats.082 = rbind(EpidemicStats.082.1,EpidemicStats.082.2)

Prevalence.Aspire = data.frame(rbind(Prevalence.2014.Harare.1549m,Prevalence.2014.Kampala.1549m,
                          Prevalence.2016.Blantyre.1549m,Prevalence.2016.Lilongwe.1549m,
                          Prevalence.Capetown.2014.1549m,Prevalence.Johannesburg.2014.1549m,Prevalence.Durban.2014.1549m))

Incidence.Aspire = data.frame(rbind(Incidence.2014.Harare.1549m,Incidence.2014.Kampala.1549m,
                                    Incidence.2016.Blantyre.1549m,Incidence.2016.Lilongwe.1549m,
                                    Incidence.Capetown.2014.1549m,Incidence.Johannesburg.2014.1549m,Incidence.Durban.2014.1549m))


#ART.Aspire = data.frame(rbind(ART.Aspire.Zimbabwe.AM,ART.Aspire.Uganda.AM,
#                              ART.Aspire.Malawi.AM,ART.Aspire.Malawi.AM,
#                              ART.Aspire.SouthAfrica.AM,ART.Aspire.SouthAfrica.AM,ART.ASPIRE.SouthAfrica.AM))

#Prevalence.Aspire$site = c("Harare","Kampala","Blantyre","Lilongwe","Capetown","Johannesburg","Durban")
#Incidence.Aspire$site = c("Harare","Kampala","Blantyre","Lilongwe","Capetown","Johannesburg","Durban")
#ART.Aspire$site = c("Harare","Kampala","Blantyre","Lilongwe","Capetown","Johannesburg","Durban")

#EpidemicStats.Aspire = merge(Prevalence.Aspire,Incidence.Aspire,by='site',suffixes = c(".prev",".inc"))
#EpidemicStats.Aspire = merge(EpidemicStats.Aspire,ART.Aspire,by='site')

#EpidemicStats.Aspire$country = c("Zimbabwe","Uganda","Malawi","Malawi","South Africa","South Africa","South Africa")[order(Prevalence.Aspire$site)]

#EpidemicStats.Aspire$study = "Aspire"
#EpidemicStats.Aspire$VLS = F
#EpidemicStats.Aspire$N = ASPIRE.makeup[EpidemicStats.Aspire$site]

#Prevalence.Echo = data.frame(rbind(Prevalence.2017.Manzini.1549m, Prevalence.2017.Kisumu.1549m, Prevalence.2016.Lusaka.1549m,
#                        Prevalence.Kwazulu.2017.1549m, Prevalence.Buffalo.2017.1549m,Prevalence.Tshwane.2017.1549m,
#                        Prevalence.Klerksdorp.2017.1549m,Prevalence.Capetown.2017.1549m,Prevalence.Johannesburg.2017.1549m,
#                        Prevalence.Durban.2017.1549m))



#Incidence.Echo = data.frame(rbind(Incidence.2017.Manzini.1549m, Incidence.2017.Kisumu.1549m, Incidence.2016.Lusaka.1549m,
#                                  Incidence.Kwazulu.2017.1549m, Incidence.Buffalo.2017.1549m,Incidence.Tshwane.2017.1549m,
#                                  Incidence.Klerksdorp.2017.1549m,Incidence.Capetown.2017.1549m,Incidence.Johannesburg.2017.1549m,
#                                  Incidence.Durban.2017.1549m))

#VLS.Echo = data.frame(rbind(VLS.Esw.15pm.2017,ART.Kisumu.2017.am, VLS.Zam.1559m.2016, NEWVLS.Kwazulu.AM,NEWVLS.EastCape.AM,
#                            NEWVLS.Guateng.AM, NEWVLS.Northwest.AM, NEWVLS.WestCape.AM, NEWVLS.Guateng.AM, NEWVLS.Kwazulu.AM))

# Prevalence.Echo$site = c("Manzini","Kisumu","Lusaka","Kwazulu","Buffalo","Tshwane","Klerksdorp","Capetown","Johannesburg","Durban")
# Incidence.Echo$site = c("Manzini","Kisumu","Lusaka","Kwazulu","Buffalo","Tshwane","Klerksdorp","Capetown","Johannesburg","Durban")
# VLS.Echo$site = c("Manzini","Kisumu","Lusaka","Kwazulu","Buffalo","Tshwane","Klerksdorp","Capetown","Johannesburg","Durban")
# 
# EpidemicStats.ECHO = merge(Prevalence.Echo,Incidence.Echo,by='site',suffixes = c(".prev",".inc"))
# EpidemicStats.ECHO = merge(EpidemicStats.ECHO,VLS.Echo,by='site')
# 
# EpidemicStats.ECHO$country = c("Eswatini","Kenya","Zambia","South Africa","South Africa","South Africa","South Africa","South Africa","South Africa","South Africa")[order(Prevalence.Echo$site)]
# 
# EpidemicStats.ECHO$study="ECHO"
# EpidemicStats.ECHO$VLS=c(T,F,T,T,T,T,T,T,T,T)[order(Prevalence.Echo$site)]
# EpidemicStats.ECHO$N = ECHO.makeup[EpidemicStats.ECHO$site]
# 
# 
# Prevalence.FEMPREP = data.frame(rbind(Prevalence.2010.Siaya.1549m,Prevalence.Tshwane.2012.1549m,Prevalence.Manguang.2012.1549m))
# Incidence.FEMPREP = data.frame(rbind(Incidence.2010.Siaya.1549m,Incidence.Tshwane.2012.1549m,Incidence.Manguang.2012.1549m))
# #ART.FEMPREP = data.frame(rbind(ART.FEMPREP.Kenya.AM, OLDART.SouthAfrica.AM, OLDART.SouthAfrica.AM))
# 
# Prevalence.FEMPREP$site = c("Bondo","Tshwane","Manguang")
# Incidence.FEMPREP$site = c("Bondo","Tshwane","Manguang")
# ART.FEMPREP$site = c("Bondo","Tshwane","Manguang")
# 
# EpidemicStats.FEMPREP = merge(Prevalence.FEMPREP,Incidence.FEMPREP,by = 'site',suffixes = c(".prev",".inc"))
# EpidemicStats.FEMPREP = merge(EpidemicStats.FEMPREP,ART.FEMPREP,by='site')
# 
# EpidemicStats.FEMPREP$site = c("Bondo","Tshwane","Manguang")
# EpidemicStats.FEMPREP$country = c("Kenya","South Africa","South Africa")[order(Prevalence.FEMPREP$site)]
# EpidemicStats.FEMPREP$N = FEMPREP.makeup[EpidemicStats.FEMPREP$site]
# EpidemicStats.FEMPREP$study = "FEMPREP"
# EpidemicStats.FEMPREP$VLS = F
# 
# Prevalence.035 = data.frame(rbind(Prevalence.2007.Blantyre.1549m,Prevalence.2007.Lilongwe.1549m,
#                                   Prevalence.2007.Lusaka.1549m,Prevalence.2007.Harare.1549m,
#                                   Prevalence.Hlabisa.2008.1549m,Prevalence.Durban.2008.1549m))
# Incidence.035 = data.frame(rbind(Incidence.2007.Blantyre.1549m,Incidence.2007.Lilongwe.1549m,
#                                   Incidence.2007.Lusaka.1549m,Incidence.2007.Harare.1549m,
#                                   Incidence.Hlabisa.2008.1549m,Incidence.Durban.2008.1549m))
# 
# #ART.035 = data.frame(rbind(ART.035.Malawi.AM, ART.035.Malawi.AM, ART.035.Zambia.AM, ART.035.Zimbabwe.AM,
# #                           ART.035.SouthAfrica.AM, ART.035.SouthAfrica.AM))
# 
# 
# Prevalence.035$site = c("Blantyre","Lilongwe","Lusaka","Harare","Hlabisa","Durban")
# Incidence.035$site = c("Blantyre","Lilongwe","Lusaka","Harare","Hlabisa","Durban")
# ART.035$site = c("Blantyre","Lilongwe","Lusaka","Harare","Hlabisa","Durban")
# 
# EpidemicStats.035 = merge(Prevalence.035,Incidence.035,by = 'site',suffixes = c(".prev",".inc"))
# EpidemicStats.035 = merge(EpidemicStats.035,ART.035,by = 'site')
# EpidemicStats.035$country = c("Malawi","Malawi","Zambia","Zimbabwe","South Africa","South Africa")[order(Prevalence.035$site)]
# EpidemicStats.035$N = HPTN035.makeup[EpidemicStats.035$site]
# EpidemicStats.035$study = "HPTN 035"
# EpidemicStats.035$VLS = F




combinelogistic.by.site=function(DF,Nruns=1000){
  K=ncol(DF)
  f1 = function(x){
    with(as.list(x),{
      samplelogistic(c(lo=lo.prev,hi=hi.prev),N = Nruns)
    })
  }
  f2 = function(x){
    with(as.list(x),{
      samplelogistic(c(lo=lo.inc,hi=hi.inc),N = Nruns)
    })
  }
  f3 = function(x){
    with(as.list(x),{
      samplelogistic(c(lo=lo,hi=hi),N = Nruns)
    })
  }
  DF.sample1=adply(DF,.margins = 1,.fun=f1)
  DF.sample2=adply(DF,.margins = 1,.fun=f2)
  DF.sample3=adply(DF,.margins = 1,.fun=f3)
  w=as.numeric(DF$N)
  w=w/sum(w)
  samples1=w%*%as.matrix(DF.sample1[,K+1:Nruns])
  samples2=w%*%as.matrix(DF.sample2[,K+1:Nruns])
  samples3=w%*%as.matrix(DF.sample3[,K+1:Nruns])
  out1=quantile(samples1,probs=c(.025,0.5,.975))
  out2=quantile(samples2,probs=c(.025,0.5,.975))
  out3=quantile(samples3,probs=c(.025,0.5,.975))
  out=c(out1,out2,out3)
  names(out)=c('lo.prev','med.prev','hi.prev','lo.inc','med.inc','hi.inc','lo.tr','med.tr','hi.tr')
  out
}

combinelogistic.by.site.prevonly=function(DF,Nruns=1000){
  K=ncol(DF)
  f1 = function(x){
    with(as.list(x),{
      samplelogistic(c(lo=lo,hi=hi),N = Nruns)
    })
  }

  DF.sample1=adply(DF,.margins = 1,.fun=f1)
  w=as.numeric(DF$N)
  w=w/sum(w)
  samples1=w%*%as.matrix(DF.sample1[,K+1:Nruns])

  out=quantile(samples1,probs=c(.025,.975))
  names(out)=c('lo','hi')
  out
}

combinelogistic.by.siteinc=function(DF,DFinc,N=1000){
  K=ncol(DF)
  DF.sample=adply(DF,.margins = 1,.fun=samplelogistic, N=N)
  DFinc.sample=adply(DFinc,.margins = 1,.fun=samplelogistic, N=N)
  w=as.numeric(W[DF.sample$site])
  M1=w*as.matrix(DFinc.sample[,K+1:N])
  samples=colSums((M1)*as.matrix(DF.sample[,K+1:N]))/colSums(M1)
  out=quantile(samples,probs=c(.025,.975))
  names(out)=c('lo','hi')
  out
}

combinelogistic.by.age=function(DF,W,N=1000){
  K=ncol(DF)
  DF.sample=adply(DF,.margins = 1,.fun=samplelogistic, N=N)
  w=W/sum(W)
  samples=w%*%as.matrix(DF.sample[,K+1:N])
  out=quantile(samples,probs=c(.025,.975))
  names(out)=c('lo','hi')
  out
}


#combine.by.age=function(DF,W){
#  DF$N=W
#  ddply(DF,.variables = 'ages',combinelogistic.by.site.prevonly)
#}


#EpidemicStats = rbind(EpidemicStats.Voice, EpidemicStats.035,EpidemicStats.082,EpidemicStats.ECHO,EpidemicStats.Aspire,EpidemicStats.FEMPREP)

#Incidence.Voice.Combine=combinelogistic.by.site(Incidence.Voice)

#Prevalence.Aspire.Combine = combinelogistic.by.site(Prevalence.Aspire,ASPIRE.makeup)
#Incidence.Aspire.Combine = combinelogistic.by.site(Incidence.Aspire,ASPIRE.makeup)

#Prevalence.Echo.Combine = combinelogistic.by.site(Prevalence.Echo,ECHO.makeup)
#Incidence.Echo.Combine = combinelogistic.by.site(Incidence.Echo,ECHO.makeup)

#Prevalence.Voice.Combine.by.age=combine.by.age(Prevalence.Voice.by.age,VOICE.makeup)
#Incidence.Voice.Combine.by.age=combine.by.age(Incidence.Voice.by.age,VOICE.makeup)

#Prevalence.HPTN082.Combine=combinelogistic.by.site(Prevalence.082,HPTN082.makeup)
#Incidence.HPTN082.Combine=combinelogistic.by.site(Incidence.082,HPTN082.makeup)

#Prevalence.HPTN082.Combine.by.age=combine.by.age(Prevalence.082.by.age,HPTN082.makeup)
#Incidence.HPTN082.Combine.by.age=combine.by.age(Incidence.082.by.age,HPTN082.makeup)

#EpidemicStats.VOICE.combine = combinelogistic.by.site(subset(EpidemicStats, study=="VOICE"))
#EpidemicStats.082.combine = combinelogistic.by.site(subset(EpidemicStats, study=="HPTN 082"))

#Prevalence.HPTN082.mainpartners=c(EpidemicStats.082.combine['lo.prev'],EpidemicStats.082.combine['hi.prev'])
#names(Prevalence.HPTN082.mainpartners) = c('lo','hi')
#Prevalence.HPTN082.casualpartners=AddEffect(Prevalence.HPTN082.mainpartners,multipartnereffect)
#Incidence.HPTN082.mainpartners=c(EpidemicStats.082.combine['lo.inc'],EpidemicStats.082.combine['hi.inc'])
#names(Incidence.HPTN082.mainpartners) = c('lo','hi')

#Prevalence.HPTN082.mainpartners.ageadjust=combinelogistic.by.age(Prevalence.HPTN082.Combine.by.age,N.age.082)
#Incidence.HPTN082.mainpartners.ageadjust=combinelogistic.by.age(Incidence.HPTN082.Combine.by.age,N.age.082)

#Prevalence.VOICE.mainpartners=c(EpidemicStats.VOICE.combine['lo.prev'],EpidemicStats.VOICE.combine['hi.prev'])
#names(Prevalence.VOICE.mainpartners) = c('lo','hi')
#Prevalence.VOICE.casualpartners=AddEffect(Prevalence.VOICE.mainpartners,multipartnereffect)
#Incidence.VOICE.mainpartners=c(EpidemicStats.VOICE.combine['lo.inc'],EpidemicStats.VOICE.combine['hi.inc'])
#names(Incidence.VOICE.mainpartners) = c('lo','hi')

#EpidemicStats.VOICE.SA = combinelogistic.by.site(subset(EpidemicStats, study=="VOICE"&country=="South Africa"))
#Prevalence.VOICE.SA = c(EpidemicStats.VOICE.SA['lo.prev'],EpidemicStats.VOICE.SA['hi.prev'])
#names(Prevalence.VOICE.SA) = c('lo','hi')
#Prevalence.VOICE.SA.casualpartners=AddEffect(Prevalence.VOICE.SA,multipartnereffect)

#write.csv(EpidemicStats, 'EpidemicStats.csv', row.names = FALSE)