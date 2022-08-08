#These are all the functions necessary to run the markov model (as it stands now)+ some unnecessary ones


#Naming convention
#A=transition matrix. Describes the transitions between different numbers of partners, and their HIV status
#B=sexual frequency matrix. It's a matrix because it depends on today's status and yesterday's status
#R=per act risk matrix

#These matrices are carried around in a giant list called 'fT'. The 'with' command adds the components of fT to the name-space

#This first function simulates infection and transition each day. It takes a shortcut as the process is time-independent
InfectionPredict4=function(Ntimes,fT,eff,params){
  #Include the model parameters (params) and the matrices (fT) in the namespace
  with(as.list(c(params,fT)),{ 
    a=a0 #Initial condition
    P=A*B*R*(1-condomfrequency*condomefficacy)*(1-eff) #Daily risk of infection
    Q=t(A-P) #Update matrix is transition minus infection
    Qpow=Q%^%(Ntimes)
    eye=diag(rep(1,length(diag(A)))) #Identity matrix the same size as A
    
    #This is equivalent to sum_i^Ntimes Q^i a0
    y=solve(Q-eye, Qpow%*%a0-a0)
    numfollowups=sum(y)
    numinfections = sum(t(P)%*%y)
    #numsexA=sum(t(B*AP)%*%y)
    #numsexV=sum(t(B*(1-AP))%*%y)
    list(numinfections,numfollowups)
  })
}

#In this version, the efficacy can change every day. Marginally slower, but more readable
EfficacyPredict=function(Ntimes,fT,effA,effV,params){
  with(as.list(c(params,fT)),{
    a=a0
    PA=A*B*RA*AP*(1-condomfrequency*condomefficacy) #Daily infection rate due to anal sex
    PV=A*B*RV*(1-AP)*(1-condomfrequency*condomefficacy) #Daily infection rate due to vaginal sex
    numinfections=0
    numfollowups=0
    numsexA=0
    numsexV=0
    for(i in seq(1,Ntimes)){
      P.i=PA*(1-effA[i])+PV*(1-effV[i])
      numfollowups=numfollowups+1-numinfections
      
      newa=a%*%A
      
      infections = a%*%P.i
      
      numsexA=a%*%(B*AP)
      numsexV=a%*%(B*(1-AP))
      a=newa-infections #Remove infected individuals from trial population
      numinfections=numinfections+sum(infections)
    }
    list(numinfections,numfollowups)
  })
}

#This version is currently not used. But it uses the sex acts to inform which state an individual is in (ie from 067 data).
#It uses the 'forward algorithm' from markov theory
ForwardProb=function(Ntimes,fT,effA,effV,params,sexA,sexV){
  with(as.list(c(params,fT)),{
  a=a0
  for(i in seq(Days)){
    sex.i=subset(S,day==i)
    anal.i=sex.i$anal
    kv=sexV[i]
    ka=sexA[i]
    #Get sex rate
    Psex=B^k*exp(-B)/factorial(kv+ka)
    
    #Vector of probabilities of having that combination of anal/vaginal sex depending on state of partnerships
    Prob.anal=analprob^(ka)*(1-AP)^(kv)*factorial(k)/(factorial(ka)*factorial(kv))
    Psex=Psex*Prob.anal
    
    I=(1-(1-effV[i]*RA)^(ka)*(1-effA[i]*RV)^(kv))*(1-condomfrequency*condomefficacy)
    
    #Get bayesian update matrix
    C=Psex*A
    p=p%*%C
  }
  list(numinfections,numfollowups)})
}

#This is a fast version of the infection probability (equivalent to what Daniel did)

InfectionPredict2alt=function(fT,eff,params){
  with(as.list(params),{
    A=fT$A
    B=fT$B
    R=fT$R
    a0=fT$a0
    AP=fT$AP
    #a=a0
    M=(A*B*R)*(1-condomfrequency*condomefficacy)*(1-eff)
    #Q=A-M
    #Qpow=Q%^%(Ntimes)
    #eye=diag(rep(1,length(diag(A)))) #Identity matrix the same size as A
    #y=solve(Q-eye, Qpow%*%a0-a0)
    numfollowups=1
    numinfections=sum(a0%*%M)
    sexrateA=sum(a0%*%(A*B*AP))
    sexrateV=sum(a0%*%(A*B*(1-AP)))
    list(numinfections,numfollowups,sexrateA,sexrateV)
  })
}



#A slightly simpler/slower version of InfectionPredict4 (include to make sure they get the same answer)
InfectionPredict3=function(Ntimes,fT,eff,params){
  with(as.list(c(params,fT)),{
    a=a0
    P=A*B*R*(1-condomfrequency*condomefficacy)*(1-eff)
    numinfections=0
    numfollowups=0
    numpills=0
    numsexA=0
    numsexV=0
    for(i in seq(1,Ntimes)){
      numfollowups=numfollowups+1-numinfections
      
      newa=a%*%A
      
      infections = a%*%P
      
      numsexA=a%*%(B*AP)
      numsexV=a%*%(B*(1-AP))
      a=newa-infections #Remove infected individuals from trial population
      numinfections=numinfections+sum(infections)
      numpills=numpills+sum(numpills)
    }
    list(numinfections,numfollowups)
  })
}

#Wrapper function for all transition matrices. This generates the list fT and needs to be called prior to simulation
fullTransition=function(params){
  with(as.list(params),{
    
    
    #condomfrequency=.45
    
    A=getA(params) #Transition matrix (same for low/casual risk)
    B.main=getB(params) #Sex rate with main partner
    ap.main=getanalprob(params) #Anal fraction with main partner
    
    BA.main=B.main*ap.main #Anal with main
    BV.main=B.main*ap.main #Vaginal with main
    
    R.main=getR(params) #Per act risk with main partners
    riskcasual=getRiskCasual(params) #Per act risk with casual
    a0.full=geta0(params) #Initial condition
    
    B=B.main+sexcasual #Total sex rates
    BA=BA.main+sexcasual*anal.casual
    BV=BV.main+sexcasual*(1-anal.casual)
    
    AP=ifelse(B>0,BA/B,0) #Anal frequency
    RA=ifelse(BA>0,(R.main*BA.main+sexcasual*anal.casual*riskcasual)/BA,0) #Per act risk from anal
    RV=ifelse(BV>0,(R.main*BA.main+sexcasual*(1-anal.casual)*riskcasual)/BV,0) #Per act risk from vaginal
    
    #RV=(R.main*BV.main+sexcasual*(1-anal.casual)*riskcasual)/BV
    #RA=(R.main*BA.main+sexcasual*anal.casual*riskcasual)/BA
    #RA[is.infinite(RA)]=0
    
    R.casual=AP*RA*sexRiskA+(1-AP)*RV*sexRiskV #Per act risk total
    R.low=(sexRiskA*ap.main+sexRiskV*(1-ap.main))*R.main
    
    #E.main=getE()
    #EA=E.main*BA.main+anal.casual*sexcasual
    #EV=E.main*BV.main+(1-anal.casual)*sexcasual
    
    fT.L=list(A=A,B=B.main,R=R.low,AP=ap.main,a0=a0.full,RA=R.main*sexRiskA,RV=R.main*sexRiskV)#,EA=E.main,EV=E.main)
    
    fT.C=list(A=A,B=B,R=R.casual,AP=AP,a0=a0.full,RA=RA*sexRiskA,RV=RV*sexRiskV)#,EA=EA,EV=EV)
    list(fT.L,fT.C)
  })
}
#Below are the functions that actually build all the matrices. 
#These are internal functions that should not necessarily be interacted with.

#Helper functions for markov model

NHIV=10 #Number of HIV statuses (useful to have as global variable)
#1=uninfected, 2=Acute,3=Asymptomatic, 4=Asymptomatic/ART early,5=Asymptomatic/ART Suppressed, 6=Asymptomatic/ART Unsuppressed 7=Late, 8=Late/ART early,   9=Late/ART Suppressed 10=Late/ART Unsuppressed,

Nstatus=NHIV*3 #Multiply by three to get all partnership types (anal only, vaginal only, anal/vaginal)
Nstates=2*Nstatus+1 #Total number of states a particpant can be in, accounting for  short/long term and also no partner

index.ONART=c(4,5,6,8,9,10)-1 #Offset by one because uninfected doesn't count
index.VLS=c(5,9)-1 #Offset by one because uninfected doesn't count


#Proportion of HIV statuses (including negative) in main/casual partners
getS0=function(params,main=TRUE){
  with(as.list(params),{
    hiv.prop=getHIVprop(params)
    s0 = c(1-prevalence.casual,prevalence.casual*hiv.prop)
    if(main){
      s0 = c(1-prevalence,prevalence*hiv.prop)
    }
    s0
  })
}

#Proportion of main partners in each possible class (HIV status/sex type)
getfullS0=function(params){
  with(as.list(params),{
    s0=getS0(params)
    c(s0*analonly,s0*vaginalonly,s0*analvaginal)
  })
}

#Progression rates from one HIV phase to the next (length 10)
#Progression is both the transition from the acute->asymptomatic->late, but alsothe rate the
#individuals become unsuppressed while on ART.
getHIVprogression=function(params){
  with(as.list(params),{
    #artunsuppressionrate=(1-1/artsuppressed)*artsuppressionrate
    hiv.progression=c(1/lengthacute,1/lengthasymptomatic,0,artunsuppressionrate,0,0,0,artunsuppressionrate,0)
    c(incidence,hiv.progression)
  })
}

#Art initiation rates in each HIV phase (length 10)
getARTinitiation=function(params){
  with(as.list(params),{
    c(0,0,artinitasymptomatic,artsuppressionrate,0,artsuppressionrate,artinitlate,artsuppressionrate,0,artsuppressionrate)
  })
}

#Death rates of partners
getDeath=function(params){
  with(as.list(params),{
    c(0,0,0,0,0,0,1/lengthlate,1/lengthlate,0,1/lengthlate)
  })
}

#Rate of a partner quitting ART
getARTQuit=function(params){
  with(as.list(params),{
    c(0,0,0,artquitrateearly,artquitrate,artquitrate,0,artquitrateearly,artquitrate,artquitrate)
  })
}

#Generate the risk matrix (per act risk from main partners)
getR=function(params){
  with(as.list(params),{
    HIVRisk=getHIVRisk(params)
    rep(1,Nstates)%*%t(c(0,rep(HIVRisk,6)))
  })
}

#Transition probabilities between each state
getA=function(params){
  with(as.list(params),{
    
    #New partners sampled from vector S0
    s0=getfullS0(params)
    
    #transition matrix HIV.mat says how HIV status of partners changes
    HIV.mat=getHIVMat(params)
    
    #Initialize matrix with zeros
    A = matrix(0,nrow = Nstates,ncol=Nstates)
    
    #Set diagonal entries to one
    diag(A)=1
    
    #Breakups
    dissolution=rep(c(dis.s,dis.l),each=Nstatus)
    
    #Transition from short to long
    transition=rep(c(trans,0),each=Nstatus)
    
    #Fillout entries
    #A[start,end]=rate, A[start,start]=-rate
    
    A[1,1+1:Nstatus]=acq*s0 #Acquisition of partners (always start at 1)
    diag(A)=diag(A)-c(acq,dissolution+transition)
    A[,1]=A[,1]+c(0,dissolution)
    
    #Transition from short to long represents a jump of Nstatus indices
    for(i in seq(Nstatus)){
      A[i+1,i+1+Nstatus]=A[i+1,i+1+Nstatus]+transition[i]
    }
    
    #Multiply each 32 by 32 block by the HIV transition matrix
    for(i in seq(0,1)){
      for(j in seq(0,1)){
        A[1+Nstatus*i+1:Nstatus,1+Nstatus*j+1:Nstatus]=HIV.mat%*%A[1+Nstatus*i+1:Nstatus,1+Nstatus*j+1:Nstatus]   
      }
    }
    
    #Enforce rows adding to one (take into account death of partner)
    A[,1]=A[,1]+1-rowSums(A)
    A
  })
}

#Return the per act risk by HIV status
getHIVRisk=function(params){
  with(as.list(params),{
    suppressed=c(0,0,0,0,1,0,0,0,1,0)
    unsuppressed=c(0,0,0,1,0,1,0,1,0,1)
    artrisk=(1-artefficacy*unsuppressed-suppressed)
    phaserisk=c(0,relriskacute,1,1,1,1,relrisklate,relrisklate,relrisklate,relrisklate)
    phaserisk*artrisk
  })
}

#HIV transition matrix
getHIVMat=function(params,simple=FALSE){
  with(as.list(params),{
    #Rate of progression
    progression = getHIVprogression(params)
    
    #Rate of art initiation
    art.initiation = getARTinitiation(params)
    
    art.quit=getARTQuit(params)
    
    aids.death=getDeath(params)
    
    HIV.mat=matrix(0,nrow=NHIV,ncol=NHIV)
    diag(HIV.mat)=1-progression-art.initiation-aids.death-art.quit
    
    progressdest = c(2,3,7,4,6,6,7,8,10,10)
    artdest=c(1,2,4,5,5,5,8,9,9,9)
    quitdest=c(1,2,3,3,3,3,7,7,7,7)
    
    for(i in seq(NHIV)){
      HIV.mat[i,progressdest[i]]=HIV.mat[i,progressdest[i]]+progression[i]
      HIV.mat[i,artdest[i]]=HIV.mat[i,artdest[i]]+art.initiation[i]
      HIV.mat[i,quitdest[i]]=HIV.mat[i,quitdest[i]]+art.quit[i]
    }
    
    if(simple){return(HIV.mat)}
    zeros=matrix(0,nrow=NHIV,ncol=NHIV)
    
    rbind(cbind(HIV.mat,zeros,zeros),cbind(zeros,HIV.mat,zeros),cbind(zeros,zeros,HIV.mat))
  })
}


getHIVprop=function(p){
  with(as.list(p),{
    k1=artinitasymptomatic/(artsuppressionrate+artquitrateearly)
    k2=artunsuppressionrate/(artsuppressionrate+artquitrate)
    k3=artinitlate/(artsuppressionrate+artquitrateearly+1/lengthlate)
    k4=artunsuppressionrate/(artsuppressionrate+artquitrate+1/lengthlate)
    k5=artsuppressionrate/((1+k4)*artquitrate+k4/lengthlate)
    
    
    
    N1=lengthacute
    N2=lengthasymptomatic
    
    T2=k1*N2
    S2=artsuppressionrate*T2/((1+k2)*artquitrate)
    U2=k2*S2
    
    N3=1/(1/lengthlate+artinitlate-artquitrate*k3*k5*(1+k4)-artquitrateearly*k3)
    T3=k3*N3
    
    
    S3=k5*T3
    U3=k4*S3
    
    y=c(N1,N2,T2,S2,U2,N3,T3,S3,U3)
    y/sum(y)
  })
}

getARTprop=function(params){
  hivprop=getHIVprop(params)
  sum(hivprop[index.ONART])
}

getVLSprop=function(params){
  hivprop=getHIVprop(params)
  sum(hivprop[index.VLS])
}

#Sex rates
getB=function(params){
  with(as.list(params),{
    sex=c(0,sexrate,sexrate)
    
    #relative sex rate for each concurrency
    sexadj=rep(1,3*NHIV)
    
    sex.full=c(sex[1],rep(sex[2:3],each=Nstatus)*rep(sexadj,2))
    
    #account for lower sex rates with extended partners
    #sexmatrix (take into account always having sex during initiation)
    B = rep(1,Nstates)%*%t(sex.full)
    
    B[1,1+1:Nstatus]=1
    B
  })
}

#Fraction of acts which are anal
getanalprob=function(params){
  with(as.list(params),{
    analprob=c(0,rep(c(1,0,analfrac),2,each=NHIV))
    rep(1,Nstates)%*%t(analprob)
  })
}

#Initial condition (steady state)
geta0=function(params){
  with(as.list(params),{
    s0 = getS0(params)
    k1=acq/(trans+dis.s)
    k2=trans/dis.l
    none=1/(1+k1+k1*k2)
    short=k1*none
    long=k2*short
    if(cohabit){
      tmp=short+long
      short=short/tmp
      long=long/tmp
      none=0
    }
    a0=c(none,short*c(analonly,vaginalonly,analvaginal),long*c(analonly,vaginalonly,analvaginal))
    c(a0[1],rep(s0,6)*rep(a0[2:7],each=length(s0)))
  })
}

#Risk from casual partners
getRiskCasual=function(params){
  with(as.list(params),{
    HIVRisk=getHIVRisk(params)
    s0=getS0(params,main=FALSE)
    sum(s0*HIVRisk)
  })
}

#Risk from main partner
getRiskMain=function(params){
  with(as.list(params),{
    HIVRisk=getHIVRisk(params)
    s0=getS0(params)
    sum(s0*HIVRisk)+incidence*relriskacute
  })
}

#Per exposure risk
getPerExposureRisk=function(params){
  with(as.list(params),{
    HIVRisk=getHIVRisk(params)
    HIVprop=getHIVprop(params)
    sum(HIVprop*HIVRisk[2:NHIV])
  })
}

#Rate of exposure
getE=function(){
  exposure=rep(c(0,rep(1,NHIV-1)),6)
  rep(1,Nstates)%*%t(c(0,exposure))
}

