#Risk ratios from VOICE data
riskratio.low=c(AG=1.28,MP=1.18,DR=1.07,FN=1.04,SE=1.03,ST=1.14)
riskratio.med=c(AG=1.70,MP=1.80, DR=1.41,FN=1.38, SE=1.63,ST=1.49)
riskratio.high=c(AG=2.27,MP=2.75,DR=1.87,FN=1.83,SE=2.58,ST=1.93)

#Fractions with each risk factor in VOICE cohort
riskfrac=c(AG=0.51,MP=0.68,DR=0.26,FN=0.17,SE=0.75,ST=0.2)
target.incidence=6.05
target.incidence.stderr=0.78

riskscore=seq(0,9)
HIVinfections=c(1,1,4,3,15,27,45,72,46,49)
TotalWomen=c(125,27,532,196,879,548,911,787,443,386)

HIVRisk.df=data.frame(riskscore=riskscore,HIVinfection=HIVinfections,TotalWomen=TotalWomen)

hrsums5=colSums(subset(HIVRisk.df,riskscore>5))

VOICESERO=263
VOICEPERSONYEARS=4384

VOICESERO.SA=254
VOICEPERSONYEARS.SA=3458

VOICESERO.ZIM=3
VOICEPERSONYEARS.ZIM=600

VOICESERO.UG=6
VOICEPERSONYEARS.UG=290

VOICEMINPROB=dpois(VOICESERO,lambda=VOICESERO,log=TRUE)

#This actually computes the likelihood (the main work of the calibration function)

riskratio.se=(log(riskratio.med)-log(riskratio.low))/qnorm(.975)


calibrationtarget.value=log(riskratio.med)[components]
calibrationtarget.stderr=riskratio.se[components]