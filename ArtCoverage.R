#UNAIDS 2018

#2017 numbers
NEWART.SouthAfrica.AM=c(lo=0.48,hi=.58)
NEWART.Uganda.AM=c(lo=.58,hi=0.68)
NEWART.Zimbabwe.AM=c(lo=.66,hi=.86)
NEWART.Malawi.AM=c(lo=0.56,hi=0.68)

ART.Ken.2010.AM=c(lo=.18,hi=.29)#VOICE should be average of 2011 and 2010
ART.Ken.2011.AM=c(lo=.21,hi=.34)





ART.Zam.2010.OA=c(lo=.31,hi=.39)#VOICE should be average of 2011 and 2010
ART.Zam.2011.OA=c(lo=.36,hi=.46)


#########################Zimbabwe###################################
#UNAIDS
ART.Zim.2010.OA=c(lo=.26,hi=.34)#VOICE should be average of 2011 and 2010
ART.Zim.2011.OA=c(lo=.34,hi=.45)
ART.Zim.2012.OA=c(lo=.40,hi=.51)
ART.Zim.2013.OA=c(lo=.46,hi=.59)
ART.Zim.2014.OA=c(lo=.53,hi=.69)
ART.Zim.2015.OA=c(lo=.59,hi=.76)
ART.Zim.2016.OA=c(lo=.65,hi=.84)
ART.Zim.2017.OA=c(lo=.74,hi=.95)#HPTN082 should use 2017

ART.Zim.by.year=rbind(ART.Zim.2010.OA,ART.Zim.2011.OA,ART.Zim.2012.OA,
                      ART.Zim.2013.OA,ART.Zim.2014.OA,ART.Zim.2015.OA,ART.Zim.2016.OA,ART.Zim.2017.OA)


VLS.Zim.2019.AM = c(lo = .59, hi = 0.79)

#ZimPhia 2016
VLS.Harare.1564.2016=c(lo=.518,hi=.649)
VLS.Zimbabwe.1564.2016=c(lo=.583,hi=.625)
VLS.Zimbabwe.1564m.2016=c(lo=.511,hi=.576)

Harare.Effect.2016 = CombineCI(VLS.Harare.1564.2016, VLS.Zimbabwe.1564.2016)

VLS.Harare.2019.AM = AddEffect(VLS.Zim.2019.AM, Harare.effect.2016)


OLDART.Zimbabwe.OA=(ART.Zim.2010.OA+ART.Zim.2011.OA)/2 #2011+2010
NEWART.Zimbabwe.OA=c(lo=.74,hi=.95) #2017

ART.Aspire.Zimbabwe.OA=(ART.Zim.2012.OA+ART.Zim.2013.OA+ART.Zim.2014.OA+ART.Zim.2015.OA+ART.Zim.2016.OA+ART.Zim.2017.OA)/6

ART.035.Zimbabwe.OA = (combinelogistic.by.year(ART.Zim.2010.OA, ART.Zim.2011.OA,2010,2011,2006)+
                         combinelogistic.by.year(ART.Zim.2010.OA, ART.Zim.2011.OA,2010,2011,2007)+
                         combinelogistic.by.year(ART.Zim.2010.OA, ART.Zim.2011.OA,2010,2011,2008))/3

VLS.Harare.AM.2016=AddEffect(VLS.Harare.1564.2016,CombineCI(VLS.Zimbabwe.1564m.2016,VLS.Zimbabwe.1564.2016))
NEWVLS.Harare.AM=AddEffect(VLS.Harare.AM.2016,CombineCI(ART.Zim.2017.OA,ART.Zim.2016.OA))

VLS.Harare.AM.2019 = AddEffect(VLS.Zim.2019.AM, Harare.Effect.2016)
######################################################################################

#####################Botswana###############################################
VLS.Botswana.15p.2019 = c(lo = 61, hi = 74) / 100

########################################################


#############Uganda#######################################################
#UNAIDS
ART.Ug.2010.OA=c(lo=.20,hi=.22)#VOICE should be average of 2011 and 2010
ART.Ug.2011.OA=c(lo=.22,hi=.25)
ART.Ug.2012.OA=c(lo=.30,hi=.33)
ART.Ug.2013.OA=c(lo=.39,hi=.44)
ART.Ug.2014.OA=c(lo=.48,hi=.54)
ART.Ug.2015.OA=c(lo=.52,hi=.59)
ART.Ug.2016.OA=c(lo=.59,hi=.67)
ART.Ug.2017.OA=c(lo=.68,hi=.77)#HPTN082 should use 2017

ART.Ug.by.year=rbind(ART.Ug.2010.OA,ART.Ug.2011.OA,ART.Ug.2012.OA,
                     ART.Ug.2013.OA,ART.Ug.2014.OA,ART.Ug.2015.OA,ART.Ug.2016.OA,ART.Ug.2017.OA)

VLS.Ug.2019.AM = c(lo = 53, hi = 77) / 100
##########################################################################



ART.SA.2010.OA=c(lo=.22,hi=.27)#VOICE should be average of 2011 and 2010
ART.SA.2011.OA=c(lo=.29,hi=.36)
ART.SA.2012.OA=c(lo=.34,hi=.42)
ART.SA.2013.OA=c(lo=.38,hi=.47)
ART.SA.2014.OA=c(lo=.43,hi=.52)
ART.SA.2015.OA=c(lo=.47,hi=.57)
ART.SA.2016.OA=c(lo=.51,hi=.61)
ART.SA.2017.OA=c(lo=.56,hi=.66)#HPTN082 should use 2017

ART.SA.by.year=rbind(ART.SA.2010.OA,ART.SA.2011.OA,ART.SA.2012.OA,
                     ART.SA.2013.OA,ART.SA.2014.OA,ART.SA.2015.OA,ART.SA.2016.OA,ART.SA.2017.OA)




#South Africa 2017 survey
NEWVLS.EastCape.OA=c(lo=.680,hi=.680)
NEWVLS.Kwazulu.OA=c(lo=.675,hi=.675)
NEWVLS.Northwest.OA=c(lo=.579,hi=.579)
NEWVLS.Guateng.OA=c(lo=.569,hi=.569)
NEWVLS.WestCape.OA=c(lo=.547,hi=.547)
NEWVLS.SouthAfrica.1549m=c(lo=.451,hi=.564)
NEWVLS.SouthAfrica.OA=c(lo=.590,hi=.65) #Derived from 62.3 overall VLS, but add some uncertainty as 95%CI not reported.


#UNAIDS
VLS.SouthAfrica.2019.15pm = c(lo = 51, hi = 62)/ 100

NEWVLS.Johannesburg.AM=AddEffect(NEWVLS.Guateng.OA,CombineCI(NEWVLS.SouthAfrica.1549m,NEWVLS.SouthAfrica.OA))
NEWVLS.Capetown.AM=AddEffect(NEWVLS.WestCape.OA,CombineCI(NEWVLS.SouthAfrica.1549m,NEWVLS.SouthAfrica.OA))
NEWVLS.Kwazulu.AM=AddEffect(NEWVLS.Kwazulu.OA,CombineCI(NEWVLS.SouthAfrica.1549m,NEWVLS.SouthAfrica.OA))
NEWVLS.EastCape.AM=AddEffect(NEWVLS.EastCape.OA,CombineCI(NEWVLS.SouthAfrica.1549m,NEWVLS.SouthAfrica.OA))
NEWVLS.WestCape.AM=AddEffect(NEWVLS.WestCape.OA,CombineCI(NEWVLS.SouthAfrica.1549m,NEWVLS.SouthAfrica.OA))
NEWVLS.Northwest.AM=AddEffect(NEWVLS.Northwest.OA,CombineCI(NEWVLS.SouthAfrica.1549m,NEWVLS.SouthAfrica.OA))
NEWVLS.Guateng.AM=AddEffect(NEWVLS.Guateng.OA,CombineCI(NEWVLS.SouthAfrica.1549m,NEWVLS.SouthAfrica.OA))

VLS.Johannesburg.AM=AddEffect(NEWVLS.Guateng.OA,CombineCI(VLS.SouthAfrica.2019.15pm,NEWVLS.SouthAfrica.OA))
VLS.Capetown.AM=AddEffect(NEWVLS.WestCape.OA,CombineCI(VLS.SouthAfrica.2019.15pm,NEWVLS.SouthAfrica.OA))
VLS.Kwazulu.AM=AddEffect(NEWVLS.Kwazulu.OA,CombineCI(VLS.SouthAfrica.2019.15pm,NEWVLS.SouthAfrica.OA))
VLS.EastCape.AM=AddEffect(NEWVLS.EastCape.OA,CombineCI(VLS.SouthAfrica.2019.15pm,NEWVLS.SouthAfrica.OA))
VLS.WestCape.AM=AddEffect(NEWVLS.WestCape.OA,CombineCI(VLS.SouthAfrica.2019.15pm,NEWVLS.SouthAfrica.OA))
VLS.Northwest.AM=AddEffect(NEWVLS.Northwest.OA,CombineCI(VLS.SouthAfrica.2019.15pm,NEWVLS.SouthAfrica.OA))
VLS.Guateng.AM=AddEffect(NEWVLS.Guateng.OA,CombineCI(VLS.SouthAfrica.2019.15pm,NEWVLS.SouthAfrica.OA))


OLDART.Uganda.OA=(ART.Ug.2010.OA+ART.Ug.2011.OA)/2 #2011
NEWART.Uganda.OA=c(lo=.68,hi=.77) #2017

ART.Aspire.Uganda.OA=(ART.Ug.2012.OA+ART.Ug.2013.OA+ART.Ug.2014.OA+ART.Ug.2015.OA+ART.Ug.2016.OA+ART.Ug.2017.OA)/6

OLDART.SouthAfrica.OA=(ART.SA.2010.OA+ART.SA.2011.OA)/2 #2011+2010
NEWART.SouthAfrica.OA=c(lo=.56,hi=.66) #2017

ART.Aspire.SouthAfrica.OA=(ART.SA.2012.OA+ART.SA.2013.OA+ART.SA.2014.OA+ART.SA.2015.OA+ART.SA.2016.OA+ART.SA.2017.OA)/6

ART.035.SouthAfrica.OA = (combinelogistic.by.year(ART.SA.2010.OA, ART.SA.2011.OA,2010,2011,2006)+
                       combinelogistic.by.year(ART.SA.2010.OA, ART.SA.2011.OA,2010,2011,2007)+
                       combinelogistic.by.year(ART.SA.2010.OA, ART.SA.2011.OA,2010,2011,2008))/3


############################Malawi############################################################
#UNAIDS
ART.Mal.2010.OA=c(lo=.25,hi=.32)#VOICE should be average of 2011 and 2010
ART.Mal.2011.OA=c(lo=.32,hi=.40)
ART.Mal.2012.OA=c(lo=.39,hi=.49)
ART.Mal.2013.OA=c(lo=.45,hi=.55)
ART.Mal.2014.OA=c(lo=.49,hi=.61)
ART.Mal.2015.OA=c(lo=.53,hi=.65)
ART.Mal.2016.OA=c(lo=.60,hi=.73)
ART.Mal.2017.OA=c(lo=.65,hi=.79)#HPTN082 should use 2017

VLS.Mal.2019.AM = c(lo = 0.56, hi = 0.67)

#Malawi PHIA 2015
VLS.Mal.2015.OA = c(lo = 66.0, hi = 70.7) / 100
VLS.Lilongwe.2015.OA = c(lo = 59.3, hi = 70.4) / 100
VLS.Blantyre.2015.OA = c(lo = 53.9, hi = 65.0) / 100


Lilongwe.effect.2015.OA = CombineCI(VLS.Lilongwe.2015.OA, VLS.Mal.2015.OA)
Blantyre.effect.2015.OA = CombineCI(VLS.Blantyre.2015.OA, VLS.Mal.2015.OA)

VLS.Lilongwe.2019.AM = AddEffect(VLS.Mal.2019.AM, Lilongwe.effect.2015.OA)
VLS.Blantyre.2019.AM = AddEffect(VLS.Mal.2019.AM, Blantyre.effect.2015.OA)

NEWART.Malawi.OA=ART.Mal.2017.OA

ART.Aspire.Malawi.OA=(ART.Mal.2012.OA+ART.Mal.2013.OA+ART.Mal.2014.OA+ART.Mal.2015.OA+ART.Mal.2016.OA+ART.Mal.2017.OA)/6


ART.035.Malawi.OA = (combinelogistic.by.year(ART.Mal.2010.OA, ART.Mal.2011.OA,2010,2011,2006)+
                       combinelogistic.by.year(ART.Mal.2010.OA, ART.Mal.2011.OA,2010,2011,2007)+
                       combinelogistic.by.year(ART.Mal.2010.OA, ART.Mal.2011.OA,2010,2011,2008))/3

###############Eswatini##############

VLS.Esw.15pm.2019 = c(lo = 80, hi = 96)/100

VLS.Esw.15p.2017 = c(lo=71.3, hi=75.0)/100
VLS.Esw.15pm.2017 = c(lo=64.5, hi=70.6)/100

VLS.Manzini.15p.2017 = c(lo = 68.2, hi=74.7)/100
VLS.Siteki.15p.2017 = c(lo = 68.7, hi=76.5)/100
VLS.Manzini.15pm.2017 = AddEffect(VLS.Manzini.15p.2017, CombineCI(VLS.Esw.15pm.2017, VLS.Esw.15p.2017))

VLS.Siteki.15pm.2019 = AddEffect(VLS.Siteki.15p.2017, CombineCI(VLS.Esw.15pm.2019, VLS.Esw.15p.2017))

##############Zambia###############
NEWART.Zambia.OA = c(lo = 0.54, hi = 0.58) #Derived from 56.2% with 2434 respondents using binom.test
NEWART.Zambia.AM = c(lo = 0.46, hi = 0.53) #Derived from 49.3% and 693 respondents using binom.test

VLS.Zam.1559.2016 = c(lo=56.6, hi=61.7)/100
VLS.Zam.1559m.2016 = c(lo=53.1, hi=61.3)/100

VLS.Lusaka.1559.2016 = c(lo = 58.0, hi=67.3)/100
VLS.Lusaka.1559m.2016 = AddEffect(VLS.Lusaka.1559.2016, CombineCI(VLS.Zam.1559m.2016, VLS.Zam.1559.2016))


ART.035.Zambia.OA = (combinelogistic.by.year(ART.Zam.2010.OA, ART.Zam.2011.OA,2010,2011,2006)+
                         combinelogistic.by.year(ART.Zam.2010.OA, ART.Zam.2011.OA,2010,2011,2007)+
                         combinelogistic.by.year(ART.Zam.2010.OA, ART.Zam.2011.OA,2010,2011,2008))/3

#############Kenya#####################

ART.Kenya.2017.a = c(lo=0.73, hi=0.73)
ART.Kenya.2011.am = c(lo=0.31, hi=0.31)

ART.Kenya.2017.am = c(lo = 0.52, hi = 0.74)

ART.Kisumu.2017.a = c(lo=0.90,hi=0.90)

ART.Kisumu.2017.am = AddEffect(ART.Kisumu.2017.a, CombineCI(ART.Kenya.2017.am, ART.Kenya.2017.a))

#PHIA Prelimnary report (2018? appears as 2019 in UN data book)
VLS.Kenya.2019.1549m = c(lo = 53.6, hi = 67.7) / 100
VLS.Kenya.2019.1564 = c(lo = 68.8, hi = 74.4) / 100
VLS.Kisumu.2019.1564 = c(lo = 75.5, hi = 90.9) / 100

VLS.Kisumu.2019.1549m = AddEffect(VLS.Kisumu.2019.1564, CombineCI(VLS.Kenya.2019.1549m, VLS.Kenya.2019.1564))



#Harare very similar to Zimbabwe as a whole (according to Zimphia)
#Kampala very similar to Uganda as a whole (according to Uganda Phia)

OLDART.SouthAfrica.AM=AddEffect(NEWART.SouthAfrica.AM,CombineCI(OLDART.SouthAfrica.OA,NEWART.SouthAfrica.OA))
OLDART.Zimbabwe.AM=AddEffect(NEWART.Zimbabwe.AM,CombineCI(OLDART.Zimbabwe.OA,NEWART.Zimbabwe.OA))
OLDART.Uganda.AM=AddEffect(NEWART.Uganda.AM,CombineCI(OLDART.Uganda.OA,NEWART.Uganda.OA))

ART.FEMPREP.Kenya.AM = (ART.Ken.2010.AM+ART.Ken.2011.AM)/2

ART.035.SouthAfrica.AM=AddEffect(ART.035.SouthAfrica.OA,CombineCI(NEWART.SouthAfrica.AM,NEWART.SouthAfrica.OA))
ART.035.Zambia.AM=AddEffect(ART.035.Zambia.OA,CombineCI(NEWART.Zambia.AM,NEWART.Zambia.OA))
ART.035.Zimbabwe.AM=AddEffect(ART.035.Zimbabwe.OA,CombineCI(NEWART.Zimbabwe.AM,NEWART.Zimbabwe.OA))
ART.035.Malawi.AM=AddEffect(ART.035.Malawi.OA,CombineCI(NEWART.Malawi.AM,NEWART.Malawi.OA))


ART.Aspire.SouthAfrica.AM=AddEffect(ART.Aspire.SouthAfrica.OA,CombineCI(NEWART.SouthAfrica.AM,NEWART.SouthAfrica.OA))
ART.Aspire.Uganda.AM=AddEffect(ART.Aspire.Uganda.OA,CombineCI(NEWART.Uganda.AM,NEWART.Uganda.OA))
ART.Aspire.Zimbabwe.AM=AddEffect(ART.Aspire.Zimbabwe.OA,CombineCI(NEWART.Zimbabwe.AM,NEWART.Zimbabwe.OA))
ART.Aspire.Malawi.AM=AddEffect(ART.Aspire.Malawi.OA,CombineCI(NEWART.Malawi.AM,NEWART.Malawi.OA))





ART.VOICE=as.data.frame(rbind(OLDART.SouthAfrica.AM,OLDART.Zimbabwe.AM,OLDART.Uganda.AM))
ART.HPTN082=as.data.frame(rbind(NEWART.SouthAfrica.AM,NEWART.Zimbabwe.AM,NEWART.Uganda.AM))
VLS.HPTN082=as.data.frame(rbind(NEWVLS.Harare.AM,NEWVLS.Capetown.AM,NEWVLS.Johannesburg.AM))

ART.VOICE=ART.VOICE[c(1,1,1,2,3,1),]
ART.VOICE$site=names(VOICE.makeup)

VLS.HPTN082$site=names(HPTN082.makeup)

ART.HPTN082=ART.HPTN082[c(2,1,1),]
ART.HPTN082$site=names(HPTN082.makeup)



