#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(expm,plyr)
source("InfectionModel.R")
source("RiskAdjustmentFunctions.R")
source('HPTN082Import.R')
source('VOICEimport.R')
source("Prevalence.R")
source("ArtCoverage.R")
source("siteadjustments.R")
library(shiny, shinyjs)

# Define UI for application that draws a histogram
ui <- fluidPage(
   
  shinyjs::useShinyjs(),
   # Application title
   titlePanel("Predict HIV incidence in female population"),
   
   # Sidebar with a slider input for number of bins 
   sidebarLayout(sidebarPanel(
       
       
       selectInput("MalePOP", label = "HIV Epidemic Parameters", choices = c("Aspire:Lilongwe",
                                                                     "Aspire:Blantyre",
                                                                     "Aspire:Harare",
                                                                     "Aspire:Kampala",
                                                                     "Aspire:Capetown",
                                                                     "Aspire:Durban",
                                                                     "Aspire:Johannesberg",
                                                                     "VOICE:Harare",
                                                                     "VOICE:Kampala",
                                                                     "VOICE:Klerkdorp",
                                                                     "VOICE:Durban",
                                                                     "VOICE:Johannesberg",
                                                                     "HPTN 082:Capetown",
                                                                     "HPTN 082:Johannesberg",
                                                                     "HPTN 082:Harare"),
                   selected = "Aspire:Lilongwe"),
       
       checkboxInput("MalePOPcustom", "Use Custom Epidemic Setting", value = FALSE),
       sliderInput("PREV", 
                   label = "Male HIV Prevalence",
                   min = .1, max = 99.9, value = c(2.0, 9.9), round=-1, step = .1),
       
       sliderInput("INC", 
                   label = "Male HIV Incidence",
                   min = .01, max = 10, value = as.numeric(round(Incidence.2016.Lilongwe.1549m,2)), round=-2, step = .01),
       
       selectInput("VLSorART", label = "Choose Treatment Cascade Metric", choices = c("Specify Viral Suppression","Specify ART coverage"), selected = "Specify ART coverage"),
       sliderInput("VLS", 
                   label = "Male HIV Viral Suppression",
                   min = .1, max = 99.9, value = c(45, 60), round=-1, step = .1),
       
       selectInput("FEMPOP", label = "Risk Factors of Study Population", choices = c("ASPIRE Trial", "VOICE Trial", "HPTN 082"), selected = "ASPIRE Trial"),
       
       checkboxInput("FEMPOPcustom", "Use Custom Study Population", value = FALSE),
       sliderInput("MPpop", 
                   label = "Does not live with main partner",
                   min = 0, max = 100, value = 58),
       
       sliderInput("FNpop", 
                   label = "Does not receive financial assistance from partner",
                   min = 0, max = 100, value = 46),
       
       sliderInput("DRpop", 
                   label = "Alcohol use within three months",
                   min = 0, max = 100, value = 26),
       
       sliderInput("STpop", 
                   label = "STI at enrollment",
                   min = 0, max = 100, value = 21),
       
       sliderInput("SEpop", 
                   label = "Does not know if partner has other partners or not.",
                   min = 0, max = 100, value = 57),
       
       sliderInput("AGpop", 
                   label = "Age 25 and younger",
                   min = 0, max = 100, value = 39)),
       
     
     
     mainPanel(sliderInput("NRUNS", 
                                   label = "Number of Simulations",
                                   min = 10, max = 200, value = 100,step = 10),
               fluidRow(
                 actionButton("GO", label="Go"),
                actionButton("PAUSE", label="Pause"),
                actionButton("RESET", label="Reset"),
                plotOutput("progress", width = 100, height = 20)
                ),
               plotOutput("hist"))
     )

)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  NP=nrow(params.voice.orig)
  bar.aspect.ratio=20
  incidence <<- NULL
  rv=reactiveValues(computing=F,prior.runs=0, reset=F, runs.completed=0)
  shinyjs::disable("PAUSE")
  shinyjs::disable("RESET")
  shinyjs::disable("PREV")
  shinyjs::disable("INC")
  shinyjs::disable("VLS")
  shinyjs::disable("MPpop")
  shinyjs::disable("FNpop")
  shinyjs::disable("DRpop")
  shinyjs::disable("AGpop")
  shinyjs::disable("SEpop")
  shinyjs::disable("STpop")
   observeEvent(input$PAUSE,{
     computing<-isolate(rv$computing)
     rv$computing<- !computing
     shinyjs::toggleState("RESET")
     updateActionButton(session,"PAUSE",label=ifelse(computing,"Resume","Pause"))
   })
   observeEvent(input$RESET,{
     incidence <<- NULL
     shinyjs::enable("PREV")
     shinyjs::enable("INC")
     shinyjs::enable("VLS")
     shinyjs::enable("VLSorART")
     shinyjs::enable("MPpop")
     shinyjs::enable("FNpop")
     shinyjs::enable("DRpop")
     shinyjs::enable("AGpop")
     shinyjs::enable("SEpop")
     shinyjs::enable("STpop")
     shinyjs::enable("GO")
     shinyjs::enable("FEMpop")
     shinyjs::enable("Malepop")
     shinyjs::enable("FEMpopcustom")
     shinyjs::enable("Malepopcustom")
     updateActionButton(session,"PAUSE",label = "Pause")
     rv$reset=T
     rv$runs.completed=0
     rv$prior.runs=0
   })
   
   observeEvent(input$VLSorART,{
     newlabel = ifelse(input$VLSorART=="Specificy Viral Suppression","Male HIV Viral Suppression","Male ART coverage")
     updateSliderInput(session, "VLS", label = newlabel, min = 0.1)
   })
   observeEvent(input$FEMPOP,{
     riskfrac=riskfrac.ASPIRE*100
     if(input$FEMPOP=="VOICE Trial"){
       riskfrac=riskfrac.VOICE*100
     }
     
     if(input$FEMPOP=="HPTN 082"){
       riskfrac=riskfrac.082*100
     }
     
     if(input$FEMPOP=="ASPIRE Trial"){
       riskfrac=riskfrac.ASPIRE*100
     }
     
     updateSliderInput(session, "AGpop", value=as.numeric(riskfrac['AG']))
     updateSliderInput(session, "DRpop", value=as.numeric(riskfrac['DR']))
     updateSliderInput(session, "FNpop", value=as.numeric(riskfrac['FN']))
     updateSliderInput(session, "SEpop", value=as.numeric(riskfrac['SE']))
     updateSliderInput(session, "STpop", value=as.numeric(riskfrac['ST']))
     updateSliderInput(session, "MPpop", value=as.numeric(riskfrac['MP']))
   })
   
   
   observeEvent(input$MalePOP,{

     
     prev = Prevalence.2016.Lilongwe.1549m
     inc = Incidence.2016.Lilongwe.1549m
     vls = ART.Aspire.Malawi.AM
     vlsorart = F
     
     if(input$MalePOP=="Aspire:Lilongwe"){
       prev = Prevalence.2016.Lilongwe.1549m
       inc = Incidence.2016.Lilongwe.1549m
       vls = ART.Aspire.Malawi.AM
     }
     if(input$MalePOP=="Aspire:Blantyre"){
       prev = Prevalence.2016.Blantyre.1549m
       inc = Incidence.2016.Blantyre.1549m
       vls = ART.Aspire.Malawi.AM
     }
     if(input$MalePOP=="Aspire:Harare"){
       prev = Prevalence.2014.Harare.1549m
       inc = Incidence.2014.Harare.1549m
       vls = ART.Aspire.Zimbabwe.AM
     }
     if(input$MalePOP=="Aspire:Kampala"){
       prev = Prevalence.2014.Kampala.1549m
       inc = Incidence.2014.Kampala.1549m
       vls = ART.Aspire.Uganda.AM
     }
     if(input$MalePOP=="Aspire:Capetown"){
       prev = Prevalence.Capetown.2014.1549m
       inc = Incidence.Capetown.2014.1549m
       vls = ART.Aspire.SouthAfrica.AM
     }
     if(input$MalePOP=="Aspire:Durban"){
       prev = Prevalence.Durban.2014.1549m
       inc = Incidence.Durban.2014.1549m
       vls = ART.Aspire.SouthAfrica.AM
     }
     if(input$MalePOP=="Aspire:Johannesberg"){
       prev = Prevalence.Johannesburg.2014.1549m
       inc = Incidence.Johannesburg.2014.1549m
       vls = ART.Aspire.SouthAfrica.AM
     }
     if(input$MalePOP=="VOICE:Harare"){
       prev = Prevalence.2010.Harare.1549m
       inc = Incidence.2010.Harare.1549m
       vls = OLDART.Zimbabwe.AM
     }
     if(input$MalePOP=="VOICE:Kampala"){
       prev = Prevalence.2011.Kampala.1549m
       inc = Incidence.2011.Kampala.1549m
       vls = OLDART.Uganda.AM
     }
     if(input$MalePOP=="VOICE:Klerksdorp"){
       prev = Prevalence.Klerksdorp.2012.1549m
       inc = Incidence.Klerksdorp.2012.1549m
       vls = OLDART.SouthAfrica.AM
     }
     if(input$MalePOP=="VOICE:Durban"){
       prev = Prevalence.Durban.2012.1549m
       inc = Incidence.Durban.2012.1549m
       vls = OLDART.SouthAfrica.AM
     }
     if(input$MalePOP=="VOICE:Johannesberg"){
       prev = Prevalence.Johannesburg.2012.1549m
       inc = Incidence.Johannesburg.2012.1549m
       vls = OLDART.SouthAfrica.AM
     }
     if(input$MalePOP=="HPTN 082:Capetown"){
       prev = Prevalence.Capetown.2017.1549m
       inc = Incidence.Capetown.2017.1549m
       vls = NEWVLS.Capetown.AM
       vlsorart = TRUE
     }
     if(input$MalePOP=="HPTN 082:Johannesberg"){
       prev = Prevalence.Johannesburg.2017.1549m
       inc = Incidence.Johannesburg.2017.1549m
       vls = NEWVLS.Johannesburg.AM
       vlsorart = TRUE
     }
     if(input$MalePOP=="HPTN 082:Harare"){
       prev = Prevalence.2017.Harare.1549m
       inc = Incidence.2017.Harare.1549m
       vls = NEWVLS.Harare.AM
       vlsorart = TRUE
     }
     updateSliderInput(session, "PREV", value = 100*as.numeric(prev))
     updateSliderInput(session, "INC", value = 100*as.numeric(inc))
     updateSliderInput(session, "VLS", value = 100*as.numeric(vls))
     updateCheckboxInput(session, "VLSorART", value = ifelse(vlsorart,"Specify Viral Suppression","Specify ART coverage"))
   })
   
   observeEvent(input$FEMPOPcustom,{
     if(!input$FEMPOPcustom){
       shinyjs::disable("MPpop")
       shinyjs::disable("FNpop")
       shinyjs::disable("DRpop")
       shinyjs::disable("AGpop")
       shinyjs::disable("SEpop")
       shinyjs::disable("STpop")
       shinyjs::show("FEMPOP")
     }
     else{
       shinyjs::enable("MPpop")
       shinyjs::enable("FNpop")
       shinyjs::enable("DRpop")
       shinyjs::enable("AGpop")
       shinyjs::enable("SEpop")
       shinyjs::enable("STpop")
       shinyjs::hide("FEMPOP")
     }
   })
   
   observeEvent(input$MalePOPcustom,{
     if(!input$MalePOPcustom){
       shinyjs::disable("PREV")
       shinyjs::disable("INC")
       shinyjs::disable("VLS")
       shinyjs::show("MalePOP")
     }
     else{
       shinyjs::enable("PREV")
       shinyjs::enable("INC")
       shinyjs::enable("VLS")
       shinyjs::hide("MalePOP")
     }
   })
    observeEvent(input$GO,{
      #rep(NA,NP)
      rv$computing <- T
      shinyjs::disable("GO")
      shinyjs::disable("RESET")
      shinyjs::enable("PAUSE")
      shinyjs::disable("NRUNS")
      shinyjs::disable("PREV")
      shinyjs::disable("INC")
      shinyjs::disable("VLS")
      shinyjs::disable("VLSorART")
      shinyjs::disable("MPpop")
      shinyjs::disable("FNpop")
      shinyjs::disable("DRpop")
      shinyjs::disable("AGpop")
      shinyjs::disable("SEpop")
      shinyjs::disable("STpop")
      shinyjs::disable("FEMpop")
      shinyjs::disable("Malepop")
      shinyjs::disable("FEMpopcustom")
      shinyjs::disable("Malepopcustom")
      rv$prior.runs = rv$prior.runs + rv$runs.completed
      rv$runs.completed=0
      output$hist <- renderPlot({
        computing=rv$computing
        reset=rv$reset
        
        rv$reset=F
        
        
        riskfrac.POP=isolate(c(MP=input$MPpop, FN=input$FNpop, DR=input$DRpop, ST=input$STpop, SE=input$SEpop, AG=input$AGpop)/100)
        
        VLS.pop=isolate(c(lo=input$VLS[1], hi=input$VLS[2])/100)
        INC.pop=isolate(c(lo=input$INC[1], hi=input$INC[2])/100)
        PREV.pop=isolate(c(lo=input$PREV[1], hi=input$PREV[2])/100)
        VLSORART = isolate(input$VLSorART == "Specify Viral Suppression")
        
        
        if(computing){
          i=sample(NP, 1)
          out.incidence=site.adjust.params(params.voice.orig[i,],
                                         VLS.pop, 
                                         INC.pop, 
                                         PREV.pop, 
                                         riskfrac.POP,
                                         VLS = VLSORART)
        
          incidence<<-c(incidence,100*out.incidence$incidence)
          
        }
        #incidence[current_i]=out.incidence$incidence
        
        
        
        

        
        n=length(incidence) - isolate(rv$prior.runs)
        rv$runs.completed=n
        
        if(n==0){
          return()
        }
        q=round(quantile(incidence, probs = c(0.025, 0.5, 0.975)),2)
        
        
        hist(incidence, breaks=20, col = 'darkgray', border = 'white', main=paste('Median:',q[2],'95% CI: (',q[1],'-',q[3],')'))
        if(isolate(rv$computing)&n<isolate(input$NRUNS)){
          invalidateLater(100, session)
        }
        else if (n>=isolate(input$NRUNS)){
          rv$computing=FALSE
          shinyjs::enable("GO")
          shinyjs::enable("RESET")
          shinyjs::disable("PAUSE")
        }
      })
   })
    output$progress <- renderPlot({
      par(mar=c(0,0,0,0))
      w.x=0.05
      w.y=0.05
      plot(NULL,NULL,xlim=c(-w.x,bar.aspect.ratio+w.x),ylim=c(-w.y,1+w.y), xlab='', ylab='',bty='n',xaxt='n',yaxt='n')
      r=rv$runs.completed
      
      N=isolate(input$NRUNS)
      rect(-w.x,-w.y,bar.aspect.ratio+w.x,1+w.y,border="Black")
      if(r>0){
        xvals=bar.aspect.ratio*seq(1,r)/N
 
        rect(xvals-bar.aspect.ratio/N,0,xvals,1, col = 'skyblue', border=NA)
      }
    }, height = 20, width=20*bar.aspect.ratio)
}



# Run the application 
shinyApp(ui = ui, server = server)

