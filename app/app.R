library(shiny)
source('model.R',local=TRUE)
colors <- c('#0072B2','#D55E00')
ui <- fluidPage(
 tags$head(tags$style(HTML('body{background:#f5f8fa;color:#183848;font-size:16px}.container-fluid{max-width:1400px}.well,.tab-content{background:white;border-radius:8px}.tab-content{padding:18px}.note{background:#e8f3f4;border-left:4px solid #007e78;padding:14px;margin-bottom:15px}h1,h3{color:#153b43}.btn-primary{background:#007e78}pre{white-space:pre-wrap}'))),
 titlePanel('Competition: who persists, and why?'),
 p('ENSCI 3120 | Predict, run, compare, explain'),
 sidebarLayout(sidebarPanel(width=3,
 selectInput('preset','Investigation',c('A: shared resources'='A','B: change competitive effects'='B','C: reverse competitive effects'='C','D: stronger effects, start 1'='D','E: stronger effects, start 2'='E')),
 actionButton('load','Load settings'),helpText('Load a scenario, make a prediction, then press Run.'),
 numericInput('N1','Initial species 1 abundance',20,min=0,max=2000),numericInput('N2','Initial species 2 abundance',20,min=0,max=2000),
 numericInput('r1','Species 1 growth rate (r1)',.2,min=.001,max=2,step=.05),numericInput('r2','Species 2 growth rate (r2)',.2,min=.001,max=2,step=.05),
 numericInput('K1','Species 1 carrying capacity (K1)',100,min=10,max=2000),numericInput('K2','Species 2 carrying capacity (K2)',100,min=10,max=2000),
 numericInput('alpha','alpha: effect of species 2 on species 1',.5,min=0,max=3,step=.1),
 numericInput('beta','beta: effect of species 1 on species 2',.5,min=0,max=3,step=.1),
 sliderInput('time','Duration (arbitrary time units)',min=1,max=500,value=150),
 actionButton('run','Run simulation',class='btn-primary'),tags$hr(),
 actionButton('keep','Keep run for comparison'),actionButton('clear','Clear comparison'),
 downloadButton('csv','Download current data')
 ),mainPanel(width=9,
 div(class='note',textOutput('status')),
 tabsetPanel(
 tabPanel('Investigate',h3('Predict before you run'),p('Will both species persist? Would changing their starting abundances change the outcome?'),
 textAreaInput('prediction','Your prediction (optional, saved only in this session)',rows=2,width='100%'),
 checkboxInput('alone','Show each species growing without its competitor',FALSE),
 plotOutput('trajectory',height='390px'),verbatimTextOutput('settings'),
 tags$details(tags$summary('Reveal the model interpretation'),verbatimTextOutput('interpretation'))),
 tabPanel('Phase plane',p('A filled circle marks the start; a square marks the end of the displayed run. Thin lines show zero-growth isoclines. Arrows show direction, not speed.'),
 checkboxInput('field','Show direction arrows',TRUE),plotOutput('phase',height='490px'),
 p('Species 1 grows below its isocline and declines above it (when N1 > 0). Species 2 grows below its isocline and declines above it (when N2 > 0). Each species also has zero growth on its own zero-abundance axis.')),
 tabPanel('Model and meaning',h3('Logistic growth with a competitor'),
 tags$pre('dN1/dt = r1*N1*[1 - (N1 + alpha*N2)/K1]\ndN2/dt = r2*N2*[1 - (N2 + beta*N1)/K2]'),
 p('alpha converts species 2 abundance into its competitive effect in species 1 units. beta does the reverse. For example, alpha = 0.5 means two species 2 individuals have the same modeled effect on species 1 as one species 1 individual.'),
 div(class='note','Notation reminder: beta here is a competition coefficient. It is NOT the transmission coefficient used in the SIR app.'),
 tags$pre('Species 1 isocline: N1 + alpha*N2 = K1\nSpecies 2 isocline: N2 + beta*N1 = K2'),
 p('Assumptions: a constant environment, positive intrinsic growth rates, nonnegative competitive effects, continuous abundances, no immigration, and no predators or parasites. Resources are represented indirectly through carrying capacities and competition coefficients.'),
 p('With both species present and both coefficients positive, competition is a negative effect on each species relative to the same species without a competitor. A winner can still be harmed by competition. At zero coefficients, one or both modeled effects are neutral.'),
 p('Very small populations are not automatically rounded to extinction. Seeing both species at the final plotted time does not prove stable coexistence. A species starting at zero remains absent.'),
 p('Original R implementation of the carrying-capacity formulation documented by ',a('EcoEvoApps',href='https://ecoevoapps.gitlab.io/docs/articles/lotka-volterra-competition.html',target='_blank'),'. Runtime: Shiny and base R.')),
 tabPanel('Classroom activity',h3('1. Can both species persist?'),p('Load A. Predict, run, then show the curves without competitors. Why is coexistence still a negative/negative interaction?'),
 h3('2. Does a winner benefit from competition?'),p('Keep A for comparison. Load B and predict. Run and compare each species with its no-competitor curve. Load C to reverse the asymmetry.'),
 h3('3. Does starting abundance matter?'),p('Load D, predict, and run. Keep D. Load E and repeat. Only the starting abundances change between D and E. Explain why the winner changes.'),
 h3('4. Distinguish stable coexistence from a balanced start'),p('Return to D and set both starting abundances to 20. Then change just species 1 to 21 and extend the duration. Why can an exactly balanced trajectory be misleading?'),
 p('Submit explanations in Top Hat when prompted. This app does not collect or transmit your answers.'))
 ))
))
server <- function(input,output,session) {
 config <- reactive({req(input$N1,input$N2,input$r1,input$r2,input$K1,input$K2,input$alpha,input$beta,input$time)
 list(N1=input$N1,N2=input$N2,r1=input$r1,r2=input$r2,K1=input$K1,K2=input$K2,alpha=input$alpha,beta=input$beta,time=input$time)})
 current <- reactiveVal(NULL);old <- reactiveVal(NULL);error <- reactiveVal(NULL)
 observeEvent(input$load,{
  p <- preset_config(input$preset)
  for(n in setdiff(names(p),'time')) updateNumericInput(session,n,value=p[[n]])
  updateSliderInput(session,'time',value=p$time)
 })
 observeEvent(input$run,{
  r <- tryCatch(simulate_model(config()),error=function(e)e)
  if(inherits(r,'error')){error(conditionMessage(r));return()}
  current(r);error(NULL)
 })
 observeEvent(input$keep,{req(current());old(current())})
 observeEvent(input$clear,old(NULL))
 output$status <- renderText({
  if(!is.null(error())) return(paste('Unable to run:',error(),'Any plots shown are from the last successful run.'))
  if(is.null(current())) return('Make a prediction, then press Run simulation.')
  if(!identical(config(),current()$config)) return('Controls have changed. Plots still show the last completed run; press Run to update.')
  'Plots show your completed run. Keep it to compare with a second run.'
 })
 result <- reactive({validate(need(!is.null(current()),'Press Run simulation to display results.'));current()})
 label <- function(p) sprintf('N1=%g, N2=%g; r1=%g, r2=%g; K1=%g, K2=%g; alpha=%g, beta=%g; duration=%g',p$N1,p$N2,p$r1,p$r2,p$K1,p$K2,p$alpha,p$beta,p$time)
 output$settings <- renderText({r<-result();paste('Current:',label(r$config),if(!is.null(old()))paste('\nKept:',label(old()$config))else '')})
 output$interpretation <- renderText({r<-result();p<-r$config;d<-r$data
 paste(classify(p),sprintf('\nAt plotted end: N1 = %.4f; N2 = %.4f.',tail(d$N1,1),tail(d$N2,1)),
 sprintf('\nGrowth rates when rare against a resident at K:\nSpecies 1: %.4f; species 2: %.4f per time unit.',p$r1*(1-p$alpha*p$K2/p$K1),p$r2*(1-p$beta*p$K1/p$K2)),
 '\nClassification describes the model with both species initially present; it is not inferred just from the last plotted values. Exact balanced starts can remain on an unstable trajectory.')})
 output$trajectory <- renderPlot({
  r<-result();d<-r$data;o<-old();yy<-c(d$N1,d$N2)
  if(isTRUE(input$alone)) yy<-c(yy,d$alone1,d$alone2)
  if(!is.null(o)) yy<-c(yy,o$data$N1,o$data$N2)
  xmax<-max(d$time,if(!is.null(o))o$data$time else 0)
  par(mar=c(4,4,2,1));matplot(d$time,d[c('N1','N2')],type='l',lty=1,lwd=3,col=colors,xlim=c(0,xmax),ylim=c(0,max(1,yy)),xlab='Time (arbitrary units)',ylab='Population abundance')
  if(!is.null(o))matlines(o$data$time,o$data[c('N1','N2')],lty=2,lwd=2,col=colors)
  if(isTRUE(input$alone))matlines(d$time,d[c('alone1','alone2')],lty=3,lwd=2,col=colors)
  legend('topright',c('Species 1','Species 2','Current: solid','Kept: dashed','Without competitor: dotted'),col=c(colors,rep('grey40',3)),lty=c(1,1,1,2,3),lwd=2,cex=.85,bg='white')
 },res=110)
 output$phase <- renderPlot({
  r<-result();p<-r$config;d<-r$data;o<-old()
  xmax<-1.1*max(p$K1,d$N1,if(!is.null(o))o$data$N1 else 0)
  ymax<-1.1*max(p$K2,d$N2,if(!is.null(o))o$data$N2 else 0)
  par(mar=c(4,4,2,1));plot(0,0,type='n',xlim=c(0,xmax),ylim=c(0,ymax),xlab='Species 1 abundance (N1)',ylab='Species 2 abundance (N2)')
  if(isTRUE(input$field)){
   g<-expand.grid(x=seq(.05*xmax,.95*xmax,length.out=13),y=seq(.05*ymax,.95*ymax,length.out=13))
   dx<-p$r1*g$x*(1-(g$x+p$alpha*g$y)/p$K1)/xmax;dy<-p$r2*g$y*(1-(g$y+p$beta*g$x)/p$K2)/ymax
   norm<-sqrt(dx^2+dy^2);ok<-norm>1e-12
   arrows(g$x[ok],g$y[ok],g$x[ok]+.035*xmax*dx[ok]/norm[ok],g$y[ok]+.035*ymax*dy[ok]/norm[ok],length=.055,col='#C3CCCC')
  }
  if(p$alpha>0) abline(a=p$K1/p$alpha,b=-1/p$alpha,col=colors[1],lwd=2) else abline(v=p$K1,col=colors[1],lwd=2)
  abline(a=p$K2,b=-p$beta,col=colors[2],lwd=2)
  if(!is.null(o))lines(o$data$N1,o$data$N2,lty=2,lwd=2,col='#777777')
  lines(d$N1,d$N2,lwd=3,col='#183848');points(d$N1[1],d$N2[1],pch=16,cex=1.3);points(tail(d$N1,1),tail(d$N2,1),pch=15)
  legend('topright',c('Species 1 isocline','Species 2 isocline','Current trajectory','Kept trajectory'),col=c(colors,'#183848','#777777'),lty=c(1,1,1,2),lwd=2,cex=.8,bg='white')
 },res=110)
 output$csv <- downloadHandler(filename=function()'competition-run.csv',content=function(file){r<-result();d<-r$data;for(n in names(r$config))d[[paste0('setting_',n)]]<-r$config[[n]];write.csv(d,file,row.names=FALSE)})
}
shinyApp(ui,server)
