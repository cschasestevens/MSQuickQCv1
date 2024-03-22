#essential columns (normalized data) are: sampleType, time, and label
#essential columns (annotations and MSI) : Class, MSMS match, mzrt match, 
# Similar MSMS match, Similar mzrt match

options(shiny.maxRequestSize = 50 * 1024^2)
#### Set-up ####

## Load Shiny-specific libraries
library(shinyBS)
library(DT)
library(shinyFeedback)
library(shiny)
library(IMIFA)
library(ComplexHeatmap)
library(readxl)
library(ggplot2)
library(purrr)
library(ggpubr)
library(ggrepel)
library(RColorBrewer)
library(circlize)
library(corrplot)
library(magrittr)
library(viridis)
library(patchwork)
library(kableExtra)
library(ggsci)

## Source scripts

source("scripts/1_lib.R",
       local = knitr::knit_global())

## Parameter list (if present in working directory)
if(
  file.exists(
    "qc_params.txt"
    )
  ) {
  
  list.params <- read.table("qc_params.txt",header = T,sep = "\t")
  
  }

exists("list.params")




#### UI ####

## Title

ui.title <- h2(id = "title1", "MSQuickQC: A rapid solution for quality control of mass spectrometry data sets")

## Add section break

ui.section.break <- hr(id = "brsect")


#### UI - Instructions ####

ui.instructions.title <- column(12,
                                h4("Formatting Instructions"),
                                id = "tpar")

ui.instructions <- column(12,

  p("- Accepts .xlsx, .csv, and tab-delimited .txt files and can accommodate up to two data sets for each report;",
    br("  file types must match if two files are uploaded!")),
  
  p("- All data files must include the following columns (case-sensitive): ",
     strong("sampleType",style = "color: #E74C3C"),", ",
     strong("label",style = "color: #E74C3C"),", ",
     strong("Group",style = "color: #E74C3C")),
  
  p("- Optional columns include (case-sensitive): ",
     strong("time",style = "color: #E74C3C")),
  
  p("- Descriptions of each column:",
  br(" ",strong("sampleType",style = "color: #E74C3C"),": ","distinguishes experimental samples from QCs"),
  br(", ",strong("time",style = "color: #E74C3C"),": ","indicates acquisition order of each sample"),
  br(", ",strong("label",style = "color: #E74C3C"),": ","unique sample IDs"),
  br(", ",strong("Group",style = "color: #E74C3C"),": ","treatment group of each sample")),
  
  p("- Metadata columns should appear before compound intensities but do not need to be in specific order"),
  
  p("- Experimental samples in ",
     strong("sampleType",style = "color: #E74C3C"),
     " must be labeled as",
     strong("Sample",style = "color: #E74C3C")),
  
  p("- If quality control samples (i.e. pooled samples, standard reference samples, method blanks) are not present,",
    " denote by 'NA' for each QC type when specifying data parameters. Otherwise, include QC names in 'sampleType'",
    " and 'Group.'"),
  
  p("- Zeroes in data sets are acceptable and are replaced with 10% of the minimum value present in the data",
    "  to prevent generation of 'infinite/NA' values during data transformation and PCA."),
  
  p("*See example data for specific formatting help:",
  br(" (",tags$a("Example file",href = "example_file.xlsx"), ") with all QCs, internal standards, and multiple treatments"),
  br(", (",tags$a("Additional example",href = "example_file2.xlsx"), ") with all QCs, internal standards, and multiple treatments"),
  br(", (",tags$a("Example file",href = "example_file_noqc.xlsx"), ") with no QCs, internal standards, or acquisition order")),
  
  p("*You can also run the report script outside of the app using this ",tags$a(".RMD file",href = "MSQuickQC_v4_input.Rmd")),
  p("See the", tags$a("README",href = "README.md") ,"for more information and email questions to Chase Stevens at either: ",
    strong("Nathanial_Stevens@med.unc.edu",style = "color: #E74C3C"),
    "or",
    strong("cschasestevens@gmail.com",style = "color: #E74C3C")
    )
  )


#### UI - Data Upload ####

ui.upload.title <- column(12,
                          h4("Input Parameters"),
                          id = "tupl")



ifelse(exists("list.params"),
       ui.upload <- column(12,
                           
                           # Input files
                           column(6,
                                  h5("Input File(s):"),
                                  fileInput('file1',
                                            "Select data files...",
                                            multiple = T,
                                            accept = c(".xlsx",".csv",".txt")),
                                  
                                  # List parameter section
                                  h5("List Parameters:",
                                     br("- type 'NA' if not used"),
                                     br("and separate entries for 2-file inputs by a ",
                                        strong(",",style = "color: #E74C3C")
                                     )
                                  ),
                                  
                                  textInput("qc1","QC #1",
                                            value = ifelse(is.na(list.params$qc1),
                                                           "NA",list.params$qc1)),
                                  textInput("qc2","QC #2",
                                            value = ifelse(is.na(list.params$qc2),
                                                           "NA",list.params$qc2)),
                                  textInput("qc3","QC #3",
                                            value = ifelse(is.na(list.params$qc3),
                                                           "NA",list.params$qc3)),
                                  textInput("istd","Internal Standard Identifier (provide one)",
                                            value = ifelse(is.na(list.params$iSTD),
                                                           "NA",list.params$iSTD)),
                                  textInput("md","Metadata Columns",
                                            value = list.params$Metadata.col),
                                  textInput("date","Acquisition Date(s) (yyyy/mm/dd)",
                                            value = list.params$Date),
                                  textInput("study","Study Name (provide one)",
                                            value = list.params$Study.Name),
                                  textInput("plat","Analysis Platform(s)",
                                            value = list.params$Platform)),
                           column(6,
                                  textInput("polar","Polarity",
                                            value = list.params$Polarity),
                                  textInput("instr","Instrument(s)",
                                            value = list.params$Instrument),
                                  textInput("proc","Processing Software",
                                            value = list.params$Processing),
                                  textInput("norm","Normalization Method",
                                            value = list.params$Norm.meth),
                                  textInput("auth","Report Author",
                                            value = list.params$Report.Author),
                                  
                                  ## After input params are specified
                                  actionButton("submit","Generate Parameter List"),
                                  
                                  # Report generation/downloads
                                  h5('Quality Control Screening:'),
                                  downloadButton("report", 
                                                 "Generate Report"),
                                  downloadButton("plots", 
                                                 "Download Plots"),
                                  downloadButton("tables", 
                                                 "Download Tables"))),
       ui.upload <- column(12,
                           
                           # Input files
                           column(6,
                                  h5("Input File(s):"),
                                  fileInput('file1',
                                            "Select data files...",
                                            multiple = T,
                                            accept = c(".xlsx",".csv",".txt")),
                                  
                                  # List parameter section
                                  h5("List Parameters:",
                                     br("- type 'NA' if not used"),
                                     br("and separate entries for 2-file inputs by a ",
                                        strong(",",style = "color: #E74C3C")
                                     )
                                  ),
                                  
                                  textInput("qc1","QC #1",
                                            value = "ex. qc1 or qc1a,qc1b"),
                                  textInput("qc2","QC #2",
                                            value = "ex. qc2 or qc2a,qc2b"),
                                  textInput("qc3","QC #3",
                                            value = "ex. qc3 or qc3a,qc3b"),
                                  textInput("istd","Internal Standard Identifier (provide one)",
                                            value = "iSTD"),
                                  textInput("md","Metadata Columns",
                                            value = "4,4"),
                                  textInput("date","Acquisition Date(s) (yyyy/mm/dd)",
                                            value = "ex. 2024/01/01 or 2024/01/01,2024/01/01"),
                                  textInput("study","Study Name (provide one)",
                                            value = "Mouse lung lipidomics"),
                                  textInput("plat","Analysis Platform(s)",
                                            value = "ex. CSH or CSH,HILIC")),
                           column(6,
                                  textInput("polar","Polarity",
                                            value = "ex. POS or POS,NEG"),
                                  textInput("instr","Instrument(s)",
                                            value = "ex. Thermo Q-Exactive HF or QEHF,QTOF"),
                                  textInput("proc","Processing Software",
                                            value = "ex. MS-DIAL v.4.0 or MS-DIAL v.4.0,MS-DIAL v.4.1"),
                                  textInput("norm","Normalization Method",
                                            value = "ex. SERRF or Raw,SERRF"),
                                  textInput("auth","Report Author (provide one)","Insert Name"),
                                  
                                  ## After input params are specified
                                  actionButton("submit","Generate Parameter List"),
                                  
                                  # Report generation/downloads
                                  h5('Quality Control Screening:'),
                                  downloadButton("report", 
                                                 "Generate Report"),
                                  downloadButton("plots", 
                                                 "Download Plots"),
                                  downloadButton("tables", 
                                                 "Download Tables"))))


#### UI - Data Preview ####

ui.prev.title <- column(12,
                          h4("Data Preview"),
                          id = "prev")

ui.file.prev <- column(4,
                       dataTableOutput('contents1', 
                                       width = "80%"))

ui.data.prev <- column(4,
                       tabsetPanel(tabPanel("File 1",
                                            dataTableOutput('contents2.f1',
                                                            width = "80%")),
                                   tabPanel("File 2",
                                            dataTableOutput('contents2.f2',
                                                            width = "80%")
                                            )
                                   )
                       )


ui.parm.prev <- column(4,
                       dataTableOutput('contents3',
                                       width = "80%"))





#### Combined UI elements ####

ui <- fluidPage(includeCSS("www/style_app.css"),
        shinyFeedback::useShinyFeedback(),
        ui.title,
        ui.section.break,
        fluidRow(ui.instructions.title,
                 ui.instructions,
                 ui.upload.title,
                 ui.upload),
        ui.section.break,
        fluidRow(ui.prev.title,
                 ui.file.prev,
                 ui.data.prev,
                 ui.parm.prev))
                  
               

      



















#### Server ####

# Define server logic required to draw a histogram
server = function(input, output){
    
  ## Set default values for input file and parameters
  
  values <- reactiveValues(file.in = NULL,
                           list.qc.data = NULL)
  
  
  ## Save input files to app files folder
  
  observeEvent(input$file1,
               {
                 
                 ### If 1 input file is detected
                 
                 if(length(input$file1[["name"]]) == 1) {
                   
                   if(tools::file_ext(input$file1[["name"]]) == "xlsx") {
                     
                     values$file.in <- readxl::read_excel(input$file1$datapath)
                     
                   }
                   
                   if(tools::file_ext(input$file1[["name"]]) == "csv") {
                     
                     values$file.in <- read.csv(input$file1$datapath,
                              header = T,
                              sep = ",")
                     
                   }
                   
                   if(tools::file_ext(input$file1[["name"]]) == "txt") {
                     
                     values$file.in <- read.table(input$file1$datapath,
                                header = T,
                                sep = "\t",
                                quote = "")
                     
                   }
                   
                   
                   ## Table outputs
                   
                   output$contents1 <- renderDataTable({

                     datatable(input$file1,
                               extensions = c('Buttons','Scroller'),
                               options = list(scrollX = "80%",
                                              scrollY = "200px"))
                     
                     
                   })
                   
                   output$contents2.f1 <- renderDataTable({

                     datatable(head(values$file.in[1:10],5),
                               extensions = c('Buttons','Scroller'),
                               options = list(scrollX = "80%",
                                              scrollY = "200px"))
                     
                     
                   })
                   
                 }
                 
                 
                 ### If 2 input files are detected
                 
                 if(length(input$file1[["name"]]) == 2) {
                   
                   ## Read input files
                   
                   ifelse(tools::file_ext(input$file1[1,"name"]) == "xlsx",
                          values$file.in <- lapply(1:2,
                                                   function(x) readxl::read_excel(input$file1[x,"datapath"])),
                          ifelse(tools::file_ext(input$file1[1,"name"]) == "csv",
                                        values$file.in <- lapply(1:2,
                                                                 function(x) read.csv(input$file1[x,"datapath"],
                                                                   header = T,
                                                                   sep = ",",
                                                                   check.names = T)),
                                        values$file.in <- lapply(1:2,
                                                                 function(x) read.table(input$file1[x,"datapath"],
                                                                     header = T,
                                                                     sep = "\t",
                                                                     quote = "",
                                                                     check.names = T))
                                        )
                                 )
                   
                   
                   ## Table outputs
                   
                   ### Input file path(s)
                   
                   output$contents1 <- renderDataTable({
                     
                     datatable(input$file1,
                               extensions = c('Buttons','Scroller'),
                               options = list(scrollX = "80%",
                                              scrollY = "200px"))
                     
                     
                   })
                   
                   ### File 1 data preview
                   
                   output$contents2.f1 <- renderDataTable({
                     
                     datatable(
                       cbind(
                         head(
                           values$file.in[[1]][1:10],5)
                         ),
                       extensions = c('Buttons','Scroller'),
                       options = list(scrollX = "80%",
                                      scrollY = "200px")
                       )
                     
                     
                   })
                   
                   ### File 2 data preview
                   
                   output$contents2.f2 <- renderDataTable({
                     
                     datatable(
                       cbind(
                         head(
                           values$file.in[[2]][1:10],5)
                         ),
                       extensions = c('Buttons','Scroller'),
                       options = list(scrollX = "80%",
                                      scrollY = "200px")
                       )
                     
                     
                   })
                   
                   
                   
                 }

                 
                 
                 ## Save data in folder for QC screening
                 
                 ### 1 input file
                 
                 if(length(input$file1[["name"]]) == 1) {
                   
                   write.table(values$file.in,
                               file = paste(gsub("\\..*","",
                                                 input$file1[["name"]]
                                            ),
                                            ".txt",
                                            sep = ""),
                               sep = "\t",
                               row.names = F
                   )
                   
                 }
                 
                 ### 2 input files
                 
                 if(length(input$file1[["name"]]) == 2) {
                   
                   lapply(as.vector(seq(1:length(input$file1[["name"]])
                                        )
                                    ),
                   function(x) {
                     
                     write.table(values$file.in[[x]],
                                 file = paste(gsub("\\..*","",
                                                   input$file1[[x,"name"]]
                                              ),
                                              ".txt",
                                              sep = ""),
                                 sep = "\t",
                                 row.names = F
                     )
                     
                   })
                   
                 }
                 
               })


  observeEvent(input$submit,
               {
                 req(input$file1)

                 
                 if(length(input$file1[["name"]]) == 2){
                   
                   values$list.qc.data <- data.frame(File.ID = input$file1[["name"]],
                                                     Path = input$file1[["datapath"]],
                                                     qc1 = unlist(strsplit(input$qc1,",")),
                                                     qc2 = unlist(strsplit(input$qc2,",")),
                                                     qc3 = unlist(strsplit(input$qc3,",")),
                                                     iSTD = rep(input$istd,2),
                                                     Metadata.col = as.numeric(unlist(strsplit(input$md,","))),
                                                     Date = unlist(strsplit(input$date,",")),
                                                     Study.Name = rep(input$study,2),
                                                     Platform = unlist(strsplit(input$plat,",")),
                                                     Polarity = unlist(strsplit(input$polar,",")),
                                                     Instrument = unlist(strsplit(input$instr,",")),
                                                     Processing = unlist(strsplit(input$proc,",")),
                                                     Norm.meth = unlist(strsplit(input$norm,",")),
                                                     Report.Author = rep(input$auth,2))
                   
                 }
                 
                 if(length(input$file1[["name"]]) == 1){
                   
                   values$list.qc.data <- data.frame(File.ID = input$file1[["name"]],
                                                     Path = input$file1[["datapath"]],
                                                     qc1 = input$qc1,
                                                     qc2 = input$qc2,
                                                     qc3 = input$qc3,
                                                     iSTD = input$istd,
                                                     Metadata.col = as.numeric(input$md),
                                                     Date= input$date,
                                                     Study.Name = input$study,
                                                     Platform = input$plat,
                                                     Polarity = input$polar,
                                                     Instrument = input$instr,
                                                     Processing = input$proc,
                                                     Norm.meth = input$norm,
                                                     Report.Author = input$auth)
                   
                 }
                   
                   
                   
                   output$contents3 <- renderDataTable({
                     
                     datatable(
                       cbind(
                         values$list.qc.data
                       ),
                       extensions = c('Buttons','Scroller'),
                       options = list(scrollX = "80%",
                                      scrollY = "200px")
                     )
                     })
                   
                   
                   write.table(values$list.qc.data,
                               file = paste("qc_params",
                                            ".txt",
                                            sep = ""),
                               sep = "\t",
                               row.names = F
                   )
                   
                   
                 })
  
  
  
  output$report <- downloadHandler(
    filename = "report_output.html",
    content = function(file) {
      withProgress(message = 'Please wait while the
                         report is generated...', {
                           
                           # Required data inputs
                           
                           req(input$file1,input$submit,
                               values$file.in,values$list.qc.data)
                           
                           # Parameters to pass into R markdown
                           
                           ifelse(length(input$file1[["name"]]) == 2,
                                  list.data.in <-  list("File 1" = list("Datasheet" = read.table(paste(gsub("\\..*","",
                                                                                                 input$file1[[1,"name"]]
                                                                                            ),
                                                                                            ".txt",
                                                                                            sep = ""),
                                                                                      header = T)),
                                                        "File 2" = list("Datasheet" = read.table(paste(gsub("\\..*","",
                                                                                                 input$file1[[2,"name"]]
                                                                                            ),
                                                                                            ".txt",
                                                                                            sep = ""),
                                                                                      header = T))),
                                  list.data.in <-  list("File 1" = list("Datasheet" = read.table(paste(gsub("\\..*","",
                                                                                                                                  input$file1[["name"]]
                                                                                                                             ),
                                                                                                                             ".txt",
                                                                                                                             sep = ""),
                                                                                                                       header = T))))
                           
                           
                           params <- list(list.data.in = list.data.in,
                                          list.qc.data.in = read.table(paste("qc_params",
                                                                             ".txt",
                                                                             sep = ""),
                                                                       header = T),
                                          rendered_by_shiny = T)
                           
                           # Generate report
                           
                           rmarkdown::render("MSQuickQC_v1.1_input_shiny.Rmd", 
                                             output_file = file,
                                             params = params,
                                             envir = new.env(parent = globalenv()
                                                             )
                                             )
                           
                         })
      
      })
  
  output$plots <- downloadHandler(
    filename = "plots.zip",
    
    content = function(file) {
      
      file.copy("plots.zip",
                file,
                overwrite = T)
      
      unlink("plots/",recursive = T)
      
      unlink("plots.zip",recursive = T)
      
      
    }
    
  )
  
  output$tables <- downloadHandler(
    filename = "tables.zip",
    
    content = function(file) {
      
      file.copy("tables.zip",
                file,
                overwrite = T)
      
      unlink("tables/",recursive = T)
      
      unlink("tables.zip",recursive = T)
      
      
    }
    
  )
  
  
  
  
  
  
  
    
    
    
}

# Run the application 
shinyApp(ui = ui, server = server)
