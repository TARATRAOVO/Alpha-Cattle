homePanel<-tabPanel(
  title = "Home",
  icon = icon("home"),
  br(),
  column(
    width = 8,offset = 2,
  #div(class = "inlay", style = "height:200px;width:100%;background-color: #2c3e50;"),
  p("Characteristics",
    style='font-size:50px;margin-top:80px;margin-bottom:20px;'
  )),
  column(
    width = 8,offset = 2,
        div(class = "jumbotron",
            h1("HCGSP"), p("Holstein Cattle Genomic Selection and Prediction(HCGSP) provides
              a platform and software to enable multiple modes of genomic selection and Prediction in Holstein Cattles."), p(a(class = "btn btn-primary btn-lg action-button shiny-bound-input", id='tabBut', "Learn more")))


   ),
  #imageOutput("pic"),

 tags$style(HTML(
    " #pic_pca{
    text-align:center;
    }
    #pic_cor{
    text-align:center;
    }
    #pic_joint{
    text-align:center;
    }
    "
  )),
 div(class = "inlay", style = "height:650px;width:100%;background-color: white;"),

 
column(
   width = 8,offset = 2,
   tags$div(
     class = "document",

     div(tags$b("Genomic selection for combined population", style='font-size:25px;')),
     br(),

      tags$li("The accuracy of GEBV is largely influenced by the size of the training population (Goddard and Hayes 2009). 
                     However, in animal breeding practice, it is sometimes difficult to obtain a large training population. 
                     For instance, due to capital, management and other reasons, the number of Holstein cattles with genotype information measured in a ranch is often small.
                     The cost of genotyping is still quite high for large-scale individuals. 
                     In response, one solution is to carry out joint breeding to enlarge the training population and further improve the accuracy of GEBV.",
                     ),
      br(),
             imageOutput("pic_pca"),
     br(),
     tags$li("We test a combined reference population was established with animals from two different populations (TB and XX) of the Yorkshire pig. Differences in 
             genetic background (Figure A) were evident between the two populations and the accuracies of predictions using the merged population and the single population were compared (Figure B). 
             Our results confirmed that using an admixed reference population is meaningful for joint breeding.",
     ),
     br(),
     div(tags$b("Large Reference Population and Corresponding Traits", style='font-size:25px;')),
     br(),
     tags$li("We have 3847 reference populations with genome-wide information. They come from Guangming Ranch across China, and they share the same characteristics as the Chinese Holstein cattle genome. 
             At the same time, we also have matching phenotype group information, including milk yield, milk protein, DHI, somatic cell score, body type score, limb and hoof character score, birth spacing, service life 
             and many other traits you care about. With the help of such a large phenotype group and genome data, after you upload your phenotype group and genome files, you can get very accurate trait prediction results 
             to help you decide the leaving or staying issue of each bull.
             "),
      br(),
     imageOutput("pic_cor"),
     br(),
     div(tags$b("Predict CPI and select according to CPI", style='font-size:25px;')),
     br(),
     tags$li(" Our database contains all the required traits for calculating the CPI of Holstein cattle. CPI index is applicable to domestic progeny testing and verification bulls with daughter production traits 
              and daughter body type identification results that are routinely evaluated in China. At the same time, we can also calculate the GCPI index. If you want to select bulls according to an index, 
              we can calculate the CPI and GCPI of the individual bulls you uploaded (genome sequencing information is required)
             "),
     br(),
     imageOutput("pic_joint"),
     br(),
   )
 ),

column(width = 12,align="center",
       style = "height:50px;color:white;width:105%;background-color: #2c3e50;margin-left:-15px;margin-top:180px", 
       h6("Copyright 2022 College of Agriculture and Biology, ShangHaiJiaoTong University. 800# Dongchuan Road, Shanghai, China. E-mail: panyuchun1963@aliyun.com. all rights reserved."),
       
        )


)