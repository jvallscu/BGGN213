Class18: Cancer genomics and Immunoinformatics
================
Joan M. Valls Cuevas
5/31/2019

What is cancer?
---------------

Q5. What is the most frequentely mutated position associated with cancer in the KRas protein (i.e. the amino acid with the most mutations)?

Pancreas


``` r
BiocManager::install("GenomicDataCommons")
```
Load tools
``` r
library(GenomicDataCommons)

library(TCGAbiolinks)

library(maftools)
library(bio3d)
```

``` r
status()
```

    ## $commit
    ## [1] "3e22a4257d5079ae9f7193950b51ed9dfc561ed1"
    ## 
    ## $data_release
    ## [1] "Data Release 17.0 - June 05, 2019"
    ## 
    ## $status
    ## [1] "OK"
    ## 
    ## $tag
    ## [1] "1.21.0"
    ## 
    ## $version
    ## [1] 1

``` r
projects <- getGDCprojects()
head(projects)
```

    ##   dbgap_accession_number
    ## 1              phs001287
    ## 2              phs001374
    ## 3              phs001628
    ## 4              phs000466
    ## 5              phs000467
    ## 6              phs001179
    ##                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 disease_type
    ## 1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               Adenomas and Adenocarcinomas
    ## 2                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         Epithelial Neoplasms, NOS, Squamous Cell Neoplasms
    ## 3                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          Myeloid Leukemias
    ## 4                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           Clear Cell Sarcoma of the Kidney
    ## 5                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              Neuroblastoma
    ## 6 Germ Cell Neoplasms, Acinar Cell Neoplasms, Miscellaneous Tumors, Thymic Epithelial Neoplasms, Gliomas, Basal Cell Neoplasms, Neuroepitheliomatous Neoplasms, Ductal and Lobular Neoplasms, Complex Mixed and Stromal Neoplasms, Complex Epithelial Neoplasms, Adnexal and Skin Appendage Neoplasms, Mesothelial Neoplasms, Mucoepidermoid Neoplasms, Not Reported, Specialized Gonadal Neoplasms, Cystic, Mucinous and Serous Neoplasms, Adenomas and Adenocarcinomas, Epithelial Neoplasms, NOS, Squamous Cell Neoplasms, Transitional Cell Papillomas and Carcinomas, Paragangliomas and Glomus Tumors, Nevi and Melanomas, Meningiomas
    ##   releasable released state
    ## 1      FALSE     TRUE  open
    ## 2      FALSE     TRUE  open
    ## 3      FALSE     TRUE  open
    ## 4       TRUE     TRUE  open
    ## 5       TRUE     TRUE  open
    ## 6      FALSE     TRUE  open
    ##                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               primary_site
    ## 1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   Kidney, Bronchus and lung, Uterus, NOS
    ## 2                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        Bronchus and lung
    ## 3                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            Hematopoietic and reticuloendothelial systems
    ## 4                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   Kidney
    ## 5                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           Nervous System
    ## 6 Testis, Gallbladder, Unknown, Other and unspecified parts of biliary tract, Adrenal gland, Thyroid gland, Spinal cord, cranial nerves, and other parts of central nervous system, Peripheral nerves and autonomic nervous system, Stomach, Cervix uteri, Bladder, Small intestine, Breast, Prostate gland, Other and unspecified female genital organs, Other and unspecified major salivary glands, Rectum, Retroperitoneum and peritoneum, Pancreas, Heart, mediastinum, and pleura, Bronchus and lung, Liver and intrahepatic bile ducts, Other and ill-defined sites, Thymus, Penis, Nasopharynx, Ovary, Uterus, NOS, Vulva, Anus and anal canal, Other and unspecified urinary organs, Trachea, Ureter, Other endocrine glands and related structures, Not Reported, Colon, Kidney, Vagina, Skin, Esophagus, Eye and adnexa, Other and ill-defined digestive organs
    ##              project_id                    id
    ## 1               CPTAC-3               CPTAC-3
    ## 2        VAREPOP-APOLLO        VAREPOP-APOLLO
    ## 3 BEATAML1.0-CRENOLANIB BEATAML1.0-CRENOLANIB
    ## 4           TARGET-CCSK           TARGET-CCSK
    ## 5            TARGET-NBL            TARGET-NBL
    ## 6                 FM-AD                 FM-AD
    ##                                                                                              name
    ## 1                                                                                                
    ## 2                                                          VA Research Precision Oncology Program
    ## 3 Clinical Resistance to Crenolanib in Acute Myeloid Leukemia Due to Diverse Molecular Mechanisms
    ## 4                                                                Clear Cell Sarcoma of the Kidney
    ## 5                                                                                   Neuroblastoma
    ## 6                                       Foundation Medicine Adult Cancer Clinical Dataset (FM-AD)
    ##        tumor
    ## 1          3
    ## 2     APOLLO
    ## 3 CRENOLANIB
    ## 4       CCSK
    ## 5        NBL
    ## 6         AD

``` r
cases_by_project <- cases() %>%
  facet("project.project_id") %>%
  aggregations()
head(cases_by_project)
```

    ## $project.project_id
    ##                      key doc_count
    ## 1                  FM-AD     18004
    ## 2             TARGET-NBL      1120
    ## 3              TCGA-BRCA      1098
    ## 4             TARGET-AML       988
    ## 5              TARGET-WT       652
    ## 6               TCGA-GBM       617
    ## 7                TCGA-OV       608
    ## 8              TCGA-LUAD       585
    ## 9              TCGA-UCEC       560
    ## 10             TCGA-KIRC       537
    ## 11             TCGA-HNSC       528
    ## 12              TCGA-LGG       516
    ## 13             TCGA-THCA       507
    ## 14             TCGA-LUSC       504
    ## 15             TCGA-PRAD       500
    ## 16          NCICCR-DLBCL       489
    ## 17             TCGA-SKCM       470
    ## 18             TCGA-COAD       461
    ## 19             TCGA-STAD       443
    ## 20             TCGA-BLCA       412
    ## 21             TARGET-OS       381
    ## 22             TCGA-LIHC       377
    ## 23               CPTAC-3       322
    ## 24             TCGA-CESC       307
    ## 25             TCGA-KIRP       291
    ## 26             TCGA-SARC       261
    ## 27             TCGA-LAML       200
    ## 28             TCGA-ESCA       185
    ## 29             TCGA-PAAD       185
    ## 30             TCGA-PCPG       179
    ## 31             TCGA-READ       172
    ## 32             TCGA-TGCT       150
    ## 33         TARGET-ALL-P3       131
    ## 34             TCGA-THYM       124
    ## 35             TCGA-KICH       113
    ## 36              TCGA-ACC        92
    ## 37             TCGA-MESO        87
    ## 38              TCGA-UVM        80
    ## 39             TARGET-RT        75
    ## 40             TCGA-DLBC        58
    ## 41              TCGA-UCS        57
    ## 42 BEATAML1.0-CRENOLANIB        56
    ## 43             TCGA-CHOL        51
    ## 44           CTSP-DLBCL1        45
    ## 45           TARGET-CCSK        13
    ## 46             HCMI-CMDC         7
    ## 47        VAREPOP-APOLLO         7

``` r
x <- cases_by_project$project.project_id

# Make a custom color vector for our plot
colvec <- rep("lightblue", nrow(x))
colvec[x$key == "TCGA-PAAD"] <- "red"

# Plot with 'log' for y axis and rotate labels with 'las'
#par(___)  
barplot(x$doc_count, names.arg= x$key, log="y", col=colvec, las=2)
```

![](class18_files/figure-markdown_github/unnamed-chunk-6-1.png)

Section 2
---------

> P53\_wt Cellular tumor antigen p53 - Healthy Tissue MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP DEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAK SVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHE RCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNS SCMGGMNRRPILTIITLEDSSGNLLGRNSFEVRVCACPGRDRRTEEENLRKKGEPHHELP PGSTKRALPNNTSSSPQPKKKPLDGEYFTLQIRGRERFEMFRELNEALELKDAQAGKEPG GSRAHSSHLKSKKGQSTSRHKKLMFKTEGPDSD

> P53\_mutant Cellular tumor antigen p53 - Tumor Tissue MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMLDLMLSPDDIEQWFTEDPGP DEAPWMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAK SVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHE RCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFVHSVVVPYEPPEVGSDCTTIHYNYMCNS SCMGGMNRRPILTIITLEV

``` r
my_data <- read.fasta("lecture18_sequences.fa")

read.fasta("lecture18_sequences.fa")
```

    ##              1        .         .         .         .         .         60 
    ## P53_wt       MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGP
    ## P53_mutant   MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMLDLMLSPDDIEQWFTEDPGP
    ##              **************************************** ******************* 
    ##              1        .         .         .         .         .         60 
    ## 
    ##             61        .         .         .         .         .         120 
    ## P53_wt       DEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAK
    ## P53_mutant   DEAPWMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAK
    ##              **** ******************************************************* 
    ##             61        .         .         .         .         .         120 
    ## 
    ##            121        .         .         .         .         .         180 
    ## P53_wt       SVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHE
    ## P53_mutant   SVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHE
    ##              ************************************************************ 
    ##            121        .         .         .         .         .         180 
    ## 
    ##            181        .         .         .         .         .         240 
    ## P53_wt       RCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNS
    ## P53_mutant   RCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFVHSVVVPYEPPEVGSDCTTIHYNYMCNS
    ##              ******************************** *************************** 
    ##            181        .         .         .         .         .         240 
    ## 
    ##            241        .         .         .         .         .         300 
    ## P53_wt       SCMGGMNRRPILTIITLEDSSGNLLGRNSFEVRVCACPGRDRRTEEENLRKKGEPHHELP
    ## P53_mutant   SCMGGMNRRPILTIITLEV-----------------------------------------
    ##              ******************                                           
    ##            241        .         .         .         .         .         300 
    ## 
    ##            301        .         .         .         .         .         360 
    ## P53_wt       PGSTKRALPNNTSSSPQPKKKPLDGEYFTLQIRGRERFEMFRELNEALELKDAQAGKEPG
    ## P53_mutant   ------------------------------------------------------------
    ##                                                                           
    ##            301        .         .         .         .         .         360 
    ## 
    ##            361        .         .         .  393 
    ## P53_wt       GSRAHSSHLKSKKGQSTSRHKKLMFKTEGPDSD
    ## P53_mutant   ---------------------------------
    ##                                                
    ##            361        .         .         .  393 
    ## 
    ## Call:
    ##   read.fasta(file = "lecture18_sequences.fa")
    ## 
    ## Class:
    ##   fasta
    ## 
    ## Alignment dimensions:
    ##   2 sequence rows; 393 position columns (259 non-gap, 134 gap) 
    ## 
    ## + attr: id, ali, call

``` r
(my_data$ali[1,] == my_data$ali[2,])
```

    ##   [1]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
    ##  [12]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
    ##  [23]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
    ##  [34]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE FALSE  TRUE  TRUE  TRUE
    ##  [45]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
    ##  [56]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE FALSE  TRUE
    ##  [67]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
    ##  [78]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
    ##  [89]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
    ## [100]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
    ## [111]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
    ## [122]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
    ## [133]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
    ## [144]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
    ## [155]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
    ## [166]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
    ## [177]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
    ## [188]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
    ## [199]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
    ## [210]  TRUE  TRUE  TRUE FALSE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
    ## [221]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
    ## [232]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
    ## [243]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
    ## [254]  TRUE  TRUE  TRUE  TRUE  TRUE FALSE FALSE FALSE FALSE FALSE FALSE
    ## [265] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ## [276] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ## [287] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ## [298] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ## [309] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ## [320] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ## [331] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ## [342] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ## [353] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ## [364] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ## [375] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ## [386] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE

``` r
start <- 41 + 8
end <- 41 - 8

seq1 <- my_data$ali[,end:start]


my_data$ali[,]
```

    ##            [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12]
    ## P53_wt     "M"  "E"  "E"  "P"  "Q"  "S"  "D"  "P"  "S"  "V"   "E"   "P"  
    ## P53_mutant "M"  "E"  "E"  "P"  "Q"  "S"  "D"  "P"  "S"  "V"   "E"   "P"  
    ##            [,13] [,14] [,15] [,16] [,17] [,18] [,19] [,20] [,21] [,22]
    ## P53_wt     "P"   "L"   "S"   "Q"   "E"   "T"   "F"   "S"   "D"   "L"  
    ## P53_mutant "P"   "L"   "S"   "Q"   "E"   "T"   "F"   "S"   "D"   "L"  
    ##            [,23] [,24] [,25] [,26] [,27] [,28] [,29] [,30] [,31] [,32]
    ## P53_wt     "W"   "K"   "L"   "L"   "P"   "E"   "N"   "N"   "V"   "L"  
    ## P53_mutant "W"   "K"   "L"   "L"   "P"   "E"   "N"   "N"   "V"   "L"  
    ##            [,33] [,34] [,35] [,36] [,37] [,38] [,39] [,40] [,41] [,42]
    ## P53_wt     "S"   "P"   "L"   "P"   "S"   "Q"   "A"   "M"   "D"   "D"  
    ## P53_mutant "S"   "P"   "L"   "P"   "S"   "Q"   "A"   "M"   "L"   "D"  
    ##            [,43] [,44] [,45] [,46] [,47] [,48] [,49] [,50] [,51] [,52]
    ## P53_wt     "L"   "M"   "L"   "S"   "P"   "D"   "D"   "I"   "E"   "Q"  
    ## P53_mutant "L"   "M"   "L"   "S"   "P"   "D"   "D"   "I"   "E"   "Q"  
    ##            [,53] [,54] [,55] [,56] [,57] [,58] [,59] [,60] [,61] [,62]
    ## P53_wt     "W"   "F"   "T"   "E"   "D"   "P"   "G"   "P"   "D"   "E"  
    ## P53_mutant "W"   "F"   "T"   "E"   "D"   "P"   "G"   "P"   "D"   "E"  
    ##            [,63] [,64] [,65] [,66] [,67] [,68] [,69] [,70] [,71] [,72]
    ## P53_wt     "A"   "P"   "R"   "M"   "P"   "E"   "A"   "A"   "P"   "P"  
    ## P53_mutant "A"   "P"   "W"   "M"   "P"   "E"   "A"   "A"   "P"   "P"  
    ##            [,73] [,74] [,75] [,76] [,77] [,78] [,79] [,80] [,81] [,82]
    ## P53_wt     "V"   "A"   "P"   "A"   "P"   "A"   "A"   "P"   "T"   "P"  
    ## P53_mutant "V"   "A"   "P"   "A"   "P"   "A"   "A"   "P"   "T"   "P"  
    ##            [,83] [,84] [,85] [,86] [,87] [,88] [,89] [,90] [,91] [,92]
    ## P53_wt     "A"   "A"   "P"   "A"   "P"   "A"   "P"   "S"   "W"   "P"  
    ## P53_mutant "A"   "A"   "P"   "A"   "P"   "A"   "P"   "S"   "W"   "P"  
    ##            [,93] [,94] [,95] [,96] [,97] [,98] [,99] [,100] [,101] [,102]
    ## P53_wt     "L"   "S"   "S"   "S"   "V"   "P"   "S"   "Q"    "K"    "T"   
    ## P53_mutant "L"   "S"   "S"   "S"   "V"   "P"   "S"   "Q"    "K"    "T"   
    ##            [,103] [,104] [,105] [,106] [,107] [,108] [,109] [,110] [,111]
    ## P53_wt     "Y"    "Q"    "G"    "S"    "Y"    "G"    "F"    "R"    "L"   
    ## P53_mutant "Y"    "Q"    "G"    "S"    "Y"    "G"    "F"    "R"    "L"   
    ##            [,112] [,113] [,114] [,115] [,116] [,117] [,118] [,119] [,120]
    ## P53_wt     "G"    "F"    "L"    "H"    "S"    "G"    "T"    "A"    "K"   
    ## P53_mutant "G"    "F"    "L"    "H"    "S"    "G"    "T"    "A"    "K"   
    ##            [,121] [,122] [,123] [,124] [,125] [,126] [,127] [,128] [,129]
    ## P53_wt     "S"    "V"    "T"    "C"    "T"    "Y"    "S"    "P"    "A"   
    ## P53_mutant "S"    "V"    "T"    "C"    "T"    "Y"    "S"    "P"    "A"   
    ##            [,130] [,131] [,132] [,133] [,134] [,135] [,136] [,137] [,138]
    ## P53_wt     "L"    "N"    "K"    "M"    "F"    "C"    "Q"    "L"    "A"   
    ## P53_mutant "L"    "N"    "K"    "M"    "F"    "C"    "Q"    "L"    "A"   
    ##            [,139] [,140] [,141] [,142] [,143] [,144] [,145] [,146] [,147]
    ## P53_wt     "K"    "T"    "C"    "P"    "V"    "Q"    "L"    "W"    "V"   
    ## P53_mutant "K"    "T"    "C"    "P"    "V"    "Q"    "L"    "W"    "V"   
    ##            [,148] [,149] [,150] [,151] [,152] [,153] [,154] [,155] [,156]
    ## P53_wt     "D"    "S"    "T"    "P"    "P"    "P"    "G"    "T"    "R"   
    ## P53_mutant "D"    "S"    "T"    "P"    "P"    "P"    "G"    "T"    "R"   
    ##            [,157] [,158] [,159] [,160] [,161] [,162] [,163] [,164] [,165]
    ## P53_wt     "V"    "R"    "A"    "M"    "A"    "I"    "Y"    "K"    "Q"   
    ## P53_mutant "V"    "R"    "A"    "M"    "A"    "I"    "Y"    "K"    "Q"   
    ##            [,166] [,167] [,168] [,169] [,170] [,171] [,172] [,173] [,174]
    ## P53_wt     "S"    "Q"    "H"    "M"    "T"    "E"    "V"    "V"    "R"   
    ## P53_mutant "S"    "Q"    "H"    "M"    "T"    "E"    "V"    "V"    "R"   
    ##            [,175] [,176] [,177] [,178] [,179] [,180] [,181] [,182] [,183]
    ## P53_wt     "R"    "C"    "P"    "H"    "H"    "E"    "R"    "C"    "S"   
    ## P53_mutant "R"    "C"    "P"    "H"    "H"    "E"    "R"    "C"    "S"   
    ##            [,184] [,185] [,186] [,187] [,188] [,189] [,190] [,191] [,192]
    ## P53_wt     "D"    "S"    "D"    "G"    "L"    "A"    "P"    "P"    "Q"   
    ## P53_mutant "D"    "S"    "D"    "G"    "L"    "A"    "P"    "P"    "Q"   
    ##            [,193] [,194] [,195] [,196] [,197] [,198] [,199] [,200] [,201]
    ## P53_wt     "H"    "L"    "I"    "R"    "V"    "E"    "G"    "N"    "L"   
    ## P53_mutant "H"    "L"    "I"    "R"    "V"    "E"    "G"    "N"    "L"   
    ##            [,202] [,203] [,204] [,205] [,206] [,207] [,208] [,209] [,210]
    ## P53_wt     "R"    "V"    "E"    "Y"    "L"    "D"    "D"    "R"    "N"   
    ## P53_mutant "R"    "V"    "E"    "Y"    "L"    "D"    "D"    "R"    "N"   
    ##            [,211] [,212] [,213] [,214] [,215] [,216] [,217] [,218] [,219]
    ## P53_wt     "T"    "F"    "R"    "H"    "S"    "V"    "V"    "V"    "P"   
    ## P53_mutant "T"    "F"    "V"    "H"    "S"    "V"    "V"    "V"    "P"   
    ##            [,220] [,221] [,222] [,223] [,224] [,225] [,226] [,227] [,228]
    ## P53_wt     "Y"    "E"    "P"    "P"    "E"    "V"    "G"    "S"    "D"   
    ## P53_mutant "Y"    "E"    "P"    "P"    "E"    "V"    "G"    "S"    "D"   
    ##            [,229] [,230] [,231] [,232] [,233] [,234] [,235] [,236] [,237]
    ## P53_wt     "C"    "T"    "T"    "I"    "H"    "Y"    "N"    "Y"    "M"   
    ## P53_mutant "C"    "T"    "T"    "I"    "H"    "Y"    "N"    "Y"    "M"   
    ##            [,238] [,239] [,240] [,241] [,242] [,243] [,244] [,245] [,246]
    ## P53_wt     "C"    "N"    "S"    "S"    "C"    "M"    "G"    "G"    "M"   
    ## P53_mutant "C"    "N"    "S"    "S"    "C"    "M"    "G"    "G"    "M"   
    ##            [,247] [,248] [,249] [,250] [,251] [,252] [,253] [,254] [,255]
    ## P53_wt     "N"    "R"    "R"    "P"    "I"    "L"    "T"    "I"    "I"   
    ## P53_mutant "N"    "R"    "R"    "P"    "I"    "L"    "T"    "I"    "I"   
    ##            [,256] [,257] [,258] [,259] [,260] [,261] [,262] [,263] [,264]
    ## P53_wt     "T"    "L"    "E"    "D"    "S"    "S"    "G"    "N"    "L"   
    ## P53_mutant "T"    "L"    "E"    "V"    "-"    "-"    "-"    "-"    "-"   
    ##            [,265] [,266] [,267] [,268] [,269] [,270] [,271] [,272] [,273]
    ## P53_wt     "L"    "G"    "R"    "N"    "S"    "F"    "E"    "V"    "R"   
    ## P53_mutant "-"    "-"    "-"    "-"    "-"    "-"    "-"    "-"    "-"   
    ##            [,274] [,275] [,276] [,277] [,278] [,279] [,280] [,281] [,282]
    ## P53_wt     "V"    "C"    "A"    "C"    "P"    "G"    "R"    "D"    "R"   
    ## P53_mutant "-"    "-"    "-"    "-"    "-"    "-"    "-"    "-"    "-"   
    ##            [,283] [,284] [,285] [,286] [,287] [,288] [,289] [,290] [,291]
    ## P53_wt     "R"    "T"    "E"    "E"    "E"    "N"    "L"    "R"    "K"   
    ## P53_mutant "-"    "-"    "-"    "-"    "-"    "-"    "-"    "-"    "-"   
    ##            [,292] [,293] [,294] [,295] [,296] [,297] [,298] [,299] [,300]
    ## P53_wt     "K"    "G"    "E"    "P"    "H"    "H"    "E"    "L"    "P"   
    ## P53_mutant "-"    "-"    "-"    "-"    "-"    "-"    "-"    "-"    "-"   
    ##            [,301] [,302] [,303] [,304] [,305] [,306] [,307] [,308] [,309]
    ## P53_wt     "P"    "G"    "S"    "T"    "K"    "R"    "A"    "L"    "P"   
    ## P53_mutant "-"    "-"    "-"    "-"    "-"    "-"    "-"    "-"    "-"   
    ##            [,310] [,311] [,312] [,313] [,314] [,315] [,316] [,317] [,318]
    ## P53_wt     "N"    "N"    "T"    "S"    "S"    "S"    "P"    "Q"    "P"   
    ## P53_mutant "-"    "-"    "-"    "-"    "-"    "-"    "-"    "-"    "-"   
    ##            [,319] [,320] [,321] [,322] [,323] [,324] [,325] [,326] [,327]
    ## P53_wt     "K"    "K"    "K"    "P"    "L"    "D"    "G"    "E"    "Y"   
    ## P53_mutant "-"    "-"    "-"    "-"    "-"    "-"    "-"    "-"    "-"   
    ##            [,328] [,329] [,330] [,331] [,332] [,333] [,334] [,335] [,336]
    ## P53_wt     "F"    "T"    "L"    "Q"    "I"    "R"    "G"    "R"    "E"   
    ## P53_mutant "-"    "-"    "-"    "-"    "-"    "-"    "-"    "-"    "-"   
    ##            [,337] [,338] [,339] [,340] [,341] [,342] [,343] [,344] [,345]
    ## P53_wt     "R"    "F"    "E"    "M"    "F"    "R"    "E"    "L"    "N"   
    ## P53_mutant "-"    "-"    "-"    "-"    "-"    "-"    "-"    "-"    "-"   
    ##            [,346] [,347] [,348] [,349] [,350] [,351] [,352] [,353] [,354]
    ## P53_wt     "E"    "A"    "L"    "E"    "L"    "K"    "D"    "A"    "Q"   
    ## P53_mutant "-"    "-"    "-"    "-"    "-"    "-"    "-"    "-"    "-"   
    ##            [,355] [,356] [,357] [,358] [,359] [,360] [,361] [,362] [,363]
    ## P53_wt     "A"    "G"    "K"    "E"    "P"    "G"    "G"    "S"    "R"   
    ## P53_mutant "-"    "-"    "-"    "-"    "-"    "-"    "-"    "-"    "-"   
    ##            [,364] [,365] [,366] [,367] [,368] [,369] [,370] [,371] [,372]
    ## P53_wt     "A"    "H"    "S"    "S"    "H"    "L"    "K"    "S"    "K"   
    ## P53_mutant "-"    "-"    "-"    "-"    "-"    "-"    "-"    "-"    "-"   
    ##            [,373] [,374] [,375] [,376] [,377] [,378] [,379] [,380] [,381]
    ## P53_wt     "K"    "G"    "Q"    "S"    "T"    "S"    "R"    "H"    "K"   
    ## P53_mutant "-"    "-"    "-"    "-"    "-"    "-"    "-"    "-"    "-"   
    ##            [,382] [,383] [,384] [,385] [,386] [,387] [,388] [,389] [,390]
    ## P53_wt     "K"    "L"    "M"    "F"    "K"    "T"    "E"    "G"    "P"   
    ## P53_mutant "-"    "-"    "-"    "-"    "-"    "-"    "-"    "-"    "-"   
    ##            [,391] [,392] [,393]
    ## P53_wt     "D"    "S"    "D"   
    ## P53_mutant "-"    "-"    "-"

``` r
ide <- conserv(my_data, method = "identity")

ide < 1 
```

    ##   [1] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ##  [12] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ##  [23] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ##  [34] FALSE FALSE FALSE FALSE FALSE FALSE FALSE  TRUE FALSE FALSE FALSE
    ##  [45] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ##  [56] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE  TRUE FALSE
    ##  [67] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ##  [78] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ##  [89] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ## [100] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ## [111] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ## [122] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ## [133] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ## [144] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ## [155] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ## [166] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ## [177] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ## [188] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ## [199] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ## [210] FALSE FALSE FALSE  TRUE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ## [221] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ## [232] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ## [243] FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE FALSE
    ## [254] FALSE FALSE FALSE FALSE FALSE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
    ## [265]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
    ## [276]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
    ## [287]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
    ## [298]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
    ## [309]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
    ## [320]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
    ## [331]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
    ## [342]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
    ## [353]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
    ## [364]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
    ## [375]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE
    ## [386]  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE  TRUE

``` r
gaps <- gap.inspect(my_data)
mismatch.inds <- which(ide < 1)

gaps.inds <- gaps$t.inds

tumor.sites <- mismatch.inds[!mismatch.inds %in% gaps.inds]


#template a[!a %in% b] 
#find the positions in mismatch.inds that are NOT in gap.inds
```

``` r
my_data$ali[,tumor.sites]
```

    ##            [,1] [,2] [,3] [,4]
    ## P53_wt     "D"  "R"  "R"  "D" 
    ## P53_mutant "L"  "W"  "V"  "V"

``` r
paste( my_data$ali[1,tumor.sites], tumor.sites,
       my_data$ali[2,tumor.sites], sep="")
```

    ## [1] "D41L"  "R65W"  "R213V" "D259V"

``` r
start.ind <- tumor.sites - 8
end.ind <-  tumor.sites + 8
my_data$ali[2, start.ind[1]: end.ind[1]]
```

    ##  [1] "S" "P" "L" "P" "S" "Q" "A" "M" "L" "D" "L" "M" "L" "S" "P" "D" "D"

``` r
tumor <- NULL

for( i in 1:length(start.ind)) {
  tumor <-  seqbind(tumor,
                    my_data$ali[2,start.ind[i]: end.ind[i]])
  
}
tumor
```

    ##        1        .      17 
    ## seq1   SPLPSQAMLDLMLSPDD
    ## seq2   DPGPDEAPWMPEAAPPV
    ## seq3   YLDDRNTFVHSVVVPYE
    ## seq4   ILTIITLEV--------
    ##                          
    ##        1        .      17 
    ## 
    ## Call:
    ##   seqbind(tumor, my_data$ali[2, start.ind[i]:end.ind[i]])
    ## 
    ## Class:
    ##   fasta
    ## 
    ## Alignment dimensions:
    ##   4 sequence rows; 17 position columns (9 non-gap, 8 gap) 
    ## 
    ## + attr: id, ali, call

``` r
#bio3d::write.fasta(ids=ids,alignment=tumor, file="subseq.fa")
```
