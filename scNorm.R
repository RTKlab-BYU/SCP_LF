library(SCnorm)
library(tidyverse)
library(dplyr)
library(limma)


filepaths = c("20min_January/report-first-pass.pg_matrix.tsv",
 "20min_January/report.pg_matrix.tsv",
"10min_January/report-first-pass.pg_matrix.tsv",
"10min_January/report.pg_matrix.tsv",

"20min_January/report-first-pass.pr_matrix.tsv",
"20min_January/report.pr_matrix.tsv",
"10min_January/report-first-pass.pr_matrix.tsv",
"10min_January/report.pr_matrix.tsv")

group_names = list()
group_names[[filepaths[1]]] = c("Colony 1.3", "Colony 1.6")
group_names[[filepaths[2]]] = c("Colony 1.3", "Colony 1.6")
group_names[[filepaths[3]]] = c("Colony 1.3", "Colony 1.6")
group_names[[filepaths[4]]] = c("Colony 1.3", "Colony 1.6")
group_names[[filepaths[5]]] = c("Colony 1.3", "Colony 1.6")
group_names[[filepaths[6]]] = c("Colony 1.3", "Colony 1.6")
group_names[[filepaths[7]]] = c("Colony 1.3", "Colony 1.6")
group_names[[filepaths[8]]] = c("Colony 1.3", "Colony 1.6")

#isProtein = TRUE
isProteins = c(TRUE,TRUE,TRUE,TRUE,FALSE,FALSE,FALSE,FALSE)


group_filters = list()
group_filters[[filepaths[1]]] = c("Col1-3","Col1-6")
group_filters[[filepaths[2]]] = c("Col1-3","Col1-6")
group_filters[[filepaths[3]]] = c("Col1-3","Col1-6")
group_filters[[filepaths[4]]] = c("Col1-3","Col1-6")
group_filters[[filepaths[5]]] = c("Col1-3","Col1-6")
group_filters[[filepaths[6]]] = c("Col1-3","Col1-6")
group_filters[[filepaths[7]]] = c("Col1-3","Col1-6")
group_filters[[filepaths[8]]] = c("Col1-3","Col1-6")

out_filter = c("N3","N4","H6","M3","L3","J10","J11")

################Functions################

DIANN.AbundanceMatrix = function(proteinFile, isProtein=TRUE, eachFilenameConvention = ".d") {
  if (isProtein){
    y = read_tsv(proteinFile) %>%
      filter(!grepl("contam", Protein.Group),!grepl("Keratin",First.Protein.Description)) %>%
      dplyr::select(contains(eachFilenameConvention), contains( "Protein.Ids")) %>%
      mutate(Accession = Protein.Ids) %>%
      select(-Protein.Ids)
  } else {
    y = read_tsv(proteinFile) %>%
      filter(!grepl("contam", Protein.Group)) %>%
      dplyr::select(contains(eachFilenameConvention), contains("Stripped.Sequence")) %>%
      mutate(`Annotated Sequence` = Stripped.Sequence) %>%
      select(-Stripped.Sequence) 
  }
  y
}

DIANN.FilterByName = function(abundances, filterIn, filterOut, isProtein = TRUE){
  if(isProtein){
    name = "Accession"
  } else {
    name = "Annotated Sequence"
  }
  for(eachFilter in filterIn) {
    abundances = dplyr::select_if(abundances, 
                                             grepl(paste0(name,"|",eachFilter),
                                                   colnames(abundances)))
  }
  for(eachFilter in filterOut) {
    abundances = dplyr::select_if(abundances, 
                                             !grepl(eachFilter,
                                                    colnames(abundances)))
  }
  abundances
}
################Import##############


j = 1
for (filepath in filepaths){
  y = list()
  isProtein = isProteins[j]
  x = DIANN.AbundanceMatrix(proteinFile = filepath, isProtein = isProtein)
  OUTPUT_FILENAME = paste0(gsub("\\.tsv","",filepath),"_normalized.tsv")
  
  i = 1  
  for(eachGroup in group_names[[filepath]]) {
    y[[eachGroup]] = x %>%
      DIANN.FilterByName(filterIn = group_filters[[filepath]][i], filterOut = out_filter, isProtein = isProtein)
    i = i + 1
    
  }
  


  remove(x)

################Make Matrix################
  if(isProtein){
    name = "Accession"
    new_name = "Protein.Ids"
  } else {
    name = "Annotated Sequence"
    new_name = "Stripped.Sequence"
  }
  
  Conditions = c()
  maxnumcells = 0
  i = 1
  
 
  for(eachGroup in group_names[[filepath]]) {
    mycols = colnames(y[[eachGroup]])
    mycols = mycols[mycols != name]
    myrows = y[[eachGroup]][[name]]
    current = y[[eachGroup]][mycols] %>%
      as.matrix()
    colnames(current) = mycols
    rownames(current) = myrows
    current = current[rowSums(is.na(current)) != ncol(current),]
    if(i==1){
      DIANN_report_matrix = current
    } else {
      DIANN_report_matrix = current %>%
        merge(DIANN_report_matrix,by="row.names")
    }
    numcells = length(mycols)
      maxnumcells = numcells
    Conditions = c(Conditions,rep(eachGroup,numcells))
    i = i + 1
    
  }
  
  
  newcols = DIANN_report_matrix[,1]
  DIANN_report_matrix = DIANN_report_matrix[,-1]
  rownames(DIANN_report_matrix) = newcols
  DIANN_report_matrix[is.na(DIANN_report_matrix)] = 0
  
  remove(current)
  remove(y)
 
################Normalize################


  mySCData = SingleCellExperiment::SingleCellExperiment(assays = list('counts' = DIANN_report_matrix))
  
  
  exampleGC = runif(dim(mySCData)[1], 0, 1)
  names(exampleGC) = rownames(mySCData)
  minCells = floor(maxnumcells*0.7) #30% missing values and above will be removed
  
  DataNorm = SCnorm(mySCData, Conditions,
                     PrintProgressPlots = TRUE,
                       FilterCellNum = minCells, withinSample = exampleGC,
                     NCores=1)
  
  NormedData = SingleCellExperiment::normcounts(DataNorm)
  
  NormedData[NormedData == 0] = NA
  
  proteinIDs = row.names(NormedData)
  
  NormedData =as.data.frame(NormedData)
  
  NormedData[new_name] = proteinIDs
  
  write_tsv(NormedData, OUTPUT_FILENAME)
  j = j + 1
}
