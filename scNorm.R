library(SCnorm)
library(tidyverse)
library(dplyr)
library(limma)


filepaths = c("20min_January/report.pg_matrix.tsv")
group_names = list()
group_names[[filepaths[1]]] = c("Colony 1.3", "Colony 1.6")


filepath = filepaths[1]
OUTPUT_FILENAME = paste0(gsub("\\.tsv","",filepath),"_normalized.tsv")

group_filters = list()
group_filters[[filepaths[1]]] = c("Col1-3","Col1-6")
out_filter = c("M@di")

################Functions################

DIANN.AbundanceMatrix = function(proteinFile, isProtein=TRUE, eachFilenameConvention = ".d") {
  if (isProtein){
    y = read_tsv(proteinFile) %>%
      filter(!grepl("contam", Protein.Group)) %>%
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

DIANN.FilterByName = function(abundances, filterIn, filterOut){
  name = "Accession"
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

y = list()

for (filepath in filepaths){
  x = DIANN.AbundanceMatrix(proteinFile = filepath)

  i = 1  
  for(eachGroup in group_names[[filepath]]) {
    y[[eachGroup]] = x %>%
      DIANN.FilterByName(filterIn = group_filters[[filepath]][i], filterOut = out_filter)
    i = i + 1
    
  }
  
}

remove(x)

################Make Matrix################

Conditions = c()
maxnumcells = 0
i = 1

for (filepath in filepaths){
  for(eachGroup in group_names[[filepath]]) {
    mycols = colnames(y[[eachGroup]])
    mycols = mycols[mycols != "Accession"]
    myrows = y[[eachGroup]][["Accession"]]
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
                   NCores=8)

NormedData = SingleCellExperiment::normcounts(DataNorm)

NormedData[NormedData == 0] = NA

write_tsv(NormedData, OUTPUT_FILENAME)
