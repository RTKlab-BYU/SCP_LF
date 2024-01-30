library(SCnorm)
library(tidyverse)


filepaths = c("20min_January/report.pg_matrix.tsv")
MS2_filenames = c("20min_January/report-first-pass.pg_matrix.tsv")
group_names = list()
group_names[[filepaths[1]]] = c("Colony 1.3", "Colony 1.6")



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
yms2 = list()
j = 1
for (filepath in filepaths){
  x = DIANN.AbundanceMatrix(proteinFile = filepath)
  x2 = DIANN.AbundanceMatrix(proteinFile = MS2_filenames[j])
  
  i = 1
  for(eachGroup in group_names[[filepath]]) {
    y[[eachGroup]] = x %>%
      DIANN.FilterByName(filterIn = group_filters[[filepath]][i], filterOut = out_filter)
    
    yms2[[eachGroup]] = x2 %>%
      DIANN.FilterByName(filterIn = group_filters[[filepath]][i], filterOut = out_filter)
    i = i + 1
  }
  j = j + 1
}



################Make Matrix################

Conditions = c()
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
    if(i==1){
      DIANN_report_matrix = current
    } else {
      DIANN_report_matrix = current %>%
        merge(DIANN_report_matrix)
    }
    Conditions = c(Conditions,rep(eachGroup,length(mycols)))
    i = i + 1
    
  }
}

DIANN_report_matrix = DIANN_report_matrix %>%
  t()
    
################Normalize################


mySCData = SingleCellExperiment::SingleCellExperiment(assays = list('counts' = DIANN_report_matrix))


exampleGC <- runif(dim(mySCData)[1], 0, 1)
names(exampleGC) <- rownames(mySCData)

DataNorm <- SCnorm(mySCData, Conditions,
                   PrintProgressPlots = TRUE,
                     FilterCellNum = 10, withinSample = exampleGC)