library(FlowSorted.Blood.EPIC) #change .EPIC to .450k
#http://bioconductor.org/packages/release/data/experiment/html/FlowSorted.Blood.450k.html
library(minfi)
library(wateRmelon)
library(DMRcate)
library(parallel)
library(mCSEA)
library(annotatr)
library(dplyr)
library(FlowSorted.Blood.EPIC)
library(sva)
library(corrplot)
library(mCSEA)
library(purrr)
library(clusterProfiler)
library(glue)
library(org.Hs.eg.db)
library(enrichplot)
library(readr)
library(gdata)


location.based.pc<-function(beta.obj,annocpgs){
  
  beta.obj<-as.matrix(beta.obj)
  beta.obj<-beta.obj[which(rownames(beta.obj) %in% annocpgs),]
  theresult<-princomp(x=beta.obj,cor=T)
  gc()
  return(theresult)
}

rankProbes2<-function (methData, pheno = NULL, paired = FALSE, explanatory = 1,
                       covariates = c(), pairColumn = c(), caseGroup = 1, refGroup = 2,
                       continuous = NULL, typeInput = "beta", typeAnalysis = "M", corfit=NULL)
{
  if (!any(class(methData) == "data.frame" | class(methData) ==
           "matrix" | class(methData) == "SummarizedExperiment" |
           class(methData) == "RangedSummarizedExperiment")) {
    stop("methData must be a data frame, a matrix or a SummarizedExperiment\n            object")
  }
  if (!any(class(pheno) == "data.frame" | class(pheno) == "matrix" |
           is.null(pheno))) {
    stop("pheno must be a data frame, a matrix or NULL")
  }
  if (class(paired) != "logical") {
    stop("paired must be a logical object (TRUE/FALSE)")
  }
  if (paired) {
    if (!any(class(pairColumn) != "character" | !is.numeric(pairColumn))) {
      stop("pairColumn must be a character or numeric object")
    }
  }
  if (!any(class(explanatory) != "character" | !is.numeric(explanatory))) {
    stop("explanatory must be a character or numeric object")
  }
  if (!any(class(covariates) == "character" | class(covariates) ==
           "list" | is.null(covariates))) {
    stop("covariates must be a character vector, a list or NULL")
  }
  if (!any(class(caseGroup) != "character" | class(caseGroup) !=
           "numeric")) {
    stop("caseGroup must be a character or numeric object")
  }
  if (!any(class(refGroup) != "character" | class(refGroup) !=
           "numeric")) {
    stop("refGroup must be a character or numeric object")
  }
  if (!any(class(continuous) == "character" | class(continuous) ==
           "list" | is.null(continuous))) {
    stop("continuous must be a character vector, a list or NULL")
  }
  if (class(typeInput) != "character") {
    stop("typeInput must be a character object")
  }
  if (class(typeAnalysis) != "character") {
    stop("typeAnalysis must be a character object")
  }
  if (class(methData) == "SummarizedExperiment" | class(methData) ==
      "RangedSummarizedExperiment") {
    if (is.null(pheno)) {
      pheno <- SummarizedExperiment::colData(methData)
    }
    methData <- SummarizedExperiment::assay(methData)
  }
  else {
    if (is.null(pheno)) {
      stop("If methData is not a SummarizedExperiment, you must provide\n                pheno parameter")
    }
  }
  if (is.null(continuous)) {
    continuous <- c()
    categorical <- colnames(pheno)
  }
  else {
    if (class(continuous) != "character") {
      continuous <- colnames(pheno)[continuous]
    }
    categorical <- setdiff(colnames(pheno), continuous)
  }
  for (column in colnames(pheno)) {
    if (column %in% categorical) {
      if (typeof(pheno[, column]) == "list") {
        message(paste(column, "variable skipped due to it is a list"))
        pheno[, -which(names(pheno) == column)]
      }
      else {
        pheno[, column] <- factor(pheno[, column])
      }
    }
    else {
      pheno[, column] <- as.numeric(as.character(pheno[,
                                                       column]))
    }
  }
  typeInput <- match.arg(typeInput, c("beta", "M"))
  typeAnalysis <- match.arg(typeAnalysis, c("M", "beta"))
  if (class(explanatory) == "numeric") {
    explanatory <- colnames(pheno)[explanatory]
  }
  if (is.numeric(covariates)) {
    covariates <- colnames(pheno)[covariates]
  }
  if (is.numeric(pairColumn)) {
    pairColumn <- colnames(pheno)[pairColumn]
  }
  if (class(caseGroup) == "numeric") {
    caseGroup <- levels(pheno[, explanatory])[caseGroup]
  }
  if (class(refGroup) == "numeric") {
    refGroup <- levels(pheno[, explanatory])[refGroup]
  }
  if (length(intersect(explanatory, covariates)) > 0) {
    stop("You specified some variable(s) as both explanatory and covariate")
  }
  samples2Include <- rownames(pheno)[pheno[, explanatory] %in%
                                       c(caseGroup, refGroup)]
  methData <- methData[, samples2Include]
  pheno <- data.frame(pheno[samples2Include, c(explanatory,
                                               covariates, pairColumn)])
  colnames(pheno) <- c(explanatory, covariates, pairColumn)
  pheno <- droplevels(pheno)
  if (typeInput == "beta") {
    if (any(min(methData, na.rm = TRUE) < 0 | max(methData,
                                                  na.rm = TRUE) > 1)) {
      warning("Introduced beta-values are not between 0 and 1. Are you\n                    sure these are not M-values?")
    }
    if (typeAnalysis == "beta") {
      dataLimma <- methData
    }
    else {
      message("Transforming beta-values to M-values")
      dataLimma <- log2(methData) - log2(1 - methData)
    }
  }
  else {
    if (min(methData, na.rm = TRUE) >= 0 && max(methData,
                                                na.rm = TRUE) <= 1) {
      warning("Introduced M-values are between 0 and 1. Are you sure these\n                    are not beta-values?")
    }
    if (typeAnalysis == "beta") {
      message("Transforming M-values to beta-values")
      dataLimma <- 2^(methData)/(1 + 2^(methData))
    }
    else {
      dataLimma <- methData
    }
  }
  message("Calculating linear model...")
  message(paste("\tExplanatory variable:", explanatory))
  if (is.factor(pheno[, explanatory])) {
    message(paste("\tCase group:", caseGroup))
    message(paste("\tReference group:", refGroup))
    message(paste("\tTotal samples:", ncol(methData)))
    pheno[, explanatory] <- relevel(pheno[, explanatory],
                                    ref = refGroup)
  }
  if (is.null(covariates) && !paired) {
    message("\tCovariates: None")
    message(paste("\tCategorical variables:", paste(categorical,
                                                    collapse = " ")))
    if (length(continuous) > 0) {
      message(paste("\tContinuous variables:", paste(continuous,
                                                     collapse = " ")))
    }
    else {
      message("\tContinuous variables: None")
    }
    model <- model.matrix(~get(explanatory), data = pheno)
    corfit <- duplicateCorrelation(dataLimma,design=model,block=corfit)
  }
  else {
    message(paste("\tPaired analysis using", pairColumn,
                  "column"))
    message(paste("\tCovariates:", paste(covariates, collapse = " ")))
    message(paste("\tCategorical variables:", paste(categorical,
                                                    collapse = " ")))
    message(paste("\tContinuous variables:", paste(continuous,
                                                   collapse = " ")))
    model <- model.matrix(~., data = pheno)
    corfit <- duplicateCorrelation(dataLimma,design=model,block=corfit)
  }
  linearModel <- limma::eBayes(limma::lmFit(dataLimma, model,correlation=corfit$consensus.correlation))
  tValues <- linearModel[["t"]][, 2]
  pvals=linearModel[["p.value"]][,2]
  outdf=data.frame(tval=tValues,pvals=pvals)
  return(outdf)
}

make_dmr=function(pr_asc){
  #print(pr_asc)
  #beg doesn't ACTUALLY necessarily mean the beginning of the dmr
  #its mixed up/arbitrary order
  beg=pr_asc[1]
  chr=epano[epano$Name==beg,1]
  prange=epano[epano$Name %in% pr_asc,2]
  pbeg=min(prange)
  pend=max(prange)
  df=data.frame(chr=chr,start=pbeg,end=pend)
  return(df)
}



prettify_dmrs<-function(pasc){
  dfa=map(prom_asc,~safe_mdmr(.x))
  outlist=list()
  i=1
  for(i in 1:length(dfa)){
    if(!is.null(dfa[i][[1]]$result)){
      res=dfa[i][[1]]$result
      res=cbind(res,gname=names(dfa[i]))
      outlist[[i]]=res
    }
    i=i+1
  }
  res=bind_rows(outlist)
  return(res)
}

annotatr_wrap<-function(GR,annot_out,overlaps_out_fig,overlaps_res){
  #Select annotations for intersection with regions
  # Note inclusion of custom annotation, and use of shortcuts
  annots = c('hg19_cpgs', 'hg19_basicgenes', 'hg19_genes_intergenic',
             'hg19_genes_intronexonboundaries',"hg19_Gm12878-chromatin")
  
  # Build the annotations (a single GRanges object)
  annotations = build_annotations(genome = 'hg19', annotations = annots)
  
  # Intersect the regions we read in with the annotations
  dm_annotated = annotate_regions(
    regions = GR,
    annotations = annotations,
    ignore.strand = TRUE,
    quiet = FALSE)
  # A GRanges object is returned
  #print(dm_annotated)
  
  #need to specifiy hg19 genome for randomize to work (makes sense)
  genome(GR)="hg19"
  # Randomize the input regions
  dm_random_regions = randomize_regions(
    regions = GR,
    allow.overlaps = TRUE,
    per.chromosome = TRUE)
  
  # Annotate the random regions using the same annotations as above
  # These will be used in later functions
  dm_random_annotated = annotate_regions(
    regions = dm_random_regions,
    annotations = annotations,
    ignore.strand = TRUE,
    quiet = TRUE)
  
  df_annot=data.frame(dm_annotated)
  write_csv(df_annot,annot_out)
  
  dm_vs_kg_annotations = plot_annotation(
    annotated_regions = dm_annotated,
    x_label = 'knownGene Annotations',
    plot_title="McSEA DMRs vs Random",
    annotated_random = dm_random_annotated,
    y_label = 'Count')
   dm_vs_kg_annotations=dm_vs_kg_annotations+
    scale_x_discrete(guide=guide_axis(angle=55)) 
  print(dm_vs_kg_annotations)
  ggsave(overlaps_out_fig,width=11,height=6)
  
  dm_annsum = summarize_annotations(
    annotated_regions = dm_annotated,
    quiet = TRUE)
  print(dm_annsum)
  write_csv(dm_annsum,overlaps_res)
  an_red=
  return(dm_annotated)
}



massbar=function(gene,plotname){
  enr=c("MF","CC","BP")
  for(cat in enr){
    ego=enrichGO(gene,OrgDb=org.Hs.eg.db,keyType="SYMBOL",ont=cat)
    odf=data.frame(ego)
    if(nrow(odf)==0) next
    p=barplot(ego)
    o=glue("{plotname}_{cat}.pdf")
    ggsave(o,width=10,height=10)
    o=glue("{plotname}_{cat}.csv")
    write_csv(odf,o)
  }
}

allbar=function(gene,plotname){
  ego=enrichGO(gene,OrgDb=org.Hs.eg.db,keyType='SYMBOL',ont='ALL',pool=TRUE)
  
  odf=data.frame(ego)
  if (nrow(odf)==0) return
  #ego=simplify(ego)
  p=barplot(ego)
  o=glue("{plotname}_all_ont.pdf")
  ggsave(o,width=10,height=10)
  o=glue("plotname}_all_ont.csv")
  odf$cat=rownames(odf)
  write_csv(odf,o)
}

allcnet=function(gene,foldchange,plotname){
  ego=enrichGO(gene,OrgDb=org.Hs.eg.db,keyType='SYMBOL',ont='ALL',pool=TRUE)
  
  odf=data.frame(ego)
  if (nrow(odf)==0) return
  #ego=simplify(ego)
  p=cnetplot(ego,foldChange=foldchange)
  o=glue("{plotname}_all_ont.pdf")
  ggsave(o,width=10,height=10)
  o=glue("plotname}_all_ont.csv")
  odf$cat=rownames(odf)
  write_csv(odf,o)
}

masscnet=function(gene,fchange,plotname){
  enr=c("MF","CC","BP")
  for(cat in enr){
    ego=enrichGO(gene,OrgDb=org.Hs.eg.db,keyType="SYMBOL",ont=cat)
    odf=as.data.frame(ego)
    
    if(nrow(odf)==0) next
    p=cnetplot(ego,foldChange=fchange)
    o=glue("{plotname}_{cat}.pdf")
    ggsave(o,width=10,height=10)
    o=glue("{plotname}_{cat}.csv")
    write_csv(odf,o)
  }
}

