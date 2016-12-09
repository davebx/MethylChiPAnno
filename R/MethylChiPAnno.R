##' @importMethodsFrom Biobase, channel, classVersion, experimentData, featureData, featureNames, "featureNames<-", isCurrent, pData, "pData<-", phenoData, protocolData, sampleNames, "sampleNames<-", storageMode, exprs, fData, varMetadata, "varMetadata<-", "featureData<-", AnnotatedDataFrame
##' @importFrom minfi pData
##' @importFrom minfi getBeta
##' @importFrom minfi GenomicMethylSet
##' @importFrom minfi makeGenomicRatioSetFromMatrix
##' @importFrom base rep
##' @importFrom base c
##' @importFrom stats model.matrix
##' @importFrom minfi bumphunter
##' @importFrom utils read.delim
##' @importFrom base intersect
##' @importFrom ChIPpeakAnno findOverlapsOfPeaks
##' @importFrom TxDb.Hsapiens.UCSC.hg19.knownGene TxDb.Hsapiens.UCSC.hg19.knownGene
##' @importFrom ChIPseeker annotatePeak
##' @importFrom base data.frame
##' @importFrom utils write.table
##' @import GenomicRanges
##' @import SummarizedExperiment
##' @import S4Vectors
##' @import IRanges
##' @import BiocGenerics
##' @import minfi
##' @exportClasses RGChannelSet, RGChannelSetExtended, MethylSet, RatioSet, GenomicMethylSet, GenomicRatioSet, IlluminaMethylationManifest, IlluminaMethylationAnnotation
##' @exportMethods show, initialize, getBeta, getM, getCN, getMeth, getUnmeth, getManifest, annotation, preprocessMethod, combine, sampleNames, featureNames, pData, mapToGenome, ratioConvert, bumphunter

# Import required libraries
  require("data.table",quietly = TRUE)
  require("minfi",quietly = TRUE)
  require("ChIPseeker",quietly = TRUE)
  require("ChIPpeakAnno",quietly = TRUE)
  require("TxDb.Hsapiens.UCSC.hg19.knownGene",quietly = TRUE)
  require("FDb.InfiniumMethylation.hg19",quietly = TRUE)
  require("GEOquery",quietly = TRUE)
  require("rtracklayer",quietly = TRUE)
  require("BiocGenerics",quietly = TRUE)
  #library("MethylChiPAnno", quietly=TRUE, warn.conflicts=FALSE,verbose = FALSE)

# define parameters
  clusterSize=2
  cutoff=0.2
  platform_id='HM450'
  genome_id='hg19'

# cmd arguments
#  args <- commandArgs(trailingOnly = TRUE)
#  methyl_file = args[1]
#  ChiPseq_file =args[2]
#  output_file = args[3]
  methyl_file <- "test-data/input.csv"
  ChiPseq_file <- "test-data/Galaxy3.bed"

  options(warn=-1)

# read uploaded data
  TAB=read.csv(methyl_file)
  ChiPseq=read.delim(ChiPseq_file)

# preprocess data to pass into bumphunter function
  ChiPseq <- GRanges(seqnames= ChiPseq[,1],
                  ranges=IRanges
                  (start=ChiPseq[,2],
                  end= ChiPseq[,3]))
  if(is.null(TAB)){
    stop("Must specify input files")
  }else{
    mysamples <- lapply(TAB$ID,function(x) getGEO(x))
  }
    s0 <- lapply(mysamples,Table)
  id_ref<-lapply(s0,function(x)x$ID_REF)

  if(length(unique(id_ref)) != 1) {
    stop("Error different ID_REF for samples")
  } else if (is.null((unlist(unique(id_ref))))) {
    stop("NO GSM data avaliable")
  } else {
    values<-do.call("cbind",lapply(s0,function(x)x$VALUE))
    colnames(values)=TAB$ID
    rownames(values)=id_ref[[1]]
    cg <-  rownames(values)
    probe <- c(cg)
    hm450.hg19 <- getPlatform(platform=platform_id, genome=genome_id)
    probe.info <- hm450.hg19[probe]
    f <- data.table(probe=names(probe.info),CHR=as.data.frame(probe.info@seqnames)$value,
                    BP=as.numeric(probe.info@elementMetadata$probeStart))
    designMatrix <- model.matrix(~ TAB$Phenotype)

# bumphunter Run with processed data
    DMR <- bumphunter(values, design = designMatrix,
                      pos=f$BP,cutoff=cutoff,chr=f$CHR)

# choose DMR's of a certain length threshold
    MAT <- DMR$table[which(DMR$table$L>=clusterSize),]
    METH <- GRanges(seqnames=MAT$chr,
                    ranges=IRanges
                    (start=MAT$start,
                    end=MAT$end),
                    value_pos=MAT$value)

#finding/counting overlaps
    peaks<-findOverlaps(METH,ChiPseq,maxgap=5000)

# save result, which contains overlaps between Methyl and ChiPseq data
    p <- as.data.frame(ChiPseq[subjectHits(peaks),])
#    write.csv(p,  output_file, row.names=FALSE)
  }

