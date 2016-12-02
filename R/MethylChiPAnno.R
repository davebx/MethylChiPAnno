# Send R errors to stderr
options(show.error.messages = F, error = function(){cat(geterrmessage(), file = stderr()); q("no", 1, F)})

# Avoid crashing Galaxy with an UTF8 error on German LC settings
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

# Import library
library("getopt")

options(stringAsfactors = FALSE, useFancyQuotes = FALSE)

  options(warn=-1)
  library("data.table",quietly = TRUE)
  library("minfi",quietly = TRUE)
  library("ChIPseeker",quietly = TRUE)
  library("ChIPpeakAnno",quietly = TRUE)
  library("TxDb.Hsapiens.UCSC.hg19.knownGene",quietly = TRUE)
  library("FDb.InfiniumMethylation.hg19",quietly = TRUE)
  library("GEOquery",quietly = TRUE)
  library("rtracklayer",quietly = TRUE)
  library("BiocGenerics",quietly = TRUE)
  
 
  
  clusterSize=2
  cutoff=0.2
  platform_id='HM450'
  genome_id='hg19'
  
  args <- commandArgs(trailingOnly = TRUE)
  methyl_file = args[1]
  ChiPseq_file =args[2]
  output_file = args[3]
  options(warn=-1)
  TAB=read.csv(methyl_file)
  ChiPseq=read.delim(ChiPseq_file)
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
    
    DMR <- bumphunter(values, design = designMatrix,
                      pos=f$BP,cutoff=cutoff,chr=f$CHR)
    
    MAT <- DMR$table[which(DMR$table$L>=clusterSize),]
    METH <- GRanges(seqnames=MAT$chr,
                    ranges=IRanges
                    (start=MAT$start,
                    end=MAT$end),
                    value_pos=MAT$value)
    ?findOverlapsOfPeaks
    peaks<-findOverlaps(METH,ChiPseq,maxgap=5000)
    p <- as.data.frame(ChiPseq[subjectHits(peaks),])
    write.csv(p,  output_file, row.names=FALSE)
  }
 


