library(data.table)
library(dplyr)
library(tidyr)
library(emdbook)


cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

options(stringsAsFactors=F)

data.dir <- "/Users/db2175/NonDData/Data/Carmi/"

array.genotypes <- c("AA", "AB", "BB")

cnv.genotypes <- c("N", "A", "B", "AA", "AB", "BB", "AAA", "AAB", "ABB", "BBB", "AAAA", "AAAB", "AABB", "ABBB", "BBBB")   

array.error.matrix <- matrix(
  c(0.999, 0.00006, 0.00094,
    0.00006,  0.999, 0.00094,
    0.0105,  0.0105,  0.979),
  byrow=T, nrow=3, ncol=3,
  dimnames=list(array.genotypes, array.genotypes)
)

ref.allele.probability <- c(N=0.54, A=0.998, B=0.011, AA=0.998, BB=0.011, 
                            AB=0.54, AAA=0.998, BBB=0.011, AAB=0.8, ABB=0.4, 
                       AAAA=0.998, BBBB=0.011, ABBB=0.3, AABB=0.54, AAAB=0.8)

#ref.allele.probability <- c(AA=0.998, AB=0.54, BB=0.011)

bb.params <- c(prob=0.554, theta=0.473)

# pd should be a data.table or data.frame with columns Type (Sequencing or Array)
# and GenData
# returns 
util.GenerateDataProbabilities <- function(pd, make.long=F){
  if (make.long){
    pd <- gather(pd, Sample, GenData, -CHROM, -POS) %>%
      mutate(Type=ifelse(GenData %in% array.genotypes, "Array", "Sequencing"))
  }
  pd <- as.data.table(pd)
  seq.data <- pd[Type=="Sequencing" & !is.na(GenData), ]
  array.data <- pd[Type=="Array" & !is.na(GenData), ]
  na.data <- pd[is.na(GenData), ]
  # want to calculate probability of observed data given
  # each genotype
  for (true.genotype in array.genotypes){
    array.genotype.probs <- array.error.matrix[true.genotype, ]
    array.data[, (true.genotype):=array.genotype.probs[array.data$GenData]]
  }
  for (true.genotype in cnv.genotypes){
    spl <- strsplit(seq.data$GenData, split="/", fixed=T)
    ref <- as.numeric(unlist(spl)[c(T,F)])
    alt <- as.numeric(unlist(spl)[c(F,T)])
    depth <- ref + alt
    
    seq.data[, (true.genotype):=dbetabinom(x=ref, size=depth, 
                                             #prob=bb.params['prob'], 
                                             prob=ref.allele.probability[true.genotype],
                                             theta=bb.params['theta'])]
  }
  #array.data[is.na(AA), `:=`(AA=1, AB=1, BB=1)]
  #na.data[, `:=`(AA=1, AB=1, BB=1)]
  pd <- bind_rows(seq.data, array.data, na.data)
}

util.GetRecombinationProbability <- function(cM){
  cm.diffs <- cM[2:length(cM)] - 
    cM[1:(length(cM) - 1)]
  probs <- (1 - exp(-2*cm.diffs/100))/2
}

# positions should be a numeric vector
# returns probabilities for males and female of recombination between
# each successive pair of positions
# in list with elements male.recomb.probs and female.recomb.probs
util.GetRecombinationProbabilities <- function(chr, positions){
  gmap.male <- read.table(paste0("/Users/db2175/Dropbox/Daniel/Biostatistics/Carmi/Refined_genetic_map_b37/male_", chr, ".txt"), header=T, stringsAsFactors=F)
  gmap.female <- read.table(paste0("/Users/db2175/Dropbox/Daniel/Biostatistics/Carmi/Refined_genetic_map_b37/female_", chr, ".txt"), header=T, stringsAsFactors=F)
  male.interp <- approx(gmap.male$pos, gmap.male$cM, xout=positions)
  female.interp <- approx(gmap.female$pos, gmap.female$cM, xout=positions)
  male.recomb.probs <- util.GetRecombinationProbability(male.interp$y)
  female.recomb.probs <- util.GetRecombinationProbability(female.interp$y)
  male.recomb.probs[is.na(male.recomb.probs)] <- min(male.recomb.probs, na.rm=T)
  female.recomb.probs[is.na(female.recomb.probs)] <- min(female.recomb.probs, na.rm=T)
  return(list(male.recomb.probs=male.recomb.probs, 
              female.recomb.probs=female.recomb.probs))
}

# from in rows
# chr is 1,2,3...
# pos is numeric

GetChromosomeSizes <- function(){
  chrom.sizes <- read.table("hg19.chrom.sizes", 
                            col.names=c("chr", "length")) %>%
    filter(chr %in% paste0("chr", c(1:22, "X", "Y") )) %>% 
    mutate(chr=factor(chr, levels=c(paste0("chr", c(1:22, "X", "Y"))))) %>%
    #mutate(chr=as.numeric(gsub("chr", "", chr))) %>%
    arrange(chr) %>% mutate(chr=as.character(chr))
}

GetGenomicCoordinates <- function(chr, pos, buffer=25e+06){
  chrom.sizes <- GetChromosomeSizes()
  chrom.sizes$start <- c(0, cumsum(chrom.sizes$length[1:23] + buffer))
  start.dict <- chrom.sizes$start
  names(start.dict) <- chrom.sizes$chr
  coord <- start.dict[chr] + pos
}

GetMidChromosomalCoordinates <- function(buffer=25e+06){
  chrom.sizes <- GetChromosomeSizes()
  GetGenomicCoordinates(chr=paste0("chr", c(1:22, "X", "Y")), pos=chrom.sizes$length / 2, buffer=buffer)
}