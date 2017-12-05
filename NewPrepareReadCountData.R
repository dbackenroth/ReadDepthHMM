options(stringsAsFactors=F)
options(scipen=30)

library(data.table)
library(dplyr)

read.locs.dir <- "/Users/db2175/NonDData/Data/CarmiBAMs/PGDlympho_data/PGDlympho_data/Readstarts"
sequence.stats.dir <- "/Users/db2175/NonDData/Data/CarmiBAMs/PGDlympho_data/PGDlympho_data/SequenceStatFiles/"

all.chrs <- paste0("chr", c(1:22, "X", "Y"))

# returns three column data.table
# with columns chr (chr1, ...), pos and MapQ
# for sample.id
# all reads mapping to chrs 1:22, X and Y

GetReadStartsFilename <- function(sample.id){
  readcount.files <- list.files(read.locs.dir, full.names=T)
  raw.f <- readcount.files[grepl(sample.id, readcount.files, fixed=T)]
  raw.f <- raw.f[endsWith(raw.f, ".gz")]
}

GetReadStarts <- function(sample.id="2M-3"){
  library(data.table)
  raw.f <- GetReadStartsFilename(sample.id=sample.id)
  res <- fread(paste0("gunzip -c ", raw.f), header=F)
  ProcessReadStarts(res)
}

ProcessReadStarts <- function(res){
  setnames(res, c("chr", "pos", "MapQ"))
  res[, `:=`(chr=gsub("chr", "", chr))]
  res[, `:=`(chr=paste0("chr", chr))]
  res <- res[chr %in% paste0("chr", c(1:22, "X", "Y", "MT", "M")), ]
  setkey(res, chr, pos)
  res
}

# chr should have "chr" prefix
GetRegionWithTabix <- function(sample="2M-3", chr="chr1", start=1, end=1000000){
  filename <- GetReadStartsFilename(sample)
  fi <- FamilyInfoWithSex() %>% as.data.table()
  if (!fi[Sample==sample, chr.prefix]) chr <- gsub("chr", "", chr)
  res <- system(paste0("tabix ", filename, " ", chr, ":", start, "-", end), intern=T)
  tab <- fread(paste0(res, collapse="\n"), col.names=)
  ProcessReadStarts(tab)
}

GetRegionWithTabixSamples <- function(samples, chr, start, end, 
                                      bins=NULL){
  l <- list()
  for (s in samples){
    res <- GetRegionWithTabix(s, chr, start, end) %>% mutate(Sample=s)
    res <- as.data.table(res)
    if (!is.null(bins)){
      bins$BIN <- 1:nrow(bins)
      breaks <- sort(unique(c(bins$start, bins$end)))
      breaks[1] <- breaks[1] - 100000
      breaks[length(breaks)] <- breaks[length(breaks)] + 100000
      res[, `:=`(BIN=cut(pos, breaks=breaks, labels=F, include.lowest=T))]
      summ <- mutate(res, MapQabove9=MapQ > 9) %>%
        group_by(BIN, MapQabove9) %>%
        summarize(NumReads=n(), NumUniqueStarts=length(unique(pos))) %>%
        ungroup()
      summ <- merge(summ, bins) %>% select(-BIN)
      l[[s]] <- summ 
    } else {
      l[[s]] <- res
    }
    l[[s]] <- l[[s]] %>% mutate(Sample=s)
  }
  bind_rows(l)
}

# returns 5 column data.table
# with columns chr start end NumGC NumAT
# chr is chr1, ..., chrY
GetSequenceStats <- function(chr="chr22", bin.size=1000000){
  stats <- fread(paste0("gunzip -c ", sequence.stats.dir, "/", chr, ".", bin.size, ".stats.txt.gz"))
  stats <- stats[, `:=`(NumGC=`7_num_C`+`8_num_G`, NumAT=`6_num_A`+`9_num_T`)]
  stats <- stats[, c(1,2,3,9,10)]
  setnames(stats, c("chr", "start", "end", "NumGC", "NumAT"))
  stats[, `:=`(chr=paste0("chr", chr))]
  stats
}

# returns sequence stats for that bin size for entire genome
GetAllSequenceStats <- function(bin.size=1000000){
  l <- list()
  for (i in all.chrs){
    l[[i]] <- GetSequenceStats(chr=i, bin.size=bin.size)
  }
  bind_rows(l)
}

GetSampleOfWindows <- function(bin.size=100000, 
                               num.windows=30000){
  stats <- GetAllSequenceStats(bin.size)
  stats[, `:=`(NumBases=NumGC+NumAT)]
  stats <- stats[NumBases==bin.size, ]
  set.seed(1)
  if (nrow(stats) < num.windows){
    return(stats)
  } else {
    sampled <- sample(1:nrow(stats), num.windows)
    return(stats[sampled])
  }
}

# bins should have columns chr, start, end
# returns data.frame with column
# BIN MapQabove9 NumReads NumUniqueStarts chr start end (and add'l columns in bins, like NumGC NumAT)
GetBinnedReadStarts <- function(sample.id="2M-3", bins){
  chrs <- unique(bins$chr)
  l <- list()
  read.starts <- GetReadStarts(sample.id)
  for (i in chrs){
    chr.bins <- bins[chr==i, ]
    chr.read.starts <- read.starts[chr==i, ]
    breaks <- sort(unique(c(chr.bins$start, chr.bins$end)))
    breaks[1] <- breaks[1]-100
    breaks[length(breaks)] <- breaks[length(breaks)]+100
    chr.read.starts[, `:=`(BIN=cut(pos, breaks=breaks, labels=F))]
    chr.summ <- mutate(chr.read.starts, MapQabove9=MapQ>9) %>%
      group_by(BIN, MapQabove9) %>%
      summarize(NumReads=n(), NumUniqueStarts=length(unique(pos))) %>%
      ungroup()
    chr.bins <- arrange(chr.bins, start) %>% mutate(BIN=1:n())
    l[[i]] <- merge(chr.summ, chr.bins)
  }
  bind_rows(l)
}

gzipped <- function(a){paste0(a, ".gz")}

LoadHighMapQ <- function(bin.size){
  f <- paste0("/Users/db2175/NonDData/Data/CarmiBAMs/PGDlympho_data/PGDlympho_data/BinnedReadStarts/", bin.size, "highmapq.bins.txt")
  if (!file.exists(f)){
    ll <- LoadBinnedReadStarts(bin.size=bin.size, high.mapq=T)
    write.table(ll, f, row.names=F)
  }
  res <- read.table(f, header=T) %>% 
    rename(PctGC=PercGC, ReadCount=NumReads) %>% select(-NumGC, -NumAT)
}

LoadBinnedReadStarts <- function(samples=FamilyInfoWithSex()$Sample, bin.size=1000000, high.mapq=T){
  prefix <- paste0("/Users/db2175/NonDData/Data/CarmiBAMs/PGDlympho_data/PGDlympho_data/BinnedReadStarts/", bin.size, "/")
  l <- list()
  for (s in samples){
    f <- paste0(prefix, s, bin.size, ".binned.txt")
  r <- read.table(f, header=T) %>% 
      mutate(Sample=s) %>% 
      filter(MapQabove9) %>% 
      select(-MapQabove9)
    l[[s]] <- r
  }
  all <- bind_rows(l) %>% mutate(BinID=paste(chr, BIN))
  bin.dat <- select(all, chr, BIN, BinID, start, end, NumGC, NumAT) %>% unique()
  dat1 <- select(all, BinID, NumReads, Sample) %>% spread(Sample, NumReads, fill=0) %>% gather(Sample, NumReads, -BinID)
  dat2 <- select(all, BinID, NumUniqueStarts, Sample) %>% spread(Sample, NumUniqueStarts, fill=0) %>% gather(Sample, NumUniqueStarts, -BinID)
  merged <- merge(dat1, dat2)
  merged <- merge(merged, bin.dat)
  merged$PercGC <- merged$NumGC / (merged$NumGC + merged$NumAT)
  merged$NumBases <- merged$NumGC + merged$NumAT
  merged
}

SaveAllBinnedReadStarts <- function(samples=FamilyInfoWithSex()$Sample, bin.size=1000000){
  prefix <- paste0("/Users/db2175/NonDData/Data/CarmiBAMs/PGDlympho_data/PGDlympho_data/BinnedReadStarts/", bin.size, "/")
  bins <- GetAllSequenceStats(bin.size=bin.size)
  l <- list()
  for (s in samples){
    print(s)
    f <- paste0(prefix, s, bin.size, ".binned.txt")
    r <- GetBinnedReadStarts(sample.id=s, bins=bins)
    write.table(r, file=f, row.names=F, quote=F)
  }
}

FamilyInfoWithSex <- function(){
  csv <- read.csv("FamilyInfoWithSex.csv")
  csv$chr.prefix <- T
  csv$chr.prefix[csv$Family=="Ethiopian"] <- F
  csv
}

# returns three column data table with columns Sample, Family and Status
# Status is either Normal (good sequencing quality, no large CNVs)
#                  Failed (terrible sequencing quality) or
#                  Presumed aneuploid (poor sequencing quality or large CNVs)
FamilyInfo <- function(){
  failed <- c("R1_46L_3blk", "R1_75L_3", "R1_71L_1", "R2_126L_27", "R1_79L_7")
  
  ethiopian <- c("33866A", "33876A", "33882A", "33889A")
  
  presumed_normal <- c("R2_106L_1", "R2_108L_5", "R2_109L_9", 
                       "R2_110L_10", 
                       "R2_118L_12", "R2_120L_17", "R2_122L_22", 
                       "R2_125L_26", 
                       "R2_2M_3", "R1_55L_6", "R1_56L_11", 
                       "R1_57L_13", "R1_58L_15", 
                       "R1_68L_11", "R1_76L_4", 
                       "R1_XX_135F", "R2_XX_135F", 
                       "R1_XY_98F", "R2_XY_98F")
  
  csv1 <- fread("/Users/db2175/NonDData/Data/CarmiBAMs/PGDlympho_data/PGDlympho_data/PGSblind_BAMs/transloseq-blindrun1.csv", skip="Sample_ID", header=T)
  csv2 <- fread("/Users/db2175/NonDData/Data/CarmiBAMs/PGDlympho_data/PGDlympho_data/PGSblind_BAMs/pipeline_sample_sheet_PGS_blindRun2.csv", skip="Sample_ID", header=T)
  csv <- rbind(csv1, csv2)[!Sample_ID=="", .(Sample_ID, Sample_Plate)] %>%
    unique(by=NULL)
  setnames(csv, c("Sample", "Family"))
  csv <- rbind(csv, data.table(Sample=ethiopian, Family="Ethiopian"))
  normal_samples <- gsub("_", "-", c(ethiopian, substring(presumed_normal, 4)), fixed=T)
  failed_samples <- gsub("_", "-", substring(failed, 4), fixed=T)
  # not all Ethiopian samples should be low-noise
  # maybe normal should be low-noise or something? 
  csv$Status <- "Presumed aneuploid"
  csv$Status[csv$Sample %in% normal_samples] <- "Normal"
  csv$Status[csv$Sample %in% failed_samples] <- "Failed"
  csv$Status[csv$Sample=="33882A"] <- "Presumed aneuploid"
  browser()
  csv
}

InferSex <- function(){
  fi <- FamilyInfo()
  rd <- LoadHighMapQ(bin.size)
  gc <- GCPredictions(rd)
  gc <- mutate(gc, Ratio=ReadCount / Predicted)
  gc <- merge(gc, fi)
  summary <- filter(gc, !Status=="Failed") %>% 
    group_by(Sample, Status, chr) %>% 
    summarize(MedianRatio=median(Ratio))
  sex.chrs <- filter(summary, chr %in% c("chrX", "chrY")) %>% spread(chr, MedianRatio)
  ggplot(sex.chrs, aes(x=chrX, y=chrY)) + geom_point()
  fi$Sex <- ifelse(fi$Sample %in% (filter(sex.chrs, chrY < 0.1))$Sample, "Female", "Male")
  fi$Sex[fi$Status=="Failed"] <- NA
  write.csv(fi, "FamilyInfoWithSex.csv", row.names=F)
}


# GetAllReadDepth1000000 <- function(){
#   fread("/Users/db2175/NonDData/Data/CarmiBAMs/PGDlympho_data/PGDlympho_data/AllReadDepth1000000.txt") %>% as.data.frame()
# }
# 
# AllReadDepth1000000 <- function(){
#   rd <- ReferenceSamplesReadDepth()
#   cd <- CombineReadDepth()
#   all <- bind_rows(rd, cd)
#   write.table(all, "/Users/db2175/NonDData/Data/CarmiBAMs/PGDlympho_data/PGDlympho_data/AllReadDepth1000000.txt", quote=T, sep="\t", row.names=F)
# }
# 
# 
# 

#d <- "/Users/db2175/NonDData/Data/CarmiBAMs/PGDlympho_data/PGDlympho_data/PGSblind_BAMs/"

# runs.dirs <- data.frame(
#   RunName=c("R1", "R2", "R1", "R2"), 
#   FileList=c(paste0(d, "references_run1/samples.txt"), 
#              paste0(d, "references_run2/samples.txt"), 
#              paste0(d, "samples_run1/samples.txt"), 
#              paste0(d, "samples_run2/samples.txt")), 
#   RunCode=c("PGSRef1", "PGSRef2", "PGSSamples1", "PGSSamples2"))




# GetGCCorrect <- function(){
#   f <- "/Users/db2175/NonDData/Data/CarmiBAMs/PGDlympho_data/PGDlympho_data/AllReadDepth1000000.GC.Rdata"
#   if (!file.exists(f)){
#     norm <- GCCorrect()
#     save(norm, file=f)
#   } else {
#     load(f)
#   }
#   norm
# }

# # chr start end PctGC Sample ReadCount Family Normal Failed
# CombineReadDepth <- function(){
#   rd <- list()
#   for (chr in 1:22){
#     r <- read.table(paste0(read.depth.dir, "/chr", chr, ".txt"))
#     colnames(r) <- c("chr", "start", "end", "33866", 
#                      "33876", "33882", "33889")
#     colnames(r) <- c("chr", "start", "end", paste0("Child", 1:4))
#     s <- read.table(paste0(sequence.stats.dir, "/chr", chr, ".stats.txt"), comment.char="", header=T)
#     s <- s[, c(1:3, 5, 10)]
#     colnames(s) <- c("chr", "start", "end", "PctGC", "NumN")
#     m <- merge(r, s)
#     m$start <- m$start - 1
#     rd[[chr]] <- m
#   }
#   rd <- bind_rows(rd) %>% 
#     filter(NumN==0) %>% 
#     select(-NumN) %>% 
#     gather(Sample, ReadCount, -chr, -start, -end, -PctGC)
#   rd$Family <- "Eth"
#   rd$Normal <- rd$Sample %in% c("Child1", "Child2", "Child4")
#   rd$Failed <- F
#   rd
# }
# 
# # chr start end PctGC Sample ReadCount Family Normal Failed
# ReferenceSamplesReadDepth <- function(){
#   rd <- list()
#   for (i in 1:nrow(runs.dirs)){
#     sample.names <- read.table(runs.dirs$FileList[i], header=F)$V1
#     c.sample.names <- paste0(runs.dirs$RunName[i], "_", sample.names)
#     c.sample.names <- gsub("-", "_", c.sample.names)
#     sample.dict <- sample.names
#     names(sample.dict) <- c.sample.names
#     sg <- runs.dirs$RunCode[i]
#     for (chr in 1:22){
#       r <- read.table(paste0(read.depth.dir, "/", sg, "chr", 
#                              chr, ".1000000.txt"))
#       colnames(r) <- c("chr", "start", "end", c.sample.names)
#       r$chr <- gsub("chr", "", r$chr) %>% as.numeric()
#       s <- read.table(
#         paste0(sequence.stats.dir, "/chr", chr, ".stats.txt"), 
#         comment.char="", header=T)
#       s <- s[, c(1:3, 5, 10)]
#       colnames(s) <- c("chr", "start", "end", "PctGC", "NumN")
#       s$start <- s$start - 1
#       m <- merge(r, s) %>% gather(Sample, ReadCount, -chr, -start, -end, -PctGC, -NumN)
#       m$Sample <- sample.dict[m$Sample]
#       rd[[paste0(sg, chr)]] <- m
#     }
#   }
#   all <- bind_rows(rd) %>% 
#     filter(NumN==0) %>% select(-NumN)
#   csv <- FamilyInfo()
#   all <- merge(all, csv) %>% select(-Sample)
#   all$Normal <- all$Sample %in% presumed_normal
#   all$Failed <- all$Sample %in% failed
#   all
# }
