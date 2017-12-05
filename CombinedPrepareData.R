source("CombinedUtilities.R")

dict <- c(F33886="Father", M33873="Mother", Ch19758="Child1", Ch33878="Child2", 
          Ch33879="Child3", Ch33890="Child4")

MI <- function(child){
  l <- list()
  for (chr in paste0("chr", 1:22)){
    l[[chr]] <- LoadGroundTruth(chr)
  }
  gt <- bind_rows(l)
  code <- paste0(gt$Father, ":", gt$Mother, ":", gt[, paste0("Child", child)])
  tt <- table(code)
  incons <- c("AA:AA:AB", "AA:AA:BB", "AA:AB:BB", "AA:BB:AA", "AA:BB:BB", 
              "AB:AA:BB", "AB:BB:AA", "BB:AB:AA", "BB:BB:AA", 
              "AA:AA:BB", "BB:BB:AB", "BB:AA:AA")
  print(mean(code %in% incons))
}

# returns array data in the form of a data.frame
# with columns CHROM POS Sample and Gen_Call (AA, AB or AA)
# (long format)
PrepareArrayData <- function(){
  array.file <- paste0(data.dir, "/snpinfo.txt")
  vcf.file <- paste0(data.dir, "vcfinfo.txt")
  mother <- "M33873"
  father <- "F33886"
  children <- c("Ch19758", "Ch33878", "Ch33879", "Ch33890")
  array <- fread(array.file, header=T)[, c("CHROM", "POS", 
                                           paste0(c(mother, father, children), 
                                                  "BaseCalls")), with=F]
  array[, `:=`(CHROM=paste0("chr", CHROM))]
  array <- gather(array, Sample, Genotype, -CHROM, -POS) %>% 
    mutate(Sample=gsub("BaseCalls", "", Sample)) %>% as.data.table()
  seq.pos.dat <- fread(vcf.file, header=T)[REF %in% c("A", "C", "G", "T"), 
                            .(CHROM, POS, REF)]
  seq.pos.dat <- unique(seq.pos.dat, by=NULL)
  setkey(seq.pos.dat, CHROM, POS)
  setkey(array, CHROM, POS)
  merged <- merge(array, seq.pos.dat)   # to get REF column for samples with no sequencing data
  setkey(merged, CHROM, POS, Sample)
  merged$Gen_Call <- "AB"
  merged$Gen_Call[merged$Genotype %in% c("CC", "TT", "AA", "GG") & 
                    substr(merged$Genotype, 1, 1)==merged$REF] <- "AA" 
  merged$Gen_Call[merged$Genotype %in% c("CC", "TT", "AA", "GG") & 
                    !substr(merged$Genotype, 1, 1)==merged$REF] <- "BB"
  # convert sample names to child names
  merged$Sample <- dict[merged$Sample]
  merged <- as.data.frame(merged[, .(CHROM, POS, Sample, Gen_Call)]) %>% 
    arrange(CHROM, POS)
}

PrepareGroundTruthAllChromosomes <- function(){
  array.dat <- PrepareArrayData()
  for (chr in paste0("chr", 1:22)){
    chr.f <- paste0(data.dir, "/GroundTruth/", chr, ".ground.truth.csv")
    gt <- PrepareGroundTruth(chrom=chr, array.dat=array.dat)
    write.csv(gt, file=chr.f, row.names=F)
  }
}

GroundTruthAllChromosomes <- function(){
  l <- list(22)
  for (i in 1:22){
    l[[i]] <- LoadGroundTruth(paste0("chr", i))
  }
  bind_rows(l)
}

LoadGroundTruth <- function(chr){
  chr.f <- paste0(data.dir, "/GroundTruth/", chr, ".ground.truth.csv")
  read.csv(chr.f)
}

# returns data frame with columns CHROM POS and Child1, Child2, etc...
# wide format
PrepareGroundTruth <- function(chrom="chr1", array.dat=NULL){
  if (is.null(array.dat)) array.dat <- PrepareArrayData()
  ground.truth <- filter(array.dat, CHROM==chrom) %>% 
    select(CHROM, POS, Sample, Gen_Call) %>%
    spread(Sample, Gen_Call) %>%
    arrange(CHROM, POS)
}

# return data 
# observed.data is a data.frame with columns CHROM POS and Father Mother 
# (all samples)
# data for arrays is genotypes (AA, AB or BB)
# for sequencing is Ref/Alt  (e.g. 7/1)

SaveRealData <- function(chrom="chr1"){
  f <- paste0(data.dir, "/RealData/", chrom, ".real.data")
  
  array.dat <- PrepareArrayData()
  seq.dat <- fread(paste0(data.dir, "vcfinfo.txt"), 
                   header=T)[REF %in% c("A", "C", "G", "T"), ]
  #if (!"CHROM" %in% colnames(seq)){
  #  cat("Fix hack for 1KG vcf without chrom column\n")
  #  seq$CHROM <- "chr1"
  #}
  seq.dat <- seq.dat[, .(CHROM, POS, REF, Sample, Ref, Alt)]
  seq.dat <- seq.dat[CHROM==chrom, ]
  seq.dat$Gen_Call <- paste0(seq.dat$Ref, "/", seq.dat$Alt)
  seq.dat <- mutate(seq.dat, Sample=dict[Sample])
  array.dat <- filter(array.dat, CHROM==chrom)
  save(seq.dat, array.dat, file=f)
}

# CHROM POS Sample REF Ref Alt Ground_Truth
CombineDataChildren <- function(){
  cd <- CombineSequencingData() %>% select(-Gen_Call)
  ad <- GroundTruthAllChromosomes() %>% gather(Sample, Ground_Truth, -CHROM, -POS)
  merged <- merge(cd, ad)
  write.table(merged, paste0(data.dir, "/CombinedDataChildren.txt"), row.names=F)
}

GetCombinedDataChildren <- function(){
  fread(paste0(data.dir, "/CombinedDataChildren.txt")) %>% as.data.table()
}

# columns CHROM POS REF Sample Ref Alt Gen_Call
CombineSequencingData <- function(){
  l <- list()
  for (chr in paste0("chr", 1:22)){
    ll <- LoadRealData(chr)
    l[[chr]] <- ll$seq.dat
  }
  all <- bind_rows(l)
}

LoadRealData <- function(chrom="chr1"){
  f <- paste0(data.dir, "/RealData/", chrom, ".real.data")
  load(f)
  return(list(seq.dat=seq.dat, array.dat=array.dat))
}

GetAllRealData <- function(chrom="chr1"){
  rd <- LoadRealData(chrom=chrom)
  array.dat <- rd$array.dat
  seq.dat <- rd$seq.dat %>% select(CHROM, POS, Sample, Gen_Call)
  bind_rows(array.dat, seq.dat)
}

PrepareRealData <- function(subject.data, chrom="chr1"){
  rd <- LoadRealData(chrom=chrom)
  array.dat <- rd$array.dat
  seq.dat <- rd$seq.dat
  sequencing.samples <- subject.data$Sample[subject.data$Type=="Sequencing"]
  array.samples <- subject.data$Sample[subject.data$Type=="Array"]
  combined <- bind_rows(filter(seq.dat, Sample %in% sequencing.samples), 
                      filter(array.dat, Sample %in% array.samples))
  merged <- merge(combined, subject.data) %>%
    rename(Dat=Gen_Call) %>%
    select(CHROM, POS, Sample, Dat) %>%
    spread(Sample, Dat) %>%
    arrange(CHROM, POS)
  observed.data <- merged[, c("CHROM", "POS", subject.data$Sample)]
}

SubjectDatas <- function(){
  # one sibling with SNP, three with sequencing
  subject.data <- data.frame(Sample=c("Father", "Mother", 
                                      "Child1", "Child2", "Child3", "Child4"), 
                             Type=c(rep("Array", 3), rep("Sequencing", 3)), 
                             stringsAsFactors=F)
  # four siblings with sequencing
  subject.data <- data.frame(Sample=c("Father", "Mother", "Child1", "Child2", 
                                      "Child3", "Child4"), 
                             Type=c(rep("Array", 2), rep("Sequencing", 4)), 
                             stringsAsFactors=F)
  # all siblings with array data
  subject.data <- data.frame(Sample=c("Father", "Mother", "Child1", "Child2", 
                                      "Child3", "Child4"), 
                             Type=rep("Array", 6), 
                             stringsAsFactors=F)
}
