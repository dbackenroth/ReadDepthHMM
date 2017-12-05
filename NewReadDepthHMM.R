library(dplyr)
library(ggplot2)
library(gridExtra)
library(mgcv)

source("NewPrepareReadCountData.R")
source("CombinedUtilities.R")
source("CNVPlots.R")

# runs the HMM twice
# the second time, calculates an adjustment factor based on the median
# ratio of read count for inferred normal regions to the expected 
# read count for those regions
RDRunTwice <- function(bin.size=1000000, sample="33882A", override.factor=NA){
  if (!is.na(override.factor)){
    mult.factor <- override.factor
  } else mult.factor <- 1
  f <- "temp.gc.uniques.Rdata"
  #f <- "temp.gc.Rdata"
  load(f)
  browser()
  #rd <- LoadHighMapQ(bin.size)
  #rd$ReadCount <- rd$NumUniqueStarts
  #gc <- GCPredictions(rd)
  #save(gc, file=f)
  #browser()
  gc <- filter(gc, NumBases==1000000) #, !chr%in%c("chrX", "chrY"))   # Sample, ReadCount, NumUniqueStarts, chr, BIN, start, end, PctGC, NumBases, Predicted
  fi <- FamilyInfoWithSex()
  gc <- filter(gc, !Sample %in% fi$Sample[fi$Status=="Failed"]) #, !chr %in% c("chrX", "chrY"))
  normal.samples <- filter(fi, Status=="Normal")$Sample
  norm <- Normalize(rd=gc, family.info=fi)
  norm <- filter(norm, Sample==sample) %>% arrange(chr, start)
  r1 <- RDRunHMM(dat=norm, mult.factor=mult.factor)
  print(table(r1$l$State))
  normal <- r1$l %>% filter(State=="TwoCopies")
  adj.factor <- median(normal$NormGC/(normal$WindowMedian/mult.factor))
  print(adj.factor)
  r2 <- RDRunHMM(dat=norm, mult.factor=adj.factor)
  num.copies <- c(ZeroCopies=0, OneCopy=1, TwoCopies=2, ThreeCopies=3, 
                  FourCopies=4)
  st <- num.copies[r2$l$State]
  r2$segments <- GetSegments(r2$l) %>% mutate(Sample=sample)
  short.segments <- filter(r2$segments, end-start < 5000000)
  source("CNVPlots.R")
  for (i in 1:nrow(short.segments)){
    print(i)
    chr <- short.segments$chr[i]
    start <- short.segments$start[i]
    end <- short.segments$end[i]
    ref.samples <- setdiff(normal.samples, sample)
    PlotCNVWithRibbon(chr=chr, cnv.start=start, cnv.end=end, 
                       plot.start=start-1000000, plot.end=end+1000000, 
                       cnv.sample=sample, ref.samples=ref.samples)
  browser()
    ggsave(paste0("ShortCNVsNOTUNIQUE/", i, ".pdf"))
  }
  #p <- PlotHMMResults(r2$l)
  r2
  #list(hmm.res=r2, plot=p)
}

RDRunHMM <- function(dat, mult.factor){
  dat <- select(dat, chr, start, end, NormGC, WindowMedian, WindowSD)
  stopifnot(is.numeric(mult.factor) & length(mult.factor)==1)
  dat$WindowMedian <- dat$WindowMedian * mult.factor
  dat <- dat %>% 
    mutate(Deviation=abs(NormGC-WindowMedian))
  # probably better without intercept
  sd.mean.fit <- lm(Deviation ~ WindowMedian+0, data=dat)
  #library(nnls)
  
  #sd.mean.fit <- nnlm(Deviation)
  #gaulss.fit1 <- gam(list(NormGC~s(WindowMedian), ~s(WindowMedian)+s(WindowSD)), data=dat, family=gaulss())
  #gaulss.fit2 <- gam(list(NormGC~s(WindowMedian), ~s(WindowMedian)), data=dat, family=gaulss())
  #b <- gam(list(accel~s(times,k=20,bs="ad"),~s(times)),
  #         data=mcycle,family=gaulss())
  #bb1 <- predict(gaulss.fit1, newdata=dat, type="response")
  #bb2 <- predict(gaulss.fit2, newdata=dat, type="response")
  #sum(mydnbinomlog(x=dat$NormGC, mean=dat$WindowMedian, sd=bb1))
  #sd.pred <- predict(sd.mean.fit, newdata=dat)
  #plot(sd.pred, 1/bb[, 2])
  
  #browser()
  emps <- GetEmissionProbabilities(dat, sd.mean.fit)
  tm <- RDTransitionMatrix()
  ip <- RDInitialStateProb()
  l <- list()
  all.emps <- list()
  for (Chr in unique(dat$chr)){
    chr.dat <- filter(emps, chr==Chr) %>% arrange(start)
    chr.emps <- select(chr.dat, Del2, Del1, Normal, Dup1, Dup2) %>% as.matrix()
    vit <- rd.DoViterbi(emission.probs.m=chr.emps, initial.state.prob=ip, transition.m=tm)
    State <- c("ZeroCopies", "OneCopy", "TwoCopies", "ThreeCopies", "FourCopies")[vit$viterbi]
    l[[Chr]] <- cbind(chr.dat, State)
    all.emps[[Chr]] <- 
      select(chr.dat, chr, start, end, Del2, Del1, Normal, Dup1, Dup2)
  }
  list(l=bind_rows(l), emps=bind_rows(all.emps))
}

GetEmissionProbabilities <- function(dat, sd.mean.fit){
  dat <- select(dat, chr, start, end, NormGC, WindowMedian)
  mydnbinomlog <- function(x, mean, sd){
    #dnorm(x, mean=mean, sd=sd, log=T)
    size <- mean^2 - (sd^2 - mean)
    dnbinom(x=x, size=size, mu=mean, log=T)
  }
  CalcLogProbs <- function(x, mean, fit){
    mydnbinomlog(x=round(x), mean=mean, sd=predict(fit, newdata=data.frame(WindowMedian=mean)))
  }
  cond <- dat$NormGC > 0 & dat$WindowMedian==0
  good.rows <- which(!cond)
  
  dat$Normal <- 0
  dat$Del1 <- 0
  dat$Del2 <- 0
  dat$Dup1 <- 0
  dat$Dup2 <- 0
  dat$Normal[good.rows] <- CalcLogProbs(x=dat$NormGC[good.rows], mean=dat$WindowMedian[good.rows], fit=sd.mean.fit)
  dat$Del1[good.rows] <- CalcLogProbs(x=dat$NormGC[good.rows], mean=dat$WindowMedian[good.rows] / 2, fit=sd.mean.fit)
  dat$Del2[good.rows] <- CalcLogProbs(x=dat$NormGC[good.rows], mean=dat$WindowMedian[good.rows]/10, fit=sd.mean.fit)
  dat$Dup1[good.rows] <- CalcLogProbs(x=dat$NormGC[good.rows], mean=dat$WindowMedian[good.rows] * 3 / 2, fit=sd.mean.fit)
  dat$Dup2[good.rows] <- CalcLogProbs(x=dat$NormGC[good.rows], mean=dat$WindowMedian[good.rows] * 2, fit=sd.mean.fit)
  
  dat
}

rd.DoViterbi <- function(emission.probs.m,
                      initial.state.prob, 
                      transition.m){
  NUM.STATES <- length(initial.state.prob)
  N <- nrow(emission.probs.m)
  viterbi.matrix <- matrix(NA, nrow=N, ncol=NUM.STATES)
  viterbi.pointers <- matrix(NA, nrow=N, ncol=NUM.STATES)
  viterbi.matrix[1, ] <- initial.state.prob + emission.probs.m[1, ]
  for (i in 2:N) {
    temp.matrix <- viterbi.matrix[i - 1, ] + 
      transition.m
    viterbi.matrix[i, ] <- apply(temp.matrix, 2, max)
    emission.probs <- c(emission.probs.m[i,])
    dim(emission.probs) <- c(NUM.STATES, 1)
    viterbi.matrix[i, ] <- viterbi.matrix[i, ] + emission.probs
    viterbi.pointers[i, ] <- apply(temp.matrix, 2, which.max)
  }
  viterbi.states = vector(length = N)
  viterbi.states[N] = which.max(viterbi.matrix[N, ])
  max.prob <- max(viterbi.matrix[N, ])
  for (i in (N - 1):1) {
    viterbi.states[i] <- viterbi.pointers[i + 1, viterbi.states[i + 1]]
  }
  return(list(mp=max.prob, viterbi=viterbi.states))
}

RDTransitionMatrix <- function(){
  transition.m <- matrix(c(NA, 1/300, 1/300, 0, 0,
                           1/8000, NA, 1/300, 0, 0, 
                           0, 1/8000, NA, 1/8000, 0, 
                           0, 0, 1/300, NA, 1/8000,
                           0, 0, 1/300, 1/300, NA), 
                         byrow=T, nrow=5)
  for (i in 1:5){
    transition.m[i,i] <- 1 - sum(transition.m[i, ], na.rm=T)
  } 
  log(transition.m)
}

RDInitialStateProb <- function(){
  log(c(1/4000000, 1/2000, 1998/2000, 1/2000, 1/4000000))
}

GetSegments <- function(l){
  l <- l %>% mutate(gen.start=GetGenomicCoordinates(chr, start), 
                    gen.end=GetGenomicCoordinates(chr, end))
  cnvs.combined <- group_by(l, chr) %>%
    select(chr, gen.start, gen.end, State, start, end) %>%
    arrange(gen.start) %>%
    mutate(DiffState=c(0, !State[2:n()]==State[1:(n()-1)])) %>%
    mutate(Segment=cumsum(DiffState)) %>%
    group_by(chr, Segment) %>%
    summarize(seg.start=min(gen.start), 
              start=min(start),
              end=max(end),
              seg.end=max(gen.end), 
              numwindows=n(),
              State=State[1], 
              NumStates=length(unique(State))) %>%
    filter(!State=="TwoCopies") %>% as.data.frame()
}

Normalize <- function(rd, family.info){
  normal.samples <- family.info$Sample[family.info$Status=="Normal"]
  male.samples <- filter(family.info, Sex=="Male")$Sample
  female.samples <- filter(family.info, Sex=="Female")$Sample
  rd$Normal <- rd$Sample %in% normal.samples
  rd$MaleNormal <- rd$Sample %in% 
    intersect(normal.samples, male.samples)
  rd$FemaleNormal <- rd$Sample %in% 
    intersect(normal.samples, female.samples)
  
  autosomes <- paste0("chr", 1:22)
  autosomal.rd <- filter(rd, chr %in% autosomes)
  
  norm.autosomal <- group_by(autosomal.rd, chr, start) %>%
    mutate(NormGC=ReadCount * mean(Predicted[Normal]) / Predicted) %>%
    mutate(WindowMedian=median(NormGC[Normal])) %>%
    mutate(WindowSD=sd(NormGC[Normal])) %>%
    ungroup() %>% select(-Predicted)
  
  male.sex.rd <- filter(rd, chr %in% c("chrX", "chrY"), Sample %in% male.samples)
  norm.male.sex <- group_by(male.sex.rd, chr, start) %>%
    mutate(NormGC=ReadCount * mean(Predicted[MaleNormal]) / Predicted) %>%
    mutate(WindowMedian=median(NormGC[MaleNormal])) %>%
    mutate(WindowSD=sd(NormGC[MaleNormal])) %>%
    ungroup() %>% select(-Predicted)
  
  female.sex.rd <- filter(rd, chr %in% c("chrX", "chrY"), Sample %in% female.samples)
  norm.female.sex <- group_by(female.sex.rd, chr, start) %>%
    mutate(NormGC=ReadCount * mean(Predicted[FemaleNormal]) / Predicted) %>%
    mutate(WindowMedian=median(NormGC[FemaleNormal])) %>%
    mutate(WindowSD=sd(NormGC[FemaleNormal])) %>%
    ungroup() %>% select(-Predicted)
  
  norm <- bind_rows(norm.autosomal, norm.male.sex, norm.female.sex)
  norm
}

# MeanNormals, SDNormals, MedianNormals
GCPredictions <- function(rd){
  l <- list()
  for (sample in unique(rd$Sample)){
    print(sample)
    sub <- filter(rd, Sample==sample)
    mod.dat <- filter(sub, !chr %in% c("chrX", "chrY"), !ReadCount==0, NumBases==max(NumBases))
    mod <- gam(ReadCount ~ s(PctGC), family=Gamma(link="log"), data=mod.dat)
    sub$Predicted <- predict(mod, newdata=data.frame(PctGC=sub$PctGC), type='response')
    l[[sample]] <- sub
  }
  rd <- bind_rows(l)
}