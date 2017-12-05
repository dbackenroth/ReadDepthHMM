PlotHMMResults <- function(l){
  l <- select(l, chr, start, end, State, NormGC, WindowMedian, Sample)
  
  l <- l %>% mutate(gen.start=GetGenomicCoordinates(chr, start), 
                    gen.end=GetGenomicCoordinates(chr, end))
  midcoords <- GetMidChromosomalCoordinates()
  cnvs.combined <- group_by(l, chr) %>%
    select(chr, gen.start, gen.end, State) %>%
    arrange(gen.start) %>%
    mutate(DiffState=c(0, !State[2:n()]==State[1:(n()-1)])) %>%
    mutate(Segment=cumsum(DiffState)) %>%
    group_by(chr, Segment) %>%
    summarize(seg.start=min(gen.start), 
              seg.end=max(gen.end), 
              State=State[1], 
              NumStates=length(unique(State))) %>%
    filter(!State=="TwoCopies")
  print(as.data.frame(cnvs.combined))
  p <- ggplot(NULL) + 
    geom_point(data=l, 
               aes(x=gen.start, y=NormGC/WindowMedian, group=chr, col=State), size=0.25) +
    theme_bw() + 
    geom_rect(data=cnvs.combined, aes(xmin=seg.start, xmax=seg.end, ymin=0, ymax=2, fill=State), alpha=0.25) + 
    geom_hline(yintercept=1, alpha=0.1) + 
    scale_y_continuous(breaks=c(0.5, 1, 1.5), limit=c(0, 2)) + 
    ylab("RC/median") + 
    scale_x_continuous(breaks=midcoords, labels=c(1:22, "X", "Y")) + 
    facet_wrap(~Sample, ncol=1) + 
    scale_fill_manual("", guide=F, values=c(ZeroCopies="yellow", 
                                            OneCopy="red", 
                                            ThreeCopies="blue", 
                                            FourCopies="purple")) + 
    scale_color_manual("", guide=F, values=c(ZeroCopies="yellow", 
                                             OneCopy="red", 
                                             TwoCopies="black",
                                             ThreeCopies="blue", 
                                             FourCopies="purple")) + 
    xlab("") + 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())
}

PlotCNVWithRibbon <- function(chr, cnv.start, cnv.end, 
                   plot.start, plot.end, 
                   cnv.sample, ref.samples){
  bins <- data.frame(start=seq(plot.start, plot.end, by=5000))
  bins$end <- bins$start + 5000
  cnv.dat <- GetRegionWithTabixSamples(samples=c(ref.samples, cnv.sample), 
                                       chr=chr, start=plot.start, end=plot.end, 
                                       bins=bins)
  # add in structural zeros
  cnv.spread <- select(cnv.dat, -NumUniqueStarts, -end) %>% spread(Sample, NumReads, fill=0) %>% gather(Sample, NumReads, -MapQabove9, -start)
  cnv.summ <- group_by(cnv.spread, Sample) %>%
    mutate(ReadPct=NumReads/sum(NumReads)) %>% filter(MapQabove9)
  cnv.summ <- group_by(cnv.summ, MapQabove9, start) %>% 
    summarize(RefMedian=median(ReadPct[Sample %in% ref.samples]), 
              Ref90=quantile(ReadPct[Sample %in% ref.samples], 0.9), 
              Ref10=quantile(ReadPct[Sample %in% ref.samples], 0.1), 
              CNVSample=ReadPct[Sample==cnv.sample])
  cnv.summ <- mutate(cnv.summ, RatioSample=CNVSample/RefMedian, 
                     Ratio10=Ref10/RefMedian, 
                     Ratio90=Ref90/RefMedian)
  p <- ggplot(filter(cnv.summ, MapQabove9), aes(x=start, y=CNVSample)) + geom_point()  + geom_ribbon(aes(x=start, ymin=Ref10, ymax=Ref90), alpha=0.3) + geom_vline(xintercept=c(cnv.start, cnv.end)) + scale_y_sqrt()
}


PlotOneCNV <- function(cnvs.combined=NULL, sample=NULL, dat=NULL, margin=10000000){
  #save(cnvs.combined, sample, dat, margin, file="temp.Rdata")
  #load("temp.Rdata")
  margin <- 10000000
  #dat <- filter(dat, Sample %in% paste0("Child", 1:4))
  cnv.dat <- filter(dat, chr==cnvs.combined$chr, start >= cnvs.combined$start - margin, 
                    end <= cnvs.combined$end + margin) 
  cnv.dat$SampleColor <- ifelse(cnv.dat$Sample==sample, "Sample", ifelse(cnv.dat$Sample %in% paste0("Child", 1:4), "EthiopianFamily", "Other"))
  alpha.dict <- c(Sample=1, EthiopianFamily=0.9, Other=0.1)
  cnv.dat$Alpha <- alpha.dict[cnv.dat$SampleColor]
  r <- range(cnv.dat$NormGC/cnv.dat$MedianNormals)
  p <- ggplot(NULL) + 
    geom_line(data=cnv.dat, aes(x=start, y=NormGC/MedianNormals, color=SampleColor, group=Sample, alpha=I(Alpha))) + 
    theme_bw() + 
    geom_rect(data=cnvs.combined, aes(xmin=start-500000, xmax=end-500000, ymin=r[1], ymax=r[2]), alpha=0.25) + 
    ggtitle(paste0(sample, " chr", cnvs.combined$chr))
  p
}