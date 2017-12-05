MT <- function(){
  samples <- c("33866A", "33876A" , "33882A", "33889A")
  gg <- GetRegionWithTabixSamples(samples, "chrMT", 1, 20000) %>% filter(MapQ>=10)
  m <- matrix(0, nrow=17000, ncol=4)
  for (i in 1:4){
    gg.s <- filter(gg, Sample==samples[i])
    for (j in 1:nrow(gg.s)){
      p <- gg.s$pos[j]
      m[p:(p+50), i] <- m[p:(p+50), i] + 1
    }
  }
  m <- m[1:16571, ]
  colnames(m) <- samples
  pdf("MitochondrialCoverage.pdf")
  plot(1:16571, m[, 1], col="black", ylim=range(m), type='l', xlab="MT Position", 
       ylab="Depth of coverage")
  lines(m[, 2], col="blue")
  lines(m[, 3], col="red")
  lines(m[, 4], col="green")
  dev.off()
  print(colSums(m))
  print(cor(m))
}