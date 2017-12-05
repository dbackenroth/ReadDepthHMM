library(R.utils)
library(dplyr)
library(tidyr)
library(stringi)
options(stringsAsFactors=F)

# add emissions probabilities over windows
# fix recombination probabilities over windows
# fix genotype estimation over windows

array.genotypes <- c("AA", "AB", "BB")

source("CombinedUtilities.R")
source("CombinedPrepareData.R")
source("PrepareReadCountData.R")
source("../HMMCode/ReadDepthHMM.R")
source("TM.R")
#source("Simulate.R")

del.prob <- (1 / 30000) / 10
out.of.del.prob <- 1 / 3000
dup.prob <- (1 / 30000) / 10
out.of.dup.prob <- 1 / 3000

childrenDF <- data.frame(Name=c("Child1", "Child4"), 
                         Reference=c(T,F), 
                         Array=c(T,F), 
                         InferCNVs=c(F,T))



# children.df has columns Name, Reference, Array, InferCNVs
#      Reference, Array and InferCNVs are booleans
CombinedHMMRealData <- function(children.df=childrenDF, 
                        chr="chr3", 
                        forward.backward=F){
  all.data <- GetAllRealData(chr=chr)     # CHROM POS Sample Gen_Call (AA, AB, BB or Ref/Alt, e.g., 5/4)
  children <- children.df$Name
  array.children <- filter(children.df, Array)
  sequencing.children <- filter(children.df, !Array)
  if (length(sequencing.children$Name) > 1) stop("Only one child with sequencing allowed\n")
  cat("Combine with code for RD alone")
  rd.res <- RDRunTwice(sample=sequencing.children$Name, override.factor=NA)
  # returns list with two components
  #     emps is the emission probabilities, has columns chr, start, end, Sample, Del2, Del1, Normal, Dup1, Dup2
  chr.numeric <- as.numeric(gsub("chr", "", chr))
  rd.eps <- filter(rd.res$emps, chr==chr.numeric)
  # keep the right children, label data as sequencing or array
  all.data <- filter(all.data, Sample %in% c("Mother", "Father", children)) %>%
    mutate(Type=ifelse(Gen_Call %in% array.genotypes, "Array", "Sequencing"))
  # get rid of sequencing data for Child1
  all.data <- filter(all.data, 
                     !(Sample %in% array.children & Type=="Sequencing"))
  # get rid of array data for other children
  all.data <- filter(all.data, 
                     !(Sample %in% sequencing.children & Type=="Array"))
  
  
  array.positions <- filter(all.data, Sample=="Father")$POS
  all.data <- filter(all.data, POS %in% array.positions)
  # observed.data has columns CHROM POS Child1 Child3 Father Mother
  # sample columns are either AB, BB, AA or Ref/Alt
  observed.data <- all.data %>% select(-Type) %>% spread(Sample, Gen_Call) %>% arrange(POS)
  
  hmm.res <- DoHMM(observed.data, 
                   forward.backward=forward.backward, 
                   children.df=children.df, rd.eps=rd.eps)
  hmm.res
}

# observations should be a data frame with columns
#   CHROM, POS, Father, Mother, Child1, ChildX
#   CHROM      POS Child1 ChildX Father Mother
#   chr21 15206086     AB    0/0     AB     AB
#   chr21 15291929     AB    1/4     AA     AB
#   1 is the number of reference alleles, 4 is the number of alternate alleles
#   there should be data for only one chromosome in observations

DoHMM <- function(observations, forward.backward=F, 
                  children.df, rd.eps=NULL){
  # recomb is a list of two vectors of length 1 less than rd.eps$start or observations$POD
  # the names are male.recomb.probs and female.recomb.probs
  if (is.null(rd.eps)) positions <- observations$POS else positions <- rd.eps$start
  recomb <- util.GetRecombinationProbabilities(observations$CHROM[1], positions)
  # get the probabilities of each genotype state given the data
  # CHROM POS Sample GenData Type N A B AA AB BB AA AAB ABB BBB AAAA AAAB AABB ABBB BBBB
  obs.prob <- util.GenerateDataProbabilities(observations, make.long=T)
  prob.data <- MergeProbData(obs.prob, cnv.children=children.df$Name[children.df$InferCNVs], 
                             samples=c("Mother", "Father", children.df$Name))
  hmm.results <- RunHMM(probabilities=prob.data, 
                        male.recomb.probs=recomb$male.recomb.probs, 
                        female.recomb.probs=recomb$female.recomb.probs, 
                        children.df=children.df, 
                        forward.backward=forward.backward, rd.eps=rd.eps)
  
  
  gen <- hmm.results$best.genotype
  rr <- data.frame(BIN=1:length(hmm.results$best.state), 
                   STATE=hmm.results$best.state)
  p <- ggplot(rr, aes(x=BIN, y=as.factor(STATE))) + geom_point() + theme_bw() + xlab("Bin number") + ylab("State")
  browser()
  print(p)
  PredState <- hmm.results$best.state
  res <- list(observations=observations, 
              predictions=cbind(observations[, c("CHROM", "POS")], gen, PredState))
  if (forward.backward){
    marginal.probs <- hmm.results$forward.m + hmm.results$backward.m
    marginal.probs <- t(apply(marginal.probs, 1, function(x){x-SumProbabilities(x)}))
    marginal.probs <- exp(marginal.probs)
    mp <- as.data.frame(marginal.probs)
    cns <- colnames(mp)
    for (child in 2:n.children){
      child.state <- unlist(lapply(as.list(cns), function(x){
        unlist(strsplit(x, ";"))[child - 1]}))
      pats <- which(child.state %in% c("id", "hp"))
      mp[, paste0("Child", child, "SamePat")] <- apply(mp[, pats], 1, sum)
      mats <- which(child.state %in% c("id", "hm"))
      mp[, paste0("Child", child, "SameMat")] <- apply(mp[, mats], 1, sum)
    }
    res$marginal <- cbind(observations[, c("CHROM", "POS")], mp)
  }
  res
}

GetSNPEmissionProbs <- function(children.df, probabilities){
  # returns a data frame with columns State Child1 Child3 n CondProb Father Mother
  emissions.prob.dict <- 
    MakeEmissionsProbabilityDictionary(children.df=children.df)
  # there are 36 states
  states <- unique(emissions.prob.dict$State)
  n.states <- length(states)
  
  n.variants <- nrow(probabilities)
  emissions.probs <- matrix(NA, nrow=n.variants, ncol=n.states)
  colnames(emissions.probs) <- states
  
  n.children <- length(children.df$Name)
  
  combinations <- unique(emissions.prob.dict$GenCode)
  n.combinations <- length(combinations)
  comb.probabilities <- matrix(NA, nrow=n.variants, ncol=n.combinations)
  colnames(comb.probabilities) <- combinations
  # for each combination multiply the probabilities of all the
  # component genotypes
  # now we know for each combination of genotypes that probability of the
  # data given that genotype
  # comb.probabilities has 1 row for each SNP and one column for each combination
  # of genotypes for parents and both children
  for (comb in combinations){
    probs.to.multiply <- unlist(strsplit(comb, ":", fixed=T))
    ptm <- probabilities[, probs.to.multiply]
    comb.probabilities[, comb] <- apply(ptm, 1, prod)
  }
  # comb.probabilities for each SNP (in row) gives probability of that genotype, 
  # taking into account all samples
  emissions.prob.dict <- as.data.table(emissions.prob.dict)
  setkey(emissions.prob.dict, State)
  # for each state, add probabilities of all genotypes consistent with that state
  for (i in 1:n.states){
    state <- states[i]
    state.emissions.prob.dict <- emissions.prob.dict[State==state & 
                                                       CondProb>0, ]
    state.probs <- rep(0, n.variants)
    for (j in 1:nrow(state.emissions.prob.dict)){
      genotype.code <- state.emissions.prob.dict$GenCode[j]
      genotype.probs <- comb.probabilities[, genotype.code]
      state.probs <- state.probs + genotype.probs * state.emissions.prob.dict$CondProb[j]
    }
    emissions.probs[, i] <- log(state.probs)
  }
  list(emissions.probs=emissions.probs, comb.probabilities=comb.probabilities)
}

CombineEmissions <- function(positions, emissions.probs, rd.eps, children.df){
  states <- colnames(emissions.probs)
  n.children <- length(children.df$Name)
  
  # add up emissions.probs for each bin
  ep.df <- data.frame(emissions.probs)
  ep.df$POS <- positions
  boundaries <- c(rd.eps$start, rd.eps$end[nrow(rd.eps)])
  ep.df$BIN <- cut(ep.df$POS, breaks=boundaries, labels=F)
  ep.df.save <- ep.df
  ep.for.later <- filter(ep.df, !is.na(BIN))
  ep.df <- filter(ep.df, !is.na(BIN)) %>% select(-POS)
  ep.binned <- group_by(ep.df, BIN) %>% summarize_all(funs(sum)) %>% arrange(BIN)
  colnames(ep.binned) <- gsub("X", "", colnames(ep.binned), fixed=T)
  colnames(ep.binned) <- gsub(".", ":", colnames(ep.binned), fixed=T)
  colnames(ep.for.later) <- gsub("X", "", colnames(ep.for.later), fixed=T)
  colnames(ep.for.later) <- gsub(".", ":", colnames(ep.for.later), fixed=T)
  cat("Need to fix this when have more than one sample with CNVs\n")
  rd.eps$BIN <- 1:nrow(rd.eps)
  state.rdstate.dict <- data.frame(STATE=colnames(emissions.probs), 
                                   RD=NA)
  ddd <- c(0, 1,  1, 2, 2, 2)
  for (i in 1:nrow(state.rdstate.dict)){
    statenums <- as.numeric(unlist(strsplit(state.rdstate.dict$STATE[i], ":", fixed=T)))
    num.copies <- sum(ddd[statenums+1])
    state.rdstate.dict$RD[i] <- c("Del2", "Del1", "Normal", "Dup1", "Dup2")[num.copies+1]
  }
  # to fill in missing rows in ep.binned
  ep.fill <- matrix(0, ncol=ncol(emissions.probs), nrow=nrow(rd.eps))
  colnames(ep.fill) <- colnames(emissions.probs)
  ep.fill <- as.data.frame(ep.fill) %>% mutate(BIN=1:nrow(rd.eps)) %>% filter(!BIN %in% ep.binned$BIN)
  ep.binned <- bind_rows(ep.binned, ep.fill) %>% arrange(BIN)
  rd.eps <- arrange(rd.eps, BIN)
  # need to fill in missing rows in ep.binned
  for (i in 1:nrow(state.rdstate.dict)){
    cc <- state.rdstate.dict$STATE[i]
    ep.binned[, cc] <- ep.binned[, cc] + rd.eps[, state.rdstate.dict$RD[i]]
  }
  combined.emissions.probs <- select(ep.binned, -BIN) %>% as.matrix()
}

# probabilities has columns POS Mother_AA Mother_BB Mother_AB 
#                               Father_AA Father_BB, Father_AB etc...
RunHMM <- function(probabilities, male.recomb.probs, 
                   female.recomb.probs, children.df,
                   forward.backward=F, rd.eps=NULL){
  
  e.list <- GetSNPEmissionProbs(children.df, probabilities)
  comb.probabilities <- e.list$comb.probabilities
  
  combined.emissions.probs <- CombineEmissions(positions=probabilities$POS, 
                                               emissions.probs=e.list$emissions.probs, 
                                               rd.eps=rd.eps, 
                                               children.df=children.df)
  initial.state.prob <- log(rep(1/ncol(combined.emissions.probs), ncol(combined.emissions.probs)))
  cat("Starting Viterbi\n")
  
  viterbi <- DoViterbi(combined.emissions.probs, male.recomb.probs=male.recomb.probs,
                       female.recomb.probs=female.recomb.probs,
                       initial.state.prob=initial.state.prob, 
                       n.children=nrow(children.df))
  browser()
  
  if (forward.backward){
    cat("Starting forward pass\n")
    forward.m <- GetForwardMatrix(emission.probs.m=emissions.probs, 
                                  male.recomb.probs=male.recomb.probs, 
                                  female.recomb.probs=female.recomb.probs, 
                                  initial.state.prob=initial.state.prob, 
                                  n.children=n.children)
    cat("Starting backward pass\n")
    backward.m <- GetBackwardMatrix(emission.probs.m=emissions.probs, 
                                    male.recomb.probs=male.recomb.probs, 
                                    female.recomb.probs=female.recomb.probs, 
                                    initial.state.prob=initial.state.prob, 
                                    n.children=n.children)
    colnames(forward.m) <- colnames(backward.m) <- colnames(emissions.probs)
  } else {
    forward.m <- NULL
    backward.m <- NULL
  }
  state <- states[viterbi$viterbi]
  # get most likely genotypes
  # use comb.probabilities
  # 
  u.states <- unique(state)
  
  state.df <- data.frame(BIN=1:nrow(rd.eps), STATE=state)
  comb.probs <- as.data.frame(e.list$comb.probabilities) %>% mutate(POS=probabilities$POS) %>%
    filter(POS %in% ep.for.later$POS)
  bin.pos <- ep.for.later[, c("POS", "BIN")]
  comb.probs <- merge(merge(comb.probs, bin.pos), state.df)
  best.genotype <- rep(NA, nrow(comb.probs))
  # need to rewrite to take into account binning procedure
  for (s in u.states){
    which.variants <- which(comb.probs$STATE==s)
    state.genotypes <- emissions.prob.dict[State==s, GenCode]
    s.comb.probabilities <- comb.probabilities[which.variants, state.genotypes]
    s.comb.probs.colnames <- colnames(s.comb.probabilities)
    temp.max <- apply(s.comb.probabilities, 1, which.max)
    max.state <- s.comb.probs.colnames[temp.max]
    best.genotype[which.variants] <- max.state
  }
  sample.names <- unlist(strsplit(unlist(strsplit(best.genotype[1], ":")), "_"))[c(T,F)]
  best.genotype <- data.frame(Temp=best.genotype) %>% separate_("Temp", sample.names, sep=":")
  best.genotype <- apply(best.genotype, 2, function(x){unlist(strsplit(x, "_"))[c(F,T)]}) %>% as.data.frame()
  best.genotype$POS <- comb.probs$POS
  return(list(best.genotype=best.genotype, best.state=state, 
              forward.m=forward.m, backward.m=backward.m, 
              ep=ep.df.save))
}

MakeEmissionsProbabilityDictionary <- function(children.df){
  grid <- expand.grid(array.genotypes, array.genotypes, stringsAsFactors=F)
  colnames(grid) <- c("Father", "Mother")
  r <- vector('list', nrow(grid))
  for (i in 1:nrow(grid)){
    r[[i]] <- MakeCombinations(father=grid$Father[i], mother=grid$Mother[i], 
                               children.df=children.df)
  }
  all <- bind_rows(r)
  samples <- c("Father", "Mother", children.df$Name)
  
  # construct a genotype code for emissions.prob.dict
  # of form Father_AA:Mother_AA:Child1_AA:Child3_A
  c.df <- all[, samples] %>% as.data.frame()
  for (s in samples){
    c.df[, s] <- paste(s, c.df[, s], sep="_")
  }
  all$GenCode <- apply(c.df, 1, paste, collapse=":")
  all
}

# given father's and mother's genotype and number of children
# generates table of all possible inheritance states
# and child genotypes
# along with conditional probability (probability of observing
# given child genotypes conditional on parent genotypes and inheritance state)
# 0 is no allele passed (deletion)
# 1 is first haplotype passed
# 2 is second haplotype passed
# 3 is two copies of first haplotype passed
# 4 is two copies of second haplotype passed
# 5 is both haplotypes passed

# the state vector is in the form 01:15 where on each side of the colon
# there are a number of digits equal to the number of non-reference children
# (there should be only one reference child)
# the left hand side is for the father, the right hand side is for the mother
MakeCombinations <- function(father="AB", mother="AB", 
                             children.df=childrenDF){
  # for children with array data we assume diploid genotypes
  # for other children we allow six different states with respect to each parent
  children.names <- children.df$Name
  reference.child.name <- children.df$Name[children.df$Reference]
  seg.list <- list()
  non.ref.children <- filter(children.df, !Reference)
  # for array children, only states 1 and 2 are possible
  # for sequencing children, 6 different states are possible
  # this gets the cartesian product of all possible states
  # across the children for one parent
  for (i in 1:nrow(non.ref.children)){
    if (non.ref.children$Array[i] | !non.ref.children$InferCNVs[i]){
      seg.list[[non.ref.children$Name[i]]] <- c(1,2)
    } else {
      seg.list[[non.ref.children$Name[i]]] <- 0:5
    }
  }
  segregations <- expand.grid(seg.list, stringsAsFactors=F) %>% 
    as.matrix() %>%
    apply(1, function(x){paste(x, collapse="")})
  # now generate all possibilities for father AND mother
  grid <- expand.grid(segregations, segregations, stringsAsFactors=F)
  colnames(grid) <- c("FatherCode", "MotherCode")
  # add observed genotypes
  # first genotype is the one transmitted to the reference child
  # second is the one not transmitted
  # code is in 0:5
  TransmittedGenotype <- function(parental.genotype, code){
    transmitted <- rep(NA, length(code))
    first <- substr(parental.genotype, 1, 1)
    second <- substr(parental.genotype, 2, 2)
    transmitted[code=="0"] <- ""
    transmitted[code=="1"] <- first
    transmitted[code=="2"] <- second
    transmitted[code=="3"] <- paste0(first, first)
    transmitted[code=="4"] <- paste0(second, second)
    transmitted[code=="5"] <- paste0(first, second)
    transmitted
  }
  grid.l <- list()
  k <- 0
  str_rev <- function(x){
    paste0(rev(unlist(strsplit(x, ""))), collapse="")
  }
  # cycle through all combinations of haplotypes transmitted to reference child
  for (mother.rev in c(F,T)){
    for (father.rev in c(F,T)){
      use.grid <- grid
      mother.gen <- ifelse(mother.rev, str_rev(mother), mother)
      father.gen <- ifelse(father.rev, str_rev(father), father)
      # construct reference child genotype
      use.grid[, paste0(reference.child.name, "FromFather")] <- substr(father.gen, 1, 1)
      use.grid[, paste0(reference.child.name, "FromMother")] <- substr(mother.gen, 1, 1)
      
      # construct genotypes for other children
      for (parent in c("Father", "Mother")){
        for (i in 1:nrow(non.ref.children)){
          child.name <- non.ref.children$Name[i]
          code <- substr(use.grid[, paste0(parent, "Code")], i, i)
          parental.genotype <- ifelse(parent=="Father", father.gen, mother.gen)
          use.grid[, paste0(child.name, "From", parent)] <- 
            TransmittedGenotype(parental.genotype, code)
        }
      }
      SortStrings <- function(X){
        sapply(X, function(x){paste0(sort(unlist(strsplit(x, ""))), collapse="")})
      }
      for (child in children.names){
        use.grid[, child] <- 
          SortStrings(paste0(use.grid[, paste0(child, "FromFather")], 
                             use.grid[, paste0(child, "FromMother")]))
      }
      use.grid$State <- paste0(use.grid$FatherCode, ":", use.grid$MotherCode)
      use.grid <- use.grid[, c("State", children.df$Name)]
      k <- k + 1
      grid.l[[k]] <- use.grid
    }
  }
  grid <- bind_rows(grid.l)
  # Columns you want to group by
  grp_cols <- c("State", children.names)
  dots <- lapply(grp_cols, as.symbol)
  u.grid <- group_by_(grid, .dots=dots) %>%
    summarize(n=n()) %>% ungroup() %>% group_by(State) %>%
    mutate(CondProb=n/sum(n))
  u.grid$Father <- father
  u.grid$Mother <- mother
  u.grid
}

DoViterbi <- function(emission.probs.m, male.recomb.probs, female.recomb.probs, 
                      initial.state.prob, n.children){
  states <- colnames(emission.probs.m)
  NUM.STATES <- length(initial.state.prob)
  N <- nrow(emission.probs.m)
  viterbi.matrix <- matrix(NA, nrow=N, ncol=NUM.STATES)
  viterbi.pointers <- matrix(NA, nrow=N, ncol=NUM.STATES)
  viterbi.matrix[1, ] <- initial.state.prob + emission.probs.m[1, ]
  for (i in 2:N) {
    tm <- log(TM(states, male.recomb.probs[i-1], female.recomb.probs[i-1], 
             n.children))
    temp.matrix <- viterbi.matrix[i - 1, ] + 
      tm
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
  print(max.prob)
  return(list(mp=max.prob, viterbi=viterbi.states))
}

# adds two log-space probabilities using the identity
# log (p1 + p2) = log p1 + log(1 + exp(log p2 - log p1))
AddTwoProbabilities <- function(x, y){
  if (is.infinite(x)) return (y)
  if (is.infinite(y)) return (x)
  sum.probs <- max(x, y) + log1p(exp(-abs(x - y)))
}

# adds multiple log-space probabilities
SumProbabilities <- function(x){
  sum.probs <- x[1]
  for (i in 2:length(x)){
    sum.probs <- AddTwoProbabilities(sum.probs, x[i])
  }
  return(sum.probs)
}

# get the forward probabilities
GetForwardMatrix <- function(emission.probs.m, male.recomb.probs, 
                             female.recomb.probs, 
                             initial.state.prob, n.children){
  states <- colnames(emission.probs.m)
  stop("Need to implement transition matrices")
  NUM.STATES <- length(initial.state.prob)
  N <- nrow(emission.probs.m)
  forward.matrix <- matrix(NA, nrow=N, ncol=NUM.STATES)   # matrix to hold forward probabilities
  forward.matrix[1, ] <- initial.state.prob + emission.probs.m[1, ]
  for (i in 2:N){
    # compute matrix with probability we were in state j and are now in state i
    # in temp.matrix[j, i] (ignoring emission of current token)
    temp.matrix <- forward.matrix[i - 1, ] + 
      TransitionM2(father.matrix=tm.matrices$father, 
                   mother.matrix=tm.matrices$mother,
                   male.recomb.prob=male.recomb.probs[i-1], 
                   female.recomb.prob=female.recomb.probs[i-1], 
                   n.children=n.children, 
                   cnv.multiplier=cnv.multiplier.matrix)
    # find the probability that we are in each of the three states
    sum.probs <- apply(temp.matrix, 2, SumProbabilities)
    forward.matrix[i, ] <- sum.probs + emission.probs.m[i, ]
  }  
  return(forward.matrix)  
}

# get the backward probabilities
GetBackwardMatrix <- function(emission.probs.m, male.recomb.probs, 
                              female.recomb.probs, 
                              initial.state.prob, n.children){
  states <- colnames(emission.probs.m)
  stop("Need to implement transition matrices")
  NUM.STATES <- length(initial.state.prob)
  N <- nrow(emission.probs.m)
  
  backward.matrix <- matrix(NA, nrow=N, ncol=NUM.STATES)   # matrix to hold backward probabilities
  
  backward.matrix[N, ] <- rep(0, NUM.STATES)
  for (i in (N - 1):1){
    temp.matrix <- TransitionM2(father.matrix=tm.matrices$father, 
                                mother.matrix=tm.matrices$mother,
                                male.recomb.prob=male.recomb.probs[i], 
                                female.recomb.prob=female.recomb.probs[i], 
                                n.children=n.children, 
                                cnv.multiplier=cnv.multiplier.matrix) + 
      matrix(backward.matrix[i + 1, ], NUM.STATES, NUM.STATES, byrow=T) +
      matrix(emission.probs.m[i+1, ], NUM.STATES, NUM.STATES, byrow=T)
    backward.matrix[i, ] <- apply(temp.matrix, 1, SumProbabilities)
  }  
  final.prob <- backward.matrix[1, ] + emission.probs.m[1, ] + initial.state.prob
  return(backward.matrix)  
}

MergeProbData <- function(obs.prob, cnv.children, samples){
  prob.data <- list()
  for (sample in samples){
    s.data <- obs.prob[Sample==sample, ]
    if (sample %in% cnv.children){
      sample.genotypes <- cnv.genotypes
    } else {
      sample.genotypes <- array.genotypes
    }
    setnames(s.data, sample.genotypes, paste0(sample, "_", sample.genotypes))
    s.data <- as.data.frame(s.data) %>% arrange(POS)
    prob.data[[sample]] <- s.data[, c("POS", paste0(sample, "_", sample.genotypes))]
    N_column <- paste0(sample, "_N")
    if (N_column %in% colnames(prob.data[[sample]])){
      setnames(prob.data[[sample]], N_column, paste0(sample, "_"))
    }
  }
  prob.data <- Reduce(merge, prob.data)
}

DoHMMErrorChecking <- function(observations, children.df){
  # check that all samples have data in observations
  samples <- c("Mother", "Father", children.df$Name)
  if (!all(samples %in% colnames(observations))){
    stop("DoHMM error:All samples", samples, "should be observations data.frame\n")
  }
  # check that only one chromosome has data
  if (length(unique(observations$CHROM))>1){
    stop("DoHMM error: Only one chromosome allowed at a time")
  }
  if (!all(observations$POS==sort(observations$POS))) stop("Position must be sorted")
}
