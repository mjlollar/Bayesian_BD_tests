#!/usr/bin/env Rscript
library("DescTools")
#### Running as bash script
# Run: Rscript bd_get_pvalues.r <your_input_file> <outname>
args = commandArgs(trailingOnly=TRUE)
if (length(args) != 2){
  stop("Zero arguments provided (two expected)", call.=FALSE)
}

#df_bd <- read.csv("2_aa_rep14_1.txt", sep='\t', skip=1,header=FALSE)
#df_probs <- read.csv("2_aa_rep14_1.txt", sep='\t', nrows=1,header=FALSE)
#out_file <- file("passedrep_2_aa_14.txt")

out_name <- paste("accepted_rep_", as.character(args[4]), ".txt", sep='')
out_file <- file(out_name)
df_bd <- read.csv(args[1], sep='\t', skip=1,header=FALSE)
df_probs <- read.csv(args[1], sep='\t', nrows=1,header=FALSE)
prob_s_if_ster = as.character(df_probs[1,1])
prob_s_if_fert = as.character(df_probs[1,2])

##Add BD pvalue column
bd_tester <- function(x, output){
  ## Cell values with Fisher adjustment to prevent divide by zero occurrences
  bd_1 <- x[3] + 0.5
  bd_2 <- x[5] + 0.5
  bd_3 <- x[4] + 0.5
  bd_4 <- x[6] + 0.5
  bd_5 <- x[7] + 0.5
  bd_6 <- x[9] + 0.5
  bd_7 <- x[8] + 0.5
  bd_8 <- x[10] + 0.5
  ## Calculate Odds 
  odds_one <- ((bd_1) / (bd_1 + bd_3))
  odds_two <- ((bd_2) / (bd_2 + bd_4))
  odds_three <- ((bd_5) / (bd_5 + bd_7))
  odds_four <- ((bd_6) / (bd_6 + bd_8))
  max_odds <- max(c(odds_one, odds_two, odds_three, odds_four))
  if (max_odds %in% c(odds_two, odds_three, odds_four) == TRUE){
    #Skip calculation if max odds are not odds 1
    y <- 999
    return(y)
  } else {
    ## Calculate Breslow Day test
    bd_table <- xtabs(freq ~ ., cbind(expand.grid(phenotype=c("sterile", "fertile"),
                                                  window_one=c("focal", "non-focal"),
                                                  window_two=c("focal", "non-focal")),
                                      freq=c(bd_1, bd_2, bd_3, bd_4, bd_5, bd_6, bd_7, bd_8)))
    y <- BreslowDayTest(bd_table)$p.value
    return(y)
  }
}
pvalues <- apply(df_bd,1,bd_tester)
df_bd <- cbind(df_bd,pvalue = pvalues)

##Get lowest pvalues
lowest_p <- min(df_bd[,"pvalue"])
lowest_mins <- df_bd$V11[df_bd$pvalue==lowest_p]
lowest_winAs <- df_bd$V1[df_bd$pvalue==lowest_p]
lowest_winAs <- sort(unique(lowest_winAs))
winAL_start <- lowest_winAs[1]
winAR_start <- lowest_winAs[length(lowest_winAs)]
lowest_winBs <- df_bd$V2[df_bd$pvalue==lowest_p]
lowest_winBs <- sort(unique(lowest_winBs))
winBL_start <- lowest_winBs[1]
winBR_start <- lowest_winBs[length(lowest_winBs)]

##Add to replicate file if lowest bd has thresh=1
if (1 %in% lowest_mins){
  ##Get intervals

  far_left <- numeric(0)
  far_right <- numeric(0)
  far_up <- numeric(0)
  far_down <- numeric(0)

  ##Left
  for (windowL in lowest_winBs){
    bool_checkL <- TRUE
    sanity_breakL <- 1
    l = winAL_start
    while (isTRUE(bool_checkL)){
      if (df_bd$V11[df_bd$V1==l & df_bd$V2==windowL]==1){
        l = l - 1
        sanity_breakL = sanity_breakL + 1
        if (sanity_breakL==2000000){
          break
        }
      } else {
        far_left[length(far_left)+1] <- as.character(l)
        bool_checkL <- FALSE
        }
    }
  }
  ##Right
  for (windowR in lowest_winBs){
    bool_checkR <- TRUE
    sanity_breakR <- 1
    r = winAR_start
    while (isTRUE(bool_checkR)){
      if (df_bd$V11[df_bd$V1==r & df_bd$V2==windowR]==1){
        r = r + 1
        sanity_breakR = sanity_breakR + 1
        if (sanity_breakR==2000000){
          break
        }
      } else {
        far_right[length(far_right)+1] <- as.character(r)
        bool_checkR <- FALSE
        }
    }
  }
  ##Up
  for (windowU in lowest_winAs){
    bool_checkU <- TRUE
    sanity_breakU <- 1
    u = winBL_start
    while (isTRUE(bool_checkU)){
      if (df_bd$V11[df_bd$V1==windowU & df_bd$V2==u]==1){
        u = u - 1
        sanity_breakU = sanity_breakU + 1
        if (sanity_breakU==2000000){
          break
        }
      } else {
        far_up[length(far_up)+1] <- as.character(u)
        bool_checkU <- FALSE
        }
    }
  }
  ##Down
  for (windowD in lowest_winAs){
    bool_checkD <- TRUE
    sanity_breakD <- 1
    d = winBR_start
    while (isTRUE(bool_checkD)){
      if (df_bd$V11[df_bd$V1==windowD & df_bd$V2==d]==1){
        d = d + 1
        sanity_breakD = sanity_breakD + 1
        if (sanity_breakD==2000000){
          break
        }
      } else {
        far_down[length(far_down)+1] <- as.character(d)
        bool_checkD <- FALSE
        }
    }
  }
  ##Get maximum CI values
  winAL <- as.character(min(far_left))
  winAR <- as.character(max(far_right))
  winBL <- as.character(min(far_up))
  winBR <- as.character(max(far_down))
  #Create header if file does not exist or append to new file
  #if (file.exists(out_file)){
  writeLines(c(prob_s_if_ster,prob_s_if_fert,winAL,winAR,winBL,winBR),out_file)
  #} else{
    #writeLines(c("probS", "probF", "winAL", "winAR", "winBL", "winBR"),file=out_file,sep="\n")
    #writeLines(c(prob_s_if_ster,prob_s_if_fert,winAL,winAR,winBL,winBR),file=out_file,append=TRUE,sep="\n")
  #}
} else{
  print('Replicate did not pass lowest BD pvalue == threshold 1')
}
