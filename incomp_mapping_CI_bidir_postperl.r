Last login: Wed Feb  5 16:49:41 on ttys000

The default interactive shell is now zsh.
To update your account to use zsh, please run `chsh -s /bin/zsh`.
For more details, please visit https://support.apple.com/kb/HT208050.
(base) jeps-MBP-3:~ matt$ ssh mlollar@ap2001.chtc.wisc.edu
Password: 
Password: 
Duo two-factor login for mlollar

Enter a passcode or select one of the following options:

 1. Duo Push to XXX-XXX-6700

Passcode or option (1-1): 1
Success. Logging you in...
Success. Logging you in...
Last login: Fri Jan  3 13:33:41 2025 from 10.134.100.1
_____________________________________________________________________
 #####  #     # #######  #####  Issues?  Email chtc@cs.wisc.edu
#     # #     #    #    #     # Unauthorized use prohibited by:
#       #     #    #    #       WI Statutes: s. 947.0125
#       #######    #    #       U.S. Code: 18 USC 1030
#       #     #    #    #       U.S. Code: 18 USC 2510-2522
#     # #     #    #    #     # U.S. Code: 18 USC 2701-2712
 #####  #     #    #     #####  U.S. Code: 18 USC § 1831
For off campus ssh access use https://www.doit.wisc.edu/network/vpn/
_____________________________________________________________________

         Online office hours are available twice a week:
            Tuesdays, 10:30am - 12pm (Central time)
            Thursdays, 3:00 - 4:30pm (Central time)
        Join via this link: go.wisc.edu/chtc-officehours
     Sign in via this link: go.wisc.edu/chtc-officehours-signin

     IMPORTANT ANNOUNCEMENT REGARDING RECOVERED DATA:
     We will be removing recovered data on Feb 17 to improve performance of our /staging file
     system. Please move (mv) or delete your recovered data (located at /recovery) before this date.
     For full details see: https://chtc.cs.wisc.edu/uw-research-computing/data-recovery-fall2024
Filesystem quota report
Storage             Used (GB)    Limit (GB)    Files (#)    File Cap (#)    Quota (%)
----------------  -----------  ------------  -----------  --------------  -----------
/home/mlollar          416.95          1020         2103               0        40.88
/staging/mlollar            0           100            1            1000            0

[mlollar@ap2001 ~]$ ls
apptainer_r_build.log  ci_2024              old_submit1_files
apptainer_r_maker.sub  fake_fly_replicates  old_submit3_files
[mlollar@ap2001 ~]$ cd ci_2024/
[mlollar@ap2001 ci_2024]$ ls
CI_inputs_test.txt                  incomp_mapping_CI_bidir.r
CI_inputs.txt                       incomp_mapping_CI_bidir_v4.pl
ci_secondstep_inputs.txt            inputs
ci_testing_24.py                    inputs_postperl
errors                              logs
errors_perl                         logs_perl
focal_test.py                       outputs
incomp_mapping_CI_bidir_minto1.pl   power_test_windows_4.txt
incomp_mapping_CI_bidir.pl          subs
incomp_mapping_CI_bidir_postperl.r
[mlollar@ap2001 ci_2024]$ cd subs/
[mlollar@ap2001 subs]$ ls
ci_perl_minto1.sh  ci_perl_test_minto1.sub  ci_postperl_aa8025pen.sh  ci_postperl_aa8025pen.sub  first_step
[mlollar@ap2001 subs]$ vi ci_perl_test_minto1.sub 
[mlollar@ap2001 subs]$ v ci_postperl_aa8025pen.sub
-bash: v: command not found
[mlollar@ap2001 subs]$ vi ci_postperl_aa8025pen.sub
[mlollar@ap2001 subs]$ cd ..
[mlollar@ap2001 ci_2024]$ ls
CI_inputs_test.txt        errors_perl                         incomp_mapping_CI_bidir.r      logs_perl
CI_inputs.txt             focal_test.py                       incomp_mapping_CI_bidir_v4.pl  outputs
ci_secondstep_inputs.txt  incomp_mapping_CI_bidir_minto1.pl   inputs                         power_test_windows_4.txt
ci_testing_24.py          incomp_mapping_CI_bidir.pl          inputs_postperl                subs
errors                    incomp_mapping_CI_bidir_postperl.r  logs
[mlollar@ap2001 ci_2024]$ vi incomp_mapping_CI_bidir_postperl.r 

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

"incomp_mapping_CI_bidir_postperl.r" 162L, 5049B                                                      14,54         Top
