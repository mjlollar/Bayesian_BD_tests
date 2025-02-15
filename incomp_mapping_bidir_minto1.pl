#!/usr/bin/perl -w
use strict;
use Getopt::Std;
use List::Util qw/shuffle/;

#Get command line arguments
our($opt_f, $opt_g, $opt_n, $opt_r, $opt_e, $opt_a, $opt_b, $opt_c, $opt_d, $opt_s, $opt_t, $opt_u, $opt_v, $opt_w, $opt_x, $opt_y, $opt_z);
getopt('abcdefgnrstuvwxyz');

unless (defined($opt_f) && defined($opt_g) && defined($opt_e) && defined($opt_n) && defined($opt_r) && defined($opt_a) && defined($opt_b) && defined($opt_c) && defined($opt_d) && defined($opt_s) && defined($opt_t) && defined($opt_u) && defined($opt_v) && defined($opt_w) && defined($opt_x) && defined($opt_y) && defined($opt_z)){
    print <<EOF;

This is a forward simulation program designed to test power for epistasis mapping.

It simulates without causative loci, and returns genotypes for males for each window.

It requires the following arguments:

-g [number of generations of crosses before sampling (g>=2)]
-r [number of accepted replicates to simulate]
-f [Output file name prefix, don't include extension]

-a [window position (0-based, genome-wide with X first) of the left boundary of the interval where focal genotype involves ancestry '0']
-c [window position (0-based, genome-wide with X first) of the right boundary of the interval where focal genotype involves ancestry '0']
-b [window position (0-based, genome-wide with X first) of the left boundary of the interval where focal genotype involves ancestry '2']
-d [window position (0-based, genome-wide with X first) of the right boundary of the interval where focal genotype involves ancestry '2']

-e [tolerance - the proportion a simulated count can be off from the empircal count - suggest 0.1]

-n [number of males to initially simulate - should be much larger than number actually used - suggest 2000]

-s [# of steriles that are focal at locus A and focal at locus B]
-t [# of steriles that are focal at locus A and non-focal at locus B]
-u [# of steriles that are non-focal at locus A and focal at locus B]
-v [# of steriles that are non-focal at locus A and non-focal at locus B]
-w [# of fertiles that are focal at locus A and focal at locus B]
-x [# of fertiles that are focal at locus A and non-focal at locus B]
-y [# of fertiles that are non-focal at locus A and focal at locus B]
-z [# of fertiles that are non-focal at locus A and non-focal at locus B]

EOF
exit 1;
}

#Get command line arguments
#my $NLoci = "";
my $gens = $opt_g;
#my $SelectionProp = "";
#my $BatchNum = "";
#my $NCrosses = 1;
my $LastNf = $opt_n;
my $replicates = $opt_r;
my $FileStem = $opt_f;
my $LocusAPosLeft = $opt_a;
my $LocusAPosRight = $opt_c;
my $LocusBPosLeft = $opt_b;
my $LocusBPosRight = $opt_d;
my $LocusAPos = 0;
my $LocusBPos = 0;
my $LocusARange = ($LocusAPosRight - $LocusAPosLeft) + 1;
my $LocusBRange = ($LocusBPosRight - $LocusBPosLeft) + 1;
my $tolerance = $opt_e;
#my $FileStem = $opt_o;
my $EmpSterilesFF = $opt_s;
my $EmpSterilesFN = $opt_t;
my $EmpSterilesNF = $opt_u;
my $EmpSterilesNN = $opt_v;
my $EmpFertilesFF = $opt_w;
my $EmpFertilesFN = $opt_x;
my $EmpFertilesNF = $opt_y;
my $EmpFertilesNN = $opt_z;

my $EmpSterileCount = $EmpSterilesFF + $EmpSterilesFN + $EmpSterilesNF + $EmpSterilesNN;
my $EmpFertileCount = $EmpFertilesFF + $EmpFertilesFN + $EmpFertilesNF + $EmpFertilesNN;

my $EmpSterilesFFUpper = $EmpSterilesFF * (1 + $tolerance);
if ($EmpSterilesFFUpper < ($EmpSterilesFF + 1.01)){
  $EmpSterilesFFUpper = $EmpSterilesFF + 1.01;
}
my $EmpSterilesFFLower = $EmpSterilesFF * (1 - $tolerance);
if ($EmpSterilesFFLower > ($EmpSterilesFF - 1.01)){
  $EmpSterilesFFLower = $EmpSterilesFF - 1.01;
}
my $EmpSterilesFNUpper = $EmpSterilesFN * (1 + $tolerance);
if ($EmpSterilesFNUpper < ($EmpSterilesFN + 1.01)){
  $EmpSterilesFNUpper = $EmpSterilesFN + 1.01;
}
my $EmpSterilesFNLower = $EmpSterilesFN * (1 - $tolerance);
if ($EmpSterilesFNLower > ($EmpSterilesFN - 1.01)){
  $EmpSterilesFNLower = $EmpSterilesFN - 1.01;
}
my $EmpSterilesNFUpper = $EmpSterilesNF * (1 + $tolerance);
if ($EmpSterilesNFUpper < ($EmpSterilesNF + 1.01)){
  $EmpSterilesNFUpper = $EmpSterilesNF + 1.01;
}
my $EmpSterilesNFLower = $EmpSterilesNF * (1 - $tolerance);
if ($EmpSterilesNFLower > ($EmpSterilesNF - 1.01)){
  $EmpSterilesNFLower = $EmpSterilesNF - 1.01;
}
my $EmpSterilesNNUpper = $EmpSterilesNN * (1 + $tolerance);
if ($EmpSterilesNNUpper < ($EmpSterilesNN + 1.01)){
  $EmpSterilesNNUpper = $EmpSterilesNN + 1.01;
}
my $EmpSterilesNNLower = $EmpSterilesNN * (1 - $tolerance);
if ($EmpSterilesNNLower > ($EmpSterilesNN - 1.01)){
  $EmpSterilesNNLower = $EmpSterilesNN - 1.01;
}
my $EmpFertilesFFUpper = $EmpFertilesFF * (1 + $tolerance);
if ($EmpFertilesFFUpper < ($EmpFertilesFF + 1.01)){
  $EmpFertilesFFUpper = $EmpFertilesFF + 1.01;
}
my $EmpFertilesFFLower = $EmpFertilesFF * (1 - $tolerance);
if ($EmpFertilesFFLower > ($EmpFertilesFF - 1.01)){
  $EmpFertilesFFLower = $EmpFertilesFF - 1.01;
}
my $EmpFertilesFNUpper = $EmpFertilesFN * (1 + $tolerance);
if ($EmpFertilesFNUpper < ($EmpFertilesFN + 1.01)){
  $EmpFertilesFNUpper = $EmpFertilesFN + 1.01;
}
my $EmpFertilesFNLower = $EmpFertilesFN * (1 - $tolerance);
if ($EmpFertilesFNLower > ($EmpFertilesFN - 1.01)){
  $EmpFertilesFNLower = $EmpFertilesFN - 1.01;
}
my $EmpFertilesNFUpper = $EmpFertilesNF * (1 + $tolerance);
if ($EmpFertilesNFUpper < ($EmpFertilesNF + 1.01)){
  $EmpFertilesNFUpper = $EmpFertilesNF + 1.01;
}
my $EmpFertilesNFLower = $EmpFertilesNF * (1 - $tolerance);
if ($EmpFertilesNFLower > ($EmpFertilesNF - 1.01)){
  $EmpFertilesNFLower = $EmpFertilesNF - 1.01;
}
my $EmpFertilesNNUpper = $EmpFertilesNN * (1 + $tolerance);
if ($EmpFertilesNNUpper < ($EmpFertilesNN + 1.01)){
  $EmpFertilesNNUpper = $EmpFertilesNN + 1.01;
}
my $EmpFertilesNNLower = $EmpFertilesNN * (1 - $tolerance);
if ($EmpFertilesNNLower > ($EmpFertilesNN - 1.01)){
  $EmpFertilesNNLower = $EmpFertilesNN - 1.01;
}

my $F2N = $LastNf * 2;

my $ScanWidth = 50;

my $WinFile = 'power_test_windows_4.txt';
#my $OutputFile = $FileStem . '_nullsims' . $BatchNum . '.txt';
my $OutputFile = '';
#open O, ">$OutputFile";
#print O "ProbSterileIfFocal\tProbSterileIfNonFocal\tLocusAPos\tLocusBPos\n";

my $EnvVarSD = 1;  #Environmental variance:  value of 1 is analogous to 20% of the total genetic effort
my $F1N = $F2N;  #Not currently varying the number of F1's independently
my $LastN = $LastNf * 2;
#my $replicates = $TrueReps * $NCrosses; #how many replicates to actually simulate
my $multiple = 1;  #combine signals from this many independent crosses
#my $OutputFile = 'SibSamNull_F' . $gens . '_N' . $F2N . '_LN' . $LastN . '_S' . $SelectionProp . '_C' . $NCrosses . '_R'  . $TrueReps . '.txt';

my $SmoothEachSide = 4;  #include this many windows on each side of the focal window in a weighted moving average (e.g. if = 4, weights are 1,2,3,4,5,4,3,2,1)
my $PrimThresh = 0.1;

#Set up locus positions and effects (current code assumes all loci are of equal magnitude and are distributed evenly across genome
my @LocusPosX = ();
my @LocusPosA = ();
my @LocusStrengthX = ();
my @LocusStrengthA = ();
my $i = 0;
my $pos = 0;
my $depth = 0;

#Get cM and depth information from window files
my $CMsX = 0;
my $CMsA = 0;
my @line = ();
my @DepthsXHigh = ();
my @DepthsAHigh = ();
my @DepthsXLow = ();
my @DepthsALow = ();
my @CMStopsX = ();
my @CMStopsA = ();
my @AutoCMBreaks = ();
my @WinChrsX = ();
my @WinChrsA = ();
my @WinChrs = ();

open W, "<$WinFile" or die "can not open $WinFile\n";;
scalar (<W>);
while (<W>){
  chomp;
  last if m/^$/;
  @line = split;
  if ($line[0] =~ m/X/){
    push @WinChrsX, $line[0];
    push @CMStopsX, $line[3];
    $depth = int($line[4]);
    push @DepthsXHigh, $depth;
    $depth = int($line[5]);
    push @DepthsXLow, $depth;
  }
  else{
    push @WinChrsA, $line[0];
    push @CMStopsA, $line[3];
    $depth = int($line[4]);
    push @DepthsAHigh, $depth;
    $depth = int($line[5]);
    push @DepthsALow, $depth;
  }
}
close W;
$CMsX = $CMStopsX[-1];
@WinChrs = (@WinChrsX, @WinChrsA);
for ($i = 0; $i < @WinChrs; $i++){
  $WinChrs[$i] =~ s/R//;
  $WinChrs[$i] =~ s/L//;
}

#Look for chromosome breaks in autosomal window data (cases where this window's CM stop position is smaller than the last window's)
for ($i = 0; $i < (@CMStopsA - 1); $i++){
  if ($CMStopsA[$i+1] < $CMStopsA[$i]){
    $CMStopsA[$i] += $CMsA;
    $CMsA = $CMStopsA[$i];
    push @AutoCMBreaks, $CMsA;
  }
  else{
    $CMStopsA[$i] += $CMsA;
  }
}
$CMStopsA[-1] += $CMsA;
$CMsA = $CMStopsA[-1];
$i = @AutoCMBreaks + 1;
print "Found window data for X chromosome ($CMsX cMs) and $i autosomes ($CMsA total cMs)\n";

#set up arrays
my @PositionsX = ();
my @AncDiffsX = ();
my @PositionsA = ();
my @AncDiffsA = ();
for ($i = 0; $i < @CMStopsX; $i++){
  $pos = $CMStopsX[$i] / $CMsX;
  push @PositionsX, $pos;
  push @AncDiffsX, 0;
}
for ($i = 0; $i < @CMStopsA; $i++){
  $pos = $CMStopsA[$i] / $CMsA;
  push @PositionsA, $pos;
  push @AncDiffsA, 0;
}

#Factorial subroutine
my $A = 0;
my $B = 0;
my $FactIn = 0;
sub factorial{
    $A = 1;
    $B = 1;
    while ($A <= $FactIn) {
	$A *= $B;
	$B++;
    }
    return $A;
}

#Probability of a given number of recombination events (Poisson) for each chromosome
#(allowing up to 10 recombination events per chromosome per generation)
#Cumulative probabilities are given for each number, for easy use with random numbers later.
my $TotalLength = $CMsX;
my $Pr = 0;
my $PrSum = 0;
my @RecombProbsX = ();
my $ExpRecEvents = 0.01 * $TotalLength;
for ($i = 0; $i < 11; $i++){
  $FactIn = $i;
  $Pr = ((exp(-1 * $ExpRecEvents)) * ($ExpRecEvents ** $i)) / &factorial;
  $PrSum = $PrSum + $Pr;
  push @RecombProbsX, $PrSum;
}
my @RecombProbsA = ();
$TotalLength = $CMsA;
$ExpRecEvents = 0.01 * $TotalLength;
$PrSum = 0;
for ($i = 0; $i < 11; $i++){
  $FactIn = $i;
  $Pr = ((exp(-1 * $ExpRecEvents)) * ($ExpRecEvents ** $i)) / &factorial;
  $PrSum = $PrSum + $Pr;
  push @RecombProbsA, $PrSum;
}

for ($i = 0; $i < @AutoCMBreaks; $i++){
  $AutoCMBreaks[$i] = $AutoCMBreaks[$i] / $CMsA;
}

#build breakpoint matrix for F1 individuals
my @OldFemaleXAoA = ();
my @OldFemaleAAoA = ();
my @OldMaleXAoA = ();
my @OldMaleAAoA = ();
my @NewFemaleXAoA = ();
my @NewFemaleAAoA = ();
my @NewMaleXAoA = ();
my @NewMaleAAoA = ();

my @IndBreaks1 = (); #chromosomes from from 1st parental line
my @IndBreaks2 = (0); #chromosomes from from 2nd parental line

#normal subroutine
sub gaussian_rand {
    my ($u1, $u2);  # uniformly distributed random numbers
    my $w;          # variance, then a weight
    my ($g1, $g2);  # gaussian-distributed numbers

    do {
        $u1 = 2 * rand() - 1;
        $u2 = 2 * rand() - 1;
        $w = $u1*$u1 + $u2*$u2;
    } while ( $w >= 1 );

    $w = sqrt( (-2 * log($w))  / $w );
    $g1 = $u2 * $w;
    return $g1;
}


    
#more declarations
my $g = 0;
my $j = 0;
my $k = 0;
my $r = 0;
my $NEachSex = 0;
my $random = 0;
my $random2 = 0;
my $parents = 0;
my $PAllele1 = 0;
my $PAllele2 = 0;
my $AncBeforeBreak = 0;
my $AncAfterBreak = 0;
my $CurrentAllele = 0;
my $LastBreak = 0;
my $AncProp = 0;
my $EnvVar = 0;
my $anc = 0;
my $allele = 0;
my $break = 0;
my $LowerThresh = 0;
my $UpperThresh = 0;
my $AncDiff = 0;
my $MaxDiff = 0;
my $MaxDiffA = 0;
my $MaxDiffX = 0;
my $MaxPos = 0;
my $RegionStart = 0;
my $RegionStop = 0;
my $PowerCount = 0;
my $SmoothNum = 0;
my $SmoothDenom = 0;
my $weight = 0;
my $PeakNow = 0;
my $ProbSterileIfFocal = 0;
my $ProbSterileIfNonFocal = 0;
my $CountsMatch = 0;
my $SimSterileCount = 0;
my $SimFertileCount = 0;
my $IndividualsAssembled = 0;
my $LocusAStart = 0;
my $LocusAStop = 0;
my $LocusBStart = 0;
my $LocusBStop = 0;
my $SimSterilesFF = 0;
my $SimSterilesFN = 0;
my $SimSterilesNF = 0;
my $SimSterilesNN = 0;
my $SimFertilesFF = 0;
my $SimFertilesFN = 0;
my $SimFertilesNF = 0;
my $SimFertilesNN = 0;
my $attempts = 0;
my @RecombPos = ();
my @AncBreaks = ();
my @A1Breaks = ();
my @A2Breaks = ();
my @ancestries = ();
my @BreakpointNumbers = ();
my @AncPropsX = ();
my @AncPropsA = ();
my @AncestryTracker = ();
my @BreakpointTracker = ();
my @AncXAoA = ();
my @AncAAoA = ();
my @blank = ();
my @phenotypes = ();
my @sorted = ();
my @array = ();
my @AncXHighAoA = ();
my @AncAHighAoA = ();
my @AncXLowAoA = ();
my @AncALowAoA = ();
my @AncPropsXHigh = ();
my @AncPropsXLow = ();
my @AncPropsAHigh = ();
my @AncPropsALow = ();
my @PeakDist = ();
my @PeakDistAoA = ();
my @AllPeakDists = ();
my @MaxDiffs = ();
my @MaxDiffAs = ();
my @MaxDiffXs = ();
my @AncDiffsAoA = ();
my @AllPeaks = ();
my @AllPeaksAoA = ();
my @SmoothDiffs = ();
my @FullIndAoA = ();
my @OutputAoA = ();
my @BDCountAoA = ();

#Replicate loop (not indented)
for ($r = 1; $r <= $replicates; $r++){

@OldFemaleXAoA = ();
@OldFemaleAAoA = ();
@OldMaleXAoA = ();
@OldMaleAAoA = ();

for ($i = 0; $i < ($F1N / 2); $i++){
  push @OldFemaleXAoA, [ @IndBreaks1 ];
  push @OldFemaleXAoA, [ @IndBreaks2 ];
  push @OldFemaleAAoA, [ @IndBreaks1 ];
  push @OldFemaleAAoA, [ @IndBreaks2 ];
  push @OldMaleAAoA, [ @IndBreaks1 ];
  push @OldMaleAAoA, [ @IndBreaks2 ];
}
for ($i = 0; $i < ($F1N / 4); $i++){
  push @OldMaleXAoA, [ @IndBreaks1 ];
  push @OldMaleXAoA, [ @IndBreaks2 ];
}
$NEachSex = $F2N / 2;

#For each individual in the next generation, choose parents, enact recombination from mother, transfer breakpoints
for ($g = 2; $g <= $gens; $g++){
  $parents = @OldFemaleXAoA / 2;
  if ($g == $gens){
    $NEachSex = $LastN / 2;
  }

#for females in the next generation, choose mother
  for ($i = 0; $i < $NEachSex; $i++){
    $random = int(rand($parents));
    $random2 = rand;
    if ($random2 < 0.5){
      $PAllele1 = $random * 2;
      $PAllele2 = $PAllele1 + 1;   
    }
    else{
      $PAllele2 = $random * 2;
      $PAllele1 = $PAllele2 + 1;   
    }

#for X chromosome, determine the number of recombination events (store that many random breakpoints) from maternal alleles
    $random = rand;
    @RecombPos = ();
    for ($j = 0; $j < @RecombProbsX; $j++){
      if ($random > $RecombProbsX[$j]){
	push @RecombPos, rand;
      }
      else{
	last;
      }
    }
    if (@RecombPos > 1){
      @RecombPos = sort @RecombPos;
    }
#in light of recombination positions, determine new ancestry along recombinant chromosome
    if (@RecombPos == 0){
      @AncBreaks = @{$OldFemaleXAoA[$PAllele1]};
      push @NewFemaleXAoA, [ @AncBreaks ];
    }
    else{
      @AncBreaks = ();
      @A1Breaks = @{$OldFemaleXAoA[$PAllele1]};
      @A2Breaks = @{$OldFemaleXAoA[$PAllele2]};
      $CurrentAllele = 1;
      $LastBreak = -1;
      for ($j = 0; $j < @RecombPos; $j++){
	$AncBeforeBreak = 1;
	$AncAfterBreak = 1;
	if ($CurrentAllele % 2 == 1){
	  for ($k = 0; $k < @A1Breaks; $k++){
	    last if ($A1Breaks[$k] > $RecombPos[$j]);
	    if ($A1Breaks[$k] > $LastBreak){
	      push @AncBreaks, $A1Breaks[$k];
	      $LastBreak = $A1Breaks[$k];
	    }
	    $AncBeforeBreak++;
	  }
	  for ($k = 0; $k < @A2Breaks; $k++){
	    last if ($A2Breaks[$k] > $RecombPos[$j]);
	    $AncAfterBreak++;
	  }
	}
	else{
	  for ($k = 0; $k < @A2Breaks; $k++){
	    last if ($A2Breaks[$k] > $RecombPos[$j]);
	    if ($A2Breaks[$k] > $LastBreak){
	      push @AncBreaks, $A2Breaks[$k];
	      $LastBreak = $A2Breaks[$k];
	    }
	    $AncBeforeBreak++;
	  }
	  for ($k = 0; $k < @A1Breaks; $k++){
	    last if ($A1Breaks[$k] > $RecombPos[$j]);
	    $AncAfterBreak++;
	  }
	}
	if ( ($AncBeforeBreak % 2) != ($AncAfterBreak % 2) ){
	  push @AncBreaks, $RecombPos[$j];
	  $LastBreak = $RecombPos[$j];
	}
	$CurrentAllele++;
      }
      if ($CurrentAllele % 2 == 1){
	for ($k = 0; $k < @A1Breaks; $k++){
	  if ($A1Breaks[$k] > $LastBreak){
	      push @AncBreaks, $A1Breaks[$k];
	    }
	}
      }
      else{
	for ($k = 0; $k < @A2Breaks; $k++){
	  if ($A2Breaks[$k] > $LastBreak){
	      push @AncBreaks, $A2Breaks[$k];
	    }
	}
      }
      push @NewFemaleXAoA, [ @AncBreaks ];
    }

#for the merged autosome, determine the number of recombination events (store that many random breakpoints) from maternal alleles
    $random = rand;
    @RecombPos = ();
    for ($j = 0; $j < @RecombProbsA; $j++){
      if ($random > $RecombProbsA[$j]){
	push @RecombPos, rand;
      }
      else{
	last;
      }
    }
#each time there's a chromosome break, switch to the parent's other copy with 50% probability
    for ($j = 0; $j < @AutoCMBreaks; $j++){
      $random = rand;
      if ($random > 0.5){
	push @RecombPos, $AutoCMBreaks[$j];
      }
    }
    if (@RecombPos > 1){
      @RecombPos = sort @RecombPos;
    }
#in light of recombination positions, determine new ancestry along recombinant chromosome
    if (@RecombPos == 0){
      @AncBreaks = @{$OldFemaleAAoA[$PAllele1]};
      push @NewFemaleAAoA, [ @AncBreaks ];
    }
    else{
      @AncBreaks = ();
      @A1Breaks = @{$OldFemaleAAoA[$PAllele1]};
      @A2Breaks = @{$OldFemaleAAoA[$PAllele2]};
      $CurrentAllele = 1;
      $LastBreak = -1;
      for ($j = 0; $j < @RecombPos; $j++){
	$AncBeforeBreak = 1;
	$AncAfterBreak = 1;
	if ($CurrentAllele % 2 == 1){
	  for ($k = 0; $k < @A1Breaks; $k++){
	    last if ($A1Breaks[$k] > $RecombPos[$j]);
	    if ($A1Breaks[$k] > $LastBreak){
	      push @AncBreaks, $A1Breaks[$k];
	      $LastBreak = $A1Breaks[$k];
	    }
	    $AncBeforeBreak++;
	  }
	  for ($k = 0; $k < @A2Breaks; $k++){
	    last if ($A2Breaks[$k] > $RecombPos[$j]);
	    $AncAfterBreak++;
	  }
	}
	else{
	  for ($k = 0; $k < @A2Breaks; $k++){
	    last if ($A2Breaks[$k] > $RecombPos[$j]);
	    if ($A2Breaks[$k] > $LastBreak){
	      push @AncBreaks, $A2Breaks[$k];
	      $LastBreak = $A2Breaks[$k];
	    }
	    $AncBeforeBreak++;
	  }
	  for ($k = 0; $k < @A1Breaks; $k++){
	    last if ($A1Breaks[$k] > $RecombPos[$j]);
	    $AncAfterBreak++;
	  }
	}
	if ( ($AncBeforeBreak % 2) != ($AncAfterBreak % 2) ){
	  push @AncBreaks, $RecombPos[$j];
	  $LastBreak = $RecombPos[$j];
	}
	$CurrentAllele++;
      }
      if ($CurrentAllele % 2 == 1){
	for ($k = 0; $k < @A1Breaks; $k++){
	  if ($A1Breaks[$k] > $LastBreak){
	      push @AncBreaks, $A1Breaks[$k];
	    }
	}
      }
      else{
	for ($k = 0; $k < @A2Breaks; $k++){
	  if ($A2Breaks[$k] > $LastBreak){
	      push @AncBreaks, $A2Breaks[$k];
	    }
	}
      }
      for ($k = 0; $k < (@AncBreaks-1); $k++){
	if ($AncBreaks[$k] == $AncBreaks[$k+1]){
	  splice (@AncBreaks, $k, 2);
	  $k--;
	}
      }
      push @NewFemaleAAoA, [ @AncBreaks ];
    }

#for females in the next generation, choose father
#transfer paternal chromosomes without recombination (except allow autosomal independent assortment)
    $PAllele1 = int(rand($parents));
    @AncBreaks = @{$OldMaleXAoA[$PAllele1]};
    push @NewFemaleXAoA, [ @AncBreaks ];
    $random = rand;
    if ($random > 0.5){
      $PAllele1 = $PAllele1 * 2;
      $PAllele2 = $PAllele1 + 1;
    }
    else{
      $PAllele2 = $PAllele1 * 2;
      $PAllele1 = $PAllele2 + 1;
    }
#each time there's a chromosome break, switch to the parent's other copy with 50% probability
    @RecombPos = ();
    for ($j = 0; $j < @AutoCMBreaks; $j++){
      $random = rand;
      if ($random > 0.5){
	push @RecombPos, $AutoCMBreaks[$j];
      }
    }
#in light of recombination positions, determine new ancestry along recombinant chromosome
    if (@RecombPos == 0){
      @AncBreaks = @{$OldMaleAAoA[$PAllele1]};
      push @NewFemaleAAoA, [ @AncBreaks ];
    }
    else{
      @AncBreaks = ();
      @A1Breaks = @{$OldMaleAAoA[$PAllele1]};
      @A2Breaks = @{$OldMaleAAoA[$PAllele2]};
      $CurrentAllele = 1;
      $LastBreak = -1;
      for ($j = 0; $j < @RecombPos; $j++){
	$AncBeforeBreak = 1;
	$AncAfterBreak = 1;
	if ($CurrentAllele % 2 == 1){
	  for ($k = 0; $k < @A1Breaks; $k++){
	    last if ($A1Breaks[$k] > $RecombPos[$j]);
	    if ($A1Breaks[$k] > $LastBreak){
	      push @AncBreaks, $A1Breaks[$k];
	      $LastBreak = $A1Breaks[$k];
	    }
	    $AncBeforeBreak++;
	  }
	  for ($k = 0; $k < @A2Breaks; $k++){
	    last if ($A2Breaks[$k] > $RecombPos[$j]);
	    $AncAfterBreak++;
	  }
	}
	else{
	  for ($k = 0; $k < @A2Breaks; $k++){
	    last if ($A2Breaks[$k] > $RecombPos[$j]);
	    if ($A2Breaks[$k] > $LastBreak){
	      push @AncBreaks, $A2Breaks[$k];
	      $LastBreak = $A2Breaks[$k];
	    }
	    $AncBeforeBreak++;
	  }
	  for ($k = 0; $k < @A1Breaks; $k++){
	    last if ($A1Breaks[$k] > $RecombPos[$j]);
	    $AncAfterBreak++;
	  }
	}
	if ( ($AncBeforeBreak % 2) != ($AncAfterBreak % 2) ){
	  push @AncBreaks, $RecombPos[$j];
	  $LastBreak = $RecombPos[$j];
	}
	$CurrentAllele++;
      }
      if ($CurrentAllele % 2 == 1){
	for ($k = 0; $k < @A1Breaks; $k++){
	  if ($A1Breaks[$k] > $LastBreak){
	      push @AncBreaks, $A1Breaks[$k];
	    }
	}
      }
      else{
	for ($k = 0; $k < @A2Breaks; $k++){
	  if ($A2Breaks[$k] > $LastBreak){
	      push @AncBreaks, $A2Breaks[$k];
	    }
	}
      }
      for ($k = 0; $k < (@AncBreaks-1); $k++){
	if ($AncBreaks[$k] == $AncBreaks[$k+1]){
	  splice (@AncBreaks, $k, 2);
	  $k--;
	}
      }
      push @NewFemaleAAoA, [ @AncBreaks ];
    }
  }

#for males in the next generation, choose mother
  for ($i = 0; $i < $NEachSex; $i++){
    $random = int(rand($parents));
    $random2 = rand;
    if ($random2 < 0.5){
      $PAllele1 = $random * 2;
      $PAllele2 = $PAllele1 + 1;   
    }
    else{
      $PAllele2 = $random * 2;
      $PAllele1 = $PAllele2 + 1;   
    }

#for X chromosome, determine the number of recombination events (store that many random breakpoints) from maternal alleles
    $random = rand;
    @RecombPos = ();
    for ($j = 0; $j < @RecombProbsX; $j++){
      if ($random > $RecombProbsX[$j]){
	push @RecombPos, rand;
      }
      else{
	last;
      }
    }
    if (@RecombPos > 1){
      @RecombPos = sort @RecombPos;
    }
#in light of recombination positions, determine new ancestry along recombinant chromosome
    if (@RecombPos == 0){
      @AncBreaks = @{$OldFemaleXAoA[$PAllele1]};
      push @NewMaleXAoA, [ @AncBreaks ];
    }
    else{
      @AncBreaks = ();
      @A1Breaks = @{$OldFemaleXAoA[$PAllele1]};
      @A2Breaks = @{$OldFemaleXAoA[$PAllele2]};
      $CurrentAllele = 1;
      $LastBreak = -1;
      for ($j = 0; $j < @RecombPos; $j++){
	$AncBeforeBreak = 1;
	$AncAfterBreak = 1;
	if ($CurrentAllele % 2 == 1){
	  for ($k = 0; $k < @A1Breaks; $k++){
	    last if ($A1Breaks[$k] > $RecombPos[$j]);
	    if ($A1Breaks[$k] > $LastBreak){
	      push @AncBreaks, $A1Breaks[$k];
	      $LastBreak = $A1Breaks[$k];
	    }
	    $AncBeforeBreak++;
	  }
	  for ($k = 0; $k < @A2Breaks; $k++){
	    last if ($A2Breaks[$k] > $RecombPos[$j]);
	    $AncAfterBreak++;
	  }
	}
	else{
	  for ($k = 0; $k < @A2Breaks; $k++){
	    last if ($A2Breaks[$k] > $RecombPos[$j]);
	    if ($A2Breaks[$k] > $LastBreak){
	      push @AncBreaks, $A2Breaks[$k];
	      $LastBreak = $A2Breaks[$k];
	    }
	    $AncBeforeBreak++;
	  }
	  for ($k = 0; $k < @A1Breaks; $k++){
	    last if ($A1Breaks[$k] > $RecombPos[$j]);
	    $AncAfterBreak++;
	  }
	}
	if ( ($AncBeforeBreak % 2) != ($AncAfterBreak % 2) ){
	  push @AncBreaks, $RecombPos[$j];
	  $LastBreak = $RecombPos[$j];
	}
	$CurrentAllele++;
      }
      if ($CurrentAllele % 2 == 1){
	for ($k = 0; $k < @A1Breaks; $k++){
	  if ($A1Breaks[$k] > $LastBreak){
	      push @AncBreaks, $A1Breaks[$k];
	    }
	}
      }
      else{
	for ($k = 0; $k < @A2Breaks; $k++){
	  if ($A2Breaks[$k] > $LastBreak){
	      push @AncBreaks, $A2Breaks[$k];
	    }
	}
      }
      push @NewMaleXAoA, [ @AncBreaks ];
    }

#for the merged autosome, determine the number of recombination events (store that many random breakpoints) from maternal alleles
    $random = rand;
    @RecombPos = ();
    for ($j = 0; $j < @RecombProbsA; $j++){
      if ($random > $RecombProbsA[$j]){
	push @RecombPos, rand;
      }
      else{
	last;
      }
    }
#each time there's a chromosome break, switch to the parent's other copy with 50% probability
    for ($j = 0; $j < @AutoCMBreaks; $j++){
      $random = rand;
      if ($random > 0.5){
	push @RecombPos, $AutoCMBreaks[$j];
      }
    }
    if (@RecombPos > 1){
      @RecombPos = sort @RecombPos;
    }
#in light of recombination positions, determine new ancestry along recombinant chromosome
    if (@RecombPos == 0){
      @AncBreaks = @{$OldFemaleAAoA[$PAllele1]};
      push @NewMaleAAoA, [ @AncBreaks ];
    }
    else{
      @AncBreaks = ();
      @A1Breaks = @{$OldFemaleAAoA[$PAllele1]};
      @A2Breaks = @{$OldFemaleAAoA[$PAllele2]};
      $CurrentAllele = 1;
      $LastBreak = -1;
      for ($j = 0; $j < @RecombPos; $j++){
	$AncBeforeBreak = 1;
	$AncAfterBreak = 1;
	if ($CurrentAllele % 2 == 1){
	  for ($k = 0; $k < @A1Breaks; $k++){
	    last if ($A1Breaks[$k] > $RecombPos[$j]);
	    if ($A1Breaks[$k] > $LastBreak){
	      push @AncBreaks, $A1Breaks[$k];
	      $LastBreak = $A1Breaks[$k];
	    }
	    $AncBeforeBreak++;
	  }
	  for ($k = 0; $k < @A2Breaks; $k++){
	    last if ($A2Breaks[$k] > $RecombPos[$j]);
	    $AncAfterBreak++;
	  }
	}
	else{
	  for ($k = 0; $k < @A2Breaks; $k++){
	    last if ($A2Breaks[$k] > $RecombPos[$j]);
	    if ($A2Breaks[$k] > $LastBreak){
	      push @AncBreaks, $A2Breaks[$k];
	      $LastBreak = $A2Breaks[$k];
	    }
	    $AncBeforeBreak++;
	  }
	  for ($k = 0; $k < @A1Breaks; $k++){
	    last if ($A1Breaks[$k] > $RecombPos[$j]);
	    $AncAfterBreak++;
	  }
	}
	if ( ($AncBeforeBreak % 2) != ($AncAfterBreak % 2) ){
	  push @AncBreaks, $RecombPos[$j];
	  $LastBreak = $RecombPos[$j];
	}
	$CurrentAllele++;
      }
      if ($CurrentAllele % 2 == 1){
	for ($k = 0; $k < @A1Breaks; $k++){
	  if ($A1Breaks[$k] > $LastBreak){
	      push @AncBreaks, $A1Breaks[$k];
	    }
	}
      }
      else{
	for ($k = 0; $k < @A2Breaks; $k++){
	  if ($A2Breaks[$k] > $LastBreak){
	      push @AncBreaks, $A2Breaks[$k];
	    }
	}
      }
      for ($k = 0; $k < (@AncBreaks-1); $k++){
	if ($AncBreaks[$k] == $AncBreaks[$k+1]){
	  splice (@AncBreaks, $k, 2);
	  $k--;
	}
      }
      push @NewMaleAAoA, [ @AncBreaks ];
    }

#for males in the next generation, choose father
#transfer paternal chromosomes without recombination (except allow autosomal independent assortment)
    $PAllele1 = int(rand($parents));
    $random = rand;
    if ($random > 0.5){
      $PAllele1 = $PAllele1 * 2;
      $PAllele2 = $PAllele1 + 1;
    }
    else{
      $PAllele2 = $PAllele1 * 2;
      $PAllele1 = $PAllele2 + 1;
    }
#each time there's a chromosome break, switch to the parent's other copy with 50% probability
    @RecombPos = ();
    for ($j = 0; $j < @AutoCMBreaks; $j++){
      $random = rand;
      if ($random > 0.5){
	push @RecombPos, $AutoCMBreaks[$j];
      }
    }
#in light of recombination positions, determine new ancestry along recombinant chromosome
    if (@RecombPos == 0){
      @AncBreaks = @{$OldMaleAAoA[$PAllele1]};
      push @NewMaleAAoA, [ @AncBreaks ];
    }
    else{
      @AncBreaks = ();
      @A1Breaks = @{$OldMaleAAoA[$PAllele1]};
      @A2Breaks = @{$OldMaleAAoA[$PAllele2]};
      $CurrentAllele = 1;
      $LastBreak = -1;
      for ($j = 0; $j < @RecombPos; $j++){
	$AncBeforeBreak = 1;
	$AncAfterBreak = 1;
	if ($CurrentAllele % 2 == 1){
	  for ($k = 0; $k < @A1Breaks; $k++){
	    last if ($A1Breaks[$k] > $RecombPos[$j]);
	    if ($A1Breaks[$k] > $LastBreak){
	      push @AncBreaks, $A1Breaks[$k];
	      $LastBreak = $A1Breaks[$k];
	    }
	    $AncBeforeBreak++;
	  }
	  for ($k = 0; $k < @A2Breaks; $k++){
	    last if ($A2Breaks[$k] > $RecombPos[$j]);
	    $AncAfterBreak++;
	  }
	}
	else{
	  for ($k = 0; $k < @A2Breaks; $k++){
	    last if ($A2Breaks[$k] > $RecombPos[$j]);
	    if ($A2Breaks[$k] > $LastBreak){
	      push @AncBreaks, $A2Breaks[$k];
	      $LastBreak = $A2Breaks[$k];
	    }
	    $AncBeforeBreak++;
	  }
	  for ($k = 0; $k < @A1Breaks; $k++){
	    last if ($A1Breaks[$k] > $RecombPos[$j]);
	    $AncAfterBreak++;
	  }
	}
	if ( ($AncBeforeBreak % 2) != ($AncAfterBreak % 2) ){
	  push @AncBreaks, $RecombPos[$j];
	  $LastBreak = $RecombPos[$j];
	}
	$CurrentAllele++;
      }
      if ($CurrentAllele % 2 == 1){
	for ($k = 0; $k < @A1Breaks; $k++){
	  if ($A1Breaks[$k] > $LastBreak){
	      push @AncBreaks, $A1Breaks[$k];
	    }
	}
      }
      else{
	for ($k = 0; $k < @A2Breaks; $k++){
	  if ($A2Breaks[$k] > $LastBreak){
	      push @AncBreaks, $A2Breaks[$k];
	    }
	}
      }
      for ($k = 0; $k < (@AncBreaks-1); $k++){
	if ($AncBreaks[$k] == $AncBreaks[$k+1]){
	  splice (@AncBreaks, $k, 2);
	  $k--;
	}
      }
      push @NewMaleAAoA, [ @AncBreaks ];
    }
  }

#transfer new ancestry breakpoint matrices to old matrices to prepare for next generation
  @OldFemaleXAoA = @NewFemaleXAoA;
  @NewFemaleXAoA = ();
  @OldFemaleAAoA = @NewFemaleAAoA;
  @NewFemaleAAoA = ();
  @OldMaleXAoA = @NewMaleXAoA;
  @NewMaleXAoA = ();
  @OldMaleAAoA = @NewMaleAAoA;
  @NewMaleAAoA = ();
}

#individual ancestry matrix for males for each position
@AncXAoA = ();
@AncAAoA = ();
for ($i = 0; $i < @OldMaleXAoA; $i++){
  push @AncXAoA, [ @blank ];
  push @AncAAoA, [ @blank ];
}

for ($i = 0; $i < @OldMaleXAoA; $i++){
  $anc = 1;
  $break = 0;
  for ($j = 0; $j < @PositionsX; $j++){
    while ( ($break < @{$OldMaleXAoA[$i]} ) && ($OldMaleXAoA[$i][$break] <= $PositionsX[$j])){
      $anc++;
      $break++;
    }
    if ( ($anc % 2) == 1){
      push @{$AncXAoA[$i]}, 1;
    }
    else{
      push @{$AncXAoA[$i]}, 0;
    }
  }
}
for ($i = 0; $i < @OldMaleAAoA; $i++){
  $anc = 1;
  $break = 0;
  for ($j = 0; $j < @PositionsA; $j++){
    while ( ($break < @{$OldMaleAAoA[$i]} ) && ($OldMaleAAoA[$i][$break] <= $PositionsA[$j])){
      $anc++;
      $break++;
    }
    if ( ($anc % 2) == 1){
      push @{$AncAAoA[$i]}, 1;
    }
    else{
      push @{$AncAAoA[$i]}, 0;
    }
  }
}

#convert to output format
@FullIndAoA = ();
for ($i = 0; $i < @AncXAoA; $i++){
  @line = ();
#return X genotypes (0 or 2)
  for ($k = 0; $k < @{$AncXAoA[$i]}; $k++){
    $anc = $AncXAoA[$i][$k] * 2;
    push @line, $anc;
  }
###
#    $pos = @line;
#    print "Output $pos X-linked genotypes and ";
###
  $j = $i * 2;
  for ($k = 0; $k < @{$AncAAoA[$j]}; $k++){
    $anc = $AncAAoA[$j][$k] + $AncAAoA[$j+1][$k];
    push @line, $anc;
  }
###
#    $pos = @line - $pos;
#    print "$pos autosomal genotypes.";
###
  push @FullIndAoA, [ @line ];
}

###
#for ($i = 0; $i < @FullIndAoA; $i++){
#  die if (@{$FullIndAoA[$i]} != 2579);
#}
###

# NEW IN THIS SCRIPT...
#Assemble a group of fertile and sterile males based upon focal 2-locus genotypes and randomly generated probabilities
$CountsMatch = 0;
@OutputAoA = ();
@BDCountAoA = ();
$attempts = 0;
while ($CountsMatch == 0){
  $attempts++;
  @FullIndAoA = shuffle @FullIndAoA;
  @OutputAoA = ();
  @BDCountAoA = ();  ###just added
  $ProbSterileIfFocal = rand;
  $ProbSterileIfNonFocal = rand($ProbSterileIfFocal);
  $IndividualsAssembled = 0;
  $SimFertileCount = 0;
  $SimSterileCount = 0;
  $random = rand($LocusARange);
  $LocusAPos = int($LocusAPosLeft + $random);
  $random = rand($LocusBRange);
  $LocusBPos = int($LocusBPosLeft + $random);
  for ($i = 0; $i < @FullIndAoA; $i++){
    if (($FullIndAoA[$i][$LocusAPos] == 0) && ($FullIndAoA[$i][$LocusBPos] == 2)){  #test if ind has focal 2-locus genotype
      $random = rand;
      if (($random < $ProbSterileIfFocal) && ($SimSterileCount < $EmpSterileCount)){  #test if sterile (and if we need more of those)
	@line = @{$FullIndAoA[$i]};  
	unshift @OutputAoA, [ @line ];
	$SimSterileCount++;
      }
      elsif (($random > $ProbSterileIfFocal) && ($SimFertileCount < $EmpFertileCount)){  #if fertile and we need more of those
	@line = @{$FullIndAoA[$i]};  ###just added
	push @OutputAoA, [ @line ]; 
	$SimFertileCount++;
      }
    }
    else{  #non-focal
      $random = rand;
      if (($random < $ProbSterileIfNonFocal) && ($SimSterileCount < $EmpSterileCount)){  #test if sterile (and if we need more of those)
	@line = @{$FullIndAoA[$i]};  
	unshift @OutputAoA, [ @line ];
	$SimSterileCount++;
      }
      elsif (($random > $ProbSterileIfNonFocal) && ($SimFertileCount < $EmpFertileCount)){  #if fertile and we need more of those
	@line = @{$FullIndAoA[$i]};  ###just added
	push @OutputAoA, [ @line ]; 
	$SimFertileCount++;
      }
    }
    if (($SimSterileCount == $EmpSterileCount) && ($SimFertileCount == $EmpFertileCount)){
      $IndividualsAssembled = 1;
      last;
    }
  }
  next if ($IndividualsAssembled == 0);  

#Scan all window pairs within 50 windows of the focal windows
#to see if any of them match all 8 empirical Breslow-Day counts within tolerance
#Start by identifying the start and stop positions of the scan (in case we can't go 50 windows due to chr breaks)
  for ($LocusAStart = $LocusAPosLeft; $LocusAStart > ($LocusAPosLeft - $ScanWidth); $LocusAStart--){
    last if ($WinChrs[$LocusAStart] ne $WinChrs[$LocusAStart-1]);
    last if ($LocusAStart == 0);
  }
  for ($LocusAStop = $LocusAPosRight; $LocusAStop < ($LocusAPosRight + $ScanWidth); $LocusAStop++){
    last if ($WinChrs[$LocusAStop] ne $WinChrs[$LocusAStop+1]);
    last if ($LocusAStop == (@WinChrs - 1));
  }
  for ($LocusBStart = $LocusBPosLeft; $LocusBStart > ($LocusBPosLeft - $ScanWidth); $LocusBStart--){
    last if ($WinChrs[$LocusBStart] ne $WinChrs[$LocusBStart-1]);
    last if ($LocusAStart == 0);
  }
  for ($LocusBStop = $LocusBPosRight; $LocusBStop < ($LocusBPosRight + $ScanWidth); $LocusBStop++){
    last if ($WinChrs[$LocusBStop] ne $WinChrs[$LocusBStop+1]);
    last if ($LocusBStop == (@WinChrs - 1));
  }      

###
#  $i = @OutputAoA;
#  print "At attempt $attempts, OutputAoA has $i rows. "; 
#  $i = @{$OutputAoA[0]};
#  print "First row has $i columns and ";
#  $i = @{$OutputAoA[-1]};
#  print "last row has $i columns.\n";
#  print "$LocusAStart $LocusAStop $LocusBStart $LocusBStop\n";
#  $i = @{$OutputAoA[60]};
#  print "OutputAoA row 60 has $i columns\n";
###
  
#Check if all 8 simulated Breslow-Day counts are within tolerance of empirical counts
  for ($i = $LocusAStart; $i <= $LocusAStop; $i++){
    for ($j = $LocusBStart; $j <= $LocusBStop; $j++){
      $SimSterilesFF = 0;
      $SimSterilesFN = 0;
      $SimSterilesNF = 0;
      $SimSterilesNN = 0;
      $SimFertilesFF = 0;
      $SimFertilesFN = 0;
      $SimFertilesNF = 0;
      $SimFertilesNN = 0;
      for ($k = 0; $k < @OutputAoA; $k++){
	if ($k < $EmpSterileCount){
	  if ($OutputAoA[$k][$i] == 0){ ### why never uninit here or just below? Is k never below EmpSterileCount?
	    if ($OutputAoA[$k][$j] == 2){ 
	      $SimSterilesFF++;
	    }
	    else{
	      $SimSterilesFN++;
	    }
	  }
	  else{
	    if ($OutputAoA[$k][$j] == 2){
	      $SimSterilesNF++;
	    }
	    else{
	      $SimSterilesNN++;
	    }
	  }
	}
	else{
###
#	  unless (defined($OutputAoA[$k][$i])){
#	    print "OutputAoA undef at $k $i\n";
#	      $i = @OutputAoA;
#  print "At attempt $attempts, OutputAoA has $i rows. "; 
#  $i = @{$OutputAoA[0]};
#  print "First row has $i columns and ";
#  $i = @{$OutputAoA[-1]};
#  print "last row has $i columns.\n";
  #print "$LocusAStart $LocusAStop $LocusBStart $LocusBStop\n";
#  $i = @{$OutputAoA[60]};
#  print "OutputAoA row 60 has $i columns\n";
#	    die;
#	  }
###	  
	  if ($OutputAoA[$k][$i] == 0){ ### uninit value
	    if ($OutputAoA[$k][$j] == 2){ ### uninit value
	      $SimFertilesFF++;
	    }
	    else{
	      $SimFertilesFN++;
	    }
	  }
	  else{
	    if ($OutputAoA[$k][$j] == 2){
	      $SimFertilesNF++;
	    }
	    else{
	      $SimFertilesNN++;
	    }
	  }
	}
      }
      @line = ();
      push @line, $i;
      push @line, $j;
      push @line, $SimSterilesFF;
      push @line, $SimFertilesFF;
      push @line, $SimSterilesNF;
      push @line, $SimFertilesNF;      
      push @line, $SimSterilesFN;
      push @line, $SimFertilesFN;
      push @line, $SimSterilesNN;
      push @line, $SimFertilesNN;
###
#      print "Probs: $ProbSterileIfFocal $ProbSterileIfNonFocal  SimCounts: $SimSterilesFF $SimSterilesFN $SimSterilesNF $SimSterilesNN $SimFertilesFF $SimFertilesFN $SimFertilesNF $SimFertilesNN\n";
###      
      if (($SimSterilesFF > $EmpSterilesFFLower) && ($SimSterilesFF < $EmpSterilesFFUpper) && ($SimSterilesFN > $EmpSterilesFNLower) && ($SimSterilesFN < $EmpSterilesFNUpper) && ($SimSterilesNF > $EmpSterilesNFLower) && ($SimSterilesNF < $EmpSterilesNFUpper) && ($SimSterilesNN > $EmpSterilesNNLower) && ($SimSterilesNN < $EmpSterilesNNUpper) && ($SimFertilesFF > $EmpFertilesFFLower) && ($SimFertilesFF < $EmpFertilesFFUpper) && ($SimFertilesFN > $EmpFertilesFNLower) && ($SimFertilesFN < $EmpFertilesFNUpper) && ($SimFertilesNF > $EmpFertilesNFLower) && ($SimFertilesNF < $EmpFertilesNFUpper) && ($SimFertilesNN > $EmpFertilesNNLower) && ($SimFertilesNN < $EmpFertilesNNUpper)){
	$CountsMatch = 1;
	print "Found a matching data set for replicate $r at attempt $attempts for windows $i and $j\n";
	push @line, 1;
#	last;
      }
      else{
	push @line, 0;
      }
      push @BDCountAoA, [ @line ];
    }
#    last if ($CountsMatch == 1);
  }
  if ($CountsMatch == 1){
    $OutputFile = $FileStem . '_' . $r . '.txt';
    open O, ">$OutputFile";
    print O "$ProbSterileIfFocal\t$ProbSterileIfNonFocal\t$LocusAPos\t$LocusBPos\n";
    for ($i = 0; $i < @BDCountAoA; $i++){
      for ($j = 0; $j < @{$BDCountAoA[$i]}; $j++){
        print O $BDCountAoA[$i][$j];
        if ($j == (@{$BDCountAoA[$i]} - 1)){
          print O "\n";
        }
        else{
          print O "\t";
        }
      }
    }
    close O;
    last;
  }
  elsif (($attempts % 100) == 0){
    print "Tested $attempts data sets while trying for accepted replicate $r\n";
###
#    print "This time, SimSterilesFF $SimSterilesFF SimSterilesFN $SimSterilesFN SimSterilesNF $SimSterilesNF SimSterilesNN $SimSterilesNN SimFertilesFF $SimFertilesFF SimFertilesFN $SimFertilesFN SimFertilesNF $SimFertilesNF SimFertilesNN $SimFertilesNN\n"; 
###
  }
}
}

#output
#$OutputFile = $FileStem . '_' . $r . '.txt';
#open O, ">$OutputFile";
#print O "$ProbSterileIfFocal\t$ProbSterileIfNonFocal\t$LocusAPos\t$LocusBPos\n";
#for ($i = 0; $i < @OutputAoA; $i++){
#  for ($j = 0; $j < @{$OutputAoA[$i]}; $j++){
#    print O $OutputAoA[$i][$j];
#    if ($j == (@{$OutputAoA[$i]} - 1)){
#      print O "\n";
#    }
#    else{
#      print O "\t";
#    }
#  }
#}
#close O;
