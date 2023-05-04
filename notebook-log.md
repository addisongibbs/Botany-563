## My Data
- I chose data on the complete mitochondrial genome of the freshwater bass family.
- I chose this because my research is in water toxicology, so I do a lot of work with fish in the lab. I've never gotten to work with bass, but I found some data on one of them and that got me interested.
- GenBank had only 4 different species' complete mitochondrial genome, which was a little sad but was okay because I had something I could still work with.
## MUSCLE
- I ran it through MUSCLE to get an aligned sequence. I chose MUSCLE mostly because it was the only one I could get to work with my data and on my generation of MacBook.
- Download MUSCLE:
    `muscle3.8.31_i86darwin64.tar.gz`
- Using MUSCLE:
`(base) addisongibbs@AG botany-563 % software/muscle3.8.31_i86darwin64 -in SEQUENCE.fasta -out SEQ-ALIGNED.fasta`

MUSCLE v3.8.31 by Robert C. Edgar

http://www.drive5.com/muscle
This software is donated to the public domain.
Please cite: Edgar, R.C. Nucleic Acids Res 32(5), 1792-97.

SEQUENCE 4 seqs, max length 16490, avg  length 16486
00:00:00      2 MB(0%)  Iter   1  100.00%  K-mer dist pass 1
00:00:00      2 MB(0%)  Iter   1  100.00%  K-mer dist pass 2
00:00:00     21 MB(0%)  Iter   1   33.33%  Align node       
00:00:16    329 MB(4%)  Iter   1  100.00%  Align node
00:00:23    334 MB(4%)  Iter   1  100.00%  Align node
00:00:23    334 MB(4%)  Iter   1  100.00%  Root alignment
00:00:23    334 MB(4%)  Iter   2  100.00%  Root alignment

00:01:00    345 MB(4%)  Iter   3  100.00%  Refine biparts
00:01:29    345 MB(4%)  Iter   4  100.00%  Refine biparts
00:01:37    345 MB(4%)  Iter   4  100.00%  Refine biparts

This gave me an aligned sequence called SEQ-ALIGNED.fasta.

## Distance and Parsimony
- First, we must install and load packages in terminal to perform these functions.
    install.packages("adegenet", dep=TRUE)
    install.packages("phangorn", dep=TRUE)
    *I did this directly in terminal because when I put it in my rscript code it gave me an error
#### Distance Based: Put this code into a file named myscript.R
    
   library(ape)
library(adegenet)
library(phangorn)
dna <- fasta2DNAbin("SEQ-ALIGNED.fasta")
D <- dist.dna(dna, model="TN93")
tre <- njs(D)
tre <- ladderize(tre)
plot(tre, cex=.6)
title("A simple NJ tree")

    
    
   ` % Rscript myscript.R    `  
    output:
Loading required package: ade4

   /// adegenet 2.1.10 is loaded ////////////

   > overview: '?adegenet'
   > tutorials/doc/questions: 'adegenetWeb()' 
   > bug reports/feature requests: adegenetIssues()



Attaching package: ‘phangorn’

The following object is masked from ‘package:adegenet’:

    AICc


 Converting FASTA alignment into a DNAbin object... 


 Finding the size of a single genome... 


 genome size is: 16,490 nucleotides 

( 237  lines per genome )

 Importing sequences... 
..........
 Forming final object... 

...done.

    
    to save it in the background:
    screen -S distance
        Rscript myscript.R
        control + A + D - detaches screen
    screen -r distance - brings it back

#### Parsimony Based:
- put this code into an rscript file called parsimony.R
   library(ape)
library(adegenet)
library(phangorn)
dna <- fasta2DNAbin(file="SEQ-ALIGNED.fasta")
dna2 <- as.phyDat(dna)
tre.ini <- njs(dist.dna(dna,model="raw"))
parsimony(tre.ini, dna2)
tre.pars <- optim.parsimony(tre.ini, dna2)
plot(tre.pars, cex=0.6)
title("plotted tree")

    
   `% Rscript parsimony.R`
   output:
The following object is masked from ‘package:adegenet’:

    AICc


 Converting FASTA alignment into a DNAbin object... 


 Finding the size of a single genome... 


 genome size is: 16,490 nucleotides 

( 237  lines per genome )

 Importing sequences... 
..........
 Forming final object... 

...done.

[1] 23059
Final p-score 1592 after  0 nni operations

## RAxML-NG
Checking the Version:
- `addisongibbs@AG raxml-ng_v1.1.0_macos_x86_64 % ./raxml-ng -v`

RAxML-NG v. 1.1.0 released on 29.11.2021 by The Exelixis Lab.
Developed by: Alexey M. Kozlov and Alexandros Stamatakis.
Contributors: Diego Darriba, Tomas Flouri, Benoit Morel, Sarah Lutteropp, Ben Bettisworth.
Latest version: https://github.com/amkozlov/raxml-ng
Questions/problems/suggestions? Please visit: https://groups.google.com/forum/#!forum/raxml

System: Intel(R) Core(TM) i5-8210Y CPU @ 1.60GHz, 2 cores, 8 GB RAM
-Ran it one time and obtained this ERROR:
`./raxml-ng --check --msa data/SEQ-ALIGNED.fasta --model GTR+G`

RAxML-NG v. 1.1.0 released on 29.11.2021 by The Exelixis Lab.
Developed by: Alexey M. Kozlov and Alexandros Stamatakis.
Contributors: Diego Darriba, Tomas Flouri, Benoit Morel, Sarah Lutteropp, Ben Bettisworth.
Latest version: https://github.com/amkozlov/raxml-ng
Questions/problems/suggestions? Please visit: https://groups.google.com/forum/#!forum/raxml

System: Intel(R) Core(TM) i5-8210Y CPU @ 1.60GHz, 2 cores, 8 GB RAM

RAxML-NG was called at 27-Apr-2023 13:50:51 as follows:

./raxml-ng --check --msa data/SEQ-ALIGNED.fasta --model GTR+G

Analysis options:
  run mode: Alignment validation
  start tree(s): 
  random seed: 1682621451
  SIMD kernels: AVX2
  parallelization: coarse-grained (auto), PTHREADS (auto)

[00:00:00] Reading alignment from file: data/SEQ-ALIGNED.fasta
[00:00:00] Loaded alignment with 4 taxa and 16497 sites

NOTE: Reduced alignment (with duplicates and gap-only sites/taxa removed) 
NOTE: was saved to: /Users/addisongibbs/Desktop/Botany-563/raxml-ng_v1.1.0_macos_x86_64/data/SEQ-ALIGNED.fasta.raxml.reduced.phy

ERROR: Following taxon name contains invalid characters: OK945936.1 Micropterus salmoides mitochondrion, complete genome
ERROR: Following taxon name contains invalid characters: HQ391897.1 Micropterus floridanus mitochondrion, complete genome
ERROR: Following taxon name contains invalid characters: OL693877.1 Micropterus punctulatus mitochondrion, complete genome
ERROR: Following taxon name contains invalid characters: NC_011361.1 Micropterus dolomieu mitochondrion, complete genome

NOTE: Following symbols are not allowed in taxa names to ensure Newick compatibility:
NOTE: " " (space), ";" (semicolon), ":" (colon), "," (comma), "()" (parentheses), "'" (quote). 
NOTE: Please either correct the names manually, or use the reduced alignment file
NOTE: generated by RAxML-NG (see above).

ERROR: Alignment check failed (see details above)!

-Changed taxa names to MD, MF, MS, and MP. Ran again and this information replaced the ERROR messages in the above output:

Partition 0: noname
Model: GTR+FO+G4m
Alignment sites: 16497
Gaps: 0.07 %
Invariant sites: 90.76 %



Alignment can be successfully read by RAxML-NG.


Execution log saved to: /Users/addisongibbs/Desktop/Botany-563/raxml-ng_v1.1.0_macos_x86_64/data/SEQ-ALIGNED.fasta.raxml.log

Analysis started: 27-Apr-2023 13:53:57 / finished: 27-Apr-2023 13:53:57

Elapsed time: 0.015 seconds


## MrBayes
- Converting to .nexus file:
    records = SeqIO.parse("SEQ-ALIGNED.fasta", "fasta")
    count = SeqIO.write(records, "SEQ.nexus", "nexus")
    print("Converted %i records" % count)

- Setup:
    create text file mbblock with the following text:
    begin mrbayes;
    set autoclose=yes;
    prset brlenspr=unconstrained:exp(10.0);
    prset shapepr=exp(1.0);
    prset tratiopr=beta(1.0,1.0);
    prset statefreqpr=dirichlet(1.0,1.0,1.0,1.0);
    lset nst=2 rates=gamma ngammacat=4;
    mcmcp ngen=10000 samplefreq=10 printfreq=100 nruns=1 nchains=3 savebrlens=yes;
    outgroup Anacystis_nidulans;
    mcmc;
    sumt;   
end;

    -Attach mbblock.txt file to data:
    `cat SEQ-ALIGNED.nex mbblock.txt > SEQ-ALIGNED-MB.nex`

- Running MrBayes:

   ` mb SEQ-ALIGNED-MB.nex`
    
 MrBayes 3.2.7a x86_64

                      (Bayesian Analysis of Phylogeny)

                             (Parallel version)
                         (1 processors available)

              Distributed under the GNU General Public License


               Type "help" or "help <command>" for information
                     on the commands that are available.

                   Type "about" for authorship and general
                       information about the program.



   Executing file "SEQ-ALIGNED-MB.nex"
   UNIX line termination
   Longest line length = 93
   Parsing file
   Expecting NEXUS formatted file
   Reading data block
      Allocated taxon set
      Allocated matrix
      Defining new matrix with 4 taxa and 16497 characters
      Data is Dna
      Missing data coded as ?
      Gaps coded as -
      Data matrix is interleaved
      Taxon 1 -> MD
      Taxon 2 -> MF
      Taxon 3 -> MS
      Taxon 4 -> MP
      Successfully read matrix
      Setting default partition (does not divide up characters)
      Setting model defaults
      Seed (for generating default start values) = 1682625236
      Setting output file names to "SEQ-ALIGNED-MB.nex.run<i>.<p|t>"
   Exiting data block
   Reading mrbayes block
      Setting autoclose to yes
      Setting quitonerror to no
      Setting Brlenspr to Unconstrained:Exponential(10.00)
      Successfully set prior model parameters
      Setting Shapepr to Exponential(1.00)
      Successfully set prior model parameters
      Setting Tratiopr to Beta(1.00,1.00)
      Successfully set prior model parameters
      Setting Statefreqpr to Dirichlet(1.00,1.00,1.00,1.00)
      Successfully set prior model parameters
      Setting Nst to 2
      Setting Rates to Gamma
      Setting Ngammacat to 4
      Successfully set likelihood model parameters
      Setting number of generations to 10000
      Setting sample frequency to 10
      Setting print frequency to 100
      WARNING: Reallocation of zero size attempted. This is probably a bug. Problems may follow.
      WARNING: Reallocation of zero size attempted. This is probably a bug. Problems may follow.
      Setting number of runs to 1
      WARNING: Allocation of zero size attempted. This is probably a bug; problems may follow.
      Setting number of chains to 3
      Setting chain output file names to "SEQ-ALIGNED-MB.nex.<p/t>"
      Successfully set chain parameters
      Setting outgroup to taxon "MD"
      Running Markov chain
      MCMC stamp = 8854554416
      Seed = 21617450
      Swapseed = 1682625236
      Model settings:

         Data not partitioned --
            Datatype  = DNA
            Nucmodel  = 4by4
            Nst       = 2
                        Transition and transversion  rates, expressed
                        as proportions of the rate sum, have a
                        Beta(1.00,1.00) prior
            Covarion  = No
            # States  = 4
                        State frequencies have a Dirichlet prior
                        (1.00,1.00,1.00,1.00)
            Rates     = Gamma
                        The distribution is approximated using 4 categories.
                        Shape parameter is exponentially
                        distributed with parameter (1.00).

      Active parameters: 

         Parameters
         ---------------------
         Tratio              1
         Statefreq           2
         Shape               3
         Ratemultiplier      4
         Topology            5
         Brlens              6
         ---------------------

         1 --  Parameter  = Tratio
               Type       = Transition and transversion rates
               Prior      = Beta(1.00,1.00)

         2 --  Parameter  = Pi
               Type       = Stationary state frequencies
               Prior      = Dirichlet

         3 --  Parameter  = Alpha
               Type       = Shape of scaled gamma distribution of site rates
               Prior      = Exponential(1.00)

         4 --  Parameter  = Ratemultiplier
               Type       = Partition-specific rate multiplier
               Prior      = Fixed(1.0)

         5 --  Parameter  = Tau
               Type       = Topology
               Prior      = All topologies equally probable a priori
               Subparam.  = V

         6 --  Parameter  = V
               Type       = Branch lengths
               Prior      = Unconstrained:Exponential(10.0)


      Number of chains per MPI processor = 3

      The MCMC sampler will use the following moves:
         With prob.  Chain will use move
            2.08 %   Dirichlet(Tratio)
            1.04 %   Dirichlet(Pi)
            1.04 %   Slider(Pi)
            2.08 %   Multiplier(Alpha)
           10.42 %   ExtSPR(Tau,V)
           10.42 %   NNI(Tau,V)
           10.42 %   ParsSPR(Tau,V)
           41.67 %   Multiplier(V)
           14.58 %   Nodeslider(V)
            6.25 %   TLMultiplier(V)

      Division 1 has 101 unique site patterns
      Initializing conditional likelihoods

      Running benchmarks to automatically select fastest BEAGLE resource... A
B

      Using BEAGLE v4.0.0 (PRE-RELEASE) resource 0 for division 1:
         Rsrc Name : CPU (x86_64)
         Impl Name : CPU-4State-Single
         Flags: PROCESSOR_CPU PRECISION_SINGLE COMPUTATION_SYNCH EIGEN_REAL
                SCALING_MANUAL SCALERS_RAW VECTOR_NONE THREADING_CPP         MODEL STATES: 4

      Initial log likelihoods and log prior probs:
         Chain 1 -- -33245.493265 -- 10.206073
         Chain 2 -- -33245.493265 -- 10.206073
         Chain 3 -- -32408.396517 -- 10.206073


      Chain results (10000 generations requested):

         0 -- [-33245.493] (-33245.493) (-32408.397) (...0 remote chains...) 
       100 -- [-30574.157] (-31514.348) (-31493.071) (...0 remote chains...) -- 0:00:00
       200 -- [-30502.129] (-31426.403) (-31363.455) (...0 remote chains...) -- 0:00:00
       300 -- [-30497.743] (-31328.598) (-31271.241) (...0 remote chains...) -- 0:00:00
       400 -- [-30498.798] (-31097.301) (-31106.607) (...0 remote chains...) -- 0:00:00
       500 -- [-30433.854] (-31030.629) (-30586.731) (...0 remote chains...) -- 0:00:00
       600 -- [-30381.376] (-31034.211) (-30440.966) (...0 remote chains...) -- 0:00:00
       700 -- [-30382.833] (-30660.990) (-30443.911) (...0 remote chains...) -- 0:00:00
       800 -- [-30378.999] (-30597.996) (-30398.508) (...0 remote chains...) -- 0:00:00
       900 -- (-30373.832) (-30539.177) [-30382.084] (...0 remote chains...) -- 0:00:00
      1000 -- (-30372.743) (-30485.230) [-30375.672] (...0 remote chains...) -- 0:00:00
      1100 -- [-30364.175] (-30482.876) (-30373.687) (...0 remote chains...) -- 0:00:00
      1200 -- [-30362.105] (-30478.623) (-30377.733) (...0 remote chains...) -- 0:00:00
      1300 -- [-30363.247] (-30471.429) (-30377.862) (...0 remote chains...) -- 0:00:00
      1400 -- (-30365.586) [-30368.769] (-30373.465) (...0 remote chains...) -- 0:00:00
      1500 -- [-30348.693] (-30362.353) (-30374.050) (...0 remote chains...) -- 0:00:00
      1600 -- [-30346.301] (-30354.013) (-30371.333) (...0 remote chains...) -- 0:00:00
      1700 -- [-30346.350] (-30354.928) (-30370.900) (...0 remote chains...) -- 0:00:00
      1800 -- (-30350.440) [-30346.389] (-30370.874) (...0 remote chains...) -- 0:00:00
      1900 -- [-30350.238] (-30347.714) (-30369.555) (...0 remote chains...) -- 0:00:00
      2000 -- (-30351.098) [-30347.717] (-30370.075) (...0 remote chains...) -- 0:00:00
      2100 -- [-30346.447] (-30346.972) (-30370.269) (...0 remote chains...) -- 0:00:00
      2200 -- [-30351.456] (-30343.059) (-30366.168) (...0 remote chains...) -- 0:00:00
      2300 -- (-30347.251) [-30342.730] (-30357.691) (...0 remote chains...) -- 0:00:00
      2400 -- [-30348.080] (-30343.033) (-30361.555) (...0 remote chains...) -- 0:00:00
      2500 -- (-30348.126) [-30342.673] (-30352.976) (...0 remote chains...) -- 0:00:00
      2600 -- (-30348.207) [-30340.946] (-30353.354) (...0 remote chains...) -- 0:00:00
      2700 -- [-30342.129] (-30341.975) (-30345.593) (...0 remote chains...) -- 0:00:00
      2800 -- (-30343.278) [-30345.184] (-30349.762) (...0 remote chains...) -- 0:00:00
      2900 -- (-30343.790) [-30347.506] (-30348.727) (...0 remote chains...) -- 0:00:00
      3000 -- [-30343.962] (-30348.767) (-30347.823) (...0 remote chains...) -- 0:00:00
      3100 -- [-30341.685] (-30346.153) (-30352.217) (...0 remote chains...) -- 0:00:00
      3200 -- [-30343.118] (-30345.511) (-30356.461) (...0 remote chains...) -- 0:00:00
      3300 -- (-30343.917) [-30347.446] (-30351.640) (...0 remote chains...) -- 0:00:00
      3400 -- (-30346.596) [-30339.205] (-30350.682) (...0 remote chains...) -- 0:00:00
      3500 -- [-30340.202] (-30341.190) (-30355.712) (...0 remote chains...) -- 0:00:00
      3600 -- [-30342.479] (-30340.950) (-30351.605) (...0 remote chains...) -- 0:00:00
      3700 -- (-30339.130) [-30342.408] (-30353.638) (...0 remote chains...) -- 0:00:00
      3800 -- [-30342.209] (-30343.402) (-30340.686) (...0 remote chains...) -- 0:00:00
      3900 -- [-30338.862] (-30341.611) (-30344.425) (...0 remote chains...) -- 0:00:00
      4000 -- (-30340.339) [-30343.315] (-30340.701) (...0 remote chains...) -- 0:00:00
      4100 -- [-30339.302] (-30341.226) (-30343.133) (...0 remote chains...) -- 0:00:00
      4200 -- [-30340.923] (-30347.926) (-30342.142) (...0 remote chains...) -- 0:00:00
      4300 -- (-30341.860) [-30341.550] (-30340.374) (...0 remote chains...) -- 0:00:00
      4400 -- (-30341.802) [-30347.190] (-30342.353) (...0 remote chains...) -- 0:00:00
      4500 -- [-30339.748] (-30344.568) (-30342.229) (...0 remote chains...) -- 0:00:01
      4600 -- (-30342.603) (-30346.968) [-30343.194] (...0 remote chains...) -- 0:00:01
      4700 -- (-30339.753) (-30348.866) [-30338.306] (...0 remote chains...) -- 0:00:01
      4800 -- (-30343.637) (-30346.029) [-30340.197] (...0 remote chains...) -- 0:00:01
      4900 -- (-30338.944) (-30346.253) [-30338.791] (...0 remote chains...) -- 0:00:01
      5000 -- [-30342.981] (-30343.708) (-30341.554) (...0 remote chains...) -- 0:00:01
      5100 -- [-30339.869] (-30343.233) (-30340.883) (...0 remote chains...) -- 0:00:00
      5200 -- [-30340.460] (-30342.892) (-30339.163) (...0 remote chains...) -- 0:00:00
      5300 -- (-30339.978) (-30342.213) [-30340.832] (...0 remote chains...) -- 0:00:00
      5400 -- (-30343.073) (-30344.328) [-30339.825] (...0 remote chains...) -- 0:00:00
      5500 -- (-30341.117) (-30344.059) [-30338.220] (...0 remote chains...) -- 0:00:00
      5600 -- [-30339.571] (-30343.047) (-30337.876) (...0 remote chains...) -- 0:00:00
      5700 -- (-30342.817) [-30343.816] (-30340.152) (...0 remote chains...) -- 0:00:00
      5800 -- (-30345.714) (-30345.854) [-30342.811] (...0 remote chains...) -- 0:00:00
      5900 -- (-30339.548) (-30353.269) [-30342.465] (...0 remote chains...) -- 0:00:00
      6000 -- [-30341.459] (-30348.787) (-30342.128) (...0 remote chains...) -- 0:00:00
      6100 -- (-30340.127) [-30343.946] (-30340.774) (...0 remote chains...) -- 0:00:00
      6200 -- [-30344.733] (-30346.942) (-30342.648) (...0 remote chains...) -- 0:00:00
      6300 -- [-30338.907] (-30343.771) (-30339.861) (...0 remote chains...) -- 0:00:00
      6400 -- [-30339.709] (-30344.988) (-30338.547) (...0 remote chains...) -- 0:00:00
      6500 -- [-30341.693] (-30344.686) (-30340.217) (...0 remote chains...) -- 0:00:00
      6600 -- (-30338.872) (-30343.799) [-30340.339] (...0 remote chains...) -- 0:00:00
      6700 -- (-30339.780) [-30342.626] (-30344.485) (...0 remote chains...) -- 0:00:00
      6800 -- (-30340.132) [-30343.178] (-30344.642) (...0 remote chains...) -- 0:00:00
      6900 -- (-30339.216) (-30343.357) [-30341.135] (...0 remote chains...) -- 0:00:00
      7000 -- [-30340.391] (-30341.890) (-30344.584) (...0 remote chains...) -- 0:00:00
      7100 -- (-30340.258) [-30345.313] (-30340.289) (...0 remote chains...) -- 0:00:00
      7200 -- (-30338.674) [-30341.261] (-30341.619) (...0 remote chains...) -- 0:00:00
      7300 -- (-30337.344) (-30341.626) [-30340.537] (...0 remote chains...) -- 0:00:00
      7400 -- (-30338.820) [-30344.482] (-30340.571) (...0 remote chains...) -- 0:00:00
      7500 -- [-30340.280] (-30344.602) (-30343.185) (...0 remote chains...) -- 0:00:00
      7600 -- (-30341.731) (-30340.077) [-30342.225] (...0 remote chains...) -- 0:00:00
      7700 -- (-30341.098) (-30340.743) [-30340.160] (...0 remote chains...) -- 0:00:00
      7800 -- (-30342.786) [-30343.234] (-30339.736) (...0 remote chains...) -- 0:00:00
      7900 -- (-30341.353) (-30345.757) [-30340.918] (...0 remote chains...) -- 0:00:00
      8000 -- (-30344.145) (-30340.514) [-30339.884] (...0 remote chains...) -- 0:00:00
      8100 -- [-30341.402] (-30341.331) (-30339.047) (...0 remote chains...) -- 0:00:00
      8200 -- [-30340.394] (-30342.491) (-30341.522) (...0 remote chains...) -- 0:00:00
      8300 -- (-30342.827) (-30340.724) [-30340.903] (...0 remote chains...) -- 0:00:00
      8400 -- (-30343.657) [-30343.802] (-30340.506) (...0 remote chains...) -- 0:00:00
      8500 -- (-30345.233) [-30340.843] (-30343.745) (...0 remote chains...) -- 0:00:00
      8600 -- [-30340.335] (-30340.938) (-30339.245) (...0 remote chains...) -- 0:00:00
      8700 -- [-30340.179] (-30343.506) (-30339.883) (...0 remote chains...) -- 0:00:00
      8800 -- (-30339.265) (-30345.651) [-30338.885] (...0 remote chains...) -- 0:00:00
      8900 -- [-30339.781] (-30346.047) (-30341.482) (...0 remote chains...) -- 0:00:00
      9000 -- [-30337.946] (-30348.245) (-30343.755) (...0 remote chains...) -- 0:00:00
      9100 -- (-30338.659) (-30346.907) [-30341.201] (...0 remote chains...) -- 0:00:00
      9200 -- (-30338.787) (-30347.979) [-30341.754] (...0 remote chains...) -- 0:00:00
      9300 -- [-30340.557] (-30347.858) (-30343.602) (...0 remote chains...) -- 0:00:00
      9400 -- [-30339.744] (-30347.391) (-30342.387) (...0 remote chains...) -- 0:00:00
      9500 -- [-30344.064] (-30348.410) (-30339.898) (...0 remote chains...) -- 0:00:00
      9600 -- [-30339.095] (-30353.036) (-30340.072) (...0 remote chains...) -- 0:00:00
      9700 -- (-30340.613) (-30351.953) [-30340.951] (...0 remote chains...) -- 0:00:00
      9800 -- [-30340.730] (-30349.857) (-30340.810) (...0 remote chains...) -- 0:00:00
      9900 -- [-30343.658] (-30348.051) (-30340.968) (...0 remote chains...) -- 0:00:00
      10000 -- (-30340.070) [-30344.198] (-30342.027) (...0 remote chains...) -- 0:00:00

      Analysis completed in 1 second
      Analysis used 1.16 seconds of CPU time on processor 0
      Log likelihood of best state for "cold" chain was -30337.28

      Acceptance rates for the moves in the "cold" chain:
         With prob.   (last 100)   chain accepted proposals by move
            17.6 %     ( 16 %)     Dirichlet(Tratio)
             5.0 %     (  5 %)     Dirichlet(Pi)
            15.0 %     ( 15 %)     Slider(Pi)
            86.5 %     ( 85 %)     Multiplier(Alpha)
             0.0 %     (  0 %)     ExtSPR(Tau,V)
             0.1 %     (  0 %)     NNI(Tau,V)
             0.0 %     (  0 %)     ParsSPR(Tau,V)
            25.3 %     ( 26 %)     Multiplier(V)
             4.3 %     (  8 %)     Nodeslider(V)
             6.6 %     (  8 %)     TLMultiplier(V)

      Chain swap information:

                 1     2     3 
           --------------------
         1 |        0.74  0.56 
         2 |  3325        0.69 
         3 |  3382  3293       

      Upper diagonal: Proportion of successful state exchanges between chains
      Lower diagonal: Number of attempted state exchanges between chains

      Chain information:

        ID -- Heat 
       -----------
         1 -- 1.00  (cold chain)
         2 -- 0.91 
         3 -- 0.83 

      Heat = 1 / (1 + T * (ID - 1))
         (where T = 0.10 is the temperature and ID is the chain number)

      Summarizing trees in file "SEQ-ALIGNED-MB.nex.t"
      Using relative burnin ('relburnin=yes'), discarding the first 25 % of sampled trees
      Writing statistics to files SEQ-ALIGNED-MB.nex.<parts|tstat|vstat|trprobs|con>
      Examining file ...
      Found one tree block in file "SEQ-ALIGNED-MB.nex.t" with 1001 trees in last block

      Tree reading status:

      0      10      20      30      40      50      60      70      80      90     100
      v-------v-------v-------v-------v-------v-------v-------v-------v-------v-------v
      *********************************************************************************

      Read 1001 trees from last tree block (sampling 751 of them)
                                                                                   
      General explanation:                                                          
                                                                                   
      In an unrooted tree, a taxon bipartition (split) is specified by removing a   
      branch, thereby dividing the species into those to the left and those to the  
      right of the branch. Here, taxa to one side of the removed branch are denoted 
      '.' and those to the other side are denoted '*'. Specifically, the '.' symbol 
      is used for the taxa on the same side as the outgroup.                        
                                                                                   
      In a rooted or clock tree, the tree is rooted using the model and not by      
      reference to an outgroup. Each bipartition therefore corresponds to a clade,  
      that is, a group that includes all the descendants of a particular branch in  
      the tree.  Taxa that are included in each clade are denoted using '*', and    
      taxa that are not included are denoted using the '.' symbol.                  
                                                                                   
      The output first includes a key to all the bipartitions with frequency larger 
      or equual to (Minpartfreq) in at least one run. Minpartfreq is a parameter to 
      sumt command and currently it is set to 0.10.  This is followed by a table  
      with statistics for the informative bipartitions (those including at least    
      two taxa), sorted from highest to lowest probability. For each bipartition,   
      the table gives the number of times the partition or split was observed in all
      runs (#obs) and the posterior probability of the bipartition (Probab.), which 
      is the same as the split frequency. If several runs are summarized, this is   
      followed by the minimum split frequency (Min(s)), the maximum frequency       
      (Max(s)), and the standard deviation of frequencies (Stddev(s)) across runs.  
      The latter value should approach 0 for all bipartitions as MCMC runs converge.
                                                                                   
      This is followed by a table summarizing branch lengths, node heights (if a    
      clock model was used) and relaxed clock parameters (if a relaxed clock model  
      was used). The mean, variance, and 95 % credible interval are given for each 
      of these parameters. If several runs are summarized, the potential scale      
      reduction factor (PSRF) is also given; it should approach 1 as runs converge. 
      Node heights will take calibration points into account, if such points were   
      used in the analysis.                                                         
                                                                                    
      Note that Stddev may be unreliable if the partition is not present in all     
      runs (the last column indicates the number of runs that sampled the partition 
      if more than one run is summarized). The PSRF is not calculated at all if     
      the partition is not present in all runs.The PSRF is also sensitive to small  
      sample sizes and it should only be considered a rough guide to convergence    
      since some of the assumptions allowing one to interpret it as a true potential
      scale reduction factor are violated in MrBayes.                               
                                                                                    
      List of taxa in bipartitions:                                                 
                                                                                   
         1 -- MD
         2 -- MF
         3 -- MS
         4 -- MP

      Key to taxon bipartitions (saved to file "SEQ-ALIGNED-MB.nex.parts"):

      ID -- Partition
      ----------
       1 -- .***
       2 -- .*..
       3 -- ..*.
       4 -- ...*
       5 -- ..**
      ----------

      Summary statistics for informative taxon bipartitions
         (saved to file "SEQ-ALIGNED-MB.nex.tstat"):

      ID   #obs    Probab.
      --------------------
       5   751    1.000000
      --------------------

      Summary statistics for branch and node parameters
         (saved to file "SEQ-ALIGNED-MB.nex.vstat"):

                                             95% HPD Interval
                                           --------------------
      Parameter      Mean       Variance     Lower       Upper       Median
      --------------------------------------------------------------------
      length[1]    0.088194    0.000020    0.079580    0.096141    0.088058
      length[2]    0.017747    0.000002    0.015041    0.020562    0.017937
      length[3]    0.003583    0.000000    0.002440    0.004572    0.003536
      length[4]    0.002947    0.000000    0.002119    0.003877    0.002950
      length[5]    0.016467    0.000002    0.013661    0.019259    0.016358
      --------------------------------------------------------------------




      Clade credibility values:

      /------------------------------------------------------------------------ MD (1)
      |                                                                               
      |------------------------------------------------------------------------ MF (2)
      +                                                                               
      |                                   /------------------------------------ MS (3)
      \----------------100----------------+                                           
                                          \------------------------------------ MP (4)
                                                                                      

      Phylogram (based on average branch lengths):

      /------------------------------------------------------------------------ MD (1)
      |                                                                               
      |--------------- MF (2)
      +                                                                               
      |            /--- MS (3)
      \------------+                                                                  
                   \--- MP (4)
                                                                                      
      |-------| 0.010 expected changes per site

      Calculating tree probabilities...

      Credible sets of trees (1 tree sampled):
         99 % credible set contains 1 tree

   Exiting mrbayes block
   Reached end of file

   Tasks completed, exiting program because mode is noninteractive
   To return control to the command line after completion of file processing, 
   set mode to interactive with 'mb -i <filename>' (i is for interactive)
   or use 'set mode=interactive'

## IQTree
I ran IQTree through terminal. It seemed I had previously downloaded the packages because I did not have to do any installation for this project.
    Running IQTree on my data:
    `iqtree -s SEQ-ALIGNED.fasta -bb 1000 -nt AUTO`
    Output:
    IQ-TREE multicore version 2.1.4-beta COVID-edition for Mac OS X 64-bit built Jun 24 2021
Developed by Bui Quang Minh, James Barbetti, Nguyen Lam Tung,
Olga Chernomor, Heiko Schmidt, Dominik Schrempf, Michael Woodhams.

Host:    AG (AVX2, FMA3, 8 GB RAM)
Command: iqtree -s SEQ-ALIGNED.fasta -bb 1000 -nt AUTO
Seed:    27646 (Using SPRNG - Scalable Parallel Random Number Generator)
Time:    Mon May  1 13:22:10 2023
OMP: Info #276: omp_set_nested routine deprecated, please use omp_set_max_active_levels instead.
Kernel:  AVX+FMA - auto-detect threads (4 CPU cores detected)

Reading alignment file SEQ-ALIGNED.fasta ... Fasta format detected
Alignment most likely contains DNA/RNA sequences
Alignment has 4 sequences with 16497 columns, 101 distinct patterns
255 parsimony-informative, 1269 singleton sites, 14973 constant sites
             Gap/Ambiguity  Composition  p-value
   1  NC_011361.1    0.05%    passed     91.79%
   2  HQ391897.1     0.11%    passed     99.22%
   3  OK945936.1     0.04%    passed     99.90%
   4  OL693877.1     0.06%    passed     99.39%
****  TOTAL          0.07%  0 sequences failed composition chi2 test (p-value<5%; df=3)


Create initial parsimony tree by phylogenetic likelihood library (PLL)... 0.001 seconds
Measuring multi-threading efficiency up to 4 CPU cores
Increase to 10 rounds for branch lengths
13512 trees examined
Threads: 1 / Time: 4.171 sec / Speedup: 1.000 / Efficiency: 100% / LogL: -31485
Threads: 2 / Time: 7.548 sec / Speedup: 0.553 / Efficiency: 28% / LogL: -31485
BEST NUMBER OF THREADS: 1

Perform fast likelihood tree search using GTR+I+G model...
Estimate model parameters (epsilon = 5.000)
Perform nearest neighbor interchange...
Estimate model parameters (epsilon = 1.000)
1. Initial log-likelihood: -30317.687
Optimal log-likelihood: -30317.460
Rate parameters:  A-C: 1.36778  A-G: 23.71026  A-T: 2.73055  C-G: 1.03164  C-T: 19.40384  G-T: 1.00000
Base frequencies:  A: 0.267  C: 0.296  G: 0.170  T: 0.267
Proportion of invariable sites: 0.450
Gamma shape alpha: 0.485
Parameters optimization took 1 rounds (0.002 sec)
Time for fast ML tree search: 0.014 seconds

NOTE: ModelFinder requires 0 MB RAM!
ModelFinder will test up to 286 DNA models (sample size: 16497) ...
 No. Model         -LnL         df  AIC          AICc         BIC
  1  GTR+F         30371.060    13  60768.119    60768.141    60868.361
  2  GTR+F+I       30316.994    14  60661.988    60662.013    60769.941
  3  GTR+F+G4      30317.012    14  60662.025    60662.050    60769.978
  4  GTR+F+I+G4    30317.463    15  60664.925    60664.954    60780.589
  5  GTR+F+R2      30319.478    15  60668.957    60668.986    60784.621
  6  GTR+F+R3      30317.983    17  60669.967    60670.004    60801.053
 15  SYM+I         30654.933    11  61331.866    61331.882    61416.687
 16  SYM+G4        30655.716    11  61333.431    61333.447    61418.252
 28  TVM+F+I       30321.546    13  60669.093    60669.115    60769.335
 29  TVM+F+G4      30321.564    13  60669.127    60669.150    60769.370
 41  TVMe+I        30656.223    10  61332.445    61332.459    61409.555
 42  TVMe+G4       30656.581    10  61333.161    61333.174    61410.270
 54  TIM3+F+I      30327.068    12  60678.135    60678.154    60770.667
 55  TIM3+F+G4     30327.269    12  60678.539    60678.558    60771.070
 67  TIM3e+I       30676.948    9   61371.897    61371.907    61441.295
 68  TIM3e+G4      30677.028    9   61372.056    61372.067    61441.455
 80  TIM2+F+I      30323.293    12  60670.585    60670.604    60763.117
 81  TIM2+F+G4     30323.415    12  60670.831    60670.850    60763.362
 93  TIM2e+I       30659.943    9   61337.885    61337.896    61407.284
 94  TIM2e+G4      30660.110    9   61338.219    61338.230    61407.618
106  TIM+F+I       30327.279    12  60678.558    60678.577    60771.090
107  TIM+F+G4      30327.496    12  60678.992    60679.011    60771.523
119  TIMe+I        30676.326    9   61370.652    61370.663    61440.050
120  TIMe+G4       30676.379    9   61370.759    61370.770    61440.157
132  TPM3u+F+I     30331.970    11  60685.940    60685.956    60770.760
133  TPM3u+F+G4    30332.058    11  60686.115    60686.131    60770.936
145  TPM3+F+I      30331.970    11  60685.940    60685.956    60770.760
146  TPM3+F+G4     30332.058    11  60686.115    60686.131    60770.936
158  TPM2u+F+I     30328.023    11  60678.046    60678.062    60762.867
159  TPM2u+F+G4    30328.101    11  60678.202    60678.218    60763.022
171  TPM2+F+I      30328.010    11  60678.019    60678.035    60762.840
172  TPM2+F+G4     30328.098    11  60678.196    60678.212    60763.016
184  K3Pu+F+I      30332.182    11  60686.364    60686.380    60771.184
185  K3Pu+F+G4     30332.339    11  60686.679    60686.695    60771.499
197  K3P+I         30677.606    8   61371.211    61371.220    61432.899
198  K3P+G4        30677.732    8   61371.464    61371.472    61433.151
210  TN+F+I        30331.859    11  60685.718    60685.734    60770.539
211  TN+F+G4       30331.944    11  60685.888    60685.904    60770.708
223  TNe+I         30680.456    8   61376.912    61376.921    61438.600
224  TNe+G4        30680.465    8   61376.929    61376.938    61438.617
236  HKY+F+I       30336.863    10  60693.727    60693.740    60770.836
237  HKY+F+G4      30336.903    10  60693.806    60693.819    60770.915
249  K2P+I         30681.717    7   61377.434    61377.441    61431.411
250  K2P+G4        30681.724    7   61377.448    61377.455    61431.425
262  F81+F+I       31215.479    9   62448.958    62448.969    62518.357
263  F81+F+G4      31215.613    9   62449.227    62449.238    62518.625
275  JC+I          31564.406    6   63140.812    63140.817    63187.078
276  JC+G4         31564.507    6   63141.015    63141.020    63187.280
Akaike Information Criterion:           GTR+F+I
Corrected Akaike Information Criterion: GTR+F+I
Bayesian Information Criterion:         TPM2+F+I
Best-fit model: TPM2+F+I chosen according to BIC

All model information printed to SEQ-ALIGNED.fasta.model.gz
CPU time for ModelFinder: 10.810 seconds (0h:0m:10s)
Wall-clock time for ModelFinder: 12.071 seconds (0h:0m:12s)
Generating 1000 samples for ultrafast bootstrap (seed: 27646)...

NOTE: 0 MB RAM (0 GB) is required!
Estimate model parameters (epsilon = 0.100)
1. Initial log-likelihood: -30328.010
Optimal log-likelihood: -30328.004
Rate parameters:  A-C: 1.96594  A-G: 20.97964  A-T: 1.96594  C-G: 1.00000  C-T: 20.97964  G-T: 1.00000
Base frequencies:  A: 0.267  C: 0.296  G: 0.170  T: 0.267
Proportion of invariable sites: 0.754
Parameters optimization took 1 rounds (0.001 sec)
Computing ML distances based on estimated model parameters...
Computing ML distances took 0.001914 sec (of wall-clock time) 0.000193 sec(of CPU time)
OMP: Info #276: omp_set_nested routine deprecated, please use omp_set_max_active_levels instead.
Computing RapidNJ tree took 0.000457 sec (of wall-clock time) 0.000419 sec (of CPU time)
OMP: Info #276: omp_set_nested routine deprecated, please use omp_set_max_active_levels instead.
Log-likelihood of RapidNJ tree: -30328.004
--------------------------------------------------------------------
|             INITIALIZING CANDIDATE TREE SET                      |
--------------------------------------------------------------------
Generating 99 parsimony trees... 0.035 second
Computing log-likelihood of 2 initial trees ... 0.000 seconds
Current best score: -30328.004

Do NNI search on 3 best initial trees
Finish initializing candidate tree set (3)
Current best tree score: -30328.004 / CPU time: 0.043
Number of iterations: 3
--------------------------------------------------------------------
|               OPTIMIZING CANDIDATE TREE SET                      |
--------------------------------------------------------------------
UPDATE BEST LOG-LIKELIHOOD: -30328.004
Iteration 10 / LogL: -30328.004 / Time: 0h:0m:0s
Iteration 20 / LogL: -30424.509 / Time: 0h:0m:0s
Iteration 30 / LogL: -30328.332 / Time: 0h:0m:0s (0h:0m:0s left)
Iteration 40 / LogL: -30423.810 / Time: 0h:0m:0s (0h:0m:0s left)
Iteration 50 / LogL: -30328.265 / Time: 0h:0m:0s (0h:0m:0s left)
Log-likelihood cutoff on original alignment: -30337.121
Iteration 60 / LogL: -30328.266 / Time: 0h:0m:0s (0h:0m:0s left)
Iteration 70 / LogL: -30424.509 / Time: 0h:0m:0s (0h:0m:0s left)
Iteration 80 / LogL: -30328.266 / Time: 0h:0m:0s (0h:0m:0s left)
Iteration 90 / LogL: -30328.332 / Time: 0h:0m:0s (0h:0m:0s left)
Iteration 100 / LogL: -30423.810 / Time: 0h:0m:0s (0h:0m:0s left)
Log-likelihood cutoff on original alignment: -30337.121
NOTE: Bootstrap correlation coefficient of split occurrence frequencies: 1.000
TREE SEARCH COMPLETED AFTER 101 ITERATIONS / Time: 0h:0m:0s

--------------------------------------------------------------------
|                    FINALIZING TREE SEARCH                        |
--------------------------------------------------------------------
Performs final model parameters optimization
Estimate model parameters (epsilon = 0.010)
1. Initial log-likelihood: -30328.004
2. Current log-likelihood: -30327.990
Optimal log-likelihood: -30327.989
Rate parameters:  A-C: 1.97721  A-G: 21.10609  A-T: 1.97721  C-G: 1.00000  C-T: 21.10609  G-T: 1.00000
Base frequencies:  A: 0.267  C: 0.296  G: 0.170  T: 0.267
Proportion of invariable sites: 0.754
Parameters optimization took 2 rounds (0.001 sec)
BEST SCORE FOUND : -30327.989
Creating bootstrap support values...
Split supports printed to NEXUS file SEQ-ALIGNED.fasta.splits.nex
Total tree length: 0.131

Total number of iterations: 101
CPU time used for tree search: 0.518 sec (0h:0m:0s)
Wall-clock time used for tree search: 0.372 sec (0h:0m:0s)
Total CPU time used: 0.554 sec (0h:0m:0s)
Total wall-clock time used: 0.416 sec (0h:0m:0s)

Computing bootstrap consensus tree...
Reading input file SEQ-ALIGNED.fasta.splits.nex...
4 taxa and 5 splits.
Consensus tree written to SEQ-ALIGNED.fasta.contree
Reading input trees file SEQ-ALIGNED.fasta.contree
Log-likelihood of consensus tree: -30327.989

Analysis results written to: 
  IQ-TREE report:                SEQ-ALIGNED.fasta.iqtree
  Maximum-likelihood tree:       SEQ-ALIGNED.fasta.treefile
  Likelihood distances:          SEQ-ALIGNED.fasta.mldist

Ultrafast bootstrap approximation results written to:
  Split support values:          SEQ-ALIGNED.fasta.splits.nex
  Consensus tree:                SEQ-ALIGNED.fasta.contree
  Screen log file:               SEQ-ALIGNED.fasta.log


REFERENCES
---------- 

Subha Kalyaanamoorthy, Bui Quang Minh, Thomas KF Wong, Arndt von Haeseler,
and Lars S Jermiin (2017) ModelFinder: Fast model selection for
accurate phylogenetic estimates. Nature Methods, 14:587–589.
https://doi.org/10.1038/nmeth.4285


Lam-Tung Nguyen, Heiko A. Schmidt, Arndt von Haeseler, and Bui Quang Minh
(2015) IQ-TREE: A fast and effective stochastic algorithm for estimating
maximum likelihood phylogenies. Mol Biol Evol, 32:268-274.
https://doi.org/10.1093/molbev/msu300


Diep Thi Hoang, Olga Chernomor, Arndt von Haeseler, Bui Quang Minh,
and Le Sy Vinh (2017) UFBoot2: Improving the ultrafast bootstrap
approximation. Mol Biol Evol, in press.
https://doi.org/10.1093/molbev/msx281

SEQUENCE ALIGNMENT
------------------

Input data: 4 sequences with 16497 nucleotide sites
Number of constant sites: 14973 (= 90.762% of all sites)
Number of invariant (constant or ambiguous constant) sites: 14973 (= 90.762% of all sites)
Number of parsimony informative sites: 255
Number of distinct site patterns: 101

ModelFinder
-----------

Best-fit model according to BIC: TPM2u+F+I

List of models sorted by BIC scores: 

Model             LogL          AIC      w-AIC      AICc     w-AICc       BIC      w-BIC
TPM2u+F+I       -30327.9510  60677.9020 - 0.0001  60677.9181 - 0.0001  60762.7223 + 0.1621
TPM2+F+I        -30327.9617  60677.9233 - 0.0001  60677.9393 - 0.0001  60762.7436 + 0.1604
TPM2u+F+G4      -30327.9667  60677.9334 - 0.0001  60677.9494 - 0.0001  60762.7537 + 0.1596
TPM2+F+G4       -30327.9726  60677.9452 - 0.0001  60677.9613 - 0.0001  60762.7655 + 0.1587
TIM2+F+I        -30323.1849  60670.3699 - 0.0051  60670.3888 - 0.0052  60762.9011 + 0.1483
TIM2+F+G4       -30323.1964  60670.3929 - 0.0051  60670.4118 - 0.0051  60762.9241 + 0.1466
TVM+F+I         -30321.3816  60668.7632 - 0.0115  60668.7853 - 0.0115  60769.0054 - 0.0070
TVM+F+G4        -30321.4344  60668.8688 - 0.0109  60668.8909 - 0.0109  60769.1110 - 0.0066
GTR+F+I         -30316.7368  60661.4736 + 0.4399  60661.4991 + 0.4400  60769.4267 - 0.0057
GTR+F+G4        -30316.7768  60661.5535 + 0.4227  60661.5790 + 0.4227  60769.5066 - 0.0055
TN+F+I          -30331.9062  60685.8123 - 0.0000  60685.8284 - 0.0000  60770.6326 - 0.0031
TN+F+G4         -30331.9447  60685.8893 - 0.0000  60685.9053 - 0.0000  60770.7096 - 0.0030
TPM3u+F+I       -30331.9754  60685.9509 - 0.0000  60685.9669 - 0.0000  60770.7712 - 0.0029
TPM3+F+I        -30331.9757  60685.9514 - 0.0000  60685.9674 - 0.0000  60770.7717 - 0.0029
TIM3+F+I        -30327.1340  60678.2681 - 0.0001  60678.2870 - 0.0001  60770.7993 - 0.0029
HKY+F+I         -30336.8976  60693.7952 - 0.0000  60693.8085 - 0.0000  60770.9045 - 0.0027
TPM3+F+G4       -30332.0422  60686.0843 - 0.0000  60686.1003 - 0.0000  60770.9046 - 0.0027
HKY+F+G4        -30336.9030  60693.8059 - 0.0000  60693.8193 - 0.0000  60770.9153 - 0.0027
TPM3u+F+G4      -30332.0498  60686.0996 - 0.0000  60686.1157 - 0.0000  60770.9199 - 0.0027
TIM3+F+G4       -30327.2140  60678.4280 - 0.0001  60678.4469 - 0.0001  60770.9592 - 0.0026
TIM+F+G4        -30327.2695  60678.5391 - 0.0001  60678.5580 - 0.0001  60771.0703 - 0.0025
K3Pu+F+G4       -30332.3457  60686.6914 - 0.0000  60686.7074 - 0.0000  60771.5117 - 0.0020
K3Pu+F+I        -30332.3503  60686.7006 - 0.0000  60686.7166 - 0.0000  60771.5209 - 0.0020
TIM+F+I         -30327.5157  60679.0314 - 0.0001  60679.0504 - 0.0001  60771.5626 - 0.0020
TPM2u+F+I+G4    -30328.1326  60680.2651 - 0.0000  60680.2840 - 0.0000  60772.7963 - 0.0011
TPM2+F+I+G4     -30328.1472  60680.2944 - 0.0000  60680.3133 - 0.0000  60772.8256 - 0.0010
TIM2+F+I+G4     -30323.6551  60673.3103 - 0.0012  60673.3324 - 0.0012  60773.5524 - 0.0007
TVM+F+I+G4      -30321.5714  60671.1427 - 0.0035  60671.1682 - 0.0035  60779.0958 - 0.0000
GTR+F+I+G4      -30317.2268  60664.4535 + 0.0991  60664.4826 + 0.0990  60780.1175 - 0.0000
TPM3u+F+I+G4    -30332.1331  60688.2662 - 0.0000  60688.2851 - 0.0000  60780.7974 - 0.0000
TPM3+F+I+G4     -30332.1564  60688.3127 - 0.0000  60688.3316 - 0.0000  60780.8439 - 0.0000
HKY+F+I+G4      -30337.1433  60696.2866 - 0.0000  60696.3026 - 0.0000  60781.1069 - 0.0000
TN+F+I+G4       -30332.4734  60688.9468 - 0.0000  60688.9657 - 0.0000  60781.4780 - 0.0000
TIM3+F+I+G4     -30327.6652  60681.3304 - 0.0000  60681.3525 - 0.0000  60781.5725 - 0.0000
K3Pu+F+I+G4     -30332.6585  60689.3170 - 0.0000  60689.3359 - 0.0000  60781.8482 - 0.0000
TIM+F+I+G4      -30327.8441  60681.6881 - 0.0000  60681.7102 - 0.0000  60781.9303 - 0.0000
TIM2+F          -30377.1448  60776.2896 - 0.0000  60776.3056 - 0.0000  60861.1099 - 0.0000
TPM2+F          -30382.1854  60784.3708 - 0.0000  60784.3842 - 0.0000  60861.4802 - 0.0000
TPM2u+F         -30382.1952  60784.3903 - 0.0000  60784.4037 - 0.0000  60861.4997 - 0.0000
GTR+F           -30371.0576  60768.1153 - 0.0000  60768.1373 - 0.0000  60868.3574 - 0.0000
TVM+F           -30376.1043  60776.2086 - 0.0000  60776.2275 - 0.0000  60868.7398 - 0.0000
TN+F            -30386.2125  60792.4250 - 0.0000  60792.4383 - 0.0000  60869.5343 - 0.0000
HKY+F           -30391.2375  60800.4750 - 0.0000  60800.4859 - 0.0000  60869.8734 - 0.0000
TIM+F           -30382.0615  60786.1231 - 0.0000  60786.1391 - 0.0000  60870.9433 - 0.0000
TIM3+F          -30382.0968  60786.1936 - 0.0000  60786.2096 - 0.0000  60871.0139 - 0.0000
TPM3u+F         -30387.0633  60794.1266 - 0.0000  60794.1399 - 0.0000  60871.2359 - 0.0000
TPM3+F          -30387.0669  60794.1338 - 0.0000  60794.1471 - 0.0000  60871.2431 - 0.0000
K3Pu+F          -30387.0915  60794.1830 - 0.0000  60794.1963 - 0.0000  60871.2923 - 0.0000
TIM2e+I         -30659.9377  61337.8755 - 0.0000  61337.8864 - 0.0000  61407.2739 - 0.0000
TIM2e+G4        -30659.9588  61337.9176 - 0.0000  61337.9286 - 0.0000  61407.3160 - 0.0000
TVMe+G4         -30656.2115  61332.4230 - 0.0000  61332.4364 - 0.0000  61409.5324 - 0.0000
TVMe+I          -30656.2136  61332.4271 - 0.0000  61332.4404 - 0.0000  61409.5364 - 0.0000
TIM2e+I+G4      -30660.0272  61340.0544 - 0.0000  61340.0677 - 0.0000  61417.1637 - 0.0000
SYM+I           -30655.2250  61332.4500 - 0.0000  61332.4660 - 0.0000  61417.2703 - 0.0000
SYM+G4          -30655.2351  61332.4701 - 0.0000  61332.4862 - 0.0000  61417.2904 - 0.0000
TVMe+I+G4       -30656.4686  61334.9372 - 0.0000  61334.9532 - 0.0000  61419.7575 - 0.0000
SYM+I+G4        -30655.3529  61334.7058 - 0.0000  61334.7247 - 0.0000  61427.2370 - 0.0000
K2P+I           -30681.7337  61377.4674 - 0.0000  61377.4742 - 0.0000  61431.4439 - 0.0000
K2P+G4          -30681.7377  61377.4755 - 0.0000  61377.4823 - 0.0000  61431.4520 - 0.0000
K3P+I           -30677.5991  61371.1981 - 0.0000  61371.2069 - 0.0000  61432.8856 - 0.0000
K3P+G4          -30677.6024  61371.2047 - 0.0000  61371.2134 - 0.0000  61432.8922 - 0.0000
TNe+I           -30680.4281  61376.8562 - 0.0000  61376.8649 - 0.0000  61438.5436 - 0.0000
TNe+G4          -30680.4311  61376.8622 - 0.0000  61376.8709 - 0.0000  61438.5497 - 0.0000
TIMe+I          -30676.3055  61370.6110 - 0.0000  61370.6219 - 0.0000  61440.0094 - 0.0000
TIMe+G4         -30676.3307  61370.6614 - 0.0000  61370.6723 - 0.0000  61440.0598 - 0.0000
TIM3e+I         -30676.8672  61371.7344 - 0.0000  61371.7453 - 0.0000  61441.1328 - 0.0000
TIM3e+G4        -30676.8947  61371.7893 - 0.0000  61371.8003 - 0.0000  61441.1877 - 0.0000
K2P+I+G4        -30682.0746  61380.1492 - 0.0000  61380.1580 - 0.0000  61441.8367 - 0.0000
K3P+I+G4        -30678.0266  61374.0533 - 0.0000  61374.0642 - 0.0000  61443.4517 - 0.0000
TNe+I+G4        -30680.6238  61379.2475 - 0.0000  61379.2585 - 0.0000  61448.6459 - 0.0000
TIMe+I+G4       -30676.6046  61373.2091 - 0.0000  61373.2224 - 0.0000  61450.3184 - 0.0000
TIM3e+I+G4      -30677.0170  61374.0339 - 0.0000  61374.0472 - 0.0000  61451.1432 - 0.0000
TIM2e           -30709.7781  61435.5562 - 0.0000  61435.5649 - 0.0000  61497.2436 - 0.0000
TVMe            -30707.1851  61432.3702 - 0.0000  61432.3811 - 0.0000  61501.7686 - 0.0000
SYM             -30705.2559  61430.5117 - 0.0000  61430.5251 - 0.0000  61507.6211 - 0.0000
K2P             -30731.7891  61475.5782 - 0.0000  61475.5833 - 0.0000  61521.8438 - 0.0000
K3P             -30728.1935  61470.3870 - 0.0000  61470.3938 - 0.0000  61524.3635 - 0.0000
TNe             -30729.4711  61472.9422 - 0.0000  61472.9490 - 0.0000  61526.9187 - 0.0000
TIMe            -30725.9230  61467.8461 - 0.0000  61467.8548 - 0.0000  61529.5336 - 0.0000
TIM3e           -30726.2772  61468.5544 - 0.0000  61468.5632 - 0.0000  61530.2419 - 0.0000
F81+F+I         -31215.4791  62448.9583 - 0.0000  62448.9692 - 0.0000  62518.3567 - 0.0000
F81+F+G4        -31215.6137  62449.2274 - 0.0000  62449.2383 - 0.0000  62518.6258 - 0.0000
F81+F+I+G4      -31215.5847  62451.1694 - 0.0000  62451.1828 - 0.0000  62528.2788 - 0.0000
F81+F           -31235.3094  62486.6188 - 0.0000  62486.6275 - 0.0000  62548.3062 - 0.0000
JC+I            -31564.4061  63140.8122 - 0.0000  63140.8173 - 0.0000  63187.0778 - 0.0000
JC+G4           -31564.5074  63141.0148 - 0.0000  63141.0199 - 0.0000  63187.2804 - 0.0000
JC+I+G4         -31564.4803  63142.9607 - 0.0000  63142.9675 - 0.0000  63196.9372 - 0.0000
JC              -31583.4467  63176.8933 - 0.0000  63176.8970 - 0.0000  63215.4480 - 0.0000

AIC, w-AIC   : Akaike information criterion scores and weights.
AICc, w-AICc : Corrected AIC scores and weights.
BIC, w-BIC   : Bayesian information criterion scores and weights.

Plus signs denote the 95% confidence sets.
Minus signs denote significant exclusion.

SUBSTITUTION PROCESS
--------------------

Model of substitution: TPM2u+F+I

Rate parameter R:

  A-C: 1.9307
  A-G: 20.5226
  A-T: 1.9307
  C-G: 1.0000
  C-T: 20.5226
  G-T: 1.0000

State frequencies: (empirical counts from alignment)

  pi(A) = 0.2675
  pi(C) = 0.2958
  pi(G) = 0.1697
  pi(T) = 0.267

Rate matrix Q:

  A   -0.7776   0.09716    0.5927    0.0877
  C   0.08787    -1.049   0.02888    0.9323
  G    0.9341   0.05033     -1.03   0.04543
  T   0.08787     1.033   0.02888     -1.15

Model of rate heterogeneity: Invar
Proportion of invariable sites: 0.7513

 Category  Relative_rate  Proportion
  0         0              0.7513
  1         4.021          0.2487

MAXIMUM LIKELIHOOD TREE
-----------------------

Log-likelihood of the tree: -30328.0707 (s.e. 163.3832)
Unconstrained log-likelihood (without tree): -30421.3891
Number of free parameters (#branches + #model parameters): 11
Akaike information criterion (AIC) score: 60678.1415
Corrected Akaike information criterion (AICc) score: 60678.1575
Bayesian information criterion (BIC) score: 60762.9618

Total tree length (sum of branch lengths): 0.1299
Sum of internal branch lengths: 0.0165 (12.6722% of tree length)

NOTE: Tree is UNROOTED although outgroup taxon 'MD' is drawn at root
Numbers in parentheses are SH-aLRT support (%) / ultrafast bootstrap support (%)

+-----------------------------------------------------------MD
|
+----------MF
|
|          +--MS
+----------| (100/100)
           +--MP

Tree in newick format:

(MD:0.0893201427,MF:0.0177099747,(MS:0.0035299700,MP:0.0028660135)100/100:0.0164592820);

CONSENSUS TREE
--------------

Consensus tree is constructed from 1000bootstrap trees
Log-likelihood of consensus tree: -30328.070727
Robinson-Foulds distance between ML tree and consensus tree: 0

Branches with support >0.000000% are kept (extended consensus)
Branch lengths are optimized by maximum likelihood on original alignment
Numbers in parentheses are bootstrap supports (%)

+-----------------------------------------------------------MD
|
+----------MF
|
|          +--MS
+----------| (100)
           +--MP


Consensus tree in newick format: 

(MD:0.0893281292,MF:0.0177106782,(MS:0.0035307578,MP:0.0028683415)100:0.0164587324);
