## My Data
I chose data on the complete genome of the freshwater bass family.
I chose this because my research is in water toxicology, so I do a lot of work with fish in the lab. I've never gotten to work with bass, but I found some data on one of them and that kind of got me interested.
GenBank had only 4 different species' complete mitochondrial genome, which was a little sad but was okay because I then had something I could still work with.
## MUSCLE
I ran it through MUSCLE to get an aligned sequence. I chose MUSCLE mostly because it was the only one I could get to work with my data and on my generation of MacBook.
Download MUSCLE:
    muscle3.8.31_i86darwin64.tar.gz
Using MUSCLE:
 /Users/addisongibbs/Desktop/Botany-563/muscle3.8.31_i86darwin64 -in SEQUENCE.fasta -out sequence-aligned.fasta 

MUSCLE v3.8.31 by Robert C. Edgar

http://www.drive5.com/muscle
This software is donated to the public domain.
Please cite: Edgar, R.C. Nucleic Acids Res 32(5), 1792-97.

SEQUENCE 4 seqs, max length 16490, avg  length 16486
00:00:00      2 MB(0%)  Iter   1  100.00%  K-mer dist pass 1
00:00:00      2 MB(0%)  Iter   1  100.00%  K-mer dist pass 2
00:00:20    334 MB(4%)  Iter   1  100.00%  Align node       
00:00:20    334 MB(4%)  Iter   1  100.00%  Root alignment
00:00:20    334 MB(4%)  Iter   2  100.00%  Root alignment
00:00:53    345 MB(4%)  Iter   3  100.00%  Refine biparts
00:01:28    345 MB(4%)  Iter   4  100.00%  Refine biparts

This gave me an aligned sequence.

## Distance and Parsimony
First, we must install and load packages in terminal to perform these functions.
    install.packages("adegenet", dep=TRUE)
    install.packages("phangorn", dep=TRUE)
    *I did this directly in terminal because when I put it in my rscript code it gave me an error
Distance Based: Put this code into a file named myscript.R
    
   library(ape)
library(adegenet)
library(phangorn)
dna <- fasta2DNAbin("SEQUENCE.fasta")
D <- dist.dna(dna, model="TN93")
tre <- njs(D)
tre <- ladderize(tre)
plot(tre, cex=.6)
title("A simple NJ tree")

    
    
    % Rscript myscript.R      
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

Parsimony Based:
put this code into an rscript file called parsimony.R
   library(ape)
library(adegenet)
library(phangorn)
dna <- fasta2DNAbin(file="SEQUENCE.fasta")
dna2 <- as.phyDat(dna)
tre.ini <- njs(dist.dna(dna,model="raw"))
parsimony(tre.ini, dna2)
tre.pars <- optim.parsimony(tre.ini, dna2)
plot(tre.pars, cex=0.6)
title("plotted tree")

    
   % Rscript parsimony.R
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
Final p-score 23059 after  0 nni operations

## OrthoFinder
    Install OrthoFinder:
conda install -c bioconda orthofinder
Collecting package metadata (current_repodata.json): done
Solving environment: done
    Assure Orthofinder is installed correctly
orthofinder -h 

OrthoFinder version 2.5.4 Copyright (C) 2014 David Emms

SIMPLE USAGE:
Run full OrthoFinder analysis on FASTA format proteomes in <dir>
  orthofinder [options] -f <dir>

Add new species in <dir1> to previous run in <dir2> and run new analysis
  orthofinder [options] -f <dir1> -b <dir2>

OPTIONS:
 -t <int>        Number of parallel sequence search threads [Default = 4]
 -a <int>        Number of parallel analysis threads
 -d              Input is DNA sequences
 -M <txt>        Method for gene tree inference. Options 'dendroblast' & 'msa'
                 [Default = dendroblast]
 -S <txt>        Sequence search program [Default = diamond]
                 Options: blast, diamond, diamond_ultra_sens, blast_gz, mmseqs, blast_nucl
 -A <txt>        MSA program, requires '-M msa' [Default = mafft]
                 Options: mafft, muscle
 -T <txt>        Tree inference method, requires '-M msa' [Default = fasttree]
                 Options: fasttree, raxml, raxml-ng, iqtree
 -s <file>       User-specified rooted species tree
 -I <int>        MCL inflation parameter [Default = 1.5]
 -x <file>       Info for outputting results in OrthoXML format
 -p <dir>        Write the temporary pickle files to <dir>
 -1              Only perform one-way sequence search
 -X              Don't add species names to sequence IDs
 -y              Split paralogous clades below root of a HOG into separate HOGs
 -z              Don't trim MSAs (columns>=90% gap, min. alignment length 500)
 -n <txt>        Name to append to the results directory
 -o <txt>        Non-default results directory
 -h              Print this help text

WORKFLOW STOPPING OPTIONS:
 -op             Stop after preparing input files for BLAST
 -og             Stop after inferring orthogroups
 -os             Stop after writing sequence files for orthogroups
                 (requires '-M msa')
 -oa             Stop after inferring alignments for orthogroups
                 (requires '-M msa')
 -ot             Stop after inferring gene trees for orthogroups 

WORKFLOW RESTART COMMANDS:
 -b  <dir>         Start OrthoFinder from pre-computed BLAST results in <dir>
 -fg <dir>         Start OrthoFinder from pre-computed orthogroups in <dir>
 -ft <dir>         Start OrthoFinder from pre-computed gene trees in <dir>

LICENSE:
 Distributed under the GNU General Public License (GPLv3). See License.md

CITATION:
 When publishing work that uses OrthoFinder please cite:
 Emms D.M. & Kelly S. (2019), Genome Biology 20:238

 If you use the species tree in your work then please also cite:
 Emms D.M. & Kelly S. (2017), MBE 34(12): 3267-3278
 Emms D.M. & Kelly S. (2018), bioRxiv https://doi.org/10.1101/267914
 
    Running OrthoFinder:

 
## RAxML-NG
Checking the Version:
addisongibbs@AG raxml-ng_v1.1.0_macos_x86_64 % ./raxml-ng -v

RAxML-NG v. 1.1.0 released on 29.11.2021 by The Exelixis Lab.
Developed by: Alexey M. Kozlov and Alexandros Stamatakis.
Contributors: Diego Darriba, Tomas Flouri, Benoit Morel, Sarah Lutteropp, Ben Bettisworth.
Latest version: https://github.com/amkozlov/raxml-ng
Questions/problems/suggestions? Please visit: https://groups.google.com/forum/#!forum/raxml

System: Intel(R) Core(TM) i5-8210Y CPU @ 1.60GHz, 2 cores, 8 GB RAM
## MrBayes
## IQTree
## Coalescent (ASTRAL)
Some code:
java -jar astral.5.7.8.jar -i test_data/song_mammals.424.gene.tre

================== ASTRAL ===================== 

This is ASTRAL version 5.7.8
Gene trees are treated as unrooted
424 trees read from test_data/song_mammals.424.gene.tre
index0
All output trees will be *arbitrarily* rooted at Chicken

======== Running the main analysis
Number of taxa: 37 (37 species)
Taxa: [Chicken, Marmoset, Orangutan, Human, Chimpanzee, Gorilla, Macaque, Galagos, Mouse_Lemur, Tree_Shrew, Mouse, Rat, Kangaroo_Rat, Guinea_Pig, Squirrel, Tarsier, Rabbit, Pika, Microbat, Megabat, Horse, Dolphin, Cow, Alpaca, Pig, Dog, Cat, Shrew, Hedgehog, Lesser_Hedgehog_Tenrec, Hyrax, Elephant, Sloth, Armadillos, Platypus, Opossum, Wallaby]
Taxon occupancy: {Rat=424, Tarsier=424, Dolphin=424, Rabbit=424, Macaque=424, Pika=424, Alpaca=424, Shrew=424, Sloth=424, Tree_Shrew=424, Kangaroo_Rat=424, Armadillos=424, Chimpanzee=424, Horse=424, Dog=424, Human=424, Lesser_Hedgehog_Tenrec=424, Microbat=424, Platypus=424, Wallaby=424, Cow=424, Pig=424, Marmoset=424, Megabat=424, Hedgehog=424, Mouse=424, Guinea_Pig=424, Mouse_Lemur=424, Cat=424, Hyrax=424, Elephant=424, Chicken=424, Orangutan=424, Opossum=424, Galagos=424, Squirrel=424, Gorilla=424}
Number of gene trees: 424
0 trees have missing taxa
Calculating quartet distance matrix (for completion of X)
Species tree distances calculated ...
Building set of clusters (X) from gene trees 
------------------------------
gradient0: 1933
Number of Clusters after addition by distance: 1933
calculating extra bipartitions to be added at level 1 ...
Adding to X using resolutions of greedy consensus ...
Limit for sigma of degrees:975
polytomy size limit : 4
discarded polytomies:  [3, 3, 4, 4]
Threshold 0.0:
Threshold 0.01:
Threshold 0.02:
Threshold 0.05:
Threshold 0.1:
Threshold 0.2:
polytomy of size 4; rounds with additions with at least 5 support: 0; clusters: 1933
Threshold 0.3333333333333333:
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 1933
polytomy of size 3; rounds with additions with at least 5 support: 0; clusters: 1933
polytomy of size 4; rounds with additions with at least 5 support: 0; clusters: 1933
max k is :0
Number of Clusters after addition by greedy: 1933
gradient0 in heuristiic: 1933
partitions formed in 1.2 secs
Dynamic Programming starting after 1.2 secs
Using tree-based weight calculation.
Using polytree-based weight calculation.
Polytree max score: 28003080
Polytree building time: 0.224 seconds.
Number of quartet trees in the gene trees: 28003080
Size of largest cluster: 37
Greedy score: 24862814
estimationFactor: 1.1263037241078182
Sub-optimal score: 25489533
Total Number of elements weighted: 3372
Normalized score (portion of input quartet trees satisfied before correcting for multiple individuals): 0.9115752624354179
Optimization score: 25526915
Optimal tree inferred in 2.299 secs.
(Chimpanzee,(Human,(Gorilla,(Orangutan,(Macaque,(Marmoset,(Tarsier,((Galagos,Mouse_Lemur),((Tree_Shrew,((Rabbit,Pika),(Squirrel,(Guinea_Pig,(Kangaroo_Rat,(Mouse,Rat)))))),((((Chicken,Platypus),(Opossum,Wallaby)),((Sloth,Armadillos),(Lesser_Hedgehog_Tenrec,(Hyrax,Elephant)))),((Shrew,Hedgehog),((Microbat,Megabat),((Horse,(Dog,Cat)),(Alpaca,(Pig,(Dolphin,Cow))))))))))))))));
Final quartet score is: 25526915
Final normalized quartet score is: 0.9115752624354179
Extended species tree:
(Chicken,(Platypus,((Opossum,Wallaby)1:4.9534768802562965,(((Sloth,Armadillos)1:4.3332364705044455,(Lesser_Hedgehog_Tenrec,(Hyrax,Elephant)1:1.5072311535707152)1:3.0324267685203896)1:0.11752571582479492,(((Shrew,Hedgehog)1:1.0128082767511968,((Microbat,Megabat)1:1.5529600021930556,((Horse,(Dog,Cat)1:2.9057840368910477)0.9:0.0641150813890718,(Alpaca,(Pig,(Dolphin,Cow)1:1.3549566469016843)1:0.6392263251934307)1:3.5727535641858488)1:0.11521870556469156)1:0.43083761402314513)1:1.9714179086025256,((Tree_Shrew,((Rabbit,Pika)1:3.449399483480025,(Squirrel,(Guinea_Pig,(Kangaroo_Rat,(Mouse,Rat)1:5.045850200387328)1:0.6153942169908233)1:0.14915704136865612)1:1.6583346999346553)1:0.843253312273408)0.91:0.08610422555073367,((Galagos,Mouse_Lemur)1:2.4173938646853443,(Tarsier,(Marmoset,(Macaque,(Orangutan,(Gorilla,(Human,Chimpanzee)1:0.6434676443273446)1:2.3363022618905545)1:2.5809965452562977)1:2.687856981495734)1:4.347341076686)1:0.6511829814609134)1:2.1216694610241857)1:2.3363739622649566)1:1.5631449093243042)1:3.3544501191225256)1:0.9262330413691204));
(Chicken,(Platypus,((Opossum,Wallaby)1:4.9534768802562965,(((Sloth,Armadillos)1:4.3332364705044455,(Lesser_Hedgehog_Tenrec,(Hyrax,Elephant)1:1.5072311535707152)1:3.0324267685203896)1:0.11752571582479492,(((Shrew,Hedgehog)1:1.0128082767511968,((Microbat,Megabat)1:1.5529600021930556,((Horse,(Dog,Cat)1:2.9057840368910477)0.9:0.0641150813890718,(Alpaca,(Pig,(Dolphin,Cow)1:1.3549566469016843)1:0.6392263251934307)1:3.5727535641858488)1:0.11521870556469156)1:0.43083761402314513)1:1.9714179086025256,((Tree_Shrew,((Rabbit,Pika)1:3.449399483480025,(Squirrel,(Guinea_Pig,(Kangaroo_Rat,(Mouse,Rat)1:5.045850200387328)1:0.6153942169908233)1:0.14915704136865612)1:1.6583346999346553)1:0.843253312273408)0.91:0.08610422555073367,((Galagos,Mouse_Lemur)1:2.4173938646853443,(Tarsier,(Marmoset,(Macaque,(Orangutan,(Gorilla,(Human,Chimpanzee)1:0.6434676443273446)1:2.3363022618905545)1:2.5809965452562977)1:2.687856981495734)1:4.347341076686)1:0.6511829814609134)1:2.1216694610241857)1:2.3363739622649566)1:1.5631449093243042)1:3.3544501191225256)1:0.9262330413691204):0.0); 
Weight calculation took 0.68930444 secs
ASTRAL finished in 2.933 secs
