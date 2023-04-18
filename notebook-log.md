## My Data
I chose data on the complete genome of the freshwater bass family.
I chose this because my research is in water toxicology, so I do a lot of work with fish in the lab. I've never gotten to work with bass, but I found some data on one of them and that kind of got me interested.
GenBank had only 4 different species' complete mitochondrial genome, which was a little sad but was okay because I then had something I could still work with.
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

## RAxML-NG

