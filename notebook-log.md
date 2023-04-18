## My Data
I chose ">MT150995.1 Perca flavescens voucher none44 cytochrome c oxidase subunit I (COX1) gene, partial cds; mitochondrial"
I chose this because my lab is on fish and I have dealt with yellow perch in the lab.
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
First, we must install and load packages in R to perform these functions.
    install.packages("adegenet", dep=TRUE)
    install.packages("phangorn", dep=TRUE)
    library(ape)
    library(adegenet)
    library(phangorn)
Distance Based: Put this code into a file named myscript.R
    
    dna <- fasta2DNAbin("sequence-aligned-muscle.fasta")
*this took a very long time*
    D <- dist.dna(dna, model="TN93")
    tre <- nj(D)
    tre <- ladderize(tre)
    plot(tre, cex=.6)
    title("A simple NJ tree")
    "Rscript myscript.R" in terminal
    
    to save it in the background:
    screen -S distance
        Rscript myscript.R
        control + A + D - detaches screen
    screen -r distance - brings it back

Parsimony Based:
    dna2 <- as.phyDat(dna)
    tre.ini <- nj(dist.dna(dna,model="raw"))
    parsimony(tre.ini, dna2)
    tre.pars <- optim.parsimony(tre.ini, dna2)
    plot(tre.pars, cex=0.6)
    

