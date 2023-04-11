## My Data
I chose ">MT150995.1 Perca flavescens voucher none44 cytochrome c oxidase subunit I (COX1) gene, partial cds; mitochondrial"
I chose this because my lab is on fish and I have dealt with yellow perch in the lab.
I ran it through MUSCLE to get an aligned sequence. I chose MUSCLE mostly because it was the only one I could get to work with my data and on my generation of MacBook.
Using MUSCLE
    in terminal: muscle -align sequence.fasta -output sequence-aligned-muscle.fasta
This gave me an aligned sequence.

## Distance and Parsimony
First, we must install and load packages in R to perform these functions.
    install.packages("adegenet", dep=TRUE)
    install.packages("phangorn", dep=TRUE)
    library(ape)
    library(adegenet)
    library(phangorn)
Distance Based:
    dna <- fasta2DNAbin("sequence-aligned-muscle.fasta")
*this took a very long time*
    D <- dist.dna(dna, model="TN93")
    tre <- nj(D)
    tre <- ladderize(tre)
    plot(tre, cex=.6)
    title("A simple NJ tree")

Parsimony Based:
    dna2 <- as.phyDat(dna)
    tre.ini <- nj(dist.dna(dna,model="raw"))
    parsimony(tre.ini, dna2)
    tre.pars <- optim.parsimony(tre.ini, dna2)
    plot(tre.pars, cex=0.6)
    

