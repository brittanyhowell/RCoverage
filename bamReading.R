setwd("~/Documents/University/Honours_2016/Project/Ranalysis/")

library(Rsamtools)
library(GenomicRanges)

# Read in both ORF L1s
L1sB <- read.table("../Data/L1Location/L1_Mouse_bothorf.txt")
colnames(L1sB) <- c("chr", "start", "end", "strand")

# Read in orf1 only L1s
L1s1 <- read.table("../Data/L1Location/L1_Mouse_orf1only.txt")
colnames(L1s1) <- c("chr", "start", "end", "strand")

# Read in orf2 only L1s
L1s2 <- read.table("../Data/L1Location/L1_Mouse_orf2only.txt")
colnames(L1s2) <- c("chr", "start", "end", "strand")

# Read in no ORF L1s
L1sN <- read.table("../Data/L1Location/L1_Mouse_noorfs.txt")
colnames(L1sN) <- c("chr", "start", "end", "strand")

names(L1s)
L1s <- c(L1sB, L1s1)

# remove all L1s on unplaced scafolds 
L1s <- L1s[grep(pattern = "chr", x = L1s[,"chr"]),]

L1chr <- as.character(unique(L1s$chr))
allRangesL1 <- NULL
for(i in 1:length(L1chr)){
  dat1 <- IRanges(start = L1s[L1s$chr == L1chr[i], "start"], end = L1s[L1s$chr == L1chr[i], "end"], names = L1s[L1s$chr == L1chr[i], "strand"])
  allRangesL1 <- c(allRangesL1, dat1)
}
names(allRangesL1) = L1chr

allRangesL1 <- RangesList(allRangesL1)


p1 <- ScanBamParam(what=scanBamWhat(), which=allRangesL1)

reads <- scanBam(file = "../Genomes/STAR/Output/mm10Aligned.sortedByCoord.out.bam",
                 index = "../Genomes/STAR/Output/mm10Aligned.sortedByCoord.out.bam.bai" ,
                 param = p1)

readNo <-NULL
for(i in 1:length(reads)){
  readNo <- c(readNo, length(reads[[i]]$pos))
}
hist(log10(readNo), breaks =100)


count <- countBam(file = "../Genomes/STAR/Output/mm10Aligned.sortedByCoord.out.bam",
                  index = "../Genomes/STAR/Output/mm10Aligned.sortedByCoord.out.bam.bai" ,
                  param = p1)



hist(count$records, xlim = c(0,500), breaks = 200000)
plot(count$width, count$nucleotides/count$width)

hist(count$nucleotides/count$width, xlim = c(0,1), breaks = 2000)

plot(count$width, count$records, ylim = c(0, 200))

#plot(coverage(IRanges(start = reads$`chr1:100754214-100762297`$pos, width = 101)))


# Reading coverage table

coverage <- read.table(file = "/Users/brittanyhowell/Documents/University/Honours_2016/Project/ReadCoverage/coverage-NoHist-INTERSECTED.out")
colnames(coverage) <- c("chromosome", "start", "stop", "read", "Unique", "Strand", "Depth", "# bases at depth", "length of read", "Fraction of read at depth")

hist(coverage$`length of read`)

coverage <- read.table(file = "/Users/brittanyhowell/Documents/University/Honours_2016/Project/ReadCoverage/coverage-Rev-INTERSECTED.out")
colnames(coverage) <- c("chromosome", "start", "stop", "read", "Unique", "Strand", "Depth", "Number bases at depth", "length of read", "Fraction of read at depth")

hist(coverage$`Number bases at depth`,  xlim = c(0,600), ylim = c(0,200), breaks = 1000, main = 'Bases in L1 elements at non zero depth', xlab = 'Number of bases at depth', ylab = 'Frequency')
median(coverage$`Number bases at depth`)
mean(coverage$Depth)

hist(coverage$Depth,  xlim = c(0,50), ylim = c(0,300), breaks = 10000, main = 'Read coverage depth of L1 bases which overlap reads', xlab = 'Depth of read coverage for covered bases in an L1 element', ylab = 'Frequency of read coverage')

plot(median(coverage$`Fraction of read at depth`), median(coverage$Depth), main = 'Depth of coverage')
