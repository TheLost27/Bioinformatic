# Bioconductor yükleme

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Biostrings")


library(Biostrings)


# İkili Dizi Hizalama
?pairwiseAlignment

pairwiseAlignment(pattern = "succeed", subject = "supersede")


A1<-pairwiseAlignment(pattern = "AGTA", subject = "AACTAACTA")

pairwiseAlignment(pattern = "AACTA", subject = "AACTAACTA")

pairwiseAlignment(pattern = "biyoinformatik", subject = "istatistik",type="local")

A2<-pairwiseAlignment(pattern = "aatgctcgta", subject = "acatctgaa",type="local")

score(A2)

nmatch(A2) # eşleşme sayısı

nmismatch(A2) # eşleşmeme sayısı

nmismatch(A1)

deletion(A1) # boşluk sayısı,pozisyonu

pid(A1) # benzerlik oranı

type(A1) # hizalama türü



?nucleotideSubstitutionMatrix

mat <- nucleotideSubstitutionMatrix(match = 1, mismatch = -3)

globalAlign <-
  pairwiseAlignment(pattern="AATGCTA", subject="ACGTCAAG", substitutionMatrix = mat,
                    gapOpening = 5, gapExtension = 2)
globalAlign

localAlign <-
  pairwiseAlignment(pattern="AATGCTA", subject="ACGTCAAG", type = "local", substitutionMatrix = mat,
                    gapOpening = 5, gapExtension = 2)
localAlign

# GenBank (NCBI) dan veri okuma işlemi

install.packages("ape")
library(ape)

?read.GenBank


Delta <- read.GenBank("OK091006.1", as.character = "TRUE")
Omicron <- read.GenBank("OM287553.1", as.character = "TRUE")


# GenBank'dan veri çekip ikili dizi hizalama yapalım.

BRCA1_musmusculus <- read.GenBank("EU349657.1", as.character = "TRUE")
BRCA1_homosapiens <- read.GenBank("MF945608.1", as.character = "TRUE")

seq1 = paste(BRCA1_musmusculus, collapse="")
seq2 = paste(BRCA1_homosapiens, collapse="")
pairwiseAlignment(pattern=seq1, subject=seq2)

al<-pairwiseAlignment(pattern=seq1, subject=seq2, type="local")
al

pairwiseAlignment(pattern=seq1, subject=seq2, type="local",scoreOnly=TRUE)

# Aminoasit dizilerini hizalama

data(package="Biostrings")
data("BLOSUM50")
data("PAM250")
BLOSUM50[1:10,1:10]
AAalign<-pairwiseAlignment(AAString("PAWHEAE"), AAString("HEAGAWGHEE"), 
                           substitutionMatrix = BLOSUM50, gapOpening = 0, 
                           gapExtension = 8)

AAalign

summary(AAalign)


insulin_hs <-"MALWMRLLPLLALLALWAPAPTRAFVNQHLCGSHLVEALYLVCGERGFFYTPKARREVEDLQVRDVELAGAPGEGGLQPLALEGALQKRGIVEQCCTSICSLYQLENYCN"

insulin_cl <-"MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN"


ins_align <-pairwiseAlignment(AAString(insulin_hs),AAString(insulin_cl), substitutionMatrix = PAM250, gapOpening = 0, gapExtension = 8)
ins_align

summary(ins_align)



################ .fasta uzantılı dosyadaki veriyi kullanma ################

#Nükleotit verisi

library(seqinr)

?read.fasta

prokaryotes <- read.fasta(file = "prok.fasta", seqtype = "DNA")
seq1 <-as.character(prokaryotes[[1]]) 
seq1 = paste(seq1, collapse="")
seq2 <-as.character(prokaryotes[[2]])
seq2 = paste(seq2, collapse="")

pairalign<-pairwiseAlignment(pattern=seq2, subject=seq1)
summary(pairalign)

# Aminoasit verisi

coxgenes <- read.fasta(file = "cox1multi.fasta", seqtype="AA")


# Nokta Matrisi çizimi
?dotPlot

cox1<-coxgenes[[1]]
cox2<-coxgenes[[2]]

dotPlot(cox1, cox2, main = "Human vs Mouse Cox1 Dotplot")

dotPlot(cox1, cox2, wsize = 3, wstep = 3, nmatch = 3, main = "Human vs Mouse Cox1 Dotplot\nwsize = 3, wstep = 3, nmatch = 3")

dotPlot(cox1[1:100], cox2[1:100], wsize = 3, wstep = 3, nmatch = 3, main = "Human vs Mouse Cox1 first 100 AA Dotplot\nwsize = 3, wstep = 3, nmatch = 3")


