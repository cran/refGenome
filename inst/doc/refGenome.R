### R code from vignette source 'refGenome.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: refGenome.Rnw:36-37
###################################################
options(width=60)


###################################################
### code chunk number 2: refGenome.Rnw:55-58
###################################################
library(refGenome)
beg<-ensemblGenome()
basedir(beg)<-system.file("extdata",package="refGenome")


###################################################
### code chunk number 3: refGenome.Rnw:87-90
###################################################
ens_gtf<-"hs.ensembl.62.small.gtf"
read.gtf(beg,ens_gtf)
beg


###################################################
### code chunk number 4: refGenome.Rnw:98-100
###################################################
tableAttributeTypes(beg)
moveAttributes(beg,c("gene_name","transcript_name","exon_number"))


###################################################
### code chunk number 5: refGenome.Rnw:131-137 (eval = FALSE)
###################################################
## uc<-ucscGenome()
## basedir(uc)<-"/my/ucsc/basedir"
## read.gtf(uc,"ucsc_knownGene.gtf")
## addXref(uc,"kgXref.csv")
## addEnsembl(uc,"knownToEnsembl.csv")
## addIsoforms(uc,"ucsc_knownisoforms.csv")


###################################################
### code chunk number 6: refGenome.Rnw:144-148
###################################################
ucfile<-system.file("extdata", "hs.ucsc.small.RData", package="refGenome")
uc<-loadGenome(ucfile)
ensfile<-system.file("extdata", "hs.ensembl.62.small.RData", package="refGenome")
ens<-loadGenome(ensfile)


###################################################
### code chunk number 7: refGenome.Rnw:161-162
###################################################
tableSeqids(ens)


###################################################
### code chunk number 8: refGenome.Rnw:167-169
###################################################
en1<-extractSeqids(ens,"^1$")
en1


###################################################
### code chunk number 9: refGenome.Rnw:178-180
###################################################
ensPrimAssembly()
ucPrimAssembly()


###################################################
### code chunk number 10: refGenome.Rnw:184-188
###################################################
enpa<-extractSeqids(ens,ensPrimAssembly())
tableSeqids(enpa)
ucpa<-extractSeqids(uc,ucPrimAssembly())
tableSeqids(ucpa)


###################################################
### code chunk number 11: refGenome.Rnw:194-197
###################################################
tableFeatures(enpa)
enpf<-extractFeature(enpa,"exon")
enpf


###################################################
### code chunk number 12: refGenome.Rnw:206-208
###################################################
dxe<-extractByGeneName(enpa,"DDX11L1")
dxu<-extractByGeneName(ucpa,"DDX11L1")


###################################################
### code chunk number 13: refGenome.Rnw:213-215
###################################################
tableTranscript.id(enpa)
tableTranscript.id(ucpa)


###################################################
### code chunk number 14: refGenome.Rnw:219-221
###################################################
extractTranscript(ens,"ENST00000456328")
extractTranscript(uc,"uc010nxr.1")


###################################################
### code chunk number 15: refGenome.Rnw:230-234
###################################################
gpe<-getGenePositions(ens)
gpe
gpu<-getGenePositions(uc)
gpu


###################################################
### code chunk number 16: refGenome.Rnw:246-248
###################################################
enex<-refExons(ens)
ucex<-refExons(uc)


###################################################
### code chunk number 17: refGenome.Rnw:251-252
###################################################
enex


###################################################
### code chunk number 18: refGenome.Rnw:262-266
###################################################
jens<-getSpliceTable(ens)
jens
juc<-getSpliceTable(uc)
juc


###################################################
### code chunk number 19: refGenome.Rnw:273-279
###################################################
ujens<-unifyJuncs(jens)
ujuc<-unifyJuncs(juc)
jeg<-getGenePositions(jens)
jug<-getGenePositions(juc)
head(ujens)
head(jug)


###################################################
### code chunk number 20: refGenome.Rnw:293-302
###################################################
qry<-data.frame(
                  id=1:6,
                  start=c(10,18,61,78,82,110),
                  end=c(15,22,63,87,90,120))
ref<-data.frame(
                  id=1:5,
                  start=c(20,40,60,80,100),
                  end=c(25,45,65,85,105))
overlap(qry,ref)


