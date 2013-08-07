### R code from vignette source 'refGenome.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: refGenome.Rnw:36-37
###################################################
options(width=60)


###################################################
### code chunk number 2: refGenome.Rnw:51-54
###################################################
library(refGenome)
beg<-ensemblGenome()
basedir(beg)<-system.file("extdata",package="refGenome")


###################################################
### code chunk number 3: refGenome.Rnw:83-86
###################################################
ens_gtf<-"hs.ensembl.62.small.gtf"
read.gtf(beg,ens_gtf)
beg


###################################################
### code chunk number 4: refGenome.Rnw:94-96
###################################################
tableAttributeTypes(beg)
moveAttributes(beg,c("gene_name","transcript_name","exon_number"))


###################################################
### code chunk number 5: refGenome.Rnw:127-133 (eval = FALSE)
###################################################
## uc<-ucscGenome()
## basedir(uc)<-"/my/ucsc/basedir"
## read.gtf(uc,"ucsc_knownGene.gtf")
## addXref(uc,"kgXref.csv")
## addEnsembl(uc,"knownToEnsembl.csv")
## addIsoforms(uc,"ucsc_knownisoforms.csv")


###################################################
### code chunk number 6: refGenome.Rnw:140-144
###################################################
ucfile<-system.file("extdata", "hs.ucsc.small.RData", package="refGenome")
uc<-loadGenome(ucfile)
ensfile<-system.file("extdata", "hs.ensembl.62.small.RData", package="refGenome")
ens<-loadGenome(ensfile)


###################################################
### code chunk number 7: refGenome.Rnw:157-158
###################################################
tableSeqids(ens)


###################################################
### code chunk number 8: refGenome.Rnw:163-165
###################################################
en1<-extractSeqids(ens,"^1$")
en1


###################################################
### code chunk number 9: refGenome.Rnw:174-176
###################################################
ensPrimAssembly()
ucPrimAssembly()


###################################################
### code chunk number 10: refGenome.Rnw:180-184
###################################################
enpa<-extractSeqids(ens,ensPrimAssembly())
tableSeqids(enpa)
ucpa<-extractSeqids(uc,ucPrimAssembly())
tableSeqids(ucpa)


###################################################
### code chunk number 11: refGenome.Rnw:190-193
###################################################
tableFeatures(enpa)
enpf<-extractFeature(enpa,"exon")
enpf


###################################################
### code chunk number 12: refGenome.Rnw:202-204
###################################################
dxe<-extractByGeneName(enpa,"DDX11L1")
dxu<-extractByGeneName(ucpa,"DDX11L1")


###################################################
### code chunk number 13: refGenome.Rnw:211-213
###################################################
tableTranscript.id(enpa)
tableTranscript.id(ucpa)


###################################################
### code chunk number 14: refGenome.Rnw:217-219
###################################################
extractTranscript(ens,"ENST00000456328")
extractTranscript(uc,"uc010nxr.1")


###################################################
### code chunk number 15: refGenome.Rnw:228-232
###################################################
gpe<-getGenePositions(ens)
gpe
gpu<-getGenePositions(ucpa)
gpu


###################################################
### code chunk number 16: refGenome.Rnw:244-246
###################################################
enex<-refExons(ens)
ucex<-refExons(uc)


###################################################
### code chunk number 17: refGenome.Rnw:256-260
###################################################
jens<-getSpliceTable(ens)
jens
juc<-getSpliceTable(uc)
juc


###################################################
### code chunk number 18: refGenome.Rnw:267-271
###################################################
ujens<-unifyJuncs(jens)
ujuc<-unifyJuncs(juc)
jeg<-getGenePositions(jens)
jug<-getGenePositions(juc)


###################################################
### code chunk number 19: refGenome.Rnw:285-294
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


