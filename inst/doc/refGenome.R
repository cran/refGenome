### R code from vignette source 'refGenome.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: refGenome.Rnw:23-24
###################################################
options(width=60)


###################################################
### code chunk number 2: refGenome.Rnw:39-42
###################################################
library(refGenome)
ens<-ensemblGenome()
basedir(ens)<-system.file("extdata",package="refGenome")


###################################################
### code chunk number 3: refGenome.Rnw:72-75
###################################################
ens_gtf<-"hs.ensembl.62.small.gtf"
read.gtf(ens,ens_gtf)
ens


###################################################
### code chunk number 4: refGenome.Rnw:83-85
###################################################
tableAttributeTypes(ens)
moveAttributes(ens,c("gene_name","transcript_name","exon_number"))


###################################################
### code chunk number 5: refGenome.Rnw:117-123 (eval = FALSE)
###################################################
## uc<-ucscGenome()
## basedir(uc)<-"/my/ucsc/basedir"
## read.gtf(uc,"ucsc_knownGene.gtf")
## addXref(uc,"kgXref.csv")
## addEnsembl(uc,"knownToEnsembl.csv")
## addIsoforms(uc,"ucsc_knownisoforms.csv")


###################################################
### code chunk number 6: refGenome.Rnw:130-132
###################################################
ucfile<-system.file("extdata", "hs.ucsc.small.RData", package="refGenome")
uc<-load.ucsc(ucfile)


###################################################
### code chunk number 7: refGenome.Rnw:147-148
###################################################
tableSeqids(ens)


###################################################
### code chunk number 8: refGenome.Rnw:153-155
###################################################
en1<-extractSeqids(ens,"^1$")
en1


###################################################
### code chunk number 9: refGenome.Rnw:164-166
###################################################
ensPrimAssembly()
ucPrimAssembly()


###################################################
### code chunk number 10: refGenome.Rnw:170-174
###################################################
enpa<-extractSeqids(ens,ensPrimAssembly())
tableSeqids(enpa)
ucpa<-extractSeqids(uc,ucPrimAssembly())
tableSeqids(ucpa)


###################################################
### code chunk number 11: refGenome.Rnw:180-183
###################################################
tableFeatures(enpa)
enpf<-extractFeature(enpa,"exon")
enpf


###################################################
### code chunk number 12: refGenome.Rnw:192-194
###################################################
dxe<-extractByGeneName(enpa,"DDX11L1")
dxu<-extractByGeneName(ucpa,"DDX11L1")


###################################################
### code chunk number 13: refGenome.Rnw:198-200
###################################################
tableTranscript.id(dxe)
tableTranscript.id(dxu)


###################################################
### code chunk number 14: refGenome.Rnw:204-206
###################################################
extractTranscript(dxe,"ENST00000456328")
extractTranscript(dxu,"uc010nxr.1")


###################################################
### code chunk number 15: refGenome.Rnw:218-222
###################################################
gpe<-getGenePositions(enpa)
gpe
gpu<-getGenePositions(ucpa)
gpu


###################################################
### code chunk number 16: refGenome.Rnw:234-243
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


