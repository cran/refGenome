
##################################################################################################
##                                                                                              ##
##  Project   :   refGenome                                                                     ##
##  Created   :   19.03.2012                                                                    ##
##  Author    :   W. Kaisers                                                                    ##
##  Content   :   Funktionality for importing and managing genomic reference data               ##
##                (Ucsc, Ensembl, Genbank)                                                      ##
##                for usage in R                                                                ##
##  Version   :   1.0.0                                                                         ##
##                                                                                              ##
##  Changelog :                                                                                 ##
##  05.Jun.12 :   addIsoforms and addEnsemble delete qualifier rows before (re-) inserting      ##
##  06.Jun.12 :   Added correction in get_ens_attribute_df which removed a memory leak          ##
##                (Also the token_list module is now valgrind tested)                           ##
##  08.Jun.12 :   Added class refdb             (encapsulates database access)                  ##
##  09.Jun.12 :   Added class refDataLocations  (encapuslates directory and file management)    ##
##  13.Jun.12 :   Implementation updates for refGenome,ensembl,ucsc and genbank finished        ##
##  26.Nov.12 :   Added function 'extractSeqids' and 'tableSeqids'                              ##
##  09.Mar.13 :   Added class and function ensemblExons
##                                                                                              ##
##################################################################################################

.onUnload<-function(libpath) { library.dynam.unload("refGenome",libpath) }

###################################################################################################
##                                                                                               ##
## refGenome      Functionality for processing Sequence and Annotation data in existing          ##
##                File (and dir) structures.                                                     ##
##                                                                                               ##
##                refdir/(input data)                                                            ##
##                refdir/seqs                                                                    ##
##                refdir/dbname.db3                                                              ##
##                                                                                               ##
##                                                                                               ##
###################################################################################################

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
#                                                                                                  #
# Split gtf Attribute column data                                                                  #
# and return list with two data.frames                                                             #
#                                                                                                  #
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #

# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
# A: The Brent Lab (Washington University, St.Louis)                                               #
# http:#mblab.wustl.edu/GTF2.html                                                                  #
# [attributes] All four features have the same two mandatory attributes at the end of the record:  #
#                                                                                                  #
# gene_id value;       A globally unique identifier for the genomic source of the transcript       #
# transcript_id value; A globally unique identifier for the predicted transcript.                  #
#                                                                                                  #
# These attributes are designed for handling multiple transcripts from the same genomic region.    #
# Any other attributes or comments must appear after these two and will be ignored.                #
#                                                                                                  #
# Attributes must end in a semicolon which must then be separated from the start of any subsequent #
# attribute by exactly one space character (NOT a tab character). Textual attributes *should* be   #
# surrounded by doublequotes.                                                                      #
#                                                                                                  #
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
# B: Wellcome Trust Sanger Institute                                                               #
# http:#www.sanger.ac.uk/resources/software/gff/spec.html                                          #
#                                                                                                  #
# Free text values *must* be quoted with double quotes.                                            #
# Note: all non-printing characters in such free text value strings (e.g. newlines, tabs, control  #
# characters, etc) must be explicitly represented by their C (UNIX) style backslash-escaped        #
# representation (e.g. newlines as '\n', tabs as '\t').                                            #
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
setGeneric("extractByGeneName",function(object,geneNames)standardGeneric("extractByGeneName"))
setGeneric("getGenePositions",function(object)standardGeneric("getGenePositions"))

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
#                                                                                                   #
#  refGenome:                                                                                       #
#   Virtual base class for ucscGenome and ensemblGenome                                             #
#                                                                                                   #
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #

# con can also be NULL
setClass("refGenome",
            representation(
                  ev="environment",
                  basedir="character",
                  "VIRTUAL"),
            prototype=prototype(
              ev=new.env(),
              basedir="."
              )
         )

setMethod("show","refGenome",function(object)
{
  if(!exists("gtf",where=object@ev,inherits=FALSE))
    cat("An empty object of class '",class(object),"'.\n",sep="")
  else
  {
    n<-min(dim(object@ev$gtf)[1],6L)
    bm<-Sys.localeconv()[7]
    cat("Object of class '",class(object),"' with ",format(dim(object@ev$gtf)[1],big.mark=bm)," rows and ",dim(object@ev$gtf)[2]," columns.\n",sep="")
    print(object@ev$gtf[1:n,])    
  }
})

setGeneric("basedir",function(object)standardGeneric("basedir"))
setMethod("basedir","refGenome",function(object) {return(object@basedir)})

setGeneric("basedir<-",function(object,value)standardGeneric("basedir<-"))
setReplaceMethod("basedir","refGenome",function(object,value)
{
  if(!file.exists(value))
    cat("[basedir.refGenome] Directory '",value,"' does not exist!\n",sep="")
  object@basedir<-value
  return(object)
})



setGeneric("getGtf",function(object)standardGeneric("getGtf"))
setMethod("getGtf","refGenome",function(object)
{
  if(!exists("gtf",where=object@ev,inherits=FALSE))  
    return(NULL)
  return(object@ev$gtf)
})

setGeneric("setGtf",function(object,value)standardGeneric("setGtf"))
setMethod("setGtf","refGenome",function(object,value)
{
  assign("gtf",value,envir=object@ev)
  return(invisible())
})

setGeneric("getAttr",function(object)standardGeneric("getAttr"))
setMethod("getAttr","refGenome",function(object)
{
 if(!exists("gtfattributes",where=object@ev,inherits=FALSE))
   return(NULL)
 return(object@ev$gtfattributes)
})

setGeneric("setAttr",function(object,value)standardGeneric("setAttr"))
setMethod("setAttr","refGenome",function(object,value)
{
  assign("gtfattributes",value,envir=object@ev)
  return(invisible())
})

setGeneric("tableAttributeTypes",function(object)standarGeneric("tableAttributeTypes"))
setMethod("tableAttributeTypes","refGenome",function(object)
{
  if(!exists("gtfattributes",where=object@ev,inherits=FALSE))
  {
    cat("[tableAttributeTypes.refGenome] gtfattributes table does not exist.\n")
  }
  else
  {   
    bm<-Sys.localeconv()[7]
    cat("[tableAttributeTypes.refGenome] Row number in gtf-table: ",format(dim(object@ev$gtf)[1],big.mark=bm),".\n",sep="")
    return(table(object@ev$gtfattributes$type))
  }
})


setGeneric("read.gtf",function(object,filename="transcripts.gtf",sep="\t",useBasedir=TRUE,...)standardGeneric("read.gtf"))
setMethod("read.gtf","refGenome",function(object,filename="transcripts.gtf",sep="\t",useBasedir=TRUE,...)
{
  if(useBasedir && (length(object@basedir)>0))
    filename<-file.path(object@basedir,filename)
  if(!file.exists(filename))
    stop("[read.gtf.refGenome] file '",filename,"' does not exist!\n",sep="")
  cat("[read.gtf.refGenome] Reading file '",basename(filename),"'.\n",sep="")
  
  # Import gtf via read.table
  tbl<-read.table(filename,sep=sep,quote="",comment.char="#")
  if(dim(tbl)[2]!=9)
    stop("[read.gtf.refGenome] Wrong number of columns: ",dim(tbl)[2],"\n")
  names(tbl)<-c("seqid","source","feature","start","end","score","strand","frame","attributes")
  n<-dim(tbl)[1]
  tbl$id<-as.integer(1:n)
  # (for seqid class factor is adequate)
  tbl$attributes<-as.character(tbl$attributes)
  
  cat("[read.gtf.refGenome] Parsing attributes.\n")
  l<-.Call("split_gtf_attr",tbl$id,tbl$attributes,PACKAGE="refGenome")
  tbl$gene_id<-l$fixed$gene_id
  tbl$transcript_id<-l$fixed$transcript_id
  if(dim(l$variable)[1]>0)
  {  
    gtfattributes<-l$variable
    assign("gtfattributes",gtfattributes,envir=object@ev)    
  }
  else
    gtfattributes<-NULL
  # Remove text after successful parsing
  tbl$attributes<-NULL
  # sort tbl
  tbl<-tbl[order(tbl$seqid,tbl$start),]  
  assign("gtf",tbl,envir=object@ev)

  bm<-Sys.localeconv()[7]
  if(!is.null(gtfattributes))
    cat("[read.gtf.refGenome] Finished reading",format(dim(tbl)[1],big.mark=bm),"gtf lines and",format(dim(gtfattributes)[1],big.mark=bm),"gtfattributes lines.\n")
  else
    cat("[read.gtf.refGenome] Finished reading",format(dim(tbl)[1],big.mark=bm),"gtf lines (no attributes).\n")
  return(invisible())
})


setGeneric("saveGenome",function(object,filename,...)standardGeneric("saveGenome"))
setMethod("saveGenome","refGenome",function(object,filename,useBasedir=FALSE,...)
{
  if(!is.character(filename))
    stop("[saveGenome.refGenome] filename must be character!")
  
  if(useBasedir && length(object@basedir)>0)
    save(file=file.path(object@basedir,filename),list=ls(envir=object@ev),envir=object@ev,...)
  else
    save(file=filename,list=ls(envir=object@ev),envir=object@ev,...)
  return(invisible())
})
       
setGeneric("writeDB",function(object,filename,useBasedir=TRUE,...)standardGeneric("writeDB"))
setMethod("writeDB","refGenome",function(object,filename,useBasedir=TRUE,...)
{
  bm<-Sys.localeconv()[7]
  
  if(missing(filename))
    stop("filename argument is not optional ")     
  
  require("RSQLite")
  if(useBasedir && (length(object@basedir)>0))
    filename<-file.path(object@basedir,filename)  
  
  drv<-dbDriver("SQLite")
  con<-dbConnect(drv,filename)
  
  # No open database connection
  if(exists("gtf",where=object@ev,inherits=FALSE))
  {
    dbWriteTable(con,"gtf",object@ev$gtf,append=FALSE,overwrite=TRUE)
    cat("[writeDB.refGenome]",format(dim(object@ev$gtf)[1],big.mark=bm),"rows written to table 'gtf'.\n")
  }
  if(exists("gtfattributes",where=object@ev,inherits=FALSE))
  {
    dbWriteTable(con,"gtfattributes",object@ev$gtfattributes,append=FALSE,overwrite=TRUE)
    cat("[writeDB.refGenome]",format(dim(object@ev$gtfattributes)[1],big.mark=bm),"rows written to table 'gtfattributes'.\n")  
  }
  if(exists("xref",where=object@ev,inherits=FALSE))
  {
    if(dim(object@ev$xref)[1]>0)
    {
      dbWriteTable(con,"xref",object@ev$xref,append=FALSE,overwrite=TRUE)
      cat("[writeDB.refGenome]",format(dim(object@ev$xref)[1],big.mark=bm),"rows written to table 'xref'.\n")
    }
  }
  dbDisconnect(con)
  return(invisible())
})


# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
#                                                                                                   #
#  ucscGenome                                                                                       #
#                                                                                                   #
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #

setClass("ucscGenome",contains="refGenome")

ucscGenome<-function(dbfile)
{
  uc<-new("ucscGenome",ev=new.env())
  if(missing(dbfile))
    return(uc)
  else
  {
    drv<-dbDriver("SQLite")
    uc@basedir<-dirname(dbfile)
  }
  return(uc)
}
           
load.ucsc<-function(filename)           
{
   uc<-new("ucscGenome",ev=new.env(),basedir=dirname(filename))
   load(filename,envir=uc@ev)
   return(invisible(uc))
}

setGeneric("addIsoforms",function(object,filename="ucsc_knownisoforms.csv",sep="\t",useBasedir=TRUE,...)standardGeneric("addIsoforms"))
setMethod("addIsoforms","ucscGenome",function(object,filename="ucsc_knownisoforms.csv",sep="\t",useBasedir=TRUE,...)
{
  if(useBasedir && (length(object@basedir)>0))
    filename<-file.path(object@basedir,filename)
  if(!exists("gtf",where=object@ev,inherits=FALSE))
    stop("[addIsoforms.ucscGenome] No gtf table present! Use 'read.gtf' or 'setGtf'!\n")
  if(!file.exists(filename))
    stop("[addIsoforms.ucscGenome] file '",filename,"' does not exist!\n",sep="")
  cat("[addIsoforms.ucscGenome] Reading file '",basename(filename),"'.\n",sep="")  
  dt<-read.table(filename,sep=sep,quote="",stringsAsFactors=FALSE,comment.char="#")
  if(dim(dt)[2]!=2)
    stop("[addIsoforms.ucscGenome] Wrong number of columns: ",dim(df)[2],"\n")
  names(dt)<-c("clusterId","transcript_id")

  mtc<-match(object@ev$gtf$transcript_id,dt$transcript_id)
  if(any(is.na(mtc)))
    cat("[addIsoforms.ucscGenome] Warning: gtf$transcript_id misses matches in isoforms table!\n")
  object@ev$gtf$clusterId<-dt$clusterId[mtc]            
  cat("[addIsoforms.ucscGenome] Finished.\n")
  return(invisible())
})


setGeneric("addEnsembl",function(object,filename="ucsc_knownToEnsembl.csv",sep="\t",useBasedir=TRUE,...)standardGeneric("addEnsembl"))
setMethod("addEnsembl","ucscGenome",function(object,filename="ucsc_knownToEnsembl.csv",sep="\t",useBasedir=TRUE,...)
{ 
  if(useBasedir && (length(object@basedir)>0))
    filename<-file.path(object@basedir,filename)
  if(!exists("gtf",where=object@ev,inherits=FALSE))
    stop("[addEnsembl.ucscGenome] No gtf table present! Use 'read.gtf' or 'setGtf'!\n")
  if(!file.exists(filename))
    stop("[addEnsembl.ucscGenome] file '",filename,"' does not exist!\n",sep="")  
  cat("[addEnsembl.ucscGenome] Reading file '",basename(filename),"'.\n",sep="")
  dt<-read.table(filename,sep=sep,quote="",stringsAsFactors=FALSE,comment.char="#")
  if(dim(dt)[2]!=2)
    stop("[addEnsembl.ucscGenome] Wrong number of columns: ",dim(df)[2],"\n")
  names(dt)<-c("transcript_id","ensembl")
  mtc<-match(object@ev$gtf$transcript_id,dt$transcript_id)
  object@ev$gtf$ensembl<-dt$ensembl[mtc]
  return(invisible())
})


setGeneric("addXref",function(object,filename="kgXref.csv",sep="\t",useBasedir=TRUE,...)standardGeneric("addXref"))
setMethod("addXref","ucscGenome",function(object,filename="kgXref.csv",sep="\t",useBasedir=TRUE,...)
{
  if(useBasedir && (length(object@basedir)>0))
    filename<-file.path(object@basedir,filename)
  if(!exists("gtf",where=object@ev,inherits=FALSE))
    stop("[addXref.ucscGenome] No gtf table present! Use 'read.gtf' or 'setGtf'!\n")
  if(!file.exists(filename))
    stop("[addXref.ucscGenome] file '",filename,"' does not exist!\n",sep="")
  cat("[addXref.ucscGenome] Reading file '",basename(filename),"'.\n",sep="")  
  dt<-read.table(filename,sep=sep,quote="",stringsAsFactors=FALSE,comment.char="#")
  if(dim(dt)[2]!=10)
    stop("[addXref.ucscGenome] Wrong number of columns: ",dim(df)[2],"\n")
  names(dt)<-c("kgID","mRNA","spID","spDisplayID","gene_name","refseq","protAcc","description","rfamAcc","tRnaName")
  mtc<-match(object@ev$gtf$gene_id,dt$kgID)
  object@ev$gtf$gene_name<-dt$gene_name[mtc]
  assign("xref",dt,envir=object@ev)
  return(invisible())
})

setGeneric("getXref",function(object)standardGeneric("getXref"))
setMethod("getXref","ucscGenome",function(object)
{
  if(!exists("xref",where=object@ev,inherits=FALSE))
    stop("[getXref.ucscGenome] No xref table present! Use 'addXref'!\n")
  return(invisible(object@ev$xref))
})

load.ucsc.db<-function(filename)
{
  if(!file.exists(filename))
    stop("File '",basename(filename),"' does not exist")
  dbdrv<-dbDriver("SQLite")
  dbcon<-dbConnect(dbdrv,filename)
  bm<-Sys.localeconv()[7]
  
  uc<-new("ucscGenome",basedir=dirname(filename))

  copy_table<-function(con,tablename)
  {
    if(!dbExistsTable(con,tablename))
    {
      cat("[loadUcsc.db] Table '",tablename,"' does not exist.\n",sep="")
    }
    else
    {
      assign(tablename,dbReadTable(con,tablename),envir=uc@ev)
      cat("[loadUcsc.db] ",format(dim(uc@ev$gtf)[1],big.mark=bm)," rows written to table '",tablename,"'.\n",sep="")
    }
  }  
  copy_table(dbcon,"gtf")
  copy_table(dbcon,"gtfattributes")
  copy_table(dbcon,"xref")  
  return(uc)  
}

setMethod("extractByGeneName","ucscGenome",function(object,geneNames)
{
  if(!exists("gtf",where=object@ev,inherits=FALSE))
    stop("[extractByGeneName.ucscGenome] gtf table does not exist! Use 'read.gtf'.")
  if(is.na(match("gene_name",names(object@ev$gtf))))
    stop("[extractByGeneName.ucscGenome] 'gene_name' column does not exist in 'gtf' table!")
  if(!is.character(geneNames))
    stop("[extractByGeneName.ucscGenome] geneNames must be character!")
  
  mtc<-match(geneNames,object@ev$gtf$gene_name)
  if(any(is.na(mtc)))
  {
    cat("[extractByGeneName.ucscGenome] Missing matches for ",sum(is.na(mtc))," gene_name(s):\n",sep="")
    print(geneNames[is.na(mtc)])
  }
  mtc<-mtc[!is.na(mtc)]
  # Retrieving gene_id's for geneNames
  dtb<-data.frame(gene_name=object@ev$gtf$gene_name[mtc])
  # Returning all rows that match with found gene_name's
  ev<-new.env()
  gtf<-merge(object@ev$gtf,dtb)
  gtf<-gtf[order(gtf$gene_name,gtf$start),]
  #c("id","seqid","start","end","strand","frame","feature","source","gene_id","transcript_id")]
  assign("gtf",gtf[order(gtf$seqid,gtf$start),],envir=ev)
  
  return(new("ucscGenome",ev=ev,basedir=object@basedir))
})


#setGeneric("getGenePositions",function(object)standardGeneric("getGenePositions"))
setMethod("getGenePositions","ucscGenome",function(object)
{
  if(!exists("gtf",where=object@ev,inherits=FALSE))
    return(NULL)
  if(!is.data.frame(object@ev$gtf))
    stop("[getGenePositions.refGenome] gtf-table must be data.frame!")
  first<-function(x) return(x[1])
  dtb<-object@ev$gtf
  gene<-aggregate(dtb[,c("seqid","source","gene_name","strand")],by=list(gene_id=dtb$gene_id),first)
  gstart<-aggregate(data.frame(start=dtb$start),by=list(gene_id=dtb$gene_id),min)
  gend<-aggregate(data.frame(end=dtb$end),by=list(gene_id=dtb$gene_id),max)
  res<-merge(gene,gstart,by="gene_id")
  res<-merge(res,gend,by="gene_id")
  res<-res[order(res$seqid,res$start),c(2,6,7,5,1,4,3)]
  return(res)
})

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
#                                                                                                   #
#  ensemblGenome                                                                                    #
#                                                                                                   #
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #

setClass("ensemblGenome",contains="refGenome")

ensemblGenome<-function(basedir)
{
  ens<-new("ensemblGenome",ev=new.env())  
  if(missing(basedir))
    return(ens)
  else
    ens@basedir<-basedir
  return(ens)
}
           
load.ensembl<-function(src,basedir)           
{
  if(is.character(src))
  {
    if(missing(basedir))    
      basedir<-dirname(src)
    ens<-new("ensemblGenome",ev=new.env(),basedir=basedir)
    load(src,envir=ens@ev)
    return(invisible(ens))    
  }
  if(is(src,"url"))
  {
    if(missing(basedir))
      basedir<-getwd()
    ens<-new("ensemblGenome",ev=new.env(),basedir=basedir)
    load(src,envir=ens@ev)
    close(src)
    return(invisible(ens))    
  }
  stop("[load.ensembl] Invalid 'src' argument!")
}           

setGeneric("moveAttributes",function(object,names)standardGeneric("moveAttributes"))
setMethod("moveAttributes","ensemblGenome",function(object,names)
{
  if(!exists("gtf",where=object@ev,inherits=FALSE))
    stop("[moveAttributes.ensemblGenome] No gtf table present! Use 'read.gtf' or 'setGtf'!\n")
  if(!exists("gtfattributes",where=object@ev,inherits=FALSE))
    stop("[moveAttributes.ensemblGenome] No gtfattributes table present!\n") 
  if(!is(names,"character"))
    stop("[moveAttributes.ensemblGenome] names must be character!")

  initSize<-dim(object@ev$gtfattributes)[1]
  bm<-Sys.localeconv()[7]
  n<-length(names)
  for(i in 1:n)
  {
    mvtype<-object@ev$gtfattributes$type==names[i]
    smv<-sum(mvtype)
    if(smv==0)
      cat("[moveAttributes.ensemblGenome] Attribute type '",names[i],"' not found in 'gtfattributes' table!\n",sep="")
    else
    {
      mvatt<-object@ev$gtfattributes[mvtype,]
      cat("[moveAttributes.ensemblGenome] Moving ",format(smv,big.mark=bm,width=5)," attributes of type '",names[i],"' from 'gtfattributes' to 'gtf' table.\n",sep="")
      mtc<-match(object@ev$gtf$id,mvatt$id)
      object@ev$gtf[,names[i]]<-mvatt$value[mtc]
      object@ev$gtfattributes<-object@ev$gtfattributes[!mvtype,]    
    }
  }
  cat("[moveAttributes.ensemblGenome] Finished. Reduced attributes table size from ",format(initSize,big.mark=bm)," to ",format(dim(object@ev$gtfattributes)[1],big.mark=bm)," rows. \n",sep="")
  return(invisible())
})

load.ensembl.db<-function(filename)
{
  if(!file.exists(filename))
    stop("File '",basename(filename),"' does not exist!")
  dbdrv<-dbDriver("SQLite")
  dbcon<-dbConnect(dbdrv,filename)
  ev<-new.env()
  bm<-Sys.localeconv()[7]
  ens<-new("ensemblGenome",basedir=dirname(filename))
  
  copy_table<-function(con,tablename)
  {
    if(!dbExistsTable(con,tablename))
      cat("[load.ensembl.db] Table '",tablename,"' does not exist!\n",sep="")
    else
    {
      assign(tablename,dbReadTable(con,tablename),envir=ens@ev)
      cat("[load.ensembl.db] ",format(dim(ens@ev$gtf)[1],big.mark=bm)," rows loaded to table '",tablename,"'.\n",sep="")
    }
  }  
  copy_table(dbcon,"gtf")
  copy_table(dbcon,"gtfattributes")
  return(ens)  
}

#setGeneric("extractByGeneName",function(object,geneNames)standardGeneric("extractByGeneName"))
setMethod("extractByGeneName","ensemblGenome",function(object,geneNames)
{
  if(!exists("gtf",where=object@ev,inherits=FALSE))
    stop("[extractByGeneName.ensemblGenome] gtf table does not exist! Use 'read.gtf'.")
  if(is.na(match("gene_name",names(object@ev$gtf))))
    stop("[extractByGeneName.ensemblGenome] 'gene_id' column does not exist in 'gtf' table! Use 'moveAttributes'.")
  if(!is.character(geneNames))
    stop("[extractByGeneName.ensemblGenome] geneNames must be character!")
  
  mtc<-match(geneNames,object@ev$gtf$gene_name)
  if(any(is.na(mtc)))
  {
    cat("[extractByGeneName.ensemblGenome] Missing matches for ",sum(is.na(mtc))," gene_name(s):\n",sep="")
    print(geneNames[is.na(mtc)])
  }
  mtc<-mtc[!is.na(mtc)]
  # Retrieving gene_id's for geneNames
  dtb<-data.frame(gene_name=object@ev$gtf$gene_name[mtc])

  # Returning all rows that match with found gene_id's
  ev<-new.env()
  assign("gtf",merge(object@ev$gtf,dtb),envir=ev)
  # reorder?
  return(new("ensemblGenome",ev=ev,basedir=object@basedir))
})

setGeneric("tableFeatures",function(object)standardGeneric("tableFeatures"))
setMethod("tableFeatures","ensemblGenome",function(object)
{return(table(object@ev$gtf$feature))})
setMethod("tableFeatures","ucscGenome",function(object)
{return(table(object@ev$gtf$feature))})

setGeneric("extractFeature",function(object,feature)standardGeneric("extractFeature"))
setMethod("extractFeature","ensemblGenome",function(object,feature="exon")
{
  if(!exists("gtf",where=object@ev,inherits=FALSE))
    return(NULL)
  if(!is.data.frame(object@ev$gtf))
    stop("[extractFeature.ensemblGenome] gtf-table must be data.frame!")
  if(!is.character(feature))
    stop("[extractFeature.ensemblGenome] feature must be character")
  if(length(feature)>1)
    stop("[extractFeature.ensemblGenome] feature must have length 1!")
  # Returning all rows that has "CDS" feature
  ev<-new.env()
  gtf<-object@ev$gtf[object@ev$gtf$feature==feature,]
  gtf<-gtf[order(gtf$seqid,gtf$start),]
  assign("gtf",gtf,envir=ev)  
  if(exists("gtfattributes",where=object@ev,inherits=FALSE))
  {
    # Should only be present in ensemblGenome
    tbl<-data.frame(id=ev$gtf$id)
    assign("gtfattributes",merge(object@ev$gtfattributes,tbl),envir=ev)
  }
  return(new("ensemblGenome",ev=ev,basedir=object@basedir))
})

setGeneric("extractTranscript",function(object,transcripts)standardGeneric("extractTranscript"))
setMethod("extractTranscript","ensemblGenome",function(object,transcripts)
{
  if(!exists("gtf",where=object@ev,inherits=FALSE))
    stop("[extractTranscript.ensemblGenome] gtf table does not exist! Use 'read.gtf'.")
  if(is.na(match("transcript_id",names(object@ev$gtf))))
    stop("[extractTranscript.ensemblGenome] 'transcript_id' column does not exist in 'gtf' table! Use 'moveAttributes'.")
  if(!is.character(transcripts))
    stop("[extractTranscript.ensemblGenome] transcripts must be character!")
  mtc<-match(transcripts,object@ev$gtf$transcript_id)
  if(any(is.na(mtc)))
  {
    cat("[extractTranscript.ensemblGenome] Missing matches for ",sum(is.na(mtc))," transcripts:\n",sep="")
    print(transcripts[is.na(mtc)])
  }
  # reorder (transcrpt_id)
  dtb<-data.frame(transcript_id=transcripts)

  ev<-new.env()
  assign("gtf",merge(object@ev$gtf,dtb),envir=ev)
  if(exists("gtfattributes",where=object@ev,inherits=FALSE))
  {
    # Should only be present in ensemblGenome
    tbl<-data.frame(id=ev$gtf$id)
    assign("gtfattributes",merge(object@ev$gtfattributes,tbl),envir=ev)
  }
  return(new("ensemblGenome",ev=ev,basedir=object@basedir))
})

setMethod("extractTranscript","ucscGenome",function(object,transcripts)
{
  if(!exists("gtf",where=object@ev,inherits=FALSE))
    stop("[extractTranscript.ucscGenome] gtf table does not exist! Use 'read.gtf'.")
  if(is.na(match("transcript_id",names(object@ev$gtf))))
    stop("[extractTranscript.ucscGenome] 'transcript_id' column does not exist in 'gtf' table!.")
  if(!is.character(transcripts))
    stop("[extractTranscript.ucscGenome] transcripts must be character!")
  mtc<-match(transcripts,object@ev$gtf$transcript_id)
  if(any(is.na(mtc)))
  {
    cat("[extractTranscript.ucscGenome] Missing matches for ",sum(is.na(mtc))," transcripts:\n",sep="")
    print(transcripts[is.na(mtc)])
  }
  # reorder (transcrpt_id)
  dtb<-data.frame(transcript_id=transcripts)
  
  ev<-new.env()
  assign("gtf",merge(object@ev$gtf,dtb),envir=ev)
  return(new("ucscGenome",ev=ev,basedir=object@basedir))
})


setGeneric("tableTranscript.id",function(object)standardGeneric("tableTranscript.id"))
setMethod("tableTranscript.id","ensemblGenome",function(object)
{
  if(!exists("gtf",where=object@ev,inherits=FALSE))
    stop("[tableTranscript.id.ensemblGenome] gtf table does not exist! Use 'read.gtf'.")
  if(is.na(match("transcript_id",names(object@ev$gtf))))
    stop("[tableTranscript.id.ensemblGenome] gtf does not contain column 'transcript_id'. Use 'moveAttributes'")
  return(table(object@ev$gtf$transcript_id))
})
setMethod("tableTranscript.id","ucscGenome",function(object)
{
  if(!exists("gtf",where=object@ev,inherits=FALSE))
    stop("[tableTranscript.id.ucscGenome] gtf table does not exist! Use 'read.gtf'.")
  if(is.na(match("transcript_id",names(object@ev$gtf))))
    stop("[tableTranscript.id.ensemblGenome] gtf does not contain column 'transcript_id'.")
  return(table(object@ev$gtf$transcript_id))
})


setGeneric("tableTranscript.name",function(object)standardGeneric("tableTranscript.name"))
setMethod("tableTranscript.name","ensemblGenome",function(object)
{
  if(!exists("gtf",where=object@ev,inherits=FALSE))
    stop("[tableTranscript.name.ensemblGenome] gtf table does not exist! Use 'read.gtf'.")
  if(is.na(match("transcript_name",names(object@ev$gtf))))
    stop("[tableTranscript.name.ensemblGenome] gtf table does not contain column 'transcript_name'. Use 'moveAttributes'!")
  return(table(object@ev$gtf$transcript_name))
})

setGeneric("extractPaGenes",function(object)standardGeneric("extractPaGenes"))
setMethod("extractPaGenes","ensemblGenome",function(object)
{
  if(is.na(match("gene_name",names(object@ev$gtf))))
    stop("[extractPaGenes.ensemblGenome] gtf table does not contain column 'gene_name'. Use 'moveAttributes'!")
  if(is.na(match("exon_number",names(object@ev$gtf))))
    stop("[extractPaGenes.ensemblGenome] gtf table does not contain column 'exon_number'. Use 'moveAttributes'!")
  
  enpa<-extractSeqids(object,ensPrimAssembly())
  first<-function(x) return(x[1])
  dtb<-enpa@ev$gtf
  gene<-aggregate(dtb[,c("seqid","source","gene_name")],by=list(gene_id=dtb$gene_id),first)
  gstart<-aggregate(data.frame(start=dtb$start),by=list(gene_id=dtb$gene_id),min)
  gend<-aggregate(dtb[,c("end","exon_number")],by=list(gene_id=dtb$gene_id),max)
  res<-merge(gene,gstart,by="gene_id")
  res<-merge(res,gend,by="gene_id")
  res<-res[order(res$seqid,res$start),c(2,5,6,4,1,7)]
  return(res)
})

setMethod("getGenePositions","ensemblGenome",function(object)
{
  if(!exists("gtf",where=object@ev,inherits=FALSE))
    return(NULL)
  if(!is.data.frame(object@ev$gtf))
    stop("[getGenePositions.refGenome] gtf-table must be data.frame!")
  if(is.na(match("exon_number",names(object@ev$gtf))))
    stop("[getGenePositions.ensemblGenome] gtf table does not contain column 'exon_number'. Use 'moveAttributes'!")
  if(is.na(match("gene_name",names(object@ev$gtf))))
    stop("[getGenePositions.ensemblGenome] gtf table does not contain column 'gene_name'. Use 'moveAttributes'!")
  
  first<-function(x) return(x[1])
  dtb<-object@ev$gtf
  gene<-aggregate(dtb[,c("seqid","source","gene_name")],by=list(gene_id=dtb$gene_id),first)
  gstart<-aggregate(data.frame(start=dtb$start),by=list(gene_id=dtb$gene_id),min)
  gend<-aggregate(dtb[,c("end","exon_number")],by=list(gene_id=dtb$gene_id),max)
  res<-merge(gene,gstart,by="gene_id")
  res<-merge(res,gend,by="gene_id")
  res<-res[order(res$seqid,res$start),c(2,5,6,4,1,7)]
  return(res)
})

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
#  tableSeqids and extract seqids, targeted  on regex:
#  Primary intention is to retreave the primary assembly and remove haplotypes
ucPrimAssembly<-function() {return("^chr[0-9XYM]{1,2}$")}
ensPrimAssembly<-function(){return("^([0-9]{1,2})$|^[XY]|MT$")}

setGeneric("tableSeqids",function(object,regex)standardGeneric("tableSeqids"))
setMethod("tableSeqids","refGenome",function(object,regex)
{
  if(missing(regex))
    return(table(object@ev$gtf$seqid))
  else
  {
    if(!is.character(regex))
      stop("[tableSeqids.refGenome] regex must be character!")
    tbl<-table(object@ev$gtf$seqid)
    seqNames<-names(tbl)
    return(data.frame(seqid=seqNames,nAnnotations=as.numeric(tbl),lgl=grepl(regex,seqNames)))
  }
})

setGeneric("extractSeqids",function(object,regex)standardGeneric("extractSeqids"))
setMethod("extractSeqids","refGenome",function(object,regex)
{
  if(!is.character(regex))
    stop("[extractSeqids.refGenome] regex must be character!")
  lgl<-grepl(regex,object@ev$gtf$seqid)
  ans<-new(class(object),ev=new.env(),basedir=object@basedir)
  ans@ev$gtf<-object@ev$gtf[lgl,]
  # Remove unused factor levels
  ans@ev$gtf$seqid<-factor(ans@ev$gtf$seqid)
  if(exists("xref",where=object@ev,inherits=FALSE))
  {
    # Should only be present in uscsGenome
    tbl<-data.frame(kgID=unique(object@ev$gtf$gene_id))
    ans@ev$xref<-merge(object@ev$xref,tbl)
  }
  if(exists("gtfattributes",where=object@ev,inherits=FALSE))
  {
    # Should only be present in ensemblGenome
    tbl<-data.frame(id=ans@ev$gtf$id)
    ans@ev$gtfattributes<-merge(object@ev$gtfattributes,tbl)
    # Delete empty rows?
  }
  return(ans)
})

# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
# Ensembl-Exons

setClass("ensemblExons",contains="ensemblGenome")

load.ensembl.exons<-function(src,basedir)           
{
  if(is.character(src))
  {
    if(missing(basedir))    
      basedir<-dirname(src)
    ens<-new("ensemblExons",basedir=basedir)
    load(src,envir=ens@ev)
    return(invisible(ens))    
  }
  if(is(src,"url"))
  {
    if(missing(basedir))
      basedir<-getwd()
    ens<-new("ensemblExons",ev=new.env(),basedir=basedir)
    load(src,envir=ens@ev)
    close(src)
    return(invisible(ens))    
  }
  stop("[loadEnsembl] Invalid 'src' argument!")
}           


setGeneric("ensemblExons",function(object,regex)standardGeneric("ensemblExons"))
setMethod("ensemblExons","ensemblGenome",function(object)
{
  cat("[ensemblExons.ensemblGenome] Extracting tables.\n")
  cds<-extractFeature(object,"CDS")@ev$gtf
  exons<-extractFeature(object,"exon")@ev$gtf
  stacod<-extractFeature(object,"start_codon")
  stocod<-extractFeature(object,"stop_codon")
  
  # +++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
  # CDS
  cat("[ensemblExons.ensemblGenome] Adding 'CDS'.\n")
  columns<-c("start","end","transcript_id","exon_number")
  cds<-cds[,columns]
  names(cds)[1:2]<-c("cds_start","cds_end")
  ex_cds<-merge(exons,cds,by=c("transcript_id","exon_number"),all=T)
  ex_cds$cds_start<-ex_cds$cds_start-ex_cds$start
  ex_cds$cds_end<-ex_cds$end-ex_cds$cds_end
  
  # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
  # start_codon
  cat("[ensemblExons.ensemblGenome] Adding 'start_codon'.\n")
  columns<-c("start","transcript_id","exon_number")
  start_codons<-as.data.frame(stacod@ev$gtf)[,columns]
  names(start_codons)[1]<-c("start_codon")
  ex_cds_st<-merge(ex_cds,start_codons,by=c("transcript_id","exon_number"),all=T)
  ex_cds_st$start_codon<-ex_cds_st$start_codon-ex_cds_st$start
  
  # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++ #
  # stop_codon
  cat("[ensemblExons.ensemblGenome] Adding 'stop_codon'.\n")  
  columns<-c("start","transcript_id","exon_number")
  stop_codons<-as.data.frame(stocod@ev$gtf)[,columns]
  names(stop_codons)[1]<-c("stop_codon")
  ex_cds_ststp<-merge(ex_cds_st,stop_codons,by=c("transcript_id","exon_number"),all=T)
  ex_cds_ststp$stop_codon<-ex_cds_ststp$stop_codon-ex_cds_ststp$start
  
  nCol<-dim(exons)[2]
  mtc<-match(names(exons),names(ex_cds_ststp)[1:nCol])
  ex_cds_ststp<-ex_cds_ststp[,c(mtc,nCol+1:4)]
  ex_cds_ststp<-ex_cds_ststp[order(ex_cds_ststp$seqid,ex_cds_ststp$start),]
  ex_cds_ststp$exon_number<-as.numeric(ex_cds_ststp$exon_number)
  
  res<-new("ensemblExons",ev=new.env(),basedir=object@basedir)
  assign("gtf",ex_cds_ststp,envir=res@ev)
  cat("[ensemblExons.ensemblGenome] Finished.\n")
  return(res)
})

setGeneric("getSpliceTable",function(object)standardGeneric("getSpliceTable"))
setMethod("getSpliceTable","ensemblExons",function(object)
{
  gtf<-object@ev$gtf
  maxex<-max(gtf$exon_number)
  cat("[getSpliceTable.ensemblExons] Extracting junctions for maximal exon_number",maxex,".\n")  
  
  l<-list()  
  for(i in 2:maxex)
  {
    cat("\r[getSpliceTable.ensemblExons] Extracting junctions for exon_number ",format(i,width=3),"/",maxex,sep="")
    lex<-gtf[gtf$exon_number==(i-1),]
    names(lex)<-paste("left_",names(lex),sep="")
    
    rex<-gtf[gtf$exon_number==i,]
    names(rex)<-paste("right_",names(rex),sep="")
    rex$left_transcript_id<-rex$right_transcript_id
    rex$left_exon_number<-(i-1)
    
    l[[(i-1)]]<-merge(lex,rex)
    gtf<-gtf[gtf$exon_number!=(i-1),]
  }
  cat("\n[getSpliceTable.ensemblExons] Merging tables.\n")
  
  #rec<-c(10,2,14,5,6,16,17,18,11,1,12,13,7,8,9)
  #junc<-junc[,rec]
  junc<-do.call(rbind,l)
  n<-dim(junc)[1]
  row.names(junc)<-1:n
  
  cat("[getSpliceTable] Finished.\n")
  return(junc)
})



###################################################################################################
##                                                                                               ##
## overlap :  Takes tables of query and reference ranges,                                        ##
##            overlaps both and reports overlaps (+ misses) for every query range ()             ##
##                                                                                               ##
##            Results encoding                                                                   ##
##            no      OVERLAP_NO_OVER     1     query misses reference                           ##
##            r       OVERLAP_R_OVER      2     query overhangs on right side                    ##
##            b       OVERLAP_B_OVER      3     query overhangs on both sides (right and left)   ##
##            n       OVERLAP_N_OVER      4     query overhangs on neither side                  ##
##            l       OVERLAP_L_OVER      5     query overhangs on left side                     ##
##                                                                                               ##
###################################################################################################

overlap<-function(qry,ref)
{
  res<-.Call("overlap_ranges",as.integer(qry$id),as.integer(qry$start),as.integer(qry$end),as.integer(ref$id),as.integer(ref$start),as.integer(ref$end),PACKAGE="refGenome")
  return(res)
}
