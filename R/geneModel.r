
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Declaration of generics for geneModel.r
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

setGeneric("getTranscript",
            function(object, i) standardGeneric("getTranscript"))

setGeneric("getExonData",
            function(object) standardGeneric("getExonData"))

setGeneric("getCdsData",
            function(object) standardGeneric("getCdsData"))

setGeneric("geneModel",
           function(object, gene_id, interior=TRUE) standardGeneric("geneModel"))

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Generate consecutive identifer from two data.frame columns
# Only for internal use (not exported)
# Used in geneModel - function.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
uid <- function(dfr, f, s)
{
  # - - - - - - - - - - - - - - #
  # dfr : data.frame
  # f   : name of first column
  # s   : name of second column
  # - - - - - - - - - - - - - - #

  f <- f[1]
  s <- s[1]
  dfr <- dfr[order(dfr[, f], dfr[, s]), ]

  fn <- c(dfr[-1, f], dfr[nrow(dfr), f])
  sn <- c(dfr[-1, s], dfr[nrow(dfr), s])
  fne <- dfr[, f] == fn
  fse <- dfr[, s] == sn
  inc <- as.numeric(!(fne & fse))
  # Leading 1 means that uid starts with 1
  cinc <- c(1, inc[-length(inc)])
  uid <- cumsum(cinc)

  dfr$uid <- cumsum(cinc)
  return(invisible(dfr))
}


# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
# CLASS transcriptModel
# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #

# ident  : transcript_id, transcript_name, gene_id, gene_name, seq_name
# coords : start, end

.transcriptModel <- setClass("transcriptModel",
    slots=c(
        id="character",         # transcript_id
        name="character",       # transcript_name
        gene_id="character",
        gene_name="character",
        seq_name="character",
        strand="character",
        biotype="character",
        coords="integer",       # start, end of transcript
        exons="data.frame",     # start, end of exons
        cds="data.frame",       # start, end of CDS
        stcodon="integer",      # start_codon, stop_codon
        prime_utr="integer",    # five, three
        #introns="data.frame",   # begin <(?) end (Position of 1st nucleotde)
        version="integer"
    )
)

setMethod("initialize", "transcriptModel", function(.Object)
{
    .Object@id <- ""
    .Object@name <- ""
    #
    .Object@coords <- rep(0L, 2)
    names(.Object@coords) <- c("start", "end")
    #
    .Object@exons <- data.frame(start=0L, end=0L)
    .Object@cds <- data.frame(start=0L, end=0L)
    #
    .Object@stcodon <- rep(0L, 2)
    names(.Object@stcodon) <- c("start", "stop")
    #
    .Object@prime_utr <- rep(0L, 2)
    names(.Object@prime_utr) <- c("five", "three")
    #
    return(.Object)
})

setMethod("show", "transcriptModel", function(object){
    bm<-Sys.localeconv()[7]
    cat("An object of class '", class(object), "'.\n", sep="")
    cat("ID          : ", object@id, "\n")
    cat("Name        : ", object@name, "\n")
    cat("Gene ID     : ", object@gene_id, "\n")
    cat("Gene Name   : ", object@gene_name, "\n")
    cat("Start       : ", format(object@coords[1], big.mark=bm), "\t\t")
    cat("End         : ", format(object@coords[2], big.mark=bm), "\n")
    cat("Start codon : ", format(object@stcodon[1], big.mark=bm), "\t\t")
    cat("Stop  codon : ", format(object@stcodon[2], big.mark=bm), "\n")
    cat("5' prime utr: ", format(object@prime_utr[1], big.mark=bm), "\t\t")
    cat("3' prime utr: ", format(object@prime_utr[2], big.mark=bm), "\n")
    cat("Seq Name    : ", object@seq_name, "\n")
    cat("Strand      : ", object@strand, "\n")
})


plot.transcriptModel <- function(x, ylim, col, lwd=2, ...)
{
    # adjustcolor: Package grDevices
    # cols=border, background = alpha(border, 0.5)
    for(i in 1:nrow(x@exons))
        rect(
            x@exons$start[i], ylim[1],
            x@exons$end[i],   ylim[2],
            col=adjustcolor(col[1], alpha.f=0.5), border=col[1],
            lwd=lwd, ...)
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Extract TranscriptModel Object from GTF table
# For internal use only (not exported)
# Used in geneModel function
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

getTranscriptFromGtf <- function(gtf, transcript_id)
{
    transcript_id <- as.character(transcript_id[1])
    gtf <- gtf[gtf$transcript_id==transcript_id, ]
    if(nrow(gtf)==0)
        stop("Transcript_id '", transcript_id, "' not found", sep="")

    res <- new("transcriptModel") # .transcriptModel()
    res@id <- transcript_id
    res@name <- as.character(gtf$transcript_name[1])
    res@gene_id <- as.character(gtf$gene_id[1])
    res@gene_name <- as.character(gtf$gene_name[1])
    res@strand <- as.character(gtf$strand[1])
    res@seq_name <- as.character(gtf$seqid[1])
    res@coords[1] <- min(gtf$start)
    res@coords[2] <- max(gtf$end)
    res@version <- as.integer(gtf$transcript_version[1])

    wc <- which(gtf$feature=="start_codon")
    res@stcodon[1] <- ifelse(length(wc)>0, gtf$start[wc[1]], NA)

    wc <- which(gtf$feature=="stop_codon")
    res@stcodon[2] <- ifelse(length(wc)>0, gtf$start[wc[1]], NA)

    wc <- which(gtf$feature=="five_prime_utr")
    res@prime_utr[1] <- ifelse(length(wc)>0, gtf$start[wc[1]], NA)

    wc <- which(gtf$feature=="three_prime_utr")
    res@prime_utr[2] <- ifelse(length(wc)>0, gtf$start[wc[1]], NA)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - #
    # Exon data.frame
    # - - - - - - - - - - - - - - - - - - - - - - - - - - #
    exons <- gtf[gtf$feature=="exon", ]

    if(nrow(exons) > 0)
    {
        # Missing exon
        if(is.na(match("exon_number", names(exons))))
        {
            # + strand: increasing, - strand: decreasing
            exons <- exons[order(exons$start), ]
            if(exons$strand[1]=="+"){
                exons$exon_number <- 1:nrow(exons)
            }else{
                exons$exon_number <- nrow(exons):1
            }
        }else{
            exons$exon_number <- as.numeric(exons$exon_number)
        }

        exn <- c("start", "end",
                 "exon_id", "exon_number", "exon_version",
                 "seqid", "strand")

        # Missing columns are removed from column names
        mtc <- match(exn, names(exons))
        exn <- exn[!is.na(mtc)]

        res@exons <- exons[order(exons$exon_number), exn]
        rownames(res@exons) <- as.character(res@exons$exon_number)
    }else{
        cat("No exons found for transcript '", transcript_id, "'\n", sep="")

    }

    # - - - - - - - - - - - - - - - - - - - - - - - - - - #
    # CDS data.frame
    # - - - - - - - - - - - - - - - - - - - - - - - - - - #
    cds <- gtf[gtf$feature=="CDS", ]
    res@cds <- cds

    return(res)
}



setMethod("getExonData", "transcriptModel", function(object){
    return(object@exons)
})

setMethod("getCdsData", "transcriptModel", function(object){
    return(object@exons)
})


# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
# CLASS geneModel
# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #

.geneModel <- setClass("geneModel",
        slots=c(gene_id="character",
                gene_name="character",
                seq_name="character",
                strand="character",
                transcripts="character",
                coords="integer"),
        contains="refGenome")

setMethod("initialize", "geneModel", function(.Object)
{
    # ev needs to be done explicitly assigned here
    # because otherwise obscure copies appear
    .Object@ev <- new.env()
    .Object@gene_id <- ""
    .Object@gene_name <- ""
    .Object@seq_name <- ""
    .Object@strand <- "*"
    .Object@transcripts <- ""
    .Object@coords <- rep(0L, 2)

    return(.Object)
})

setGeneric("geneName", function(object) standardGeneric("geneName"))
setMethod("geneName", "geneModel", function(object) {
    return(object@gene_name)
})

setGeneric("geneName<-", function(object, value) standardGeneric("geneName<-"))
setReplaceMethod("geneName", c("geneModel", "character"), function(object, value){
    object@gene_name <- value[1]
    return(object)
})

setReplaceMethod("geneName", c("geneModel", "factor"), function(object, value){
    geneName(object) <- as.character(value[1])
    return(object)
})

setGeneric("geneId", function(object) standardGeneric("geneId"))
setMethod("geneId", "geneModel", function(object) {
    return(object@gene_name)
})

setGeneric("geneId<-", function(object, value) standardGeneric("geneId<-"))
setReplaceMethod("geneId", c("geneModel", "character"), function(object, value){
    object@gene_name <- value[1]
    return(object)
})

setReplaceMethod("geneId", c("geneModel", "factor"), function(object, value){
    geneId(object) <- as.character(value[1])
    return(object)
})




setMethod("show", "geneModel", function(object)
{
    bm<-Sys.localeconv()[7]
    cat("Object of class '", class(object), "'\n", sep="")
    cat("Gene id     : ", object@gene_id[1]  , "\n")
    cat("Gene name   : ", object@gene_name[1], "\n")
    cat("Seqid       : ", object@seq_name[1], "\n")
    cat("Strand      : ", object@strand[1], "\n")
    cat("Start       : ", format(object@coords[1], big.mark=bm), "\n")
    cat("End         : ", format(object@coords[2], big.mark=bm), "\n")
    cat("Transcripts : ", length(object@transcripts), "\n")



    if(!exists("exons", where=object@ev, inherits=FALSE))
        cat("(No exon table present)\n")
    else
    {
        n<-min(nrow(object@ev$exons), 6L)
        cat("Exon Nr     : ",
            format(nrow(object@ev$exons), big.mark = bm),
            "\n")

        print(head(object@transcripts))
    }
})

plot.geneModel <- function(x, cols=c("firebrick2", "gray40"), ...)
{
    gene_text <- paste("Gene id :", x@gene_id,
                       " Seqid : "  , x@seq_name,
                       " Strand: "  , x@strand)

    ylim=c(0,10)
    op <- par(mar=c(5,8,4,2) + 0.1)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    # Do main plot
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    plot(x@coords, ylim, type="n", bty="n", yaxt="n",
         main=paste("Gene model :", x@gene_name),
         xlab=paste("Position on seqid", x@seq_name),
         ylab="", ...)

    #bm<-Sys.localeconv()[7]
    mtext(gene_text, side=3, line=0)

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    # Strand arrow
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    if(x@strand == "+"){
        arrows(x0=x@coords[1], y0=ylim[2],
               x1=x@coords[2], y1=ylim[2],
               code=2, length=0.1)
    } else {
        arrows(x0=x@coords[1], y0=ylim[2],
               x1=x@coords[2], y1=ylim[2],
               code=1, length=0.1)
    }


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    # Draw "All exons" line
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    if(!exists("exons", envir=x@ev))
    {
        cat("[plot.geneModel] No exon data found!\n")
    }else{
        # Draw all exon boxes for gene using second color value
        # adjustcolor: Package grDevices
        # cols=border, background = alpha(border, 0.5)
        exons <- x@ev$exons
        for(i in 1:nrow(exons))
            rect(exons$start[i], 0.2, exons$end[i], 0.8,
                 col=adjustcolor(cols[2], alpha.f=0.5), border=cols[2])
    }

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    # Draw exon boxes for transcripts
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

    if(!exists("transcripts", envir=x@ev))
    {
        cat("[plot.geneModel] No transcript data found!\n")
    }else{
        trans <- x@ev$transcripts
        ntrans <- length(trans)

        ylim_t <- c(1, 8)
        ywid <- ylim_t[2] - ylim_t[1]
        allrects <- ywid * 2/3
        allgaps  <- ywid * 1/3
        rectwid  <- allrects / ntrans       # vertical number of transcript boxes
        gapwid   <- allgaps / (ntrans - 1)  # number of vertical gaps

        ylolim <- (0:(ntrans-1) * (rectwid + gapwid)) + ylim_t[1]

        for(i in 1:ntrans)
            plot(trans[[i]], ylim=c(ylolim[i], ylolim[i] + rectwid),
                 col=cols[1])
    }


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    # Draw labels on y axis
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
    axis(side=2,
            at=c(0.5, ylolim + rectwid/2),
            tick=FALSE,
            line=NA,
            labels=c("All exons",names(x@transcripts)),
            las=1, cex.axis=0.6)
    par(op)
}



setMethod("getTranscript", c("geneModel", "numeric"), function(object, i)
{
        return(object@ev$transcripts[[i]])
})

setMethod("getTranscript", c("geneModel", "character"), function(object, i)
{
    mtc <- match(i[1], names(object@ev$transcripts))
    if(!is.na(mtc))
        return(object@ev$transcripts[[mtc]])
    mtc <- match(i[1], object@transcripts)
    if(!is.na(mtc))
        return(object@ev$transcripts[[mtc]])
    stop("No Match for transcript name '", i[1], "'")
})


setMethod("geneModel", c("ensemblGenome", "character"),
    function(object, gene_id, interior=TRUE)
    {

        # Object only contains data for one single gene_id
        gene_id <- as.character(gene_id[1])
        gg <- extractByGeneId(object, gene_id)

        if(!exists("genes", envir=object@ev))
            stop("genes table missing in ensemblGenome object")

        gtb <- gg@ev$gtf # shortcut


        # Create output object
        res <- .geneModel()
        res@gene_id <- gene_id
        res@gene_name <- as.character(gg@ev$gtf$gene_name[1])

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
        # Extract gene data
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
        mtc <- match(gene_id, gg@ev$genes$gene_id)
        assign("genes", gg@ev$genes[mtc, ], envir=res@ev)

        genes<-res@ev$genes

        res@seq_name <- as.character(genes$seqid[1])
        res@strand <- as.character(genes$strand[1])
        res@coords <- c(genes$start[1], genes$end[1])
        names(res@coords) <- c("start", "end")


        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
        # Exon and transcript data eventually is skipped
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
        if(!interior)
            return(res)

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
        # Extract exon data
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

        exn <- c("seqid", "start", "end",
                 "exon_id", "exon_number", "exon_version",
                 "transcript_version", "transcript_id","transcript_name")

        # Missing columns are removed from column names
        mtc <- match(exn, names(gtb))
        exn <- exn[which(!is.na(mtc))]

        exons <- gtb[gtb$feature=="exon", exn]
        assign("exons", exons, envir=res@ev)


        # Generate coordinates of unique exons
        exons <- uid(exons, "start", "end")
        exid <- sort(unique(exons$uid))

        mtc <- match(exid, exons$uid)
        assign("uexons", exons[mtc, ], envir=res@ev)

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
        # Extract transcript data
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

        tr <- sort(unique(gtb$transcript_name))
        mtc <- match(tr, gtb$transcript_name)

        res@transcripts <- as.character(gtb$transcript_id[mtc])
        names(res@transcripts) <- tr

        l <- lapply(rep("transcriptModel", length(res@transcripts)), new)
        for(i in 1:length(l))
        {
            tro <- getTranscriptFromGtf(gtb, res@transcripts[i])
            etr <- extractTranscript(gg, res@transcripts[i])

            if(nrow(etr@ev$gtf) < 2)
            {
                # Empty splice table
                tdf <- data.frame(begin=character(0), end=character(0))
            }else{
                etj <- getSpliceTable(etr)
                utj <- unifyJuncs(etj)
                tdf <- data.frame(begin=utj@ev$gtf$lend + 1,
                                end=utj@ev$gtf$rstart - 1)
            }
            #tro@introns <- tdf

            l[[i]] <- tro
        }
        names(l) <- tr
        assign("transcripts", l, envir=res@ev)

        return(res)
    }
)

setMethod("geneModel", c("ucscGenome", "character"),
    function(object, gene_id, interior=TRUE)
{

        # Object only contains data for one single gene_id
        gene_id <- as.character(gene_id[1])
        gg <- extractByGeneId(object, gene_id)

        if(!exists("genes", envir=object@ev))
            stop("genes table missing in ensemblGenome object")

        gtb <- gg@ev$gtf # shortcut


        # Create output object
        res <- .geneModel()
        res@gene_id <- gene_id
        res@gene_name <- as.character(gg@ev$gtf$gene_name[1])

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
        # Extract gene data
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
        mtc <- match(gene_id, gg@ev$genes$gene_id)
        assign("genes", gg@ev$genes[mtc, ], envir=res@ev)

        genes<-res@ev$genes

        res@seq_name <- as.character(genes$seqid[1])
        res@strand <- as.character(genes$strand[1])
        res@coords <- c(genes$start[1], genes$end[1])
        names(res@coords) <- c("start", "end")

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
        # Exon and transcript data eventually is skipped
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
        if(!interior)
            return(res)

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
        # Extract exon data
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

        exn <- c("seqid", "start", "end",
                 "exon_id", "exon_number", "exon_version",
                 "transcript_version", "transcript_id","transcript_name")

        # Missing columns are removed from column names
        mtc <- match(exn, names(gtb))
        exn <- exn[which(!is.na(mtc))]

        exons <- gtb[gtb$feature=="exon", exn]
        assign("exons", exons, envir=res@ev)


        # Generate coordinates of unique exons
        exons <- uid(exons, "start", "end")
        exid <- sort(unique(exons$uid))

        mtc <- match(exid, exons$uid)
        assign("uexons", exons[mtc, ], envir=res@ev)

        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
        # Extract transcript data
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

        tr <- sort(unique(gtb$transcript_id))
        mtc <- match(tr, gtb$transcript_id)

        res@transcripts <- as.character(gtb$transcript_id[mtc])
        names(res@transcripts) <- tr

        l <- lapply(rep("transcriptModel", length(res@transcripts)), new)
        for(i in 1:length(l))
        {
            l[[i]] <- getTranscriptFromGtf(gtb, res@transcripts[i])
        }
        names(l) <- tr
        assign("transcripts", l, envir=res@ev)


        return(res)
    }
)

# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
# CLASS geneList
# + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + #



.geneList <- setClass("geneList",
                             slots=c(
                                l="list"
                             )
)


setMethod("initialize", "geneList", function(.Object){

    .Object@l <- list()
    return(.Object)
})

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# S3 generics
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
setMethod("length", "geneList", function(x){ return(length(x@l)) })
setMethod("names", "geneList", function(x) { return(names(x@l))})

setMethod("names<-", c("geneList", "character"),
    function(x, value)
    {
        names(x@l) <- value
        return(x)
    }
)

setMethod("names<-", c("geneList", "numeric"),
          function(x, value)
          {
              names(x@l) <- value
              return(x)
          }
)

setMethod("show", "geneList", function(object)
{
    bm<-Sys.localeconv()[7]
    cat("Object of class '", class(object), "'\n", sep="")
    cat("Length      : ", length(object)  , "\n")
    cat("Names:\n")
    print(names(object))
})

setMethod("[", signature="geneList", function(x, i)
{
    if(length(i) == 1)
        return(x@l[[i]])

    res <- .geneList()
    res@l <- x@l[i]
    return(res)
})

setMethod("+", signature=c("geneList", "geneList"), function(e1, e2){
    res <- .geneList()
    res@l <- c(e1@l, e2@l)
    return(res)
})

setMethod("+", c("geneModel", "geneModel"), function(e1, e2){
    res <- .geneList()
    res@l <- list(e1, e2)
    names(res@l) <- c(e1@gene_id, e2@gene_id)
    return(res)
})


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# Creation of geneList objects
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
setGeneric("geneList", function(ref, genes, interior=TRUE)
                                            standardGeneric("geneList"))

setMethod("geneList", c("ensemblGenome", "character"),
    function(ref, genes, interior=TRUE)
{
    ng <- length(genes)
    # convert
    genl <- split(genes, 1:ng)
    names(genl) <- genes
    getGeneModel <- function(x) { return(geneModel(ref, x, interior))}

    gl <- .geneList()
    gl@l <- lapply(genl, getGeneModel)
    return(gl)
})

setMethod("geneList", c("ensemblGenome", "factor"),
    function(ref, genes, interior=TRUE)
{
    return(geneList(ref, as.character(genes), interior))
})

setMethod("geneList", c("ucscGenome", "character"),
    function(ref, genes, interior=TRUE)
{
    ng <- length(genes)
    # convert
    genl <- split(genes, 1:ng)
    names(genl) <- genes
    getGeneModel <- function(x) { return(geneModel(ref, x, interior))}

    gl <- .geneList()
    gl@l <- lapply(genl, getGeneModel)
    return(gl)
})

setMethod("geneList", c("ucscGenome", "factor"),
    function(ref, genes, interior=TRUE)
{
    return(geneList(ref, as.character(genes), interior))
})


