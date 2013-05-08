/**************************************************************************************************
 **************************************************************************************************
 *
 * Project	:	refGenome
 * Created	:	19.03.2012
 * Author	:	W. Kaisers
 *
 * Content	:	Managing genomic reference data for usage in R
 *
 * Version	:	1.0.0
 *
 * Changelog	:
 * 04.Jun.12	:	get_ucsc_exon_bound_df	now has unique id-values for use as PRIM KEY
 * 06.Jun.12	:	Removed memory leak 	token_list module is now valgrind tested
 * 13.Jun.12	:	get_ens_attribute_df 	additionally has id row for use as PRIM KEY
 * 17.Okt.12	:	split_gtf_attr			Added function for splitting gtf attribute column (valgrind tested)
 *
 **************************************************************************************************
 **************************************************************************************************/

#ifndef REFGENOME_C_
#define REFGENOME_C_
#include "refGenome.h"

static const int buf_size=2048; // buffer size for printing ints into chars

// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +
//
// Calculate exon_number from subsequent transcript and start values
// Expects ordering by transcript, seqid, start, end
//
// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +

SEXP get_splice_juncs(SEXP pTranscript,SEXP pId,SEXP pStart,SEXP pEnd)
{
	// + + + + + + + + + + + + + + + + + + + + + + + + + +
	// Expects that pTranscript is ordered
	// so that consecutive equal pTranscript values
	// represent junctions
	// + + + + + + + + + + + + + + + + + + + + + + + + + +

	// + + + + + + + + + + + + + + + + + + + + + + + + + +
	// Check incoming
	if(TYPEOF(pTranscript)!=INTSXP)
		error("[get_splice_juncs] pTranscript must be INT!");
	if(TYPEOF(pId)!=INTSXP)
		error("[get_splice_juncs] pId must be INT!");
	if(TYPEOF(pStart)!=INTSXP)
		error("[get_splice_juncs] pStart must be INT!");
	if(TYPEOF(pEnd)!=INTSXP)
		error("[get_splice_juncs] pEnd must be INT!");

	int inRow=LENGTH(pTranscript);
	if((LENGTH(pId)!=inRow) | (LENGTH(pStart)!=inRow) | (LENGTH(pEnd)!=inRow))
		error("[get_splice_juncs] All arguments must have same length!");


	int *tr=INTEGER(pTranscript), *id=INTEGER(pId);
	int *start=INTEGER(pStart), *end=INTEGER(pEnd);

	// + + + + + + + + + + + + + + + + + + + + + + + + + +
	// Count number of splice junctions
	unsigned nJunc=0, i,j;
	for(i=1,j=0;i<inRow;++i,++j)
	{
		if(tr[i]==tr[j])
			++nJunc;
	}

	// + + + + + + + + + + + + + + + + + + + + + + + + + +
	// Create output vectors
	SEXP pLexid;
	PROTECT(pLexid=allocVector(INTSXP,nJunc));
	SEXP pRexid;
	PROTECT(pRexid=allocVector(INTSXP,nJunc));
	SEXP pLstart;
	PROTECT(pLstart=allocVector(INTSXP,nJunc));
	SEXP pLend;
	PROTECT(pLend=allocVector(INTSXP,nJunc));
	SEXP pRstart;
	PROTECT(pRstart=allocVector(INTSXP,nJunc));
	SEXP pRend;
	PROTECT(pRend=allocVector(INTSXP,nJunc));

	int *lexid=INTEGER(pLexid);
	int *rexid=INTEGER(pRexid);
	int *lstart=INTEGER(pLstart);
	int *lend=INTEGER(pLend);
	int *rstart=INTEGER(pRstart);
	int *rend=INTEGER(pRend);

	// + + + + + + + + + + + + + + + + + + + + + + + + + +
	// Calculate splice pairs
	unsigned iJunc=0;
	for(i=1,j=0;i<inRow;++i,++j)
	{
		if(iJunc>=nJunc)
			error("[get_splice_juncs] iJunc error: i=%i\tnJunc=%i\n",i,nJunc);

		if(tr[i]==tr[j])
		{
			lexid[iJunc]=id[j];
			lstart[iJunc]=start[j];
			lend[iJunc]=end[j];

			rexid[iJunc]=id[i];
			rstart[iJunc]=start[i];
			rend[iJunc]=end[i];
			++iJunc;
		}
	}

	// + + + + + + + + + + + + + + + + + + + + + + + + + +
	// Assemble result
	SEXP dflist;
	PROTECT(dflist=allocVector(VECSXP,6));
	SET_VECTOR_ELT(dflist,0,pLexid);
	SET_VECTOR_ELT(dflist,1,pRexid);
	SET_VECTOR_ELT(dflist,2,pLstart);
	SET_VECTOR_ELT(dflist,3,pLend);
	SET_VECTOR_ELT(dflist,4,pRstart);
	SET_VECTOR_ELT(dflist,5,pRend);

	// + + + + + + + + + + + + + + + + + + + + + + + + + +
	// Column Names
	SEXP col_names;
	PROTECT(col_names=allocVector(STRSXP,6));
	SET_STRING_ELT(col_names,0,mkChar("lexid"));
	SET_STRING_ELT(col_names,1,mkChar("rexid"));
	SET_STRING_ELT(col_names,2,mkChar("lstart"));
	SET_STRING_ELT(col_names,3,mkChar("lend"));
	SET_STRING_ELT(col_names,4,mkChar("rstart"));
	SET_STRING_ELT(col_names,5,mkChar("rend"));
	setAttrib(dflist,R_NamesSymbol,col_names);

	// + + + + + + + + + + + + + + + + + + + + + + + + + +
	// Row Names
	SEXP row_names;
    PROTECT(row_names=allocVector(STRSXP,nJunc));

	char *buf=(char*) calloc(buf_size,sizeof(char));
    for(i=0;i<nJunc;++i)
    {
    	sprintf(buf,"%i",i+1);
    	SET_STRING_ELT(row_names,i,mkChar(buf));
    }
    free(buf);
    setAttrib(dflist,R_RowNamesSymbol,row_names);

    // + + + + + + + + + + + + + + + + + + + + + + + + + +
	setAttrib(dflist,R_ClassSymbol,mkString("data.frame"));

	UNPROTECT(9);
	return dflist;
}

SEXP unify_splice_juncs(SEXP pSeqid,SEXP pLstart,SEXP pLend,SEXP pRstart,SEXP pRend, SEXP pId, SEXP pGeneId, SEXP pStrand)
{
	// + + + + + + + + + + + + + + + + + + + + + + + + + +
	// Calculates unique junction coordinate values (uJunc)
	// from a sorted junction list
	//
	// Expects that pSeqid, pLend and pRstart are ordered
	// so that identical junctions are
	// represented by consecutive equal values
	// + + + + + + + + + + + + + + + + + + + + + + + + + +

	// + + + + + + + + + + + + + + + + + + + + + + + + + +
	// Check incoming

	if(TYPEOF(pSeqid)!=INTSXP)
		error("[unify_splice_juncs] pSeqid must be INT!");
	if(TYPEOF(pLstart)!=INTSXP)
		error("[unify_splice_juncs] pLstart must be INT!");
	if(TYPEOF(pLend)!=INTSXP)
		error("[unify_splice_juncs] pLend must be INT!");
	if(TYPEOF(pRstart)!=INTSXP)
		error("[unify_splice_juncs] pRstart must be INT!");
	if(TYPEOF(pRend)!=INTSXP)
		error("[unify_splice_juncs] pRend must be INT!");
	if(TYPEOF(pId)!=INTSXP)
		error("[unify_splice_juncs] pId must be INT!");
	if(TYPEOF(pGeneId)!=INTSXP)
		error("[unify_splice_juncs] pGeneId must be INT!");
	if(TYPEOF(pStrand)!=INTSXP)
		error("[unify_splice_juncs] pStrand must be INT!");

	unsigned i,j,k, nJunc, n=LENGTH(pId), nSites;

	int 	*seqid=INTEGER(pSeqid),
			*lstart=INTEGER(pLstart), *lend=INTEGER(pLend),
			*rstart=INTEGER(pRstart), *rend=INTEGER(pRend),
			*id=INTEGER(pId),
			*gid=INTEGER(pGeneId),    *strand=INTEGER(pStrand);

	int usq,ule,urs; // unified splice coordinates

	// First row identifies first nJunc site
	usq=seqid[0];
	ule=lend[0];
	urs=rstart[0];
	nJunc=1;
	i=0;
	//Rprintf("[unify_splice_juncs] Start junc %i at i= %i.\n",nJunc,i);

	// Check all subsequent rows for position equality
	for(i=1;i<n;++i)
	{
		if((seqid[i]!=usq) | (lend[i]!=ule) | (rstart[i]!=urs))
		{
			++nJunc;
			//Rprintf("[unify_splice_juncs] Start junc %i at i= %i.\n",nJunc,i);
			usq=seqid[i];
			ule=lend[i];
			urs=rstart[i];
		}
	}
	//Rprintf("[unify_splice_juncs] Found %i juncs.\n",nJunc);

	if(nJunc==0)
		return R_NilValue;

	// Create output vectors
	SEXP puId;
	PROTECT(puId=allocVector(INTSXP,nJunc));

	SEXP puSeq;
	PROTECT(puSeq=allocVector(INTSXP,nJunc));

	SEXP puLstart;
	PROTECT(puLstart=allocVector(INTSXP,nJunc));

	SEXP puLend;
	PROTECT(puLend=allocVector(INTSXP,nJunc));

	SEXP puRstart;
	PROTECT(puRstart=allocVector(INTSXP,nJunc));

	SEXP puRend;
	PROTECT(puRend=allocVector(INTSXP,nJunc));

	SEXP puNsites; // Number of junction sites per uJunc
	PROTECT(puNsites=allocVector(INTSXP,nJunc));

	SEXP puGene; // Gene-id associated with exon table
	PROTECT(puGene=allocVector(INTSXP,nJunc));

	SEXP puStrand;
	PROTECT(puStrand=allocVector(INTSXP,nJunc));

	SEXP pFexid; // First exon id (points into exon table)
	PROTECT(pFexid=allocVector(INTSXP,nJunc));


	// Calculate values for sites
	int     *uid=INTEGER(puId),         *useq=INTEGER(puSeq),
			*ulstart=INTEGER(puLstart), *ulend=INTEGER(puLend),
			*urstart=INTEGER(puRstart), *urend=INTEGER(puRend),
			*uNsites=INTEGER(puNsites), *uId=INTEGER(pFexid),
			*uGene=INTEGER(puGene),     *uStrand=INTEGER(puStrand);

	int min_lstart,max_rend;
	const unsigned nGenes=10;
	unsigned geneId[nGenes];  // geneId
	unsigned geneCt[nGenes];  // gene-count
	unsigned geneStr[nGenes]; // gene-strand
	unsigned mxGct,mxGid,mxStr;    // maxGeneCount, maxGeneId, maxStrand



	j=0; // write index of actual uJunc
	i=0; // read  index of actual junc

	// Start reading in row 0
	// Row 0 identifies first uJunc
	uid[j]=j+1;
	useq[j]=seqid[i];
	ulend[j]=lend[i];
	urstart[j]=rstart[i];
	uId[j]=id[i];
	min_lstart=lstart[i];
	max_rend=rend[i];
	nSites=1;
	//Rprintf("[unify_splice_juncs] Start junc %i at i= %i.\n",j,i);

	geneId[0]=gid[i];
	geneStr[0]=strand[i];
	geneCt[0]=1;
	for(k=1;k<nGenes;++k)
		geneId[k]=0;

	// Check all subsequent rows for position equality
	for(i=1;i<n;++i)
	{
		if((seqid[i]!=useq[j]) | (lend[i]!=ulend[j]) | (rstart[i]!=urstart[j]))
		{
			// Row i contains new uJunc-site

			// + + + + + + + + + + + + + + + + + + + + + +
			// Complete values for last uJunc-site
			uNsites[j]=nSites;
			ulstart[j]=min_lstart;
			urend[j]=max_rend;

			// Get some gene-id with "maximal" gene-id count
			mxGct=0;
			mxGid=0;
			mxStr=0;
			for(k=0;k<nGenes;++k)
			{
				if(geneId[k]==0)
					break;

				if(mxGct<geneCt[k])
				{
					mxGct=geneCt[k];
					mxGid=geneId[k];
					mxStr=geneStr[k];
				}
			}
			uGene[j]=mxGid;
			uStrand[j]=mxStr;

			//Rprintf("[unify_splice_juncs] Start junc %i at i= %i.\n",j,i);

			// + + + + + + + + + + + + + + + + + + + + + +
			// Set values for next uJunc-site
			++j;
			if(j>=nJunc)
				error("[unify_splice_juncs] Write index exceeding nJunc limit!");

			uid[j]=j+1;
			useq[j]=seqid[i];
			ulend[j]=lend[i];
			urstart[j]=rstart[i];
			uId[j]=id[i];
			min_lstart=lstart[i];
			max_rend=rend[i];
			nSites=1;

			geneId[0]=gid[i];
			geneStr[0]=strand[i];
			geneCt[0]=1;

			for(k=1;k<nGenes;++k)
			{
				geneId[k]=0;
				geneCt[k]=0;
			}

		}
		else
		{
			// + + + + + + + + + + + + + + + + + + + + + +
			// Row i is part of actual uJunc-site
			++nSites;

			// lstart and rend
			if(lstart[i]<min_lstart)
				min_lstart=lstart[i];
			if(rend[i]>max_rend)
				max_rend=rend[i];

			// Count geneId's
			for(k=0;k<nGenes;++k)
			{
				if(gid[i]==geneId[k])
				{
					++(geneCt[k]);
					break;
				}
				if(geneId[k]==0)
				{
					geneId[k]=gid[i];
					geneStr[k]=strand[i];
					geneCt[k]=1;
					break;
				}
			}
		}
	}

	// + + + + + + + + + + + + + + + + + + + + + +
	// Complete values for last uJunc-site
	uNsites[j]=nSites;
	ulstart[j]=min_lstart;
	urend[j]=max_rend;

	// Get some gene-id with "maximal" gene-id count
	mxGct=0;
	mxGid=0;
	mxStr=0;
	for(k=0;k<nGenes;++k)
	{
		if(geneId[k]==0)
			break;

		if(mxGct<geneCt[k])
		{
			mxGct=geneCt[k];
			mxGid=geneId[k];
			mxStr=geneStr[k];
		}
	}
	uGene[j]=mxGid;
	uStrand[j]=mxStr;


	// Assemble output
	SEXP dflist;
	PROTECT(dflist=allocVector(VECSXP,10));

	SET_VECTOR_ELT(dflist,0,puId);
	SET_VECTOR_ELT(dflist,1,puSeq);
	SET_VECTOR_ELT(dflist,2,puLstart);
	SET_VECTOR_ELT(dflist,3,puLend);
	SET_VECTOR_ELT(dflist,4,puRstart);
	SET_VECTOR_ELT(dflist,5,puRend);
	SET_VECTOR_ELT(dflist,6,puNsites);
	SET_VECTOR_ELT(dflist,7,puGene);
	SET_VECTOR_ELT(dflist,8,puStrand);
	SET_VECTOR_ELT(dflist,9,pFexid);

	// + + + + + + + + + + + + + + + + + + + + + + + + + +
	// Column Names
	SEXP col_names;
	PROTECT(col_names=allocVector(STRSXP,10));
	SET_STRING_ELT(col_names,0,mkChar("id"));
	SET_STRING_ELT(col_names,1,mkChar("seqid"));
	SET_STRING_ELT(col_names,2,mkChar("lstart"));
	SET_STRING_ELT(col_names,3,mkChar("lend"));
	SET_STRING_ELT(col_names,4,mkChar("rstart"));
	SET_STRING_ELT(col_names,5,mkChar("rend"));
	SET_STRING_ELT(col_names,6,mkChar("nSites"));
	SET_STRING_ELT(col_names,7,mkChar("gene_id"));
	SET_STRING_ELT(col_names,8,mkChar("strand"));
	SET_STRING_ELT(col_names,9,mkChar("fexid"));
	setAttrib(dflist,R_NamesSymbol,col_names);

	// + + + + + + + + + + + + + + + + + + + + + + + + + +
	// Row Names
	SEXP row_names;
    PROTECT(row_names=allocVector(STRSXP,nJunc));

	char *buf=(char*) calloc(buf_size,sizeof(char));
    for(i=0;i<nJunc;++i)
    {
    	sprintf(buf,"%i",i+1);
    	SET_STRING_ELT(row_names,i,mkChar(buf));
    }
    free(buf);
    setAttrib(dflist,R_RowNamesSymbol,row_names);

    // + + + + + + + + + + + + + + + + + + + + + + + + + +
	setAttrib(dflist,R_ClassSymbol,mkString("data.frame"));

	UNPROTECT(13);
	return dflist;
}


SEXP get_exon_number(SEXP pTranscript,SEXP pSeqid, SEXP pStart, SEXP pEnd)
{
	if(TYPEOF(pTranscript)!=INTSXP)
		error("[get_exon_number] pTranscript must be INT!");
	if(TYPEOF(pSeqid)!=INTSXP)
		error("[get_exon_number] pSeqid must be INT!");
	if(TYPEOF(pStart)!=INTSXP)
		error("[get_exon_number] pStart must be INT!");
	if(TYPEOF(pEnd)!=INTSXP)
		error("[get_exon_number] pEnd must be INT!");

	int n=LENGTH(pTranscript);
	if(LENGTH(pSeqid)!=n || LENGTH(pSeqid)!=n || LENGTH(pStart)!=n || LENGTH(pEnd)!=n)
		error("[get_exon_number] All args must have same length!");

	SEXP res;
	PROTECT(res=allocVector(INTSXP,n));
	int i,j, exon_number=1, nSeqMm=0, nStartMm=0;

	INTEGER(res)[0]=exon_number;

	for(i=1,j=0;i<n;++i,++j)
	{
		if(INTEGER(pTranscript)[j]==INTEGER(pTranscript)[i])
		{
			// same transcript
			++exon_number;

			// Security checks
			if(INTEGER(pSeqid)[j]!=INTEGER(pSeqid)[i])
				++nSeqMm;
			if(INTEGER(pEnd)[j]>=INTEGER(pStart)[i])
				++nStartMm;
		}
		else
			exon_number=1; // new transcript

		INTEGER(res)[i]=exon_number;
	}

	if(nSeqMm>0)
		Rprintf("[get_exon_number] Found %i sequence mismatches!\n",nSeqMm);
	if(nStartMm>0)
		Rprintf("[get_exon_number] Found %i start-end mismatches!\n",nStartMm);

	UNPROTECT(1);
	return res;
}



// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +
//
// overlap: Compares list of query and reference ranges
// and reports overlaps in data.frame
//
// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +

SEXP overlap_ranges(SEXP qryid, SEXP qrystart, SEXP qryend, SEXP refid,SEXP refstart,SEXP refend)
{
	// expects qryend >= qrystart, refend >= refstart
	// expects qrystart and refstart ascending sorted

	// Check type of incoming args
	if(TYPEOF(qryid)!=INTSXP)
		error("[overlap_ranges] qryid is no INT!\n");
	if(TYPEOF(qrystart)!=INTSXP)
		error("[overlap_ranges] qrystart is no INT!\n");
	if(TYPEOF(qryend)!=INTSXP)
		error("[overlap_ranges] qryend is no INT!\n");
	if(TYPEOF(refid)!=INTSXP)
		error("[overlap_ranges] refid is no INT!\n");
	if(TYPEOF(refstart)!=INTSXP)
		error("[overlap_ranges] refstart is no INT!\n");
	if(TYPEOF(refend)!=INTSXP)
		error("[overlap_ranges] refend is no INT!\n");


	// Check size of incoming args
	// Size of qry determines size of output parameters
	unsigned nRows=LENGTH(qryid);
	if(nRows==0)
		error("[overlap_ranges] qryid had length zero!");
	if(LENGTH(qrystart)!=nRows)
		error("[overlap_ranges] length(qrystart)!=length(qryid)!");
	if(LENGTH(qryend)!=nRows)
		error("[overlap_ranges] length(qryend)!=length(qryid)!");

	unsigned nRef=LENGTH(refid);
	if(nRef==0)
		error("[overlap_ranges] refid has length zero!");
	if(LENGTH(refstart)!=nRef)
		error("[overlap_ranges] length(refstart)!=length(refid)!");
	if(LENGTH(refend)!=nRef)
		error("[overlap_ranges] length(refstart)!=length(refend)!");

	// read args
	unsigned *qid   =(unsigned*)INTEGER(qryid);
	unsigned *qstart=(unsigned*)INTEGER(qrystart);
	unsigned *qend  =(unsigned*)INTEGER(qryend);
	unsigned *rid   =(unsigned*)INTEGER(refid);
	unsigned *rstart=(unsigned*)INTEGER(refstart);
	unsigned *rend  =(unsigned*)INTEGER(refend);

	unsigned nProtected=0;
	unsigned qidx=0, ridx=0; 	// query and ref indices + max idices

	// Column 0: overlap code
	SEXP vov;
	PROTECT(vov=allocVector(INTSXP,nRows));
	++nProtected;

	// Column 1: difference to annotated position
	SEXP vldiff;
	PROTECT(vldiff=allocVector(INTSXP,nRows));
	++nProtected;

	// Column 2: difference to annotated position
	SEXP vrdiff;
	PROTECT(vrdiff=allocVector(INTSXP,nRows));
	++nProtected;

	// Column 3: id of query item
	SEXP vqid;
	PROTECT(vqid=allocVector(INTSXP,nRows));
	++nProtected;

	// Column 4: id of ref item
	SEXP vrid;
	PROTECT(vrid=allocVector(INTSXP,nRows));
	++nProtected;

	while((qidx<nRows) & (ridx<nRef))
	{
		// qry misses right
		if(qstart[qidx] > rend[ridx])
		{
			++ridx;
			continue;
		}

		if(qend[qidx] >rend[ridx])
		{
			if(qstart[qidx] >= rstart[ridx])
			{
				// qry overhang right: R_OVER
				INTEGER(vqid)  [qidx]=qid[qidx];
				INTEGER(vrid)  [qidx]=rid[ridx];
				INTEGER(vov)   [qidx]=OVERLAP_R_OVER;
				INTEGER(vldiff)[qidx]=qstart[qidx]-rstart[ridx];
				INTEGER(vrdiff)[qidx]=qend[qidx]-rend[ridx];
				++qidx;
				continue;
			}
			else
			{
				// qry overhang on both sides: B_OVER
				INTEGER(vqid)  [qidx]=qid[qidx];
				INTEGER(vrid)  [qidx]=rid[ridx];
				INTEGER(vov)   [qidx]=OVERLAP_B_OVER;
				INTEGER(vldiff)[qidx]=rstart[ridx]-qstart[qidx];
				INTEGER(vrdiff)[qidx]=qend[qidx]-rend[ridx];
				++qidx;
				continue;
			}
		}
		else if(qend[qidx] >= rstart[ridx])
		{
			if(qstart[qidx] >= rstart[ridx])
			{
				// no qry-overhang: N_OVER
				INTEGER(vqid)  [qidx]=qid[qidx];
				INTEGER(vrid)  [qidx]=rid[ridx];
				INTEGER(vov)   [qidx]=OVERLAP_N_OVER;
				INTEGER(vldiff)[qidx]=qstart[qidx]-rstart[ridx];
				INTEGER(vrdiff)[qidx]=rend[ridx]-qend[qidx];
				++qidx;
				continue;
			}
			else
			{
				// qry overhang left: L_OVER
				INTEGER(vqid)  [qidx]=qid[qidx];
				INTEGER(vrid)  [qidx]=rid[ridx];
				INTEGER(vov)   [qidx]=OVERLAP_L_OVER;
				INTEGER(vldiff)[qidx]=rstart[ridx]-qstart[qidx];
				INTEGER(vrdiff)[qidx]=rend[ridx]-qend[qidx];
				++qidx;
				continue;
			}
		}
		else
		{
			// No overlap
			INTEGER(vqid) [qidx]=qid[qidx];
			INTEGER(vrid) [qidx]=0;
			INTEGER(vov)  [qidx]=OVERLAP_NO_OVER;
			// rdiff gives distance to next ref on right side
			INTEGER(vrdiff)[qidx]=rstart[ridx]-qend[qidx];

			// ldiff gives distance to next ref on left side (or 0 if not exists)
			if(ridx>0)
				INTEGER(vldiff)[qidx]=qstart[qidx]-rend[ridx-1];
			else
				INTEGER(vldiff)[qidx]=0;

			++qidx;
			continue;
		}
	}
	// Process remaining qry rows when rightmost ref ranges are passed
	if((qidx<nRows) & (qstart[qidx] > rend[nRef-1]))
	{
		unsigned refend=rend[nRef-1];
		while(qidx<nRows)
		{
			INTEGER(vqid)  [qidx]=qid[qidx];
			INTEGER(vrid)  [qidx]=0;
			INTEGER(vov)   [qidx]=OVERLAP_NO_OVER;
			INTEGER(vrdiff)[qidx]=0;
			INTEGER(vldiff)[qidx]=qstart[qidx]-refend;
			++qidx;
		}
	}

	// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +
	// Convert vov to factor
	SEXP levs;
	int nLevels=5;
	PROTECT(levs=allocVector(STRSXP,nLevels));
	++nProtected;

	SET_STRING_ELT(levs,0,mkChar("no"));
	SET_STRING_ELT(levs,1,mkChar("r"));
	SET_STRING_ELT(levs,2,mkChar("b"));
	SET_STRING_ELT(levs,3,mkChar("n"));
	SET_STRING_ELT(levs,4,mkChar("l"));
	setAttrib(vov,R_LevelsSymbol,levs);

	SEXP csymb;
	PROTECT(csymb=mkString("factor"));
	++nProtected;
	setAttrib(vov,R_ClassSymbol,csymb);

	// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +
	// Create data.frame
	unsigned nCols=5;
	SEXP dflist;
	PROTECT(dflist=allocVector(VECSXP,nCols));
	++nProtected;

	SET_VECTOR_ELT(dflist,0,vov);
	SET_VECTOR_ELT(dflist,1,vldiff);
	SET_VECTOR_ELT(dflist,2,vrdiff);
	SET_VECTOR_ELT(dflist,3,vqid);
	SET_VECTOR_ELT(dflist,4,vrid);

	// Column Names
	SEXP col_names;
	PROTECT(col_names=allocVector(STRSXP,nCols));
	++nProtected;

	SET_STRING_ELT(col_names,0,mkChar("overlap"));
	SET_STRING_ELT(col_names,1,mkChar("leftDiff"));
	SET_STRING_ELT(col_names,2,mkChar("rightDiff"));
	SET_STRING_ELT(col_names,3,mkChar("queryid"));
	SET_STRING_ELT(col_names,4,mkChar("refid"));
	setAttrib(dflist,R_NamesSymbol,col_names);

	// Create row names for data.frame
	SEXP row_names;
    PROTECT(row_names=allocVector(STRSXP,nRows));
    ++nProtected;

	int buf_size=1024;
	char *buf=(char*) calloc(buf_size,sizeof(char));
    for(qidx=0;qidx<nRows;++qidx)
    {
    	sprintf(buf,"%i",qidx);
    	SET_STRING_ELT(row_names,qidx,mkChar(buf));
    }
    free(buf);
    setAttrib(dflist,R_RowNamesSymbol,row_names);

    // Make list to data.frame
	setAttrib(dflist,R_ClassSymbol,mkString("data.frame"));
	UNPROTECT(nProtected);

	return dflist;
}

// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +#
//                                                                                                  #
// Split gtf Attribute column data                                                                  #
// and return list with two data.frames                                                             #
//                                                                                                  #
// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +#

// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +#
// A: The Brent Lab (Washington University, St.Louis)                                               #
// http://mblab.wustl.edu/GTF2.html                                                                 #
// [attributes] All four features have the same two mandatory attributes at the end of the record:  #
//                                                                                                  #
// gene_id value;       A globally unique identifier for the genomic source of the transcript       #
// transcript_id value; A globally unique identifier for the predicted transcript.                  #
//                                                                                                  #
// These attributes are designed for handling multiple transcripts from the same genomic region.    #
// Any other attributes or comments must appear after these two and will be ignored.                #
//                                                                                                  #
// Attributes must end in a semicolon which must then be separated from the start of any subsequent #
// attribute by exactly one space character (NOT a tab character). Textual attributes *should* be   #
// surrounded by doublequotes.                                                                      #
//                                                                                                  #
// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +#
// B: Wellcome Trust Sanger Institute                                                               #
// http://www.sanger.ac.uk/resources/software/gff/spec.html                                         #
//                                                                                                  #
// Free text values *must* be quoted with double quotes.                                            #
// Note: all non-printing characters in such free text value strings (e.g. newlines, tabs, control  #
// characters, etc) must be explicitly represented by their C (UNIX) style backslash-escaped        #
// representation (e.g. newlines as '\n', tabs as '\t').                                            #
// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +#


SEXP split_gtf_attr(SEXP id_vec,SEXP attr_vec)
{
	// Check type of incoming args
	if(TYPEOF(id_vec)!=INTSXP)
		error("[split_gtf_attr] id_vec is no INT!\n");
	if(TYPEOF(attr_vec)!=STRSXP)
		error("[split_gtf_attr] attr_vec is no STR: %i!\n",TYPEOF(attr_vec));

	unsigned long int i,n;
	n=LENGTH(id_vec);
	if(LENGTH(attr_vec)!=n)
		error("[split_gtf_attr] id_vec and attr_vec must have same length!\n");

	unsigned nProtected=0;
	int *id=INTEGER_POINTER(id_vec);

	// Column 1: id vector
	SEXP idvec;
	PROTECT(idvec=allocVector(INTSXP,n));
	++nProtected;

	// Column 2: gene_id
	SEXP geneIdVec;
	PROTECT(geneIdVec=allocVector(STRSXP,n));
	++nProtected;

	// Column 3: transcript_id
	SEXP transIdVec;
	PROTECT(transIdVec=allocVector(STRSXP,n));
	++nProtected;

	// flag characters
	const char space=' ', delim=';', quote='"', zero='\0';
	// iterator positions
	const char *token_first,*token_second,*iter;
	// string length
	unsigned long first_len,second_len;
	// attribute index
	unsigned attr_index;

	ptr_pair_list *l=ptr_pair_list_init();
	for(i=0;i<n;++i)
	{
		INTEGER(idvec)[i]=id[i];
		iter=CHAR(STRING_ELT(attr_vec,i));
		attr_index=1;

		// Skip leading spaces
		while((*iter==space))
			++iter;

		while(*iter!=zero)
		{
			// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +
			// First token
			// Take start position
			token_first=iter;
			// Rprintf("[split_gtf_attr] token_first: '%s'\n",token_first);

			if((*iter!='g') && (attr_index==1))
				error("[split_gtf_attr] First item must be 'gene_id': '%s'!",iter);
			if((*iter!='t') && (attr_index==2))
				error("[split_gtf_attr] Second item must be 'transcript_id': '%s'!",iter);

			// proceed until space and take length
			while((*iter!=space) && (*iter!=zero))
				++iter;
			if(*iter==zero)
				error("[split_gtf_attr] Found end of string in first token in line %lu: '%s'!\n",i+1,iter);
			if(iter==token_first)
				error("[split_gtf_attr] First token ist empty: '%s'!\n",iter);

			first_len=iter-token_first;

			// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +
			// second token:
			// skip spaces and quotes, then take start position
			while((*iter==space) || (*iter==quote))
				++iter;
			token_second=iter;
			// Rprintf("[split_gtf_attr] token_second: '%s'\n",token_second);

			// proceed until space or quote
			while((*iter!=space) && (*iter!=quote) && (*iter!=delim) && (*iter!=zero))
				++iter;
			second_len=iter-token_second;

			// second token may be closed by quote
			if(*iter==quote)
				++iter;

			// second token must be terminated by delim
			if(*iter!=delim)
				error("[split_gtf_attr] Second token must end on ';': '%s'!",token_second);
			++iter;

			// There may be terminating spaces
			while(*iter==space)
				++iter;

			// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +
			//  Process collected pointer positions
			//  First token must be gene_id
			if(attr_index==1)
				SET_STRING_ELT(geneIdVec,i,mkCharLen(token_second,second_len));
			// Second token must be transcript_id
			else if(attr_index==2)
				SET_STRING_ELT(transIdVec,i,mkCharLen(token_second,second_len));
			// Probably some more tokens
			else
				ptr_pair_list_push_back(l,token_first,first_len,token_second,second_len,id[i]);
			++attr_index;
		}
	}

	// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +
	//  Create output data.frames
	int nCols, nRows;
	char buf[20];

	// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +
	// Merge id,gene_id and transcript_id into idflist data.frame
	nCols=3;
	nRows=n;
	SEXP idflist;
	PROTECT(idflist=allocVector(VECSXP,nCols));
	++nProtected;

	SET_VECTOR_ELT(idflist,0,idvec);
	SET_VECTOR_ELT(idflist,1,geneIdVec);
	SET_VECTOR_ELT(idflist,2,transIdVec);

	// Column Names
	SEXP icol_names;
	PROTECT(icol_names=allocVector(STRSXP,nCols));
	++nProtected;

	SET_STRING_ELT(icol_names,0,mkChar("id"));
	SET_STRING_ELT(icol_names,1,mkChar("gene_id"));
	SET_STRING_ELT(icol_names,2,mkChar("transcript_id"));
	setAttrib(idflist,R_NamesSymbol,icol_names);

	// Row Names
	SEXP irow_names;
	PROTECT(irow_names=allocVector(STRSXP,nRows));
	++nProtected;

	for(i=0;i<nRows;++i)
	{
		sprintf(buf,"%lu",i+1);
		SET_STRING_ELT(irow_names,i,mkChar(buf));
	}
	setAttrib(idflist,R_RowNamesSymbol,irow_names);
	setAttrib(idflist,R_ClassSymbol,mkString("data.frame"));

	// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +
	// Merge id, type, value into adflist data.frame

	nCols=3;
	nRows=l->size;

	// Column 0: id
	SEXP aid_vector;
	PROTECT(aid_vector=allocVector(INTSXP,nRows));
	++nProtected;

	// Column 1: type
	SEXP atype_vector;
	PROTECT(atype_vector=allocVector(STRSXP,nRows));
	++nProtected;

	// Column 2: value
	SEXP aval_vector;
	PROTECT(aval_vector=allocVector(STRSXP,nRows));
	++nProtected;

	ptr_pair_list_rewind(l);
	const ptr_pair_element *e;

	for(i=0;i<nRows;++i)
	{
		e=ptr_pair_list_get_next_element(l);
		INTEGER(aid_vector)[i]=e->id;
		SET_STRING_ELT(atype_vector,i,Rf_mkCharLen(e->first,e->first_len));
		SET_STRING_ELT(aval_vector,i,Rf_mkCharLen(e->second,e->second_len));
	}

	SEXP adflist;
	PROTECT(adflist=allocVector(VECSXP,nCols));
	++nProtected;

	SET_VECTOR_ELT(adflist,0,aid_vector);
	SET_VECTOR_ELT(adflist,1,atype_vector);
	SET_VECTOR_ELT(adflist,2,aval_vector);

	// Column names
	SEXP acol_names;
	PROTECT(acol_names=allocVector(STRSXP,nCols));
	++nProtected;
	SET_STRING_ELT(acol_names,0,mkChar("id"));
	SET_STRING_ELT(acol_names,1,mkChar("type"));
	SET_STRING_ELT(acol_names,2,mkChar("value"));
	setAttrib(adflist,R_NamesSymbol,acol_names);

	// Row names
	SEXP arow_names;
	PROTECT(arow_names=allocVector(STRSXP,nRows));
	++nProtected;
	for(i=0;i<nRows;++i)
	{
		sprintf(buf,"%lu",i+1);
		SET_STRING_ELT(arow_names,i,mkChar(buf));
	}
	setAttrib(adflist,R_RowNamesSymbol,arow_names);

	// Class symbol
	setAttrib(adflist,R_ClassSymbol,mkString("data.frame"));

	// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +
	//  ans result list
	int list_len=2;
	SEXP ans;
	PROTECT(ans=allocVector(VECSXP,list_len));
	++nProtected;

	SET_VECTOR_ELT(ans,0,idflist);
	SET_VECTOR_ELT(ans,1,adflist);

	// names
	SEXP ans_names;
	PROTECT(ans_names=allocVector(STRSXP,list_len));
	++nProtected;
	SET_STRING_ELT(ans_names,0,mkChar("fixed"));
	SET_STRING_ELT(ans_names,1,mkChar("variable"));
	setAttrib(ans,R_NamesSymbol,ans_names);

	// Class symbol
	setAttrib(ans,R_ClassSymbol,mkString("list"));

	// Return
	ptr_pair_list_destroy(l);
	UNPROTECT(nProtected);
	return ans;
}

#endif	/* REFGENOME_C_ */

