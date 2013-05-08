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
