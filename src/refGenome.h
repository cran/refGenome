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
 * 07.May.13    :   Valgrind tested package examples
 * 08.May.13    :   Corecced inline to static R_INLINE (in ptr_pair_list.h)
 *
 **************************************************************************************************
 **************************************************************************************************/

#ifndef REFGENOME_H_
#define REFGENOME_H_
#include <string.h>
#include <unistd.h>
#include <stdio.h>
#include <ctype.h>
#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/PrtUtil.h>
#include "ptr_pair_list.h"

// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +
//
// Split Ensembl gtf Attribute column data
// and return data.frame
//
// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +

SEXP split_gtf_attr(SEXP id_vec,SEXP attr_vec);

// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +
//
// Calculate exon_number from subsequent transcript and start values
// Expects ordering by transcript, seqid, start, end
//
// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +

SEXP get_exon_number(SEXP pTranscript,SEXP pSeqid, SEXP pStart, SEXP pEnd);
SEXP get_splice_juncs(SEXP pTranscript,SEXP pId,SEXP pStart,SEXP pEnd);
SEXP unify_splice_juncs(SEXP pSeqid,SEXP pLstart,SEXP pLend,SEXP pRstart,SEXP pRend, SEXP pId, SEXP pGeneId, SEXP pStrand);



// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +
//
// overlap: Compares list of query and reference ranges
// and reports overlaps in data.frame
//
// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +

#define OVERLAP_NO_OVER 1
#define OVERLAP_R_OVER  2
#define OVERLAP_B_OVER  3
#define OVERLAP_N_OVER  4
#define OVERLAP_L_OVER  5

SEXP overlap_ranges(SEXP qryid, SEXP qrystart, SEXP qryend, SEXP refid,SEXP refstart,SEXP refend);


#endif /* REFGENOME_H_ */
