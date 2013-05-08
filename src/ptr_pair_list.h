/*
 * ptr_ptr_pair_list.h
 *
 *  Created on: 16.10.2012
 *      Author: kaisers
 */

#ifndef PTR_ptr_pair_list_H_
#define PTR_ptr_pair_list_H_

#include <R.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define PTRPL_FIRST_TOKEN_EOF			-1
#define PTRPL_FIRST_TOKEN_EMPTY			-2
#define PTRPL_SECOND_TOKEN_DELIM_MISS	-3

///////////////////////////////////////////////////////////////////////////////////////////////////
// basic definitions

typedef struct ptr_pair_element
{
	unsigned long id;
	const char *first;
	unsigned first_len;
	const char *second;
	unsigned second_len;
	struct ptr_pair_element *next_el;
} ptr_pair_element;

typedef struct ptr_pair_list
{
	ptr_pair_element *first_el;
	ptr_pair_element *last_el;
	ptr_pair_element *curr_el;
	unsigned long size;
} ptr_pair_list;

///////////////////////////////////////////////////////////////////////////////////////////////////
// basic functions

static R_INLINE ptr_pair_list * ptr_pair_list_init()	{ return (ptr_pair_list*) calloc(sizeof(ptr_pair_list),1); }

static R_INLINE ptr_pair_element* ptr_pair_elem_init(const char *first, const unsigned first_len, const char *second, const unsigned second_len, unsigned long id)
{
	// Nothing is copied.
	// The list just contains pointer into target string
	ptr_pair_element *e =malloc(sizeof(ptr_pair_element));
	e->first=first;
	e->first_len=first_len;
	e->second=second;
	e->second_len=second_len;
	e->id=id;
	return e;
}

// Nothing has been copied, so nothing has to be freed
static R_INLINE void ptr_pair_elem_destroy(ptr_pair_element *e) {free(e);}

///////////////////////////////////////////////////////////////////////////////////////////////////
// list generic accessor functions

static R_INLINE void ptr_pair_list_push_back(ptr_pair_list *l, const char* first, const unsigned first_len, const char* second,const unsigned second_len, unsigned long id)
{
	ptr_pair_element *el=ptr_pair_elem_init(first,first_len,second,second_len,id);
	if(l->size==0)
	{
		l->first_el=el;
		l->last_el=el;
		l->size=1;
	}
	else
	{
		l->last_el->next_el=el;
		l->last_el=el;
		++(l->size);
	}
}

static R_INLINE void ptr_pair_list_push_front(ptr_pair_list *l, const char *first,const unsigned first_len, const char *second, const unsigned second_len, unsigned long id)
{
	ptr_pair_element *el=ptr_pair_elem_init(first,first_len,second,second_len,id);
	if(l->size==0)
	{
		l->first_el=el;
		l->last_el=el;
		l->size=1;
	}
	else
	{
		el->next_el=l->first_el;
		l->first_el=el;
		++(l->size);
	}
}

void ptr_pair_list_pop_front(ptr_pair_list *l)
{
	ptr_pair_element *e=l->first_el;
	if(l->size>1)
		l->first_el=l->first_el->next_el;
	else
		l->first_el=0;
	ptr_pair_elem_destroy(e);
	--(l->size);
}

static R_INLINE void ptr_pair_list_rewind(ptr_pair_list *l)	{ l->curr_el=0; }


///////////////////////////////////////////////////////////////////////////////////////////////////
// higher level functions

void ptr_pair_list_destroy(ptr_pair_list *l)
{
	while(l->size>0)
		ptr_pair_list_pop_front(l);
	free(l);
}

static R_INLINE const ptr_pair_element * ptr_pair_list_get_next_element(ptr_pair_list *l)
{
	if(l->first_el==0)
		return (ptr_pair_element*) 0;
	if(l->curr_el==0)
	{
		l->curr_el=l->first_el;
		return l->curr_el;
	}
	if(l->curr_el->next_el==0)
	{
		l->curr_el=0;
		return (ptr_pair_element*)0;
	}
	l->curr_el=l->curr_el->next_el;
	return l->curr_el;
}


///////////////////////////////////////////////////////////////////////////////////////////////////
// String splitting functions
///////////////////////////////////////////////////////////////////////////////////////////////////


int ptr_pair_list_push_back_gtf_attr(ptr_pair_list *l, const char* gtf_attr, unsigned long id)
{
	// flag characters
	char space=' ',delim=';',quote='"',zero ='\0';
	// iterator positions
	const char *token_first,*token_second,*iter;
	// string length
	unsigned long first_len,second_len;

	iter=gtf_attr;
	while(*iter!=zero)
	{
		// ++++++++++++++++++++++++++++++++++++++++++++++++
		// First token
		// skip spaces, then take start position
		while((*iter==space))
			++iter;
		token_first=iter;

		// proceed until space and take length
		while((*iter!=space) && (*iter!=zero))
			++iter;
		if(*iter==zero)
			return PTRPL_FIRST_TOKEN_EOF;
		if(iter==token_first)
			return PTRPL_FIRST_TOKEN_EMPTY;
		first_len=iter-token_first;

		// ++++++++++++++++++++++++++++++++++++++++++++++++
		// second token:
		// skip spaces and quotes, then take start position
		while((*iter==space) || (*iter==quote))
			++iter;
		token_second=iter;

		// proceed until space or quote
		while((*iter!=space) && (*iter!=quote) && (*iter!=delim) && (*iter!=zero))
			++iter;
		second_len=iter-token_second;

		// second token may be closed by quote
		if(*iter==quote)
			++iter;

		if(*iter!=delim)
			return PTRPL_SECOND_TOKEN_DELIM_MISS;
		++iter;
		ptr_pair_list_push_back(l,token_first,first_len,token_second,second_len,id);
	}
	return 0;
}

#endif /* PTR_ptr_pair_list_H_ */
