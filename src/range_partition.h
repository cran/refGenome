/*
 * tiled_schema.hpp
 *
 *  Created on: 27.01.2015
 *      Author: kaisers
 */

#ifndef RANGE_PARTITION_HPP_
#define RANGE_PARTITION_HPP_

//============================================================================
// Name        : range_partition.hpp
// Author      : W. Kaisers
// Version     : Segmenting schema a genomic range
// Date        : 12.01.2015
//============================================================================


#include <iostream>
#include <vector>
#include <list>
#include <iomanip>
using namespace std;



// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Represents a tiling schema on a (genomic) region
// (on single seqid's (i.e. chromosomes))
//
// General rule (optional ??): No inserts outside segment range.


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //

typedef unsigned position_type; // position type
typedef long int cat_type; // category type

struct range_node;
struct range_element;
class range_list;
class partition;
class range_inserter;



// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Data structs:
// range_node
// range_list
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //

struct range_node;
struct range_element;

struct range_node
{
	range_node(): id(1), position(0), category(0), source(0), source_id(0) {}
	cat_type id;
	position_type position;
	cat_type category;
	cat_type source;
	cat_type source_id;

	bool in(const range_element &elem) const;
	bool operator < (const range_element &elem) const;
	bool operator > (const range_element &elem) const;
	void operator++() { ++id; ++source_id; }

	// range_node
	bool operator < (const range_node &rhs) const { return position < rhs.position; }
	bool operator > (const range_node &rhs) const { return position > rhs.position; }
	bool operator == (const range_node &rhs) const { return position == rhs.position; }
	bool operator != (const range_node &rhs) const { return position != rhs.position; }

	// position_type
	bool operator < (const position_type pos) const { return position < pos; }
	bool operator > (const position_type pos) const { return position > pos; }
	bool operator == (const position_type pos) const { return position == pos; }
	bool operator != (const position_type pos) const { return position != pos; }
	operator position_type() const { return position; }
};

struct range_element
{
	cat_type id;
	position_type begin;
	position_type end;
	cat_type category;
	cat_type source;
	cat_type source_id;

	bool operator > (const range_node &node) const
	{ return node.position > end; }

	bool operator < (const range_node &node) const
	{ return node.position < begin; }

	bool operator > (const partition & scheme) const;
	bool operator < (const partition & scheme) const;
	bool in(const partition &scheme) const;

	operator range_node () const
	{
		range_node node;
		node.id = id;
		node.position = begin;
		node.category = category;
		node.source = source;
		node.source_id = id;
		return node;
	}

};

// ToDo: Shift back to partition_stream.hpp
ostream & operator<< (ostream &os, const range_node &node)
{
	os << "[node] id: " << setw(3) << node.id;
	os << "\tpos: " << setw(4) << node.position;
	os << "\tsrc_id: " << setw(3) << node.source_id;
	return os;
}
ostream & operator<<(ostream &os, const range_element &obj)
{
	os << "[range_elem] id: "	<< setw(3) << obj.id;
	os << "\tbegin: " 			<< setw(4) << obj.begin;
	os << "\tend: " 			<< setw(4) << obj.end;
	return os;
}


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Comparison operators between both element classes
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //

bool range_node::operator < (const range_element &elem) const
{ return position < elem.begin; }

bool range_node::operator > (const range_element &elem) const
{ return position >= elem.end; }

bool range_node::in(const range_element &elem) const
{ return (position >= elem.begin) && (position < elem.end); }




// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// partition class
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //

class partition
{
public:

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
	// Constructors
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
	partition(): begin_(0), end_(0), last_id_(0),
		category_(0), distance_(1), source_(0) {}
	partition(position_type begin, position_type end, position_type dist,
			cat_type category, position_type source):
		begin_(begin), end_(end), last_id_(0),
		category_(category), distance_(dist), source_(source)
	{
		range_node node;
		while(getNextNode(node))
		{
			node.id = ++last_id_;
			l_.push_back(node);
		}
	}

	partition(const range_list &rlist, position_type dist=1000, position_type begin=0, position_type end=0);

	virtual ~partition() {}

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
	// Accessors and friends
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //

	size_t size() const { return l_.size(); }
	position_type begin() const { return begin_; }
	position_type end() const { return end_; }
	bool setRange (position_type begin, position_type end);

	friend ostream & operator<< (ostream &os, const partition & s);
	friend class range_list;

	// ToDo: Remove this DEBUG - only function
	list<range_node>::iterator begin_iter() { return l_.begin(); }
	list<range_node>::iterator end_iter() { return l_.end(); }


public:
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
	// Functions for pushing and popping list elements
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
	void push_front(range_node & node, list<range_node>::iterator &iter);
	void push_front(range_node & node);

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
	// push_front +
	// additionally front-pushes node
	// so that fist node points to position
	//
	// Intended as preparation for adding
	// nodes at the beginning of the list
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
	bool set_start_node(position_type position);
	void push_front(list<range_element>::const_iterator &rit, list<range_node>::iterator &wit);

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
	// Removes all nodes < position
	// and up-shifts first node to position
	//
	// Does not ensure that first node
	// points to 'position' (i.e. when l_.begin() > position)
	//
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
	bool pop_front(position_type position);


	void push_back(range_node &node);
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
	// pop_back + additionally back-pushes node
	// so that last node points to position
	//
	// Intended as preparation for adding nodes at the end of the list
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
	bool set_end_node(position_type position);

	// ToDo: Does not check for back position...
	void push_back(list<range_element>::const_iterator &iter);

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
	// Removes all nodes > position
	// and down-shifts first node to position
	//
	// Does not ensure that last node
	// points to 'position' (i.e. when l_.back() < position)
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
	bool pop_back(position_type position);

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
	// Insertion of a range_element inside
	// Pre - condition:
	// list<range_node>::iterator begin
	// at position of or behind last insertion
	// (e.g. returned by push_front function)
	//
	// Returns l_.end() when range_element has to be
	// push_back'ed
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
	void prepare_inner_insert(const range_element &elem, list<range_node>::iterator &iter);
	void inner_insert(list<range_element>::const_iterator &rit, list<range_node>::iterator &wit);
	void insert(range_list &l);


	bool getNextNode(range_node &node);


private:
	// Border of range
	position_type begin_;
	position_type end_;
	cat_type last_id_;
	cat_type category_;
	position_type distance_;
	position_type source_;

	// List containing position nodes
	list<range_node> l_;

}; // partition


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Accessors and friends
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //

bool partition::setRange (position_type begin, position_type end)
{
	if((begin > l_.begin()->position) || (end < l_.back().position))
		return false;

	begin_ = begin;
	end_ = end;
	return true;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Functions for pushing and popping list elements
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //

void partition::push_front(range_node & node, list<range_node>::iterator &iter)
{
	iter = l_.begin();
	if(l_.size())
	{
		if(node.position >= l_.begin()->position)
			return;
	}

	node.id = ++last_id_;
	l_.push_front(node);
	++iter;
	return;
}

void partition::push_front(range_node & node)
{
	list<range_node>::iterator iter = l_.begin();
	if(l_.size())
	{
		if(node.position >= l_.begin()->position)
			return;
	}

	node.id = ++last_id_;
	l_.push_front(node);
	++iter;
	return;
}


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// push_front +
// additionally front-pushes node
// so that fist node points to position
//
// Intended as preparation for adding
// nodes at the beginning of the list
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //

bool partition::set_start_node(position_type position)
{
	if(!l_.size())
	{
		range_node node;
		node.position = position;
		push_back(node);
		return true;
	}

	pop_front(position);

	if(l_.begin()->position > position)
	{
		range_node node = *l_.begin();
		node.position = position;

		// Sets node.id
		push_front(node);
	}
	return true;
}


void partition::push_front(list<range_element>::const_iterator &rit, list<range_node>::iterator &wit)
{
	//cout << "[push_front] beg: " << setw(5) <<  rit->begin << " end: " << setw(5) << rit->end << "\n";

	wit = l_.begin();
	if(rit->begin > wit->position)
		return;

	// Eventually removes or shifts nodes at front:
	set_start_node(rit->end);

	// Translate into node and use node insertion (manages node-id)
	range_node node = *rit;
	push_front(node, wit);

	// Advance read iter
	++rit;
	return;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Removes all nodes < position
// and up-shifts first node to position
//
// Does not ensure that first node
// points to 'position' (i.e. when l_.begin() > position)
//
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
bool partition::pop_front(position_type position)
{
	if(!l_.size())
		return true;

	// There is nothing to remove
	if(l_.begin()->position >= position)
		return true;

	// Just remove a single element
	if(l_.size() == 1)
		l_.pop_front();

	// Do truncate using two iterators:
	// second points to second element in list
	list<range_node>::const_iterator second = l_.begin();
	++second;

	while( (second->position < position) && (second != l_.end()))
	{

		l_.pop_front();
		++second;
	}

	// Eventually up-shift first element to position
	l_.begin()->position = position;

	return true;
}


void partition::push_back(range_node &node)
{
	if(l_.size())
	{
		if(node.position <= l_.back().position)
			return;
	}

	// ToDo: Check for begin_ and end_; Throw exception?
	node.id = ++last_id_;
	l_.push_back(node);
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// pop_back +
// additionally back-pushes node
// so that last node points to position
//
// Intended as preparation for adding
// nodes at the end of the list
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //

bool partition::set_end_node(position_type position)
{
	if(!l_.size())
	{
		range_node node;
		node.position = position;
		push_back(node);
		return true;
	}

	pop_back(position);

	if(l_.back().position < position)
	{
		range_node node = l_.back();
		node.position = position;

		// Sets node.id
		push_back(node);
	}
	return true;
}


// ToDo: Does not check for back position...
void partition::push_back(list<range_element>::const_iterator &iter)
{
	// Removes or shifts nodes.position > iter->begin
	pop_back(iter->begin);

	range_node node = *iter;
	push_back(node);
	range_node empty;

	// ToDo: elem.end + 1 ??
	empty.position = iter->end;
	push_back(empty);

	// Advance iter behind insert
	++iter;

	//cout << "[partition] push_back id: " << node.id << "\n";
}


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Removes all nodes > position
// and down-shifts first node to position
//
// Does not ensure that last node
// points to 'position' (i.e. when l_.back() < position)
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
bool partition::pop_back(position_type position)
{
	if(!l_.size())
		return true;

	// There is nothing to remove
	if(l_.back().position < position)
		return true;

	// Just remove a single element
	if(l_.size() == 1)
		l_.pop_back();

	// Do truncate using two iterators:
	// second points to second to last element
	list<range_node>::const_iterator second = l_.end();
	--second;
	--second;
	while( (second->position > position) && (second != l_.begin()) )
	{
		l_.pop_back();
		--second;
	}

	// Eventually up-shift first element to position
	l_.back().position= position;
	return true;
}

// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Insertion of a range_element inside
// Pre - condition:
// list<range_node>::iterator begin
// at position of or behind last insertion
// (e.g. returned by push_front function)
//
// Returns l_.end() when range_element has to be
// push_back'ed
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
void partition::prepare_inner_insert(const range_element &elem, list<range_node>::iterator &iter)
{
	range_node node = *iter;

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
	// A) Advance iter until elem.begin is passed.
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
	while(iter->position < elem.begin && iter != l_.end())
		++iter;

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
	// B) Advance iter until elem.end is passed
	// B.1) Remember position			: begin
	// B.2) Count nodes inside range 	: nnodes
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
	list<range_node>::iterator begin;
	begin = iter;
	// Number of iteration steps between begin and current iter position
	unsigned nnodes=0;

	while(iter->position < elem.end && iter != l_.end())
	{
		++iter;
		++nnodes;
	}

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
	// C) Eventually add post-insert node when
	// C.1) elem.begin past all node.position behind iter -> iter advances to end
	// C.2) range is completely inside one node-interspace (nnodes = 0)
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
	if(iter == l_.end() || nnodes == 0)
	{
		// begin.pos > elem.end
		// -> Add post-insert node
		node.position = elem.end;
		node.id = ++last_id_;
		l_.insert(iter, node);
		// iter now is positioned on new node
		--iter;
	}

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
	// D) Eventually remove or shift nodes which will be hidden by range
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
	//cout << "[prep_ins] iter 2: " << *iter << "\tnnodes: " << nnodes << "\n";

	if(iter->position > elem.end && nnodes > 0)
	{
		--iter;
		iter->position = elem.end;
	}

	if(nnodes > 0)
	{
		//cout << "[ERASE] beg: " << begin->position << "\titer: " << iter->position << "\n";
		if(iter != begin)
			l_.erase(begin, iter);
	}


	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
	// Terminate
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
	return;
}

void partition::inner_insert(list<range_element>::const_iterator &rit, list<range_node>::iterator &wit)
{
	//cout << "[in_ins] pre  prep iter: " << *wit << "\n";
	prepare_inner_insert(*rit, wit);
	//cout << "[in_ins] post prep iter: " << *wit << "\n";

	if(wit != l_.end())
	{
		// ToDo: eventually check for empty list ??
		range_node node = *rit;
		node.id = ++last_id_;
		l_.insert(wit, node);

		// Set iter to new inserted node
		--wit;
		//cout << "[partition] inner_insert: id: " << node.id << "\n";
	}
	else
	{
		push_back(rit);
	}
	++rit;

	return;
}

bool partition::getNextNode(range_node &node)
{
	if(begin_ < end_)
	{
		node.id = ++last_id_;
		node.position = begin_;
		node.category = category_;
		node.source = source_;
		node.source_id = last_id_;

		begin_ += distance_;
		return true;
	}
	return false;
}


///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////



// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// range_list container
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //


class range_list
{
public:
	range_list(): begin_(0), end_(0) {}
	virtual ~range_list() {}

	size_t size() const { return l_.size(); }
	operator partition() const;

	// ToDo: Throw exception ?
	bool push_back(range_element &elem);
	bool push_font(range_element &elem);

	void setRange(position_type begin, position_type end) { begin_= begin; end_ = end; }

public:
	list<range_element>::iterator erase(list<range_element>::iterator first, list<range_element>::iterator last)
	{ return l_.erase(first, last);	}

	friend ostream & operator<< (ostream &os, const range_list & obj);
	friend class partition;

private:
	position_type begin_;
	position_type end_;

	list<range_element> l_;
};

bool range_list::push_back(range_element &elem)
{
	if(l_.size())
	{
		// Create an ordered list of ranges
		// ToDo: Throw exception?
		if(elem.begin < l_.back().end)
			return false;

		if(elem.end <= elem.begin)
			return false;
	}
	else
	{
		begin_ = elem.begin;
	}
	l_.push_back(elem);
	end_ = elem.end;
	return true;
}

bool range_list::push_font(range_element &elem)
{
	if(l_.size())
	{
		range_element &front = *l_.begin();
		if(elem.begin < elem.end && elem.end < front.begin)
		{
			l_.push_front(elem);
			begin_ = elem.begin;

			return true;
		}
		return false;
	}else
	{
		l_.push_front(elem);
		end_ = elem.end;
		return true;
	}
	return false;
}

range_list::operator partition() const
{
	partition s;
	range_node range_start, empty;

	if(!l_.size()) return s;


	s.setRange(begin_, end_);
	// left and right iter:
	list<range_element>::const_iterator left = l_.begin(), right;
	range_start.category = left->category;
	range_start.position = left->begin;
	range_start.source = left->source;
	range_start.source_id = left->id;

	// ToDo: throw exception?
	s.push_back(range_start);


	if(l_.size() == 1)
	{
		range_start.position = left->end;
		++range_start.source_id;
		s.push_back(range_start);
		return s;
	}

	// size >= 2:
	// >= 2 ranges: insert empty segments
	right = left;
	++right;

	for(; right!=l_.end(); ++right, ++left)
	{
		// eventually insert first positon after end as dummy.
		////cout << "right.begin: " << right->begin << "\tleft.end: " << left->end << "\n";
		if(right->begin > (left->end + 1))
		{
			empty.position = left->end + 1;
			//empty.source_id = left->id;
			s.push_back(empty);
		}

		// insert next start
		range_start.position = right->begin;
		range_start.source_id = right->id;
		s.push_back(range_start);
	}
	empty.position = left->end + 1;
	//empty.source_id = left->id;
	s.push_back(empty);


	return s;
}


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// range_element comparing operators
// (Used in subsequent function implementations)
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
bool range_element::in(const partition &scheme) const
{ return (this->begin >= scheme.begin() && this->end <= scheme.end()); }

bool range_element::operator > (const partition & scheme) const
{ return this->begin > scheme.end(); }

bool range_element::operator < (const partition & scheme) const
{ return this->end < scheme.begin(); }


void partition::insert(range_list &l)
{
	list<range_element>::const_iterator rit = l.l_.begin();
	list<range_node>::iterator wit;

	push_front(rit, wit);
	while(rit != l.l_.end() && wit != l_.end())
		inner_insert(rit, wit);


	while(rit != l.l_.end())
		push_back(rit);

	return;
}


// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
// Constructor: Uses range_list internals
// So has to be positioned behind range_list class declaration
// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - //
partition::partition(const range_list &rlist, position_type dist, position_type begin, position_type end):
			begin_(begin), end_(end), last_id_(0), category_(0), distance_(dist), source_(0)
{
	// Prevent infinite loop
	if(!distance_)
		return;

	if(!end)
	{
		begin_ = rlist.begin_;
		end_ = rlist.end_;
	}

	cat_type last_id = 1;
	position_type n_dist=0;

	// pnode = partition node
	range_node pnode;
	pnode.id = last_id;
	pnode.category = 0;
	pnode.source = 0;
	pnode.source_id = 0;
	pnode.position = begin_;

	list<range_element>::const_iterator iter = rlist.l_.begin();


	while(iter != rlist.l_.end())
	{
		// Pre insert leading nodes
		while(pnode.position < iter->begin)
		{
			l_.push_back(pnode);
			pnode.position += distance_;
			pnode.id = ++last_id;
		}

		// Insert node with begin position
		l_.push_back(*iter);
		l_.back().id = ++last_id;

		// Insert "end" node
		// End = last position of range + 1
		//     = begin of subsequent range
		pnode.position = iter->end;
		pnode.id = ++last_id;
		l_.push_back(pnode);
		++iter;

		// Calculate position for next partition node in schema
		n_dist = (pnode.position - begin_) / distance_ + 1;
		pnode.position = begin_ + distance_ * n_dist;
		pnode.id = ++last_id;
	}


	while(pnode.position <= end_)
	{
		l_.push_back(pnode);
		pnode.position += distance_;
		pnode.id = ++last_id;

	}
}


#endif /* RANGE_PARTITION_HPP_ */
