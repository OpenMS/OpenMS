// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Clemens Groepl $
// $Authors: $
// --------------------------------------------------------------------------
//

#include <algorithm>
#include <functional>
#include <string>
#include <utility>
#include <vector>

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/KERNEL/ComparatorUtils.h>


///////////////////////////

using namespace std;
using namespace OpenMS;

///////////////////////////

struct IntStringPair
{
	IntStringPair ( int _first, string _second, int _id = 0 )
		: first(_first), second(_second), id(_id) {}
	int first;
	string second;
	int id;
};

struct
StrangeMixedLessOfIntAndString : std::binary_function < int, string, bool >
{
	bool
	operator () ( int const & i, string const & s ) const
	{
		return i < strtol ( s.c_str(), 0, 10 );
	}
};

struct IntStringPairLessFirst : binary_function < IntStringPair, IntStringPair, bool >
{
	bool operator () ( IntStringPair const & left, IntStringPair const & right ) const
	{ return left.first < right.first; }
};

struct IntStringPairLessSecond : binary_function < IntStringPair, IntStringPair, bool >
{
	bool operator () ( IntStringPair const & left, IntStringPair const & right ) const
	{ return left.second < right.second; }
};

/////////////////////////////////////////////////////////////

START_TEST(ComparatorUtils.h, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION(PointerComparator < less<int> >)
	int i = 1000, j = 1000000;
  TEST_EQUAL(less<int>().operator()(i,j),true)
	TEST_EQUAL(less<int>().operator()(j,i),false)
	TEST_EQUAL(less<int>().operator()(i,i),false)
	TEST_EQUAL(PointerComparator < less<int> >().operator()(&i,&j),true)
	TEST_EQUAL(PointerComparator < less<int> >().operator()(&j,&i),false)
	TEST_EQUAL(PointerComparator < less<int> >().operator()(&i,&i),false)
	swap(i,j);
  TEST_EQUAL(less<int>().operator()(j,i),true)
	TEST_EQUAL(less<int>().operator()(i,j),false)
	TEST_EQUAL(less<int>().operator()(i,i),false)
	TEST_EQUAL(PointerComparator < less<int> >().operator()(&j,&i),true)
	TEST_EQUAL(PointerComparator < less<int> >().operator()(&i,&j),false)
	TEST_EQUAL(PointerComparator < less<int> >().operator()(&i,&i),false)
END_SECTION

START_SECTION(PointerComparator < ReverseComparator < less<int> > >)
	int i = 1000, j = 1000000;
	TEST_EQUAL(PointerComparator < ReverseComparator < less<int> > >().operator()(&i,&j),false)
	TEST_EQUAL(PointerComparator < ReverseComparator < less<int> > >().operator()(&i,&i),false)
	TEST_EQUAL(PointerComparator < ReverseComparator < less<int> > >().operator()(&j,&i),true)
	swap(i,j);
	TEST_EQUAL(PointerComparator < ReverseComparator < less<int> > >().operator()(&j,&i),false)
	TEST_EQUAL(PointerComparator < ReverseComparator < less<int> > >().operator()(&i,&i),false)
	TEST_EQUAL(PointerComparator < ReverseComparator < less<int> > >().operator()(&i,&j),true)
END_SECTION

START_SECTION(pointerComparator(std::less<int>()))
	int i = 88, j = 99;
  TEST_EQUAL(pointerComparator(std::less<int>())(&i,&j),true)
  TEST_EQUAL(pointerComparator(std::less<int>())(&i,&i),false)
  TEST_EQUAL(pointerComparator(std::less<int>())(&j,&i),false)
END_SECTION

START_SECTION(ReverseComparator < PointerComparator < less<int> > >)
	int i = 1000, j = 1000000;
	TEST_EQUAL(ReverseComparator < PointerComparator < less<int> > >().operator()(&i,&j),false)
	TEST_EQUAL(ReverseComparator < PointerComparator < less<int> > >().operator()(&i,&i),false)
	TEST_EQUAL(ReverseComparator < PointerComparator < less<int> > >().operator()(&j,&i),true)
	swap(i,j);
	TEST_EQUAL(ReverseComparator < PointerComparator < less<int> > >().operator()(&j,&i),false)
	TEST_EQUAL(ReverseComparator < PointerComparator < less<int> > >().operator()(&i,&i),false)
	TEST_EQUAL(ReverseComparator < PointerComparator < less<int> > >().operator()(&i,&j),true)
END_SECTION

START_SECTION(ReverseComparator < PointerComparator < less<int> > >)
	int i = 1000, j = 1000000;
	TEST_EQUAL(ReverseComparator < PointerComparator < less<int> > >().operator()(&i,&j),false)
	TEST_EQUAL(ReverseComparator < PointerComparator < less<int> > >().operator()(&i,&i),false)
	TEST_EQUAL(ReverseComparator < PointerComparator < less<int> > >().operator()(&j,&i),true)
	swap(i,j);
	TEST_EQUAL(ReverseComparator < PointerComparator < less<int> > >().operator()(&j,&i),false)
	TEST_EQUAL(ReverseComparator < PointerComparator < less<int> > >().operator()(&i,&i),false)
	TEST_EQUAL(ReverseComparator < PointerComparator < less<int> > >().operator()(&i,&j),true)
END_SECTION

START_SECTION(typedefs in ReverseComparator<> and PointerComparator<>)
  int i = 1000;
  string s = "1000000";
  TEST_EQUAL(StrangeMixedLessOfIntAndString().operator()(i,s),true)
  TEST_EQUAL(StrangeMixedLessOfIntAndString().operator()(2000000,s),false)
  TEST_EQUAL(StrangeMixedLessOfIntAndString().operator()(i,"2000"),true)
  TEST_EQUAL(StrangeMixedLessOfIntAndString().operator()(2,"10"),true)
  TEST_EQUAL(StrangeMixedLessOfIntAndString().operator()(20,"10"),false)
  TEST_EQUAL(StrangeMixedLessOfIntAndString().operator()(20,"20"),false)
  TEST_EQUAL(typeid(StrangeMixedLessOfIntAndString::first_argument_type)==typeid(int),true)
  TEST_EQUAL(typeid(StrangeMixedLessOfIntAndString::second_argument_type)==typeid(string),true)
  TEST_EQUAL(typeid(ReverseComparator<StrangeMixedLessOfIntAndString>::first_argument_type)==typeid(string),true)
  TEST_EQUAL(typeid(ReverseComparator<StrangeMixedLessOfIntAndString>::second_argument_type)==typeid(int),true)
  TEST_EQUAL(typeid(PointerComparator<StrangeMixedLessOfIntAndString>::first_argument_type)==typeid(int*),true)
  TEST_EQUAL(typeid(PointerComparator<StrangeMixedLessOfIntAndString>::second_argument_type)==typeid(string*),true)
  TEST_EQUAL(typeid(ReverseComparator<PointerComparator<StrangeMixedLessOfIntAndString> >::first_argument_type)==typeid(string*),true)
  TEST_EQUAL(typeid(ReverseComparator<PointerComparator<StrangeMixedLessOfIntAndString> >::second_argument_type)==typeid(int*),true)
  TEST_EQUAL(typeid(PointerComparator<ReverseComparator<StrangeMixedLessOfIntAndString> >::first_argument_type)==typeid(string*),true)
  TEST_EQUAL(typeid(PointerComparator<ReverseComparator<StrangeMixedLessOfIntAndString> >::second_argument_type)==typeid(int*),true)
END_SECTION

START_SECTION(reverseComparator(std::less<int>()))
	int i = 88, j = 99;
  TEST_EQUAL(reverseComparator(std::less<int>())(i,j),false)
  TEST_EQUAL(reverseComparator(std::less<int>())(i,i),false)
  TEST_EQUAL(reverseComparator(std::less<int>())(j,i),true)
END_SECTION

START_SECTION(reverseComparator(pointerComparator(std::less<int>())))
	int i = 88, j = 99;
  TEST_EQUAL(reverseComparator(pointerComparator(std::less<int>()))(&i,&j),false)
  TEST_EQUAL(reverseComparator(pointerComparator(std::less<int>()))(&i,&i),false)
  TEST_EQUAL(reverseComparator(pointerComparator(std::less<int>()))(&j,&i),true)
END_SECTION

START_SECTION(pointerComparator(reverseComparator(std::less<int>())))
	int i = 88, j = 99;
  TEST_EQUAL(pointerComparator(reverseComparator(std::less<int>()))(&i,&j),false)
  TEST_EQUAL(pointerComparator(reverseComparator(std::less<int>()))(&i,&i),false)
  TEST_EQUAL(pointerComparator(reverseComparator(std::less<int>()))(&j,&i),true)
END_SECTION

START_SECTION(LexicographicComparator<>)
  vector < IntStringPair > seq;
  seq.push_back( IntStringPair ( 1, "a", 1 ) );
  seq.push_back( IntStringPair ( 1, "b", 2 ) );
  seq.push_back( IntStringPair ( 1, "c", 3 ) );
  seq.push_back( IntStringPair ( 2, "a", 4 ) );
  seq.push_back( IntStringPair ( 2, "b", 5 ) );
  seq.push_back( IntStringPair ( 2, "c", 6 ) );
  seq.push_back( IntStringPair ( 3, "a", 7 ) );
  seq.push_back( IntStringPair ( 3, "b", 8 ) );
  seq.push_back( IntStringPair ( 3, "c", 9 ) );
  for ( vector < IntStringPair >::iterator p = seq.begin(); p != seq.end(); ++p ) STATUS( p-> id << ":  " << p->first << ' ' << p->second );

  vector < IntStringPair * > seqp;
  for ( vector < IntStringPair >::iterator pos = seq.begin(); pos != seq.end(); seqp.push_back(&(*pos++)) ) ;
  for ( vector < IntStringPair * >::iterator p = seqp.begin(); p != seqp.end(); ++p ) STATUS( (*p)->id << ":  " << (*p)->first << ' ' << (*p)->second ) ;

  sort ( seq.begin(), seq.end(), IntStringPairLessSecond() );
  sort ( seq.begin(), seq.end(), IntStringPairLessFirst() );
  sort ( seq.begin(), seq.end(), ReverseComparator<IntStringPairLessSecond>() );
  sort ( seq.begin(), seq.end(), ReverseComparator<IntStringPairLessFirst>() );

  random_shuffle(seq.begin(),seq.end());
  STATUS("after random_shuffle:");
  for ( vector < IntStringPair >::iterator p = seq.begin(); p != seq.end(); ++p ) STATUS( p-> id << ":  " << p->first << ' ' << p->second );
  STATUS("Okay!");

  for ( unsigned loops = 3; loops; --loops )
  {
  	STATUS("\n\nremaining loops: " << loops<<"\n");

  	random_shuffle(seq.begin(),seq.end());
  	sort ( seq.begin(), seq.end(), LexicographicComparator<IntStringPairLessFirst,IntStringPairLessSecond>() );
  	{
  		int order [] = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
  		for ( unsigned i = 0; i < 9; ++i )
  			TEST_EQUAL(seq[i].id,order[i]);
  	}
  	STATUS("Okay!");

  	random_shuffle(seq.begin(),seq.end());
  	sort ( seq.begin(), seq.end(), LexicographicComparator<IntStringPairLessFirst,ReverseComparator<IntStringPairLessSecond> >() );
  	{
  		int order [] = { 3, 2, 1, 6, 5, 4, 9, 8, 7 };
  		for ( unsigned i = 0; i < 9; ++i )
  			TEST_EQUAL(seq[i].id,order[i]);
  	}
  	STATUS("Okay!");

  	random_shuffle(seq.begin(),seq.end());
  	sort ( seq.begin(), seq.end(), LexicographicComparator<ReverseComparator<IntStringPairLessFirst>,IntStringPairLessSecond>() );
  	{
  		int order [] = { 7, 8, 9, 4, 5, 6, 1, 2, 3 };
  		for ( unsigned i = 0; i < 9; ++i )
  			TEST_EQUAL(seq[i].id,order[i]);
  	}
  	STATUS("Okay!");

  	random_shuffle(seq.begin(),seq.end());
  	sort ( seq.begin(), seq.end(), LexicographicComparator<ReverseComparator<IntStringPairLessFirst>,ReverseComparator<IntStringPairLessSecond> >() );
  	{
  		int order [] = { 9, 8, 7, 6, 5, 4, 3, 2, 1 };
  		for ( unsigned i = 0; i < 9; ++i )
  			TEST_EQUAL(seq[i].id,order[i]);
  	}
  	STATUS("Okay!");

  	random_shuffle(seq.begin(),seq.end());
  	sort ( seq.begin(), seq.end(), ReverseComparator<LexicographicComparator<ReverseComparator<IntStringPairLessFirst>,ReverseComparator<IntStringPairLessSecond> > >() );
  	{
  		int order [] = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
  		for ( unsigned i = 0; i < 9; ++i )
  			TEST_EQUAL(seq[i].id,order[i]);
  	}
  	STATUS("Okay!");

  	random_shuffle(seqp.begin(),seqp.end());
  	sort ( seqp.begin(), seqp.end(), PointerComparator<ReverseComparator<LexicographicComparator<ReverseComparator<IntStringPairLessFirst>,ReverseComparator<IntStringPairLessSecond> > > >() );
  	{
  		int order [] = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
  		for ( unsigned i = 0; i < 9; ++i )
  			TEST_EQUAL(seqp[i]->id,order[i]);
  	}
  	STATUS("Okay!");

  	random_shuffle(seqp.begin(),seqp.end());
  	sort ( seqp.begin(), seqp.end(), ReverseComparator<LexicographicComparator<PointerComparator<ReverseComparator<IntStringPairLessFirst> >,ReverseComparator<PointerComparator<IntStringPairLessSecond> > > >() );
  	{
  		int order [] = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
  		for ( unsigned i = 0; i < 9; ++i )
  			TEST_EQUAL(seqp[i]->id,order[i]);
  	}
  	STATUS("Okay!");

  	random_shuffle(seq.begin(),seq.end());
  	sort ( seq.begin(), seq.end(), LexicographicComparator<IntStringPairLessSecond,IntStringPairLessFirst>() );
  	{
  		int order [] = { 1, 4, 7, 2, 5, 8, 3, 6, 9 };
  		for ( unsigned i = 0; i < 9; ++i )
  			TEST_EQUAL(seq[i].id,order[i]);
  	}
  	STATUS("Okay!");

  	random_shuffle(seq.begin(),seq.end());
  	sort ( seq.begin(), seq.end(), LexicographicComparator<IntStringPairLessSecond,ReverseComparator<IntStringPairLessFirst> >() );
  	{
  		int order [] = { 7, 4, 1, 8, 5, 2, 9, 6, 3 };
  		for ( unsigned i = 0; i < 9; ++i )
  			TEST_EQUAL(seq[i].id,order[i]);
  	}
  	STATUS("Okay!");

  	random_shuffle(seq.begin(),seq.end());
  	sort ( seq.begin(), seq.end(), LexicographicComparator<ReverseComparator<IntStringPairLessSecond>,IntStringPairLessFirst>() );
  	{
  		int order [] = { 3, 6, 9, 2, 5, 8, 1, 4, 7 };
  		for ( unsigned i = 0; i < 9; ++i )
  			TEST_EQUAL(seq[i].id,order[i]);
  	}
  	STATUS("Okay!");

  	random_shuffle(seq.begin(),seq.end());
  	sort ( seq.begin(), seq.end(), LexicographicComparator<ReverseComparator<IntStringPairLessSecond>,ReverseComparator<IntStringPairLessFirst> >() );
  	{
  		int order [] = { 9, 6, 3, 8, 5, 2, 7, 4, 1 };
  		for ( unsigned i = 0; i < 9; ++i )
  			TEST_EQUAL(seq[i].id,order[i]);
  	}
  	STATUS("Okay!");

  }

END_SECTION

START_SECTION(lexicographicComparator())

	// Note: the correctness has been extensively in the preceding test, see
	// START_SECTION(LexicographicComparator<>)
	// Here we only check if the template instantiation works.
	// The different line is commented below.

  vector < IntStringPair > seq;
  seq.push_back( IntStringPair ( 1, "a", 1 ) );
  seq.push_back( IntStringPair ( 1, "b", 2 ) );
  seq.push_back( IntStringPair ( 1, "c", 3 ) );
  seq.push_back( IntStringPair ( 2, "a", 4 ) );
  seq.push_back( IntStringPair ( 2, "b", 5 ) );
  seq.push_back( IntStringPair ( 2, "c", 6 ) );
  seq.push_back( IntStringPair ( 3, "a", 7 ) );
  seq.push_back( IntStringPair ( 3, "b", 8 ) );
  seq.push_back( IntStringPair ( 3, "c", 9 ) );
  for ( vector < IntStringPair >::iterator p = seq.begin(); p != seq.end(); ++p ) STATUS( p-> id << ":  " << p->first << ' ' << p->second );

  vector < IntStringPair * > seqp;
  for ( vector < IntStringPair >::iterator pos = seq.begin(); pos != seq.end(); seqp.push_back(&(*pos++)) ) ;
  for ( vector < IntStringPair * >::iterator p = seqp.begin(); p != seqp.end(); ++p ) STATUS( (*p)->id << ":  " << (*p)->first << ' ' << (*p)->second ) ;

  sort ( seq.begin(), seq.end(), IntStringPairLessSecond() );
  sort ( seq.begin(), seq.end(), IntStringPairLessFirst() );
  sort ( seq.begin(), seq.end(), ReverseComparator<IntStringPairLessSecond>() );
  sort ( seq.begin(), seq.end(), ReverseComparator<IntStringPairLessFirst>() );

  random_shuffle(seq.begin(),seq.end());
  STATUS("after random_shuffle:");
  for ( vector < IntStringPair >::iterator p = seq.begin(); p != seq.end(); ++p ) STATUS( p-> id << ":  " << p->first << ' ' << p->second );
  STATUS("Okay!");

  for ( unsigned loops = 1; loops; --loops )
  {
  	STATUS("remaining loops: " << loops);

  	random_shuffle(seq.begin(),seq.end());
		// Note how the next line differs from the preceding test...
  	sort ( seq.begin(), seq.end(), lexicographicComparator(IntStringPairLessFirst(),IntStringPairLessSecond()) );
  	{
  		int order [] = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
  		for ( unsigned i = 0; i < 9; ++i )
  			TEST_EQUAL(seq[i].id,order[i]);
  	}
  }

END_SECTION

START_SECTION(PairComaratorFirstElement())
  std::pair<int,int> i(4,6), j(-3,7);
  PairComparatorFirstElement<std::pair<int,int> > testcomp;
  TEST_EQUAL(testcomp(i,j),false);
  TEST_EQUAL(testcomp(j,i),true);
  j = std::make_pair(4,7);
  TEST_EQUAL(testcomp(i,j),false);
  TEST_EQUAL(testcomp(j,i),false);

END_SECTION

START_SECTION(PairComaratorSecondElement())
  std::pair<int,int> i(4,6), j(-3,7);
  PairComparatorSecondElement<std::pair<int,int> > testcomp;
  TEST_EQUAL(testcomp(i,j),true);
  j = std::make_pair(4,6);
  TEST_EQUAL(testcomp(i,j),false);
END_SECTION

START_SECTION(PairComparatorFirstElementMore())
  std::pair<int,int> i(4,6), j(-3,7);
  PairComparatorFirstElementMore<std::pair<int,int> > testcomp;
  TEST_EQUAL(testcomp(i,j),true);
  TEST_EQUAL(testcomp(j,i),false);
  j = std::make_pair(4,7);
  TEST_EQUAL(testcomp(i,j),false);
  TEST_EQUAL(testcomp(j,i),false);

END_SECTION

START_SECTION(PairComparatorSecondElementMore())
  std::pair<int,int> i(4,6), j(-3,7);
  PairComparatorSecondElementMore<std::pair<int,int> > testcomp;
  TEST_EQUAL(testcomp(i,j),false);
  j = std::make_pair(4,6);
  TEST_EQUAL(testcomp(i,j),false);
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
