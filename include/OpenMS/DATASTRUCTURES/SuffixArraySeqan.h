// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Clemens Groepl,Andreas Bertsch$
// --------------------------------------------------------------------------



#ifndef OPENMS_DATASTRUCTURES_SUFFIXARRAYSEQAN_H
#define OPENMS_DATASTRUCTURES_SUFFIXARRAYSEQAN_H

#include <vector>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <seqan/index.h>
#include <OpenMS/DATASTRUCTURES/SuffixArray.h>

namespace OpenMS
{

	/**
	@brief Class that uses SEQAN library for a suffix array. It can be used to find peptide Candidates for a MS spectrum

	This class uses SEQAN suffix array. It can just be used for finding peptide Candidates for a given MS Spectrum within a certain mass tolerance. The suffix array can be saved to disc for reused so it has to be build just once.

	*/

	class OPENMS_DLLAPI SuffixArraySeqan : public SuffixArray
	{

		typedef seqan::TopDown<seqan::ParentLinks<> > TIterSpec;
		typedef seqan::Index<seqan::String<char>, seqan::Index_ESA<TIterSpec> > TIndex;
		typedef seqan::Iter<TIndex, seqan::VSTree<TIterSpec> > TIter;

		// TODO ??? was: typedef seqan::Index<seqan::String<char>, seqan::Index_ESA<seqan::TopDown<seqan::ParentLinks<seqan::Preorder> > > > TIndex;

	 public:

		/**
		@brief constructor
		@param st const string reference with the string for which the suffix array should be build
		@param filename const string reference with filename for opening or saving the suffix array
		@throw FileNotFound is thrown if the given file is not found
		@throw InvalidValue if the given suffix array string is invalid
		*/
		SuffixArraySeqan(const String& st, const String& filename);

		/**
		@brief copy constructor
		*/
		SuffixArraySeqan(const SuffixArraySeqan& source);

		/**
		@brief destructor
		*/
		virtual ~SuffixArraySeqan();

		/**
		@brief converts suffix array to a printable string
		*/
		String toString();

		/**
		@brief the function that will find all peptide candidates for a given spectrum
		@param spec const reference of double vector describing the spectrum
		@param candidates output parameters which holds the candidates of the masses given in spec after call
		@return a vector of int pairs.

		for every mass within the spectrum all candidates described by as pairs of ints are returned. All masses are searched for the same time in just one suffix array traversal. In order to accelerate the traversal the skip and lcp table are used. The mass wont be calculated for each entry but it will be updated during traversal using a stack datastructure
		*/
		void findSpec(std::vector<std::vector<std::pair<std::pair<int, int>, double > > >& candidates, const std::vector<double> & spec);

		/**
		@brief saves the suffix array to disc
		@param filename const reference string describing the filename
		@return bool if operation was succesful
		@throw UnableToCreateFile is thrown if the output files could not be created
		*/
		bool save(const String& filename);

		/**
		@brief opens the suffix array

		@param filename const reference string describing the filename
		@return bool if operation was succesful
		@throw FileNotFound is thrown if the given file could not be found
		*/
		bool open(const String& filename);

		/**
		@brief setter for tolerance
		@param t double with tolerance, only 0 or greater is allowed
		@throw InvalidValue is thrown if given tolerance is negative
		*/
		void setTolerance(double t);

		/**
		@brief getter for tolerance
		@return double with tolerance
		*/
		double getTolerance() const;

		/**
		@brief returns if an enzyme will cut after first character
		@param aa1 const char as first aminoacid
		@param aa2 const char as second aminoacid
		@return bool descibing if it is a digesting site
		*/
		bool isDigestingEnd(const char aa1, const char aa2) const;

		/**
		@brief setter for tags
		@param tags reference to vector of strings with tags
		@note sets use_tags = true
		*/
		void setTags(const std::vector<OpenMS::String>& tags);

		/**
		@brief getter for tags
		@return const reference to vector of strings
		*/
		const std::vector<OpenMS::String>& getTags ();

		/**
		@brief setter for use_tags
		@param use_tags indicating whether tags should be used or not
		*/
		void setUseTags(bool use_tags);

		/**
		@brief getter for use_tags
		@return bool indicating whether tags are used or not
		*/
		bool getUseTags();

		/**
		@brief setter for number of modifications
		@param number_of_mods
		*/
		void setNumberOfModifications(unsigned int number_of_mods);

		/**
		@brief getter for number of modifications
		@return number of modifications
		*/
		unsigned int getNumberOfModifications();

		void printStatistic ();

	 protected:

		/**
		@brief overwriting goNextSubTree_ from seqan index_esa_stree.h for mass update during suffix array traversal

		the suffix array is treated as a suffix tree. this function skips the subtree under the actual node and goes directly to the next subtree that has not been visited yet. During this traversal the mass will be updated using the stack with edge masses.

		@param it reference to the suffix array iterator
		@param m reference to actual mass
		@param allm reference to the stack with history of traversal
		@param mod_map input parameters which specifies the modification massen allowed in the candidates

		@see goNext
		*/
		inline void goNextSubTree_(TIter& it, double& m, std::stack<double>& allm, std::stack<std::map<double, int> >& mod_map)
		{
			// preorder dfs
			if (!goRight(it))
			{
				while (true)
				{
					if (goUp(it))
					{
						m -= allm.top();
						allm.pop();
						mod_map.pop();
					}
					else
					{
						break;
					}

					if (goRight(it))
					{
						m -= allm.top();
						allm.pop();
						mod_map.pop();
						break;
					}
				}
			}
			else
			{
				m -= allm.top();
				allm.pop();
				mod_map.pop();
			}
			if (isRoot(it))
			{
				clear(it);
			}
		}

		/**
		@brief goes to the next sub tree
		@param it reference to the suffix array iterator
		@see goNext
		*/
		inline void goNextSubTree_(TIter& it)
		{
			// preorder dfs
			if (!goRight(it))
			{
				while (true)
				{
					if (!goUp(it))
					{
						break;
					}
					if (goRight(it))
					{
						break;
					}
				}
			}
			if (isRoot(it))
			{
				clear(it);
			}
		}

		/**
		@brief overwriting goNext from seqan index_esa_stree.h for mass update during suffix array traversal

		the suffix array is treated as a suffix tree. this function goes to the next node that has not been visited yet. During this traversal the mass will be updated using the stack with edge masses.

		@param it reference to the suffix array iterator
		@param m reference to actual mass
		@param allm reference to the stack with history of traversal
		@param mod_map input parameters which specifies the modification masses allowed in the candidates

		@see goNextSubTree_
		*/
		inline void goNext_(TIter& it, double& m, std::stack<double>& allm, std::stack<std::map<double, int> >& mod_map)
		{
			// preorder dfs
			if (!goDown(it))
			{
				goNextSubTree_(it, m, allm, mod_map);
			}
		}



		inline void parseTree_(TIter& it, std::vector<std::pair<int, int> >& out_number, std::vector<std::pair<int, int> >& edge_length, std::vector<int>& leafe_depth)
		{
			int depth = 1;
			while (!atEnd(it))
			{
				int le = 0;
				bool isLeaf = false;
				if (length(parentEdgeLabel(it))>0){
					if (countChildren(it)>0)
					{
						edge_length.push_back(std::pair<int,int>(depth,length(parentEdgeLabel(it))));
					} else
					{
						//le <- length(representative(it));
						//isLeaf = true;
					}
				}
				if (countChildren(it)>0) {
					out_number.push_back(std::pair<int,int> (depth,countChildren(it)));
				} else {
					leafe_depth.push_back(depth);
				}
				if (goDown(it)){
					depth++;
				} else if (!goRight(it)) {
					while(!goRight(it)) {
						goUp(it);
						if (isLeaf) {
							edge_length.push_back(std::pair<int,int>(depth,le - length(parentEdgeLabel(it))));
							isLeaf = false;
						}
						depth--;
						if (isRoot(it)) return;
					}
				}
				else
				{
				}
			}
		}


		TIndex index_; ///< seqan suffix array

		TIter* it_; ///< seqan suffix array iterator

		/**
		@brief binary search for finding the index of the first element of the spectrum that matches the desired mass within the tolerance.
		@param spec const reference to spectrum
		@param m mass
		@return int with the index of the first occurence
		@note requires that there is at least one occurence
		*/
		int findFirst_ (const std::vector<double> & spec, double & m);

		/**
		@brief binary search for finding the index of the first element of the spectrum that matches the desired mass within the tolerance. it searches recursivly.
		@param spec const reference to spectrum
		@param m mass
		@param start start index
		@param end end index
		@return int with the index of the first occurence
		@note requires that there is at least one occurence
		*/
		int findFirst_ (const std::vector<double> & spec, double & m,int start, int  end);

		const String& s_; ///< reference to strings for which the suffix array is build

		double masse_[255]; ///< amino acid masses

		int number_of_modifications_; ///< number of allowed modifications

		std::vector<String> tags_; ///< all tags

		bool use_tags_; ///< if tags are used

		double tol_; ///< tolerance
	};
}

#endif //OPENMS_DATASTRUCTURES_SUFFIXARRAYSEQAN_H
