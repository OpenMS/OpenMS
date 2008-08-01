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


#ifndef OPENMS_DATASTRUCTURES_SUFFIXARRAYTRYPTICCOMPRESSED_H
#define OPENMS_DATASTRUCTURES_SUFFIXARRAYTRYPTICCOMPRESSED_H

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/SuffixArray.h>


namespace OpenMS {
	class String;

/**
	@brief Class that implements a suffix array for a String. It can be used to find peptide Candidates for a MS spectrum

	This class implements a suffix array. It can just be used for finding peptide Candidates for a given MS Spectrum within a certain mass tolerance. The suffix array can be saved to disc for reused so it has to be build just once. The suffix array consits of a vector of pair of ints for every suffix, a vector of LCP values and a so called skip vector.
	Only the sufices that are matching the function isDigestingEnd are created. Besides a suffix will not reach till the end of the string but till the next occurence of the seperator ($). So only the interessting sufices will be saved. This will reduce the used space.
*/

class SuffixArrayTrypticCompressed : public SuffixArray {
	
public:

	/**
	@brief constructor taking the string and the filename for writing or reading
	@param st the string as const reference with which the suffix array will be build
	@param filename the filename for writing or reading the suffix array
	@throw Exception::InvalidValue if string does not start with empty string ($)
	@throw FileNotFound is thrown if the given file was not found

	The constructor checks if a suffix array with given filename (without file extension) exists or not. In the first case it will simple be loaded and otherwise it will be build. Bulding the suffix array consists of several steps. At first all indices for a digesting enzyme (defined by using function isDigestingEnd) are created as an vector of int pairs. After creating all relevant indices they are sorted and the lcp and skip vectors are created.
	*/
	SuffixArrayTrypticCompressed(const String& st, const String& filename);

	/**
	@brief copy constructor
	*/
	SuffixArrayTrypticCompressed(const SuffixArrayTrypticCompressed & sa);

	/**
	@brief destructor
	*/
	virtual ~SuffixArrayTrypticCompressed();

	/**
	@brief transforms suffix array to a printable String
	*/
	String toString();

	/**
	@brief the function that will find all peptide candidates for a given spectrum
	@param spec const reference of double vector describing the spectrum
	@param candidates output parameter which contains the candidates of the masses given in spec
	@return a vector of int pairs.
	@throw InvalidValue if the spectrum is not sorted ascendingly
	
	for every mass within the spectrum all candidates described by as pairs of ints are returned. All masses are searched for the same time in just one suffix array traversal. In order to accelerate the traversal the skip and lcp table are used. The mass wont be calculated for each entry but it will be updated during traversal using a stack datastructure 
	*/
	void findSpec(std::vector<std::vector<std::pair<std::pair<int, int>, float > > >& candidates, const std::vector<double> & spec);

	/**
	@brief saves the suffix array to disc
	@param file_name const reference string describing the filename
	@return bool if operation was succesful
	@throw Exception::UnableToCreateFile if file could not be created (e.x. if you have no rigths)
	*/
	bool save(const String& file_name);
	/**
	@brief opens the suffix array
	@param file_name const reference string describing the filename
	@return bool if operation was succesful
	@throw FileNotFound
	*/
	bool open(const String& file_name);

	/**
	@brief setter for tolerance
	@param t double with tolerance
	@throw Exception::InvalidValue if tolerance is negative
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
	@param tags const vector of strings with tags with length 3 each
	@throw InvalidValue if at least one tag does not have size of 3
	*/
	void setTags(const std::vector<String>& tags);

	/**
	@brief getter for tags
	@return const vector of string with tags
	*/
	const std::vector<String>& getTags ();

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
	@return unsigned int describing number of modifications
	*/
	unsigned int getNumberOfModifications ();

	/**
	@brief output for statistic
	*/
	void printStatistic ();

protected:

	/**
	@brief constructor
	*/
	SuffixArrayTrypticCompressed();

	/**
	@brief gets the index of the next sperator for a given index
	@param p const int describing a position in the string
	@return int with the index of the next occurence of the sperator or -1 if there is no more seperator
	*/
	int getNextSep_(const int p) const;

	/**
	@brief gets the lcp for two strings described as pairs of ints
	@param last_point const pair of ints describing a substring
	@param current_point const pair of ints describing a substring
	@return int with the length of the lowest common prefix
	*/
	int getLCP_(const std::pair<int, int>& last_point, const std::pair<int, int>& current_point);

	/**
	@brief binary search for finding the index of the first element of the spectrum that matches the desired mass within the tolerance.
	@param spec const reference to spectrum
	@param m mass
	@return int with the index of the first occurence
	@note requires that there is at least one occurence
	*/
	int findFirst_(const std::vector<double>& spec, double& m);

	/**
	@brief binary search for finding the index of the first element of the spectrum that matches the desired mass within the tolerance. it searches recursivly.
	@param spec const reference to spectrum
	@param m mass
	@param start start index
	@param end end index
	@return int with the index of the first occurence
	@note requires that there is at least one occurence
	*/
	int findFirst_(const std::vector<double>& spec, double& m, int start, int end);

	/**
	@brief treats the suffix array as a tree and parses the tree using postorder traversion. This is realised by a recursive algorithm.
	@param start_index int describing the start index in indices_ vector
	@param stop_index int describing the end index in indices_ vector
	@param depth at with depth the traversion is at the actual position
	@param walked_in how many characters we have seen from root to actual position
	@param edge_len how many characters we have seen from last node to actual position
	@param out_number reference to vector of pairs of ints. For every node it will be filled with how many outgoing edge a node has in dependece of its depth
	@param edge_length will be filled with the edge_length in dependence of its depth
	@param leafe_depth will be filled with the depth of every leafe
	@note intialize: walked_in=0, depth=1, edge_len=1
	*/
	void parseTree_(int start_index, int stop_index, int depth, int walked_in, int edge_len, std::vector<std::pair<int,int> >& out_number, std::vector<std::pair<int,int> >& edge_length, std::vector<int>& leafe_depth);
	
	/**
	@brief indicates if a node during traversal has more outgoings
	@param start_index int describing the start index in indices_ vector
	@param stop_index int describing the end index in indices_ vector
	@param walked_in how many characters we have seen from root to actual position
	*/
	bool hasMoreOutgoings_(int start_index, int stop_index, int walked_in);

	const String& s_; ///< the string with which the suffix array is build

	double tol_; ///< mass tolerance for finding candidates
	
	std::vector<std::pair<int,int> > indices_; ///< vector of pairs of ints describing all relevant sufices

	std::vector<int> lcp_; ///< vector of ints with lcp values

	std::vector<int> skip_; ///< vector of ints with skip values

	//const int getIndex_ (const String & s);

	double masse_[256]; ///< mass table

	int number_of_modifications_; ///< number of allowed modifications

	std::vector<String> tags_; ///< all given tags

	bool use_tags_ ; ///< indicates whether tags are used or not

	int progress_;
};
}

#endif //OPENMS_DATASTRUCTURES_SUFFIXARRAYTRYPTICCOMPRESSED_H
