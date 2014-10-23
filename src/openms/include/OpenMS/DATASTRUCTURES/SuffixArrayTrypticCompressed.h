// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2014.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Clemens Groepl,Andreas Bertsch$
// $Authors: Chris Bauer $
// --------------------------------------------------------------------------


#ifndef OPENMS_DATASTRUCTURES_SUFFIXARRAYTRYPTICCOMPRESSED_H
#define OPENMS_DATASTRUCTURES_SUFFIXARRAYTRYPTICCOMPRESSED_H

#include <OpenMS/CHEMISTRY/WeightWrapper.h>
#include <OpenMS/DATASTRUCTURES/SuffixArray.h>

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/OpenMSConfig.h>

namespace OpenMS
{
  class String;

  /**
    @brief Class that implements a suffix array for a String. It can be used to find peptide Candidates for a MS spectrum

    This class implements a suffix array. It can just be used for finding
    peptide Candidates for a given MS Spectrum within a certain mass tolerance.
    The suffix array can be saved to disc for reuse so it has to be built just
    once. The suffix array consists of a vector of pair of ints for every
    suffix, a vector of LCP values and a so called skip vector.

    Only the suffixes that are matching the function isDigestingEnd are
    created. Besides a suffix will not reach till the end of the string but
    till the next occurrence of the separator ($). So only the interesting
    suffixes will be saved. This will reduce the used space.
  */

  class OPENMS_DLLAPI SuffixArrayTrypticCompressed :
    public SuffixArray
    , public WeightWrapper
  {

public:

    /**
    @brief constructor taking the string and the filename for writing or reading
    @param st the string as const reference with which the suffix array will be build
    @param filename the filename for writing or reading the suffix array
    @param weight_mode if not monoisotopic weight should be used, this parameters can be set to AVERAGE
    @throw Exception::InvalidValue if string does not start with empty string ($)
    @throw FileNotFound is thrown if the given file was not found

    The constructor checks if a suffix array with given filename (without file
    extension) exists or not. In the first case it will simple be loaded and
    otherwise it will be build. Building the suffix array consists of several
    steps. At first all indices for a digesting enzyme (defined by using
    function isDigestingEnd) are created as an vector of SignedSize pairs.
    After creating all relevant indices they are sorted and the lcp and skip
    vectors are created.
    */
    SuffixArrayTrypticCompressed(const String & st, const String & filename, const WeightWrapper::WEIGHTMODE weight_mode = WeightWrapper::MONO);

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
    @return a vector of SignedSize pairs.
    @throw InvalidValue if the spectrum is not sorted ascendingly

    For every mass within the spectrum all candidates described by as pairs of
    ints are returned. All masses are searched for the same time in just one
    suffix array traversal. In order to accelerate the traversal the skip and
    lcp table are used. The mass wont be calculated for each entry but it will
    be updated during traversal using a stack data structure
    */
    void findSpec(std::vector<std::vector<std::pair<std::pair<SignedSize, SignedSize>, double> > > & candidates, const std::vector<double> & spec);

    /**
    @brief saves the suffix array to disc
    @param file_name const reference string describing the filename
    @return bool if operation was successful
    @throw Exception::UnableToCreateFile if file could not be created (e.g. if you have no rights)
    */
    bool save(const String & file_name);
    /**
    @brief opens the suffix array
    @param file_name const reference string describing the filename
    @return bool if operation was successful
    @throw FileNotFound
    */
    bool open(const String & file_name);

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
    @param aa1 const char as first amino acid
    @param aa2 const char as second amino acid
    @return bool describing if it is a digesting site
    */
    bool isDigestingEnd(const char aa1, const char aa2) const;

    /**
    @brief setter for tags
    @param tags const vector of strings with tags with length 3 each
    @throw InvalidValue if at least one tag does not have size of 3
    */
    void setTags(const std::vector<String> & tags);

    /**
    @brief getter for tags
    @return const vector of string with tags
    */
    const std::vector<String> & getTags();

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
    void setNumberOfModifications(Size number_of_mods);

    /**
    @brief getter for number of modifications
    @return unsigned SignedSize describing number of modifications
    */
    Size getNumberOfModifications();

    /**
    @brief output for statistic
    */
    void printStatistic();

protected:

    /**
    @brief constructor
    */
    SuffixArrayTrypticCompressed();

    /**
    @brief gets the index of the next separator for a given index
    @param p const SignedSize describing a position in the string
    @return SignedSize with the index of the next occurrence of the separator or -1 if there is no more separator
    */
    SignedSize getNextSep_(const SignedSize p) const;

    /**
    @brief gets the lcp for two strings described as pairs of ints
    @param last_point const pair of ints describing a substring
    @param current_point const pair of ints describing a substring
    @return SignedSize with the length of the lowest common prefix
    */
    SignedSize getLCP_(const std::pair<SignedSize, SignedSize> & last_point, const std::pair<SignedSize, SignedSize> & current_point);

    /**
    @brief binary search for finding the index of the first element of the spectrum that matches the desired mass within the tolerance.
    @param spec const reference to spectrum
    @param m mass
    @return SignedSize with the index of the first occurrence
    @note requires that there is at least one occurrence
    */
    SignedSize findFirst_(const std::vector<double> & spec, double & m);

    /**
    @brief binary search for finding the index of the first element of the
    spectrum that matches the desired mass within the tolerance. It searches
    recursively.
    @param spec const reference to spectrum
    @param m mass
    @param start start index
    @param end end index
    @return SignedSize with the index of the first occurrence
    @note requires that there is at least one occurrence
    */
    SignedSize findFirst_(const std::vector<double> & spec, double & m, SignedSize start, SignedSize end);

    /**
    @brief treats the suffix array as a tree and parses the tree using postorder traversal. This is realised by a recursive algorithm.
    @param start_index SignedSize describing the start index in indices_ vector
    @param stop_index SignedSize describing the end index in indices_ vector
    @param depth at with depth the traversal is at the actual position
    @param walked_in how many characters we have seen from root to actual position
    @param edge_len how many characters we have seen from last node to actual position
    @param out_number reference to vector of pairs of ints. For every node it will be filled with how many outgoing edge a node has in dependence of its depth
    @param edge_length will be filled with the edge_length in dependence of its depth
    @param leafe_depth will be filled with the depth of every leaf
    @note initialize: walked_in=0, depth=1, edge_len=1
    */
    void parseTree_(SignedSize start_index, SignedSize stop_index, SignedSize depth, SignedSize walked_in, SignedSize edge_len, std::vector<std::pair<SignedSize, SignedSize> > & out_number, std::vector<std::pair<SignedSize, SignedSize> > & edge_length, std::vector<SignedSize> & leafe_depth);

    /**
    @brief indicates if a node during traversal has more outgoings
    @param start_index SignedSize describing the start index in indices_ vector
    @param stop_index SignedSize describing the end index in indices_ vector
    @param walked_in how many characters we have seen from root to actual position
    */
    bool hasMoreOutgoings_(SignedSize start_index, SignedSize stop_index, SignedSize walked_in);

    const String & s_; ///< the string with which the suffix array is build

    double tol_; ///< mass tolerance for finding candidates

    std::vector<std::pair<SignedSize, SignedSize> > indices_; ///< vector of pairs of ints describing all relevant suffixes

    std::vector<SignedSize> lcp_; ///< vector of ints with lcp values

    std::vector<SignedSize> skip_; ///< vector of ints with skip values

    //const SignedSize getIndex_ (const String & s);

    double masse_[256]; ///< mass table

    Size number_of_modifications_; ///< number of allowed modifications

    std::vector<String> tags_; ///< all given tags

    bool use_tags_;  ///< indicates whether tags are used or not

    SignedSize progress_;
  };
}

#endif //OPENMS_DATASTRUCTURES_SUFFIXARRAYTRYPTICCOMPRESSED_H
