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



#ifndef OPENMS_DATASTRUCTURES_SUFFIXARRAYSEQAN_H
#define OPENMS_DATASTRUCTURES_SUFFIXARRAYSEQAN_H

#include <OpenMS/CHEMISTRY/WeightWrapper.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/SuffixArray.h>
#include <OpenMS/DATASTRUCTURES/SeqanIncludeWrapper.h>

#include <functional>
#include <vector>

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/OpenMSConfig.h>

namespace OpenMS
{
    /**
  @brief comparator for two doubles with a tolerance value
  */
  struct FloatsWithTolLess :
	public std::binary_function<double, double, bool>
  {
    /**
    @brief constructor
    @param t const reference to the tolerance
    */
    explicit FloatsWithTolLess(const double & t) :
      tol_(t) {}
    /**
    @brief copy constructor
    */
    FloatsWithTolLess(const FloatsWithTolLess & rhs) :
      tol_(rhs.tol_) {}

    /**
    @brief implementation of the '<' operator for two doubles with the tolerance value
    @param f1 first double
    @param f2 second double
    @return true if first double '<' second double-tolerance
    */
    bool operator()(double f1, double f2) const
    {
      return f1 < (f2 - tol_);
    }

protected:
    double const & tol_;     ///< tolerance value
  };

  /**
  @brief comparator for two doubles with a tolerance value
  @todo Think about that this does and if it is really necessary (why double, double????) (Andreas, Clemens)
  */
  struct IntsInRangeLess :
	  public std::binary_function<double, double, bool>
  {
    /**
    @brief constructor
    @param t const reference to the tolerance
    */
    IntsInRangeLess(const int & s, const int & e) :
      start_(s), end_(e) {}
    /**
    @brief copy constructor
    */
    IntsInRangeLess(const IntsInRangeLess & source) :
      start_(source.start_), end_(source.end_) {}

    /**
    @brief implementation of the '<' operator for two doubles with the tolerance value
    @param f1 first double
    @param f2 second double
    @return true if first double '<' second double-tolerance
    */
    bool operator()(int f1, int f2) const
    {
      //cout<<"f1:"<<f1<<" f2:"<<f2<<" start:"<<start_<< " end:" << end_<<endl;
      return (f2 == end_) ? f1 <= f2 - start_ : f1 < f2;
    }

protected:
    int const & start_;     ///< start index
    int const & end_;     ///< end index
  };


  /**
  @brief Class that uses SEQAN library for a suffix array. It can be used to find peptide Candidates for a MS spectrum

  This class uses SEQAN suffix array. It can just be used for finding peptide Candidates for a given MS Spectrum within a certain mass tolerance. The suffix array can be saved to disc for reused so it has to be build just once.

  */

  class OPENMS_DLLAPI SuffixArraySeqan :
    public SuffixArray
    , public WeightWrapper
  {

    typedef seqan::TopDown<seqan::ParentLinks<> > TIterSpec;
    typedef seqan::Index<seqan::String<char>, seqan::IndexEsa<TIterSpec> > TIndex;
    typedef seqan::Iter<TIndex, seqan::VSTree<TIterSpec> > TIter;

    // TODO ??? was: typedef seqan::Index<seqan::String<char>, seqan::Index_ESA<seqan::TopDown<seqan::ParentLinks<seqan::Preorder> > > > TIndex;

public:

    /**
    @brief constructor
    @param st const string reference with the string for which the suffix array should be build
    @param filename const string reference with filename for opening or saving the suffix array
    @param weight_mode if not monoisotopic weight should be used, this parameters can be set to AVERAGE
    @throw FileNotFound is thrown if the given file is not found
    @throw InvalidValue if the given suffix array string is invalid
    */
    SuffixArraySeqan(const String & st, const String & filename, const WeightWrapper::WEIGHTMODE weight_mode = WeightWrapper::MONO);

    /**
    @brief copy constructor
    */
    SuffixArraySeqan(const SuffixArraySeqan & source);

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
    @return a vector of SignedSize pairs.

    For every mass within the spectrum all candidates described by as pairs of
    ints are returned. All masses are searched for the same time in just one
    suffix array traversal. In order to accelerate the traversal the skip and
    lcp table are used. The mass wont be calculated for each entry but it will
    be updated during traversal using a stack data structure
    */
    void findSpec(std::vector<std::vector<std::pair<std::pair<SignedSize, SignedSize>, double> > > & candidates, const std::vector<double> & spec);

    /**
    @brief saves the suffix array to disc
    @param filename const reference string describing the filename
    @return bool if operation was successful
    @throw UnableToCreateFile is thrown if the output files could not be created
    */
    bool save(const String & filename);

    /**
    @brief opens the suffix array

    @param filename const reference string describing the filename
    @return bool if operation was successful
    @throw FileNotFound is thrown if the given file could not be found
    */
    bool open(const String & filename);

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
    @param aa1 const char as first amino acid
    @param aa2 const char as second amino acid
    @return bool describing if it is a digesting site
    */
    bool isDigestingEnd(const char aa1, const char aa2) const;

    /**
    @brief setter for tags
    @param tags reference to vector of strings with tags
    @note sets use_tags = true
    */
    void setTags(const std::vector<OpenMS::String> & tags);

    /**
    @brief getter for tags
    @return const reference to vector of strings
    */
    const std::vector<OpenMS::String> & getTags();

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
    @return number of modifications
    */
    Size getNumberOfModifications();

    void printStatistic();

protected:

    /**
    @brief overwriting goNextSubTree_ from seqan index_esa_stree.h for mass update during suffix array traversal

    the suffix array is treated as a suffix tree. this function skips the subtree under the actual node and goes directly to the next subtree that has not been visited yet. During this traversal the mass will be updated using the stack with edge masses.

    @param it reference to the suffix array iterator
    @param m reference to actual mass
    @param allm reference to the stack with history of traversal
    @param mod_map input parameters which specifies the modification masses allowed in the candidates

    @see goNext
    */
    inline void goNextSubTree_(TIter & it, double & m, std::stack<double> & allm, std::stack<std::map<double, SignedSize> > & mod_map)
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
    inline void goNextSubTree_(TIter & it)
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
    inline void goNext_(TIter & it, double & m, std::stack<double> & allm, std::stack<std::map<double, SignedSize> > & mod_map)
    {
      // preorder dfs
      if (!goDown(it))
      {
        goNextSubTree_(it, m, allm, mod_map);
      }
    }

    inline void parseTree_(TIter & it, std::vector<std::pair<SignedSize, SignedSize> > & out_number, std::vector<std::pair<SignedSize, SignedSize> > & edge_length, std::vector<SignedSize> & leafe_depth)
    {
      SignedSize depth = 1;
      while (!atEnd(it))
      {
        SignedSize le = 0;
        bool isLeaf = false;
        if (length(parentEdgeLabel(it)) > 0)
        {
          if (countChildren(it) > 0)
          {
            edge_length.push_back(std::pair<SignedSize, SignedSize>(depth, length(parentEdgeLabel(it))));
          }
          else
          {
            //le <- length(representative(it));
            //isLeaf = true;
          }
        }
        if (countChildren(it) > 0)
        {
          out_number.push_back(std::pair<SignedSize, SignedSize>(depth, countChildren(it)));
        }
        else
        {
          leafe_depth.push_back(depth);
        }
        if (goDown(it))
        {
          depth++;
        }
        else if (!goRight(it))
        {
          while (!goRight(it))
          {
            goUp(it);
            if (isLeaf)
            {
              edge_length.push_back(std::pair<SignedSize, SignedSize>(depth, le - length(parentEdgeLabel(it))));
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

    TIndex index_;     ///< seqan suffix array

    TIter * it_;    ///< seqan suffix array iterator

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
    SignedSize findFirst_(const std::vector<double> & spec, double & m, SignedSize start, SignedSize  end);

    const String & s_;    ///< reference to strings for which the suffix array is build

    double masse_[255];     ///< amino acid masses

    SignedSize number_of_modifications_;     ///< number of allowed modifications

    std::vector<String> tags_;     ///< all tags

    bool use_tags_;     ///< if tags are used

    double tol_;     ///< tolerance
  };
}

#endif //OPENMS_DATASTRUCTURES_SUFFIXARRAYSEQAN_H
