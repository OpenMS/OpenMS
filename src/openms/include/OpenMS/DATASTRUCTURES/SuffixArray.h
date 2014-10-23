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


#ifndef OPENMS_DATASTRUCTURES_SUFFIXARRAY_H
#define OPENMS_DATASTRUCTURES_SUFFIXARRAY_H

#include <vector>

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/OpenMSConfig.h>


namespace OpenMS
{
  class String;

  /**
    @brief abstract class for suffix array
  */
  class OPENMS_DLLAPI SuffixArray
  {

public:

    /**
    @brief constructor taking the string and the filename for writing or reading
    @param st the string as const reference with which the suffix array will be build
    @param filename the filename for writing or reading the suffix array
    @throw Exception::InvalidValue if string does not start with empty string ($)
    @throw Exception::FileNotFound is thrown if the given filename is not found
    */
    SuffixArray(const String & st, const String & filename);

    /**
    @brief copy constructor
    */
    SuffixArray(const SuffixArray & sa);

    /**
    @brief destructor
    */
    virtual ~SuffixArray() = 0;

    /**
    @brief transforms suffix array to a printable String
    */
    virtual String toString() = 0;

    /**
    @brief the function that will find all peptide candidates for a given spectrum
    @param spec const reference of double vector describing the spectrum
    @param candidates the candidates which are returned for the masses given in spec
    @return a vector of SignedSize pairs.
    @throw InvalidValue if the spectrum is not sorted ascendingly

    */
    virtual void findSpec(std::vector<std::vector<std::pair<std::pair<SignedSize, SignedSize>, double> > > & candidates, const std::vector<double> & spec) = 0;

    /**
    @brief saves the suffix array to disc
    @param filename const reference string describing the filename
    @return bool if operation was successful
    @throw UnableToCreateFile if file could not be created (e.g. if you have no rights)
    */
    virtual bool save(const String & filename) = 0;
    /**
    @brief opens the suffix array
    @param filename const reference string describing the filename
    @return bool if operation was successful
    @throw FileNotFound
    */
    virtual bool open(const String & filename) = 0;

    /**
    @brief setter for tolerance
    @param t double with tolerance
    @throw InvalidValue if tolerance is negative
    */
    virtual void setTolerance(double t) = 0;

    /**
    @brief getter for tolerance
    @return double with tolerance
    */
    virtual double getTolerance() const  = 0;

    /**
    @brief returns if an enzyme will cut after first character
    @param aa1 const char as first aminoacid
    @param aa2 const char as second aminoacid
    @return bool describing if it is a digesting site
    */
    virtual bool isDigestingEnd(const char aa1, const char aa2) const = 0;

    /**
    @brief setter for tags
    @param tags const vector of strings with tags with length 3 each
    @throw Exception::InvalidValue if at least one tag does not have size of 3
    */
    virtual void setTags(const std::vector<String> & tags) = 0;

    /**
    @brief getter for tags
    @return const vector of string with tags
    */
    virtual const std::vector<String> & getTags() = 0;

    /**
    @brief setter for use_tags
    @param use_tags indicating whether tags should be used or not
    */
    virtual void setUseTags(bool use_tags) = 0;

    /**
    @brief getter for use_tags
    @return bool indicating whether tags are used or not
    */
    virtual bool getUseTags() = 0;

    /**
    @brief setter for number of modifications
    @param number_of_mods
    */
    virtual void setNumberOfModifications(Size number_of_mods) = 0;

    /**
    @brief getter for number of modifications
    @return Size describing number of modifications
    */
    virtual Size getNumberOfModifications() = 0;

    /**
    @brief output for statistic
    */
    virtual void printStatistic() = 0;

    /**
    @brief constructor
    */
    SuffixArray();


  };
}

#endif //OPENMS_DATASTRUCTURES_SUFARRAY_H
