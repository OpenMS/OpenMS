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


#ifndef OPENMS_DATASTRUCTURES_SUFFIXARRAYPEPTIDEFINDER_H
#define OPENMS_DATASTRUCTURES_SUFFIXARRAYPEPTIDEFINDER_H

#include <OpenMS/CHEMISTRY/WeightWrapper.h>
#include <OpenMS/DATASTRUCTURES/BigString.h>
#include <vector>

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/OpenMSConfig.h>

namespace OpenMS
{
  class SuffixArray;

  /**
      @brief wrapper for easy use of sufArray
  */
  class OPENMS_DLLAPI SuffixArrayPeptideFinder :
    public WeightWrapper
  {

public:

    /**
    @brief
    */
    typedef std::pair<String, String> FASTAEntry;

    /**
    @brief constructor
    @param filename FASTA File name
    @param method Name of the method used (trypticCompressed, seqan, trypticSeqan)
    @param weight_mode if not monoisotopic weight should be used, this parameters can be set to AVERAGE
    @throw FileNotFound is thrown if the filename is not found
    @throw ParseError is thrown if a error in parsing of the fasta file occurs
    @throw InvalidValue is thrown if an unknown method is supplied
    */
    SuffixArrayPeptideFinder(const String & filename, const String & method, const WeightWrapper::WEIGHTMODE weight_mode = WeightWrapper::MONO);

    /**
    @brief copy constructor
    */
    SuffixArrayPeptideFinder(const SuffixArrayPeptideFinder & source);

    /**
    @brief destructor
    */
    virtual ~SuffixArrayPeptideFinder();

    /**
    @brief finds all candidates for given spectrum in the suffix array
    @param spec vector holding the mass values to query
    @param candidates Output holding the candidates for input masses (one vector per mass)
                 FASTAEntry contains the FASTA header and the peptide sequence
                 The String contains the modification (if any) in the format specified by getModificationOutputMethod()
    @see sufArray.h
    */
    void getCandidates(std::vector<std::vector<std::pair<FASTAEntry, String> > > & candidates, const std::vector<double> & spec);

    /**
    @brief finds all candidate for given DTA file
    @param DTA_file DTA file location
    @param candidates Output parameters which holds the candidates suitable for the mass given in the dta file
                 FASTAEntry contains the FASTA header and the peptide sequence
                 The String contains the modification (if any) in the format specified by getModificationOutputMethod()
    @throw FileNotFound if DTA file does not exists
    @throw ParseError is thrown if the dta file could not be parsed
    @see sufArray.h
    */
    void getCandidates(std::vector<std::vector<std::pair<FASTAEntry, String> > > & candidates, const String & DTA_file);

    /**
    @brief allowed tolerance for mass match
    @param t Tolerance in u
    */
    void setTolerance(const double t);

    /**
    @brief allowed tolerance for mass match
    @return Tolerance in u
    */
    double getTolerance() const;

    /**
    @brief setter for number of modifications
    @param number_of_mods
    */
    void setNumberOfModifications(Size number_of_mods) const;

    /**
    @brief getter for number of modifications
    @return number of modifications
    */
    Size getNumberOfModifications() const;

    /**
    @brief setter for tags
    @param tags reference to vector of strings with tags
    @note sets use_tags = true
    */
    void setTags(const std::vector<String> & tags);

    /**
    @brief getter for tags
    @return const reference to vector of strings
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
    @brief set modification output method (valid are: "mass", "stringUnchecked", "stringChecked")
    @param s describing how modifications should be given back
    @throw InvalidValue is thrown if method s is not known
    */
    void setModificationOutputMethod(const String & s);

    /**
    @brief getter for modification output method
    @return String
    */
    String getModificationOutputMethod();

protected:

    String vToString_(std::vector<String> v);

    BigString big_string_;  ///< bigString object holding all peptides of fasta file

    SuffixArray * sa_;   ///< pointer to suffixarray

    String modification_output_method_; ///< output method for modifications

  };

}

#endif //OPENMS_EXAMPLES_SuffixArrayPeptideFinder_H
