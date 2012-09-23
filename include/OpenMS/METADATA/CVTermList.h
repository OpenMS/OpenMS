// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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
// $Maintainer: Andreas Bertsch $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#ifndef OPENMS_METADATA_CVTERMLIST_H
#define OPENMS_METADATA_CVTERMLIST_H

#include <OpenMS/METADATA/CVTerm.h>
#include <OpenMS/METADATA/MetaInfoInterface.h>
#include <OpenMS/DATASTRUCTURES/CVMappingRule.h>
#include <OpenMS/FORMAT/ControlledVocabulary.h>

namespace OpenMS
{
  /**
      @brief Representation of controlled vocabulary term list

      This class should be used to inherit from, to allow to add
      an arbitrary number of CV terms to the inherited class

      @ingroup Metadata
  */
  ///Represenation of a CV term used by CVMappings
  class OPENMS_DLLAPI CVTermList :
    public MetaInfoInterface
  {
public:

    /// Defaults constructor
    CVTermList();

    /// Copy constructor
    CVTermList(const CVTermList & rhs);

    /// Destructor
    virtual ~CVTermList();

    /// Assignment operator
    CVTermList & operator=(const CVTermList & rhs);

    /** @name Accessors
    */
    //@{
    /// sets the CV terms
    void setCVTerms(const std::vector<CVTerm> & terms);

    /// replaces the specified CV term
    void replaceCVTerm(const CVTerm & cv_term);

    /// replaces the specified CV terms using the given accession number
    void replaceCVTerms(const std::vector<CVTerm> & cv_terms, const String & accession);

    /// replaces all cv terms with a map (can be obtained via getCVTerms)
    void replaceCVTerms(const Map<String, std::vector<CVTerm> > & cv_term_map);

    /// returns the accession string of the term
    const Map<String, std::vector<CVTerm> > & getCVTerms() const;

    /// adds a CV term
    void addCVTerm(const CVTerm & term);

    /// checks whether the spellings of the CV terms stored are correct
    //bool checkCVTerms(const ControlledVocabulary& cv) const;

    /// corrects the CVTerm names, according to the loaded CV
    //void correctCVTermNames();
    //@}

    /** @name Predicates
    */
    //@{
    /// equality operator
    bool operator==(const CVTermList & cv_term_list) const;

    /// inequality operator
    bool operator!=(const CVTermList & cv_term_list) const;

    /// checks whether the term has a value
    bool hasCVTerm(const String & accession) const;

    /// checks whether the stored terms fullfil a given CVMappingRule
    bool checkCVTerms(const CVMappingRule & rule, const ControlledVocabulary & cv) const;

    /// return true if no terms are available
    bool empty() const;
    //}

protected:

    Map<String, std::vector<CVTerm> > cv_terms_;

  };

} // namespace OpenMS

#endif // OPENMS_DATASTRUCTURES_CVTERMLIST_H
