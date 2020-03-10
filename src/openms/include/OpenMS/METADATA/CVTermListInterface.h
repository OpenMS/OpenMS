// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/METADATA/MetaInfoInterface.h>
#include <OpenMS/METADATA/CVTerm.h>
#include <OpenMS/DATASTRUCTURES/Map.h>

namespace OpenMS
{

  class CVTermList;

  /**
      @brief Interface to the controlled vocabulary term list

      This class is an interface to CVTermList, providing the same interface
      and it can be used instead (direct plug-in). It can be used to inherit
      from instead of CVTermList. The advantage is that this class only stores
      a single pointer to a CVTermList which may be more efficient for an empty
      list than storing the whole object.

      @ingroup Metadata
  */
  class OPENMS_DLLAPI CVTermListInterface :
    public MetaInfoInterface
  {

  public:

    /** @name Constructors and Assignment
    */
    //@{
    // Constructor
    CVTermListInterface();
    /// Copy constructor
    CVTermListInterface(const CVTermListInterface & rhs);
    /// Move constructor
    CVTermListInterface(CVTermListInterface&&) noexcept;
    // Destructor (non virtual)
    ~CVTermListInterface();

    /// Assignment operator
    CVTermListInterface & operator=(const CVTermListInterface & rhs);
    /// Move assignment operator
    CVTermListInterface& operator=(CVTermListInterface&&) noexcept;
    //@}

    /// equality operator
    bool operator==(const CVTermListInterface& rhs) const;

    /// inequality operator
    bool operator!=(const CVTermListInterface& rhs) const;

    void replaceCVTerms(Map<String, std::vector<CVTerm> > & cv_terms);

    /// sets the CV terms
    void setCVTerms(const std::vector<CVTerm>& terms);

    /// replaces the specified CV term
    void replaceCVTerm(const CVTerm& cv_term);

    /// replaces the specified CV terms using the given accession number
    void replaceCVTerms(const std::vector<CVTerm>& cv_terms, const String& accession);

    /// replaces all cv terms with a map (can be obtained via getCVTerms)
    void replaceCVTerms(const Map<String, std::vector<CVTerm> >& cv_term_map);

    /// merges the given map into the member map, no duplicate checking
    void consumeCVTerms(const Map<String, std::vector<CVTerm> >& cv_term_map);

    /// returns the accession string of the term
    const Map<String, std::vector<CVTerm> >& getCVTerms() const;

    /// adds a CV term
    void addCVTerm(const CVTerm& term);

    /// checks whether the term has a value
    bool hasCVTerm(const String& accession) const;

    bool empty() const;

  private:

    void createIfNotExists_(); 

    CVTermList* cvt_ptr_;
  };

} // namespace OpenMS

