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
// $Maintainer: Andreas Bertsch $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#ifndef OPENMS_DATASTRUCTURES_CVMAPPINGTERM_H
#define OPENMS_DATASTRUCTURES_CVMAPPINGTERM_H

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/OpenMSConfig.h>

namespace OpenMS
{
  /**
      @brief Representation of controlled vocabulary term

      This class simply stores CV terms read from e.g. an OBO-file

      @ingroup Datastructures
  */
  ///Representation of a CV term used by CVMappings
  class OPENMS_DLLAPI CVMappingTerm
  {
public:

    /// Defaults constructor
    CVMappingTerm();

    /// Copy constructor
    CVMappingTerm(const CVMappingTerm & rhs);

    /// Destructor
    virtual ~CVMappingTerm();

    /// Assignment operator
    CVMappingTerm & operator=(const CVMappingTerm & rhs);

    /** @name Accessors
    */
    //@{
    /// sets the accession string of the term
    void setAccession(const String & accession);

    /// returns the accession string of the term
    const String & getAccession() const;

    /// sets whether the term name should be used, instead of the accession
    void setUseTermName(bool use_term_name);

    /// returns whether the term name should be used, instead of the accession
    bool getUseTermName() const;

    /// sets whether the term itself can be used (or only its children)
    void setUseTerm(bool use_term);

    /// returns true if the term can be used, false if only children are allowed
    bool getUseTerm() const;

    /// sets the name of the term
    void setTermName(const String & term_name);

    /// returns the name of the term
    const String & getTermName() const;

    /// sets whether this term can be repeated
    void setIsRepeatable(bool is_repeatable);

    /// returns true if this term can be repeated, false otherwise
    bool getIsRepeatable() const;

    /// sets whether children of this term are allowed
    void setAllowChildren(bool allow_children);

    /// returns true if the children of this term are allowed to be used
    bool getAllowChildren() const;

    /// sets the cv identifier reference string, e.g. UO for unit obo
    void setCVIdentifierRef(const String & cv_identifier_ref);

    /// returns the cv identifier reference string
    const String & getCVIdentifierRef() const;
    //@}

    /** @name Predicates
    */
    //@{
    /// equality operator
    bool operator==(const CVMappingTerm & rhs) const;

    /// inequality operator
    bool operator!=(const CVMappingTerm & rhs) const;
    //}

protected:

    String accession_;

    bool use_term_name_;

    bool use_term_;

    String term_name_;

    bool is_repeatable_;

    bool allow_children_;

    String cv_identifier_ref_;
  };

} // namespace OpenMS

#endif // OPENMS_DATASTRUCTURES_CVMAPPINGTERM_H
