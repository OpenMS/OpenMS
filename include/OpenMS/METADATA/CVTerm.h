// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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

#ifndef OPENMS_METADATA_CVTERM_H
#define OPENMS_METADATA_CVTERM_H

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/DataValue.h>

namespace OpenMS
{
  /**
      @brief Representation of controlled vocabulary term

      This class simply stores a CV term, its value and unit if necessary.

      @ingroup Metadata
  */
  ///Represenation of a CV term used by CVMappings
  class OPENMS_DLLAPI CVTerm
  {
public:


    struct Unit
    {
      Unit()
      {
      }

      Unit(const String & p_accession, const String & p_name, const String & p_cv_ref) :
        accession(p_accession),
        name(p_name),
        cv_ref(p_cv_ref)
      {
      }

      Unit(const Unit & rhs) :
        accession(rhs.accession),
        name(rhs.name),
        cv_ref(rhs.cv_ref)
      {
      }

      virtual ~Unit()
      {
      }

      Unit & operator=(const Unit & rhs)
      {
        if (this != &rhs)
        {
          accession = rhs.accession;
          name = rhs.name;
          cv_ref = rhs.cv_ref;
        }
        return *this;
      }

      bool operator==(const Unit & rhs) const
      {
        return accession == rhs.accession &&
               name == rhs.name &&
               cv_ref == rhs.cv_ref;
      }

      bool operator!=(const Unit & rhs) const
      {
        return !(*this == rhs);
      }

      String accession;
      String name;
      String cv_ref;
    };


    /// Default constructor
    CVTerm();

    /// Detailed constructor
    CVTerm(const String & accession, const String & name, const String & cv_identifier_ref, const String & value, const Unit & unit);

    /// Copy constructor
    CVTerm(const CVTerm & rhs);

    /// Destructor
    virtual ~CVTerm();

    /// Assignment operator
    CVTerm & operator=(const CVTerm & rhs);

    /** @name Accessors
    */
    //@{
    /// sets the accession string of the term
    void setAccession(const String & accession);

    /// returns the accession string of the term
    const String & getAccession() const;

    /// sets the name of the term
    void setName(const String & name);

    /// returns the name of the term
    const String & getName() const;

    /// sets the cv identifier reference string, e.g. UO for unit obo
    void setCVIdentifierRef(const String & cv_identifier_ref);

    /// returns the cv identifier reference string
    const String & getCVIdentifierRef() const;

    /// set the value of the term
    void setValue(const DataValue & value);

    /// returns the value of the term
    const DataValue & getValue() const;

    /// sets the unit of the term
    void setUnit(const Unit & unit);

    /// returns the unit
    const Unit & getUnit() const;
    //@}

    /** @name Predicates
    */
    //@{
    /// equality operator
    bool operator==(const CVTerm & rhs) const;

    /// inequality operator
    bool operator!=(const CVTerm & rhs) const;

    /// checks whether the term has a value
    bool hasValue() const;

    /// checks whether the term has a unit
    bool hasUnit() const;
    //}

protected:

    String accession_;

    String name_;

    String cv_identifier_ref_;

    Unit unit_;

    DataValue value_;
  };

} // namespace OpenMS

#endif // OPENMS_DATASTRUCTURES_CVTERM_H
