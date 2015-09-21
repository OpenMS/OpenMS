// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2015.
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
// $Maintainer: Xiao Liang $
// $Authors: Xiao Liang $
// --------------------------------------------------------------------------
//

#ifndef OPENMS_CHEMISTRY_ENZYME_H
#define OPENMS_CHEMISTRY_ENZYME_H

#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <iosfwd>
#include <set>
#include <vector>

namespace OpenMS
{

  /**
      @ingroup Chemistry

      @brief Representation of an enzyme

      This class represents enzymes.
  */
  class OPENMS_DLLAPI Enzyme
  {
public:

    /** @name Constructors
    */
    //@{
    /// copy constructor
    Enzyme(const Enzyme & enzyme);

    /// detailed constructor
    explicit Enzyme(const String & name,
                    const String & cleavage_regex,
                    const std::set<String> & synonyms = std::set<String>(),
                    String regex_description = "",
                    EmpiricalFormula n_term_gain = EmpiricalFormula("H"),
                    EmpiricalFormula c_term_gain = EmpiricalFormula("OH"),
                    String psi_id = "",
                    String xtandem_id = "",
                    UInt omssa_id = 0);

    /// destructor
    virtual ~Enzyme();
    //@}

    /** @name Assignment
     */
    //@{
    /// assignment operator
    Enzyme & operator=(const Enzyme & enzyme);
    //@}

    /** Accessors
    */
    //@{
    /// sets the name of the enzyme
    void setName(const String & name);

    /// returns the name of the enzyme
    const String & getName() const;

    /// sets the synonyms
    void setSynonyms(const std::set<String> & synonyms);

    /// adds a synonym
    void addSynonym(const String & synonym);

    /// returns the synonyms
    const std::set<String> & getSynonyms() const;

    /// sets the name as regex
    void setRegEx(const String & cleavage_regex);

    /// returns the name as regex
    const String & getRegEx() const;

    /// sets the regex description
    void setRegExDescription(String value);

    /// returns the regex description
    String getRegExDescription() const;

    /// sets the N-terminal gain
    void setNTermGain(EmpiricalFormula value);

    /// returns N-terminal gain
    EmpiricalFormula getNTermGain() const;

    /// sets the C-terminal gain
    void setCTermGain(EmpiricalFormula value);

    /// returns C-terminal gain
    EmpiricalFormula getCTermGain() const;

    /// sets the PSI id
    void setPSIid(String value);

    /// returns the PSI id
    String getPSIid() const;

    /// sets the XTANDEM enzyme id
    void setXTANDEMid(String value);

    /// returns the XTANDEM enzyme id
    String getXTANDEMid() const;

    /// sets the OMSSA enzyme id
    void setOMSSAid(UInt value);

    /// returns the OMSSA enzyme id
    UInt getOMSSAid() const;
    //@}

    /** @name Predicates
    */
    //@{
    /// equality operator
    bool operator==(const Enzyme & enzyme) const;

    /// inequality operator
    bool operator!=(const Enzyme & enzyme) const;

    /// equality operator for regex
    bool operator==(String cleavage_regex) const;

    /// equality operator for regex
    bool operator!=(String cleavage_regex) const;

    /// order operator
    bool operator<(const Enzyme & enzyme) const;
    //@}

    /// ostream iterator to write the enzyme to a stream
    friend OPENMS_DLLAPI std::ostream & operator<<(std::ostream & os, const Enzyme & enzyme);

protected:
    /// default constructor
    Enzyme();

    // basic
    String name_;

    String cleavage_regex_;

    std::set<String> synonyms_;

    String regex_description_;

    EmpiricalFormula n_term_gain_;

    EmpiricalFormula c_term_gain_;

    String psi_id_;

    String xtandem_id_;

    UInt omssa_id_;
  };

  OPENMS_DLLAPI std::ostream & operator<<(std::ostream & os, const Enzyme & enzyme);

}

#endif
