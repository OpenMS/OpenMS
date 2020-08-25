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
// $Maintainer: Xiao Liang $
// $Authors: Xiao Liang $
// --------------------------------------------------------------------------
//

#pragma once

#include <OpenMS/CHEMISTRY/DigestionEnzyme.h>

namespace OpenMS
{

  /**
      @ingroup Chemistry

      @brief Representation of a digestion enzyme for proteins (protease)
  */
  class OPENMS_DLLAPI DigestionEnzymeProtein :
    public DigestionEnzyme
  {
  public:

    /** @name Constructors
    */
    //@{

    /// Default constructor
    DigestionEnzymeProtein();

    /// Constructor from base class (adding defaults for the missing stuff)
    explicit DigestionEnzymeProtein(const DigestionEnzyme& d);

    /// Copy constructor
    DigestionEnzymeProtein(const DigestionEnzymeProtein&) = default;

    /// Move constructor
    DigestionEnzymeProtein(DigestionEnzymeProtein&&) = default;

    /// Detailed constructor
    explicit DigestionEnzymeProtein(const String& name,
                                    const String& cleavage_regex,
                                    const std::set<String>& synonyms = std::set<String>(),
                                    String regex_description = "",
                                    EmpiricalFormula n_term_gain = EmpiricalFormula("H"),
                                    EmpiricalFormula c_term_gain = EmpiricalFormula("OH"),
                                    String psi_id = "",
                                    String xtandem_id = "",
                                    Int comet_id = -1,
                                    String crux_id = "",
                                    Int msgf_id = -1,
                                    Int omssa_id = -1);

    /// Destructor
    ~DigestionEnzymeProtein() override;
    //@}

    /** @name Assignment
     */
    //@{
    /// Assignment operator
    DigestionEnzymeProtein& operator=(const DigestionEnzymeProtein&) = default;

    /// Move assignment operator
    DigestionEnzymeProtein& operator=(DigestionEnzymeProtein&&) & = default;
    //@}

    /** Accessors
    */
    //@{
    /// sets the N-terminal gain
    void setNTermGain(const EmpiricalFormula& value);

    /// returns N-terminal gain
    EmpiricalFormula getNTermGain() const;

    /// sets the C-terminal gain
    void setCTermGain(const EmpiricalFormula& value);

    /// returns C-terminal gain
    EmpiricalFormula getCTermGain() const;

    /// sets the PSI ID
    void setPSIID(const String& value);

    /// returns the PSI ID
    String getPSIID() const;

    /// sets the X! Tandem enzyme ID
    void setXTandemID(const String& value);

    /// returns the X! Tandem enzyme ID
    String getXTandemID() const;

    /// returns the Comet enzyme ID
    Int getCometID() const;

    /// sets the Comet enzyme ID
    void setCometID(Int value);

    /// returns the Crux enzyme ID
    String getCruxID() const;

    /// sets the Crux enzyme ID
    void setCruxID(const String& value);

    /// sets the MSGFPlus enzyme id
    void setMSGFID(Int value);

    /// returns the MSGFPlus enzyme id
    Int getMSGFID() const;

    /// sets the OMSSA enzyme ID
    void setOMSSAID(Int value);

    /// returns the OMSSA enzyme ID
    Int getOMSSAID() const;

    //@}

    /** @name Predicates
    */
    //@{
    /// equality operator
    bool operator==(const DigestionEnzymeProtein& enzyme) const;

    /// inequality operator
    bool operator!=(const DigestionEnzymeProtein& enzyme) const;

    // Note: comparison operator is not inherited. TODO rename and make virtual
    /// equality operator for regex
    bool operator==(const String& cleavage_regex) const;

    /// equality operator for regex
    bool operator!=(const String& cleavage_regex) const;

    /// order operator
    bool operator<(const DigestionEnzymeProtein& enzyme) const;
    //@}

    /**
       @brief Set the value of a member variable based on an entry from an input file

       Returns whether the key was recognized and the value set successfully.
    */
    bool setValueFromFile(const String& key, const String& value) override;

    /// ostream iterator to write the enzyme to a stream
    friend OPENMS_DLLAPI std::ostream& operator<<(std::ostream& os, const DigestionEnzymeProtein& enzyme);

  protected:
    EmpiricalFormula n_term_gain_;

    EmpiricalFormula c_term_gain_;

    String psi_id_;

    String xtandem_id_;

    Int comet_id_;

    String crux_id_;

    Int msgf_id_;

    Int omssa_id_;

  };

  OPENMS_DLLAPI std::ostream& operator<<(std::ostream& os, const DigestionEnzymeProtein& enzyme);

  typedef DigestionEnzymeProtein Protease;
}

