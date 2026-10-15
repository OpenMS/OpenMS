// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
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
    explicit DigestionEnzymeProtein(const std::string& name,
                                    const std::string& cleavage_regex,
                                    const std::set<std::string>& synonyms = std::set<std::string>(),
                                    std::string regex_description = "",
                                    EmpiricalFormula n_term_gain = EmpiricalFormula("H"),
                                    EmpiricalFormula c_term_gain = EmpiricalFormula("OH"),
                                    std::string psi_id = "",
                                    std::string xtandem_id = "",
                                    Int comet_id = -1,
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
    void setPSIID(const std::string& value);

    /// returns the PSI ID
    std::string getPSIID() const;

    /// sets the X! Tandem enzyme ID
    void setXTandemID(const std::string& value);

    /// returns the X! Tandem enzyme ID
    std::string getXTandemID() const;

    /// sets the Comet enzyme ID
    void setCometID(Int value);

    /// returns the Comet enzyme ID
    Int getCometID() const;

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
    bool operator==(const std::string& cleavage_regex) const;

    /// equality operator for regex
    bool operator!=(const std::string& cleavage_regex) const;

    /// order operator
    bool operator<(const DigestionEnzymeProtein& enzyme) const;
    //@}

    /**
       @brief Set the value of a member variable based on an entry from an input file

       Returns whether the key was recognized and the value set successfully.
    */
    bool setValueFromFile(const std::string& key, const std::string& value) override;

    /// ostream iterator to write the enzyme to a stream
    friend OPENMS_DLLAPI std::ostream& operator<<(std::ostream& os, const DigestionEnzymeProtein& enzyme);

  protected:
    EmpiricalFormula n_term_gain_;

    EmpiricalFormula c_term_gain_;

    std::string psi_id_;

    std::string xtandem_id_;

    Int comet_id_;

    Int msgf_id_;

    Int omssa_id_;

  };

  OPENMS_DLLAPI std::ostream& operator<<(std::ostream& os, const DigestionEnzymeProtein& enzyme);

  typedef DigestionEnzymeProtein Protease;
}

