// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2021.
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
// $Maintainer: Hannes Rost $
// $Authors: Hannes Rost $
// --------------------------------------------------------------------------
//

#pragma once

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CHEMISTRY/Element.h>

#include <string>

namespace OpenMS
{
  /** @ingroup Chemistry

      @brief Representation of an isotope

      This class represents a single isotope of an Element. Stable isotopes are
      expected to contain abundances while unstable isotopes generally do not
      have natural abundances but contain information about half life and decay
      mode.
  */
  class OPENMS_DLLAPI Isotope :
    public Element
  {
    public:

    enum DecayMode
    {
      NONE = 0,     ///< No decay (stable isotope)
      UNKNOWN,      ///< Unknown / Unspecified decay mode
      ALPHA,        ///< Alpha decay
      BETA_PLUS,    ///< Beta plus decay
      BETA_MINUS,   ///< Beta minus decay
      PROTON,       ///< Proton emission
      SIZE_OF_DECAYMODE
    };

    using Element::Element;

    /// detailed constructor
    Isotope(const std::string & name,
            const std::string & symbol,
            unsigned int atomic_number,
            unsigned int neutrons,
            double mono_weight,
            double abundance,
            double half_life,
            Isotope::DecayMode dm);

    Isotope & operator=(const Isotope & element) = default;
    Isotope(const Isotope & element) = default;

    /// destructor
    virtual ~Isotope();

    /// Get corresponding element
    const Element* getElement() const;

    /// set isotope half life in seconds
    void setHalfLife(double hl);
    /// get isotope half life in seconds
    double getHalfLife() const;
    /// set isotope natural abundance
    void setAbundance(double ab);
    /// get isotope natural abundance
    double getAbundance() const;
    /// set number of neutrons
    void setNeutrons(int ne);
    /// get number of neutrons
    int getNeutrons() const;
    /// set primary decay mode (for unstable isotopes)
    void setDecayMode(DecayMode dm);
    /// get primary decay mode (for unstable isotopes)
    DecayMode getDecayMode() const;

    /// Whether this is an Isotope or an Element (for casting)
    virtual bool isIsotope() const override {return true;}

    /// Whether this is a stable isotope
    bool isStable() const {return half_life_ < 0;}

    /// writes the isotope to an output stream
    friend OPENMS_DLLAPI std::ostream & operator<<(std::ostream & os, const Isotope & isotope);

  protected:

    int neutrons_ = -1;
    double abundance_ = -1; // has to consistent with getIsotopeDistribution()
    double half_life_ = -1; ///< half life in seconds
    DecayMode decay_mode_ = DecayMode::NONE;
  };

  OPENMS_DLLAPI std::ostream & operator<<(std::ostream &, const Isotope &);

} // namespace OpenMS

