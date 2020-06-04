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
// $Maintainer: Timo Sachsenberg$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>
#include <string>

/**
    @brief Main %OpenMS namespace.

    In this namespace all the main %OpenMS classes are located.
*/
namespace OpenMS
{
  /**
      @brief Mathematical and physical constants namespace.

      This namespace contains definitions for some basic mathematical and physical constants.
      All constants are double precision.
      <BR>
      There are basically two ways of accessing these constants:
      <UL>
          <LI> specify all namespaces:
          <BR>
          <tt>double my_pi = OpenMS::Constants::PI</tt>
          <BR>
          <LI>shortcut via the <tt>using directive</tt>:
          <BR>
          <tt>using namespace OpenMS::Constants;
          <BR>
          double my_pi = PI;</tt>
      </UL>

      @ingroup Concept
  */
  namespace Constants
  {
    /**	@name	Mathematical constants.
    */
    //@{

    /// PI
    extern OPENMS_DLLAPI const double  PI;

    /// Euler's number - base of the natural logarithm
    extern OPENMS_DLLAPI const double  E;

    /**	Internal threshold for equality comparisons.
            Default value is 1e-6.
    */
    extern OPENMS_DLLAPI double EPSILON;
    //@}

    /**	@name Chemical/physical constants.
    */
    //@{

    /**	Elementary charge.
          In units of C (\f$1.60217738 \cdot 10^{-19} C\f$).
    */
    extern OPENMS_DLLAPI const double   ELEMENTARY_CHARGE;           // C

    /// Elementary charge (alias)
    extern OPENMS_DLLAPI const double   e0;

    /** Electron mass.
            In units of kg (\f$9.1093897 \cdot 10^{-31}\f$ kg).
    */
    extern OPENMS_DLLAPI const double   ELECTRON_MASS;                   // kg

    /** Electron mass
            In units (\f$1,822.88850204(77)^{-1}\f$u).
    */
    extern OPENMS_DLLAPI const double ELECTRON_MASS_U;             // u

    /** Proton mass.
            In units of kg (\f$1.6726230 \cdot 10^{-27}\f$ kg).
    */
    extern OPENMS_DLLAPI const double   PROTON_MASS;                     // kg

    /** Proton mass.
            In units (\f$1.00727646677(10)\f$u)
    */
    extern OPENMS_DLLAPI const double PROTON_MASS_U;               // u

    /** C13C12 mass difference.
        In units (\f$1.0033548\f$u)
    */
    extern OPENMS_DLLAPI const double C13C12_MASSDIFF_U;     // u

    /** Neutron mass.
            In units of kg (\f$1.6749286 \cdot 10^{-27}\f$ kg).
    */
    extern OPENMS_DLLAPI const double   NEUTRON_MASS;                    // kg

    /** Neutron mass.
            In units (\f$1.0086649156(6)\f$u)
    */
    extern OPENMS_DLLAPI const double NEUTRON_MASS_U;              // u

    /** Avogadro constant.
            In units of \f$mol^{-1}\f$ (\f$6.0221367 \cdot 10^{23} mol^{-1}\f$).
    */
    extern OPENMS_DLLAPI const double   AVOGADRO;

    /** Avogadro constant (alias)
    */
    extern OPENMS_DLLAPI const double   NA;

    /** Avogadro constant (alias)
    */
    extern OPENMS_DLLAPI const double   MOL;

    /** Boltzmann constant.
            In units of J/K (\f$1.380657 \cdot 10^{-23}\f$ J/K).
    */
    extern OPENMS_DLLAPI const double   BOLTZMANN;

    /** Boltzmann constant (alias)
    */
    extern OPENMS_DLLAPI const double   k;

    /** Planck constant.
            In units of Js (\f$6.6260754 \cdot 10^{-34}\f$ Js).
    */
    extern OPENMS_DLLAPI const double   PLANCK;

    /** Planck constant (alias)
    */
    extern OPENMS_DLLAPI const double   h;

    /** Gas constant (= NA * k)
    */
    extern OPENMS_DLLAPI const double   GAS_CONSTANT;

    /** Gas constant (alias)
    */
    extern OPENMS_DLLAPI const double R;

    /** Faraday constant (= NA * e0)
    */
    extern OPENMS_DLLAPI const double   FARADAY;

    /** Faraday constant (alias)
    */
    extern OPENMS_DLLAPI const double   F;

    /** Bohr radius.
            In units m (\f$5.29177249 \cdot 10^{-11}\f$ m).
    */
    extern OPENMS_DLLAPI const double   BOHR_RADIUS;

    /** Bohr radius (alias)
    */
    extern OPENMS_DLLAPI const double   a0;

    //  the following values from:
    //  P.W.Atkins: Physical Chemistry, 5th ed., Oxford University Press, 1995

    /** Vacuum permittivity.
            In units of \f$C^2J^{-1}m^{-1}\f$ (\f$8.85419 \cdot 10^{-12} C^2J^{-1}m^{-1}\f$).
    */
    extern OPENMS_DLLAPI const double   VACUUM_PERMITTIVITY;

    /** Vacuum permeability.
            In units of \f$Js^2C^{-2}m^{-1}\f$ (\f$4\pi \cdot 10^{-7} Js^2C^{-2}m^{-1}\f$).
    */
    extern OPENMS_DLLAPI const double   VACUUM_PERMEABILITY;

    /** Speed of light.
            In units of m/s (\f$2.99792458 \cdot 10^8 ms^{-1}\f$).
    */
    extern OPENMS_DLLAPI const double   SPEED_OF_LIGHT;

    /** Speed of Light (alias)
    */
    extern OPENMS_DLLAPI const double   c;

    /** Gravitational constant.
            In units of \f$Nm^2kg^{-2}\f$ (\f$6.67259 \cdot 10^{-11} Nm^2kg^{-2}\f$).
    */
    extern OPENMS_DLLAPI const double   GRAVITATIONAL_CONSTANT;

    /** Fine structure constant.
            Without unit (\f$7.29735 \cdot 10^{-3}\f$).
    */
    extern OPENMS_DLLAPI const double   FINE_STRUCTURE_CONSTANT;
    //@}

    /**	@name	Conversion factors
    */
    //@{

    /** Degree per rad.
            57.2957795130823209
    */
    extern OPENMS_DLLAPI const double   DEG_PER_RAD;

    /** Rad per degree.
            0.0174532925199432957
    */
    extern OPENMS_DLLAPI const double   RAD_PER_DEG;

    /** mm per inch.
            25.4
    */
    extern OPENMS_DLLAPI const double   MM_PER_INCH;

    /** m per foot.
            3.048
    */
    extern OPENMS_DLLAPI const double   M_PER_FOOT;

    /** Joules per calorie.
            4.184
    */
    extern OPENMS_DLLAPI const double   JOULE_PER_CAL;

    /** Calories per Joule.
            1/JOULE_PER_CAL
    */
    extern OPENMS_DLLAPI const double   CAL_PER_JOULE;

    namespace UserParam
    {
      /** User parameter name for the M/Z of other chromatograms which have been merged into this one
              String
       */
      extern OPENMS_DLLAPI const std::string   MERGED_CHROMATOGRAM_MZS;

      /** User parameter name for precursor mz error in ppm
              String
      */
      extern OPENMS_DLLAPI const std::string   PRECURSOR_ERROR_PPM_USERPARAM;

      /** User parameter name for fragment mz error in ppm
              String
      */
      extern OPENMS_DLLAPI const std::string   FRAGMENT_ERROR_MEDIAN_PPM_USERPARAM;

      /** User parameter name for fragment annotations
              String
      */
      extern OPENMS_DLLAPI const std::string   FRAGMENT_ANNOTATION_USERPARAM;

      /** User parameter name for the spectrum reference in PeptideIdentification (is is not yet treated as a class attribute)
              String
      */
      extern OPENMS_DLLAPI const std::string   SPECTRUM_REFERENCE;

      /** User parameter name for target/decoy annotation of a PeptideHit, e.g. as annotated by PeptideIndexer. One of: target, decoy, target+decoy
              String
      */
      extern OPENMS_DLLAPI const std::string   TARGET_DECOY;

      /** User parameter name for a delta score: a score ratio between a rank x hit and the rank x+1 hit
              String
      */
      extern OPENMS_DLLAPI const std::string   DELTA_SCORE;

      /** User parameter name to indicate a monoisotopic peak misassignment. Used for precursor correction. (usually an integer x with the correction being -x times C13C12_MASSDIFF_U)
              String
      */
      extern OPENMS_DLLAPI const std::string   ISOTOPE_ERROR;

      // Cross-Linking Mass Spectrometry user parameters
      /** Name of OpenPepXL main score (PSI CV term)
              String
      */
      extern OPENMS_DLLAPI const std::string   OPENPEPXL_SCORE;

      /** User parameter name for the sequence of the second peptide in a cross-link
              String
      */
      extern OPENMS_DLLAPI const std::string   OPENPEPXL_BETA_SEQUENCE;

      /** User parameter name for the protein accessions of the second peptide in a cross-link
              String
      */
      extern OPENMS_DLLAPI const std::string   OPENPEPXL_BETA_ACCESSIONS;

      /** User parameter name for the 1st position of cross-link (alpha peptide position in a real cross-link, 1st of two positions in a loop-link, modified position in a mono-link)
              String
      */
      extern OPENMS_DLLAPI const std::string   OPENPEPXL_XL_POS1;

      /** User parameter name for the 2nd position of cross-link (beta peptide position in a real cross-link, 2nd of two positions in a loop-link, "-" in a mono-link)
              String
      */
      extern OPENMS_DLLAPI const std::string   OPENPEPXL_XL_POS2;

      /** User parameter name for the 1st cross-link position on the protein
              String
      */
      extern OPENMS_DLLAPI const std::string   OPENPEPXL_XL_POS1_PROT;

      /** User parameter name for the 2nd cross-link position on the protein
              String
      */
      extern OPENMS_DLLAPI const std::string   OPENPEPXL_XL_POS2_PROT;

      /** User parameter name for the cross-link type, one of: cross-link, loop-link, mono-link
              String
      */
      extern OPENMS_DLLAPI const std::string   OPENPEPXL_XL_TYPE;

      /** User parameter name for the cross-link rank (ranks of PeptideHits across different PeptideIdentifications)
              String
      */
      extern OPENMS_DLLAPI const std::string   OPENPEPXL_XL_RANK;

      /** User parameter name for the name of a cross-link
              String
      */
      extern OPENMS_DLLAPI const std::string   OPENPEPXL_XL_MOD;

      /** User parameter name for the mass of a cross-link
              String
      */
      extern OPENMS_DLLAPI const std::string   OPENPEPXL_XL_MASS;

      /** User parameter name for the terminal specificity of a cross-link on the alpha peptide (to distinguish a link to the first or last residue side chain from a terminal link)
              String
      */
      extern OPENMS_DLLAPI const std::string   OPENPEPXL_XL_TERM_SPEC_ALPHA;

      /** User parameter name for the terminal specificity of a cross-link on the beta peptide (to distinguish a link to the first or last residue side chain from a terminal link)
              String
      */
      extern OPENMS_DLLAPI const std::string   OPENPEPXL_XL_TERM_SPEC_BETA;

      /** User parameter name for the RT of the heavy spectrum precursor in a labeled cross-linking experiment
              String
      */
      extern OPENMS_DLLAPI const std::string   OPENPEPXL_HEAVY_SPEC_RT;

      /** User parameter name for the m/z of the heavy spectrum precursor in a labeled cross-linking experiment
              String
      */
      extern OPENMS_DLLAPI const std::string   OPENPEPXL_HEAVY_SPEC_MZ;

      /** User parameter name for the spectrum reference of the heavy spectrum in a labeled cross-linking experiment
              String
      */
      extern OPENMS_DLLAPI const std::string   OPENPEPXL_HEAVY_SPEC_REF;

      /** User parameter name for target/decoy annotation of alpha peptides
              String
      */
      extern OPENMS_DLLAPI const std::string   OPENPEPXL_TARGET_DECOY_ALPHA;

      /** User parameter name for target/decoy annotation of beta peptides
              String
      */
      extern OPENMS_DLLAPI const std::string   OPENPEPXL_TARGET_DECOY_BETA;

      /** User parameter name for PeptideEvidence info for the beta/acceptor peptide: pre
              String
      */
      extern OPENMS_DLLAPI const std::string   OPENPEPXL_BETA_PEPEV_PRE;

      /** User parameter name for PeptideEvidence info for the beta/acceptor peptide: post
              String
      */
      extern OPENMS_DLLAPI const std::string   OPENPEPXL_BETA_PEPEV_POST;

      /** User parameter name for PeptideEvidence info for the beta/acceptor peptide: start
              String
      */
      extern OPENMS_DLLAPI const std::string   OPENPEPXL_BETA_PEPEV_START;

      /** User parameter name for PeptideEvidence info for the beta/acceptor peptide: end
              String
      */
      extern OPENMS_DLLAPI const std::string   OPENPEPXL_BETA_PEPEV_END;

      /** User parameter name for XL-MS FDR values
              String
      */
      extern OPENMS_DLLAPI const std::string   XFDR_FDR;
    }

    //@}
  }
}
