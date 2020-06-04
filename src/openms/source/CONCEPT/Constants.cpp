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

#include <OpenMS/CONCEPT/Constants.h>

namespace OpenMS
{
  namespace Constants
  {
    // EPSILON (used for comparisons)
    double EPSILON = 1e-6;

    // PI
    const double PI = 3.14159265358979323846;

    // Euler's number - base of the natural logarithm
    const double E = 2.718281828459045235;

    // Elementary charge
    const double ELEMENTARY_CHARGE = 1.60217738E-19; // C

    // Elementary charge (alias)
    const double e0 = ELEMENTARY_CHARGE;

    // Electron mass
    const double ELECTRON_MASS = 9.1093897E-31; // kg

    // Electron mass in units
    const double ELECTRON_MASS_U = 1.0 / 1822.8885020477; // u

    // Proton mass
    const double PROTON_MASS = 1.6726230E-27; // kg

    // Proton mass in units
    const double  PROTON_MASS_U = 1.0072764667710;// u

    // Mass difference between Carbon-13 and Carbon-12 in units
    const double C13C12_MASSDIFF_U = 1.0033548378; // u

    // Neutron mass
    const double NEUTRON_MASS = 1.6749286E-27; // kg

    // Neutron mass in units
    const double NEUTRON_MASS_U = 1.00866491566; // u

    // Avogadro constant
    const double AVOGADRO = 6.0221367E+23; // 1 / mol

    // Avogadro constant (alias)
    const double NA = AVOGADRO;

    // Avogadro constant (alias)
    const double MOL = AVOGADRO;

    // Boltzmann constant
    const double BOLTZMANN = 1.380657E-23; // J / K

    // Boltzmann constant (alias)
    const double k = BOLTZMANN;

    // Planck constant
    const double PLANCK = 6.6260754E-34; // J * sec

    // Planck constant (alias)
    const double h = PLANCK;

    // Gas constant (= NA * k)
    const double GAS_CONSTANT = NA * k;

    // Gas constant (alias)
    const double R = GAS_CONSTANT;

    // Faraday constant (= NA * e0)
    const double FARADAY = NA * e0;

    // Faraday constant (alias)
    const double F = FARADAY;

    // Bohr radius
    const double BOHR_RADIUS = 5.29177249E-11; // m

    // Bohr radius (alias)
    const double a0 = BOHR_RADIUS;

    //  the following values from:
    //  P.W.Atkins: Physical Chemistry, 5th ed., Oxford University Press, 1995

    // Vacuum permittivity
    const double VACUUM_PERMITTIVITY = 8.85419E-12; // C^2 / (J * m)

    // Vacuum permeability
    const double VACUUM_PERMEABILITY = (4 * PI * 1E-7); // J s^2 / (C^2 * m)

    // Speed of light
    const double SPEED_OF_LIGHT = 2.99792458E+8; // m / s

    // Speed of Light (alias)
    const double c = SPEED_OF_LIGHT;

    // Gravitational constant
    const double GRAVITATIONAL_CONSTANT  = 6.67259E-11; // N m^2 / kg^2

    // Fine structure constant
    const double FINE_STRUCTURE_CONSTANT = 7.29735E-3; // 1

    // Degree per rad
    const double DEG_PER_RAD = 57.2957795130823209;

    // Rad per degree
    const double RAD_PER_DEG = 0.0174532925199432957;

    // mm per inch
    const double MM_PER_INCH = 25.4;

    // m per foot
    const double M_PER_FOOT = 3.048;

    // Joule per calorie
    const double JOULE_PER_CAL  = 4.184;

    // Calories per Joule
    const double CAL_PER_JOULE = (1 / 4.184);

    namespace UserParam
    {
      // User parameter name for the M/Z of other chromatograms which have been merged into this one
      const std::string MERGED_CHROMATOGRAM_MZS = "merged_chromatogram_mzs";

      // User parameter name for precursor mz error in ppm
      const std::string PRECURSOR_ERROR_PPM_USERPARAM = "precursor_mz_error_ppm";

      // User parameter name for fragment mz error in ppm
      const std::string FRAGMENT_ERROR_MEDIAN_PPM_USERPARAM = "fragment_mz_error_median_ppm";

      // User parameter name for precursor mz error in ppm
      const std::string FRAGMENT_ANNOTATION_USERPARAM = "fragment_annotation";

      // User parameter name for the spectrum reference in PeptideIdentification (is is not yet treated as a class attribute)
      const std::string SPECTRUM_REFERENCE = "spectrum_reference";

      // User parameter name for target/decoy annotation of a PeptideHit, e.g. as annotated by PeptideIndexer. One of: target, decoy, target+decoy
      const std::string TARGET_DECOY = "target_decoy";

      // User parameter name for a delta score: a score ratio between a rank x hit and the rank x+1 hit
      const std::string DELTA_SCORE = "delta_score";

      // User parameter name for monoisotopic peak misassignment correction of a precursor (usually an integer x with the correction being -x times C13C12_MASSDIFF_U)
      const std::string ISOTOPE_ERROR = "isotope_error";

      // Cross-Linking Mass Spectrometry user parameters
      // Name of OpenPepXL main score (MzIdentML CV term)
      const std::string OPENPEPXL_SCORE = "OpenPepXL:score";

      // User parameter name for the sequence of the second peptide in a cross-link
      const std::string OPENPEPXL_BETA_SEQUENCE = "sequence_beta";

      // User parameter name for the protein accessions of the second peptide in a cross-link
      const std::string OPENPEPXL_BETA_ACCESSIONS = "accessions_beta";

      // User parameter name for the 1st position of cross-link (alpha peptide position in a real cross-link, 1st of two positions in a loop-link, modified position in a mono-link)
      const std::string OPENPEPXL_XL_POS1 = "xl_pos1";

      // User parameter name for the 2nd position of cross-link (beta peptide position in a real cross-link, 2nd of two positions in a loop-link, "-" in a mono-link)
      const std::string OPENPEPXL_XL_POS2 = "xl_pos2";

      // User parameter name for the 1st cross-link position on the protein
      const std::string OPENPEPXL_XL_POS1_PROT = "xl_pos1_protein";

      // User parameter name for the 1st cross-link position on the protein
      const std::string OPENPEPXL_XL_POS2_PROT = "xl_pos2_protein";

      // User parameter name for the cross-link type, one of: cross-link, loop-link, mono-link
      const std::string OPENPEPXL_XL_TYPE = "xl_type";

      // User parameter name for the cross-link rank (ranks of PeptideHits across different PeptideIdentifications)
      const std::string OPENPEPXL_XL_RANK = "xl_rank";

      // User parameter name for the name of a cross-link
      const std::string OPENPEPXL_XL_MOD = "xl_mod";

      // User parameter name for the mass of a cross-link
      const std::string OPENPEPXL_XL_MASS = "xl_mass";

      // User parameter name for the terminal specificity of a cross-link on the alpha peptide (to distinguish a link to the first or last residue side chain from a terminal link)
      const std::string OPENPEPXL_XL_TERM_SPEC_ALPHA = "xl_term_spec_alpha";

      // User parameter name for the terminal specificity of a cross-link on the beta peptide (to distinguish a link to the first or last residue side chain from a terminal link)
      const std::string OPENPEPXL_XL_TERM_SPEC_BETA = "xl_term_spec_beta";

      // User parameter name for the RT of the heavy spectrum precursor in a labeled cross-linking experiment
      const std::string OPENPEPXL_HEAVY_SPEC_RT = "spec_heavy_RT";

      // User parameter name for the m/z of the heavy spectrum precursor in a labeled cross-linking experiment
      const std::string OPENPEPXL_HEAVY_SPEC_MZ = "spec_heavy_MZ";

      // User parameter name for the spectrum reference of the heavy spectrum in a labeled cross-linking experiment
      const std::string OPENPEPXL_HEAVY_SPEC_REF = "spectrum_reference_heavy";

      // User parameter name for target/decoy annotation of alpha peptides
      const std::string OPENPEPXL_TARGET_DECOY_ALPHA = "xl_target_decoy_alpha";

      // User parameter name for target/decoy annotation of beta peptides
      const std::string OPENPEPXL_TARGET_DECOY_BETA = "xl_target_decoy_beta";

      // User parameter name for PeptideEvidence info for the beta/acceptor peptide: pre
      const std::string OPENPEPXL_BETA_PEPEV_PRE = "BetaPepEv:pre";

      // User parameter name for PeptideEvidence info for the beta/acceptor peptide: post
      const std::string OPENPEPXL_BETA_PEPEV_POST = "BetaPepEv:post";

      // User parameter name for PeptideEvidence info for the beta/acceptor peptide: start
      const std::string OPENPEPXL_BETA_PEPEV_START = "BetaPepEv:start";

      // User parameter name for PeptideEvidence info for the beta/acceptor peptide: end
      const std::string OPENPEPXL_BETA_PEPEV_END = "BetaPepEv:end";

      // User parameter name for XL-MS FDR values
      const std::string XFDR_FDR = "XFDR:FDR";
    }
  }
}
