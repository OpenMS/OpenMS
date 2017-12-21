// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------


#ifndef OPENMS_ANALYSIS_ID_PROTONDISTRIBUTIONMODEL_H
#define OPENMS_ANALYSIS_ID_PROTONDISTRIBUTIONMODEL_H

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <vector>

namespace OpenMS
{
  class AASequence;

  /**
      @brief A proton distribution model to calculate the proton distribution over charged peptides

      The model uses proton affinity values of backbone nitrogens and sidechains to calculate the
      proton distribution of charged peptide among these sites. The possible sites are the peptide
      bonds between the amino acids, the side chains and the C-terminus and N-terminus. The calculation
      is done calculating a Boltzmann distribution of the sites.

      Details and the proton affinities can be found in
      Z. Zhang, Prediction of Low-Energy Collision-Induced Dissociation Spectra of Peptides,
  Anal. Chem., 76 (14), 3908 - 3922, 2004

      A proton distribution can be calculated using the getProtonDistribution method. The backbone
      probabilities are reported in the first parameter (index 0 for the N-terminus, index 1 for the
      first peptide bond...), the site chain probabilities are reported in the second parameter
      (index 0, for the first amino acid...). The peptide and the number of protons as well as type
      of peptide (can be Reside::YIon for peptides and y-ions and any other ion type).

      Charge state intensities of differently charged equal (e.g. y7+ and y7++) ions can be calculated
      using the getChargeStateIntensities function.

      @htmlinclude OpenMS_ProtonDistributionModel.parameters
  */
  class OPENMS_DLLAPI ProtonDistributionModel :
    public DefaultParamHandler
  {
public:

    /** Constructor and destructors
    */
    //@{
    /// default constructor
    ProtonDistributionModel();

    /// copy constructor
    ProtonDistributionModel(const ProtonDistributionModel & model);

    /// destructor
    ~ProtonDistributionModel() override;
    //@}

    /// assignment operator
    ProtonDistributionModel & operator=(const ProtonDistributionModel & pdm);

    /** @name Enumerations
    */
    //@{
    /// the type of fragmentation
    enum FragmentationType
    {
      ChargeDirected = 0,
      ChargeRemote,
      SideChain
    };
    //@}

    /** @brief calculates a proton distribution of the given charged peptide

            @param bb_charges the calculated probabilities of the backbone sites (including N-terminus and C-terminus)
            @param sc_charges the calculated probabilities of the side chain sites
            @param peptide the peptide
            @param charge the charge
            @param res_type the type of the ion given in peptide. Peptides are handled as y-ions, i.e. Residue::YIon
    */
    void getProtonDistribution(std::vector<double> & bb_charges, std::vector<double> & sc_charges, const AASequence & peptide, Int charge, Residue::ResidueType res_type = Residue::YIon);

    /** @brief calculates the charge state intensities of different charge states of the same ion

            @param peptide the peptide
            @param n_term_ion the prefix ion sequence
            @param c_term_ion the suffix ion sequence
            @param charge the charge
            @param n_term_type the ion type of the N-terminal ion; valid values are Residue::AIon, Residue::BIon
            @param n_term_intensities the probability of seeing a charged prefix ions (first index corresponds to ion of charge 1)
            @param c_term_intensities the probability of seeing a charged suffix ions (first index corresponds to ion of charge 2)
            @param type the type of fragmentation (charge-directed, charge-remote of side chain)
    */
    //void getChargeStateIntensities(const AASequence& peptide, const AASequence& n_term_ion, const AASequence& c_term_ion, Int charge, Residue::ResidueType n_term_type, double& n_term1,  double& c_term1, double& n_term2, double& c_term2, FragmentationType type);


    void getChargeStateIntensities(const AASequence & peptide, const AASequence & n_term_ion, const AASequence & c_term_ion, Int charge, Residue::ResidueType n_term_type,
                                   std::vector<double> & n_term_intensities, std::vector<double> & c_term_intensities, FragmentationType type);

    /// sets the proton distributions of the whole peptide, they are needed for the getChargeStateIntensities_ method and need to be recalculated each time if not given
    void setPeptideProtonDistribution(const std::vector<double> & bb_charge, const std::vector<double> & sc_charge);

protected:

    // calculates the proton distribution
    void calculateProtonDistribution_(const AASequence & peptide, Int charge, Residue::ResidueType res_type = Residue::YIon, bool fixed_proton = false, Size cleavage_site = 0, bool use_most_basic_site = false);

    void calculateProtonDistributionCharge1_(const AASequence & peptide, Residue::ResidueType res_type);
    void calculateProtonDistributionCharge2_(const AASequence & peptide, Residue::ResidueType res_type, bool fixed_proton, Size cleavage_site, bool use_most_basic_site);
    void calculateProtonDistributionGreater2_(const AASequence & peptide, Int charge, Residue::ResidueType res_type);

    void calculateProtonDistributionIonPair_(const AASequence & peptide, Residue::ResidueType type, Size cleavage_site);

    // returns the proton affinity of the peptide with the given charge and ion type
    //double getProtonAffinity_(const AASequence& ion, Int charge, Residue::ResidueType res_type);

    // returns the (relative) Intensities of the possible charge states of the ion from peptide
    //std::vector<double> getChargeStateIntensities_(const AASequence& peptide, const AASequence& ion, Int charge, Residue::ResidueType res_type);


    void calcChargeStateIntensities_(const AASequence & peptide, const AASequence & n_term_ion, const AASequence & c_term_ion, Int charge,
                                     Residue::ResidueType n_term_type, std::vector<double> & n_term_intensities, std::vector<double> & c_term_intensities, FragmentationType type);

    // calculates the intensities of the different possible charge states
    void calcChargeStateIntensities_(const AASequence & peptide, const AASequence & n_term_ion, const AASequence & c_term_ion, Int charge, Residue::ResidueType n_term_type, double & n_term1, double & c_term1, double & n_term2, double & c_term2, FragmentationType type);

    // returns the left and right GB values, NH2 and COOH if at terminus
    void getLeftAndRightGBValues_(const AASequence & peptide, double & left_gb, double & right_gb, Size position);

    std::vector<double> sc_charge_;
    std::vector<double> bb_charge_;
    std::vector<double> sc_charge_full_;
    std::vector<double> bb_charge_full_;
    std::vector<double> sc_charge_ion_n_term_;
    std::vector<double> bb_charge_ion_n_term_;
    std::vector<double> sc_charge_ion_c_term_;
    std::vector<double> bb_charge_ion_c_term_;
    double E_;
    double E_c_term_;
    double E_n_term_;

  };
}
#endif // OPENMS_ANALYSIS_ID_PROTONDISTRIBUTIONMODEL_H
