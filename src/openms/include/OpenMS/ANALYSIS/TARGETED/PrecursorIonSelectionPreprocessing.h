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
// $Authors: $
// --------------------------------------------------------------------------
//

#ifndef OPENMS_ANALYSIS_TARGETED_PRECURSORIONSELECTIONPREPROCESSING_H
#define OPENMS_ANALYSIS_TARGETED_PRECURSORIONSELECTIONPREPROCESSING_H

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/KERNEL/FeatureMap.h>

#include <algorithm>
#include <cmath>
#include <set>
#include <fstream>
namespace OpenMS
{

  /**
    @brief This class implements the database preprocessing needing for precursor ion selection.

    @htmlinclude OpenMS_PrecursorIonSelectionPreprocessing.parameters
   */
  class OPENMS_DLLAPI PrecursorIonSelectionPreprocessing :
    public DefaultParamHandler
  {
public:
    PrecursorIonSelectionPreprocessing();
    PrecursorIonSelectionPreprocessing(const PrecursorIonSelectionPreprocessing & source);
    ~PrecursorIonSelectionPreprocessing() override;

    PrecursorIonSelectionPreprocessing & operator=(const PrecursorIonSelectionPreprocessing & source);


    const std::map<String, std::vector<double> > & getProtMasses() const;


    const std::vector<double> & getMasses(String acc) const;

    const std::map<String, std::vector<double> > & getProteinRTMap() const;
    const std::map<String, std::vector<double> > & getProteinPTMap() const;
    const std::map<String, std::vector<String> > & getProteinPeptideSequenceMap() const;


    /**
      @brief Calculates tryptic peptide masses of a given database and stores masses and peptide sequences

      @param db_path Path to database file (fasta)
      @param save Flag if preprocessing should be stored.

      @throws Exception::FileNotFound if file with preprocessing or db can't be found
      @throws Exception::UnableToCreateFile if preprocessing file can't be written
    */
    void dbPreprocessing(String db_path, bool save = true);

    /**
      @brief Calculates tryptic peptide masses of a given database and stores masses and peptide sequences

      @param db_path Path to database file (fasta)
      @param rt_model_path
      @param dt_model_path
      @param save Flag if preprocessing should be stored.

      @throws Exception::FileNotFound if file with preprocessing or db can't be found
      @throws Exception::UnableToCreateFile if preprocessing file can't be written
    */
    void dbPreprocessing(String db_path, String rt_model_path, String dt_model_path, bool save = true);


    /**
      @brief Loads tryptic peptide masses of a given database.

      @throws Exception::FileNotFound if file with preprocessing can't be found
      @throws Exception::InvalidParameter if precursor_mass_tolerance_unit is ppm and
      file containing bin borders can't be found
    */
    void loadPreprocessing();

    /// get the weighted frequency of a mass
    double getWeight(double mass);

    double getRT(String prot_id, Size peptide_index);

    double getPT(String prot_id, Size peptide_index);

    void setFixedModifications(StringList & modifications);
    const std::map<char, std::vector<String> > & getFixedModifications()
    {
      return fixed_modifications_;
    }

    void setGaussianParameters(double mu, double sigma);
    double getGaussMu()
    {
      return mu_;
    }

    double getGaussSigma()
    {
      return sigma_;
    }

    double getRTProbability(String prot_id, Size peptide_index, Feature & feature);
    double getRTProbability(double pred_rt, Feature & feature);

protected:
    /// saves the preprocessed db
    void savePreprocessedDB_(String db_path, String path);
    void savePreprocessedDBWithRT_(String db_path, String path);
    /// loads the preprocessed db
    void loadPreprocessedDB_(String path);
    /// pre-process fasta identifier
    void filterTaxonomyIdentifier_(FASTAFile::FASTAEntry & entry);
    Int getScanNumber_(double rt);
    double getRTProbability_(double min_obs_rt, double max_obs_rt, double pred_rt);
    /// update members method from DefaultParamHandler to update the members
    void updateMembers_() override;

    /// all tryptic masses of the distinct peptides in the database
    std::vector<double> masses_;
    /// the sequences of the tryptic peptides
    std::set<AASequence> sequences_;
    /// stores masses of tryptic peptides for proteins, key is the accession number
    std::map<String, std::vector<double> > prot_masses_;
    /// the masses of the bins used for preprocessing (only used if bins are not equidistant, i.e. with ppm)
    std::vector<double> bin_masses_;
    /// counter for the bins
    std::vector<UInt> counter_;
    /// maximal relative frequency of a mass
    UInt f_max_;

    bool fixed_mods_;
    std::map<String, std::vector<double> > rt_prot_map_;
    std::map<String, std::vector<double> > pt_prot_map_;
    std::map<String, std::vector<String> > prot_peptide_seq_map_;
    std::map<char, std::vector<String> > fixed_modifications_;
    double sigma_;
    double mu_;


  };
}

#endif //#ifndef OPENMS_ANALYSIS_ID_PRECURSORIONSELECTIONPREPROCESSING_H
