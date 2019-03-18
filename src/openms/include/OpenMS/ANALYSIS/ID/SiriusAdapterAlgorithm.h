// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
// $Maintainer: Oliver Alka $
// $Authors: Oliver Alka $
// --------------------------------------------------------------------------

#pragma once 

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/ANALYSIS/QUANTITATION/KDTreeFeatureMaps.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureMapping.h>

#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>

#include <QString>

namespace OpenMS
{
  class OPENMS_DLLAPI SiriusAdapterAlgorithm : public DefaultParamHandler
    {
    public:
      /// default constructor
      SiriusAdapterAlgorithm();

      /// Struct for temporary folder structure
      struct SiriusTmpStruct
      {
        String tmp_dir;
        String tmp_ms_file;
        String tmp_out_dir;
      };

      /**
        @brief Construct temporary folder structure for SIRIUS (SiriusTmpStruct)

        @return SiriusTmpStruct
      */
      static SiriusAdapterAlgorithm::SiriusTmpStruct constructSiriusTmpStruct();

      /**
      @brief Checks if executable was povided 

      @return Pair "path to executable" and "path to the working directory"

      @param executable Path to the executable
      */
      static std::pair<String, String> checkSiriusExecutablePath(String& executable);


      /**
      @brief Preprocessing needed for SIRIUS

      @return FeatureToMS2Indices
            
      Filter number of masstraces and perform feature mapping.

      @param featureinfo: Path to featureXML
      @param spectra: Input of MSExperiment with spectra information
      @param fp_map_kd: KDTree used for query and match spectra with features
      @param sirius_algo: Parameters needed for preprocessing
      @param feature_mapping: Empty FeatureToMs2Indices
      */
      static void preprocessingSirius(const String& featureinfo,
                                      const MSExperiment& spectra,
                                      std::vector<FeatureMap>& v_fp,
                                      KDTreeFeatureMaps& fp_map_kd,
                                      const SiriusAdapterAlgorithm& sirius_algo,
                                      FeatureMapping::FeatureToMs2Indices& feature_mapping);

      /**
      @brief logs number of features and spectra used

      Prints the number of features and spectra used (LOG_INFO)

      @param featureinfo: Path to featureXML
      @param feature_mapping: FeatureToMs2Indices with feature mapping
      @param spectra: Input of MSExperiment with spectra information
      @param sirius_algo: Parameters needed for preprocessing
      */
      static void checkFeatureSpectraNumber(const String& featureinfo,
                                            const FeatureMapping::FeatureToMs2Indices& feature_mapping, 
                                            const MSExperiment& spectra, 
                                            const SiriusAdapterAlgorithm& sirius_algo);

      
      /**
      @brief Call SIRIUS with QProcess

      @return Vector with paths to a compound 

      @param tmp_ms_file: path to temporary .ms file
      @param tmp_out_dir: path to temporary output folder
      @param executable: path to executable
      @param out_csifingerid: path to CSI:FingerID output (can be empty).
      @param sirius_algo: Parameters needed for SIRIUS

      */
      const static std::vector<String> callSiriusQProcess(const String& tmp_ms_file,
                                                          const String& tmp_out_dir,
                                                          String& executable,
                                                          const String& out_csifingerid,
                                                          const SiriusAdapterAlgorithm& sirius_algo);

      // getter (used to call functions from SiriusMSConverter, SiriusMzTabWriter, CsiFingerIDMzTabWriter)
      String getFeatureOnly();  
      String getNoMasstraceInfoIsotopePattern();
      int getIsotopePatternIterations();
      int getCandidates();
      int getTopNHits();

    protected:

      // adapter parameters (preprocessing)
      unsigned int filter_by_num_masstraces_;
      double precursor_mz_tolerance_;
      String precursor_mz_tolerance_unit_;
      double precursor_rt_tolerance_;
      int isotope_pattern_iterations_;
      // flags
      String feature_only_;
      String no_masstrace_info_isotope_pattern_;
      // parameters for SIRIUS (sirius)
      String profile_;
      int candidates_;
      String database_;
      int noise_;
      int ppm_max_;
      String isotope_;
      String elements_;
      int compound_timeout_;
      int tree_timeout_;
      int top_n_hits_;
      int cores_;
      // flags
      String auto_charge_;
      String ion_tree_;
      String no_recalibration_;
      String most_intense_ms2_;

      void updateMembers_() override;
    };
} // namespace OpenMS
