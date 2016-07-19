// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2016.
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
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

// Interfaces
#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/ISpectrumAccess.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/DataStructures.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>
#include <OpenMS/INTERFACES/IMSDataConsumer.h>

// Consumers
#include <OpenMS/FORMAT/DATAACCESS/MSDataCachedConsumer.h>
#include <OpenMS/FORMAT/DATAACCESS/MSDataTransformingConsumer.h>
#include <OpenMS/FORMAT/DATAACCESS/MSDataWritingConsumer.h>

// Files
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/TraMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/TransformationXMLFile.h>
#include <OpenMS/ANALYSIS/OPENSWATH/TransitionTSVReader.h>
#include <OpenMS/FORMAT/CachedMzML.h>
#include <OpenMS/FORMAT/SwathFile.h>
#include <OpenMS/ANALYSIS/OPENSWATH/SwathWindowLoader.h>

// Kernel and implementations
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessOpenMS.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessTransforming.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessOpenMSInMemory.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/SwathMap.h>

// Helpers
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathHelper.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/DataAccessHelper.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SimpleOpenMSSpectraAccessFactory.h>

// Algorithms
#include <OpenMS/ANALYSIS/OPENSWATH/MRMRTNormalizer.h>
#include <OpenMS/ANALYSIS/OPENSWATH/ChromatogramExtractor.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMFeatureFinderScoring.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMTransitionGroupPicker.h>
#include <OpenMS/ANALYSIS/OPENSWATH/SwathMapMassCorrection.h>

#include <assert.h>

#define OPENSWATH_WORKFLOW_DEBUG

using namespace OpenMS;

static bool SortPairDoubleByFirst(const std::pair<double,double> & left, const std::pair<double,double> & right)
{
  return left.first < right.first;
}

// The workflow class and the TSV writer
namespace OpenMS
{

  /**
   * @brief Class to write out an OpenSwath TSV output (mProphet input)
   *
   */
  class OpenSwathTSVWriter
  {
    std::ofstream ofs;
    String input_filename_;
    bool doWrite_;
    bool use_ms1_traces_;
    bool enable_uis_scoring_;

  public:

    OpenSwathTSVWriter(String output_filename, String input_filename = "inputfile", bool ms1_scores = false, bool uis_scores = false) :
      ofs(output_filename.c_str()),
      input_filename_(input_filename),
      doWrite_(!output_filename.empty()),
      use_ms1_traces_(ms1_scores),
      enable_uis_scoring_(uis_scores)
      {}

    bool isActive() 
    {
      return doWrite_;
    }

    void writeHeader()
    {
      ofs << "transition_group_id\tpeptide_group_label\trun_id\tfilename\tRT\tid\tSequence\tFullPeptideName" <<
        "\tCharge\tm/z\tIntensity\tProteinName\tdecoy\tassay_rt\tdelta_rt\tleftWidth" <<
        "\tmain_var_xx_swath_prelim_score\tnorm_RT\tnr_peaks\tpeak_apices_sum\tpotentialOutlier" <<
        "\trightWidth\trt_score\tsn_ratio\ttotal_xic\tvar_bseries_score\tvar_dotprod_score" <<
        "\tvar_intensity_score\tvar_isotope_correlation_score\tvar_isotope_overlap_score" <<
        "\tvar_library_corr\tvar_library_dotprod\tvar_library_manhattan\tvar_library_rmsd" <<
        "\tvar_library_rootmeansquare\tvar_library_sangle\tvar_log_sn_score\tvar_manhatt_score" <<
        "\tvar_massdev_score\tvar_massdev_score_weighted\tvar_norm_rt_score\tvar_xcorr_coelution" <<
        "\tvar_xcorr_coelution_weighted\tvar_xcorr_shape\tvar_xcorr_shape_weighted" <<
        "\tvar_yseries_score\tvar_elution_model_fit_score";
      if (use_ms1_traces_)
      {
        ofs << "\tvar_ms1_ppm_diff\tvar_ms1_isotope_corr\tvar_ms1_isotope_overlap\tvar_ms1_xcorr_coelution\tvar_ms1_xcorr_shape";
      }
      ofs << "\txx_lda_prelim_score\txx_swath_prelim_score";
      if (use_ms1_traces_)
      {
        ofs << "\taggr_prec_Peak_Area\taggr_prec_Peak_Apex\taggr_prec_Fragment_Annotation";
      }
      ofs << "\taggr_Peak_Area\taggr_Peak_Apex\taggr_Fragment_Annotation";
      if (enable_uis_scoring_)
      {
        ofs << "\tuis_target_transition_names"
            << "\tuis_target_var_ind_log_intensity"
            << "\tuis_target_num_transitions"
            << "\tuis_target_var_ind_xcorr_coelution"
            << "\tuis_target_main_var_ind_xcorr_shape"
            << "\tuis_target_var_ind_log_sn_score"
            << "\tuis_target_var_ind_massdev_score"
            << "\tuis_target_var_ind_isotope_correlation"
            << "\tuis_target_var_ind_isotope_overlap"
            << "\tuis_decoy_transition_names"
            << "\tuis_decoy_var_ind_log_intensity"
            << "\tuis_decoy_num_transitions"
            << "\tuis_decoy_var_ind_xcorr_coelution"
            << "\tuis_decoy_main_var_ind_xcorr_shape"
            << "\tuis_decoy_var_ind_log_sn_score"
            << "\tuis_decoy_var_ind_massdev_score"
            << "\tuis_decoy_var_ind_isotope_correlation"
            << "\tuis_decoy_var_ind_isotope_overlap";
      }
      ofs << "\n";
    }

    String prepareLine(const OpenSwath::LightCompound& pep,
        const OpenSwath::LightTransition* transition,
        FeatureMap& output, String id)
    {
        String result = "";
        String decoy = "0"; // 0 = false
        if (transition->decoy) 
        {
          decoy = "1";
        }

        for (FeatureMap::iterator feature_it = output.begin(); feature_it != output.end(); ++feature_it)
        {

          char intensity_char[40];
          char intensity_apex_char[40];
          String aggr_Peak_Area = "";
          String aggr_Peak_Apex = "";
          String aggr_Fragment_Annotation = "";
          String aggr_prec_Peak_Area = "";
          String aggr_prec_Peak_Apex = "";
          String aggr_prec_Fragment_Annotation = "";
          for (std::vector<Feature>::iterator sub_it = feature_it->getSubordinates().begin(); sub_it != feature_it->getSubordinates().end(); ++sub_it)
          {
            sprintf(intensity_char, "%f", sub_it->getIntensity());
            sprintf(intensity_apex_char, "%f", (double)sub_it->getMetaValue("peak_apex_int"));
            if (sub_it->metaValueExists("FeatureLevel") && sub_it->getMetaValue("FeatureLevel") == "MS2")
            {
              aggr_Peak_Area += (String)intensity_char + ";";
              aggr_Peak_Apex += (String)intensity_apex_char + ";";
              aggr_Fragment_Annotation += (String)sub_it->getMetaValue("native_id") + ";";
            }
            else if (sub_it->metaValueExists("FeatureLevel") && sub_it->getMetaValue("FeatureLevel") == "MS1")
            {
              aggr_prec_Peak_Area += (String)intensity_char + ";";
              aggr_Peak_Apex += (String)intensity_apex_char + ";";
              aggr_prec_Fragment_Annotation += (String)sub_it->getMetaValue("native_id") + ";";
            }
          }
          if (!feature_it->getSubordinates().empty())
          {
            aggr_Peak_Area = aggr_Peak_Area.substr(0, aggr_Peak_Area.size() - 1);
            aggr_Peak_Apex = aggr_Peak_Apex.substr(0, aggr_Peak_Apex.size() - 1);
            aggr_Fragment_Annotation = aggr_Fragment_Annotation.substr(0, aggr_Fragment_Annotation.size() - 1);
            aggr_prec_Peak_Area = aggr_prec_Peak_Area.substr(0, aggr_prec_Peak_Area.size() - 1);
            aggr_prec_Peak_Apex = aggr_prec_Peak_Apex.substr(0, aggr_prec_Peak_Apex.size() - 1);
            aggr_prec_Fragment_Annotation = aggr_prec_Fragment_Annotation.substr(0, aggr_prec_Fragment_Annotation.size() - 1);
          }

          String full_peptide_name = "";
          for (int loc = -1; loc <= (int)pep.sequence.size(); loc++)
          {
            if (loc > -1 && loc < (int)pep.sequence.size())
            {
              full_peptide_name += pep.sequence[loc];
            }
            // C-terminal and N-terminal modifications may be at positions -1 or pep.sequence
            for (Size modloc = 0; modloc < pep.modifications.size(); modloc++)
            {
              if (pep.modifications[modloc].location == loc)
              {
                full_peptide_name += "(" + pep.modifications[modloc].unimod_id + ")";
              }
            }
          }

          // Compute peptide group label (use the provided label or use the
          // transition group).
          String group_label = pep.peptide_group_label;
          if (group_label.empty()) group_label = id;
          if (group_label == "light") group_label = id; // legacy fix since there are many TraMLs floating around which have "light" in there

          // If a protein is present, take the first one
          String protein_name = "";
          if (!pep.protein_refs.empty() )
          {
            protein_name = pep.protein_refs[0];
          }

          String line = "";
          line += id + "_run0"
            + "\t" + group_label
            + "\t" + "0"
            + "\t" + input_filename_
            + "\t" + (String)feature_it->getRT()
            + "\t" + "f_" + feature_it->getUniqueId()  // TODO might not be unique!!!
            + "\t" + (String)pep.sequence
            + "\t" + full_peptide_name
            + "\t" + (String)pep.charge
            + "\t" + (String)transition->precursor_mz
            + "\t" + (String)feature_it->getIntensity()
            + "\t" + protein_name 
            + "\t" + decoy
            // Note: missing MetaValues will just produce a DataValue::EMPTY which lead to an empty column
            + "\t" + (String)feature_it->getMetaValue("assay_rt")
            + "\t" + (String)feature_it->getMetaValue("delta_rt")
            + "\t" + (String)feature_it->getMetaValue("leftWidth")
            + "\t" + (String)feature_it->getMetaValue("main_var_xx_swath_prelim_score")
            + "\t" + (String)feature_it->getMetaValue("norm_RT")
            + "\t" + (String)feature_it->getMetaValue("nr_peaks")
            + "\t" + (String)feature_it->getMetaValue("peak_apices_sum")
            + "\t" + (String)feature_it->getMetaValue("potentialOutlier")
            + "\t" + (String)feature_it->getMetaValue("rightWidth")
            + "\t" + (String)feature_it->getMetaValue("rt_score")
            + "\t" + (String)feature_it->getMetaValue("sn_ratio")
            + "\t" + (String)feature_it->getMetaValue("total_xic")
            + "\t" + (String)feature_it->getMetaValue("var_bseries_score")
            + "\t" + (String)feature_it->getMetaValue("var_dotprod_score")
            + "\t" + (String)feature_it->getMetaValue("var_intensity_score")
            + "\t" + (String)feature_it->getMetaValue("var_isotope_correlation_score")
            + "\t" + (String)feature_it->getMetaValue("var_isotope_overlap_score")
            + "\t" + (String)feature_it->getMetaValue("var_library_corr")
            + "\t" + (String)feature_it->getMetaValue("var_library_dotprod")
            + "\t" + (String)feature_it->getMetaValue("var_library_manhattan")
            + "\t" + (String)feature_it->getMetaValue("var_library_rmsd")
            + "\t" + (String)feature_it->getMetaValue("var_library_rootmeansquare")
            + "\t" + (String)feature_it->getMetaValue("var_library_sangle")
            + "\t" + (String)feature_it->getMetaValue("var_log_sn_score")
            + "\t" + (String)feature_it->getMetaValue("var_manhatt_score")
            + "\t" + (String)feature_it->getMetaValue("var_massdev_score")
            + "\t" + (String)feature_it->getMetaValue("var_massdev_score_weighted")
            + "\t" + (String)feature_it->getMetaValue("var_norm_rt_score")
            + "\t" + (String)feature_it->getMetaValue("var_xcorr_coelution")
            + "\t" + (String)feature_it->getMetaValue("var_xcorr_coelution_weighted")
            + "\t" + (String)feature_it->getMetaValue("var_xcorr_shape")
            + "\t" + (String)feature_it->getMetaValue("var_xcorr_shape_weighted")
            + "\t" + (String)feature_it->getMetaValue("var_yseries_score")
            + "\t" + (String)feature_it->getMetaValue("var_elution_model_fit_score");

            if (use_ms1_traces_) 
            {
              line += "\t" + (String)feature_it->getMetaValue("var_ms1_ppm_diff")
              + "\t" + (String)feature_it->getMetaValue("var_ms1_isotope_correlation")
              + "\t" + (String)feature_it->getMetaValue("var_ms1_isotope_overlap")
              + "\t" + (String)feature_it->getMetaValue("var_ms1_xcorr_coelution")
              + "\t" + (String)feature_it->getMetaValue("var_ms1_xcorr_shape");
            }

            line += "\t" + (String)feature_it->getMetaValue("xx_lda_prelim_score")
            + "\t" + (String)feature_it->getMetaValue("xx_swath_prelim_score");
            if (use_ms1_traces_) 
            {
              line += "\t" + aggr_prec_Peak_Area + "\t" + aggr_prec_Peak_Apex + "\t" + aggr_prec_Fragment_Annotation;
            }
            line += "\t" + aggr_Peak_Area + "\t" + aggr_Peak_Apex + "\t" + aggr_Fragment_Annotation;
            if (enable_uis_scoring_)
            {
              line += "\t" + (String)feature_it->getMetaValue("id_target_transition_names")
              + "\t" + (String)feature_it->getMetaValue("id_target_ind_log_intensity")
              + "\t" + (String)feature_it->getMetaValue("id_target_num_transitions")
              + "\t" + (String)feature_it->getMetaValue("id_target_ind_xcorr_coelution")
              + "\t" + (String)feature_it->getMetaValue("id_target_ind_xcorr_shape")
              + "\t" + (String)feature_it->getMetaValue("id_target_ind_log_sn_score")
              + "\t" + (String)feature_it->getMetaValue("id_target_ind_massdev_score")
              + "\t" + (String)feature_it->getMetaValue("id_target_ind_isotope_correlation")
              + "\t" + (String)feature_it->getMetaValue("id_target_ind_isotope_overlap")
              + "\t" + (String)feature_it->getMetaValue("id_decoy_transition_names")
              + "\t" + (String)feature_it->getMetaValue("id_decoy_ind_log_intensity")
              + "\t" + (String)feature_it->getMetaValue("id_decoy_num_transitions")
              + "\t" + (String)feature_it->getMetaValue("id_decoy_ind_xcorr_coelution")
              + "\t" + (String)feature_it->getMetaValue("id_decoy_ind_xcorr_shape")
              + "\t" + (String)feature_it->getMetaValue("id_decoy_ind_log_sn_score")
              + "\t" + (String)feature_it->getMetaValue("id_decoy_ind_massdev_score")
              + "\t" + (String)feature_it->getMetaValue("id_decoy_ind_isotope_correlation")
              + "\t" + (String)feature_it->getMetaValue("id_decoy_ind_isotope_overlap");
            }
            line += "\n";          result += line;
        } // end of iteration
      return result;
    }

    void writeLines(std::vector<String> to_output)
    {
      for (Size i = 0; i < to_output.size(); i++) { ofs << to_output[i]; }
    }

  };

  /**
   * @brief Class to execute an OpenSwath Workflow
   *
   * performExtraction will perform the OpenSWATH analysis. Optionally, an RT
   * transformation (mapping peptides to normalized space) can be obtained
   * beforehand using performRTNormalization.
   *
   */
  class OpenSwathWorkflow :
    public ProgressLogger
  {

  public:

    explicit OpenSwathWorkflow(bool use_ms1_traces) :
      use_ms1_traces_(use_ms1_traces)
    {}

    /** @brief ChromatogramExtractor parameters
     *
     * A small helper struct to pass the parameters for the chromatogram
     * extraction through to the actual algorithm.
     *
    */
    struct ChromExtractParams
    {
      /// Whether to not extract anything closer than this (in Da) from the upper edge
      double min_upper_edge_dist;
      /// Extraction window in Da or ppm (e.g. 50ppm means extraction +/- 25ppm)
      double mz_extraction_window;
      /// Whether the extraction window is given in ppm or Da
      bool ppm;
      /// The extraction function in mass space
      String extraction_function;
      /// The retention time extraction window
      double rt_extraction_window;
      /// Whether to extract some extra in the retention time (can be useful if one wants to look at the chromatogram outside the window)
      double extra_rt_extract;
    };

    /** @brief Compute the alignment against a set of RT-normalization peptides
     *
     * This function extracts the RT normalization chromatograms
     * (simpleExtractChromatograms) and then uses the chromatograms to find
     * features (in RTNormalization).
     *
     * @param irt_transitions A set of transitions used for the RT normalization peptides
     * @param swath_maps The raw data (swath maps)
     * @param min_rsq Minimal R^2 value that is expected for the RT regression
     * @param min_coverage Minimal coverage of the chromatographic space that needs to be achieved
     * @param feature_finder_param Parameter set for the feature finding in chromatographic dimension 
     * @param cp_irt Parameter set for the chromatogram extraction
     * @param irt_detection_param Parameter set for the detection of the iRTs (outlier detection, peptides per bin etc)
     * @param mz_correction_function If correction in m/z is desired, which function should be used
     * @param debug_level Debug level (writes out the RT normalization chromatograms if larger than 1)
     *
    */
    TransformationDescription performRTNormalization(const OpenMS::TargetedExperiment & irt_transitions,
            std::vector< OpenSwath::SwathMap > & swath_maps, double min_rsq, double min_coverage,
            const Param & feature_finder_param, const ChromExtractParams & cp_irt, const Param& irt_detection_param, 
            const String & mz_correction_function, Size debug_level)
    {
      LOG_DEBUG << "performRTNormalization method starting" << std::endl;
      std::vector< OpenMS::MSChromatogram<> > irt_chromatograms;
      simpleExtractChromatograms(swath_maps, irt_transitions, irt_chromatograms, cp_irt);

      // debug output of the iRT chromatograms
      if (debug_level > 1)
      {
        try
        {
          MSExperiment<> exp;
          exp.setChromatograms(irt_chromatograms);
          MzMLFile().store("debug_irts.mzML", exp);
        }
        catch (OpenMS::Exception::UnableToCreateFile& /*e*/)
        {
          LOG_DEBUG << "Error creating file 'debug_irts.mzML', not writing out iRT chromatogram file"  << std::endl;
        }
        catch (OpenMS::Exception::BaseException& /*e*/)
        {
          LOG_DEBUG << "Error writint to file 'debug_irts.mzML', not writing out iRT chromatogram file"  << std::endl;
        }
      }
      LOG_DEBUG << "Extracted number of chromatograms from iRT files: " << irt_chromatograms.size() <<  std::endl;

      // get RT normalization from data
      TransformationDescription tr = RTNormalization(irt_transitions,
              irt_chromatograms, min_rsq, min_coverage, feature_finder_param, irt_detection_param, swath_maps, mz_correction_function);
      return tr;
    }

    /** @brief Execute the OpenSWATH workflow on a set of SwathMaps and transitions.
     *
     * Executes the following operations on the given input:
     *
     * 1. OpenSwathHelper::selectSwathTransitions
     * 2. ChromatogramExtractor prepare, extract
     * 3. scoreAllChromatograms
     * 4. Write out chromatograms and found features
     *
     * @param swath_maps The raw data (swath maps)
     * @param trafo Transformation description (translating this runs' RT to normalized RT space)
     * @param cp Parameter set for the chromatogram extraction
     * @param feature_finder_param Parameter set for the feature finding in chromatographic dimension 
     * @param transition_exp The set of assays to be extracted and scored
     * @param out_featureFile Output feature map to store identified features
     * @param store_features Whether features should be appended to the output feature map
     * @param tsv_writer TSV Writer object to store identified features in csv format
     * @param chromConsumer Chromatogram consumer object to store the extracted chromatograms 
     * @param batchSize Size of the batches which should be extracted and scored
     *
    */
    void performExtraction(const std::vector< OpenSwath::SwathMap > & swath_maps,
      const TransformationDescription trafo,
      const ChromExtractParams & cp, const Param & feature_finder_param,
      const OpenSwath::LightTargetedExperiment& transition_exp,
      FeatureMap& out_featureFile, bool store_features,
      OpenSwathTSVWriter & tsv_writer, Interfaces::IMSDataConsumer<> * chromConsumer,
      int batchSize, bool load_into_memory)
    {
      tsv_writer.writeHeader();

      // Compute inversion of the transformation
      TransformationDescription trafo_inverse = trafo;
      trafo_inverse.invert();

      std::cout << "Will analyze " << transition_exp.transitions.size() << " transitions in total." << std::endl;
      int progress = 0;
      this->startProgress(0, swath_maps.size(), "Extracting and scoring transitions");

      // (i) Obtain precursor chromatograms (MS1) if precursor extraction is enabled
      std::map< std::string, OpenSwath::ChromatogramPtr > ms1_chromatograms;
      for (SignedSize i = 0; i < boost::numeric_cast<SignedSize>(swath_maps.size()); ++i)
      {
        if (swath_maps[i].ms1 && use_ms1_traces_) 
        {
          // store reference to MS1 map for later -> note that this is *not* threadsafe!
          ms1_map_ = swath_maps[i].sptr;

          if (load_into_memory)
          {
            // This creates an InMemory object that keeps all data in memory
            // but provides the same access functionality to the raw data as
            // any object implementing ISpectrumAccess 
            ms1_map_ = boost::shared_ptr<SpectrumAccessOpenMSInMemory>( new SpectrumAccessOpenMSInMemory(*ms1_map_) );
          }

          std::vector< OpenSwath::ChromatogramPtr > chrom_list;
          std::vector< ChromatogramExtractor::ExtractionCoordinates > coordinates;
          OpenSwath::LightTargetedExperiment transition_exp_used = transition_exp; // copy for const correctness
          ChromatogramExtractor extractor;

          // prepare the extraction coordinates & extract chromatogram
          prepare_coordinates_wrap(chrom_list, coordinates, transition_exp_used, true, trafo_inverse, cp);
          extractor.extractChromatograms(ms1_map_, chrom_list, coordinates, cp.mz_extraction_window,
              cp.ppm, cp.extraction_function);

          std::vector< OpenMS::MSChromatogram<> > chromatograms;
          extractor.return_chromatogram(chrom_list, coordinates, transition_exp_used,  SpectrumSettings(), chromatograms, true);

          for (Size j = 0; j < coordinates.size(); j++)
          {
            if ( chromatograms[j].empty() ) {continue;} // skip empty chromatograms
            ms1_chromatograms [ coordinates[j].id ] = chrom_list[j];
            // write MS1 chromatograms to disk
            chromConsumer->consumeChromatogram( chromatograms[j] );
          }
        }
      }

      // (ii) Perform extraction and scoring of fragment ion chromatograms (MS2)
      // We set dynamic scheduling such that the maps are worked on in the order
      // in which they were given to the program / acquired. This gives much
      // better load balancing than static allocation.
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1)
#endif
      for (SignedSize i = 0; i < boost::numeric_cast<SignedSize>(swath_maps.size()); ++i)
      {
        if (!swath_maps[i].ms1) // skip MS1
        {

          OpenSwath::SpectrumAccessPtr current_swath_map = swath_maps[i].sptr;

          if (load_into_memory)
          {
            // This creates an InMemory object that keeps all data in memory
            // but provides the same access functionality to the raw data as
            // any object implementing ISpectrumAccess 
            current_swath_map = boost::shared_ptr<SpectrumAccessOpenMSInMemory>( new SpectrumAccessOpenMSInMemory(*current_swath_map) );
            std::cout << " loading all data completely into memory !!! " << std::endl;
          }

          // Step 1: select which transitions to extract (proceed in batches)
          OpenSwath::LightTargetedExperiment transition_exp_used_all;
          OpenSwathHelper::selectSwathTransitions(transition_exp, transition_exp_used_all,
              cp.min_upper_edge_dist, swath_maps[i].lower, swath_maps[i].upper);
          if (transition_exp_used_all.getTransitions().size() > 0) // skip if no transitions found
          {

            int batch_size;
            if (batchSize <= 0 || batchSize >= (int)transition_exp_used_all.getCompounds().size())
            {
              batch_size = transition_exp_used_all.getCompounds().size();
            }
            else
            {
              batch_size = batchSize;
            }

#ifdef _OPENMP
#pragma omp critical (featureFinder)
#endif
            {
              std::cout << "Thread " <<
#ifdef _OPENMP
              omp_get_thread_num() << " " <<
#endif
              "will analyze " << transition_exp_used_all.getCompounds().size() <<  " compounds and "
              << transition_exp_used_all.getTransitions().size() <<  " transitions "
              "from SWATH " << i << " in batches of " << batch_size << std::endl;
            }

            for (size_t pep_idx = 0; pep_idx <= (transition_exp_used_all.getCompounds().size() / batch_size); pep_idx++)
            {
              // Create the new, batch-size transition experiment
              OpenSwath::LightTargetedExperiment transition_exp_used;
              selectCompoundsForBatch_(transition_exp_used_all, transition_exp_used, batch_size, pep_idx);

              // Step 2.1: extract these transitions
              ChromatogramExtractor extractor;
              boost::shared_ptr<MSExperiment<Peak1D> > chrom_exp(new MSExperiment<Peak1D>);
              std::vector< OpenSwath::ChromatogramPtr > chrom_list;
              std::vector< ChromatogramExtractor::ExtractionCoordinates > coordinates;

              // Step 2.2: prepare the extraction coordinates & extract chromatograms
              prepare_coordinates_wrap(chrom_list, coordinates, transition_exp_used, false, trafo_inverse, cp);
              extractor.extractChromatograms(current_swath_map, chrom_list, coordinates, cp.mz_extraction_window,
                  cp.ppm, cp.extraction_function);

              // Step 2.3: convert chromatograms back and write to output
              std::vector< OpenMS::MSChromatogram<> > chromatograms;
              extractor.return_chromatogram(chrom_list, coordinates, transition_exp_used,  SpectrumSettings(), chromatograms, false);
              chrom_exp->setChromatograms(chromatograms);
              OpenSwath::SpectrumAccessPtr chromatogram_ptr = OpenSwath::SpectrumAccessPtr(new OpenMS::SpectrumAccessOpenMS(chrom_exp));

              // Step 3: score these extracted transitions
              FeatureMap featureFile;
              scoreAllChromatograms(chromatogram_ptr, current_swath_map, transition_exp_used,
                  feature_finder_param, trafo, cp.rt_extraction_window, featureFile, tsv_writer, 
                  ms1_chromatograms);

              // Step 4: write all chromatograms and features out into an output object / file
              // (this needs to be done in a critical section since we only have one
              // output file and one output map).
#ifdef _OPENMP
#pragma omp critical (featureFinder)
#endif
              {
                // write chromatograms to output if so desired
                for (Size chrom_idx = 0; chrom_idx < chromatograms.size(); ++chrom_idx)
                {
                  chromConsumer->consumeChromatogram(chromatograms[chrom_idx]);
                }

                // write features to output if so desired
                if (store_features)
                {
                  for (FeatureMap::iterator feature_it = featureFile.begin();
                       feature_it != featureFile.end(); ++feature_it)
                  {
                    out_featureFile.push_back(*feature_it);
                  }
                  for (std::vector<ProteinIdentification>::iterator protid_it =
                         featureFile.getProteinIdentifications().begin();
                       protid_it != featureFile.getProteinIdentifications().end();
                       ++protid_it)
                  {
                    out_featureFile.getProteinIdentifications().push_back(*protid_it);
                  }
                  this->setProgress(progress++);
                }
              }
            }

          } // continue 2 (no continue due to OpenMP)
        } // continue 1 (no continue due to OpenMP)
      }
      this->endProgress();
    }

  private:

    /** @brief Select which compounds to analyze in the next batch (and copy to output)
     *
     * This function will select which compounds to analyze in the next batch j
     * and will copy the the corresponding compounds and transitions into the
     * output structure. The output will contain batch_size compounds.
     *
     * @param transition_exp_used_all input (all transitions for this swath)
     * @param transition_exp_used output (will contain only transitions for the next batch)
     * @param batch_size how many compounds per batch
     * @param j batch number (peptides from j*batch_size to j*batch_size+batch_size will be copied)
     *
     * @note The proteins will be copied completely without checking for a match
     *
    */
    void selectCompoundsForBatch_(const OpenSwath::LightTargetedExperiment& transition_exp_used_all,
      OpenSwath::LightTargetedExperiment& transition_exp_used, int batch_size, size_t j)
    {
      // compute batch start/end
      size_t start = j * batch_size;
      size_t end = j * batch_size + batch_size;
      if (end > transition_exp_used_all.compounds.size())
      {
        end = transition_exp_used_all.compounds.size();
      }

      // Create the new, batch-size transition experiment
      transition_exp_used.proteins = transition_exp_used_all.proteins;
      transition_exp_used.compounds.insert(transition_exp_used.compounds.end(),
          transition_exp_used_all.compounds.begin() + start, transition_exp_used_all.compounds.begin() + end);
      copyBatchTransitions_(transition_exp_used.compounds, transition_exp_used_all.transitions, transition_exp_used.transitions);
    }

    /** @brief Copy the required transitions to output
     *
     * Copy all transitions matching to one of the compounds in the selected
     * peptide vector from all_transitions to the output.
     *
     * @param used_compounds Which peptides or metabolites to be used
     * @param all_transitions Transitions vector from which to select transitions
     * @param output Output vector containing matching transitions (taken from all_transitions)
     *
    */
    void copyBatchTransitions_(const std::vector<OpenSwath::LightCompound>& used_compounds,
      const std::vector<OpenSwath::LightTransition>& all_transitions,
      std::vector<OpenSwath::LightTransition>& output)
    {
      std::set<std::string> selected_compounds;
      for (Size i = 0; i < used_compounds.size(); i++)
      {
        selected_compounds.insert(used_compounds[i].id);
      }

      for (Size i = 0; i < all_transitions.size(); i++)
      {
        if (selected_compounds.find(all_transitions[i].peptide_ref) != selected_compounds.end())
        {
          output.push_back(all_transitions[i]);
        }
      }
    }

    /// Simple method to extract chromatograms (for the RT-normalization peptides)
    void simpleExtractChromatograms(const std::vector< OpenSwath::SwathMap > & swath_maps,
      const OpenMS::TargetedExperiment & irt_transitions,
      std::vector< OpenMS::MSChromatogram<> > & chromatograms, const ChromExtractParams & cp)
    {
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for (SignedSize map_idx = 0; map_idx < boost::numeric_cast<SignedSize>(swath_maps.size()); ++map_idx)
      {
        std::vector< OpenMS::MSChromatogram<> > tmp_chromatograms;
        if (!swath_maps[map_idx].ms1) // skip MS1
        {

          TargetedExperiment transition_exp_used;
          OpenSwathHelper::selectSwathTransitions(irt_transitions, transition_exp_used,
              cp.min_upper_edge_dist, swath_maps[map_idx].lower, swath_maps[map_idx].upper);
          if (transition_exp_used.getTransitions().size() > 0) // skip if no transitions found
          {

            std::vector< OpenSwath::ChromatogramPtr > tmp_out;
            std::vector< ChromatogramExtractor::ExtractionCoordinates > coordinates;
            ChromatogramExtractor extractor;
            extractor.prepare_coordinates(tmp_out, coordinates, transition_exp_used, cp.rt_extraction_window, false);
            extractor.extractChromatograms(swath_maps[map_idx].sptr, tmp_out, coordinates, cp.mz_extraction_window,
                cp.ppm, cp.extraction_function);
            extractor.return_chromatogram(tmp_out, coordinates,
                transition_exp_used, SpectrumSettings(), tmp_chromatograms, false);

#ifdef _OPENMP
#pragma omp critical (featureFinder)
#endif
            {
              LOG_DEBUG << "Extracted "  << tmp_chromatograms.size() << " chromatograms from SWATH map " <<
                map_idx << " with m/z " << swath_maps[map_idx].lower << " to " << swath_maps[map_idx].upper << ":" << std::endl;
              for (Size chrom_idx = 0; chrom_idx < tmp_chromatograms.size(); chrom_idx++)
              {
                // Check TIC and remove empty chromatograms (can happen if the
                // extraction window is outside the mass spectrometric acquisition
                // window).
                double tic = std::accumulate(tmp_out[chrom_idx]->getIntensityArray()->data.begin(),
                                             tmp_out[chrom_idx]->getIntensityArray()->data.end(),0.0);
                LOG_DEBUG << "Chromatogram "  << coordinates[chrom_idx].id << " with size "
                  << tmp_out[chrom_idx]->getIntensityArray()->data.size() << " and TIC " << tic  << std::endl;
                if (tic > 0.0)
                {
                  // add the chromatogram to the output
                  chromatograms.push_back(tmp_chromatograms[chrom_idx]);
                }
                else
                {
                  std::cerr << " - Warning: Empty chromatogram " << coordinates[chrom_idx].id << 
                    " detected. Will skip it!" << std::endl;
                }
              }
            }
          }
          else
          {
            LOG_DEBUG << "Extracted no transitions from SWATH map " << map_idx << " with m/z " <<
                swath_maps[map_idx].lower << " to " << swath_maps[map_idx].upper << ":" << std::endl;
          }
        }
      }
    }

    /** @brief Perform RT normalization using the MRMFeatureFinderScoring 
     *
     * @param transition_exp_ The transitions for the normalization peptides
     * @param chromatograms The extracted chromatograms
     * @param min_rsq Minimal R^2 value that is expected for the RT regression
     * @param min_coverage Minimal coverage of the chromatographic space that needs to be achieved
     * @param feature_finder_param Parameter set for the feature finding in chromatographic dimension 
     * @param irt_detection_param Parameter set for the detection of the iRTs (outlier detection, peptides per bin etc)
     * @param swath_maps The raw data for the m/z correction
     * @param mz_correction_function If correction in m/z is desired, which function should be used
     *
     * @note: feature_finder_param are copied because they are changed here.
     *
    */
    TransformationDescription RTNormalization(TargetedExperiment transition_exp_,
            std::vector< OpenMS::MSChromatogram<> > chromatograms, double min_rsq, double min_coverage,
            Param feature_finder_param, const Param& irt_detection_param, 
            std::vector< OpenSwath::SwathMap > & swath_maps, const String & mz_correction_function)
    {
      LOG_DEBUG << "Start of RTNormalization method" << std::endl;
      this->startProgress(0, 1, "Retention time normalization");

      OpenSwath::LightTargetedExperiment targeted_exp;
      OpenSwathDataAccessHelper::convertTargetedExp(transition_exp_, targeted_exp);

      bool estimateBestPeptides = irt_detection_param.getValue("estimateBestPeptides").toBool();
      if (estimateBestPeptides)
      {
        LOG_DEBUG << "Activated the 'estimateBestPeptides' option." << std::endl;
      }

      // 1. Estimate the retention time range of the whole experiment
      std::pair<double,double> RTRange = OpenSwathHelper::estimateRTRange(targeted_exp);
      LOG_DEBUG << "Detected retention time range from " << RTRange.first << " to " << RTRange.second << std::endl;

      // 2. Store the peptide retention times in an intermediate map
      std::map<OpenMS::String, double> PeptideRTMap;
      for (Size i = 0; i < targeted_exp.getCompounds().size(); i++)
      {
        PeptideRTMap[targeted_exp.getCompounds()[i].id] = targeted_exp.getCompounds()[i].rt;
      }

      // 3. Extract the RT pairs from the input data
      OpenSwath::LightTargetedExperiment transition_exp_used = targeted_exp;

      // Change the feature finding parameters:
      //  - no RT score (since we don't know the correct retention time)
      //  - no RT window
      //  - no elution model score
      //  - no peak quality (use all peaks)
      //  - if best peptides should be used, use peak quality
      MRMFeatureFinderScoring featureFinder;
      feature_finder_param.setValue("Scores:use_rt_score", "false");
      feature_finder_param.setValue("Scores:use_elution_model_score", "false");
      feature_finder_param.setValue("rt_extraction_window", -1.0);
      feature_finder_param.setValue("TransitionGroupPicker:PeakPickerMRM:signal_to_noise", 1.0); // set to 1.0 in all cases
      feature_finder_param.setValue("TransitionGroupPicker:compute_peak_quality", "false"); // no peak quality -> take all peaks!
      if (estimateBestPeptides)
      {
        feature_finder_param.setValue("TransitionGroupPicker:compute_peak_quality", "true");
        feature_finder_param.setValue("TransitionGroupPicker:minimal_quality", irt_detection_param.getValue("InitialQualityCutoff"));
      }
      featureFinder.setParameters(feature_finder_param);

      FeatureMap featureFile; // for results
      OpenMS::MRMFeatureFinderScoring::TransitionGroupMapType transition_group_map; // for results
      boost::shared_ptr<MSExperiment<Peak1D> > empty_swath_map(new MSExperiment<Peak1D>); // empty map 
      OpenSwath::SpectrumAccessPtr empty_swath_ptr = SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(empty_swath_map);
      TransformationDescription empty_trafo; // empty transformation

      // Prepare the data with the chromatograms
      boost::shared_ptr<MSExperiment<Peak1D> > xic_map(new MSExperiment<Peak1D>); 
      xic_map->setChromatograms(chromatograms);
      OpenSwath::SpectrumAccessPtr chromatogram_ptr = OpenSwath::SpectrumAccessPtr(new OpenMS::SpectrumAccessOpenMS(xic_map));

      featureFinder.setStrictFlag(false); // TODO remove this, it should be strict (e.g. all transitions need to be present for RT norm)
      featureFinder.pickExperiment(chromatogram_ptr, featureFile, transition_exp_used, empty_trafo, empty_swath_ptr, transition_group_map);

      // find most likely correct feature for each group and add it to the
      // "pairs" vector by computing pairs of iRT and real RT. 
      // Note that the quality threshold will only be applied if
      // estimateBestPeptides is true
      std::vector<std::pair<double, double> > pairs; // store the RT pairs to write the output trafoXML
      std::map<std::string, double> res = OpenSwathHelper::simpleFindBestFeature(transition_group_map, 
        estimateBestPeptides, irt_detection_param.getValue("OverallQualityCutoff"));

      for (std::map<std::string, double>::iterator it = res.begin(); it != res.end(); ++it)
      {
        pairs.push_back(std::make_pair(it->second, PeptideRTMap[it->first])); // pair<exp_rt, theor_rt>
      }

      // 4. Correct mz deviations using SwathMapMassCorrection
      SwathMapMassCorrection::correctMZ(transition_group_map, swath_maps, mz_correction_function);

      // 5. Perform the outlier detection
      std::vector<std::pair<double, double> > pairs_corrected;
      String outlier_method = irt_detection_param.getValue("outlierMethod");
      if (outlier_method == "iter_residual" || outlier_method == "iter_jackknife")
      {
        pairs_corrected = MRMRTNormalizer::removeOutliersIterative(pairs, min_rsq, min_coverage,
        irt_detection_param.getValue("useIterativeChauvenet").toBool(), outlier_method);
      }
      else if (outlier_method == "ransac")
      {
        // First, estimate of the maximum deviation from RT that is tolerated:
        //   Because 120 min gradient can have around 4 min elution shift, we use
        //   a default value of 3 % of the gradient to find upper RT threshold (3.6 min).
        double pcnt_rt_threshold = irt_detection_param.getValue("RANSACMaxPercentRTThreshold");
        double max_rt_threshold = (RTRange.second - RTRange.first) * pcnt_rt_threshold / 100.0;

        pairs_corrected = MRMRTNormalizer::removeOutliersRANSAC(pairs, min_rsq, min_coverage,
          irt_detection_param.getValue("RANSACMaxIterations"), max_rt_threshold,
          irt_detection_param.getValue("RANSACSamplingSize"));
      }
      else if (outlier_method == "none") 
      {
        pairs_corrected = pairs;
      }
      else 
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__,
          String("Illegal argument '") + outlier_method + "' used for outlierMethod (valid: 'iter_residual', 'iter_jackknife', 'ransac', 'none').");
      }

      // 5. Check whether the found peptides fulfill the binned coverage criteria
      // set by the user.
      if (estimateBestPeptides)
      {
        bool enoughPeptides = computeBinnedCoverage(RTRange, pairs_corrected,
          irt_detection_param.getValue("NrRTBins"),
          irt_detection_param.getValue("MinPeptidesPerBin"),
          irt_detection_param.getValue("MinBinsFilled") );

        if (!enoughPeptides)
        {
          throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__,
            "There were not enough bins with the minimal number of peptides");
        }
      }

      // store transformation, using a linear model as default
      TransformationDescription trafo_out;
      trafo_out.setDataPoints(pairs_corrected);
      Param model_params;
      model_params.setValue("symmetric_regression", "false");
      String model_type = "linear";
      trafo_out.fitModel(model_type, model_params);

      this->endProgress();
      return trafo_out;
    }

    /// Function to compute the coverage of the iRT peptides across the gradient
    ///   Cmp with RTNormalizer
    bool computeBinnedCoverage(const std::pair<double,double> & rtRange, 
        const std::vector<std::pair<double, double> > & pairs, int nrBins, 
        int minPeptidesPerBin, int minBinsFilled)
    {
      std::vector<int> binCounter(nrBins, 0);
      for (std::vector<std::pair<double, double> >::const_iterator pair_it = pairs.begin(); pair_it != pairs.end(); ++pair_it)
      {
        double normRT = (pair_it->second - rtRange.first) / (rtRange.second - rtRange.first); // compute a value between [0,1)
        normRT *= nrBins;
        int bin = (int)normRT;
        if (bin >= nrBins)
        {
          // this should never happen, but just to make sure
          std::cerr << "OpenSwathWorkflow::countPeptidesInBins : computed bin was too large (" << 
            bin << "), setting it to the maximum of " << nrBins << std::endl;
          bin = nrBins - 1;
        }
        binCounter[ bin ]++;
      }

      int binsFilled = 0;
      for (Size i = 0; i < binCounter.size(); i++)
      {
        LOG_DEBUG <<" In bin " << i << " out of " << binCounter.size() << 
          " we have " << binCounter[i] << " peptides " << std::endl;
        if (binCounter[i] >= minPeptidesPerBin) 
        {
          binsFilled++;
        }
      }

      return (binsFilled >= minBinsFilled);
    }

    /// Helper function to score a set of chromatograms
    /// Will iterate over all assays contained in transition_exp and for each
    /// assay fetch the corresponding chromatograms and find peakgroups.
    ///
    void scoreAllChromatograms(
        const OpenSwath::SpectrumAccessPtr input,
        const OpenSwath::SpectrumAccessPtr swath_map,
        OpenSwath::LightTargetedExperiment& transition_exp,
        const Param& feature_finder_param,
        TransformationDescription trafo, const double rt_extraction_window,
        FeatureMap& output, OpenSwathTSVWriter & tsv_writer, 
        std::map< std::string, OpenSwath::ChromatogramPtr > & ms1_chromatograms
        )
    {
      typedef OpenSwath::LightTransition TransitionType;
      // a transition group holds the MSSpectra with the Chromatogram peaks from above
      typedef MRMTransitionGroup<MSSpectrum <ChromatogramPeak>, TransitionType> MRMTransitionGroupType;
      // this is the type in which we store the chromatograms for this analysis
      typedef MSSpectrum<ChromatogramPeak> RichPeakChromatogram;

      TransformationDescription trafo_inv = trafo;
      trafo_inv.invert();

      MRMFeatureFinderScoring featureFinder;

      // To ensure multi-threading safe access to the individual spectra, we
      // need to use a light clone of the spectrum access (if multiple threads
      // share a single filestream and call seek on it, chaos will ensue).
      if (use_ms1_traces_)
      {
        OpenSwath::SpectrumAccessPtr threadsafe_ms1 = ms1_map_->lightClone();
        featureFinder.setMS1Map( threadsafe_ms1 );
      }

      MRMTransitionGroupPicker trgroup_picker;

      trgroup_picker.setParameters(feature_finder_param.copy("TransitionGroupPicker:", true));
      featureFinder.setParameters(feature_finder_param);
      featureFinder.prepareProteinPeptideMaps_(transition_exp);

      // Map chromatogram id to sequence number
      std::map<String, int> chromatogram_map;
      for (Size i = 0; i < input->getNrChromatograms(); i++)
      {
        chromatogram_map[input->getChromatogramNativeID(i)] = boost::numeric_cast<int>(i);
      }
      // Map peptide id to sequence number
      std::map<String, int> assay_peptide_map;
      for (Size i = 0; i < transition_exp.getCompounds().size(); i++)
      {
        assay_peptide_map[transition_exp.getCompounds()[i].id] = boost::numeric_cast<int>(i);
      }
      // Map peptide id to corresponding transitions
      typedef std::map<String, std::vector< const TransitionType* > > AssayMapT;
      AssayMapT assay_map;
      for (Size i = 0; i < transition_exp.getTransitions().size(); i++)
      {
        assay_map[transition_exp.getTransitions()[i].getPeptideRef()].push_back(&transition_exp.getTransitions()[i]);
      }

      std::vector<String> to_output;
      // Iterating over all the assays
      for (AssayMapT::iterator assay_it = assay_map.begin(); assay_it != assay_map.end(); ++assay_it)
      {
        // Create new MRMTransitionGroup
        String id = assay_it->first;
        MRMTransitionGroupType transition_group;
        transition_group.setTransitionGroupID(id);
        double expected_rt = transition_exp.getCompounds()[ assay_peptide_map[id] ].rt;
        double precursor_mz = -1;

        // Go through all transitions, for each transition get chromatogram and
        // the chromatogram and the assay to the MRMTransitionGroup
        for (Size i = 0; i < assay_it->second.size(); i++)
        {
          const TransitionType* transition = assay_it->second[i];

          if (chromatogram_map.find(transition->getNativeID()) == chromatogram_map.end())
          {
            throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__,
                "Error, did not find chromatogram for transitions" + transition->getNativeID() );
          }

          // Convert chromatogram to MSChromatogram
          OpenSwath::ChromatogramPtr cptr = input->getChromatogramById(chromatogram_map[transition->getNativeID()]);
          MSChromatogram<ChromatogramPeak> chromatogram_old;
          OpenSwathDataAccessHelper::convertToOpenMSChromatogram(chromatogram_old, cptr);
          RichPeakChromatogram chromatogram;

          // Extract and convert chromatogram to input chromatogram
          precursor_mz = transition->getPrecursorMZ();
          chromatogram.setMetaValue("product_mz", transition->getProductMZ());
          chromatogram.setMetaValue("precursor_mz", transition->getPrecursorMZ());
          chromatogram.setNativeID(transition->getNativeID());
          double de_normalized_experimental_rt = trafo_inv.apply(expected_rt);
          selectChrom_(chromatogram_old, chromatogram, rt_extraction_window, de_normalized_experimental_rt);

          // Now add the transition and the chromatogram to the MRMTransitionGroup
          transition_group.addTransition(*transition, transition->getNativeID());
          transition_group.addChromatogram(chromatogram, chromatogram.getNativeID());
        }

        // currently .tsv and .featureXML are mutually exclusive
        if (tsv_writer.isActive()) { output.clear(); }

        // Set the MS1 chromatogram if available
        if (!ms1_chromatograms.empty() && ms1_chromatograms.find(transition_group.getTransitionGroupID()) != ms1_chromatograms.end())
        {
          OpenSwath::ChromatogramPtr cptr = ms1_chromatograms[ transition_group.getTransitionGroupID() ];
          MSChromatogram<ChromatogramPeak> chromatogram_old;
          OpenSwathDataAccessHelper::convertToOpenMSChromatogram(chromatogram_old, cptr);
          RichPeakChromatogram chromatogram;
          selectChrom_(chromatogram_old, chromatogram, -1, -1);
          chromatogram.setMetaValue("precursor_mz", precursor_mz);
          chromatogram.setNativeID(transition_group.getTransitionGroupID() + "_" + "Precursor_i0");
          transition_group.addPrecursorChromatogram(chromatogram, "Precursor_i0");
        }

        // Process the MRMTransitionGroup: find peakgroups and score them
        trgroup_picker.pickTransitionGroup(transition_group);
        featureFinder.scorePeakgroups(transition_group, trafo, swath_map, output);

        // Add to the output tsv if given
        if (tsv_writer.isActive())
        {
          const OpenSwath::LightCompound pep = transition_exp.getCompounds()[ assay_peptide_map[id] ];
          const TransitionType* transition = assay_it->second[0];
          to_output.push_back(tsv_writer.prepareLine(pep, transition, output, id));
        }
      }

      // Only write at the very end since this is a step needs a barrier
      if (tsv_writer.isActive())
      {
#ifdef _OPENMP
#pragma omp critical (scoreAll)
#endif
        {
          tsv_writer.writeLines(to_output);
        }
      }
    }

    /// Wrapper function for prepare_coordinates that also correctly handles transformations
    void prepare_coordinates_wrap(std::vector< OpenSwath::ChromatogramPtr > & chrom_list,
      std::vector< ChromatogramExtractorAlgorithm::ExtractionCoordinates > & coordinates,
      OpenSwath::LightTargetedExperiment & transition_exp_used,
      const bool ms1, const TransformationDescription trafo_inverse,
      const ChromExtractParams & cp) const
    {
      if (cp.rt_extraction_window < 0)
      {
        prepare_coordinates(chrom_list, coordinates, transition_exp_used, cp.rt_extraction_window, ms1);
      }
      else
      {
        // Use an rt extraction window of 0.0 which will just write the retention time in start / end positions
        // Then correct the start/end positions and add the extra_rt_extract parameter
        prepare_coordinates(chrom_list, coordinates, transition_exp_used, 0.0, ms1);
        for (std::vector< ChromatogramExtractor::ExtractionCoordinates >::iterator it = coordinates.begin(); it != coordinates.end(); ++it)
        {
          it->rt_start = trafo_inverse.apply(it->rt_start) - (cp.rt_extraction_window + cp.extra_rt_extract)/ 2.0;
          it->rt_end = trafo_inverse.apply(it->rt_end) + (cp.rt_extraction_window + cp.extra_rt_extract)/ 2.0;
        }
      }
    }

    void prepare_coordinates(std::vector< OpenSwath::ChromatogramPtr > & output_chromatograms,
      std::vector< ChromatogramExtractorAlgorithm::ExtractionCoordinates > & coordinates,
      OpenSwath::LightTargetedExperiment & transition_exp_used,
      const double rt_extraction_window, const bool ms1) const
    {
      // hash of the peptide reference containing all transitions
      std::map<String, std::vector<OpenSwath::LightTransition*> > peptide_trans_map;
      for (Size i = 0; i < transition_exp_used.getTransitions().size(); i++)
      {
        peptide_trans_map[transition_exp_used.getTransitions()[i].getPeptideRef()].push_back(&transition_exp_used.getTransitions()[i]);
      }
      std::map<String, OpenSwath::LightCompound*> trans_peptide_map;
      for (Size i = 0; i < transition_exp_used.getCompounds().size(); i++)
      {
        trans_peptide_map[transition_exp_used.getCompounds()[i].id] = &transition_exp_used.getCompounds()[i];
      }

      // Determine iteration size: 
      // When extracting MS1/precursor transitions, we iterate over compounds.
      // Otherwise (for SWATH/fragment ions), we iterate over the transitions.
      Size itersize;
      if (ms1) {itersize = transition_exp_used.getCompounds().size();}
      else     {itersize = transition_exp_used.getTransitions().size();}

      for (Size i = 0; i < itersize; i++)
      {
        OpenSwath::ChromatogramPtr s(new OpenSwath::Chromatogram);
        output_chromatograms.push_back(s);

        ChromatogramExtractor::ExtractionCoordinates coord;
        OpenSwath::LightCompound pep; // TargetedExperiment::Peptide pep;
        OpenSwath::LightTransition transition;

        if (ms1)
        {
          pep = transition_exp_used.getCompounds()[i];

          // Catch cases where a peptide has no transitions
          if (peptide_trans_map.count(pep.id) == 0 )
          {
            LOG_INFO << "Warning: no transitions found for peptide " << pep.id << std::endl;
            coord.rt_start = -1;
            coord.rt_end = -2; // create negative range
            coord.id = pep.id;
            coordinates.push_back(coord);
            continue;
          }

          // This is slightly awkward but the m/z of the precursor is *not*
          // stored in the precursor object but only in the transition object
          // itself. So we have to get the first transition to look it up.
          transition = (*peptide_trans_map[pep.id][0]);
          coord.mz = transition.getPrecursorMZ();
          coord.id = pep.id;
        }
        else
        {
          transition = transition_exp_used.getTransitions()[i];
          pep = (*trans_peptide_map[transition.getPeptideRef()]);
          coord.mz = transition.getProductMZ();
          coord.id = transition.getNativeID();
        }

        double rt = pep.rt;
        coord.rt_start = rt - rt_extraction_window / 2.0;
        coord.rt_end = rt + rt_extraction_window / 2.0;
        coordinates.push_back(coord);
      }

      // sort result
      std::sort(coordinates.begin(), coordinates.end(), ChromatogramExtractor::ExtractionCoordinates::SortExtractionCoordinatesByMZ);
    }

    void selectChrom_(const MSChromatogram<ChromatogramPeak>& chromatogram_old,
      MSSpectrum<ChromatogramPeak>& chromatogram, double rt_extraction_window, double center_rt)
    {
      double rt_max = center_rt + rt_extraction_window;
      double rt_min = center_rt - rt_extraction_window;
      for (MSChromatogram<ChromatogramPeak>::const_iterator it = chromatogram_old.begin(); it != chromatogram_old.end(); ++it)
      {
        if (rt_extraction_window >= 0 && (it->getRT() < rt_min || it->getRT() > rt_max))
        {
          continue;
        }
        ChromatogramPeak peak;
        peak.setMZ(it->getRT());
        peak.setIntensity(it->getIntensity());
        chromatogram.push_back(peak);
      }
    }

  private:
    /**
     * @brief Spectrum Access to the MS1 map (note that this is *not* threadsafe!)
     *
     * @note This pointer is not threadsafe, please use the lightClone() function to create a copy for each thread
     * @note This pointer may be NULL if use_ms1_traces_ is set to false
     *
     */
    OpenSwath::SpectrumAccessPtr ms1_map_;

    /// Whether to use the MS1 traces
    bool use_ms1_traces_;

  };

}

// OpenMS base classes
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
  @page UTILS_OpenSwathWorkflow OpenSwathWorkflow

  @brief Complete workflow to run OpenSWATH

  This implements the OpenSwath workflow as described in Rost and Rosenberger
  et al. (Nature Biotechnology, 2014) and provides a complete, integrated
  analysis tool without the need to run multiple tools consecutively.

  It executes the following steps in order:

  <ul>
    <li>Reading of input files, which can be provided as one single mzML or multiple "split" mzMLs (one per SWATH)</li>
    <li>Computing the retention time transformation using RT-normalization peptides</li>
    <li>Reading of the transition list</li>
    <li>Extracting the specified transitions</li>
    <li>Scoring the peak groups in the extracted ion chromatograms (XIC)</li>
    <li>Reporting the peak groups and the chromatograms</li>
  </ul>

  See below or have a look at the INI file (via "OpenSwathWorkflow -write_ini myini.ini") for available parameters and more functionality.

  <h3>Input: SWATH maps and transition list </h3>
  SWATH maps can be provided as mzML files, either as single file directly from
  the machine (this assumes that the SWATH method has 1 MS1 and then n MS2
  spectra which are ordered the same way for each cycle). E.g. a valid method
  would be MS1, MS2 [400-425], MS2 [425-450], MS1, MS2 [400-425], MS2 [425-450]
  while an invalid method would be MS1, MS2 [400-425], MS2 [425-450], MS1, MS2
  [425-450], MS2 [400-425] where MS2 [xx-yy] indicates an MS2 scan with an
  isolation window starting at xx and ending at yy. OpenSwathWorkflow will try
  to read the SWATH windows from the data, if this is not possible please
  provide a tab-separated list with the correct windows using the
  -swath_windows_file parameter.

  Alternatively, a set of split files (n+1 mzML files) can be provided, each
  containing one SWATH map (or MS1 map).

  Since the file size can become rather large, it is recommended to not load the
  whole file into memory but rather cache it somewhere on the disk using a
  fast-access data format. This can be specified using the -readOptions cache
  parameter (this is recommended!).

  <h3>Output: Feature list and chromatograms </h3>
  The output of the OpenSwathWorkflow is a feature list, either as FeatureXML
  or as tsv (use -out_features or -out_tsv) while the latter is more memory
  friendly. If you analyze large dataset, it is recommended to only use
  -out_tsv and not -out_features. For downstream analysis (e.g. using mProphet or pyProphet)
  also the -out_tsv format is recommended.

  The feature list generated by -out_tsv is a tab-separated file. It can be
  used directly as input to the mProphet or pyProphet (a Python
  re-implementation of mProphet) software tool, see Reiter et al (2011, Nature
  Methods).

  In addition, the extracted chromatograms can be written out using the
  -out_chrom parameter.

  <h4> Feature list output format </h4>
  
  The tab-separated feature output contains the following information:

<CENTER>
  <table>
    <tr>
      <td ALIGN = "left" BGCOLOR="#EBEBEB"> Header row </td>
      <td ALIGN = "left" BGCOLOR="#EBEBEB"> Format </td>
      <td ALIGN = "left" BGCOLOR="#EBEBEB"> Description </td>
    </tr>
    <tr>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> transition_group_id </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> String </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> A unique id for the transition group (all chromatographic traces that are analyzed together)</td>
    </tr>

    <tr>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> peptide_group_label </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> String </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> A unique id for the peptide group (will be the same for each charge state and heavy/light status) </td>
    </tr>

    <tr>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> run_id </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> String </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> An identifier for the run (currently always 0)</td>
    </tr>

    <tr>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> filename </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> String </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> The input filename </td>
    </tr>

    <tr>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> RT </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Float </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Peak group retention time </td>
    </tr>

    <tr>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> id </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> String </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> A unique identifier for the peak group</td>
    </tr>

    <tr>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Sequence </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> String </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Peptide sequence (no modifications) </td>
    </tr>

    <tr>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> FullPeptideName </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> String </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Full peptide sequence including modifications in Unimod format</td>
    </tr>

    <tr>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Charge </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Int </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Assumed charge state</td>
    </tr>

    <tr>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> m/z </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Float </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Precursor m/z</td>
    </tr>

    <tr>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Intensity </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Float </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Peak group intensity (sum of all transitions)</td>
    </tr>

    <tr>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> ProteinName </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> String </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Name of the associated protein</td>
    </tr>

    <tr>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> decoy </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> String </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Whether the transition is decoy or not (0 = false, 1 = true) </td>
    </tr>

    <tr>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> assay_rt </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Float </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> The expected RT in seconds (based on normalized iRT value) </td>
    </tr>

    <tr>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> delta_rt </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Float </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> The difference between the expected RT and the peak group RT in seconds </td>
    </tr>

    <tr>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> leftWidth </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Float </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> The start of the peak group (left side) in seconds </td>
    </tr>

    <tr>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> main_var_xx_swath_prelim_score </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Float </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Initial score </td>
    </tr>

    <tr>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> norm_RT </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Float </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> The peak group retention time in normalized (iRT) space </td>
    </tr>

    <tr>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> nr_peaks </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Int </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> The number of transitions used </td>
    </tr>

    <tr>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> peak_apices_sum </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Float </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> The sum of all peak apices (may be used as alternative intensity) </td>
    </tr>

    <tr>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> potentialOutlier </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> String </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Potential outlier transitions (or "none" if none was detected)</td>
    </tr>

    <tr>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> rightWidth </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Float </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> The end of the peak group (left side) in seconds </td>
    </tr>

    <tr>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> rt_score </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Float </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> The raw RT score (unnormalized) </td>
    </tr>

    <tr>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> sn_ratio </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Float </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> The raw S/N ratio </td>
    </tr>

    <tr>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> total_xic </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Float </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> The total XIC of the chromatogram </td>
    </tr>

    <tr>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> var_... </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Float </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> One of multiple sub-scores used by OpenSWATH to describe the peak group </td>
    </tr>

    <tr>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> aggr_prec_Peak_Area </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> String </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Intensity (peak area) of MS1 traces separated by semicolon </td>
    </tr>

    <tr>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> aggr_prec_Peak_Apex </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> String </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Intensity (peak apex) of MS1 traces separated by semicolon </td>
    </tr>
    <tr>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> aggr_prec_Fragment_Annotation </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> String </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Annotation of MS1 traces separated by semicolon </td>
    </tr>
    <tr>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> aggr_Peak_Area </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> String </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Intensity (peak area) of fragment ion traces separated by semicolon </td>
    </tr>
    <tr>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> aggr_Peak_Apex </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> String </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> Intensity (peak apex) of fragment ion traces separated by semicolon </td>
    </tr>
    <tr>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> aggr_Fragment_Annotation </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=1> String </td>
      <td VALIGN="middle" ALIGN = "left" ROWSPAN=2> Annotation of fragment ion traces separated by semicolon </td>
    </tr>


  </table>
</CENTER>

  <B>The command line parameters of this tool are:</B>
  @verbinclude UTILS_OpenSwathWorkflow.cli
  <B>INI file documentation of this tool:</B>
  @htmlinclude UTILS_OpenSwathWorkflow.html

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES
class TOPPOpenSwathWorkflow
  : public TOPPBase
{
public:

  TOPPOpenSwathWorkflow()
    : TOPPBase("OpenSwathWorkflow", "Complete workflow to run OpenSWATH", false)
  {
  }

protected:

  void registerOptionsAndFlags_()
  {
    registerInputFileList_("in", "<files>", StringList(), "Input files separated by blank");
    setValidFormats_("in", ListUtils::create<String>("mzML,mzXML"));

    registerInputFile_("tr", "<file>", "", "transition file ('TraML','tsv' or 'csv')");
    setValidFormats_("tr", ListUtils::create<String>("traML,tsv,csv"));
    registerStringOption_("tr_type", "<type>", "", "input file type -- default: determined from file extension or content\n", false);
    setValidStrings_("tr_type", ListUtils::create<String>("traML,tsv,csv"));

    // one of the following two needs to be set
    registerInputFile_("tr_irt", "<file>", "", "transition file ('TraML')", false);
    setValidFormats_("tr_irt", ListUtils::create<String>("traML"));

    registerInputFile_("rt_norm", "<file>", "", "RT normalization file (how to map the RTs of this run to the ones stored in the library). If set, tr_irt may be omitted.", false, true);
    setValidFormats_("rt_norm", ListUtils::create<String>("trafoXML"));

    registerInputFile_("swath_windows_file", "<file>", "", "Optional, tab separated file containing the SWATH windows for extraction: lower_offset upper_offset \\newline 400 425 \\newline ... Note that the first line is a header and will be skipped.", false, true);
    registerFlag_("sort_swath_maps", "Sort of input SWATH files when matching to SWATH windows from swath_windows_file", true);

    registerFlag_("use_ms1_traces", "Extract the precursor ion trace(s) and use for scoring", true);
    registerFlag_("enable_uis_scoring", "Enable additional scoring of identification assays", true);

    // one of the following two needs to be set
    registerOutputFile_("out_features", "<file>", "", "output file", false);
    setValidFormats_("out_features", ListUtils::create<String>("featureXML"));

    registerStringOption_("out_tsv", "<file>", "", "TSV output file (mProphet compatible)", false);

    registerOutputFile_("out_chrom", "<file>", "", "Also output all computed chromatograms (chrom.mzML) output", false, true);
    setValidFormats_("out_chrom", ListUtils::create<String>("mzML"));

    registerDoubleOption_("min_upper_edge_dist", "<double>", 0.0, "Minimal distance to the edge to still consider a precursor, in Thomson", false, true);
    registerDoubleOption_("rt_extraction_window", "<double>", 600.0, "Only extract RT around this value (-1 means extract over the whole range, a value of 600 means to extract around +/- 300 s of the expected elution).", false);
    registerDoubleOption_("extra_rt_extraction_window", "<double>", 0.0, "Output an XIC with a RT-window that by this much larger (e.g. to visually inspect a larger area of the chromatogram)", false, true);
    registerDoubleOption_("mz_extraction_window", "<double>", 0.05, "Extraction window used (in Thomson, to use ppm see -ppm flag)", false);
    setMinFloat_("mz_extraction_window", 0.0);
    setMinFloat_("extra_rt_extraction_window", 0.0);
    registerFlag_("ppm", "m/z extraction_window is in ppm");

    registerDoubleOption_("min_rsq", "<double>", 0.95, "Minimum r-squared of RT peptides regression", false, true);
    registerDoubleOption_("min_coverage", "<double>", 0.6, "Minimum relative amount of RT peptides to keep", false, true);

    registerFlag_("split_file_input", "The input files each contain one single SWATH (alternatively: all SWATH are in separate files)", true);
    registerFlag_("use_elution_model_score", "Turn on elution model score (EMG fit to peak)", true);

    registerStringOption_("readOptions", "<name>", "normal", "Whether to run OpenSWATH directly on the input data, cache data to disk first or to perform a datareduction step first. If you choose cache, make sure to also set tempDirectory", false, true);
    setValidStrings_("readOptions", ListUtils::create<String>("normal,cache,cacheWorkingInMemory"));

    registerStringOption_("mz_correction_function", "<name>", "none", "Use the retention time normalization peptide MS2 masses to perform a mass correction (linear, weighted by intensity linear or quadratic) of all spectra.", false, true);
    setValidStrings_("mz_correction_function", ListUtils::create<String>("none,unweighted_regression,weighted_regression,quadratic_regression,weighted_quadratic_regression,weighted_quadratic_regression_delta_ppm,quadratic_regression_delta_ppm"));

    // TODO terminal slash !
    registerStringOption_("tempDirectory", "<tmp>", "/tmp/", "Temporary directory to store cached files for example", false, true);

    registerStringOption_("extraction_function", "<name>", "tophat", "Function used to extract the signal", false, true);
    setValidStrings_("extraction_function", ListUtils::create<String>("tophat,bartlett"));

    registerIntOption_("batchSize", "<number>", 0, "The batch size of chromatograms to process (0 means to only have one batch, sensible values are around 500-1000)", false, true);
    setMinInt_("batchSize", 0);

    registerSubsection_("Scoring", "Scoring parameters section");

    registerSubsection_("outlierDetection", "Parameters for the outlierDetection for iRT petides. Outlier detection can be done iteratively (by default) which removes one outlier per iteration or using the RANSAC algorithm.");
  }

  Param getSubsectionDefaults_(const String& name) const
  {
    if (name == "Scoring")
    {
      // set sensible default parameters
      Param feature_finder_param = MRMFeatureFinderScoring().getDefaults();
      feature_finder_param.remove("rt_extraction_window");
      feature_finder_param.setValue("rt_normalization_factor", 100.0); // for iRT peptides between 0 and 100 (more or less)

      feature_finder_param.setValue("TransitionGroupPicker:min_peak_width", 14.0);
      feature_finder_param.setValue("TransitionGroupPicker:recalculate_peaks", "true");
      feature_finder_param.setValue("TransitionGroupPicker:compute_peak_quality", "true");
      feature_finder_param.setValue("TransitionGroupPicker:minimal_quality", -1.5);
      feature_finder_param.setValue("TransitionGroupPicker:background_subtraction", "none");
      feature_finder_param.remove("TransitionGroupPicker:stop_after_intensity_ratio");

      // Peak Picker
      feature_finder_param.setValue("TransitionGroupPicker:PeakPickerMRM:use_gauss", "false");
      feature_finder_param.setValue("TransitionGroupPicker:PeakPickerMRM:sgolay_polynomial_order", 3);
      feature_finder_param.setValue("TransitionGroupPicker:PeakPickerMRM:sgolay_frame_length", 11);
      feature_finder_param.setValue("TransitionGroupPicker:PeakPickerMRM:peak_width", -1.0);
      feature_finder_param.setValue("TransitionGroupPicker:PeakPickerMRM:remove_overlapping_peaks", "true");
      feature_finder_param.setValue("TransitionGroupPicker:PeakPickerMRM:write_sn_log_messages", "false"); // no log messages
      // TODO it seems that the legacy method produces slightly larger peaks, e.g. it will not cut off peaks too early
      // however the same can be achieved by using a relatively low SN cutoff in the -Scoring:TransitionGroupPicker:PeakPickerMRM:signal_to_noise 0.5
      feature_finder_param.setValue("TransitionGroupPicker:recalculate_peaks_max_z", 0.75);
      feature_finder_param.setValue("TransitionGroupPicker:PeakPickerMRM:method", "corrected");
      feature_finder_param.setValue("TransitionGroupPicker:PeakPickerMRM:signal_to_noise", 0.1);
      feature_finder_param.setValue("TransitionGroupPicker:PeakPickerMRM:gauss_width", 30.0);
      feature_finder_param.setValue("uis_threshold_sn",-1);
      feature_finder_param.setValue("uis_threshold_peak_area",0);
      feature_finder_param.remove("TransitionGroupPicker:PeakPickerMRM:sn_win_len");
      feature_finder_param.remove("TransitionGroupPicker:PeakPickerMRM:sn_bin_count");
      feature_finder_param.remove("TransitionGroupPicker:PeakPickerMRM:stop_after_feature");

      // EMG Scoring - turn off by default since it is very CPU-intensive
      feature_finder_param.remove("Scores:use_elution_model_score");
      feature_finder_param.setValue("EMGScoring:max_iteration", 10);
      feature_finder_param.remove("EMGScoring:interpolation_step");
      feature_finder_param.remove("EMGScoring:tolerance_stdev_bounding_box");
      feature_finder_param.remove("EMGScoring:deltaAbsError");

      // remove these parameters
      feature_finder_param.remove("add_up_spectra");
      feature_finder_param.remove("spacing_for_spectra_resampling");
      feature_finder_param.remove("EMGScoring:statistics:mean");
      feature_finder_param.remove("EMGScoring:statistics:variance");
      return feature_finder_param;
    }
    else if (name == "outlierDetection")
    {
      Param p;
      p.setValue("outlierMethod", "iter_residual", "Which outlier detection method to use (valid: 'iter_residual', 'iter_jackknife', 'ransac', 'none'). Iterative methods remove one outlier at a time. Jackknife approach optimizes for maximum r-squared improvement while 'iter_residual' removes the datapoint with the largest residual error (removal by residual is computationally cheaper, use this with lots of peptides).");
      p.setValidStrings("outlierMethod", ListUtils::create<String>("iter_residual,iter_jackknife,ransac,none"));

      p.setValue("useIterativeChauvenet", "false", "Whether to use Chauvenet's criterion when using iterative methods. This should be used if the algorithm removes too many datapoints but it may lead to true outliers being retained.");
      p.setValidStrings("useIterativeChauvenet", ListUtils::create<String>("true,false"));

      p.setValue("RANSACMaxIterations", 1000, "Maximum iterations for the RANSAC outlier detection algorithm.");
      p.setValue("RANSACMaxPercentRTThreshold", 3, "Maximum threshold in RT dimension for the RANSAC outlier detection algorithm (in percent of the total gradient). Default is set to 3% which is around +/- 4 minutes on a 120 gradient.");
      p.setValue("RANSACSamplingSize", 10, "Sampling size of data points per iteration for the RANSAC outlier detection algorithm.");

      p.setValue("estimateBestPeptides", "false", "Whether the algorithms should try to choose the best peptides based on their peak shape for normalization. Use this option you do not expect all your peptides to be detected in a sample and too many 'bad' peptides enter the outlier removal step (e.g. due to them being endogenous peptides or using a less curated list of peptides).");
      p.setValidStrings("estimateBestPeptides", ListUtils::create<String>("true,false"));

      p.setValue("InitialQualityCutoff", 0.5, "The initial overall quality cutoff for a peak to be scored (range ca. -2 to 2)");
      p.setValue("OverallQualityCutoff", 5.5, "The overall quality cutoff for a peak to go into the retention time estimation (range ca. 0 to 10)");
      p.setValue("NrRTBins", 10, "Number of RT bins to use to compute coverage. This option should be used to ensure that there is a complete coverage of the RT space (this should detect cases where only a part of the RT gradient is actually covered by normalization peptides)");
      p.setValue("MinPeptidesPerBin", 1, "Minimal number of peptides that are required for a bin to counted as 'covered'");
      p.setValue("MinBinsFilled", 8, "Minimal number of bins required to be covered");
      return p;
    }
    else
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Unknown subsection", name);
    }
  }

  void loadSwathFiles(StringList& file_list, bool split_file, String tmp, String readoptions,
    boost::shared_ptr<ExperimentalSettings > & exp_meta,
    std::vector< OpenSwath::SwathMap > & swath_maps)
  {
    SwathFile swath_file;
    swath_file.setLogType(log_type_);

    if (split_file || file_list.size() > 1)
    {
      // TODO cannot use data reduction here any more ...
      swath_maps = swath_file.loadSplit(file_list, tmp, exp_meta, readoptions);
    }
    else
    {
      FileTypes::Type in_file_type = FileTypes::nameToType(file_list[0]);
      if (in_file_type == FileTypes::MZML || file_list[0].suffix(4).toLower() == "mzml"
        || file_list[0].suffix(7).toLower() == "mzml.gz"  )
      {
        swath_maps = swath_file.loadMzML(file_list[0], tmp, exp_meta, readoptions);
      }
      else if (in_file_type == FileTypes::MZXML || file_list[0].suffix(5).toLower() == "mzxml"
        || file_list[0].suffix(8).toLower() == "mzxml.gz"  )
      {
        swath_maps = swath_file.loadMzXML(file_list[0], tmp, exp_meta, readoptions);
      }
      else
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__,
            "Input file needs to have ending mzML or mzXML");
      }
    }
  }

  /**
   * @brief Load the retention time transformation file
   *
   * This function will create the retention time transformation either by
   * loading a provided .trafoXML file or determine it from the data itself by
   * extracting the transitions specified in the irt_tr_file TraML file.
   *
   * @param trafo_in Input trafoXML file (if not empty, transformation will be
   *                 loaded from this file)
   * @param irt_tr_file  Input TraML file containing transitions (if trafo_in
   *                     is empty, this file will be loaded and transitions
   *                     will be extracted)
   * @param swath_maps The raw data (swath maps)
   * @param min_rsq Minimal R^2 value that is expected for the RT regression
   * @param min_coverage Minimal coverage of the chromatographic space that needs to be achieved
   * @param feature_finder_param Parameter set for the feature finding in chromatographic dimension 
   * @param cp_irt Parameter set for the chromatogram extraction
   * @param irt_detection_param Parameter set for the detection of the iRTs (outlier detection, peptides per bin etc)
   * @param mz_correction_function If correction in m/z is desired, which function should be used
   * @param debug_level Debug level (writes out the RT normalization chromatograms if larger than 1)
   *
   */
  TransformationDescription loadTrafoFile(String trafo_in, String irt_tr_file,
    std::vector< OpenSwath::SwathMap > & swath_maps, double min_rsq, double min_coverage,
    const Param& feature_finder_param, const OpenSwathWorkflow::ChromExtractParams& cp_irt,
    const Param& irt_detection_param, const String & mz_correction_function, Size debug_level)
  {
    TransformationDescription trafo_rtnorm;
    if (!trafo_in.empty())
    {
      // get read RT normalization file
      TransformationXMLFile trafoxml;
      trafoxml.load(trafo_in, trafo_rtnorm, false);
      Param model_params = getParam_().copy("model:", true);
      model_params.setValue("symmetric_regression", "false");
      String model_type = "linear";
      trafo_rtnorm.fitModel(model_type, model_params);
    }
    else if (!irt_tr_file.empty())
    {
      OpenSwathWorkflow wf(false);
      wf.setLogType(log_type_);
      // Loading iRT file
      std::cout << "Will load iRT transitions and try to find iRT peptides" << std::endl;
      TraMLFile traml;
      OpenMS::TargetedExperiment irt_transitions;
      traml.load(irt_tr_file, irt_transitions);
      trafo_rtnorm = wf.performRTNormalization(irt_transitions, swath_maps, min_rsq, min_coverage,
          feature_finder_param, cp_irt, irt_detection_param, mz_correction_function, debug_level);
    }
    return trafo_rtnorm;
  }

  ExitCodes main_(int, const char **)
  {
    ///////////////////////////////////
    // Prepare Parameters
    ///////////////////////////////////
    StringList file_list = getStringList_("in");
    String tr_file = getStringOption_("tr");

    Param irt_detection_param = getParam_().copy("outlierDetection:", true);

    //tr_file input file type
    FileHandler fh_tr_type;
    FileTypes::Type tr_type = FileTypes::nameToType(getStringOption_("tr_type"));

    if (tr_type == FileTypes::UNKNOWN)
    {
      tr_type = fh_tr_type.getType(tr_file);
      writeDebug_(String("Input file type: ") + FileTypes::typeToName(tr_type), 2);
    }

    if (tr_type == FileTypes::UNKNOWN)
    {
      writeLog_("Error: Could not determine input file type!");
      return PARSE_ERROR;
    }

    String out = getStringOption_("out_features");
    String out_tsv = getStringOption_("out_tsv");

    String irt_tr_file = getStringOption_("tr_irt");
    String trafo_in = getStringOption_("rt_norm");

    String out_chrom = getStringOption_("out_chrom");
    bool ppm = getFlag_("ppm");
    bool split_file = getFlag_("split_file_input");
    bool use_emg_score = getFlag_("use_elution_model_score");
    bool force = getFlag_("force");
    bool sort_swath_maps = getFlag_("sort_swath_maps");
    bool use_ms1_traces = getFlag_("use_ms1_traces");
    bool enable_uis_scoring = getFlag_("enable_uis_scoring");
    double min_upper_edge_dist = getDoubleOption_("min_upper_edge_dist");
    double mz_extraction_window = getDoubleOption_("mz_extraction_window");
    double rt_extraction_window = getDoubleOption_("rt_extraction_window");
    double extra_rt_extract = getDoubleOption_("extra_rt_extraction_window");
    String extraction_function = getStringOption_("extraction_function");
    String swath_windows_file = getStringOption_("swath_windows_file");
    int batchSize = (int)getIntOption_("batchSize");
    Size debug_level = (Size)getIntOption_("debug");

    double min_rsq = getDoubleOption_("min_rsq");
    double min_coverage = getDoubleOption_("min_coverage");

    String readoptions = getStringOption_("readOptions");
    String mz_correction_function = getStringOption_("mz_correction_function");
    String tmp = getStringOption_("tempDirectory");

    ///////////////////////////////////
    // Parameter validation
    ///////////////////////////////////

    bool load_into_memory = false;
    if (readoptions == "cacheWorkingInMemory")
    {
      readoptions = "cache";
      load_into_memory = true;
    }

    if (trafo_in.empty() && irt_tr_file.empty())
    {
      std::cout << "Since neither rt_norm nor tr_irt is set, OpenSWATH will " <<
        "not use RT-transformation (rather a null transformation will be applied)" << std::endl;
    }
    if ( (out.empty() && out_tsv.empty()) || (!out.empty() && !out_tsv.empty()) )
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__,
          "Either out_features or out_tsv needs to be set (but not both)");
    }

    // Check swath window input
    if (!swath_windows_file.empty())
    {
      LOG_INFO << "Validate provided Swath windows file:" << std::endl;
      std::vector<double> swath_prec_lower;
      std::vector<double> swath_prec_upper;
      SwathWindowLoader::readSwathWindows(swath_windows_file, swath_prec_lower, swath_prec_upper);

      LOG_INFO << "Read Swath maps file with " << swath_prec_lower.size() << " windows." << std::endl;
      for (Size i = 0; i < swath_prec_lower.size(); i++)
      {
        LOG_DEBUG << "Read lower swath window " << swath_prec_lower[i] << " and upper window " << swath_prec_upper[i] << std::endl;
      }
    }

    OpenSwathWorkflow::ChromExtractParams cp;
    cp.min_upper_edge_dist   = min_upper_edge_dist;
    cp.mz_extraction_window  = mz_extraction_window;
    cp.ppm                   = ppm;
    cp.rt_extraction_window  = rt_extraction_window,
    cp.extraction_function   = extraction_function;
    cp.extra_rt_extract      = extra_rt_extract;

    OpenSwathWorkflow::ChromExtractParams cp_irt = cp;
    cp_irt.rt_extraction_window = -1; // extract the whole RT range

    Param feature_finder_param = getParam_().copy("Scoring:", true);
    if (use_emg_score)
    {
      feature_finder_param.setValue("Scores:use_elution_model_score", "true");
    }
    else
    {
      feature_finder_param.setValue("Scores:use_elution_model_score", "false");
    }
    if (use_ms1_traces)
    {
      feature_finder_param.setValue("Scores:use_ms1_correlation", "true");
      feature_finder_param.setValue("Scores:use_ms1_fullscan", "true");
    }
    if (enable_uis_scoring)
    {
      feature_finder_param.setValue("Scores:use_uis_scores", "true");
    }


    ///////////////////////////////////
    // Load the SWATH files
    ///////////////////////////////////

    // (i) Load files
    boost::shared_ptr<ExperimentalSettings> exp_meta(new ExperimentalSettings);
    std::vector< OpenSwath::SwathMap > swath_maps;
    loadSwathFiles(file_list, split_file, tmp, readoptions, exp_meta, swath_maps);

    // (ii) Allow the user to specify the SWATH windows
    if (!swath_windows_file.empty())
    {
      SwathWindowLoader::annotateSwathMapsFromFile(swath_windows_file, swath_maps, sort_swath_maps);
    }

    for (Size i = 0; i < swath_maps.size(); i++)
    {
      LOG_DEBUG << "Found swath map " << i << " with lower " << swath_maps[i].lower
        << " and upper " << swath_maps[i].upper << " and " << swath_maps[i].sptr->getNrSpectra()
        << " spectra." << std::endl;
    }

    // (iii) Sanity check: there should be no overlap between the windows:
    std::vector<std::pair<double, double> > sw_windows;
    for (Size i = 0; i < swath_maps.size(); i++)
    {
      if (!swath_maps[i].ms1)
      {
        sw_windows.push_back(std::make_pair(swath_maps[i].lower, swath_maps[i].upper));
      }
    }
    std::sort(sw_windows.begin(), sw_windows.end(), SortPairDoubleByFirst);

    for (Size i = 1; i < sw_windows.size(); i++)
    {
      double lower_map_end = sw_windows[i-1].second - min_upper_edge_dist;
      double upper_map_start = sw_windows[i].first;
      LOG_DEBUG << "Extraction will go up to " << lower_map_end << " and continue at " << upper_map_start << std::endl;

      if (upper_map_start - lower_map_end > 0.01)
      {
        LOG_WARN << "Extraction will have a gap between " << lower_map_end << " and " << upper_map_start << std::endl;
        if (!force)
        {
          LOG_ERROR << "Extraction windows have a gap. Will abort (override with -force)" << std::endl;
          return PARSE_ERROR;
        }
      }

      if (lower_map_end - upper_map_start > 0.01)
      {
        LOG_WARN << "Extraction will overlap between " << lower_map_end << " and " << upper_map_start << std::endl;
        LOG_WARN << "This will lead to multiple extraction of the transitions in the overlapping region" << 
                    "which will lead to duplicated output. It is very unlikely that you want this." << std::endl;
        LOG_WARN << "Please fix this by providing an appropriate extraction file with -swath_windows_file" << std::endl;
        if (!force)
        {
          LOG_ERROR << "Extraction windows overlap. Will abort (override with -force)" << std::endl;
          return PARSE_ERROR;
        }
      }

    }

    ///////////////////////////////////
    // Get the transformation information (using iRT peptides)
    ///////////////////////////////////
    TransformationDescription trafo_rtnorm = loadTrafoFile(trafo_in, irt_tr_file,
        swath_maps, min_rsq, min_coverage, feature_finder_param, cp_irt, irt_detection_param, mz_correction_function, debug_level);

    ///////////////////////////////////
    // Load the transitions
    ///////////////////////////////////
    OpenSwath::LightTargetedExperiment transition_exp;
    ProgressLogger progresslogger;
    std::cout << " do laod traml files " << std::endl;
    progresslogger.setLogType(log_type_);
    progresslogger.startProgress(0, swath_maps.size(), "Load TraML file");
    FileTypes::Type tr_file_type = FileTypes::nameToType(tr_file);
    if (tr_file_type == FileTypes::TRAML || tr_file.suffix(5).toLower() == "traml"  )
    {
      TargetedExperiment targeted_exp;
      TraMLFile().load(tr_file, targeted_exp);
      OpenSwathDataAccessHelper::convertTargetedExp(targeted_exp, transition_exp);
    }
    else
    {
      TransitionTSVReader().convertTSVToTargetedExperiment(tr_file.c_str(), tr_type, transition_exp);
    }
    progresslogger.endProgress();

    ///////////////////////////////////
    // Set up chrom.mzML output
    ///////////////////////////////////
    MSDataWritingConsumer * chromConsumer;
    if (!out_chrom.empty())
    {
      chromConsumer = new PlainMSDataWritingConsumer(out_chrom);
      int expected_chromatograms = transition_exp.transitions.size();
      chromConsumer->setExpectedSize(0, expected_chromatograms);
      chromConsumer->setExperimentalSettings(*exp_meta);
      chromConsumer->getOptions().setWriteIndex(true);  // ensure that we write the index
      chromConsumer->addDataProcessing(getProcessingInfo_(DataProcessing::SMOOTHING));
    }
    else
    {
      chromConsumer = new NoopMSDataWritingConsumer("");
    }

    ///////////////////////////////////
    // Extract and score
    ///////////////////////////////////
    FeatureMap out_featureFile;

    OpenSwathTSVWriter tsvwriter(out_tsv, file_list[0], use_ms1_traces, enable_uis_scoring);
    OpenSwathWorkflow wf(use_ms1_traces);
    wf.setLogType(log_type_);

    wf.performExtraction(swath_maps, trafo_rtnorm, cp, feature_finder_param, transition_exp,
        out_featureFile, !out.empty(), tsvwriter, chromConsumer, batchSize, load_into_memory);
    if (!out.empty())
    {
      addDataProcessing_(out_featureFile, getProcessingInfo_(DataProcessing::QUANTITATION));
      out_featureFile.ensureUniqueId();
      FeatureXMLFile().store(out, out_featureFile);
    }

    delete chromConsumer;

    return EXECUTION_OK;
  }

};

int main(int argc, const char ** argv)
{
  TOPPOpenSwathWorkflow tool;
  return tool.main(argc, argv);
}

/// @endcond
