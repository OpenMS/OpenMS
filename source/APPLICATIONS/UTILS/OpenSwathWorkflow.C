// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
#include <OpenMS/FORMAT/TraMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/TransformationXMLFile.h>
#include <OpenMS/ANALYSIS/OPENSWATH/TransitionTSVReader.h>
#include <OpenMS/ANALYSIS/OPENSWATH/CachedmzML.h>
#ifdef OPENMS_FORMAT_SWATHFILE_MZXMLSUPPORT
#include "MSDataReader.h"
#endif
#include <OpenMS/FORMAT/SwathFile.h>

// Kernel and implementations
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessOpenMS.h>

// Helpers
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathHelper.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/DataAccessHelper.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SimpleOpenMSSpectraAccessFactory.h>

// Algorithms
#include <OpenMS/ANALYSIS/OPENSWATH/MRMRTNormalizer.h>
#include <OpenMS/ANALYSIS/OPENSWATH/ChromatogramExtractor.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMFeatureFinderScoring.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMTransitionGroupPicker.h>

#include <assert.h>

using namespace OpenMS;

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

  public:

    OpenSwathTSVWriter(String output_filename, String input_filename = "inputfile") :
      ofs(output_filename.c_str()), 
      input_filename_(input_filename),
      doWrite_(!output_filename.empty())
      {}

    bool isActive() {return doWrite_;}

    void writeHeader()
    {
      ofs << "transition_group_id\trun_id\tfilename\tRT\tid\tSequence\tFullPeptideName\tCharge\tm/z\tIntensity\tProteinName\tdecoy\tassay_rt\tdelta_rt\tleftWidth\tmain_var_xx_swath_prelim_score\tnorm_RT\tnr_peaks\tpeak_apices_sum\tpotentialOutlier\trightWidth\trt_score\tsn_ratio\ttotal_xic\tvar_bseries_score\tvar_dotprod_score\tvar_intensity_score\tvar_isotope_correlation_score\tvar_isotope_overlap_score\tvar_library_corr\tvar_library_dotprod\tvar_library_manhattan\tvar_library_rmsd\tvar_library_rootmeansquare\tvar_library_sangle\tvar_log_sn_score\tvar_manhatt_score\tvar_massdev_score\tvar_massdev_score_weighted\tvat_norm_rt_score\tvar_xcorr_coelution\tvar_xcorr_coelution_weighted\tvar_xcorr_shape\tvar_xcorr_shape_weighted\tvar_yseries_score\tvar_elution_model_fit_score\txx_lda_prelim_score\txx_swath_prelim_score\taggr_Peak_Area\taggr_Peak_Apex\taggr_Fragment_Annotation\n";
    }

    String prepareLine(const OpenSwath::LightPeptide & pep,
        const OpenSwath::LightTransition* transition,
        FeatureMap<>& output, String id)
    {
        String result = "";
        String decoy = "0"; // 0 = false
        if (transition->decoy) decoy = "1";
        for (FeatureMap<>::iterator feature_it = output.begin(); feature_it != output.end(); feature_it++)
        {

          char intensity_char[40];
          String aggr_Peak_Area = "";
          String aggr_Peak_Apex = "";
          String aggr_Fragment_Annotation = "";
          for (std::vector<Feature>::iterator sub_it = feature_it->getSubordinates().begin(); sub_it != feature_it->getSubordinates().end(); ++sub_it)
          {
            sprintf(intensity_char, "%f", sub_it->getIntensity());
            aggr_Peak_Area += (String)intensity_char + ";";
            aggr_Peak_Apex +=  "NA;";
            aggr_Fragment_Annotation += (String)sub_it->getMetaValue("native_id") + ";";
          }
          if (!feature_it->getSubordinates().empty())
          {
            aggr_Peak_Area = aggr_Peak_Area.substr(0, aggr_Peak_Area.size() - 1);
            aggr_Peak_Apex = aggr_Peak_Apex.substr(0, aggr_Peak_Apex.size() - 1);
            aggr_Fragment_Annotation = aggr_Fragment_Annotation.substr(0, aggr_Fragment_Annotation.size() - 1);
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

          String line = "";
          line += id + "_run0"
            + "\t" + "0" 
            + "\t" + input_filename_
            + "\t" + (String)feature_it->getRT() 
            + "\t" + "f_" + feature_it->getUniqueId()  // TODO might not be unique!!! 
            + "\t" + pep.sequence
            + "\t" + full_peptide_name
            + "\t" + (String)pep.charge
            + "\t" + (String)transition->precursor_mz
            + "\t" + (String)feature_it->getIntensity() 
            + "\t" + pep.protein_ref
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
            + "\t" + (String)feature_it->getMetaValue("var_elution_model_fit_score") 
            + "\t" + (String)feature_it->getMetaValue("xx_lda_prelim_score") 
            + "\t" + (String)feature_it->getMetaValue("xx_swath_prelim_score") 
            + "\t" + aggr_Peak_Area + "\t" + aggr_Peak_Apex + "\t" + aggr_Fragment_Annotation + "\n";
          result += line;
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

    /** @brief ChromatogramExtractor parameters
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
    */
    TransformationDescription performRTNormalization(const OpenMS::TargetedExperiment & irt_transitions, 
            const std::vector< OpenSwath::SwathMap > & swath_maps, double min_rsq, double min_coverage, 
            const Param& feature_finder_param, const ChromExtractParams cp_irt)
    {
      LOG_DEBUG << "performRTNormalization method starting" << std::endl;
      std::vector< OpenMS::MSChromatogram<> > irt_chromatograms;
      simpleExtractChromatograms(swath_maps, irt_transitions, irt_chromatograms, cp_irt);
      LOG_DEBUG << "Extracted number of chromatograms from iRT files: " << irt_chromatograms.size() <<  std::endl;
      // get RT normalization from data
      return RTNormalization(irt_transitions,
              irt_chromatograms, min_rsq, min_coverage, feature_finder_param);
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
    */
    void performExtraction(const std::vector< OpenSwath::SwathMap > & swath_maps,
      const TransformationDescription trafo,
      const ChromExtractParams cp, const Param& feature_finder_param,
      const OpenSwath::LightTargetedExperiment& transition_exp, 
      FeatureMap<>& out_featureFile, String out, 
      OpenSwathTSVWriter & tsv_writer, Interfaces::IMSDataConsumer<> * chromConsumer, 
      int batchSize)
    {
      tsv_writer.writeHeader();

      TransformationDescription trafo_inverse = trafo;
      trafo_inverse.invert();

      std::cout << "Will analyze " << transition_exp.transitions.size() << " transitions in total." << std::endl;
      int progress = 0;
      this->startProgress(0, swath_maps.size(), "Extracting and scoring transitions");
      
      // We set dynamic scheduling such that the maps are worked on in the order
      // in which they were given to the program / acquired. This gives much
      // better load balancing than static allocation.
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic,1)
#endif
      for (SignedSize i = 0; i < boost::numeric_cast<SignedSize>(swath_maps.size()); ++i)
      {
        if (!swath_maps[i].ms1) { // continue if MS1


        // Step 1: select transitions
        OpenSwath::LightTargetedExperiment transition_exp_used_all;
        OpenSwathHelper::selectSwathTransitions(transition_exp, transition_exp_used_all,
            cp.min_upper_edge_dist, swath_maps[i].lower, swath_maps[i].upper);
        if (transition_exp_used_all.getTransitions().size() > 0) { // continue if no transitions found

        int batch_size;
        if (batchSize <= 0 || batchSize >= (int)transition_exp_used_all.getPeptides().size()) 
        {
          batch_size = transition_exp_used_all.getPeptides().size();
        }
        else {batch_size = batchSize;}
#ifdef _OPENMP
#pragma omp critical (featureFinder)
#endif
        { std::cout << "Thread " << 
#ifdef _OPENMP
          omp_get_thread_num() << " " <<
#endif
          "will analyze " << transition_exp_used_all.getPeptides().size() <<  " peptides and "
          << transition_exp_used_all.getTransitions().size() <<  " transitions "
          "from SWATH " << i << " in batches of " << batch_size << std::endl; }
        for (size_t j = 0; j <= (transition_exp_used_all.getPeptides().size() / batch_size) ; j++)
        {
          // Create the new, batch-size transition experiment
          OpenSwath::LightTargetedExperiment transition_exp_used;
          selectPeptidesForBatch_(transition_exp_used_all, transition_exp_used, batch_size, j);

          // Step 2: extract these transitions
          ChromatogramExtractor extractor;
          boost::shared_ptr<MSExperiment<Peak1D> > chrom_exp(new MSExperiment<Peak1D>);

          std::vector< OpenSwath::ChromatogramPtr > chrom_list;
          std::vector< ChromatogramExtractor::ExtractionCoordinates > coordinates;

          // Step 2.1: prepare the extraction coordinates
          if (cp.rt_extraction_window < 0)
          {
            prepare_coordinates(chrom_list, coordinates, transition_exp_used, cp.rt_extraction_window, false);
          }
          else
          {
            // Use an rt extraction window of 0.0 which will just write the retention time in start / end positions
            // Then correct the start/end positions and add the extra_rt_extract parameter
            prepare_coordinates(chrom_list, coordinates, transition_exp_used, 0.0, false);
            for (std::vector< ChromatogramExtractor::ExtractionCoordinates >::iterator it = coordinates.begin(); it != coordinates.end(); it++)
            {
              it->rt_start = trafo_inverse.apply(it->rt_start) - (cp.rt_extraction_window + cp.extra_rt_extract)/ 2.0;
              it->rt_end = trafo_inverse.apply(it->rt_end) + (cp.rt_extraction_window + cp.extra_rt_extract)/ 2.0;
            }
          }

          // Step 2.2: extract chromatograms
          extractor.extractChromatograms(swath_maps[i].sptr, chrom_list, coordinates, cp.mz_extraction_window,
              cp.ppm, cp.extraction_function);

          // Step 2.3: convert chromatograms back and write to output
          std::vector< OpenMS::MSChromatogram<> > chromatograms;
          extractor.return_chromatogram(chrom_list, coordinates, transition_exp_used,  SpectrumSettings(), chromatograms, false);
          chrom_exp->setChromatograms(chromatograms);
          OpenSwath::SpectrumAccessPtr chromatogram_ptr = OpenSwath::SpectrumAccessPtr(new OpenMS::SpectrumAccessOpenMS(chrom_exp));

          // Step 3: score these extracted transitions
          FeatureMap<> featureFile;
          scoreAllChromatograms(chromatogram_ptr, swath_maps[i].sptr, transition_exp_used, 
              feature_finder_param, trafo, cp.rt_extraction_window, featureFile, tsv_writer);

          // Step 4: write all chromatograms and features out into an output object / file 
          // (this needs to be done in a critical section since we only have one
          // output file and one output map).
#ifdef _OPENMP
#pragma omp critical (featureFinder)
#endif
          {
            // write chromatograms to output if so desired
            for (Size j = 0; j < chromatograms.size(); j++)
            {
              chromConsumer->consumeChromatogram(chromatograms[j]); 
            }
            // write features to output if so desired
            if (!out.empty())
            {
              for (FeatureMap<Feature>::iterator feature_it = featureFile.begin();
                   feature_it != featureFile.end(); feature_it++)
              {
                out_featureFile.push_back(*feature_it);
              }
              for (std::vector<ProteinIdentification>::iterator protid_it =
                     featureFile.getProteinIdentifications().begin();
                   protid_it != featureFile.getProteinIdentifications().end();
                   protid_it++)
              {
                out_featureFile.getProteinIdentifications().push_back(*protid_it);
              }
              this->setProgress(progress++);
            }
          }
        }

      } // continue 2
      } // continue 1
      }
      this->endProgress();
    }

  private:

    /** @brief Select which peptides to analyze in the next batch and copy the corresponding peptides and transitions to transition_exp_used
     *
     * @param transition_exp_used input (all transitions for this swath)
     * @param transition_exp_used output (contains only transitions for the next batch)
     * @param batch_size how many peptides per batch
     * @param j batch number (peptides from j*batch_size to j*batch_size+batch_size will be copied)
     *
    */
    void selectPeptidesForBatch_(const OpenSwath::LightTargetedExperiment& transition_exp_used_all, 
      OpenSwath::LightTargetedExperiment& transition_exp_used, int batch_size, size_t j)
    {
      // compute batch start/end
      size_t start = j*batch_size;
      size_t end = j*batch_size+batch_size;
      if (end > transition_exp_used_all.peptides.size() ) {end = transition_exp_used_all.peptides.size();}

      // Create the new, batch-size transition experiment
      transition_exp_used.proteins = transition_exp_used_all.proteins;
      transition_exp_used.peptides.insert(transition_exp_used.peptides.end(), 
          transition_exp_used_all.peptides.begin() + start, transition_exp_used_all.peptides.begin() + end);
      copyBatchTransitions_(transition_exp_used.peptides, transition_exp_used_all.transitions, transition_exp_used.transitions);
    }

    /// Copy the required transitions from all_transitions to output (e.g. those that match a peptide in the used_peptides vector)
    void copyBatchTransitions_(const std::vector<OpenSwath::LightPeptide>& used_peptides, 
        const std::vector<OpenSwath::LightTransition>& all_transitions, std::vector<OpenSwath::LightTransition>& output)
    {
      std::set<std::string> selected_peptides;
      for (Size i = 0; i < used_peptides.size(); i++)
      {
        selected_peptides.insert(used_peptides[i].id);
      }

      for (Size i = 0; i < all_transitions.size(); i++)
      {
        if (selected_peptides.find(all_transitions[i].peptide_ref) != selected_peptides.end())
        {
          output.push_back(all_transitions[i]);
        }
      }
    }

    /// Simple method to extract chromatograms (for the RT-normalization peptides)
    void simpleExtractChromatograms(const std::vector< OpenSwath::SwathMap > & swath_maps,
      const OpenMS::TargetedExperiment & irt_transitions, 
      std::vector< OpenMS::MSChromatogram<> > & chromatograms, ChromExtractParams cp)
    {
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for (SignedSize i = 0; i < boost::numeric_cast<SignedSize>(swath_maps.size()); ++i)
      {
        if (!swath_maps[i].ms1) { // continue if MS1
        TargetedExperiment transition_exp_used;
        OpenSwathHelper::selectSwathTransitions(irt_transitions, transition_exp_used,
            cp.min_upper_edge_dist, swath_maps[i].lower, swath_maps[i].upper);
        if (transition_exp_used.getTransitions().size() > 0) { // continue if no transitions found

        std::vector< OpenSwath::ChromatogramPtr > tmp_out;
        std::vector< ChromatogramExtractor::ExtractionCoordinates > coordinates;
        ChromatogramExtractor extractor;
        // TODO for lrage rt extraction windows!
        extractor.prepare_coordinates(tmp_out, coordinates, transition_exp_used,  cp.rt_extraction_window, false);
        extractor.extractChromatograms(swath_maps[i].sptr, tmp_out, coordinates, cp.mz_extraction_window,
            cp.ppm, cp.extraction_function);

#ifdef _OPENMP
#pragma omp critical (featureFinder)
#endif
        {
          LOG_DEBUG << "Extracted "  << tmp_out.size() << " chromatograms from SWATH map " << 
              i << " with m/z " << swath_maps[i].lower << " to " << swath_maps[i].upper << ":" << std::endl;
          for (Size i = 0; i < tmp_out.size(); i++)
          { 
            // Check TIC and remove empty chromatograms (can happen if the
            // extraction window is outside the mass spectrometric acquisition
            // window).
            double tic = std::accumulate(tmp_out[i]->getIntensityArray()->data.begin(),tmp_out[i]->getIntensityArray()->data.end(),0);
            LOG_DEBUG << "Chromatogram "  << coordinates[i].id << " with size " 
                << tmp_out[i]->getIntensityArray()->size() << " and TIC " << tic  << std::endl;
            if (tic <= 0.0) 
            {
              std::cerr  << " - Warning: Empty chromatogram " << coordinates[i].id << " detected. Will skip it!" << std::endl;
              continue;
            }

            OpenMS::MSChromatogram<> chrom;
            OpenSwathDataAccessHelper::convertToOpenMSChromatogram(chrom, tmp_out[i]);
            chrom.setNativeID(coordinates[i].id);
            chromatograms.push_back(chrom);
          }
        }
      } // continue 2
      else
      {
        LOG_DEBUG << "Extracted no transitions from SWATH map " << i << " with m/z " << 
            swath_maps[i].lower << " to " << swath_maps[i].upper << ":" << std::endl;
      }
      } // continue 1
      }
    }

    /// @note: feature_finder_param are copied because they are changed here.
    TransformationDescription RTNormalization(TargetedExperiment transition_exp_,
            std::vector< OpenMS::MSChromatogram<> > chromatograms, double min_rsq, double min_coverage, 
            Param feature_finder_param)
    {
      LOG_DEBUG << "Start of RTNormalization method" << std::endl;
      this->startProgress(0, 1, "Retention time normalization");

      OpenSwath::LightTargetedExperiment targeted_exp;
      OpenSwathDataAccessHelper::convertTargetedExp(transition_exp_, targeted_exp);

      std::vector<std::pair<double, double> > pairs; // store the RT pairs to write the output trafoXML

      // Store the peptide retention times in an intermediate map
      std::map<OpenMS::String, double> PeptideRTMap;
      for (Size i = 0; i < targeted_exp.getPeptides().size(); i++)
      {
        PeptideRTMap[targeted_exp.getPeptides()[i].id] = targeted_exp.getPeptides()[i].rt; 
      }

      OpenSwath::LightTargetedExperiment transition_exp_used = targeted_exp;

      MRMFeatureFinderScoring featureFinder;
      feature_finder_param.setValue("Scores:use_rt_score", "false");
      feature_finder_param.setValue("Scores:use_elution_model_score", "false");
      feature_finder_param.setValue("rt_extraction_window", -1.0);
      feature_finder_param.setValue("TransitionGroupPicker:PeakPickerMRM:signal_to_noise", 1.0); // set to 1.0 in all cases
      feature_finder_param.setValue("TransitionGroupPicker:compute_peak_quality", "false"); // no peak quality -> take all peaks!

      featureFinder.setParameters(feature_finder_param);
      
      FeatureMap<> featureFile; // also for results
      OpenMS::MRMFeatureFinderScoring::TransitionGroupMapType transition_group_map; // for results
      boost::shared_ptr<MSExperiment<Peak1D> > swath_map(new MSExperiment<Peak1D>);
      OpenSwath::SpectrumAccessPtr swath_ptr = SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(swath_map);

      boost::shared_ptr<MSExperiment<Peak1D> > xic_map(new MSExperiment<Peak1D>); // the map with the extracted ion chromatograms
      xic_map->setChromatograms(chromatograms);
      OpenSwath::SpectrumAccessPtr chromatogram_ptr = OpenSwath::SpectrumAccessPtr(new OpenMS::SpectrumAccessOpenMS(xic_map));
      TransformationDescription empty_trafo;

      featureFinder.setStrictFlag(false); // TODO remove this, it should be strict (e.g. all transitions need to be present for RT norm)
      featureFinder.pickExperiment(chromatogram_ptr, featureFile, transition_exp_used, empty_trafo, swath_ptr, transition_group_map);

      // find best feature, compute pairs of iRT and real RT
      simple_find_best_feature(transition_group_map, pairs, PeptideRTMap);

      std::vector<std::pair<double, double> > pairs_corrected;
      pairs_corrected = MRMRTNormalizer::rm_outliers(pairs, min_rsq, min_coverage);

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
      
    /// Simple method to find the best feature among a set of features (for the RT-normalization peptides)
    // TODO shared code!! -> OpenSwathRTNormalizer...
    void simple_find_best_feature(OpenMS::MRMFeatureFinderScoring::TransitionGroupMapType & transition_group_map, 
        std::vector<std::pair<double, double> > & pairs, std::map<OpenMS::String, double> PeptideRTMap)
    {
      for (OpenMS::MRMFeatureFinderScoring::TransitionGroupMapType::iterator trgroup_it = transition_group_map.begin();
          trgroup_it != transition_group_map.end(); trgroup_it++)
      {
        // we need at least one feature to find the best one
        OpenMS::MRMFeatureFinderScoring::MRMTransitionGroupType * transition_group = &trgroup_it->second;
        if (transition_group->getFeatures().size() == 0)
        {
          throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__,
              "RT normalization: did not find any features for group " + transition_group->getTransitionGroupID());
        }

        // Find the feature with the highest score
        double bestRT = -1;
        double highest_score = -1000;
        for (std::vector<MRMFeature>::iterator mrmfeature = transition_group->getFeaturesMuteable().begin();
             mrmfeature != transition_group->getFeaturesMuteable().end(); mrmfeature++)
        {
          if (mrmfeature->getOverallQuality() > highest_score)
          {
            bestRT = mrmfeature->getRT();
            highest_score = mrmfeature->getOverallQuality();
          }
        }
        String pepref = trgroup_it->second.getTransitions()[0].getPeptideRef();
        pairs.push_back(std::make_pair(bestRT, PeptideRTMap[pepref]));
      }
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
        FeatureMap<Feature>& output, OpenSwathTSVWriter & tsv_writer)
    {
      typedef OpenSwath::LightTransition TransitionType;
      // a transition group holds the MSSpectra with the Chromatogram peaks from above
      typedef MRMTransitionGroup<MSSpectrum <ChromatogramPeak>, TransitionType> MRMTransitionGroupType; 
      typedef std::map<String, MRMTransitionGroupType> TransitionGroupMapType;
      // this is the type in which we store the chromatograms for this analysis
      typedef MSSpectrum<ChromatogramPeak> RichPeakChromatogram; 

      double expected_rt;
      TransformationDescription trafo_inv = trafo;
      trafo_inv.invert();

      MRMFeatureFinderScoring featureFinder;
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
      for (Size i = 0; i < transition_exp.getPeptides().size(); i++)
      {
        assay_peptide_map[transition_exp.getPeptides()[i].id] = boost::numeric_cast<int>(i);
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
      for (AssayMapT::iterator assay_it = assay_map.begin(); assay_it != assay_map.end(); assay_it++)
      {
        // Create new MRMTransitionGroup
        String id = assay_it->first;
        MRMTransitionGroupType transition_group;
        transition_group.setTransitionGroupID(id);
        expected_rt = transition_exp.getPeptides()[ assay_peptide_map[id] ].rt;

        // Go through all transitions, for each transition get chromatogram and
        // the chromatogram and the assay to the MRMTransitionGroup
        for (Size i = 0; i < assay_it->second.size(); i++)
        {
          const TransitionType* transition = assay_it->second[i];

          if (chromatogram_map.find(transition->getNativeID()) == chromatogram_map.end() )
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

        // Process the MRMTransitionGroup: find peakgroups and score them
        trgroup_picker.pickTransitionGroup(transition_group);
        featureFinder.scorePeakgroups(transition_group, trafo, swath_map, output);

        // Add to the output tsv if given
        if (tsv_writer.isActive())
        {
          const OpenSwath::LightPeptide pep = transition_exp.getPeptides()[ assay_peptide_map[id] ];
          const TransitionType* transition = assay_it->second[0];
          to_output.push_back(tsv_writer.prepareLine(pep, transition, output, id));
        }
      }

      // Only write at the very end since this is a step needs a barrier
      if(tsv_writer.isActive())
      {
#ifdef _OPENMP
#pragma omp critical (scoreAll)
#endif
        {
          tsv_writer.writeLines(to_output);
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
      std::map<String, OpenSwath::LightPeptide* > trans_peptide_map;
      for (Size i = 0; i < transition_exp_used.getPeptides().size(); i++)
      {
        trans_peptide_map[transition_exp_used.getPeptides()[i].id] = &transition_exp_used.getPeptides()[i];
      }

      // Determine iteration size (nr peptides or nr transitions)
      Size itersize;
      if (ms1) {itersize = transition_exp_used.getPeptides().size();}
      else     {itersize = transition_exp_used.getTransitions().size();}

      for (Size i = 0; i < itersize; i++)
      {
        OpenSwath::ChromatogramPtr s(new OpenSwath::Chromatogram);
        output_chromatograms.push_back(s);

        ChromatogramExtractor::ExtractionCoordinates coord;
        OpenSwath::LightPeptide pep; // TargetedExperiment::Peptide pep;
        OpenSwath::LightTransition transition;

        if (ms1) 
        {
          pep = transition_exp_used.getPeptides()[i];
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

  };

}

// OpenMS base classes
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
  @page TOPP_OpenSwathWorkflow Workflow

  @brief Complete workflow to run OpenSWATH

  This implements the OpenSwath workflow as described in Roest and Rosenberger
  et al. (2013) and provides a complete, integrated analysis tool without the
  need to run multiple tools consecutively.

  It executes the following steps in order:

  <ul>
    <li>Reading of input files, can be provided as one single mzML or multiple "split" mzML (one per SWATH) </li>
    <li>Computing the retention time transformation using RT-normalization peptides</li>
    <li>Reading of the transition list</li>
    <li>Extracting the specified transitions</li>
    <li>Scoring the peakgroups in the extracted ion chromatograms (XIC)</li>
    <li>Reporting the peakgroups and the chromatograms</li>
  </ul>

  Look at the INI file (via "OpenSwathWorkflow -write_ini myini.ini") to see the available parameters and more functionality.

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
  -out_tsv and not -out_features. For downstream analysis (e.g. using mProphet)
  also the -out_tsv format is recommended.

  In addition, the extracted chromatograms can be written out using the
  -out_chrom parameter.

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES
class TOPPOpenSwathWorkflow 
  : public TOPPBase
{
public:

  TOPPOpenSwathWorkflow() 
    : TOPPBase("OpenSwathWorkflow", "Complete workflow to run OpenSWATH", true)
  {
  }

protected:

  void registerOptionsAndFlags_()
  {
    registerInputFileList_("in", "<files>", StringList(), "Input files separated by blank");
    setValidFormats_("in", StringList::create("mzML,mzXML"));

    registerInputFile_("tr", "<file>", "", "transition file ('TraML' or 'csv')");
    setValidFormats_("tr", StringList::create("csv,traML"));

    // one of the following two needs to be set
    registerInputFile_("tr_irt", "<file>", "", "transition file ('TraML' or 'csv')", false);
    setValidFormats_("tr_irt", StringList::create("csv,traML"));

    registerInputFile_("rt_norm", "<file>", "", "RT normalization file (how to map the RTs of this run to the ones stored in the library). If set, tr_irt may be omitted.", false, true);
    setValidFormats_("rt_norm", StringList::create("trafoXML"));

    registerStringOption_("swath_windows_file", "<file>", "", "Optional, tab separated file containing the SWATH windows: lower_offset upper_offset \\newline 400 425 \\newline ... ", false, true);

    // one of the following two needs to be set
    registerOutputFile_("out_features", "<file>", "", "output file", false);
    setValidFormats_("out_features", StringList::create("featureXML"));

    registerStringOption_("out_tsv", "<file>", "", "TSV output file (mProphet compatible)", false);

    registerOutputFile_("out_chrom", "<file>", "", "Also output all computed chromatograms (chrom.mzML) output", false, true);
    setValidFormats_("out_chrom", StringList::create("mzML"));

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
    setValidStrings_("readOptions", StringList::create("normal,cache"));

    // TODO terminal slash !
    registerStringOption_("tempDirectory", "<tmp>", "/tmp/", "Temporary directory to store cached files for example", false, true);

    registerStringOption_("extraction_function", "<name>", "tophat", "Function used to extract the signal", false, true);
    setValidStrings_("extraction_function", StringList::create("tophat,bartlett"));

    registerIntOption_("batchSize", "<number>", 0, "The batch size of chromatograms to process (0 means to only have one batch, sensible values are around 500-1000)", false, true);
    setMinInt_("batchSize", 0);

    registerSubsection_("Scoring", "Scoring parameters section");
  }

  Param getSubsectionDefaults_(const String & name) const
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
      feature_finder_param.remove("TransitionGroupPicker:background_subtraction");
      feature_finder_param.remove("TransitionGroupPicker:stop_after_intensity_ratio");

      // Peak Picker
      feature_finder_param.setValue("TransitionGroupPicker:PeakPickerMRM:use_gauss", "false");
      feature_finder_param.setValue("TransitionGroupPicker:PeakPickerMRM:sgolay_polynomial_order", 3);
      feature_finder_param.setValue("TransitionGroupPicker:PeakPickerMRM:sgolay_frame_length", 11);
      feature_finder_param.setValue("TransitionGroupPicker:PeakPickerMRM:peak_width", -1.0);
      feature_finder_param.setValue("TransitionGroupPicker:PeakPickerMRM:remove_overlapping_peaks", "true");
      // TODO it seems that the legacy method produces slightly larger peaks, e.g. it will not cut off peaks too early
      // however the same can be achieved by using a relatively low SN cutoff in the -Scoring:TransitionGroupPicker:PeakPickerMRM:signal_to_noise 0.5
      feature_finder_param.setValue("TransitionGroupPicker:recalculate_peaks_max_z", 0.75);
      feature_finder_param.setValue("TransitionGroupPicker:PeakPickerMRM:method", "corrected");
      feature_finder_param.setValue("TransitionGroupPicker:PeakPickerMRM:signal_to_noise", 0.1);
      feature_finder_param.setValue("TransitionGroupPicker:PeakPickerMRM:gauss_width", 30);
      feature_finder_param.remove("TransitionGroupPicker:PeakPickerMRM:gauss_width");
      feature_finder_param.remove("TransitionGroupPicker:PeakPickerMRM:sn_win_len");
      feature_finder_param.remove("TransitionGroupPicker:PeakPickerMRM:sn_bin_count");

      // EMG Scoring - turn off by default since it is very CPU-intensive
      feature_finder_param.remove("Scores:use_elution_model_score");
      feature_finder_param.setValue("EMGScoring:max_iteration", 10);
      feature_finder_param.setValue("EMGScoring:deltaRelError", 0.1);
      feature_finder_param.remove("EMGScoring:interpolation_step");
      feature_finder_param.remove("EMGScoring:tolerance_stdev_bounding_box");
      feature_finder_param.remove("EMGScoring:deltaAbsError");
         
      // remove these parameters
      feature_finder_param.remove("stop_report_after_feature");
      feature_finder_param.remove("add_up_spectra");
      feature_finder_param.remove("spacing_for_spectra_resampling");
      feature_finder_param.remove("EMGScoring:statistics:mean");
      feature_finder_param.remove("EMGScoring:statistics:variance");
      return feature_finder_param;
    }
    else
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Unknown subsection", name);
    }
  }

  void readSwathWindows(String filename, std::vector<double> & swath_prec_lower_,
    std::vector<double> & swath_prec_upper_ )
  {
    std::ifstream data(filename.c_str());
    std::string   line;
    std::string   tmp;
    std::getline(data, line); //skip header
    double lower, upper;
    while (std::getline(data, line))
    {
      std::stringstream lineStream(line);

      lineStream >> lower;
      lineStream >> upper;

      swath_prec_lower_.push_back(lower);
      swath_prec_upper_.push_back(upper);
    }
    assert(swath_prec_lower_.size() == swath_prec_upper_.size());
  }

  void annotateSwathMapsFromFile(String filename,
    std::vector< OpenSwath::SwathMap >& swath_maps)
  {
    std::vector<double> swath_prec_lower_, swath_prec_upper_;
    readSwathWindows(filename, swath_prec_lower_, swath_prec_upper_);
    assert(swath_prec_lower_.size() == swath_maps.size());
    for (Size i = 0; i < swath_maps.size(); i++)
    {
      swath_maps[i].lower = swath_prec_lower_[i];
      swath_maps[i].upper = swath_prec_upper_[i];
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

  TransformationDescription loadTrafoFile(String trafo_in, String irt_tr_file,
    const std::vector< OpenSwath::SwathMap > & swath_maps, double min_rsq, double min_coverage, 
    const Param& feature_finder_param, const OpenSwathWorkflow::ChromExtractParams& cp_irt)
  {
    TransformationDescription trafo_rtnorm;
    if (trafo_in.size() > 0) 
    {
      // get read RT normalization file
      TransformationXMLFile trafoxml;
      trafoxml.load(trafo_in, trafo_rtnorm);
      Param model_params = getParam_().copy("model:", true);
      model_params.setValue("symmetric_regression", "false");
      String model_type = "linear";
      trafo_rtnorm.fitModel(model_type, model_params);
    }
    else
    {
      OpenSwathWorkflow wf;
      wf.setLogType(log_type_);
      // Loading iRT file
      std::cout << "Will load iRT transitions and try to find iRT peptides" << std::endl;
      TraMLFile traml;
      OpenMS::TargetedExperiment irt_transitions;
      traml.load(irt_tr_file, irt_transitions);
      trafo_rtnorm = wf.performRTNormalization(irt_transitions, swath_maps, min_rsq, min_coverage, 
          feature_finder_param, cp_irt);
    }
    return trafo_rtnorm;
  }

  int computeExpectedChromatograms(const std::vector< OpenSwath::SwathMap > & swath_maps, 
    const OpenSwath::LightTargetedExperiment & transition_exp)
  {
    int expected_chromatograms = 0;
    for (Size i = 0; i < transition_exp.transitions.size(); i++)
    { 
      for (Size j = 0; j < swath_maps.size(); j++)
      {
        if (!swath_maps[j].ms1 && transition_exp.transitions[i].precursor_mz >= swath_maps[j].lower 
          && transition_exp.transitions[i].precursor_mz <= swath_maps[j].upper)
        {
          // here we just check whether there is a SWATH from which we could
          // potentially extract this transition, if we find one we abort (e.g.
          // dont consider min_upper_edge_dist here).
          expected_chromatograms++;
          break;
        }
      }
    }
    return expected_chromatograms;
  }

  ExitCodes main_(int, const char **)
  {
    ///////////////////////////////////
    // Prepare Parameters
    ///////////////////////////////////
    StringList file_list = getStringList_("in");
    String tr_file = getStringOption_("tr");

    String out = getStringOption_("out_features");
    String out_tsv = getStringOption_("out_tsv");

    String irt_tr_file = getStringOption_("tr_irt");
    String trafo_in = getStringOption_("rt_norm");

    String out_chrom = getStringOption_("out_chrom");
    bool ppm = getFlag_("ppm");
    bool split_file = getFlag_("split_file_input");
    bool use_emg_score = getFlag_("use_elution_model_score");
    DoubleReal min_upper_edge_dist = getDoubleOption_("min_upper_edge_dist");
    DoubleReal mz_extraction_window = getDoubleOption_("mz_extraction_window");
    DoubleReal rt_extraction_window = getDoubleOption_("rt_extraction_window");
    DoubleReal extra_rt_extract = getDoubleOption_("extra_rt_extraction_window");
    String extraction_function = getStringOption_("extraction_function");
    String swath_windows_file = getStringOption_("swath_windows_file");
    int batchSize = (int)getIntOption_("batchSize");

    DoubleReal min_rsq = getDoubleOption_("min_rsq");
    DoubleReal min_coverage = getDoubleOption_("min_coverage");

    String readoptions = getStringOption_("readOptions");
    String tmp = getStringOption_("tempDirectory");

    if (trafo_in.empty() && irt_tr_file.empty()) 
          throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__,
              "Either rt_norm or tr_irt needs to be set");
    if ((out.empty() && out_tsv.empty()) || (!out.empty() && !out_tsv.empty()) ) 
          throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__,
              "Either out_features or out_tsv needs to be set (but not both)");

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
    if (use_emg_score) { feature_finder_param.setValue("Scores:use_elution_model_score", "true");}
    else { feature_finder_param.setValue("Scores:use_elution_model_score", "false");}

    ///////////////////////////////////
    // Load the SWATH files
    ///////////////////////////////////
    boost::shared_ptr<ExperimentalSettings > exp_meta(new ExperimentalSettings);
    std::vector< OpenSwath::SwathMap > swath_maps;
    loadSwathFiles(file_list, split_file, tmp, readoptions, exp_meta, swath_maps);

    // Allow the user to specify the SWATH windows
    if (!swath_windows_file.empty())
      annotateSwathMapsFromFile(swath_windows_file, swath_maps);

    for (Size i = 0; i < swath_maps.size(); i++)
        LOG_DEBUG << "Found swath map " << i << " with lower " << swath_maps[i].lower << " and upper " << swath_maps[i].upper << std::endl;

    ///////////////////////////////////
    // Get the transformation information (using iRT peptides)
    ///////////////////////////////////
    TransformationDescription trafo_rtnorm = loadTrafoFile(trafo_in, irt_tr_file,
        swath_maps, min_rsq, min_coverage, feature_finder_param, cp_irt);

    ///////////////////////////////////
    // Load the transitions
    ///////////////////////////////////
    OpenSwath::LightTargetedExperiment transition_exp;
    ProgressLogger progresslogger;
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
      TransitionTSVReader().convertTSVToTargetedExperiment(tr_file.c_str(), transition_exp);
    }
    progresslogger.endProgress();

    ///////////////////////////////////
    // Set up chrom.mzML output
    ///////////////////////////////////
    MSDataWritingConsumer * chromConsumer;
    int expected_chromatograms = 0;
    if (!out_chrom.empty())
    {
      chromConsumer = new PlainMSDataWritingConsumer(out_chrom);
      expected_chromatograms = computeExpectedChromatograms(swath_maps, transition_exp); 
      chromConsumer->setExpectedSize(0, expected_chromatograms);
      chromConsumer->setExperimentalSettings(*exp_meta);
      chromConsumer->addDataProcessing(getProcessingInfo_(DataProcessing::SMOOTHING));
    }
    else
    {
      chromConsumer = new NoopMSDataWritingConsumer(out_chrom);
    }

    ///////////////////////////////////
    // Extract and score
    ///////////////////////////////////
    FeatureMap<> out_featureFile;

    OpenSwathTSVWriter tsvwriter(out_tsv, file_list[0]);
    OpenSwathWorkflow wf;
    wf.setLogType(log_type_);

    wf.performExtraction(swath_maps, trafo_rtnorm, cp, feature_finder_param, transition_exp,
        out_featureFile, out, tsvwriter, chromConsumer, batchSize);
    if (!out.empty())
    {
      addDataProcessing_(out_featureFile, getProcessingInfo_(DataProcessing::QUANTITATION));
      out_featureFile.ensureUniqueId();
      FeatureXMLFile().store(out, out_featureFile);
    }

    // Check that the number in <chromatogramList count=...> is equal to the
    // number of actually written chromatograms.
    if (!out_chrom.empty() && (int)chromConsumer->getNrChromatogramsWritten() != expected_chromatograms)
    {
      std::cerr << "Expected to extract " << transition_exp.transitions.size() << 
        " chromatograms, however " << chromConsumer->getNrChromatogramsWritten() <<  
        " were written to disk. Something is off here!" << std::endl;
      if ( chromConsumer->getNrChromatogramsWritten() < transition_exp.transitions.size() )
      {
      std::cerr << "Will try to rescue by writing " << 
        transition_exp.transitions.size() - chromConsumer->getNrChromatogramsWritten() <<  
        " extra empty chromatograms." << std::endl;
        for (Size i = 0; i < transition_exp.transitions.size() - chromConsumer->getNrChromatogramsWritten(); i++)
        {
          OpenMS::MSChromatogram<> c;
          chromConsumer->consumeChromatogram(c);
        }
      }
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

