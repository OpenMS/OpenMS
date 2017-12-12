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
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathTSVWriter.h>

namespace OpenMS
{



  OpenSwathTSVWriter::OpenSwathTSVWriter(String output_filename, 
                       String input_filename,
                       bool ms1_scores, 
                       bool sonar, 
                       bool uis_scores) :
      ofs(output_filename.c_str()),
      input_filename_(input_filename),
      doWrite_(!output_filename.empty()),
      use_ms1_traces_(ms1_scores),
      sonar_(sonar),
      enable_uis_scoring_(uis_scores)
      {}

    bool OpenSwathTSVWriter::isActive()
    {
      return doWrite_;
    }

    void OpenSwathTSVWriter::writeHeader()
    {
      ofs << "transition_group_id\tpeptide_group_label\trun_id\tfilename\tRT\tid\tSequence\tFullPeptideName" <<
        "\tCharge\tm/z\tIntensity\tProteinName\tdecoy\tassay_rt\tdelta_rt\tleftWidth" <<
        "\tmain_var_xx_swath_prelim_score\tnorm_RT\tnr_peaks\tpeak_apices_sum\tpotentialOutlier\tinitialPeakQuality" <<
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
      if (sonar_)
      {
        ofs << "\tvar_sonar_lag\tvar_sonar_shape\tvar_sonar_log_sn\tvar_sonar_log_diff\tvar_sonar_log_trend\tvar_sonar_rsq";
      }
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

    String OpenSwathTSVWriter::prepareLine(const OpenSwath::LightCompound& pep,
        const OpenSwath::LightTransition * transition,
        const FeatureMap& output, const String id)
    {
        String result = "";
        String decoy = "0"; // 0 = false
        if (transition->decoy)
        {
          decoy = "1";
        }

        for (FeatureMap::const_iterator feature_it = output.begin(); feature_it != output.end(); ++feature_it)
        {

          char intensity_char[40];
          char intensity_apex_char[40];
          String aggr_Peak_Area = "";
          String aggr_Peak_Apex = "";
          String aggr_Fragment_Annotation = "";
          String aggr_prec_Peak_Area = "";
          String aggr_prec_Peak_Apex = "";
          String aggr_prec_Fragment_Annotation = "";
          for (std::vector<Feature>::const_iterator sub_it = feature_it->getSubordinates().begin(); sub_it != feature_it->getSubordinates().end(); ++sub_it)
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
                full_peptide_name += "(UniMod:" + String(pep.modifications[modloc].unimod_id) + ")";
              }
            }
          }

          // Compute peptide group label (use the provided label or use the
          // transition group).
          String group_label = pep.peptide_group_label;
          if (group_label.empty()) group_label = id;
          if (group_label == "light") group_label = id; // legacy fix since there are many TraMLs floating around which have "light" in there
          if (group_label == "NA") group_label = id; // legacy fix since there are many TraMLs floating around which have "NA" in there

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
            + "\t" + (String)feature_it->getMetaValue("initialPeakQuality")
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
            if (sonar_)
            {
              line += "\t" + (String)feature_it->getMetaValue("var_sonar_lag")
              + "\t" + (String)feature_it->getMetaValue("var_sonar_shape")
              + "\t" + (String)feature_it->getMetaValue("var_sonar_log_sn")
              + "\t" + (String)feature_it->getMetaValue("var_sonar_log_diff")
              + "\t" + (String)feature_it->getMetaValue("var_sonar_log_trend")
              + "\t" + (String)feature_it->getMetaValue("var_sonar_rsq");

            }
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

    void OpenSwathTSVWriter::writeLines(std::vector<String> to_output)
    {
      for (Size i = 0; i < to_output.size(); i++) { ofs << to_output[i]; }
    }

}


