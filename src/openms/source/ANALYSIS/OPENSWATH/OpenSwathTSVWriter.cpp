// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathTSVWriter.h>

namespace OpenMS
{


  OpenSwathTSVWriter::OpenSwathTSVWriter(const String& output_filename,
                                         const String& input_filename,
                                         bool ms1_scores) :
    ofs(output_filename.c_str()),
    input_filename_(input_filename),
    doWrite_(!output_filename.empty()),
    use_ms1_traces_(ms1_scores)
    {
    }

    bool OpenSwathTSVWriter::isActive() const
    {
      return doWrite_;
    }

    void OpenSwathTSVWriter::writeHeader()
    {
      ofs << "transition_group_id" << "\t"
          << "peptide_group_label" << "\t"
          << "run_id" << "\t"
          << "filename"<< "\t"\
          << "RT" << "\t"
          << "id" << "\t"
          << "Sequence"<< "\t"
          << "MC" << "\t"
          << "FullPeptideName" << "\t"
          << "Charge" << "\t"
          << "m/z" << "\t"
          << "Intensity" << "\t"
          << "ProteinName" << "\t"
          << "GeneName" << "\t"
          << "decoy" << "\t"
          << "assay_rt" <<"\t"
          << "delta_rt" << "\t"
          << "leftWidth" <<
        "\tmain_var_xx_swath_prelim_score\tnorm_RT\tnr_peaks\tpeak_apices_sum\tpotentialOutlier\tinitialPeakQuality" <<
        "\trightWidth\trt_score\tsn_ratio\ttotal_xic\tvar_bseries_score\tvar_dotprod_score" <<
        "\tvar_intensity_score\tvar_isotope_correlation_score\tvar_isotope_overlap_score" <<
        "\tvar_library_corr\tvar_library_dotprod\tvar_library_manhattan\tvar_library_rmsd" <<
        "\tvar_library_rootmeansquare\tvar_library_sangle\tvar_log_sn_score\tvar_manhatt_score" <<
        "\tvar_massdev_score\tvar_massdev_score_weighted\tvar_norm_rt_score\tvar_xcorr_coelution" <<
        "\tvar_xcorr_coelution_weighted\tvar_xcorr_shape\tvar_xcorr_shape_weighted" <<
        "\tvar_im_xcorr_shape\tvar_im_xcorr_coelution\tvar_im_delta_score\tvar_im_ms1_delta_score" <<
        "\tim_drift_ms2\tim_drift_weighted_ms2\tdelta_im_ms2\tim_drift_ms1\tim_delta_ms1" <<
        "\tvar_yseries_score\tvar_elution_model_fit_score";

      if (use_ms1_traces_)
      {
        ofs << "\tvar_ms1_ppm_diff\tvar_ms1_isotope_corr\tvar_ms1_isotope_overlap\tvar_ms1_xcorr_coelution\tvar_ms1_xcorr_shape";
      }
      ofs << "\txx_lda_prelim_score\txx_swath_prelim_score";
      if (use_ms1_traces_)
      {
        ofs << "\taggr_prec_Peak_Area\taggr_prec_Peak_Apex\taggr_prec_Annotation";
      }
      ofs << "\taggr_Peak_Area\taggr_Peak_Apex\taggr_Fragment_Annotation" << "\t"
          << "rt_fwhm" << "\t"
          << "masserror_ppm";
      ofs << "\n";
    }

    String OpenSwathTSVWriter::prepareLine(const OpenSwath::LightCompound& pep,
        const OpenSwath::LightTransition * transition,
        const FeatureMap& output, const String& id) const
    {
        String result = "";
        String decoy = "0"; // 0 = false
        if (transition->decoy)
        {
          decoy = "1";
        }

        // iterator over MRMFeatures
        for (FeatureMap::const_iterator feature_it = output.begin(); feature_it != output.end(); ++feature_it)
        {
          StringList aggr_Peak_Area, aggr_Peak_Apex, aggr_Fragment_Annotation;
          StringList aggr_prec_Peak_Area, aggr_prec_Peak_Apex, aggr_prec_Fragment_Annotation;
          StringList rt_fwhm;
          String gene_name;
          for (std::vector<Feature>::const_iterator sub_it = feature_it->getSubordinates().begin(); sub_it != feature_it->getSubordinates().end(); ++sub_it)
          {
            if (sub_it->metaValueExists("FeatureLevel"))
            {
              if (sub_it->getMetaValue("FeatureLevel") == "MS2")
              {
                aggr_Peak_Area.push_back((String)sub_it->getIntensity());
                aggr_Peak_Apex.push_back(String((double)sub_it->getMetaValue("peak_apex_int")));
                aggr_Fragment_Annotation.push_back((String)sub_it->getMetaValue("native_id"));
                rt_fwhm.push_back((String)sub_it->getMetaValue("width_at_50"));
              }
              else if (sub_it->getMetaValue("FeatureLevel") == "MS1")
              {
                aggr_prec_Peak_Area.push_back((String)sub_it->getIntensity());
                aggr_prec_Peak_Apex.push_back(String((double)sub_it->getMetaValue("peak_apex_int")));
                aggr_prec_Fragment_Annotation.push_back((String)sub_it->getMetaValue("native_id"));
              }
            }
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
          // legacy fix since there are many TraMLs floating around which have "light"/"NA" in there
          if (group_label.empty() || group_label == "light" || group_label == "NA")
          {
            group_label = id;
          }

          // If a protein is present, take the first one
          String protein_name = "";
          if (!pep.protein_refs.empty() )
          {
            protein_name = pep.protein_refs[0];
          }

          String main_var = "0";
          if (feature_it->metaValueExists("main_var_xx_swath_prelim_score"))
          {
            main_var = (String)feature_it->getMetaValue("main_var_xx_swath_prelim_score");
          }
          else if (feature_it->metaValueExists("main_var_xx_lda_prelim_score"))
          {
            main_var = (String)feature_it->getMetaValue("main_var_xx_lda_prelim_score");
          }

          String line =
            id + "_run0"
            + "\t" + group_label
            + "\t" + "0"
            + "\t" + input_filename_
            + "\t" + (String)feature_it->getRT()
            + "\t" + "f_" + feature_it->getUniqueId()  // TODO might not be unique!!!
            + "\t" + (String)pep.sequence
            + "\t" + (feature_it->metaValueExists("missedCleavages") ? (String)feature_it->getMetaValue("missedCleavages") : "")
            + "\t" + full_peptide_name
            + "\t" + (String)pep.charge
            + "\t" + (String)transition->precursor_mz
            + "\t" + (String)feature_it->getIntensity()
            + "\t" + protein_name
            + "\t" + gene_name
            + "\t" + decoy
            // Note: missing MetaValues will just produce a DataValue::EMPTY which lead to an empty column
            + "\t" + (String)feature_it->getMetaValue("assay_rt")
            + "\t" + (String)feature_it->getMetaValue("delta_rt")
            + "\t" + (String)feature_it->getMetaValue("leftWidth")
            + "\t" + main_var
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

            + "\t" + (String)feature_it->getMetaValue("var_im_xcorr_shape")
            + "\t" + (String)feature_it->getMetaValue("var_im_xcorr_coelution")
            + "\t" + (String)feature_it->getMetaValue("var_im_delta_score")
            + "\t" + (String)feature_it->getMetaValue("var_im_ms1_delta_score")
            + "\t" + (String)feature_it->getMetaValue("im_drift")
            + "\t" + (String)feature_it->getMetaValue("im_drift_weighted")
            + "\t" + (String)feature_it->getMetaValue("im_delta")
            + "\t" + (String)feature_it->getMetaValue("im_ms1_drift")
            + "\t" + (String)feature_it->getMetaValue("im_ms1_delta")

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
              line += "\t" + ListUtils::concatenate(aggr_prec_Peak_Area, ";") + "\t" + ListUtils::concatenate(aggr_prec_Peak_Apex, ";") + "\t" + ListUtils::concatenate(aggr_prec_Fragment_Annotation, ";");
            }
            line += "\t" + ListUtils::concatenate(aggr_Peak_Area, ";") + "\t" + ListUtils::concatenate(aggr_Peak_Apex, ";") + "\t" + ListUtils::concatenate(aggr_Fragment_Annotation, ";");
            line += "\t" + ListUtils::concatenate(rt_fwhm, ";");
            line += "\t" + (feature_it->metaValueExists("masserror_ppm") ? ListUtils::concatenate(feature_it->getMetaValue("masserror_ppm").toDoubleList(), ";") : "");

            line += "\n";
            result += line;
        } // end of iteration
      return result;
    }

    void OpenSwathTSVWriter::writeLines(const std::vector<String>& to_output)
    {
      for (const auto& s : to_output) ofs << s;
    }

}

