// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: George Rosenberger $
// $Authors: George Rosenberger $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/MRMIonSeries.h>

#include <boost/assign.hpp>
#include <boost/lexical_cast.hpp>

namespace OpenMS
{
  MRMIonSeries::MRMIonSeries() = default;

  MRMIonSeries::~MRMIonSeries() = default;

  std::pair<String, double> MRMIonSeries::getIon(IonSeries& ionseries, const String& ionid)
  {
    if (ionseries.find(ionid) != ionseries.end())
    {
      return make_pair(ionid, ionseries[ionid]);
    }
    else
    {
      return make_pair(String("unannotated"), -1);
    }
  }

  std::pair<String, double> MRMIonSeries::annotateIon(const IonSeries& ionseries, const double ProductMZ, const double mz_threshold)
  {
    // make sure to only use annotated transitions and to use the theoretical MZ
    using namespace boost::assign;

    // Iterate over ion type and then ordinal
    std::pair<String, double> ion;
    String unannotated = "unannotated";
    ion = make_pair(unannotated, -1);
    double closest_delta = std::numeric_limits<double>::max();

    for (const auto& ordinal : ionseries)
    {
      if (std::fabs(ordinal.second - ProductMZ) <= mz_threshold && std::fabs(ordinal.second - ProductMZ) <= closest_delta)
      {
        closest_delta = std::fabs(ordinal.second - ProductMZ);
        ion = make_pair(ordinal.first, ordinal.second);
      }
    }

    return ion;
  }

  TargetedExperiment::Interpretation MRMIonSeries::annotationToCVTermList_(const String& annotation)
  {
    // CVTermList interpretation;
    TargetedExperiment::Interpretation interpretation;

    String fragment_type;
    int fragment_nr = -1;
    double fragment_loss = 0;
    // int fragment_gain = 0;

    std::vector<String> best_annotation;
    annotation.split("/", best_annotation);

    if (best_annotation[0] == "Precursor_i0" || best_annotation[0] == "MS2_Precursor_i0")
    {
      return interpretation;
    }
    else if (best_annotation[0].find("-") != std::string::npos)
    {
      std::vector<String> best_annotation_loss;
      best_annotation[0].split("-", best_annotation_loss);

      fragment_type = best_annotation_loss[0].substr(0, 1);
      fragment_nr = best_annotation_loss[0].substr(1).toInt();

      // SpectraST style neutral loss
      try
      {
        int nl = boost::lexical_cast<int>(best_annotation_loss[1]);
        fragment_loss = -1 * nl;
      }
      catch (boost::bad_lexical_cast &)
      {
        static const EmpiricalFormula nl_formula(best_annotation_loss[1]);
        fragment_loss = -1 * nl_formula.getMonoWeight();
      }
    }
    else if (best_annotation[0].find("+") != std::string::npos)
    {
      std::vector<String> best_annotation_gain;
      best_annotation[0].split("+", best_annotation_gain);

      fragment_type = best_annotation_gain[0].substr(0, 1);
      fragment_nr = best_annotation_gain[0].substr(1).toInt();
      // fragment_gain = String(best_annotation_gain[1]).toInt(); // fragment neutral gain is not implemented as CV term.
    }
    else
    {
      fragment_type = best_annotation[0].substr(0, 1);
      fragment_nr = best_annotation[0].substr(1).toInt();
    }

    if (fragment_nr != -1)
    {
      interpretation.rank = 1; // we only store the best interpretation
    }

    if (fragment_nr != -1)
    {
      interpretation.ordinal = fragment_nr;
    }

    if (fragment_loss < 0)
    {
      CVTerm frag_loss;
      frag_loss.setCVIdentifierRef("MS");
      frag_loss.setAccession("MS:1001524");
      frag_loss.setName("fragment neutral loss");
      frag_loss.setValue(fragment_loss);
      interpretation.addCVTerm(frag_loss);
    }

    // figure out which fragment it is
    if (fragment_type == "x")
    {
      interpretation.iontype = TargetedExperiment::IonType::XIon;
    }
    else if (fragment_type == "y")
    {
      interpretation.iontype = TargetedExperiment::IonType::YIon;
    }
    else if (fragment_type == "z")
    {
      interpretation.iontype = TargetedExperiment::IonType::ZIon;
    }
    else if (fragment_type == "a")
    {
      interpretation.iontype = TargetedExperiment::IonType::AIon;
    }
    else if (fragment_type == "b")
    {
      interpretation.iontype = TargetedExperiment::IonType::BIon;
    }
    else if (fragment_type == "c")
    {
      interpretation.iontype = TargetedExperiment::IonType::CIon;
    }
    else
    {
      interpretation.iontype = TargetedExperiment::IonType::NonIdentified;
    }

    return interpretation;
  }

  void MRMIonSeries::annotationToCV_(ReactionMonitoringTransition& tr)
  {
    OpenMS::ReactionMonitoringTransition::Product p = tr.getProduct();

    std::vector<String> best_annotation;
    tr.getMetaValue("annotation").toString().split("/", best_annotation);

    String annotation;
    if (best_annotation[0].find("^") != std::string::npos)
    {
      std::vector<String> best_annotation_charge;
      best_annotation[0].split("^", best_annotation_charge);
      p.setChargeState(String(best_annotation_charge[1]).toInt());
      annotation = best_annotation_charge[0];
    }
    else
    {
      p.setChargeState(1);
      annotation = best_annotation[0];
    }

    TargetedExperiment::Interpretation interpretation = annotationToCVTermList_(annotation);

    p.resetInterpretations();
    p.addInterpretation(interpretation);
    tr.setProduct(p);
  }

  void MRMIonSeries::annotateTransitionCV(ReactionMonitoringTransition& tr, const String& annotation)
  {
    tr.setMetaValue("annotation", annotation);
    annotationToCV_(tr);
  }

  void MRMIonSeries::annotateTransition(ReactionMonitoringTransition& tr, const TargetedExperiment::Peptide& peptide, const double precursor_mz_threshold, double product_mz_threshold, const bool enable_reannotation, const std::vector<String>& fragment_types, const std::vector<size_t>& fragment_charges, const bool enable_specific_losses, const bool enable_unspecific_losses, const int round_decPow)
  {
    OPENMS_PRECONDITION(peptide.hasCharge(), "Cannot annotate transition without a peptide charge state")
    // TODO: we should not have transitions without charge states
    // OPENMS_PRECONDITION(tr.isProductChargeStateSet(), "Cannot annotate transition without a charge state")

    TargetedExperiment::Interpretation interpretation;
    OpenMS::AASequence sequence = TargetedExperimentHelper::getAASequence(peptide);
    int precursor_charge = 1; // assume default to be 1 (should always be set, see precondition)
    int fragment_charge = 1; // assume default to be 1 (should always be set, see precondition)
    if (peptide.hasCharge())
    {
      precursor_charge = peptide.getChargeState();
    }
    if (tr.isProductChargeStateSet() )
    {
      fragment_charge = tr.getProductChargeState();
    }

    double prec_pos = sequence.getMZ(precursor_charge);
    bool unannotated = false;
    std::pair<String, double> target_ion = std::make_pair(String("unannotated"), -1);
    double pos = -1;
    String ionstring;

    if (!tr.getProduct().getInterpretationList().empty())
    {
      interpretation = tr.getProduct().getInterpretationList()[0];
      AASequence ion;

      if (interpretation.ordinal > 0) // if ordinal is set
      {
        int ordinal = (int)interpretation.ordinal;

        if (interpretation.iontype == TargetedExperiment::IonType::XIon)
        {
          ion = sequence.getSuffix(ordinal);
          ionstring += "x";
          pos = ion.getMZ(fragment_charge, Residue::XIon);
        }
        else if (interpretation.iontype == TargetedExperiment::IonType::YIon)
        {
          ion = sequence.getSuffix(ordinal);
          ionstring += "y";
          pos = ion.getMZ(fragment_charge, Residue::YIon);
        }
        else if (interpretation.iontype == TargetedExperiment::IonType::ZIon)
        {
          ion = sequence.getSuffix(ordinal);
          ionstring += "z";
          pos = ion.getMZ(fragment_charge, Residue::ZIon);
        }
        else if (interpretation.iontype == TargetedExperiment::IonType::AIon)
        {
          ion = sequence.getSuffix(ordinal);
          ionstring += "a";
          pos = ion.getMZ(fragment_charge, Residue::AIon);
        }
        else if (interpretation.iontype == TargetedExperiment::IonType::BIon)
        {
          ion = sequence.getSuffix(ordinal);
          ionstring += "b";
          pos = ion.getMZ(fragment_charge, Residue::BIon);
        }
        else if (interpretation.iontype == TargetedExperiment::IonType::CIon)
        {
          ion = sequence.getSuffix(ordinal);
          ionstring += "c";
          pos = ion.getMZ(fragment_charge, Residue::CIon);
        }
        else
        {
          unannotated = true;
        }
      }
      else
      {
        unannotated = true;
      }

      if (std::find(fragment_types.begin(), fragment_types.end(), ionstring) == fragment_types.end())
      {
        unannotated = true;
      }

      if (interpretation.ordinal > 0)
      {
        ionstring += String(interpretation.ordinal);
      }
      else
      {
        unannotated = true;
      }

      if (interpretation.hasCVTerm("MS:1001524") && (enable_specific_losses || enable_unspecific_losses)) // fragment ion neutral loss
      {
        double nl = interpretation.getCVTerms().at("MS:1001524")[0].getValue().toString().toDouble();
        // SpectraST style neutral losses
        if (nl == -18)
        {
          ionstring += "-H2O1";
          static const EmpiricalFormula neutralloss_h2o("H2O1"); // -18 H2O loss
          pos -= neutralloss_h2o.getMonoWeight() / fragment_charge;
        }
        else if (nl == -17)
        {
          ionstring += "-H3N1";
          static const EmpiricalFormula neutralloss_nh3("H3N1"); // -17 NH3 loss
          pos -= neutralloss_nh3.getMonoWeight() / fragment_charge;
        }
        else if (nl == -36)
        {
          ionstring += "-H4O2";
          static const EmpiricalFormula neutralloss_h2oh2o("H4O2"); // -36 2 * H2O loss
          pos -= neutralloss_h2oh2o.getMonoWeight() / fragment_charge;
        }
        else if (nl == -34)
        {
          ionstring += "-H6N2";
          static const EmpiricalFormula neutralloss_nh3nh3("H6N2"); // -34 2 * NH3 loss
          pos -= neutralloss_nh3nh3.getMonoWeight() / fragment_charge;
        }
        else if (nl == -35)
        {
          ionstring += "-H5N1O1";
          static const EmpiricalFormula neutralloss_h2onh3("H5N1O1"); // -35 H2O & NH3 loss
          pos -= neutralloss_h2onh3.getMonoWeight() / fragment_charge;
        }
        else if (nl == -64)
        {
          ionstring += "-C1H4O1S1";
          static const EmpiricalFormula neutralloss_ch4so("C1H4O1S1"); // -64 CH4SO loss
          pos -= neutralloss_ch4so.getMonoWeight() / fragment_charge;
        }
        else if (nl == -80)
        {
          ionstring += "-H1O3P1";
          static const EmpiricalFormula neutralloss_hpo3("H1O3P1"); // -80 HPO3 loss
          pos -= neutralloss_hpo3.getMonoWeight() / fragment_charge;
        }
        else if (nl == -98)
        {
          ionstring += "-H3O4P1";
          static const EmpiricalFormula neutralloss_hpo3h2o("H3O4P1"); // -98 HPO3 & H2O loss
          pos -= neutralloss_hpo3h2o.getMonoWeight() / fragment_charge;
        }
        else if (nl == -45)
        {
          ionstring += "-C1H3N1O1";
          static const EmpiricalFormula neutralloss_ch3no("C1H3N1O1"); // -45 CH3NO loss
          pos -= neutralloss_ch3no.getMonoWeight() / fragment_charge;
        }
        else if (nl == -44)
        {
          ionstring += "-C1O2";
          static const EmpiricalFormula neutralloss_co2("C1O2"); // -44 CO2 loss
          pos -= neutralloss_co2.getMonoWeight() / fragment_charge;
        }
        else if (nl == -46)
        {
          ionstring += "-C1H2O2";
          static const EmpiricalFormula neutralloss_hccoh("C1H2O2"); // -46 HCOOH loss
          pos -= neutralloss_hccoh.getMonoWeight() / fragment_charge;

        }
        // Double CV term (compatible with PSI CV terms)
        else if (nl < 0)
        {
          ionstring += String(Math::roundDecimal(nl, round_decPow));
          pos -= nl / fragment_charge;

        }
        else
        {
          unannotated = true;
        }
      }

      if (fragment_charge >= 1 && 
           std::find(fragment_charges.begin(), fragment_charges.end(), fragment_charge) != fragment_charges.end())
      {
        ionstring += "^" + String(fragment_charge);
        tr.setMetaValue("annotation", ionstring);
      }
      else
      {
        unannotated = true;
      }
    }
    else
    {
      unannotated = true;
    }

    if (enable_reannotation)
    {
      MRMIonSeries::IonSeries ionseries = getIonSeries(sequence, precursor_charge, fragment_types,
                                                       fragment_charges, enable_specific_losses, enable_unspecific_losses);
      target_ion = annotateIon(ionseries, tr.getProductMZ(), product_mz_threshold);
      ionstring = target_ion.first;
      tr.setMetaValue("annotation", ionstring);
      pos = Math::roundDecimal(target_ion.second, round_decPow);
      prec_pos = Math::roundDecimal(prec_pos, round_decPow);
      tr.setProductMZ(pos);
      tr.setPrecursorMZ(prec_pos);

      if (ionstring == "unannotated")
      {
        unannotated = true;
      }
      else
      {
        annotationToCV_(tr);
        interpretation = tr.getProduct().getInterpretationList()[0];
        unannotated = false;
      }
    }

    if (!unannotated && std::fabs(tr.getProductMZ() - pos) <= product_mz_threshold && std::fabs(tr.getPrecursorMZ() - prec_pos) <= precursor_mz_threshold)
    {
      CVTerm frag_mzdelta;
      frag_mzdelta.setCVIdentifierRef("MS");
      frag_mzdelta.setAccession("MS:1000904");
      frag_mzdelta.setName("product ion m/z delta");
      frag_mzdelta.setValue(std::fabs(Math::roundDecimal(tr.getProductMZ() - pos, round_decPow)));
      interpretation.replaceCVTerm(frag_mzdelta);
      pos = Math::roundDecimal(pos, round_decPow);
      prec_pos = Math::roundDecimal(prec_pos, round_decPow);
      tr.setProductMZ(pos);
      tr.setPrecursorMZ(prec_pos);
    }
    else
    {
      unannotated = true;
    }

    if (unannotated)
    {
      interpretation.iontype = TargetedExperiment::IonType::NonIdentified;
      tr.setMetaValue("annotation", "unannotated");
    }
    else
    {
      tr.setMetaValue("annotation", ionstring);
      annotationToCV_(tr);
    }

    OpenMS::ReactionMonitoringTransition::Product p = tr.getProduct();
    p.resetInterpretations();
    p.addInterpretation(interpretation);
    tr.setProduct(p);
  }

  std::unordered_map<String, double> MRMIonSeries::getIonSeries(const AASequence& sequence,
                                                                size_t precursor_charge,
                                                                const std::vector<String>& fragment_types,
                                                                const std::vector<size_t>& fragment_charges,
                                                                const bool enable_specific_losses,
                                                                const bool enable_unspecific_losses,
                                                                const int round_decPow)
  {
    const static EmpiricalFormula H2O = EmpiricalFormula("H2O1");
    const static EmpiricalFormula NH3 = EmpiricalFormula("H3N1");
    const static EmpiricalFormula CN2 = EmpiricalFormula("C1H2N2");
    const static EmpiricalFormula CNO = EmpiricalFormula("C1H2N1O1");

    std::unordered_map<String, double> ionseries;

    for (std::vector<String>::const_iterator ft_it = fragment_types.begin(); ft_it != fragment_types.end(); ++ft_it)
    {
      for (const auto& charge : fragment_charges)
      {
        if (charge > precursor_charge)
        {
          continue;
        }

        for (Size i = 1; i < sequence.size(); ++i)
        {
          double pos = 0;
          AASequence ion;

          if (*ft_it == "a")
          {
            ion = sequence.getPrefix(i);
            pos = ion.getMZ(charge, Residue::AIon);
          }
          else if (*ft_it == "b")
          {
            ion = sequence.getPrefix(i);
            pos = ion.getMZ(charge, Residue::BIon);
          }
          else if (*ft_it == "c")
          {
            ion = sequence.getPrefix(i);
            pos = ion.getMZ(charge, Residue::CIon);
          }
          else if (*ft_it == "x")
          {
            ion = sequence.getSuffix(i);
            pos = ion.getMZ(charge, Residue::XIon);
          }
          else if (*ft_it == "y")
          {
            ion = sequence.getSuffix(i);
            pos = ion.getMZ(charge, Residue::YIon);
          }
          else if (*ft_it == "z")
          {
            ion = sequence.getSuffix(i);
            pos = ion.getMZ(charge, Residue::ZIon);
          }
          else
          {
            throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                *ft_it + " ion series for peptide sequence \"" + sequence.toString() +
                "\" with precursor charge +" + String(precursor_charge) + " could not be generated.");
          }
          
          ionseries[*ft_it + String(i) + "^" + String(charge)] = Math::roundDecimal(pos, round_decPow);

          for (Size j = 0; j < ion.size(); ++j)
          {
            if (ion[j].hasNeutralLoss())
            {
              for (const auto& lit : ion[j].getLossFormulas())
              {
                if (enable_specific_losses && 
                    lit != H2O &&
                    lit != NH3 &&
                    lit != CN2 &&
                    lit != CNO)
                {
                  ionseries[*ft_it + String(i) + "-" + lit.toString() + "^" + String(charge)] =
                    Math::roundDecimal(pos - lit.getMonoWeight() / charge, round_decPow);
                }
                else if (enable_unspecific_losses && (
                    lit == H2O ||
                    lit == NH3 ||
                    lit == CN2 ||
                    lit == CNO))
                {
                  ionseries[*ft_it + String(i) + "-" + lit.toString() + "^" + String(charge)] =
                    Math::roundDecimal(pos - lit.getMonoWeight() / charge, round_decPow);
                }
              }
            }
          }
        }
      }
    }

    return ionseries;
  }

}
