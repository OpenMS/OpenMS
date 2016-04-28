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
// $Maintainer: George Rosenberger $
// $Authors: George Rosenberger $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/MRMIonSeries.h>

namespace OpenMS
{
  MRMIonSeries::MRMIonSeries()
  {
  }

  MRMIonSeries::~MRMIonSeries()
  {
  }

  std::pair<String, double> MRMIonSeries::getIon(IonSeries ionseries, String ionid)
  {
    std::pair<String, double> ion = make_pair(String("unannotated"), -1);

    if (ionseries.find(ionid) != ionseries.end())
    {
      ion = make_pair(ionid, ionseries[ionid]);
    }

    return ion;
  }

  std::pair<String, double> MRMIonSeries::annotateIon(IonSeries ionseries, double ProductMZ, double mz_threshold)
  {
    // make sure to only use annotated transitions and to use the theoretical MZ
    using namespace boost::assign;

    // Iterate over ion type and then ordinal
    std::pair<String, double> ion;
    String unannotated = "unannotated";
    ion = make_pair(unannotated, -1);
    double closest_delta = std::numeric_limits<double>::max();

    for (boost::unordered_map<String, double>::iterator ordinal = ionseries.begin(); ordinal != ionseries.end(); ++ordinal)
    {
      if (std::fabs(ordinal->second - ProductMZ) <= mz_threshold && std::fabs(ordinal->second - ProductMZ) <= closest_delta)
      {
        closest_delta = std::fabs(ordinal->second - ProductMZ);
        ion = make_pair(ordinal->first, ordinal->second);
      }
    }

    return ion;
  }

  CVTermList MRMIonSeries::annotationToCVTermList_(String annotation)
  {
    CVTermList interpretation;

    String fragment_type;
    int fragment_nr = -1;
    double fragment_loss = 0;
    // int fragment_gain = 0;

    std::vector<String> best_annotation;
    annotation.split("/", best_annotation);

    if (best_annotation[0].find("-") != std::string::npos)
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
      CVTerm rank;
      rank.setCVIdentifierRef("MS");
      rank.setAccession("MS:1000926");
      rank.setName("product interpretation rank");
      rank.setValue(1); // we only store the best interpretation
      interpretation.addCVTerm(rank);
    }

    if (fragment_nr != -1)
    {
      CVTerm frag_nr;
      frag_nr.setCVIdentifierRef("MS");
      frag_nr.setAccession("MS:1000903");
      frag_nr.setName("product ion series ordinal");
      frag_nr.setValue(fragment_nr);
      interpretation.addCVTerm(frag_nr);
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
      CVTerm ion;
      ion.setCVIdentifierRef("MS");
      ion.setAccession("MS:1001228");
      ion.setName("frag: x ion");
      interpretation.addCVTerm(ion);
    }
    else if (fragment_type == "y")
    {
      CVTerm ion;
      ion.setCVIdentifierRef("MS");
      ion.setAccession("MS:1001220");
      ion.setName("frag: y ion");
      interpretation.addCVTerm(ion);
    }
    else if (fragment_type == "z")
    {
      CVTerm ion;
      ion.setCVIdentifierRef("MS");
      ion.setAccession("MS:1001230");
      ion.setName("frag: z ion");
      interpretation.addCVTerm(ion);
    }
    else if (fragment_type == "a")
    {
      CVTerm ion;
      ion.setCVIdentifierRef("MS");
      ion.setAccession("MS:1001229");
      ion.setName("frag: a ion");
      interpretation.addCVTerm(ion);
    }
    else if (fragment_type == "b")
    {
      CVTerm ion;
      ion.setCVIdentifierRef("MS");
      ion.setAccession("MS:1001224");
      ion.setName("frag: b ion");
      interpretation.addCVTerm(ion);
    }
    else if (fragment_type == "c")
    {
      CVTerm ion;
      ion.setCVIdentifierRef("MS");
      ion.setAccession("MS:1001231");
      ion.setName("frag: c ion");
      interpretation.addCVTerm(ion);
    }
    else
    {
      CVTerm ion;
      ion.setCVIdentifierRef("MS");
      ion.setAccession("MS:1001240");
      ion.setName("non-identified ion");
      interpretation.addCVTerm(ion);
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

    CVTermList interpretation = annotationToCVTermList_(annotation);

    p.resetInterpretations();
    p.addInterpretation(interpretation);
    tr.setProduct(p);
  }

  void MRMIonSeries::annotateTransitionCV(ReactionMonitoringTransition& tr, String annotation)
  {
    tr.setMetaValue("annotation", annotation);
    annotationToCV_(tr);
  }

  void MRMIonSeries::annotateTransition(ReactionMonitoringTransition& tr, const TargetedExperiment::Peptide peptide, const double precursor_mz_threshold, double product_mz_threshold, bool enable_reannotation, std::vector<String> fragment_types, std::vector<size_t> fragment_charges, bool enable_specific_losses, bool enable_unspecific_losses, int round_decPow)
  {
    CVTermList interpretation;

    OpenMS::AASequence sequence = TargetedExperimentHelper::getAASequence(peptide);
    double prec_pos = sequence.getMonoWeight(Residue::Full, peptide.getChargeState()) / peptide.getChargeState();
    bool unannotated = false;
    std::pair<String, double> target_ion = std::make_pair(String("unannotated"), -1);
    double pos = -1;
    String ionstring;

    // product ion series ordinal
    if (tr.getProduct().getInterpretationList().size() > 0)
    {
      interpretation = tr.getProduct().getInterpretationList()[0];
      AASequence ion;

      if (interpretation.hasCVTerm("MS:1000903")) // if ordinal is set
      {
        int ordinal = interpretation.getCVTerms()["MS:1000903"][0].getValue().toString().toInt();

        if (interpretation.hasCVTerm("MS:1001228"))
        {
          ion = sequence.getSuffix(ordinal);
          ionstring += "x";
          pos = ion.getMonoWeight(Residue::XIon, tr.getProduct().getChargeState()) / (double) tr.getProduct().getChargeState();
        }
        else if (interpretation.hasCVTerm("MS:1001220"))
        {
          ion = sequence.getSuffix(ordinal);
          ionstring += "y";
          pos = ion.getMonoWeight(Residue::YIon, tr.getProduct().getChargeState()) / (double) tr.getProduct().getChargeState();
        }
        else if (interpretation.hasCVTerm("MS:1001230"))
        {
          ion = sequence.getSuffix(ordinal);
          ionstring += "z";
          pos = ion.getMonoWeight(Residue::ZIon, tr.getProduct().getChargeState()) / (double) tr.getProduct().getChargeState();
        }
        else if (interpretation.hasCVTerm("MS:1001229"))
        {
          ion = sequence.getSuffix(ordinal);
          ionstring += "a";
          pos = ion.getMonoWeight(Residue::AIon, tr.getProduct().getChargeState()) / (double) tr.getProduct().getChargeState();
        }
        else if (interpretation.hasCVTerm("MS:1001224"))
        {
          ion = sequence.getSuffix(ordinal);
          ionstring += "b";
          pos = ion.getMonoWeight(Residue::BIon, tr.getProduct().getChargeState()) / (double) tr.getProduct().getChargeState();
        }
        else if (interpretation.hasCVTerm("MS:1001231"))
        {
          ion = sequence.getSuffix(ordinal);
          ionstring += "c";
          pos = ion.getMonoWeight(Residue::CIon, tr.getProduct().getChargeState()) / (double) tr.getProduct().getChargeState();
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

      if (interpretation.hasCVTerm("MS:1000903"))
      {
        // product ion series ordinal
        ionstring += interpretation.getCVTerms()["MS:1000903"][0].getValue().toString();
      }
      else
      {
        unannotated = true;
      }

      if (interpretation.hasCVTerm("MS:1001524") && (enable_specific_losses || enable_unspecific_losses)) // fragment ion neutral loss
      {
        double nl = interpretation.getCVTerms()["MS:1001524"][0].getValue().toString().toDouble();
        // SpectraST style neutral losses
        if (nl == -18)
        {
          ionstring += "-H2O1";
          static const EmpiricalFormula neutralloss_h2o("H2O1"); // -18 H2O loss
          pos -= (neutralloss_h2o.getMonoWeight() / (tr.getProduct().getChargeState()));
        }
        else if (nl == -17)
        {
          ionstring += "-H3N1";
          static const EmpiricalFormula neutralloss_nh3("H3N1"); // -17 NH3 loss
          pos -= (neutralloss_nh3.getMonoWeight() / (tr.getProduct().getChargeState()));
        }
        else if (nl == -36)
        {
          ionstring += "-H4O2";
          static const EmpiricalFormula neutralloss_h2oh2o("H4O2"); // -36 2 * H2O loss
          pos -= (neutralloss_h2oh2o.getMonoWeight() / (tr.getProduct().getChargeState()));
        }
        else if (nl == -34)
        {
          ionstring += "-H6N2";
          static const EmpiricalFormula neutralloss_nh3nh3("H6N2"); // -34 2 * NH3 loss
          pos -= (neutralloss_nh3nh3.getMonoWeight() / (tr.getProduct().getChargeState()));
        }
        else if (nl == -35)
        {
          ionstring += "-H5N1O1";
          static const EmpiricalFormula neutralloss_h2onh3("H5N1O1"); // -35 H2O & NH3 loss
          pos -= (neutralloss_h2onh3.getMonoWeight() / (tr.getProduct().getChargeState()));
        }
        else if (nl == -64)
        {
          ionstring += "-C1H4O1S1";
          static const EmpiricalFormula neutralloss_ch4so("C1H4O1S1"); // -64 CH4SO loss
          pos -= (neutralloss_ch4so.getMonoWeight() / (tr.getProduct().getChargeState()));
        }
        else if (nl == -80)
        {
          ionstring += "-H1O3P1";
          static const EmpiricalFormula neutralloss_hpo3("H1O3P1"); // -80 HPO3 loss
          pos -= (neutralloss_hpo3.getMonoWeight() / (tr.getProduct().getChargeState()));
        }
        else if (nl == -98)
        {
          ionstring += "-H3O4P1";
          static const EmpiricalFormula neutralloss_hpo3h2o("H3O4P1"); // -98 HPO3 & H2O loss
          pos -= (neutralloss_hpo3h2o.getMonoWeight() / (tr.getProduct().getChargeState()));
        }
        else if (nl == -45)
        {
          ionstring += "-C1H3N1O1";
          static const EmpiricalFormula neutralloss_ch3no("C1H3N1O1"); // -45 CH3NO loss
          pos -= (neutralloss_ch3no.getMonoWeight() / (tr.getProduct().getChargeState()));
        }
        else if (nl == -44)
        {
          ionstring += "-C1O2";
          static const EmpiricalFormula neutralloss_co2("C1O2"); // -44 CO2 loss
          pos -= (neutralloss_co2.getMonoWeight() / (tr.getProduct().getChargeState()));
        }
        else if (nl == -46)
        {
          ionstring += "-C1H2O2";
          static const EmpiricalFormula neutralloss_hccoh("C1H2O2"); // -46 HCOOH loss
          pos -= (neutralloss_hccoh.getMonoWeight() / (tr.getProduct().getChargeState()));

        }
        // Double CV term (compatible with PSI CV terms)
        else if (nl < 0)
        {
          ionstring += String(Math::roundDecimal(nl, round_decPow));
          pos -= (nl / (tr.getProduct().getChargeState()));

        }
        else
        {
          unannotated = true;
        }
      }

      if (tr.getProduct().getChargeState() >= 1 && std::find(fragment_charges.begin(), fragment_charges.end(), tr.getProduct().getChargeState()) != fragment_charges.end())
      {
        ionstring += "^" + String(tr.getProduct().getChargeState());
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
      MRMIonSeries::IonSeries ionseries = getIonSeries(sequence, peptide.getChargeState(), fragment_types, fragment_charges, enable_specific_losses, enable_unspecific_losses);
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
      interpretation.addCVTerm(frag_mzdelta);
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
      CVTerm ion;
      ion.setCVIdentifierRef("MS");
      ion.setAccession("MS:1001240");
      ion.setName("non-identified ion");
      interpretation.addCVTerm(ion);
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

  boost::unordered_map<String, double> MRMIonSeries::getIonSeries(AASequence sequence, size_t precursor_charge, std::vector<String> fragment_types, std::vector<size_t> fragment_charges, bool enable_specific_losses, bool enable_unspecific_losses, int round_decPow)
  {
    boost::unordered_map<String, double> ionseries;

    for (std::vector<String>::iterator ft_it = fragment_types.begin(); ft_it != fragment_types.end(); ++ft_it)
    {
      for (std::vector<size_t>::iterator ch_it = fragment_charges.begin(); ch_it != fragment_charges.end(); ++ch_it)
      {
        size_t charge = *ch_it;

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
            pos = ion.getMonoWeight(Residue::AIon, charge) / (double) charge;
          }
          else if (*ft_it == "b")
          {
            ion = sequence.getPrefix(i);
            pos = ion.getMonoWeight(Residue::BIon, charge) / (double) charge;
          }
          else if (*ft_it == "c")
          {
            ion = sequence.getPrefix(i);
            pos = ion.getMonoWeight(Residue::CIon, charge) / (double) charge;
          }
          else if (*ft_it == "x")
          {
            ion = sequence.getSuffix(i);
            pos = ion.getMonoWeight(Residue::XIon, charge) / (double) charge;
          }
          else if (*ft_it == "y")
          {
            ion = sequence.getSuffix(i);
            pos = ion.getMonoWeight(Residue::YIon, charge) / (double) charge;
          }
          else if (*ft_it == "z")
          {
            ion = sequence.getSuffix(i);
            pos = ion.getMonoWeight(Residue::ZIon, charge) / (double) charge;
          }
          else
          {
            throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, *ft_it + " ion series for peptide sequence \"" + sequence.toString() + "\" with precursor charge +" + String(precursor_charge) + " could not be generated.");
          }

          ionseries[*ft_it + String(i) + "^" + String(charge)] = Math::roundDecimal(pos, round_decPow);

          for (Size j = 0; j < ion.size(); ++j)
          {
            if (ion[j].hasNeutralLoss())
            {
              const std::vector<EmpiricalFormula> losses = ion[j].getLossFormulas();
              for (std::vector<EmpiricalFormula>::const_iterator lit = losses.begin(); lit != losses.end(); ++lit)
              {
                if (enable_specific_losses && lit->toString() != String("H2O1") && lit->toString() != String("H3N1") && lit->toString() != String("C1H2N2") && lit->toString() != String("C1H2N1O1"))
                {
                  ionseries[*ft_it + String(i) + "-" + lit->toString() + "^" + String(charge)] = Math::roundDecimal(pos - lit->getMonoWeight() / charge, round_decPow);
                }
                else if (enable_unspecific_losses && (lit->toString() == String("H2O1") || lit->toString() == String("H3N1") || lit->toString() == String("C1H2N2") || lit->toString() == String("C1H2N1O1")))
                {
                  ionseries[*ft_it + String(i) + "-" + lit->toString() + "^" + String(charge)] = Math::roundDecimal(pos - lit->getMonoWeight() / charge, round_decPow);
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
