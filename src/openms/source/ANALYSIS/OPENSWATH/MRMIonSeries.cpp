// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2015.
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

  void MRMIonSeries::annotationToCV_(ReactionMonitoringTransition& tr)
  {
    OpenMS::ReactionMonitoringTransition::Product p = tr.getProduct();
    CVTermList interpretation;

    String fragment_type;
    int fragment_nr = -1;
    int fragment_loss = 0;
    int fragment_gain = 0;

    String annotation = tr.getMetaValue("annotation");

    std::vector<String> best_annotation;
    annotation.split("/", best_annotation);

    if (best_annotation[0].find("^") != std::string::npos)
    {
      std::vector<String> best_annotation_charge;
      best_annotation[0].split("^", best_annotation_charge);
      p.setChargeState(String(best_annotation_charge[1]).toInt());
    }
    else
    {
      p.setChargeState(1);
    }

    if (best_annotation[0].find("-") != std::string::npos)
    {
      std::vector<String> best_annotation_loss;
      best_annotation[0].split("-", best_annotation_loss);

      fragment_type = best_annotation_loss[0].substr(0, 1);
      fragment_nr = best_annotation_loss[0].substr(1).toInt();
      fragment_loss = -1 * String(best_annotation_loss[1]).toInt();
    }
    else if (best_annotation[0].find("+") != std::string::npos)
    {
      std::vector<String> best_annotation_gain;
      best_annotation[0].split("+", best_annotation_gain);

      fragment_type = best_annotation_gain[0].substr(0, 1);
      fragment_nr = best_annotation_gain[0].substr(1).toInt();
      fragment_gain = String(best_annotation_gain[1]).toInt();
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

    p.resetInterpretations();
    p.addInterpretation(interpretation);
    tr.setProduct(p);
  }

  void MRMIonSeries::annotateTransitionCV(ReactionMonitoringTransition& tr, String annotation)
  {
    tr.setMetaValue("annotation", annotation);
    annotationToCV_(tr);
  }

  void MRMIonSeries::annotateTransition(ReactionMonitoringTransition& tr, const TargetedExperiment::Peptide peptide, const double mz_threshold, bool enable_reannotation, std::vector<String> fragment_types, std::vector<size_t> fragment_charges, bool enable_losses)
  {
    CVTermList interpretation;

    OpenMS::AASequence sequence = TargetedExperimentHelper::getAASequence(peptide);
    bool unannotated = false;
    std::pair<String, double> target_ion = std::make_pair(String("unannotated"), -1);
    double pos = -1;
    String ionstring;

    // product ion series ordinal
    if (tr.getProduct().getInterpretationList().size() > 0)
    {
      interpretation = tr.getProduct().getInterpretationList()[0];
      AASequence ion;

      if (interpretation.hasCVTerm("MS:1001228") && interpretation.hasCVTerm("MS:1000903"))
      {
        ion = sequence.getSuffix(interpretation.getCVTerms()["MS:1000903"][0].getValue().toString().toInt());
        ionstring += "x";
        pos = ion.getMonoWeight(Residue::XIon, tr.getProduct().getChargeState()) / (double) tr.getProduct().getChargeState();
      }
      else if (interpretation.hasCVTerm("MS:1001220") && interpretation.hasCVTerm("MS:1000903"))
      {
        ion = sequence.getSuffix(interpretation.getCVTerms()["MS:1000903"][0].getValue().toString().toInt());
        ionstring += "y";
        pos = ion.getMonoWeight(Residue::YIon, tr.getProduct().getChargeState()) / (double) tr.getProduct().getChargeState();
      }
      else if (interpretation.hasCVTerm("MS:1001230") && interpretation.hasCVTerm("MS:1000903"))
      {
        ion = sequence.getSuffix(interpretation.getCVTerms()["MS:1000903"][0].getValue().toString().toInt());
        ionstring += "z";
        pos = ion.getMonoWeight(Residue::ZIon, tr.getProduct().getChargeState()) / (double) tr.getProduct().getChargeState();
      }
      else if (interpretation.hasCVTerm("MS:1001229") && interpretation.hasCVTerm("MS:1000903"))
      {
        ion = sequence.getPrefix(interpretation.getCVTerms()["MS:1000903"][0].getValue().toString().toInt());
        ionstring += "a";
        pos = ion.getMonoWeight(Residue::AIon, tr.getProduct().getChargeState()) / (double) tr.getProduct().getChargeState();
      }
      else if (interpretation.hasCVTerm("MS:1001224") && interpretation.hasCVTerm("MS:1000903"))
      {
        ion = sequence.getPrefix(interpretation.getCVTerms()["MS:1000903"][0].getValue().toString().toInt());
        ionstring += "b";
        pos = ion.getMonoWeight(Residue::BIon, tr.getProduct().getChargeState()) / (double) tr.getProduct().getChargeState();
      }
      else if (interpretation.hasCVTerm("MS:1001231") && interpretation.hasCVTerm("MS:1000903"))
      {
        ion = sequence.getPrefix(interpretation.getCVTerms()["MS:1000903"][0].getValue().toString().toInt());
        ionstring += "c";
        pos = ion.getMonoWeight(Residue::CIon, tr.getProduct().getChargeState()) / (double) tr.getProduct().getChargeState();
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

      if (interpretation.hasCVTerm("MS:1001524") && enable_losses) // fragment ion neutral loss
      {
        double nl = interpretation.getCVTerms()["MS:1001524"][0].getValue().toString().toDouble();
        if (nl == -18)
        {
          ionstring += "-18";
          static const EmpiricalFormula neutralloss_h2o("H2O"); // -18 H2O loss
          pos -= (neutralloss_h2o.getMonoWeight() / (tr.getProduct().getChargeState()));
        }
        else if (nl == -17)
        {
          ionstring += "-17";
          static const EmpiricalFormula neutralloss_nh3("NH3"); // -17 NH3 loss
          pos -= (neutralloss_nh3.getMonoWeight() / (tr.getProduct().getChargeState()));
        }
        else if (nl == -36)
        {
          ionstring += "-36";
          static const EmpiricalFormula neutralloss_h2oh2o("H2OH2O"); // -36 2 * H2O loss
          pos -= (neutralloss_h2oh2o.getMonoWeight() / (tr.getProduct().getChargeState()));
        }
        else if (nl == -34)
        {
          ionstring += "-34";
          static const EmpiricalFormula neutralloss_nh3nh3("NH3NH3"); // -34 2 * NH3 loss
          pos -= (neutralloss_nh3nh3.getMonoWeight() / (tr.getProduct().getChargeState()));
        }
        else if (nl == -35)
        {
          ionstring += "-35";
          static const EmpiricalFormula neutralloss_h2onh3("H2ONH3"); // -35 H2O & NH3 loss
          pos -= (neutralloss_h2onh3.getMonoWeight() / (tr.getProduct().getChargeState()));
        }
        else if (nl == -64)
        {
          ionstring += "-64";
          static const EmpiricalFormula neutralloss_ch4so("CH4SO"); // -64 CH4SO loss
          pos -= (neutralloss_ch4so.getMonoWeight() / (tr.getProduct().getChargeState()));
        }
        else if (nl == -80)
        {
          ionstring += "-80";
          static const EmpiricalFormula neutralloss_hpo3("HPO3"); // -80 HPO3 loss
          pos -= (neutralloss_hpo3.getMonoWeight() / (tr.getProduct().getChargeState()));
        }
        else if (nl == -98)
        {
          ionstring += "-98";
          static const EmpiricalFormula neutralloss_hpo3h2o("HPO3H2O"); // -98 HPO3 loss
          pos -= (neutralloss_hpo3h2o.getMonoWeight() / (tr.getProduct().getChargeState()));
        }
        else if (nl == -45)
        {
          ionstring += "-45";
          static const EmpiricalFormula neutralloss_ch3no("CH3NO"); // -45 CH3NO loss
          pos -= (neutralloss_ch3no.getMonoWeight() / (tr.getProduct().getChargeState()));
        }
        else if (nl == -44)
        {
          ionstring += "-44";
          static const EmpiricalFormula neutralloss_co2("CO2"); // -44 CO2 loss
          pos -= (neutralloss_co2.getMonoWeight() / (tr.getProduct().getChargeState()));
        }
        else if (nl == -46)
        {
          ionstring += "-46";
          static const EmpiricalFormula neutralloss_hccoh("HCOOH"); // -46 HCOOH loss
          pos -= (neutralloss_hccoh.getMonoWeight() / (tr.getProduct().getChargeState()));

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
      MRMIonSeries::IonSeries ionseries = getIonSeries(sequence, peptide.getChargeState(), fragment_types, fragment_charges, enable_losses);
      target_ion = annotateIon(ionseries, tr.getProductMZ(), mz_threshold);
      ionstring = target_ion.first;
      pos = target_ion.second;
      tr.setMetaValue("annotation", ionstring);
      tr.setProductMZ(pos);

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

    if (!unannotated && std::fabs(tr.getProductMZ() - pos) <= mz_threshold)
    {
      CVTerm frag_mzdelta;
      frag_mzdelta.setCVIdentifierRef("MS");
      frag_mzdelta.setAccession("MS:1000904");
      frag_mzdelta.setName("product ion m/z delta");
      frag_mzdelta.setValue(std::fabs(tr.getProductMZ() - pos));
      interpretation.addCVTerm(frag_mzdelta);
      tr.setProductMZ(pos);
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

  boost::unordered_map<String, double> MRMIonSeries::getIonSeries(AASequence sequence, size_t precursor_charge, std::vector<String> fragment_types, std::vector<size_t> fragment_charges, bool enable_losses)
  {
    boost::unordered_map<String, double> ionseries;

    // Neutral losses of all ion series
    static const EmpiricalFormula neutralloss_h2o("H2O"); // -18 H2O loss
    static const EmpiricalFormula neutralloss_nh3("NH3"); // -17 NH3 loss

    static const EmpiricalFormula neutralloss_h2oh2o("H2OH2O"); // -36 2 * H2O loss
    static const EmpiricalFormula neutralloss_nh3nh3("NH3NH3"); // -34 2 * NH3 loss
    static const EmpiricalFormula neutralloss_h2onh3("H2ONH3"); // -35 H2O & NH3 loss

    // Neutral loss (oxidation) of methionine only
    static const EmpiricalFormula neutralloss_ch4so("CH4SO"); // -64 CH4SO loss

    // Neutral losses (phospho) of serine and threonine only
    static const EmpiricalFormula neutralloss_hpo3("HPO3"); // -80 HPO3 loss
    static const EmpiricalFormula neutralloss_hpo3h2o("HPO3H2O"); // -98 HPO3 loss

    // Neutral loss of asparagine and glutamine only
    static const EmpiricalFormula neutralloss_ch3no("CH3NO"); // -45 CH3NO loss

    // Neutral losses of y-ions only
    static const EmpiricalFormula neutralloss_co2("CO2"); // -44 CO2 loss
    static const EmpiricalFormula neutralloss_hccoh("HCOOH"); // -46 HCOOH loss

    for (std::vector<String>::iterator ft_it = fragment_types.begin(); ft_it != fragment_types.end(); ft_it++)
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

          ionseries[*ft_it + String(i) + "^" + String(charge)] = pos;

          if (enable_losses)
          {
            ionseries[*ft_it + String(i) + "-17" + "^" + String(charge)] = pos - neutralloss_nh3.getMonoWeight() / charge;
            ionseries[*ft_it + String(i) + "-18" + "^" + String(charge)] = pos - neutralloss_h2o.getMonoWeight() / charge;
            ionseries[*ft_it + String(i) + "-34" + "^" + String(charge)] = pos - neutralloss_nh3nh3.getMonoWeight() / charge;
            ionseries[*ft_it + String(i) + "-35" + "^" + String(charge)] = pos - neutralloss_h2onh3.getMonoWeight() / charge;
            ionseries[*ft_it + String(i) + "-36" + "^" + String(charge)] = pos - neutralloss_h2oh2o.getMonoWeight() / charge;
            ionseries[*ft_it + String(i) + "-44" + "^" + String(charge)] = pos - neutralloss_co2.getMonoWeight() / charge;
            ionseries[*ft_it + String(i) + "-46" + "^" + String(charge)] = pos - neutralloss_hccoh.getMonoWeight() / charge;
            if (sequence.toString().find("N") != std::string::npos || sequence.toString().find("Q") != std::string::npos)
            // This hack is implemented to enable the annotation of residue specific modifications in the decoy fragments.
            // If the function is used for generic annotation, use ion.toString() instead of sequence.toString().
            {
              ionseries[*ft_it + String(i) + "-45" + "^" + String(charge)] = pos - neutralloss_ch3no.getMonoWeight() / charge;
            }
            if (sequence.toString().find("M(Oxidation)") != std::string::npos)
            // This hack is implemented to enable the annotation of residue specific modifications in the decoy fragments.
            // If the function is used for generic annotation, use ion.toString() instead of sequence.toString().
            {
              ionseries[*ft_it + String(i) + "-64" + "^" + String(charge)] = pos - neutralloss_ch4so.getMonoWeight() / charge;
            }
            if (sequence.toString().find("S(Phospho)") != std::string::npos || sequence.toString().find("T(Phospho)") != std::string::npos)
            // This hack is implemented to enable the annotation of residue specific modifications in the decoy fragments.
            // If the function is used for generic annotation, use ion.toString() instead of sequence.toString().
            {
              ionseries[*ft_it + String(i) + "-80" + "^" + String(charge)] = pos - neutralloss_hpo3.getMonoWeight() / charge;
              ionseries[*ft_it + String(i) + "-98" + "^" + String(charge)] = pos - neutralloss_hpo3h2o.getMonoWeight() / charge;
            }
          }
        }
      }
    }

    return ionseries;
  }

}
