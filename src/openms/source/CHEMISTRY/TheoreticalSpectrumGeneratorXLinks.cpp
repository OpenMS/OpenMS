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
// $Maintainer: Eugen Netz $
// $Authors: Eugen Netz $
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGeneratorXLinks.h>
#include <OpenMS/DATASTRUCTURES/Map.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/CHEMISTRY/IsotopeDistribution.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <limits.h>

using namespace std;

namespace OpenMS
{

  TheoreticalSpectrumGeneratorXLinks::TheoreticalSpectrumGeneratorXLinks() :
    DefaultParamHandler("TheoreticalSpectrumGeneratorXLinks")
  {
    defaults_.setValue("add_isotopes", "false", "If set to 1 isotope peaks of the product ion peaks are added");
    defaults_.setValidStrings("add_isotopes", ListUtils::create<String>("true,false"));

    defaults_.setValue("max_isotope", 2, "Defines the maximal isotopic peak which is added, add_isotopes must be set to 1");

    defaults_.setValue("add_metainfo", "true", "Adds the type of peaks as metainfo to the peaks, like y8+, [M-H2O+2H]++");
    defaults_.setValidStrings("add_metainfo", ListUtils::create<String>("true,false"));

    defaults_.setValue("add_losses", "false", "Adds common losses to those ion expect to have them, only water and ammonia loss is considered");
    defaults_.setValidStrings("add_losses", ListUtils::create<String>("true,false"));

    defaults_.setValue("add_precursor_peaks", "false", "Adds peaks of the precursor to the spectrum, which happen to occur sometimes");
    defaults_.setValidStrings("add_precursor_peaks", ListUtils::create<String>("true,false"));

    defaults_.setValue("add_abundant_immonium_ions", "false", "Add most abundant immonium ions");
    defaults_.setValidStrings("add_abundant_immonium_ions", ListUtils::create<String>("true,false"));

    defaults_.setValue("add_first_prefix_ion", "true", "If set to true e.g. b1 ions are added");
    defaults_.setValidStrings("add_first_prefix_ion", ListUtils::create<String>("true,false"));

    defaults_.setValue("add_y_ions", "true", "Add peaks of y-ions to the spectrum");
    defaults_.setValidStrings("add_y_ions", ListUtils::create<String>("true,false"));

    defaults_.setValue("add_b_ions", "true", "Add peaks of b-ions to the spectrum");
    defaults_.setValidStrings("add_b_ions", ListUtils::create<String>("true,false"));

    defaults_.setValue("add_a_ions", "false", "Add peaks of a-ions to the spectrum");
    defaults_.setValidStrings("add_a_ions", ListUtils::create<String>("true,false"));

    defaults_.setValue("add_c_ions", "false", "Add peaks of c-ions to the spectrum");
    defaults_.setValidStrings("add_c_ions", ListUtils::create<String>("true,false"));

    defaults_.setValue("add_x_ions", "false", "Add peaks of  x-ions to the spectrum");
    defaults_.setValidStrings("add_x_ions", ListUtils::create<String>("true,false"));

    defaults_.setValue("add_z_ions", "false", "Add peaks of z-ions to the spectrum");
    defaults_.setValidStrings("add_z_ions", ListUtils::create<String>("true,false"));


    // intensity options of the ions
    defaults_.setValue("y_intensity", 1.0, "Intensity of the y-ions");
    defaults_.setValue("b_intensity", 1.0, "Intensity of the b-ions");
    defaults_.setValue("a_intensity", 1.0, "Intensity of the a-ions");
    defaults_.setValue("c_intensity", 1.0, "Intensity of the c-ions");
    defaults_.setValue("x_intensity", 1.0, "Intensity of the x-ions");
    defaults_.setValue("z_intensity", 1.0, "Intensity of the z-ions");

    defaults_.setValue("relative_loss_intensity", 0.1, "Intensity of loss ions, in relation to the intact ion intensity");

    // precursor intensity
    defaults_.setValue("precursor_intensity", 1.0, "Intensity of the precursor peak");
    defaults_.setValue("precursor_H2O_intensity", 1.0, "Intensity of the H2O loss peak of the precursor");
    defaults_.setValue("precursor_NH3_intensity", 1.0, "Intensity of the NH3 loss peak of the precursor");

    defaultsToParam_();
  }

  TheoreticalSpectrumGeneratorXLinks::TheoreticalSpectrumGeneratorXLinks(const TheoreticalSpectrumGeneratorXLinks & rhs) :
    DefaultParamHandler(rhs)
  {
  }

  TheoreticalSpectrumGeneratorXLinks & TheoreticalSpectrumGeneratorXLinks::operator=(const TheoreticalSpectrumGeneratorXLinks & rhs)
  {
    if (this != &rhs)
    {
      DefaultParamHandler::operator=(rhs);
    }
    return *this;
  }

  TheoreticalSpectrumGeneratorXLinks::~TheoreticalSpectrumGeneratorXLinks()
  {
  }

  void TheoreticalSpectrumGeneratorXLinks::getCommonIonSpectrum(RichPeakSpectrum & spec, const ProteinProteinCrossLink& cross_link, Int charge, bool fragment_alpha_chain) const
  {

    for (Int z = 1; z <= charge; ++z)
    {
      if (add_b_ions_)
        addCommonPeaks(spec, cross_link, Residue::BIon, z, fragment_alpha_chain);
      if (add_y_ions_)
        addCommonPeaks(spec, cross_link, Residue::YIon, z, fragment_alpha_chain);
    }

    spec.sortByPosition();
  }

  void TheoreticalSpectrumGeneratorXLinks::getXLinkIonSpectrum(RichPeakSpectrum & spec_alpha, RichPeakSpectrum & spec_beta, const ProteinProteinCrossLink& cross_link, Int mincharge, Int maxcharge) const
  {

    for (Int z = mincharge; z <= maxcharge; ++z)
    {
      if (add_b_ions_)
        addXLinkIonPeaks(spec_alpha, spec_beta, cross_link, Residue::BIon, z);
      if (add_y_ions_)
        addXLinkIonPeaks(spec_alpha, spec_beta, cross_link, Residue::YIon, z);
    }

    spec_alpha.sortByPosition();
    spec_beta.sortByPosition();
    return;
  }

  // Function for mono- and loop-links
  void TheoreticalSpectrumGeneratorXLinks::getXLinkIonSpectrum(RichPeakSpectrum & spec_alpha, const ProteinProteinCrossLink& cross_link, Int mincharge, Int maxcharge) const
  {

    for (Int z = mincharge; z <= maxcharge; ++z)
    {
      if (add_b_ions_)
        addXLinkIonPeaks(spec_alpha, cross_link, Residue::BIon, z);
      if (add_y_ions_)
        addXLinkIonPeaks(spec_alpha, cross_link, Residue::YIon, z);
    }

    spec_alpha.sortByPosition();
    return;
  }

  void TheoreticalSpectrumGeneratorXLinks::addXLinkIonPeaks(RichPeakSpectrum & spec_alpha, RichPeakSpectrum & spec_beta, const ProteinProteinCrossLink & cross_link, Residue::ResidueType res_type, Int charge) const
  {
    const AASequence& peptideA = cross_link.alpha;
    const AASequence& peptideB = cross_link.beta;
    double cross_link_mass = cross_link.cross_linker_mass;

    if (peptideA.empty())
    {
      cout << "Warning: Attempt at creating XLink Ions Spectrum from empty string!" << endl;
      return;
    }

    const SignedSize xlink_pos_A = cross_link.cross_link_position.first;
    const SignedSize xlink_pos_B = cross_link.cross_link_position.second;

    Map<double, AASequence> ions_alpha;
    Map<double, AASequence> ions_beta;
    Map<double, String> names;
    AASequence ion;
    double intensity(1.0);
    //bool add_first_prefix_ion(param_.getValue("add_first_prefix_ion").toBool());

//    double xlink_mass = param_.getValue("cross_link_type2_mass");
    double peptideA_mass(peptideA.getMonoWeight());
    double peptideB_mass(0);
    if (xlink_pos_B != -1)
    {
      peptideB_mass = peptideB.getMonoWeight();
    }

    // Debug support output
    /*
    cout << "peptideA: " << peptideA.toString() << endl;
    cout << "peptideAX_gen_String: " << new_peptideA << endl;
    cout << "peptideAX: " << peptideA_xlink.toString() << endl;
    cout << "peptideB: " << peptideB.toString() << endl;
    cout << "peptideBX_gen_String: " << new_peptideB << endl;
    cout << "peptideBX: " << peptideB_xlink.toString() << endl;
    */
    // Generate the ion peaks:
    // Does not generate peaks of full peptide (therefore "<").
    // They are added via precursor mass (and neutral losses).
    // Could be changed in the future.
    switch (res_type)
    {
    case Residue::BIon:
    {
//      Size i = xlink_pos_A+1;
      for (Size i = xlink_pos_A+1; i < peptideA.size(); ++i)
      {
        ion = peptideA.getPrefix(i);
        double pos = (ion.getMonoWeight(Residue::BIon, charge) + (cross_link_mass + peptideB_mass)) / static_cast<double>(charge);
        // Adding a second isotopic peak, as it is the most intense one in many cases for cross-links
        double pos2 = (ion.getMonoWeight(Residue::BIon, charge) + (cross_link_mass + peptideB_mass) + Constants::NEUTRON_MASS_U) / static_cast<double>(charge);
        ions_alpha[pos] = ion;
        ions_alpha[pos2] = ion;
        names[pos] = "[alpha$b" + String(i) + "]"; // + "X" + String(charge, '+');
        names[pos2] = "[alpha$b" + String(i) + "]"; // + "X" + String(charge, '+');
//        cout << "XLink PepA, b-ion: " << ion.toString() << "\t Charge: " << charge << "\t Mono: " << (ion.getMonoWeight(Residue::BIon, charge) + (cross_link_mass + peptideB_mass)) << "\t MZ: " << pos << endl;
      }

      if (xlink_pos_B != -1) {
//      i = xlink_pos_B+1;
        for (Size i = xlink_pos_B+1; i < peptideB.size(); ++i)
        {
          ion = peptideB.getPrefix(i);
          double pos =(ion.getMonoWeight(Residue::BIon, charge) + (cross_link_mass + peptideA_mass)) / static_cast<double>(charge);
          // Adding a second isotopic peak, as it is the most intense one in many cases for cross-links
          double pos2 = (ion.getMonoWeight(Residue::BIon, charge) + (cross_link_mass + peptideA_mass) + Constants::NEUTRON_MASS_U) / static_cast<double>(charge);
          ions_beta[pos] = ion;
          ions_beta[pos2] = ion;
          names[pos] = "[beta$b" + String(i) + "]"; // + "X" + String(charge, '+');
          names[pos2] = "[beta$b" + String(i) + "]"; // + "X" + String(charge, '+');
  //        cout << "XLink PepB, b-ion: " << ion.toString() << "\t Charge: " << charge << "\t Mono: " << (ion.getMonoWeight(Residue::BIon, charge) + (cross_link_mass + peptideA_mass)) << "\t MZ: " << pos << endl;
        }
      }


      intensity = b_intensity_;
      break;
    }

    case Residue::YIon:
    {
      for (Size i = (peptideA.size() - xlink_pos_A);  i < peptideA.size(); ++i)
      {
        ion = peptideA.getSuffix(i);
        double pos = (ion.getMonoWeight(Residue::YIon, charge) + (cross_link_mass + peptideB_mass)) / static_cast<double>(charge);
        // Adding a second isotopic peak, as it is the most intense one in many cases for cross-links
        double pos2 = (ion.getMonoWeight(Residue::BIon, charge) + (cross_link_mass + peptideB_mass) + Constants::NEUTRON_MASS_U) / static_cast<double>(charge);
        ions_alpha[pos] = ion;
        ions_alpha[pos2] = ion;
        names[pos] = "[alpha$y" + String(i) + "]"; // + "X" + String(charge, '+');
        names[pos2] = "[alpha$y" + String(i) + "]"; // + "X" + String(charge, '+');
//        cout << "XLink PepA, y-ion: " << ion.toString() << "\t Charge: " << charge << "\t Mono: " << (ion.getMonoWeight(Residue::BIon, charge) + (cross_link_mass + peptideB_mass)) << "\t MZ: " << pos << endl;
      }

      if (xlink_pos_B != -1)
      {
        for (Size i = (peptideB.size() - xlink_pos_B);  i < peptideB.size(); ++i)
        {
          ion = peptideB.getSuffix(i);
          double pos = (ion.getMonoWeight(Residue::YIon, charge) + (cross_link_mass + peptideA_mass)) / static_cast<double>(charge);
          // Adding a second isotopic peak, as it is the most intense one in many cases for cross-links
          double pos2 = (ion.getMonoWeight(Residue::BIon, charge) + (cross_link_mass + peptideA_mass) + Constants::NEUTRON_MASS_U) / static_cast<double>(charge);
          ions_beta[pos] = ion;
          ions_beta[pos2] = ion;
          names[pos] = "[beta$y" + String(i) + "]"; // + "X" + String(charge, '+');
          names[pos2] = "[beta$y" + String(i) + "]"; // + "X" + String(charge, '+');
  //        cout << "XLink PepB, y-ion: " << ion.toString() << "\t Charge: " << charge << "\t Mono: " << (ion.getMonoWeight(Residue::BIon, charge) + (cross_link_mass + peptideA_mass)) << "\t MZ: " << pos << endl;
        }
      }

      intensity = y_intensity_;
      break;
    }

    default:
      cerr << "Cannot create peaks of that ion type" << endl;
    }

    RichPeak1D p;
    for (Map<double, AASequence>::ConstIterator cit = ions_alpha.begin(); cit != ions_alpha.end(); ++cit)
    {
      ion = cit->second;
      double pos = cit->first;
      String ion_name = names[pos];

        p.setMZ(pos);
        p.setIntensity(intensity);
        if (add_metainfo_)
        {
          p.setMetaValue("IonName", ion_name);
          p.setMetaValue("z", charge);
        }
        spec_alpha.push_back(p);

      if (add_losses_)
      {
        set<String> losses;
        for (AASequence::ConstIterator it = cit->second.begin(); it != cit->second.end(); ++it)
        {
          if (it->hasNeutralLoss())
          {
            vector<EmpiricalFormula> loss_formulas = it->getLossFormulas();
            for (Size i = 0; i != loss_formulas.size(); ++i)
            {
              losses.insert(loss_formulas[i].toString());
            }
          }
        }


        if (!add_isotopes_)
        {
          p.setIntensity(intensity * rel_loss_intensity_);
        }

        for (set<String>::const_iterator it = losses.begin(); it != losses.end(); ++it)
        {
          EmpiricalFormula loss_ion = ion.getFormula(res_type, charge) - EmpiricalFormula(*it);
          // thanks to Chris and Sandro
          // check for negative element frequencies (might happen if losses are not allowed for specific ions)
          bool negative_elements(false);
          for (EmpiricalFormula::ConstIterator eit = loss_ion.begin(); eit != loss_ion.end(); ++eit)
          {
            if (eit->second < 0)
            {
              negative_elements = true;
              break;
            }
          }
          if (negative_elements)
          {
            continue;
          }
          double loss_pos = (loss_ion.getMonoWeight() + (cross_link_mass + peptideB_mass)) / (double)charge;
          String loss_name = *it;

            p.setMZ(loss_pos);
            if (add_metainfo_)
            {
              p.setMetaValue("IonName", ion_name + "-" + loss_name);
              p.setMetaValue("z", charge);
            }
            spec_alpha.push_back(p);

        }
      }
    }

    for (Map<double, AASequence>::ConstIterator cit = ions_beta.begin(); cit != ions_beta.end(); ++cit)
    {
      ion = cit->second;
      double pos = cit->first;
      String ion_name = names[pos];

        p.setMZ(pos);
        p.setIntensity(intensity);
        if (add_metainfo_)
        {
          p.setMetaValue("IonName", ion_name);
          p.setMetaValue("z", charge);
        }
        spec_beta.push_back(p);


      if (add_losses_)
      {
        set<String> losses;
        for (AASequence::ConstIterator it = cit->second.begin(); it != cit->second.end(); ++it)
        {
          if (it->hasNeutralLoss())
          {
            vector<EmpiricalFormula> loss_formulas = it->getLossFormulas();
            for (Size i = 0; i != loss_formulas.size(); ++i)
            {
              losses.insert(loss_formulas[i].toString());
            }
          }
        }


        if (!add_isotopes_)
        {
          p.setIntensity(intensity * rel_loss_intensity_);
        }

        for (set<String>::const_iterator it = losses.begin(); it != losses.end(); ++it)
        {
          EmpiricalFormula loss_ion = ion.getFormula(res_type, charge) - EmpiricalFormula(*it);
          // thanks to Chris and Sandro
          // check for negative element frequencies (might happen if losses are not allowed for specific ions)
          bool negative_elements(false);
          for (EmpiricalFormula::ConstIterator eit = loss_ion.begin(); eit != loss_ion.end(); ++eit)
          {
            if (eit->second < 0)
            {
              negative_elements = true;
              break;
            }
          }
          if (negative_elements)
          {
            continue;
          }
          double loss_pos = (loss_ion.getMonoWeight() + (cross_link_mass + peptideB_mass))/ (double)charge;
          String loss_name = *it;

          if (add_isotopes_)
          {
            IsotopeDistribution dist = loss_ion.getIsotopeDistribution(max_isotope_);
            UInt j(0);
            for (IsotopeDistribution::ConstIterator iso = dist.begin(); iso != dist.end(); ++iso)
            {
              p.setMZ((double)(loss_pos + j) / (double)charge);
              p.setIntensity(intensity * rel_loss_intensity_ * iso->second);
              if (add_metainfo_ && j == 0)
              {
                p.setMetaValue("IonName", ion_name + "-" + loss_name);
              }
              spec_beta.push_back(p);
            }
          }
          else
          {
            p.setMZ(loss_pos);
            if (add_metainfo_)
            {
              p.setMetaValue("IonName", ion_name + "-" + loss_name);
              p.setMetaValue("z", charge);
            }
            spec_beta.push_back(p);
          }
        }
      }
    }

    if (add_metainfo_)
    {
      p.setMetaValue("IonName", String(""));
    }

//    spec_alpha.sortByPosition();
//    spec_beta.sortByPosition();

    return;
  }

  // MONO AND LOOP LINKS
  void TheoreticalSpectrumGeneratorXLinks::addXLinkIonPeaks(RichPeakSpectrum & spec_alpha, const ProteinProteinCrossLink & cross_link, Residue::ResidueType res_type, Int charge) const
  {
    const AASequence& peptideA = cross_link.alpha;
    double cross_link_mass = cross_link.cross_linker_mass;

    if (peptideA.empty())
    {
      cout << "Warning: Attempt at creating XLink Ions Spectrum from empty string!" << endl;
      return;
    }

    SignedSize xlink_pos_A;
    SignedSize xlink_pos_B;
    // Mono-link has only one position, which can be used for b- and y-ions
    if (cross_link.cross_link_position.second == -1)
    {
      xlink_pos_A = cross_link.cross_link_position.first;
      xlink_pos_B = cross_link.cross_link_position.first;
    }
    // Loop-link has two different positions, the smaller of the two has to be used for b-ions and the larger for y-ions
    else
    {
      if (cross_link.cross_link_position.first > cross_link.cross_link_position.second)
      {
        xlink_pos_A = cross_link.cross_link_position.first;
        xlink_pos_B = cross_link.cross_link_position.second;
      }
      else
      {
        xlink_pos_A = cross_link.cross_link_position.second;
        xlink_pos_B = cross_link.cross_link_position.first;
      }
    }


    Map<double, AASequence> ions_alpha;
    Map<double, String> names;
    AASequence ion;
    double intensity(1.0);
    //bool add_first_prefix_ion(param_.getValue("add_first_prefix_ion").toBool());

    // Debug support output
    /*
    cout << "peptideA: " << peptideA.toString() << endl;
    cout << "peptideAX_gen_String: " << new_peptideA << endl;
    cout << "peptideAX: " << peptideA_xlink.toString() << endl;
    cout << "peptideB: " << peptideB.toString() << endl;
    cout << "peptideBX_gen_String: " << new_peptideB << endl;
    cout << "peptideBX: " << peptideB_xlink.toString() << endl;
    */
    // Generate the ion peaks:
    // Does not generate peaks of full peptide (therefore "<").
    // They are added via precursor mass (and neutral losses).
    // Could be changed in the future.
    switch (res_type)
    {
    case Residue::BIon:
    {
      Size i = xlink_pos_A+1;
      for (; i < peptideA.size(); ++i)
      {
        ion = peptideA.getPrefix(i);
        double pos = (ion.getMonoWeight(Residue::BIon, charge) + cross_link_mass) / static_cast<double>(charge);
        // Adding a second isotopic peak, as it is the most intense one in many cases for cross-links
        double pos2 = (ion.getMonoWeight(Residue::BIon, charge) + cross_link_mass + Constants::NEUTRON_MASS_U) / static_cast<double>(charge);
        ions_alpha[pos] = ion;
        ions_alpha[pos2] = ion;
        names[pos] = "[alpha$b" + String(i) + "]"; // + "X" + String(charge, '+');
        names[pos2] = names[pos];
//        cout << "XLink PepA, b-ion: " << ion.toString() << "\t Charge: " << charge << "\t Mono: " << (ion.getMonoWeight(Residue::BIon, charge) + (cross_link_mass + peptideB_mass)) << "\t MZ: " << pos << endl;
      }

      intensity = b_intensity_;
      break;
    }

    case Residue::YIon:
    {
      for (Size i = (peptideA.size() - xlink_pos_B);  i < peptideA.size(); ++i)
      {
        ion = peptideA.getSuffix(i);
        double pos = (ion.getMonoWeight(Residue::YIon, charge) + cross_link_mass) / static_cast<double>(charge);
        // Adding a second isotopic peak, as it is the most intense one in many cases for cross-links
        double pos2 = (ion.getMonoWeight(Residue::BIon, charge) + cross_link_mass + Constants::NEUTRON_MASS_U) / static_cast<double>(charge);
        ions_alpha[pos] = ion;
        ions_alpha[pos2] = ion;
        names[pos] = "[alpha$y" + String(i) + "]"; // + "X" + String(charge, '+');
        names[pos2] = names[pos];
//        cout << "XLink PepA, y-ion: " << ion.toString() << "\t Charge: " << charge << "\t Mono: " << (ion.getMonoWeight(Residue::BIon, charge) + (cross_link_mass + peptideB_mass)) << "\t MZ: " << pos << endl;
      }

      intensity = y_intensity_;
      break;
    }

    default:
      cerr << "Cannot create peaks of that ion type" << endl;
    }

    RichPeak1D p;
    for (Map<double, AASequence>::ConstIterator cit = ions_alpha.begin(); cit != ions_alpha.end(); ++cit)
    {
      ion = cit->second;
      double pos = cit->first;
      String ion_name = names[pos];

        p.setMZ(pos);
        p.setIntensity(intensity);
        if (add_metainfo_)
        {
          p.setMetaValue("IonName", ion_name);
          p.setMetaValue("z", charge);
        }
        spec_alpha.push_back(p);


      if (add_losses_)
      {
        set<String> losses;
        for (AASequence::ConstIterator it = cit->second.begin(); it != cit->second.end(); ++it)
        {
          if (it->hasNeutralLoss())
          {
            vector<EmpiricalFormula> loss_formulas = it->getLossFormulas();
            for (Size i = 0; i != loss_formulas.size(); ++i)
            {
              losses.insert(loss_formulas[i].toString());
            }
          }
        }


        if (!add_isotopes_)
        {
          p.setIntensity(intensity * rel_loss_intensity_);
        }

        for (set<String>::const_iterator it = losses.begin(); it != losses.end(); ++it)
        {
          EmpiricalFormula loss_ion = ion.getFormula(res_type, charge) - EmpiricalFormula(*it);
          // thanks to Chris and Sandro
          // check for negative element frequencies (might happen if losses are not allowed for specific ions)
          bool negative_elements(false);
          for (EmpiricalFormula::ConstIterator eit = loss_ion.begin(); eit != loss_ion.end(); ++eit)
          {
            if (eit->second < 0)
            {
              negative_elements = true;
              break;
            }
          }
          if (negative_elements)
          {
            continue;
          }
          double loss_pos = (loss_ion.getMonoWeight() + cross_link_mass) / (double)charge;
          String loss_name = *it;

            p.setMZ(loss_pos);
            if (add_metainfo_)
            {
              p.setMetaValue("IonName", ion_name + "-" + loss_name);
              p.setMetaValue("z", charge);
            }
            spec_alpha.push_back(p);

        }
      }
    }

    if (add_metainfo_)
    {
      p.setMetaValue("IonName", String(""));
    }

//    spec_alpha.sortByPosition();

    return;
  }

  void TheoreticalSpectrumGeneratorXLinks::addAbundantImmoniumIons(RichPeakSpectrum & spec, const AASequence& peptide) const
  {
    //bool add_metainfo(param_.getValue("add_metainfo").toBool());

    RichPeak1D p;

    // just in case someone wants the ion names;
    p.metaRegistry().registerName("IonName", "Name of the ion");

    // Histidin immonium ion (C5H8N3)
    if (peptide.has(*ResidueDB::getInstance()->getResidue('H')))
    {
      p.setMZ(110.0718);
      p.setIntensity(1.0);
      if (add_metainfo_)
      {
        String name("iH");
        p.setMetaValue("IonName", name);
      }
      spec.push_back(p);
    }

    // Phenylalanin immonium ion (C8H10N)
    if (peptide.has(*ResidueDB::getInstance()->getResidue('F')))
    {
      p.setMZ(120.0813);
      p.setIntensity(1.0);
      if (add_metainfo_)
      {
        String name("iF");
        p.setMetaValue("IonName", name);
      }
      spec.push_back(p);
    }

    // Tyrosine immonium ion (C8H10NO)
    if (peptide.has(*ResidueDB::getInstance()->getResidue('Y')))
    {
      p.setMZ(136.0762);
      p.setIntensity(1.0);
      if (add_metainfo_)
      {
        String name("iY");
        p.setMetaValue("IonName", name);
      }
      spec.push_back(p);
    }

    // Iso/Leucin immonium ion (same mass for immonium ion)
    if (peptide.has(*ResidueDB::getInstance()->getResidue('L')))
    {
      p.setMZ(86.09698);
      p.setIntensity(1.0);
      if (add_metainfo_)
      {
        String name("iL/I");
        p.setMetaValue("IonName", name);
      }
      spec.push_back(p);
    }

    // Tryptophan immonium ion
    if (peptide.has(*ResidueDB::getInstance()->getResidue('W')))
    {
      p.setMZ(159.0922);
      p.setIntensity(1.0);
      if (add_metainfo_)
      {
        String name("iW");
        p.setMetaValue("IonName", name);
      }
      spec.push_back(p);
    }

    // Cysteine (C2H6NS)
    if (peptide.has(*ResidueDB::getInstance()->getResidue('C')))
    {
      p.setMZ(76.0221);
      p.setIntensity(1.0);
      if (add_metainfo_)
      {
        String name("iC");
        p.setMetaValue("IonName", name);
      }
      spec.push_back(p);
    }

    // Proline immonium ion (C4H8N)
    if (peptide.has(*ResidueDB::getInstance()->getResidue('P')))
    {
      p.setMZ(70.0656);
      p.setIntensity(1.0);
      if (add_metainfo_)
      {
        String name("iP");
        p.setMetaValue("IonName", name);
      }
      spec.push_back(p);
    }

    spec.sortByPosition();
  }

  void TheoreticalSpectrumGeneratorXLinks::addCommonPeaks(RichPeakSpectrum & spectrum, const ProteinProteinCrossLink & cross_link, Residue::ResidueType res_type, Int charge, bool fragment_alpha_chain) const
  {
    AASequence peptide;
    SignedSize xlink_pos_A;
    SignedSize xlink_pos_B;
    if (fragment_alpha_chain)
    {
      peptide = cross_link.alpha;
      xlink_pos_A = cross_link.cross_link_position.first;
      if (cross_link.cross_link_position.second == -1)
      {
        // Ions of alpha chain without second position: mono-link
        xlink_pos_B = cross_link.cross_link_position.first;
      }
      else
      {
        // Ions of alpha chain with second position: could be a cross-link or a loop-link
        if (cross_link.beta.size() > 0)
        {
          // Cross-link, only first position is on alpha chain
          xlink_pos_B = cross_link.cross_link_position.first;
        }
        else
        {
          // Loop-link, both positions are on alpha chain
          xlink_pos_B = cross_link.cross_link_position.second;
        }
      }
    }
    else
    {
      // Ions of beta chain, but beta is empty, or has no position for a cross-link, should never happen
      if (cross_link.cross_link_position.second == -1 || (cross_link.beta.size() == 0 ))
      {
        cout << "Warning: Attempt at creating Common Ions Spectrum from Beta chain without sequence or second cross-link position!" << endl;
        return;
      }
      // Ions of beta chain, if beta chain exists this is a cross-link, only second position is on beta chain
      peptide = cross_link.beta;
      xlink_pos_A = cross_link.cross_link_position.second;
      xlink_pos_B = cross_link.cross_link_position.second;
    }


    if (peptide.empty())
    {
      cout << "Warning: Attempt at creating Common Ions Spectrum from empty string!" << endl;
      return;
    }

    Map<double, AASequence> ions;
    Map<double, String> names;
    AASequence ion;
    double intensity(1.0);
    //bool add_first_prefix_ion(param_.getValue("add_first_prefix_ion").toBool());

    //cout << "CommonIons, peptide: " << peptide.toString() << endl;

    // Generate the ion peaks:
    // Does not generate peaks of full peptide (therefore "<").
    // They are added via precursor mass (and neutral losses).
    // Could be changed in the future.
    switch (res_type)
    {
    case Residue::BIon:
    {
      Size i = 1;
      if (!add_first_prefix_ion_)
      {
        i = 2;
      }
      for (; i < xlink_pos_A+1; ++i)
      {
        ion = peptide.getPrefix(i);
        double pos = ion.getMonoWeight(Residue::BIon, charge) / static_cast<double>(charge);
        ions[pos] = ion;
        if (fragment_alpha_chain)
        {
          names[pos] = "[alpha$b" + String(i) + "]"; // + String(charge, '+');
        }
        else
        {
          names[pos] = "[beta$b" + String(i) + "]"; // + String(charge, '+');
        }

//        cout << "CommonIons, b-ion: " << ion.toString() << "\t Charge: " << charge << "\t Mono: " << ion.getMonoWeight(Residue::BIon, charge) << "\t MZ: " << pos << endl;
      }
      intensity = b_intensity_;
      break;
    }

    case Residue::YIon:
    {
      for (Size i = 1; i < (peptide.size() - xlink_pos_B); ++i)
      {
        ion = peptide.getSuffix(i);
        double pos = ion.getMonoWeight(Residue::YIon, charge) / static_cast<double>(charge);
        ions[pos] = ion;
        if (fragment_alpha_chain)
        {
          names[pos] = "[alpha$y" + String(i) + "]"; // + String(charge, '+');
        }
        else
        {
          names[pos] = "[beta$y" + String(i) + "]"; // + String(charge, '+');
        }
//        cout << "CommonIons, y-ion: " << ion.toString() << "\t Charge: " << charge << "\t Mono: " << ion.getMonoWeight(Residue::YIon, charge) << "\t MZ: " << pos << endl;
      }
      intensity = y_intensity_;
      break;
    }

    default:
      cerr << "Cannot create peaks of that ion type" << endl;
    }

    RichPeak1D p;
    for (Map<double, AASequence>::ConstIterator cit = ions.begin(); cit != ions.end(); ++cit)
    {
      ion = cit->second;
      double pos = cit->first;
      String ion_name = names[pos];
      if (add_isotopes_)
      {
        IsotopeDistribution dist = ion.getFormula(res_type, charge).getIsotopeDistribution(max_isotope_);
        UInt j(0);
        for (IsotopeDistribution::ConstIterator it = dist.begin(); it != dist.end(); ++it, ++j)
        {
          p.setMZ((double)(pos + (double)j * Constants::NEUTRON_MASS_U) / (double)charge);
          p.setIntensity(intensity * it->second);
          if (add_metainfo_ && j == 0)
          {
            p.setMetaValue("IonName", ion_name);
            p.setMetaValue("z", charge);
          }
          spectrum.push_back(p);
        }
      }
      else
      {
        p.setMZ(pos);
        p.setIntensity(intensity);
        if (add_metainfo_)
        {
          p.setMetaValue("IonName", ion_name);
          p.setMetaValue("z", charge);
        }
        spectrum.push_back(p);
      }

      if (add_losses_)
      {
        set<String> losses;
        for (AASequence::ConstIterator it = cit->second.begin(); it != cit->second.end(); ++it)
        {
          if (it->hasNeutralLoss())
          {
            vector<EmpiricalFormula> loss_formulas = it->getLossFormulas();
            for (Size i = 0; i != loss_formulas.size(); ++i)
            {
              losses.insert(loss_formulas[i].toString());
            }
          }
        }


        if (!add_isotopes_)
        {
          p.setIntensity(intensity * rel_loss_intensity_);
        }

        for (set<String>::const_iterator it = losses.begin(); it != losses.end(); ++it)
        {
          EmpiricalFormula loss_ion = ion.getFormula(res_type, charge) - EmpiricalFormula(*it);
          // thanks to Chris and Sandro
          // check for negative element frequencies (might happen if losses are not allowed for specific ions)
          bool negative_elements(false);
          for (EmpiricalFormula::ConstIterator eit = loss_ion.begin(); eit != loss_ion.end(); ++eit)
          {
            if (eit->second < 0)
            {
              negative_elements = true;
              break;
            }
          }
          if (negative_elements)
          {
            continue;
          }
          double loss_pos = loss_ion.getMonoWeight() / (double)charge;
          String loss_name = *it;

          if (add_isotopes_)
          {
            IsotopeDistribution dist = loss_ion.getIsotopeDistribution(max_isotope_);
            UInt j(0);
            for (IsotopeDistribution::ConstIterator iso = dist.begin(); iso != dist.end(); ++iso)
            {
              p.setMZ((double)(loss_pos + j) / (double)charge);
              p.setIntensity(intensity * rel_loss_intensity_ * iso->second);
              if (add_metainfo_ && j == 0)
              {
                p.setMetaValue("IonName", ion_name + "-" + loss_name);
              }
              spectrum.push_back(p);
            }
          }
          else
          {
            p.setMZ(loss_pos);
            if (add_metainfo_)
            {
              p.setMetaValue("IonName", ion_name + "-" + loss_name);
            }
            spectrum.push_back(p);
          }
        }
      }
    }

    if (add_metainfo_)
    {
      p.setMetaValue("IonName", String(""));
    }

    spectrum.sortByPosition();

    return;
  }

  void TheoreticalSpectrumGeneratorXLinks::addPrecursorPeaks(RichPeakSpectrum & spec, const AASequence & peptide, Int charge) const
  {

    RichPeak1D p;

    // precursor peak
    double mono_pos = peptide.getMonoWeight(Residue::Full, charge) / double(charge);
    if (add_isotopes_)
    {
      IsotopeDistribution dist = peptide.getFormula(Residue::Full, charge).getIsotopeDistribution(max_isotope_);
      UInt j(0);
      for (IsotopeDistribution::ConstIterator it = dist.begin(); it != dist.end(); ++it, ++j)
      {
        p.setMZ((double)(mono_pos + j * Constants::NEUTRON_MASS_U) / (double)charge);
        p.setIntensity(pre_int_ *  it->second);
        if (add_metainfo_)
        {
          String name("[M+H]+");
          if (charge == 2)
          {
            name = "[M+2H]++";
          }
          p.setMetaValue("IonName", name);
        }
        spec.push_back(p);
      }
    }
    else
    {
      p.setMZ(mono_pos);
      p.setIntensity(pre_int_);
      if (add_metainfo_)
      {
        String name("[M+H]+");
        if (charge == 2)
        {
          name = "[M+2H]++";
        }
        p.setMetaValue("IonName", name);
      }
      spec.push_back(p);
    }
    // loss peaks of the precursor

    //loss of water
    EmpiricalFormula ion = peptide.getFormula(Residue::Full, charge) - EmpiricalFormula("H2O");
    mono_pos = ion.getMonoWeight() / double(charge);
    if (add_isotopes_)
    {
      IsotopeDistribution dist = ion.getIsotopeDistribution(max_isotope_);
      UInt j(0);
      for (IsotopeDistribution::ConstIterator it = dist.begin(); it != dist.end(); ++it, ++j)
      {
        p.setMZ((double)(mono_pos + j * Constants::NEUTRON_MASS_U) / (double)charge);
        p.setIntensity(pre_int_H2O_ *  it->second);
        if (add_metainfo_)
        {
          String name("[M+H]-H2O+");
          if (charge == 2)
          {
            name = "[M+2H]-H2O++";
          }
          p.setMetaValue("IonName", name);
        }
        spec.push_back(p);
      }
    }
    else
    {
      p.setMZ(mono_pos);
      p.setIntensity(pre_int_H2O_);
      if (add_metainfo_)
      {
        String name("[M+H]-H2O+");
        if (charge == 2)
        {
          name = "[M+2H]-H2O++";
        }
        p.setMetaValue("IonName", name);
      }
      spec.push_back(p);
    }

    //loss of ammonia
    ion = peptide.getFormula(Residue::Full, charge) - EmpiricalFormula("NH3");
    mono_pos = ion.getMonoWeight() / double(charge);
    if (add_isotopes_)
    {
      IsotopeDistribution dist = ion.getIsotopeDistribution(max_isotope_);
      UInt j(0);
      for (IsotopeDistribution::ConstIterator it = dist.begin(); it != dist.end(); ++it, ++j)
      {
        p.setMZ((double)(mono_pos + j * Constants::NEUTRON_MASS_U) / (double)charge);
        p.setIntensity(pre_int_NH3_ *  it->second);
        if (add_metainfo_)
        {
          String name("[M+H]-NH3+");
          if (charge == 2)
          {
            name = "[M+2H]-NH3++";
          }
          p.setMetaValue("IonName", name);
        }
        spec.push_back(p);
      }
    }
    else
    {
      p.setMZ(mono_pos);
      p.setIntensity(pre_int_NH3_);
      if (add_metainfo_)
      {
        String name("[M+H]-NH3+");
        if (charge == 2)
        {
          name = "[M+2H]-NH3++";
        }
        p.setMetaValue("IonName", name);
      }
      spec.push_back(p);
    }

    spec.sortByPosition();
  }

  void TheoreticalSpectrumGeneratorXLinks::updateMembers_()
  {
    add_b_ions_ = param_.getValue("add_b_ions").toBool();
    add_y_ions_ = param_.getValue("add_y_ions").toBool();
    add_a_ions_ = param_.getValue("add_a_ions").toBool();
    add_c_ions_ = param_.getValue("add_c_ions").toBool();
    add_x_ions_ = param_.getValue("add_x_ions").toBool();
    add_z_ions_ = param_.getValue("add_z_ions").toBool();
    add_first_prefix_ion_ = param_.getValue("add_first_prefix_ion").toBool();
    add_losses_ = param_.getValue("add_losses").toBool();
    add_metainfo_ = param_.getValue("add_metainfo").toBool();
    add_isotopes_ = param_.getValue("add_isotopes").toBool();
    add_precursor_peaks_ = param_.getValue("add_precursor_peaks").toBool();
    add_abundant_immonium_ions_ = param_.getValue("add_abundant_immonium_ions").toBool();
    a_intensity_ = static_cast<double>(param_.getValue("a_intensity"));
    b_intensity_ = static_cast<double>(param_.getValue("b_intensity"));
    c_intensity_ = static_cast<double>(param_.getValue("c_intensity"));
    x_intensity_ = static_cast<double>(param_.getValue("x_intensity"));
    y_intensity_ = static_cast<double>(param_.getValue("y_intensity"));
    z_intensity_ = static_cast<double>(param_.getValue("z_intensity"));
    max_isotope_ = static_cast<Int>(param_.getValue("max_isotope"));
    rel_loss_intensity_ = static_cast<double>(param_.getValue("relative_loss_intensity"));
    pre_int_ = static_cast<double>(param_.getValue("precursor_intensity"));
    pre_int_H2O_ = static_cast<double>(param_.getValue("precursor_H2O_intensity"));
    pre_int_NH3_ = static_cast<double>(param_.getValue("precursor_NH3_intensity"));
  }

}
