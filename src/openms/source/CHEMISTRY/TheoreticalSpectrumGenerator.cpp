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
// $Maintainer: Sandro Andreotti $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/DATASTRUCTURES/Map.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/CHEMISTRY/IsotopeDistribution.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CHEMISTRY/AASequence.h>

using namespace std;

namespace OpenMS
{

  TheoreticalSpectrumGenerator::TheoreticalSpectrumGenerator() :
    DefaultParamHandler("TheoreticalSpectrumGenerator")
  {
    defaults_.setValue("add_isotopes", "false", "If set to 1 isotope peaks of the product ion peaks are added");
    defaults_.setValidStrings("add_isotopes", ListUtils::create<String>("true,false"));

    defaults_.setValue("max_isotope", 2, "Defines the maximal isotopic peak which is added, add_isotopes must be set to 1");

    defaults_.setValue("add_metainfo", "false", "Adds the type of peaks as metainfo to the peaks, like y8+, [M-H2O+2H]++");
    defaults_.setValidStrings("add_metainfo", ListUtils::create<String>("true,false"));

    defaults_.setValue("add_losses", "false", "Adds common losses to those ion expect to have them, only water and ammonia loss is considered");
    defaults_.setValidStrings("add_losses", ListUtils::create<String>("true,false"));

    defaults_.setValue("add_precursor_peaks", "false", "Adds peaks of the precursor to the spectrum, which happen to occur sometimes");
    defaults_.setValidStrings("add_precursor_peaks", ListUtils::create<String>("true,false"));

    defaults_.setValue("add_abundant_immonium_ions", "false", "Add most abundant immonium ions");
    defaults_.setValidStrings("add_abundant_immonium_ions", ListUtils::create<String>("true,false"));

    defaults_.setValue("add_first_prefix_ion", "false", "If set to true e.g. b1 ions are added");
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

  TheoreticalSpectrumGenerator::TheoreticalSpectrumGenerator(const TheoreticalSpectrumGenerator & rhs) :
    DefaultParamHandler(rhs)
  {
  }

  TheoreticalSpectrumGenerator & TheoreticalSpectrumGenerator::operator=(const TheoreticalSpectrumGenerator & rhs)
  {
    if (this != &rhs)
    {
      DefaultParamHandler::operator=(rhs);
    }
    return *this;
  }

  TheoreticalSpectrumGenerator::~TheoreticalSpectrumGenerator()
  {
  }

  void TheoreticalSpectrumGenerator::getSpectrum(RichPeakSpectrum & spec, const AASequence & peptide, Int charge) const
  {
    for (Int z = 1; z <= charge; ++z)
    {
      if (add_b_ions_)
        addPeaks(spec, peptide, Residue::BIon, z);
      if (add_y_ions_)
        addPeaks(spec, peptide, Residue::YIon, z);
      if (add_a_ions_)
        addPeaks(spec, peptide, Residue::AIon, z);
      if (add_c_ions_)
        addPeaks(spec, peptide, Residue::CIon, z);
      if (add_x_ions_)
        addPeaks(spec, peptide, Residue::XIon, z);
      if (add_z_ions_)
        addPeaks(spec, peptide, Residue::ZIon, z);
    }

    if (add_precursor_peaks)
    {
      addPrecursorPeaks(spec, peptide, charge);
    }

    if (add_abundant_immonium_ions)
    {
      addAbundantImmoniumIons(spec);
    }

    return;
  }

  void TheoreticalSpectrumGenerator::addAbundantImmoniumIons(RichPeakSpectrum & spec) const
  {
    RichPeak1D p;

    // just in case someone wants the ion names;
    p.metaRegistry().registerName("IonName", "Name of the ion");

    // Histidin immonium ion
    p.setMZ(110.0718);
    p.setIntensity(1.0);
    if (add_metainfo_)
    {
      String name("iH");
      p.setMetaValue("IonName", name);
    }
    spec.push_back(p);

    // Phenylalanin immonium ion
    p.setMZ(120.0813);
    p.setIntensity(1.0);
    if (add_metainfo_)
    {
      String name("iF");
      p.setMetaValue("IonName", name);
    }
    spec.push_back(p);

    // Tyrosine immonium ion
    p.setMZ(136.0762);
    p.setIntensity(1.0);
    if (add_metainfo_)
    {
      String name("iY");
      p.setMetaValue("IonName", name);
    }
    spec.push_back(p);

    // Iso/Leucin immonium ion (same mass for immonium ion)
    p.setMZ(86.09698);
    p.setIntensity(1.0);
    if (add_metainfo_)
    {
      String name("iL/I");
      p.setMetaValue("IonName", name);
    }
    spec.push_back(p);

    // Tryptophan immonium ion
    p.setMZ(159.0922);
    p.setIntensity(1.0);
    if (add_metainfo_)
    {
      String name("iW");
      p.setMetaValue("IonName", name);
    }
    spec.push_back(p);

    spec.sortByPosition();
  }

  char TheoreticalSpectrumGenerator::residueTypeToIonLetter_(Residue::ResidueType res_type) const
  {
    switch (res_type)
    {
      case Residue::AIon: return 'a'; break;
      case Residue::BIon: return 'b'; break;
      case Residue::CIon: return 'c'; break;
      case Residue::XIon: return 'x'; break;
      case Residue::YIon: return 'y'; break;
      case Residue::ZIon: return 'z'; break;
      default:
       cerr << "Unknown residue type encountered. Can't map to ion letter." << endl;
    }
    return ' ';
  }

  void TheoreticalSpectrumGenerator::addIsotopeCluster_(RichPeakSpectrum & spectrum, const AASequence & ion, Residue::ResidueType res_type, Int charge, double intensity) const
  {
    double pos = ion.getMonoWeight(res_type, charge) / (double)charge;
    RichPeak1D p;
    IsotopeDistribution dist = ion.getFormula(res_type, charge).getIsotopeDistribution(max_isotope_);
    UInt j(0);
    for (IsotopeDistribution::ConstIterator it = dist.begin(); it != dist.end(); ++it, ++j)
    {
      // TODO: this is usually dominated by 13C-12C mass shift which deviates a bit from neutron mass
      p.setMZ((double)(pos + (double)j * Constants::NEUTRON_MASS_U) / (double)charge); 
      p.setIntensity(intensity * it->second);
      if (add_metainfo_ && j == 0)
      {
        String ion_name = String(residueTypeToIonLetter_(res_type)) + String(ion.size()) + String(charge, '+');
        p.setMetaValue("IonName", ion_name);
      }
      spectrum.push_back(p);
    }
  }

  void TheoreticalSpectrumGenerator::addPeak_(RichPeakSpectrum & spectrum, double pos, double intensity, Residue::ResidueType res_type, Size ion_index, int charge) const 
  {
    RichPeak1D p;
    p.setMZ(pos);
    p.setIntensity(intensity);
    if (add_metainfo_)
    {
      String ion_name = String(residueTypeToIonLetter_(res_type)) + String(ion_index) + String(charge, '+');
      p.setMetaValue("IonName", ion_name);
    }
    spectrum.push_back(p);    
  }

  void TheoreticalSpectrumGenerator::addLosses_(RichPeakSpectrum & spectrum, const AASequence & ion, double intensity, Residue::ResidueType res_type, int charge) const 
  {
    RichPeak1D p;

    set<String> losses;
    for (AASequence::ConstIterator it = ion.begin(); it != ion.end(); ++it)
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
      const String& loss_name = *it;

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
            // note: important to construct a string from char. If omitted it will perform pointer arithmetics on the "-" string literal
            String ion_name = String(residueTypeToIonLetter_(res_type)) + String(ion.size()) + "-" + loss_name + String(charge, '+');
            p.setMetaValue("IonName", ion_name);
          }
          spectrum.push_back(p);
        }
      }
      else
      {
        p.setMZ(loss_pos);
        if (add_metainfo_)
        {
          // note: important to construct a string from char. If omitted it will perform pointer arithmetics on the "-" string literal
          String ion_name = String(residueTypeToIonLetter_(res_type)) + String(ion.size()) + "-" + loss_name + String(charge, '+');
          p.setMetaValue("IonName", ion_name);
        }
        spectrum.push_back(p);
      }
    }

  }

  void TheoreticalSpectrumGenerator::addPeaks(RichPeakSpectrum & spectrum, const AASequence & peptide, Residue::ResidueType res_type, Int charge) const
  {
    if (peptide.empty())
    {
      return;
    }

    spectrum.reserve(peptide.size());

    // Generate the ion peaks:
    // Does not generate peaks of full peptide (therefore "<").
    // They are added via precursor mass (and neutral losses).
    // Could be changed in the future.
    switch (res_type)
    {
    case Residue::AIon:
    {
      Size i = 1;
      if (!add_first_prefix_ion_)
      {
        i = 2;
      }
      for (; i < peptide.size(); ++i)
      {
        AASequence ion = peptide.getPrefix(i);

        if (add_isotopes_) // add isotope cluster
        {
          addIsotopeCluster_(spectrum, ion, res_type, charge, a_intensity_);
        }
        else  // add single peak
        {
          double pos = ion.getMonoWeight(Residue::AIon, charge) / static_cast<double>(charge);
          addPeak_(spectrum, pos, a_intensity_, res_type, i, charge);
        }

        if (add_losses_)
        {
          addLosses_(spectrum, ion, a_intensity_, res_type, charge);
        }
      }
      break;
    }

    case Residue::BIon:
    {
      Size i = 1;
      if (!add_first_prefix_ion_)
      {
        i = 2;
      }
      for (; i < peptide.size(); ++i)
      {
        AASequence ion = peptide.getPrefix(i);
        if (add_isotopes_)
        {
          addIsotopeCluster_(spectrum, ion, res_type, charge, b_intensity_);
        }
        else
        {
          double pos = ion.getMonoWeight(Residue::BIon, charge) / static_cast<double>(charge);
          addPeak_(spectrum, pos, b_intensity_, res_type, i, charge);
        }

        if (add_losses_)
        {
          addLosses_(spectrum, ion, b_intensity_, res_type, charge);
        }
      }
      break;
    }

    case Residue::CIon:
    {
      Size i = 1;
      if (!add_first_prefix_ion_)
      {
        i = 2;
      }
      if (peptide.size() < 2)
      {
        //"Cannot create c ions of a monomer."
        throw Exception::InvalidSize(__FILE__, __LINE__, __PRETTY_FUNCTION__, 1);
      }
      for (; i < peptide.size(); ++i)
      {
        AASequence ion = peptide.getPrefix(i);
        if (add_isotopes_)
        {
          addIsotopeCluster_(spectrum, ion, res_type, charge, c_intensity_);
        }
        else
        {
          double pos = ion.getMonoWeight(Residue::CIon, charge) / static_cast<double>(charge);
          addPeak_(spectrum, pos, c_intensity_, res_type, i, charge);
        }

        if (add_losses_)
        {
          addLosses_(spectrum, ion, c_intensity_, res_type, charge);
        }
      }
      break;
    }

    case Residue::XIon:
    {
      Size i = 1;
      if (!add_first_prefix_ion_)
      {
        i = 2;
      }
      if (peptide.size() < 2)
      {
        // "Cannot create c ions of a monomer."
        throw Exception::InvalidSize(__FILE__, __LINE__, __PRETTY_FUNCTION__, 1);
      }
      for (; i < peptide.size(); ++i)
      {
        AASequence ion = peptide.getSuffix(i);
        if (add_isotopes_)
        {
          addIsotopeCluster_(spectrum, ion, res_type, charge, x_intensity_);
        }
        else
        {
          double pos = ion.getMonoWeight(Residue::XIon, charge) / static_cast<double>(charge);
          addPeak_(spectrum, pos, x_intensity_, res_type, i, charge);
        }

        if (add_losses_)
        {
          addLosses_(spectrum, ion, x_intensity_, res_type, charge);
        }
      }

      break;
    }

    case Residue::YIon:
    {
      for (Size i = 1; i < peptide.size(); ++i)
      {
        AASequence ion = peptide.getSuffix(i);
        if (add_isotopes_)
        {
          addIsotopeCluster_(spectrum, ion, res_type, charge, y_intensity_);
        }
        else
        {
          double pos = ion.getMonoWeight(Residue::YIon, charge) / static_cast<double>(charge);
          addPeak_(spectrum, pos, y_intensity_, res_type, i, charge);
        }

        if (add_losses_)
        {
          addLosses_(spectrum, ion, y_intensity_, res_type, charge);
        }
      }

      break;
    }

    case Residue::ZIon:
    {
      for (Size i = 1; i < peptide.size(); ++i)
      {
        AASequence ion = peptide.getSuffix(i);
        if (add_isotopes_)
        {
          addIsotopeCluster_(spectrum, ion, res_type, charge, z_intensity_);
        }
        else
        {
          double pos = ion.getMonoWeight(Residue::ZIon, charge) / static_cast<double>(charge);
          addPeak_(spectrum, pos, z_intensity_, res_type, i, charge);
        }

        if (add_losses_)
        {
          addLosses_(spectrum, ion, z_intensity_, res_type, charge);
        }
      }

      break;
    }

    default:
      cerr << "Cannot create peaks of that ion type" << endl;
    }

    spectrum.sortByPosition();

    return;
  }

  void TheoreticalSpectrumGenerator::addPrecursorPeaks(RichPeakSpectrum & spec, const AASequence & peptide, Int charge) const
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

  void TheoreticalSpectrumGenerator::updateMembers_()
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
    add_precursor_peaks = param_.getValue("add_precursor_peaks").toBool();
    add_abundant_immonium_ions = param_.getValue("add_abundant_immonium_ions").toBool();
    a_intensity_ = (double)param_.getValue("a_intensity");
    b_intensity_ = (double)param_.getValue("b_intensity");
    c_intensity_ = (double)param_.getValue("c_intensity");
    x_intensity_ = (double)param_.getValue("x_intensity");
    y_intensity_ = (double)param_.getValue("y_intensity");
    z_intensity_ = (double)param_.getValue("z_intensity");
    max_isotope_ = (Int)param_.getValue("max_isotope");
    rel_loss_intensity_ = (double)param_.getValue("relative_loss_intensity");
    pre_int_ = (double)param_.getValue("precursor_intensity");
    pre_int_H2O_ = (double)param_.getValue("precursor_H2O_intensity");
    pre_int_NH3_ = (double)param_.getValue("precursor_NH3_intensity");
  }

}
