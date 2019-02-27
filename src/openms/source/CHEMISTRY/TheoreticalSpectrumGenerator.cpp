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
// $Maintainer: Timo Sachsenberg, Eugen Netz $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>

#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/NASequence.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>

#include <OpenMS/KERNEL/MSSpectrum.h>

using namespace std;

namespace OpenMS
{

  TheoreticalSpectrumGenerator::TheoreticalSpectrumGenerator() :
    DefaultParamHandler("TheoreticalSpectrumGenerator")
  {
    defaults_.setValue("add_isotopes", "false", "If set to 'true' isotope peaks of the product ion peaks are added");
    defaults_.setValidStrings("add_isotopes", ListUtils::create<String>("true,false"));

    defaults_.setValue("max_isotope", 2, "Defines the maximal isotopic peak which is added if 'add_isotopes' is 'true'");

    defaults_.setValue("add_metainfo", "false", "Adds the type of peaks as metainfo to the peaks, like y8+, [M-H2O+2H]++");
    defaults_.setValidStrings("add_metainfo", ListUtils::create<String>("true,false"));

    defaults_.setValue("add_losses", "false", "Adds common losses to those ion expect to have them, only water and ammonia loss is considered");
    defaults_.setValidStrings("add_losses", ListUtils::create<String>("true,false"));

    defaults_.setValue("add_precursor_peaks", "false", "Adds peaks of the precursor to the spectrum, which happen to occur sometimes");
    defaults_.setValidStrings("add_precursor_peaks", ListUtils::create<String>("true,false"));

    defaults_.setValue("add_all_precursor_charges", "false", "Adds precursor peaks with all charges in the given range");
    defaults_.setValidStrings("add_all_precursor_charges", ListUtils::create<String>("true,false"));

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

    defaults_.setValue("add_d_ions", "false", "Add peaks of d-ions to the spectrum");
    defaults_.setValidStrings("add_d_ions", ListUtils::create<String>("true,false"));

    defaults_.setValue("add_w_ions", "false", "Add peaks of w-ions to the spectrum");
    defaults_.setValidStrings("add_w_ions", ListUtils::create<String>("true,false"));

    defaults_.setValue("add_a-B_ions", "false", "Add peaks of a-B-ions to the spectrum");
    defaults_.setValidStrings("add_a-B_ions", ListUtils::create<String>("true,false"));

    // intensity options of the ions
    defaults_.setValue("y_intensity", 1.0, "Intensity of the y-ions");
    defaults_.setValue("b_intensity", 1.0, "Intensity of the b-ions");
    defaults_.setValue("a_intensity", 1.0, "Intensity of the a-ions");
    defaults_.setValue("c_intensity", 1.0, "Intensity of the c-ions");
    defaults_.setValue("x_intensity", 1.0, "Intensity of the x-ions");
    defaults_.setValue("z_intensity", 1.0, "Intensity of the z-ions");

    defaults_.setValue("d_intensity", 1.0, "Intensity of the d-ions");
    defaults_.setValue("w_intensity", 1.0, "Intensity of the w-ions");
    defaults_.setValue("a-B_intensity", 1.0, "Intensity of the a-B-ions");

    defaults_.setValue("relative_loss_intensity", 0.1, "Intensity of loss ions, in relation to the intact ion intensity");

    // precursor intensity
    defaults_.setValue("precursor_intensity", 1.0, "Intensity of the precursor peak");
    defaults_.setValue("precursor_H2O_intensity", 1.0, "Intensity of the H2O loss peak of the precursor");
    defaults_.setValue("precursor_NH3_intensity", 1.0, "Intensity of the NH3 loss peak of the precursor");

    defaultsToParam_();
  }


  TheoreticalSpectrumGenerator::TheoreticalSpectrumGenerator(const TheoreticalSpectrumGenerator& rhs) :
    DefaultParamHandler(rhs)
  {
  }


  TheoreticalSpectrumGenerator& TheoreticalSpectrumGenerator::operator=(const TheoreticalSpectrumGenerator& rhs)
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


  void TheoreticalSpectrumGenerator::getSpectrum(PeakSpectrum& spectrum, const NASequence& nucleotide, Int min_charge, Int max_charge) const
  {
    Int sign = 1;
    if (max_charge < 0 && min_charge < 0)
    {
      sign = -1;
    }
    else if ((max_charge > 0 && min_charge < 0) ||
             (max_charge < 0 && min_charge > 0))
    {
      //Signs don't match we need to quit and thow error here to avoid messing up for loops below
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "min. and max. charge must both be either positive or negative");
    }
    else if (max_charge < min_charge)
    {
      swap(max_charge, min_charge);
    }
    PeakSpectrum::StringDataArray ion_names;
    PeakSpectrum::IntegerDataArray charges;

    if (add_metainfo_)
    {
      if (spectrum.getIntegerDataArrays().size() > 0)
      {
        charges = spectrum.getIntegerDataArrays()[0];
      }
      if (spectrum.getStringDataArrays().size() > 0)
      {
        ion_names = spectrum.getStringDataArrays()[0];
      }
      ion_names.setName("IonNames");
      charges.setName("Charges");
    }

    for (uint z = (uint)abs(min_charge); z <= (uint)abs(max_charge) && z < (uint)nucleotide.size(); ++z)
    {
      if (add_b_ions_) addPeaks_(spectrum, nucleotide, ion_names, charges, NASequence::BIon, z * sign);
      if (add_y_ions_) addPeaks_(spectrum, nucleotide, ion_names, charges, NASequence::YIon, z * sign);
      if (add_a_ions_) addPeaks_(spectrum, nucleotide, ion_names, charges, NASequence::AIon, z * sign);
      if (add_c_ions_) addPeaks_(spectrum, nucleotide, ion_names, charges, NASequence::CIon, z * sign);
      if (add_x_ions_) addPeaks_(spectrum, nucleotide, ion_names, charges, NASequence::XIon, z * sign);
      if (add_z_ions_) addPeaks_(spectrum, nucleotide, ion_names, charges, NASequence::ZIon, z * sign);
      if (add_d_ions_) addPeaks_(spectrum, nucleotide, ion_names, charges, NASequence::DIon, z * sign);
      if (add_w_ions_) addPeaks_(spectrum, nucleotide, ion_names, charges, NASequence::WIon, z * sign);
      if (add_aB_ions_) addPeaks_(spectrum, nucleotide, ion_names, charges, NASequence::AminusB, z * sign);
    }

    if (add_precursor_peaks_)
    {
      if (add_all_precursor_charges_)
      {
        for (uint z = (uint)abs(min_charge); z <= (uint)abs(max_charge); ++z)
        {
          addPrecursorPeaks_(spectrum, nucleotide, ion_names, charges, z * sign);
        }
      }
      else // add_all_precursor_charges_ = false, only add precursor with highest charge
      {
        addPrecursorPeaks_(spectrum, nucleotide, ion_names, charges, max_charge);
      }
    }

    if (add_metainfo_)
    {
      if (spectrum.getIntegerDataArrays().size() > 0)
      {
        spectrum.getIntegerDataArrays()[0] = charges;
      }
      else
      {
        spectrum.getIntegerDataArrays().push_back(charges);
      }
      if (spectrum.getStringDataArrays().size() > 0)
      {
        spectrum.getStringDataArrays()[0] = ion_names;
      }
      else
      {
        spectrum.getStringDataArrays().push_back(ion_names);
      }
    }

    spectrum.sortByPosition();
    return;
  }


  void TheoreticalSpectrumGenerator::getSpectrum(PeakSpectrum& spectrum, const AASequence& peptide, Int min_charge, Int max_charge) const
  {
    if (peptide.empty())
    {
      return;
    }

    PeakSpectrum::StringDataArray ion_names;
    PeakSpectrum::IntegerDataArray charges;

    if (add_metainfo_)
    {
      if (spectrum.getIntegerDataArrays().size() > 0)
      {
        charges = spectrum.getIntegerDataArrays()[0];
      }
      if (spectrum.getStringDataArrays().size() > 0)
      {
        ion_names = spectrum.getStringDataArrays()[0];
      }
      ion_names.setName("IonNames");
      charges.setName("Charges");
    }

    for (Int z = min_charge; z <= max_charge; ++z)
    {
      if (add_b_ions_) addPeaks_(spectrum, peptide, ion_names, charges, Residue::BIon, z);
      if (add_y_ions_) addPeaks_(spectrum, peptide, ion_names, charges, Residue::YIon, z);
      if (add_a_ions_) addPeaks_(spectrum, peptide, ion_names, charges, Residue::AIon, z);
      if (add_c_ions_) addPeaks_(spectrum, peptide, ion_names, charges, Residue::CIon, z);
      if (add_x_ions_) addPeaks_(spectrum, peptide, ion_names, charges, Residue::XIon, z);
      if (add_z_ions_) addPeaks_(spectrum, peptide, ion_names, charges, Residue::ZIon, z);
    }

    if (add_precursor_peaks_)
    {
      if (add_all_precursor_charges_)
      {
        for (Int z = min_charge; z <= max_charge; ++z)
        {
          addPrecursorPeaks_(spectrum, peptide, ion_names, charges, z);
        }
      }
      else // add_all_precursor_charges_ = false, only add precursor with highest charge
      {
        addPrecursorPeaks_(spectrum, peptide, ion_names, charges, max_charge);
      }
    }

    if (add_abundant_immonium_ions_)
    {
      addAbundantImmoniumIons_(spectrum, peptide, ion_names, charges);
    }

    if (add_metainfo_)
    {
      if (spectrum.getIntegerDataArrays().size() > 0)
      {
        spectrum.getIntegerDataArrays()[0] = charges;
      }
      else
      {
        spectrum.getIntegerDataArrays().push_back(charges);
      }
      if (spectrum.getStringDataArrays().size() > 0)
      {
        spectrum.getStringDataArrays()[0] = ion_names;
      }
      else
      {
        spectrum.getStringDataArrays().push_back(ion_names);
      }
    }

    spectrum.sortByPosition();
    return;
  }


  void TheoreticalSpectrumGenerator::addAbundantImmoniumIons_(PeakSpectrum& spectrum, const AASequence& peptide, DataArrays::StringDataArray& ion_names, DataArrays::IntegerDataArray& charges) const
  {
    Peak1D p;

    // Histidin immonium ion (C5H8N3)
    if (peptide.has(*ResidueDB::getInstance()->getResidue('H')))
    {
      p.setMZ(110.0718);
      p.setIntensity(1.0);
      if (add_metainfo_)
      {
        String ion_name("iH");
        ion_names.push_back(ion_name);
        charges.push_back(1);
      }
      spectrum.push_back(p);
    }

    // Phenylalanin immonium ion (C8H10N)
    if (peptide.has(*ResidueDB::getInstance()->getResidue('F')))
    {
      p.setMZ(120.0813);
      p.setIntensity(1.0);
      if (add_metainfo_)
      {
        String ion_name("iF");
        ion_names.push_back(ion_name);
        charges.push_back(1);
      }
      spectrum.push_back(p);
    }

    // Tyrosine immonium ion (C8H10NO)
    if (peptide.has(*ResidueDB::getInstance()->getResidue('Y')))
    {
      p.setMZ(136.0762);
      p.setIntensity(1.0);
      if (add_metainfo_)
      {
        String ion_name("iY");
        ion_names.push_back(ion_name);
        charges.push_back(1);
      }
      spectrum.push_back(p);
    }

    // Iso/Leucin immonium ion (same mass for immonium ion)
    if (peptide.has(*ResidueDB::getInstance()->getResidue('L')))
    {
      p.setMZ(86.09698);
      p.setIntensity(1.0);
      if (add_metainfo_)
      {
        String ion_name("iL/I");
        ion_names.push_back(ion_name);
        charges.push_back(1);
      }
      spectrum.push_back(p);
    }

    // Tryptophan immonium ion
    if (peptide.has(*ResidueDB::getInstance()->getResidue('W')))
    {
      p.setMZ(159.0922);
      p.setIntensity(1.0);
      if (add_metainfo_)
      {
        String ion_name("iW");
        ion_names.push_back(ion_name);
        charges.push_back(1);
      }
      spectrum.push_back(p);
    }

    // Cysteine (C2H6NS)
    if (peptide.has(*ResidueDB::getInstance()->getResidue('C')))
    {
      p.setMZ(76.0221);
      p.setIntensity(1.0);
      if (add_metainfo_)
      {
        String ion_name("iC");
        ion_names.push_back(ion_name);
        charges.push_back(1);
      }
      spectrum.push_back(p);
    }

    // Proline immonium ion (C4H8N)
    if (peptide.has(*ResidueDB::getInstance()->getResidue('P')))
    {
      p.setMZ(70.0656);
      p.setIntensity(1.0);
      if (add_metainfo_)
      {
        String ion_name("iP");
        ion_names.push_back(ion_name);
        charges.push_back(1);
      }
      spectrum.push_back(p);
    }
  }


  char TheoreticalSpectrumGenerator::residueTypeToIonLetter_(Residue::ResidueType res_type)
  {
    switch (res_type)
    {
      case Residue::AIon: return 'a';
      case Residue::BIon: return 'b';
      case Residue::CIon: return 'c';
      case Residue::XIon: return 'x';
      case Residue::YIon: return 'y';
      case Residue::ZIon: return 'z';
      default:
       cerr << "Unknown residue type encountered. Can't map to ion letter." << endl;
    }
    return ' ';
  }


  String TheoreticalSpectrumGenerator::ribonucleotideTypeToIonCode_(NASequence::NASFragmentType type, Size num)
  {
    String result;
    switch (type)
    {
    case NASequence::AIon: result = "a"; break;
    case NASequence::BIon: result = "b"; break;
    case NASequence::CIon: result = "c"; break;
    case NASequence::XIon: result = "x"; break;
    case NASequence::YIon: result = "y"; break;
    case NASequence::ZIon: result = "z"; break;
    case NASequence::DIon: result = "d"; break;
    case NASequence::WIon: result = "w"; break;
    case NASequence::AminusB: return (num > 0) ? "a" + String(num) + "-B" : "a-B";
    default:
      cerr << "Unknown ribonucleotide type encountered. Can't map to ion code." << endl;
      result = "?";
    }
    if (num > 0) result += String(num);
    return result;
  }


  void TheoreticalSpectrumGenerator::addIsotopeCluster_(PeakSpectrum& spectrum, const AASequence& ion, DataArrays::StringDataArray& ion_names, DataArrays::IntegerDataArray& charges, Residue::ResidueType res_type, Int charge, double intensity) const
  {
    double pos = ion.getMonoWeight(res_type, charge);
    Peak1D p;
    IsotopeDistribution dist = ion.getFormula(res_type, charge).getIsotopeDistribution(CoarseIsotopePatternGenerator(max_isotope_));

    String ion_name = String(Residue::residueTypeToIonLetter(res_type)) + String(ion.size()) + String((Size)abs(charge), '+');

    double j(0.0);
    for (IsotopeDistribution::ConstIterator it = dist.begin(); it != dist.end(); ++it, ++j)
    {
      // TODO: this is usually dominated by 13C-12C mass shift which deviates a bit from neutron mass
      p.setMZ((double)(pos + j * Constants::C13C12_MASSDIFF_U) / (double)charge);
      p.setIntensity(intensity * it->getIntensity());
      if (add_metainfo_) // one entry per peak
      {
        ion_names.push_back(ion_name);
        charges.push_back(charge);
      }
      spectrum.push_back(p);
    }
  }


  void TheoreticalSpectrumGenerator::addIsotopeCluster_(PeakSpectrum& spectrum, const NASequence& ion, DataArrays::StringDataArray& ion_names, DataArrays::IntegerDataArray& charges, NASequence::NASFragmentType res_type, Int charge, double intensity) const
  {
    double charge_unit = (charge < 0) ? -1.0 : 1.0; // pos. or neg. charge?
    double pos = ion.getMonoWeight(res_type, charge);
    Peak1D p;
    IsotopeDistribution dist = ion.getFormula(res_type, charge).getIsotopeDistribution(CoarseIsotopePatternGenerator(max_isotope_));

    String ion_name = ribonucleotideTypeToIonCode_(res_type, ion.size()); // + String((Size)abs(charge), charge_sign);

    double j(0.0);
    for (IsotopeDistribution::ConstIterator it = dist.begin(); it != dist.end(); ++it, ++j)
    {
      // TODO: this is usually dominated by 13C-12C mass shift which deviates a bit from neutron mass
      p.setMZ((double)(pos + j * Constants::C13C12_MASSDIFF_U * charge_unit) / (double)abs(charge));
      p.setIntensity(intensity * it->getIntensity());
      if (add_metainfo_) // one entry per peak
      {
        ion_names.push_back(ion_name);
        charges.push_back(charge);
      }
      spectrum.push_back(p);
    }
  }


  void TheoreticalSpectrumGenerator::addLosses_(PeakSpectrum& spectrum, const AASequence& ion, DataArrays::StringDataArray& ion_names, DataArrays::IntegerDataArray& charges, double intensity, Residue::ResidueType res_type, int charge) const
  {
    Peak1D p;

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
      double loss_pos = loss_ion.getMonoWeight();
      const String& loss_name = *it;

      if (add_isotopes_)
      {
        IsotopeDistribution dist = loss_ion.getIsotopeDistribution(CoarseIsotopePatternGenerator(max_isotope_));

        // note: important to construct a string from char. If omitted it will perform pointer arithmetics on the "-" string literal
        String ion_name = String(Residue::residueTypeToIonLetter(res_type)) + String(ion.size()) + "-" + loss_name + String((Size)abs(charge), '+');

        double j(0.0);
        for (IsotopeDistribution::ConstIterator iso = dist.begin(); iso != dist.end(); ++iso, ++j)
        {
          p.setMZ((double)(loss_pos + j * Constants::C13C12_MASSDIFF_U) / (double)charge);
          p.setIntensity(intensity * rel_loss_intensity_ * iso->getIntensity());
          if (add_metainfo_)
          {
            ion_names.push_back(ion_name);
            charges.push_back(charge);
          }
          spectrum.push_back(p);
        }
      }
      else
      {
        p.setMZ(loss_pos / (double)charge);
        if (add_metainfo_)
        {
          // note: important to construct a string from char. If omitted it will perform pointer arithmetics on the "-" string literal
          String ion_name = String(Residue::residueTypeToIonLetter(res_type)) + String(ion.size()) + "-" + loss_name + String((Size)abs(charge), '+');
          ion_names.push_back(ion_name);
          charges.push_back(charge);
        }
        spectrum.push_back(p);
      }
    }
  }


  void TheoreticalSpectrumGenerator::addPeaks_(PeakSpectrum& spectrum, const NASequence& nucleotide, DataArrays::StringDataArray& ion_names, DataArrays::IntegerDataArray& charges, NASequence::NASFragmentType res_type, Int charge) const
  {
    spectrum.reserve(nucleotide.size());

    // Generate the ion peaks:
    // Does not generate peaks of full sequence (therefore "<").
    // They are added via precursor mass (and neutral losses).
    // Could be changed in the future.

    double intensity(1);
    // char charge_sign = (charge > 0) ? '+' : '-';

    switch (res_type)
    {
    case NASequence::AIon:
      intensity = a_intensity_;
      break;
    case NASequence::BIon:
      intensity = b_intensity_;
      break;
    case NASequence::CIon:
      intensity = c_intensity_;
      break;
    case NASequence::XIon:
      intensity = x_intensity_;
      break;
    case NASequence::YIon:
      intensity = y_intensity_;
      break;
    case NASequence::ZIon:
      intensity = z_intensity_;
      break;
    case NASequence::DIon:
      intensity = d_intensity_;
      break;
    case NASequence::WIon:
      intensity = w_intensity_;
      break;
    case NASequence::AminusB:
      intensity = aB_intensity_;
      break;
    default:
      break;
    }

    if (res_type == NASequence::AIon || res_type == NASequence::BIon || res_type == NASequence::CIon || res_type == NASequence::AminusB || res_type == NASequence::DIon)
    {
      // @TODO: special cases for a-B ions ("NASequence::AminusB")
      // - they may not be relevant for fragments of length 1 (unless modified?)
      // - mods on the last base (that gets lost) may be completely or partially
      // retained/lost, depending on the mod (whether it's on the base or on the
      // backbone or both - we don't have that information at the moment)

      if (!add_isotopes_) // add single peak
      {
        Size length = add_first_prefix_ion_ ? 1 : 2;
        for (; length < nucleotide.size(); ++length)
        {
          NASequence ion = nucleotide.getPrefix(length);
          double mass = ion.getMonoWeight(res_type, charge);
          Peak1D p;
          p.setMZ(mass / abs(charge));
          p.setIntensity(intensity);
          spectrum.push_back(p);
          if (add_metainfo_)
          {
            String ion_name = ribonucleotideTypeToIonCode_(res_type, length); // + String((Size)abs(charge), charge_sign);
            ion_names.push_back(ion_name);
            charges.push_back(charge);
          }
        }
      }
      else // add isotope clusters (slow)
      {
        Size length = add_first_prefix_ion_ ? 1 : 2;
        for (; length < nucleotide.size(); ++length)
        {
          const NASequence ion = nucleotide.getPrefix(length);
          addIsotopeCluster_(spectrum, ion, ion_names, charges, res_type, charge, intensity); //TODO IMPLEMENT
        }
      }

//      if (add_losses_) // add loss peaks (slow)
//      {
//        Size i = add_first_prefix_ion_ ? 1 : 2;
//        for (; i < nucleotide.size(); ++i)
//        {
//          const NASequence ion = nucleotide.getPrefix(i);
//          addLosses_(spectrum, ion, ion_names, charges, intensity, res_type, charge);
//        }
//      }
    }
    else // WIon, XIon, YIon, ZIon
    {
      if (!add_isotopes_) // add single peak
      {
        for (Size length = 1; length < nucleotide.size(); ++length)
        {
          NASequence ion = nucleotide.getSuffix(length);
          double mass = ion.getMonoWeight(res_type, charge);
          Peak1D p;
          p.setMZ(mass / abs(charge));
          p.setIntensity(intensity);
          spectrum.push_back(p);
          if (add_metainfo_)
          {
            String ion_name = ribonucleotideTypeToIonCode_(res_type, length); // + String((Size)abs(charge), charge_sign);
            ion_names.push_back(ion_name);
            charges.push_back(charge);
          }
        }
      }
      else // add isotope clusters
      {
        for (Size length = 1; length < nucleotide.size(); ++length)
        {
          const NASequence ion = nucleotide.getSuffix(length);
          addIsotopeCluster_(spectrum, ion, ion_names, charges, res_type, charge, intensity); //TODO IMPLEMENT
        }
      }

      //if (add_losses_) // add loss peaks (slow)
      //{
      //  for (Size i = 1; i < nucleotide.size(); ++i)
       // {
        //  const NASequence ion = nucleotide.getSuffix(i);
        //  addLosses_(spectrum, ion, ion_names, charges, intensity, res_type, charge);
       // }
     // }
    }

    return;
  }


  void TheoreticalSpectrumGenerator::addPeaks_(PeakSpectrum& spectrum, const AASequence& peptide, DataArrays::StringDataArray& ion_names, DataArrays::IntegerDataArray& charges, Residue::ResidueType res_type, Int charge) const
  {
    int f = 1 + int(add_isotopes_) + int(add_losses_);
    spectrum.reserve(spectrum.size() + f * peptide.size());

    // Generate the ion peaks:
    // Does not generate peaks of full peptide (therefore "<").
    // They are added via precursor mass (and neutral losses).
    // Could be changed in the future.

    double intensity(1);

    switch (res_type)
    {
      case Residue::AIon: intensity = a_intensity_; break;
      case Residue::BIon: intensity = b_intensity_; break;
      case Residue::CIon: if (peptide.size() < 2) throw Exception::InvalidSize(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 1); intensity = c_intensity_; break;
      case Residue::XIon: if (peptide.size() < 2) throw Exception::InvalidSize(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, 1); intensity = x_intensity_; break;
      case Residue::YIon: intensity = y_intensity_; break;
      case Residue::ZIon: intensity = z_intensity_; break;
      default: break;
    }

    double mono_weight(Constants::PROTON_MASS_U * charge);

    if (res_type == Residue::AIon || res_type == Residue::BIon || res_type == Residue::CIon)
    {
      if (peptide.hasNTerminalModification())
      {
        mono_weight += peptide.getNTerminalModification()->getDiffMonoMass();
      }

      if (!add_isotopes_) // add single peak
      {
        Size i = add_first_prefix_ion_ ? 0 : 1;
        if (i == 1) mono_weight += peptide[0].getMonoWeight(Residue::Internal);
        for (; i < peptide.size() - 1; ++i)
        {
          mono_weight += peptide[i].getMonoWeight(Residue::Internal); // standard internal residue including named modifications
          double pos(mono_weight);
          switch (res_type)
          {
          case Residue::AIon: pos = (pos + Residue::getInternalToAIon().getMonoWeight()) / charge; break;
          case Residue::BIon: pos = (pos + Residue::getInternalToBIon().getMonoWeight()) / charge; break;
          case Residue::CIon: pos = (pos + Residue::getInternalToCIon().getMonoWeight()) / charge; break;
          default: break;
          }
          Peak1D p;
          p.setMZ(pos);
          p.setIntensity(intensity);
          spectrum.push_back(p);
          if (add_metainfo_)
          {
            String ion_name = String(Residue::residueTypeToIonLetter(res_type)) + String(i + 1) + String((Size)abs(charge), '+');
            ion_names.push_back(ion_name);
            charges.push_back(charge);
          }
        }
      }
      else // add isotope clusters (slow)
      {
        Size i = add_first_prefix_ion_ ? 1 : 2;
        for (; i < peptide.size(); ++i)
        {
          const AASequence ion = peptide.getPrefix(i);
          addIsotopeCluster_(spectrum, ion, ion_names, charges, res_type, charge, intensity);
        }
      }

      if (add_losses_) // add loss peaks (slow)
      {
        Size i = add_first_prefix_ion_ ? 1 : 2;
        for (; i < peptide.size(); ++i)
        {
          const AASequence ion = peptide.getPrefix(i);
          addLosses_(spectrum, ion, ion_names, charges, intensity, res_type, charge);
        }
      }
    }
    else // if (res_type == Residue::XIon || res_type == Residue::YIon || res_type == Residue::ZIon)
    {
      if (peptide.hasCTerminalModification())
      {
        mono_weight += peptide.getCTerminalModification()->getDiffMonoMass();
      }

      if (!add_isotopes_) // add single peak
      {
        Size i = peptide.size() - 1;

        for (; i > 0; --i)
        {
          mono_weight += peptide[i].getMonoWeight(Residue::Internal); // standard internal residue including named modifications
          double pos(mono_weight);
          switch (res_type)
          {
          case Residue::XIon: pos = (pos + Residue::getInternalToXIon().getMonoWeight()) / charge; break;
          case Residue::YIon: pos = (pos + Residue::getInternalToYIon().getMonoWeight()) / charge; break;
          case Residue::ZIon: pos = (pos + Residue::getInternalToZIon().getMonoWeight()) / charge; break;
          default: break;
          }
          Peak1D p;
          p.setMZ(pos);
          p.setIntensity(intensity);
          spectrum.push_back(p);
          if (add_metainfo_)
          {
            String ion_name = String(Residue::residueTypeToIonLetter(res_type)) + String(peptide.size() - i) + String((Size)abs(charge), '+');

            ion_names.push_back(ion_name);
            charges.push_back(charge);
          }
        }
      }
      else // add isotope clusters
      {
        for (Size i = 1; i < peptide.size(); ++i)
        {
          const AASequence ion = peptide.getSuffix(i);
          addIsotopeCluster_(spectrum, ion, ion_names, charges, res_type, charge, intensity);
        }
      }

      if (add_losses_) // add loss peaks (slow)
      {
        for (Size i = 1; i < peptide.size(); ++i)
        {
          const AASequence ion = peptide.getSuffix(i);
          addLosses_(spectrum, ion, ion_names, charges, intensity, res_type, charge);
        }
      }
    }

    return;
  }


  void TheoreticalSpectrumGenerator::addPrecursorPeaks_(PeakSpectrum& spectrum, const AASequence& peptide, DataArrays::StringDataArray& ion_names, DataArrays::IntegerDataArray& charges, Int charge) const
  {
    Peak1D p;

    String ion_name("[M+H]" + String((Size)abs(charge), '+'));

    // precursor peak
    double mono_pos = peptide.getMonoWeight(Residue::Full, charge);

    if (add_isotopes_)
    {
      IsotopeDistribution dist = peptide.getFormula(Residue::Full, charge).getIsotopeDistribution(CoarseIsotopePatternGenerator(max_isotope_));
      double j(0.0);
      for (IsotopeDistribution::ConstIterator it = dist.begin(); it != dist.end(); ++it, ++j)
      {
        p.setMZ((double)(mono_pos + j * Constants::C13C12_MASSDIFF_U) / (double)charge);
        p.setIntensity(pre_int_ * it->getIntensity());
        if (add_metainfo_)
        {
          ion_names.push_back(ion_name);
          charges.push_back(charge);
        }
        spectrum.push_back(p);
      }
    }
    else
    {
      p.setMZ(mono_pos / (double)charge);
      p.setIntensity(pre_int_);
      if (add_metainfo_)
      {
        ion_names.push_back(ion_name);
        charges.push_back(charge);
      }
      spectrum.push_back(p);
    }
    // loss peaks of the precursor

    //loss of water
    EmpiricalFormula ion = peptide.getFormula(Residue::Full, charge) - EmpiricalFormula("H2O");
    mono_pos = ion.getMonoWeight();
    if (add_isotopes_)
    {
      IsotopeDistribution dist = ion.getIsotopeDistribution(CoarseIsotopePatternGenerator(max_isotope_));
      UInt j(0);
      for (IsotopeDistribution::ConstIterator it = dist.begin(); it != dist.end(); ++it, ++j)
      {
        p.setMZ((double)(mono_pos + j * Constants::C13C12_MASSDIFF_U) / (double)charge);
        p.setIntensity(pre_int_H2O_ *  it->getIntensity());
        if (add_metainfo_)
        {
          String ion_name("[M+H]-H2O" + String((Size)abs(charge), '+'));
          ion_names.push_back(ion_name);
          charges.push_back(charge);
        }
        spectrum.push_back(p);
      }
    }
    else
    {
      p.setMZ(mono_pos / (double)charge);
      p.setIntensity(pre_int_H2O_);
      if (add_metainfo_)
      {
        String ion_name("[M+H]-H2O" + String((Size)abs(charge), '+'));
        ion_names.push_back(ion_name);
        charges.push_back(charge);
      }
      spectrum.push_back(p);
    }

    //loss of ammonia
    ion = peptide.getFormula(Residue::Full, charge) - EmpiricalFormula("NH3");
    mono_pos = ion.getMonoWeight();
    if (add_isotopes_)
    {
      IsotopeDistribution dist = ion.getIsotopeDistribution(CoarseIsotopePatternGenerator(max_isotope_));
      UInt j(0);
      for (IsotopeDistribution::ConstIterator it = dist.begin(); it != dist.end(); ++it, ++j)
      {
        p.setMZ((double)(mono_pos + j * Constants::C13C12_MASSDIFF_U) / (double)charge);
        p.setIntensity(pre_int_NH3_ *  it->getIntensity());
        if (add_metainfo_)
        {
          String ion_name("[M+H]-NH3" + String((Size)abs(charge), '+'));
          ion_names.push_back(ion_name);
          charges.push_back(charge);
        }
        spectrum.push_back(p);
      }
    }
    else
    {
      p.setMZ(mono_pos / (double)charge);
      p.setIntensity(pre_int_NH3_);
      if (add_metainfo_)
      {
        String ion_name("[M+H]-NH3" + String((Size)abs(charge), '+'));
        ion_names.push_back(ion_name);
        charges.push_back(charge);
      }
      spectrum.push_back(p);
    }
  }


  void TheoreticalSpectrumGenerator::addPrecursorPeaks_(PeakSpectrum& spectrum, const NASequence& nucleotide, DataArrays::StringDataArray& ion_names, DataArrays::IntegerDataArray& charges, Int charge) const
  {
    Peak1D p;

    char charge_sign = '+';
    double charge_unit = 1.0; //double equivalent of sign
    if (charge < 0)
    {
      charge_sign = '-';
      charge_unit = -1.0;
    }

    String ion_name("[M"+ String(charge_sign)+"H]" + String((Size)abs(charge), charge_sign));

    // precursor peak
    double mono_pos = nucleotide.getMonoWeight(NASequence::Full, charge);

    if (add_isotopes_)
    {
      IsotopeDistribution dist = nucleotide.getFormula(NASequence::Full, charge).getIsotopeDistribution(CoarseIsotopePatternGenerator(max_isotope_));
      double j(0.0);
      for (IsotopeDistribution::ConstIterator it = dist.begin(); it != dist.end(); ++it, ++j)
      {
        p.setMZ((double)(mono_pos + j * Constants::C13C12_MASSDIFF_U * charge_unit) / (double)(abs(charge)));
        p.setIntensity(pre_int_ * it->getIntensity());
        if (add_metainfo_)
        {
          ion_names.push_back(ion_name);
          charges.push_back(charge);
        }
        spectrum.push_back(p);
      }
    }
    else
    {
      p.setMZ(mono_pos / (double)(abs(charge)));
      p.setIntensity(pre_int_);
      if (add_metainfo_)
      {
        ion_names.push_back(ion_name);
        charges.push_back(charge);
      }
      spectrum.push_back(p);
    }
    // loss peaks of the precursor

    //loss of water
    EmpiricalFormula ion = nucleotide.getFormula(NASequence::Full, charge) - EmpiricalFormula("H2O");
    mono_pos = ion.getMonoWeight();
    if (add_isotopes_)
    {
      IsotopeDistribution dist = ion.getIsotopeDistribution(CoarseIsotopePatternGenerator(max_isotope_));
      UInt j(0);
      for (IsotopeDistribution::ConstIterator it = dist.begin(); it != dist.end(); ++it, ++j)
      {
        p.setMZ((double)(mono_pos + j * Constants::C13C12_MASSDIFF_U * charge_unit ) / (double)(abs(charge)));
        p.setIntensity(pre_int_H2O_ *  it->getIntensity());
        if (add_metainfo_)
        {
          String ion_name("[M"+ String(charge_sign)+"H]-H2O" + String((Size)abs(charge), charge_sign));
          ion_names.push_back(ion_name);
          charges.push_back(charge);
        }
        spectrum.push_back(p);
      }
    }
    else
    {
      p.setMZ(mono_pos / (double)(abs(charge)));
      p.setIntensity(pre_int_H2O_);
      if (add_metainfo_)
      {
        String ion_name("[M"+String(charge_sign)+"H]-H2O" + String((Size)abs(charge), charge_sign));
        ion_names.push_back(ion_name);
        charges.push_back(charge);
      }
      spectrum.push_back(p);
    }

    //loss of ammonia
    ion = nucleotide.getFormula(NASequence::Full, charge) - EmpiricalFormula("NH3");
    mono_pos = ion.getMonoWeight();
    if (add_isotopes_)
    {
      IsotopeDistribution dist = ion.getIsotopeDistribution(CoarseIsotopePatternGenerator(max_isotope_));
      UInt j(0);
      for (IsotopeDistribution::ConstIterator it = dist.begin(); it != dist.end(); ++it, ++j)
      {
        p.setMZ((double)(mono_pos + j * Constants::C13C12_MASSDIFF_U * charge_unit ) / (double)(abs(charge)));
        p.setIntensity(pre_int_NH3_ *  it->getIntensity());
        if (add_metainfo_)
        {
          String ion_name("[M"+String(charge_sign)+"H]-NH3" + String((Size)abs(charge), charge_sign));
          ion_names.push_back(ion_name);
          charges.push_back(charge);
        }
        spectrum.push_back(p);
      }
    }
    else
    {
      p.setMZ(mono_pos / (double)(abs(charge)));
      p.setIntensity(pre_int_NH3_);
      if (add_metainfo_)
      {
        String ion_name("[M"+String(charge_sign)+"H]-NH3" + String((Size)abs(charge), charge_sign));
        ion_names.push_back(ion_name);
        charges.push_back(charge);
      }
      spectrum.push_back(p);
    }
  }


  void TheoreticalSpectrumGenerator::updateMembers_()
  {
    add_b_ions_ = param_.getValue("add_b_ions").toBool();
    add_y_ions_ = param_.getValue("add_y_ions").toBool();
    add_a_ions_ = param_.getValue("add_a_ions").toBool();
    add_c_ions_ = param_.getValue("add_c_ions").toBool();
    add_x_ions_ = param_.getValue("add_x_ions").toBool();
    add_z_ions_ = param_.getValue("add_z_ions").toBool();
    add_d_ions_ = param_.getValue("add_d_ions").toBool();
    add_w_ions_ = param_.getValue("add_w_ions").toBool();
    add_aB_ions_ = param_.getValue("add_a-B_ions").toBool();
    add_first_prefix_ion_ = param_.getValue("add_first_prefix_ion").toBool();
    add_losses_ = param_.getValue("add_losses").toBool();
    add_metainfo_ = param_.getValue("add_metainfo").toBool();
    add_isotopes_ = param_.getValue("add_isotopes").toBool();
    add_precursor_peaks_ = param_.getValue("add_precursor_peaks").toBool();
    add_all_precursor_charges_ = param_.getValue("add_all_precursor_charges").toBool();
    add_abundant_immonium_ions_ = param_.getValue("add_abundant_immonium_ions").toBool();
    a_intensity_ = (double)param_.getValue("a_intensity");
    aB_intensity_ = (double)param_.getValue("a-B_intensity");
    b_intensity_ = (double)param_.getValue("b_intensity");
    c_intensity_ = (double)param_.getValue("c_intensity");
    d_intensity_ = (double)param_.getValue("d_intensity");
    w_intensity_ = (double)param_.getValue("w_intensity");
    x_intensity_ = (double)param_.getValue("x_intensity");
    y_intensity_ = (double)param_.getValue("y_intensity");
    z_intensity_ = (double)param_.getValue("z_intensity");
    max_isotope_ = (Int)param_.getValue("max_isotope");
    rel_loss_intensity_ = (double)param_.getValue("relative_loss_intensity");
    pre_int_ = (double)param_.getValue("precursor_intensity");
    pre_int_H2O_ = (double)param_.getValue("precursor_H2O_intensity");
    pre_int_NH3_ = (double)param_.getValue("precursor_NH3_intensity");
  }

} // end namespace OpenMS
