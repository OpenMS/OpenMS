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
// $Maintainer: Eugen Netz $
// $Authors: Eugen Netz $
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGeneratorXLMS.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/KERNEL/MSSpectrum.h>


using namespace std;

namespace OpenMS
{

  TheoreticalSpectrumGeneratorXLMS::TheoreticalSpectrumGeneratorXLMS() :
    DefaultParamHandler("TheoreticalSpectrumGeneratorXLMS")
  {
    // TODO only partly functional (second isotopic peak if max_isotope = 2)
    defaults_.setValue("add_isotopes", "false", "If set to 1 isotope peaks of the product ion peaks are added");
    defaults_.setValidStrings("add_isotopes", ListUtils::create<String>("true,false"));

    defaults_.setValue("max_isotope", 2, "Defines the maximal isotopic peak which is added, add_isotopes must be set to 1");

    defaults_.setValue("add_metainfo", "true", "Adds the type of peaks as metainfo to the peaks, like y8+, [M-H2O+2H]++");
    defaults_.setValidStrings("add_metainfo", ListUtils::create<String>("true,false"));

    // TODO not functional yet
    defaults_.setValue("add_losses", "false", "Adds common losses to those ion expect to have them, only water and ammonia loss is considered");
    defaults_.setValidStrings("add_losses", ListUtils::create<String>("true,false"));

    // TODO not functional yet
    defaults_.setValue("add_precursor_peaks", "false", "Adds peaks of the precursor to the spectrum, which happen to occur sometimes");
    defaults_.setValidStrings("add_precursor_peaks", ListUtils::create<String>("true,false"));

    defaults_.setValue("add_abundant_immonium_ions", "false", "Add most abundant immonium ions");
    defaults_.setValidStrings("add_abundant_immonium_ions", ListUtils::create<String>("true,false"));

    // TODO not functional yet
    defaults_.setValue("add_first_prefix_ion", "true", "If set to true e.g. b1 ions are added");
    defaults_.setValidStrings("add_first_prefix_ion", ListUtils::create<String>("true,false"));

    // TODO not functional yet
    defaults_.setValue("multiple_fragmentation_mode" , "false", "If set to true, multiple fragmentation events on the same cross-linked peptide pair are considered (HCD fragmentation)");
    defaults_.setValidStrings("multiple_fragmentation_mode", ListUtils::create<String>("true,false"));

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

  TheoreticalSpectrumGeneratorXLMS::TheoreticalSpectrumGeneratorXLMS(const TheoreticalSpectrumGeneratorXLMS & rhs) :
    DefaultParamHandler(rhs)
  {
  }

  TheoreticalSpectrumGeneratorXLMS & TheoreticalSpectrumGeneratorXLMS::operator=(const TheoreticalSpectrumGeneratorXLMS & rhs)
  {
    if (this != &rhs)
    {
      DefaultParamHandler::operator=(rhs);
    }
    return *this;
  }

  TheoreticalSpectrumGeneratorXLMS::~TheoreticalSpectrumGeneratorXLMS()
  {
  }

  void TheoreticalSpectrumGeneratorXLMS::getCommonIonSpectrum(PeakSpectrum & spectrum, AASequence peptide, Size link_pos, bool frag_alpha, int charge, Size link_pos_2) const
  {
    PeakSpectrum::IntegerDataArray charges;
    PeakSpectrum::StringDataArray ion_names;


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


    for (Int z = 1; z <= charge; ++z)
    {
      if (add_b_ions_)
      {
        addCommonPeaks_(spectrum, charges, ion_names, peptide, link_pos, frag_alpha, Residue::BIon, z, link_pos_2);
      }
      if (add_y_ions_)
      {
        addCommonPeaks_(spectrum, charges, ion_names, peptide, link_pos, frag_alpha, Residue::YIon, z, link_pos_2);
      }
      if (add_a_ions_)
      {
        addCommonPeaks_(spectrum, charges, ion_names, peptide, link_pos, frag_alpha, Residue::AIon, z, link_pos_2);
      }
      if (add_x_ions_)
      {
        addCommonPeaks_(spectrum, charges, ion_names, peptide, link_pos, frag_alpha, Residue::XIon, z, link_pos_2);
      }
      if (add_c_ions_)
      {
        addCommonPeaks_(spectrum, charges, ion_names, peptide, link_pos, frag_alpha, Residue::CIon, z, link_pos_2);
      }
      if (add_z_ions_)
      {
        addCommonPeaks_(spectrum, charges, ion_names, peptide, link_pos, frag_alpha, Residue::ZIon, z, link_pos_2);
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

void TheoreticalSpectrumGeneratorXLMS::addCommonPeaks_(PeakSpectrum & spectrum, DataArrays::IntegerDataArray & charges, DataArrays::StringDataArray & ion_names, AASequence peptide, Size link_pos, bool frag_alpha, Residue::ResidueType res_type, int charge, Size link_pos_2) const
  {
    if (peptide.empty())
    {
      cout << "Warning: Attempt at creating XLink Ions Spectrum from empty string!" << endl;
      return;
    }

    String ion_type;
    if (frag_alpha)
    {
      ion_type = "alpha|ci";
    }
    else
    {
      ion_type = "beta|ci";
    }

    // second link position, in case of a loop-link
    Size link_pos_B = link_pos_2;
    if (link_pos_2 == 0)
    {
      link_pos_B = link_pos;
    }
//    cout << "Link_pos: " << link_pos << " | Link_pos_B: " << link_pos_B << " | charge: " << static_cast<double>(charge) << endl;
//    cout << "Fragmented Peptide: " << peptide.toUnmodifiedString() << " | Size: " << peptide.size() << endl;

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

    if (res_type == Residue::AIon || res_type == Residue::BIon || res_type == Residue::CIon)
    {
//      if ((!add_isotopes_) || max_isotope_ <= 2) // add single peaks (and maybe a second isotopic peak)
//      {
        double mono_weight(Constants::PROTON_MASS_U * static_cast<double>(charge));
        if (peptide.hasNTerminalModification())
        {
          mono_weight += peptide.getNTerminalModification()->getDiffMonoMass();
        }

        switch (res_type)
        {
          case Residue::AIon: mono_weight += Residue::getInternalToAIon().getMonoWeight(); break;
          case Residue::BIon: mono_weight += Residue::getInternalToBIon().getMonoWeight(); break;
          case Residue::CIon: mono_weight += Residue::getInternalToCIon().getMonoWeight(); break;
          default: break;
        }

        Size i = 0;
        for (; i < link_pos; ++i)
        {
          mono_weight += peptide[i].getMonoWeight(Residue::Internal);
          double pos(mono_weight / static_cast<double>(charge));
          int frag_index = i+1;

          // fragment testing
//          double recalc_pos = (peptide.getPrefix(i+1).getMonoWeight(Residue::BIon) + static_cast<double>(charge) ) / static_cast<double>(charge);
//          cout << "Current residue: " << i << " = " << peptide.toUnmodifiedString()[i] << " | ion_type: " << ion_type << "$b" << frag_index << " | pos: " << pos << " | recalc_pos: " << recalc_pos << endl;

          addPeak_(spectrum, charges, ion_names, pos, intensity, res_type, frag_index, charge, ion_type);
          if (add_isotopes_ && max_isotope_ >= 2) // add second isotopic peak with fast method, if two or more peaks are asked for
          {
            pos += Constants::C13C12_MASSDIFF_U / static_cast<double>(charge);
            addPeak_(spectrum, charges, ion_names, pos, intensity, res_type, frag_index, charge, ion_type);
          }
        }
//      }
    }
    else // if (res_type == Residue::XIon || res_type == Residue::YIon || res_type == Residue::ZIon)
    {
//      if ((!add_isotopes_) || max_isotope_ <= 2) // add single peaks (and maybe a second isotopic peak)
//      {
        double mono_weight(Constants::PROTON_MASS_U * static_cast<double>(charge));
        if (peptide.hasCTerminalModification())
        {
          mono_weight += peptide.getCTerminalModification()->getDiffMonoMass();
        }

        switch (res_type)
        {
          case Residue::XIon: mono_weight += Residue::getInternalToXIon().getMonoWeight(); break;
          case Residue::YIon: mono_weight += Residue::getInternalToYIon().getMonoWeight(); break;
          case Residue::ZIon: mono_weight += Residue::getInternalToZIon().getMonoWeight(); break;
          default: break;
        }

        for (Size i = peptide.size()-1; i > link_pos_B; --i)
        {
          mono_weight += peptide[i].getMonoWeight(Residue::Internal);
          double pos(mono_weight / static_cast<double>(charge));
          int frag_index = peptide.size() - i;

//          double recalc_pos = (peptide.getSuffix(frag_index).getMonoWeight(Residue::YIon) + static_cast<double>(charge) ) / static_cast<double>(charge);
//          cout << "Current residue: " << i << " = " << peptide.toUnmodifiedString()[i] << " | ion_type: " << ion_type << "$y" << frag_index << " | pos: " << pos << " | recalc_pos: " << recalc_pos << endl;

          addPeak_(spectrum, charges, ion_names, pos, intensity, res_type, frag_index, charge, ion_type);
          if (add_isotopes_ && max_isotope_ >= 2) // add second isotopic peak with fast method, if two or more peaks are asked for
          {
            pos += Constants::C13C12_MASSDIFF_U / static_cast<double>(charge);
            addPeak_(spectrum, charges, ion_names, pos, intensity, res_type, frag_index, charge, ion_type);
          }
        }
//      }
    }

    return;
  }

  void TheoreticalSpectrumGeneratorXLMS::getXLinkIonSpectrum(PeakSpectrum & spectrum, AASequence peptide, Size link_pos, double precursor_mass, bool frag_alpha, int mincharge, int maxcharge, Size link_pos_2) const
  {
    PeakSpectrum::IntegerDataArray charges;
    PeakSpectrum::StringDataArray ion_names;

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

    for (Int z = mincharge; z <= maxcharge; ++z)
    {
      if (add_b_ions_)
      {
        addXLinkIonPeaks_(spectrum, charges, ion_names, peptide, link_pos, precursor_mass, frag_alpha, Residue::BIon, z, link_pos_2);
      }
      if (add_y_ions_)
      {
        addXLinkIonPeaks_(spectrum, charges, ion_names, peptide, link_pos, precursor_mass, frag_alpha, Residue::YIon, z, link_pos_2);
      }
      if (add_a_ions_)
      {
        addXLinkIonPeaks_(spectrum, charges, ion_names, peptide, link_pos, precursor_mass, frag_alpha, Residue::AIon, z, link_pos_2);
      }
      if (add_x_ions_)
      {
        addXLinkIonPeaks_(spectrum, charges, ion_names, peptide, link_pos, precursor_mass, frag_alpha, Residue::XIon, z, link_pos_2);
      }
      if (add_c_ions_)
      {
        addXLinkIonPeaks_(spectrum, charges, ion_names, peptide, link_pos, precursor_mass, frag_alpha, Residue::CIon, z, link_pos_2);
      }
      if (add_z_ions_)
      {
        addXLinkIonPeaks_(spectrum, charges, ion_names, peptide, link_pos, precursor_mass, frag_alpha, Residue::ZIon, z, link_pos_2);
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

  void TheoreticalSpectrumGeneratorXLMS::addXLinkIonPeaks_(PeakSpectrum & spectrum, DataArrays::IntegerDataArray & charges, DataArrays::StringDataArray & ion_names, AASequence peptide, Size link_pos, double precursor_mass, bool frag_alpha, Residue::ResidueType res_type, int charge, Size link_pos_2) const
  {
    if (peptide.empty())
    {
      cout << "Warning: Attempt at creating XLink Ions Spectrum from empty string!" << endl;
      return;
    }

    String ion_type;
    if (frag_alpha)
    {
      ion_type = "alpha|xi";
    } else
    {
      ion_type = "beta|xi";
    }

    // second link position, in case of a loop-link
    Size link_pos_B = link_pos_2;
    if (link_pos_2 == 0)
    {
      link_pos_B = link_pos;
    }
//    cout << "Link_pos: " << link_pos << " | Link_pos_B: " << link_pos_B << " | charge: " << static_cast<double>(charge) << endl;
//    cout << "Fragmented Peptide: " << peptide.toUnmodifiedString() << " | Size: " << peptide.size() << " | precursor_mass: " << precursor_mass << endl;

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

    if (res_type == Residue::AIon || res_type == Residue::BIon || res_type == Residue::CIon)
    {
//      if ((!add_isotopes_) || max_isotope_ <= 2) // add single peaks (and maybe a second isotopic peak)
//      {
        // whole mass of both peptides + cross-link (or peptide + mono-link), converted to an internal ion
        double mono_weight((Constants::PROTON_MASS_U * static_cast<double>(charge)) + precursor_mass - Residue::getInternalToFull().getMonoWeight());


        if (peptide.hasCTerminalModification())
        {
          mono_weight -= peptide.getCTerminalModification()->getDiffMonoMass();
        }

        // adjust mass to given residue type
        switch (res_type)
        {
          case Residue::AIon: mono_weight += Residue::getInternalToAIon().getMonoWeight(); break;
          case Residue::BIon: mono_weight += Residue::getInternalToBIon().getMonoWeight(); break;
          case Residue::CIon: mono_weight += Residue::getInternalToCIon().getMonoWeight(); break;
          default: break;
        }

        // subtract one residue at a time
        for (Size i = peptide.size()-1; i > link_pos_B; --i)
        {
          mono_weight -= peptide[i].getMonoWeight(Residue::Internal);
          double pos(mono_weight / static_cast<double>(charge));
          int frag_index = i;

//          double recalc_pos = (peptide.getPrefix(frag_index+1).getMonoWeight() + static_cast<double>(charge) ) / static_cast<double>(charge);
//          cout << "Current residue: " << i << " = " << peptide.toUnmodifiedString()[i] << " | ion_type: " << ion_type << "$b" << frag_index << " | pos: " << pos << "| current_residue_mass: " << peptide[i].getMonoWeight(Residue::Internal) << " | recalc_pos: " << recalc_pos << endl;

          addPeak_(spectrum, charges, ion_names, pos, intensity, res_type, frag_index, charge, ion_type);
          if (add_isotopes_ && max_isotope_ >= 2) // add second isotopic peak with fast method, if two or more peaks are asked for
          {
            pos += Constants::C13C12_MASSDIFF_U / static_cast<double>(charge);
            addPeak_(spectrum, charges, ion_names, pos, intensity, res_type, frag_index, charge, ion_type);
          }
        }
//      }
    }
    else // if (res_type == Residue::XIon || res_type == Residue::YIon || res_type == Residue::ZIon)
    {
        // not needed yet, since there is no alternative yet
//      if ((!add_isotopes_) || max_isotope_ <= 2) // add single peaks (and maybe a second isotopic peak)
//      {
        // whole mass of both peptides + cross-link (or peptide + mono-link), converted to an internal ion
        double mono_weight((Constants::PROTON_MASS_U * static_cast<double>(charge)) + precursor_mass - Residue::getInternalToFull().getMonoWeight()); // whole mass

        if (peptide.hasNTerminalModification())
        {
          mono_weight -= peptide.getNTerminalModification()->getDiffMonoMass();
        }

        // adjust mass to given residue type
        switch (res_type)
        {
          case Residue::XIon: mono_weight += Residue::getInternalToXIon().getMonoWeight(); break;
          case Residue::YIon: mono_weight += Residue::getInternalToYIon().getMonoWeight(); break;
          case Residue::ZIon: mono_weight += Residue::getInternalToZIon().getMonoWeight(); break;
          default: break;
        }

        // subtract one residue at a time
        for (Size i = 0; i < link_pos; ++i)
        {
          mono_weight -= peptide[i].getMonoWeight(Residue::Internal);
          double pos(mono_weight / static_cast<double>(charge));
          int frag_index = peptide.size() - 1 - i;

//          double recalc_pos = (peptide.getSuffix(frag_index+1).getMonoWeight() + static_cast<double>(charge) ) / static_cast<double>(charge);
//          cout << "Current residue: " << i << " = " << peptide.toUnmodifiedString()[i] << " | ion_type: " << ion_type << "$y" << frag_index << " | pos: " << pos << "| current_residue_mass: " << peptide[i].getMonoWeight(Residue::Internal) << " | recalc_pos: " << recalc_pos << endl;

          addPeak_(spectrum, charges, ion_names, pos, intensity, res_type, frag_index, charge, ion_type);
          if (add_isotopes_ && max_isotope_ >= 2) // add second isotopic peak with fast method, if two or more peaks are asked for
          {
            pos += Constants::C13C12_MASSDIFF_U / static_cast<double>(charge);
            addPeak_(spectrum, charges, ion_names, pos, intensity, res_type, frag_index, charge, ion_type);
          }
        }
//      }
    }
    return;
  }


  // helper to add a single peak to a spectrum (simple fragmentation)
  void TheoreticalSpectrumGeneratorXLMS::addPeak_(PeakSpectrum & spectrum, DataArrays::IntegerDataArray & charges, DataArrays::StringDataArray & ion_names, double pos, double intensity, Residue::ResidueType res_type, Size frag_index, int charge, String ion_type) const
  {
    Peak1D p;
    p.setMZ(pos);
    p.setIntensity(intensity);
    spectrum.push_back(p);
    if (add_metainfo_)
    {
      String ion_name = "[" + ion_type + "$" + String(residueTypeToIonLetter_(res_type)) + String(frag_index) + "]"; //+ String(charge, '+');
      ion_names.push_back(ion_name);
      charges.push_back(charge);
    }
  }

  // helper for mapping residue type to letter
  char TheoreticalSpectrumGeneratorXLMS::residueTypeToIonLetter_(Residue::ResidueType res_type) const
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

  void TheoreticalSpectrumGeneratorXLMS::updateMembers_()
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
    multiple_fragmentation_mode_ = param_.getValue("multiple_fragmentation_mode").toBool();
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
