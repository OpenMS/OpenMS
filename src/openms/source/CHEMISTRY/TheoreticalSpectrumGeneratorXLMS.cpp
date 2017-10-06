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

    defaults_.setValue("add_losses", "false", "Adds common losses to those ion expect to have them, only water and ammonia loss is considered");
    defaults_.setValidStrings("add_losses", ListUtils::create<String>("true,false"));

    defaults_.setValue("add_precursor_peaks", "true", "Adds peaks of the precursor to the spectrum, which happen to occur sometimes");
    defaults_.setValidStrings("add_precursor_peaks", ListUtils::create<String>("true,false"));

    defaults_.setValue("add_abundant_immonium_ions", "false", "Add most abundant immonium ions");
    defaults_.setValidStrings("add_abundant_immonium_ions", ListUtils::create<String>("true,false"));

    defaults_.setValue("add_k_linked_ions", "true", "Add RES-Linked ions, which are specific to XLMS");
    defaults_.setValidStrings("add_k_linked_ions", ListUtils::create<String>("true,false"));

    // TODO not functional yet
    defaults_.setValue("add_first_prefix_ion", "true", "If set to true e.g. b1 ions are added");
    defaults_.setValidStrings("add_first_prefix_ion", ListUtils::create<String>("true,false"));

    defaults_.setValue("add_y_ions", "true", "Add peaks of y-ions to the spectrum");
    defaults_.setValidStrings("add_y_ions", ListUtils::create<String>("true,false"));

    defaults_.setValue("add_b_ions", "true", "Add peaks of b-ions to the spectrum");
    defaults_.setValidStrings("add_b_ions", ListUtils::create<String>("true,false"));

    defaults_.setValue("add_a_ions", "true", "Add peaks of a-ions to the spectrum");
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

  void TheoreticalSpectrumGeneratorXLMS::getLinearIonSpectrum(PeakSpectrum & spectrum, AASequence & peptide, Size link_pos, bool frag_alpha, int charge, Size link_pos_2) const
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

    std::vector< std::set< TheoreticalSpectrumGeneratorXLMS::LossMass > > forward_losses;
    std::vector< std::set< TheoreticalSpectrumGeneratorXLMS::LossMass > > backward_losses;

    if (add_losses_)
    {
      forward_losses = getForwardLossesForLinearIons_(peptide);
      backward_losses = getBackwardLossesForLinearIons_(peptide);
    }

    for (Int z = 1; z <= charge; ++z)
    {
      if (add_b_ions_)
      {
        addLinearPeaks_(spectrum, charges, ion_names, peptide, link_pos, frag_alpha, Residue::BIon, forward_losses, backward_losses, z, link_pos_2);
      }
      if (add_y_ions_)
      {
        addLinearPeaks_(spectrum, charges, ion_names, peptide, link_pos, frag_alpha, Residue::YIon, forward_losses, backward_losses, z, link_pos_2);
      }
      if (add_a_ions_)
      {
        addLinearPeaks_(spectrum, charges, ion_names, peptide, link_pos, frag_alpha, Residue::AIon, forward_losses, backward_losses, z, link_pos_2);
      }
      if (add_x_ions_)
      {
        addLinearPeaks_(spectrum, charges, ion_names, peptide, link_pos, frag_alpha, Residue::XIon, forward_losses, backward_losses, z, link_pos_2);
      }
      if (add_c_ions_)
      {
        addLinearPeaks_(spectrum, charges, ion_names, peptide, link_pos, frag_alpha, Residue::CIon, forward_losses, backward_losses, z, link_pos_2);
      }
      if (add_z_ions_)
      {
        addLinearPeaks_(spectrum, charges, ion_names, peptide, link_pos, frag_alpha, Residue::ZIon, forward_losses, backward_losses, z, link_pos_2);
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

  void TheoreticalSpectrumGeneratorXLMS::addLinearPeaks_(PeakSpectrum & spectrum, DataArrays::IntegerDataArray & charges, DataArrays::StringDataArray & ion_names, AASequence & peptide, Size link_pos, bool frag_alpha, Residue::ResidueType res_type, std::vector< std::set< TheoreticalSpectrumGeneratorXLMS::LossMass > > & forward_losses, std::vector< std::set< TheoreticalSpectrumGeneratorXLMS::LossMass > > & backward_losses, int charge, Size link_pos_2) const
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
          if (add_losses_)
          {
            // addLinearIonLosses_(spectrum, charges, ion_names, peptide.getPrefix(i+1), res_type, frag_index, intensity, charge, ion_type);
            addLinearIonLossesFAST_(spectrum, charges, ion_names, mono_weight, res_type, frag_index, intensity, charge, ion_type, forward_losses[i]);
          }
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
          if (add_losses_)
          {
            // addLinearIonLosses_(spectrum, charges, ion_names, peptide.getSuffix(peptide.size() - i), res_type, frag_index, intensity, charge, ion_type);
            addLinearIonLossesFAST_(spectrum, charges, ion_names, pos, res_type, frag_index, intensity, charge, ion_type, backward_losses[i]);
          }
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

  void TheoreticalSpectrumGeneratorXLMS::getXLinkIonSpectrum(PeakSpectrum & spectrum, AASequence & peptide, Size link_pos, double precursor_mass, bool frag_alpha, int mincharge, int maxcharge, Size link_pos_2) const
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
      if (add_k_linked_ions_)
      {
        addKLinkedIonPeaks_(spectrum, charges, ion_names, peptide, link_pos, precursor_mass, frag_alpha, z);
      }
    }

    if (add_precursor_peaks_)
    {
      addPrecursorPeaks_(spectrum, charges, ion_names, precursor_mass, maxcharge);
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

  void TheoreticalSpectrumGeneratorXLMS::addXLinkIonPeaks_(PeakSpectrum & spectrum, DataArrays::IntegerDataArray & charges, DataArrays::StringDataArray & ion_names, AASequence & peptide, Size link_pos, double precursor_mass, bool frag_alpha, Residue::ResidueType res_type, int charge, Size link_pos_2) const
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
    add_k_linked_ions_ = param_.getValue("add_k_linked_ions").toBool();
  }

  void TheoreticalSpectrumGeneratorXLMS::getComplexXLinkIonSpectrum(PeakSpectrum & spectrum, OPXLDataStructs::ProteinProteinCrossLink & crosslink, int mincharge, int maxcharge) const
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

    std::vector< std::set< TheoreticalSpectrumGeneratorXLMS::LossMass > > forward_losses_alpha;
    std::vector< std::set< TheoreticalSpectrumGeneratorXLMS::LossMass > > backward_losses_alpha;
    std::vector< std::set< TheoreticalSpectrumGeneratorXLMS::LossMass > > forward_losses_beta;
    std::vector< std::set< TheoreticalSpectrumGeneratorXLMS::LossMass > > backward_losses_beta;

    if (add_losses_)
    {
      forward_losses_alpha = getForwardLossesForLinearIons_(crosslink.alpha);
      backward_losses_alpha = getBackwardLossesForLinearIons_(crosslink.alpha);
      forward_losses_beta = getForwardLossesForLinearIons_(crosslink.beta);
      backward_losses_beta = getBackwardLossesForLinearIons_(crosslink.beta);
    }

    for (int z = mincharge; z <= maxcharge; ++z)
    {
      addComplexXLinkIonPeaks_(spectrum, charges, ion_names, crosslink, z, forward_losses_alpha, backward_losses_alpha, forward_losses_beta, backward_losses_beta);
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

  void TheoreticalSpectrumGeneratorXLMS::addLinearIonLosses_(PeakSpectrum & spectrum, DataArrays::IntegerDataArray& charges, DataArrays::StringDataArray& ion_names, const AASequence & ion, Residue::ResidueType res_type, Size frag_index, double intensity, int charge, String ion_type) const
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

    p.setIntensity(intensity * rel_loss_intensity_);

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
      String ion_name;

      p.setMZ(loss_pos / static_cast<double>(charge));
      if (add_metainfo_)
      {
        ion_name = "[" + ion_type + "$" + String(residueTypeToIonLetter_(res_type)) + String(frag_index) + "-" + loss_name + "]";
        ion_names.push_back(ion_name);
        charges.push_back(charge);
      }
      spectrum.push_back(p);

      if (add_isotopes_ && max_isotope_ >= 2) // add second isotopic peak with fast method, if two or more peaks are asked for
      {
        p.setMZ((loss_pos + Constants::C13C12_MASSDIFF_U) / static_cast<double>(charge));
        spectrum.push_back(p);
        if (add_metainfo_)
        {
          ion_names.push_back(ion_name);
          charges.push_back(charge);
        }
      }
    }
  }

  void TheoreticalSpectrumGeneratorXLMS::addLinearIonLossesFAST_(PeakSpectrum & spectrum, DataArrays::IntegerDataArray& charges, DataArrays::StringDataArray& ion_names, double mono_weight, Residue::ResidueType res_type, Size frag_index, double intensity, int charge, String ion_type, std::set< TheoreticalSpectrumGeneratorXLMS::LossMass > & losses) const
  {
    Peak1D p;
    p.setIntensity(intensity * rel_loss_intensity_);

    for (TheoreticalSpectrumGeneratorXLMS::LossMass loss : losses)
    {
      double mass_with_loss = mono_weight - loss.mass;
      String ion_name;

      if (mass_with_loss < 0.0) { continue; }

      p.setMZ(mass_with_loss / static_cast<double>(charge));
      if (add_metainfo_)
      {
        ion_name = "[" + ion_type + "$" + String(residueTypeToIonLetter_(res_type)) + String(frag_index) + "-" + loss.name + "]";
        ion_names.push_back(ion_name);
        charges.push_back(charge);
      }
      spectrum.push_back(p);

      // if (add_isotopes_ && max_isotope_ >= 2) // add second isotopic peak with fast method, if two or more peaks are asked for
      // {
      //   p.setMZ((mass_with_loss + Constants::C13C12_MASSDIFF_U) / static_cast<double>(charge));
      //   spectrum.push_back(p);
      //   if (add_metainfo_)
      //   {
      //     ion_names.push_back(ion_name);
      //     charges.push_back(charge);
      //   }
      // }
    }
  }

  void TheoreticalSpectrumGeneratorXLMS::addPrecursorPeaks_(PeakSpectrum & spectrum, DataArrays::IntegerDataArray & charges, DataArrays::StringDataArray & ion_names, double precursor_mass, int charge) const
  {
    // TODO precursors, mostly 4+ and 5+  (derive this and all other charges from precursor charges?)
    Peak1D p;
    String ion_name("[M+H]");

    // precursor peak
    double mono_pos = precursor_mass + (Constants::PROTON_MASS_U * static_cast<double>(charge));
    p.setMZ(mono_pos / static_cast<double>(charge));
    p.setIntensity(pre_int_);
    if (add_metainfo_)
    {
      ion_names.push_back(ion_name);
      charges.push_back(charge);
    }
    spectrum.push_back(p);
    if (add_isotopes_ && max_isotope_ >= 2) // add second isotopic peak with fast method, if two or more peaks are asked for
    {
      double pos = mono_pos + (Constants::C13C12_MASSDIFF_U / static_cast<double>(charge));
      p.setMZ(pos);
      p.setIntensity(pre_int_);
      if (add_metainfo_)
      {
        ion_names.push_back(ion_name);
        charges.push_back(charge);
      }
      spectrum.push_back(p);
    }

    // loss peaks of the precursor
    // loss of water
    mono_pos = precursor_mass + (Constants::PROTON_MASS_U * static_cast<double>(charge)) - EmpiricalFormula("H2O").getMonoWeight();
    p.setMZ(mono_pos / static_cast<double>(charge));
    p.setIntensity(pre_int_H2O_);
    if (add_metainfo_)
    {
      String ion_name("[M+H]-H2O");
      ion_names.push_back(ion_name);
      charges.push_back(charge);
    }
    spectrum.push_back(p);
    if (add_isotopes_ && max_isotope_ >= 2) // add second isotopic peak with fast method, if two or more peaks are asked for
    {
      double pos = mono_pos + (Constants::C13C12_MASSDIFF_U / static_cast<double>(charge));
      p.setMZ(pos);
      p.setIntensity(pre_int_H2O_);
      if (add_metainfo_)
      {
        ion_names.push_back(ion_name);
        charges.push_back(charge);
      }
      spectrum.push_back(p);
    }

    //loss of ammonia
    mono_pos = precursor_mass + (Constants::PROTON_MASS_U * static_cast<double>(charge)) - EmpiricalFormula("NH3").getMonoWeight();
    p.setMZ(mono_pos / static_cast<double>(charge));
    p.setIntensity(pre_int_NH3_);
    if (add_metainfo_)
    {
      String ion_name("[M+H]-NH3");
      ion_names.push_back(ion_name);
      charges.push_back(charge);
    }
    spectrum.push_back(p);
    if (add_isotopes_ && max_isotope_ >= 2) // add second isotopic peak with fast method, if two or more peaks are asked for
    {
      double pos = mono_pos + (Constants::C13C12_MASSDIFF_U / static_cast<double>(charge));
      p.setMZ(pos);
      p.setIntensity(pre_int_NH3_);
      if (add_metainfo_)
      {
        ion_names.push_back(ion_name);
        charges.push_back(charge);
      }
      spectrum.push_back(p);
    }
  }

  void TheoreticalSpectrumGeneratorXLMS::addKLinkedIonPeaks_(PeakSpectrum & spectrum, DataArrays::IntegerDataArray & charges, DataArrays::StringDataArray & ion_names, AASequence & peptide, Size link_pos, double precursor_mass, bool frag_alpha, int charge) const
  {
    double mono_weight = precursor_mass;
    // link_pos can be zero, if the cross-link is N-terminal
    if (link_pos > 1)
    {
      mono_weight -= peptide.getPrefix(link_pos-1).getMonoWeight(Residue::BIon);
    }
    // same here for C-terminal links
    if (link_pos < peptide.size()-1)
    {
      mono_weight -= peptide.getSuffix(peptide.size() - link_pos).getMonoWeight(Residue::XIon);
    }

    mono_weight += Constants::PROTON_MASS_U * static_cast<double>(charge);
    double pos(mono_weight / static_cast<double>(charge));

    Peak1D p;
    p.setMZ(pos);
    p.setIntensity(1.0);
    spectrum.push_back(p);

    // here the ion type is reversed compared to other peak types,
    // because for this special ion type, it would not make sense to call it alpha$y(n)-alpha$a(n)
    // Only one residue is left of the fragmented Peptide, so we call it a RES-linked beta
    String ion_type = "alpha";
    if (frag_alpha)
    {
      ion_type = "beta";
    }
    String ion_name;

    if (add_metainfo_)
    {

      int l_pos = link_pos;
      if (l_pos < 1)
      {
        l_pos = 0;
      }
      ion_name = "[" + peptide[l_pos].getOneLetterCode() + "-linked-" + ion_type + "]";
      ion_names.push_back(ion_name);
      charges.push_back(charge);
    }

    if (add_isotopes_ && max_isotope_ >= 2) // add second isotopic peak with fast method, if two or more peaks are asked for
    {
      pos += Constants::C13C12_MASSDIFF_U / static_cast<double>(charge);
      p.setMZ(pos);
      spectrum.push_back(p);
      if (add_metainfo_)
      {
        ion_names.push_back(ion_name);
        charges.push_back(charge);
      }
    }
  }

  void TheoreticalSpectrumGeneratorXLMS::addComplexXLinkIonPeaks_(PeakSpectrum & spectrum, DataArrays::IntegerDataArray & charges, DataArrays::StringDataArray & ion_names, OPXLDataStructs::ProteinProteinCrossLink & crosslink, int charge, std::vector< std::set< TheoreticalSpectrumGeneratorXLMS::LossMass > > & forward_losses_alpha, std::vector< std::set< TheoreticalSpectrumGeneratorXLMS::LossMass > > & backward_losses_alpha, std::vector< std::set< TheoreticalSpectrumGeneratorXLMS::LossMass > > & forward_losses_beta, std::vector< std::set< TheoreticalSpectrumGeneratorXLMS::LossMass > > & backward_losses_beta) const
  {
    // TODO fragment losses for cross-links
    // TODO yb, ya, all combos except alpha-beta aa, yy, bb,   only +1

    // Notes: b+1, y+1, a+1 and ya+1 : only linear
    //        b+2, b+3, a+2, y+3 : only xlinked
    //        y+2, yb+1 : both linear and xlinked

    AASequence & alpha = crosslink.alpha;
    AASequence & beta = crosslink.beta;

    if (alpha.empty() || beta.empty())
    {
      cout << "Warning: Attempt at creating XLink Ions Spectrum from empty string!" << endl;
      return;
    }

    String ion_name;

    Size link_pos_A = crosslink.cross_link_position.first;
    Size link_pos_B = crosslink.cross_link_position.second;

    // if its a terminal mod, the masses are the same, as when it is on the first or last residue
    // have to do this otherwise link_pos_A-1 = MAX_SIZE
    if (link_pos_A == 0) link_pos_A++;
    if (link_pos_B == 0) link_pos_B++;

    double intensity((b_intensity_ + y_intensity_) / 2.0);

    double precursor_mass = alpha.getMonoWeight() + beta.getMonoWeight() + crosslink.cross_linker_mass;

    // ALPHA-B BETA-Y
    // whole mass of both peptides + cross-link (or peptide + mono-link), converted to an internal ion
    double mono_weight((Constants::PROTON_MASS_U * static_cast<double>(charge)) + precursor_mass - (Residue::getInternalToFull().getMonoWeight()*2));

    if (alpha.hasCTerminalModification())
    {
      mono_weight -= alpha.getCTerminalModification()->getDiffMonoMass();
    }
    if (beta.hasNTerminalModification())
    {
      mono_weight -= beta.getNTerminalModification()->getDiffMonoMass();
    }

    mono_weight += Residue::getInternalToBIon().getMonoWeight(); // alpha is B ion
    mono_weight += Residue::getInternalToYIon().getMonoWeight(); // beta is Y ion

    double temp_full_beta_mass(mono_weight);

    // subtract one residue at a time from alpha, link_pos starts count at 1, index at 0
    for (Size i = alpha.size()-1; i > link_pos_A-1 && i > 0; --i)
    {
      // temp_full_beta_mass keeps track of the mass before residues of beta get subtracted
      mono_weight = temp_full_beta_mass;
      mono_weight -= alpha[i].getMonoWeight(Residue::Internal);
      temp_full_beta_mass = mono_weight;

      // subtract one residue at a time from beta
      for (Size j = 0; j < link_pos_B-1 && j < beta.size()-1; ++j)
      {
        mono_weight -= beta[j].getMonoWeight(Residue::Internal);
        double pos(mono_weight / static_cast<double>(charge));

        // if (pos < 0.0) { continue; }

        Peak1D p;
        p.setMZ(pos);
        p.setIntensity(intensity);
        spectrum.push_back(p);
        if (add_metainfo_)
        {
          ion_name = "[alpha|xi$b" + String(i) + "+beta|xi$y" + String(beta.size()-j-1) + "]";
          ion_names.push_back(ion_name);
          charges.push_back(charge);
        }
        // TODO losses removed for these complex peaks for now
        // if (add_losses_)
        // {
        //   std::set< TheoreticalSpectrumGeneratorXLMS::LossMass > losses = forward_losses_alpha[i-1];
        //   losses.insert(backward_losses_beta[j+1].begin(), backward_losses_beta[j+1].end());
        //   addXLinkIonLossesFAST_(spectrum, charges, ion_names, mono_weight, intensity, charge, ion_name, losses);
        //   // addXLinkIonLosses_(spectrum, charges, ion_names, alpha.getPrefix(i-1), beta.getSuffix(beta.size()-j-1), mono_weight, intensity, charge, ion_name);
        // }

        if (add_isotopes_ && max_isotope_ >= 2) // add second isotopic peak with fast method, if two or more peaks are asked for
        {
          pos += Constants::C13C12_MASSDIFF_U / static_cast<double>(charge);
          p.setMZ(pos);
          p.setIntensity(intensity);
          spectrum.push_back(p);
          if (add_metainfo_)
          {
            ion_names.push_back(ion_name);
            charges.push_back(charge);
          }
        }
      }
    }

    // ALPHA-Y BETA-B
    // whole mass of both peptides + cross-link (or peptide + mono-link), converted to an internal ion
    mono_weight = (Constants::PROTON_MASS_U * static_cast<double>(charge)) + precursor_mass - (Residue::getInternalToFull().getMonoWeight()*2); // whole mass

    if (alpha.hasNTerminalModification())
    {
      mono_weight -= alpha.getNTerminalModification()->getDiffMonoMass();
    }
    if (beta.hasCTerminalModification())
    {
      mono_weight -= beta.getCTerminalModification()->getDiffMonoMass();
    }

    mono_weight += Residue::getInternalToBIon().getMonoWeight(); // beta is B ion
    mono_weight += Residue::getInternalToYIon().getMonoWeight(); // alpha is Y ion

    temp_full_beta_mass = mono_weight;
    // subtract one residue at a time from alpha
    for (Size i = 0; i < link_pos_A-1 && i < alpha.size()-1; ++i)
    {
      // temp_full_beta_mass keeps track of the mass before residues of beta get subtracted
      mono_weight = temp_full_beta_mass;
      mono_weight -= alpha[i].getMonoWeight(Residue::Internal);
      temp_full_beta_mass = mono_weight;

      for (Size j = beta.size()-1; j > link_pos_B-1 && j > 0; --j)
      {
        mono_weight -= beta[j].getMonoWeight(Residue::Internal);
        double pos(mono_weight / static_cast<double>(charge));

        // if (pos < 0.0) { continue; }

        Peak1D p;
        p.setMZ(pos);
        p.setIntensity(intensity);
        spectrum.push_back(p);
        if (add_metainfo_)
        {
          ion_name = "[alpha|xi$y" + String(alpha.size()-i-1) + "+beta|xi$b" + String(j) + "]";
          ion_names.push_back(ion_name);
          charges.push_back(charge);
        }
        // TODO losses removed for these complex peaks for now
        // if (add_losses_)
        // {
        //   std::set< TheoreticalSpectrumGeneratorXLMS::LossMass > losses = backward_losses_alpha[i+1];
        //   losses.insert(forward_losses_beta[j-1].begin(), forward_losses_beta[j-1].end());
        //   addXLinkIonLossesFAST_(spectrum, charges, ion_names, mono_weight, intensity, charge, ion_name, losses);
        //   // addXLinkIonLosses_(spectrum, charges, ion_names, alpha.getSuffix(alpha.size()-i-1), beta.getPrefix(j-1), mono_weight, intensity, charge, ion_name);
        // }

        if (add_isotopes_ && max_isotope_ >= 2) // add second isotopic peak with fast method, if two or more peaks are asked for
        {
          pos += Constants::C13C12_MASSDIFF_U / static_cast<double>(charge);
          p.setMZ(pos);
          p.setIntensity(intensity);
          spectrum.push_back(p);
          if (add_metainfo_)
          {
            ion_names.push_back(ion_name);
            charges.push_back(charge);
          }
        }
      }
    }
    return;
  }


  void TheoreticalSpectrumGeneratorXLMS::addXLinkIonLosses_(PeakSpectrum & spectrum, DataArrays::IntegerDataArray& charges, DataArrays::StringDataArray& ion_names, const AASequence & ion1, const AASequence & ion2, double ion_mass, double intensity, int charge, String ion_name) const
  {
    Peak1D p;

    set<String> losses;
    for (AASequence::ConstIterator it = ion1.begin(); it != ion1.end(); ++it)
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
    for (AASequence::ConstIterator it = ion2.begin(); it != ion2.end(); ++it)
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

    p.setIntensity(intensity * rel_loss_intensity_);

    for (set<String>::const_iterator it = losses.begin(); it != losses.end(); ++it)
    {
      double loss_pos = ion_mass - EmpiricalFormula(*it).getMonoWeight();
      const String& loss_name = *it;
      String loss_ion_name;

      if (loss_pos < 0.0) { continue; }

      p.setMZ(loss_pos / static_cast<double>(charge));
      if (add_metainfo_)
      {
        // remove final bracket, insert loss name and add the bracket again
        loss_ion_name = ion_name.prefix(ion_name.size()-1) + "-" + loss_name + "]";
        ion_names.push_back(loss_ion_name);
        charges.push_back(charge);
      }
      spectrum.push_back(p);

      if (add_isotopes_ && max_isotope_ >= 2) // add second isotopic peak with fast method, if two or more peaks are asked for
      {
        p.setMZ((loss_pos + Constants::C13C12_MASSDIFF_U) / static_cast<double>(charge));
        spectrum.push_back(p);
        if (add_metainfo_)
        {
          ion_names.push_back(ion_name);
          charges.push_back(charge);
        }
      }
    }
  }

  void TheoreticalSpectrumGeneratorXLMS::addXLinkIonLossesFAST_(PeakSpectrum & spectrum, DataArrays::IntegerDataArray& charges, DataArrays::StringDataArray& ion_names, double mono_weight, double intensity, int charge, String ion_name, set< TheoreticalSpectrumGeneratorXLMS::LossMass > & losses) const
  {
    Peak1D p;

    p.setIntensity(intensity * rel_loss_intensity_);

    for (TheoreticalSpectrumGeneratorXLMS::LossMass loss : losses)
    {
      double mass_with_loss = mono_weight - loss.mass;
      String loss_ion_name;

      if (mass_with_loss < 0.0) { continue; }

      p.setMZ(mass_with_loss / static_cast<double>(charge));
      if (add_metainfo_)
      {
        // remove final bracket, insert loss name and add the bracket again
        loss_ion_name = ion_name.prefix(ion_name.size()-1) + "-" + loss.name + "]";
        ion_names.push_back(loss_ion_name);
        charges.push_back(charge);
      }
      spectrum.push_back(p);

      // TODO I think these are not really necessary for loss peaks
      // if (add_isotopes_ && max_isotope_ >= 2) // add second isotopic peak with fast method, if two or more peaks are asked for
      // {
      //   p.setMZ((mass_with_loss + Constants::C13C12_MASSDIFF_U) / static_cast<double>(charge));
      //   spectrum.push_back(p);
      //   if (add_metainfo_)
      //   {
      //     ion_names.push_back(ion_name);
      //     charges.push_back(charge);
      //   }
      // }
    }
  }

  void TheoreticalSpectrumGeneratorXLMS::getXLinkIonSpectrumWithLosses(PeakSpectrum & spectrum, OPXLDataStructs::ProteinProteinCrossLink & crosslink, bool frag_alpha, int mincharge, int maxcharge) const
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

    std::vector< std::set< TheoreticalSpectrumGeneratorXLMS::LossMass > > forward_losses;
    std::vector< std::set< TheoreticalSpectrumGeneratorXLMS::LossMass > > backward_losses;

    std::set< TheoreticalSpectrumGeneratorXLMS::LossMass > losses_peptide2;
    if (frag_alpha)
    {
      losses_peptide2 = getBackwardLossesForLinearIons_(crosslink.beta)[0];
      forward_losses = getForwardLossesForLinearIons_(crosslink.alpha);
      backward_losses = getBackwardLossesForLinearIons_(crosslink.alpha);
    }
    else
    {
      losses_peptide2 = getBackwardLossesForLinearIons_(crosslink.alpha)[0];
      forward_losses = getForwardLossesForLinearIons_(crosslink.beta);
      backward_losses = getBackwardLossesForLinearIons_(crosslink.beta);
    }

    for (Int z = mincharge; z <= maxcharge; ++z)
    {
      if (add_b_ions_)
      {
        addXLinkIonPeaksWithLosses_(spectrum, charges, ion_names, crosslink, frag_alpha, Residue::BIon, forward_losses, backward_losses, losses_peptide2, z);
      }
      if (add_y_ions_)
      {
        addXLinkIonPeaksWithLosses_(spectrum, charges, ion_names, crosslink, frag_alpha, Residue::YIon, forward_losses, backward_losses, losses_peptide2, z);
      }
      if (add_a_ions_)
      {
        addXLinkIonPeaksWithLosses_(spectrum, charges, ion_names, crosslink, frag_alpha, Residue::AIon, forward_losses, backward_losses, losses_peptide2, z);
      }
      if (add_x_ions_)
      {
        addXLinkIonPeaksWithLosses_(spectrum, charges, ion_names, crosslink, frag_alpha, Residue::XIon, forward_losses, backward_losses, losses_peptide2, z);
      }
      if (add_c_ions_)
      {
        addXLinkIonPeaksWithLosses_(spectrum, charges, ion_names, crosslink, frag_alpha, Residue::CIon, forward_losses, backward_losses, losses_peptide2, z);
      }
      if (add_z_ions_)
      {
        addXLinkIonPeaksWithLosses_(spectrum, charges, ion_names, crosslink, frag_alpha, Residue::ZIon, forward_losses, backward_losses, losses_peptide2, z);
      }
      if (add_k_linked_ions_)
      {
        double precursor_mass = crosslink.alpha.getMonoWeight() + crosslink.cross_linker_mass;
        if (!crosslink.beta.empty())
        {
          precursor_mass += crosslink.beta.getMonoWeight();
        }
        AASequence peptide;
        Size link_pos;
        if (frag_alpha)
        {
          peptide = crosslink.alpha;
          link_pos = crosslink.cross_link_position.first;
        } else
        {
          peptide = crosslink.beta;
          link_pos = crosslink.cross_link_position.second;
        }
        addKLinkedIonPeaks_(spectrum, charges, ion_names, peptide, link_pos, precursor_mass, frag_alpha, z);
      }
    }

    if (add_precursor_peaks_)
    {
      double precursor_mass = crosslink.alpha.getMonoWeight() + crosslink.cross_linker_mass;
      if (!crosslink.beta.empty())
      {
        precursor_mass += crosslink.beta.getMonoWeight();
      }
      addPrecursorPeaks_(spectrum, charges, ion_names, precursor_mass, maxcharge);
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

  void TheoreticalSpectrumGeneratorXLMS::addXLinkIonPeaksWithLosses_(PeakSpectrum & spectrum, DataArrays::IntegerDataArray & charges, DataArrays::StringDataArray & ion_names, OPXLDataStructs::ProteinProteinCrossLink & crosslink, bool frag_alpha, Residue::ResidueType res_type, std::vector< std::set< TheoreticalSpectrumGeneratorXLMS::LossMass > > & forward_losses, std::vector< std::set< TheoreticalSpectrumGeneratorXLMS::LossMass > > & backward_losses, std::set< TheoreticalSpectrumGeneratorXLMS::LossMass > & losses_peptide2, int charge) const
  {
    if (crosslink.alpha.empty())
    {
      cout << "Warning: Attempt at creating XLink Ions Spectrum from empty string!" << endl;
      return;
    }


    double precursor_mass = crosslink.alpha.getMonoWeight() + crosslink.cross_linker_mass;

    if (!crosslink.beta.empty())
    {
      precursor_mass += crosslink.beta.getMonoWeight();
    }

    String ion_type;
    AASequence peptide;
    AASequence peptide2;
    Size link_pos;
    if (frag_alpha)
    {
      ion_type = "alpha|xi";
      peptide = crosslink.alpha;
      peptide2 = crosslink.beta;
      link_pos = crosslink.cross_link_position.first;
    } else
    {
      ion_type = "beta|xi";
      peptide = crosslink.beta;
      peptide2 = crosslink.alpha;
      link_pos = crosslink.cross_link_position.second;
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
        for (Size i = peptide.size()-1; i > link_pos; --i)
        {
          mono_weight -= peptide[i].getMonoWeight(Residue::Internal);
          double pos(mono_weight / static_cast<double>(charge));
          int frag_index = i;

//          double recalc_pos = (peptide.getPrefix(frag_index+1).getMonoWeight() + static_cast<double>(charge) ) / static_cast<double>(charge);
//          cout << "Current residue: " << i << " = " << peptide.toUnmodifiedString()[i] << " | ion_type: " << ion_type << "$b" << frag_index << " | pos: " << pos << "| current_residue_mass: " << peptide[i].getMonoWeight(Residue::Internal) << " | recalc_pos: " << recalc_pos << endl;

          addPeak_(spectrum, charges, ion_names, pos, intensity, res_type, frag_index, charge, ion_type);
          String ion_name = "[" + ion_type + "$" + String(residueTypeToIonLetter_(res_type)) + String(frag_index) + "]";

          if ( !forward_losses.empty() && (!forward_losses[i-1].empty() || !losses_peptide2.empty()) )
          {
            std::set< TheoreticalSpectrumGeneratorXLMS::LossMass > losses = losses_peptide2;
            losses.insert(forward_losses[i-1].begin(), forward_losses[i-1].end());

            // std::set< TheoreticalSpectrumGeneratorXLMS::LossMass > losses;
            // std::set_union(losses_peptide2.begin(), losses_peptide2.end(), forward_losses[i-1].begin(), forward_losses[i-1].end(), losses.begin());

            addXLinkIonLossesFAST_(spectrum, charges, ion_names, mono_weight, intensity, charge, ion_name, losses);
          }
          // addXLinkIonLosses_(spectrum, charges, ion_names, peptide.getPrefix(i-1), peptide2, mono_weight, intensity, charge, ion_name);
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
          String ion_name = "[" + ion_type + "$" + String(residueTypeToIonLetter_(res_type)) + String(frag_index) + "]";

          if ( !backward_losses.empty() && (!backward_losses[i+1].empty() || !losses_peptide2.empty()) )
          {
            std::set< TheoreticalSpectrumGeneratorXLMS::LossMass > losses = losses_peptide2;
            losses.insert(backward_losses[i+1].begin(), backward_losses[i+1].end());

            // std::vector< TheoreticalSpectrumGeneratorXLMS::LossMass > losses;
            // std::set_union(losses_peptide2.begin(), losses_peptide2.end(), backward_losses[i+1].begin(), backward_losses[i+1].end(), losses.begin());

            addXLinkIonLossesFAST_(spectrum, charges, ion_names, mono_weight, intensity, charge, ion_name, losses);
          }
          // addXLinkIonLosses_(spectrum, charges, ion_names, peptide.getSuffix(peptide.size()-i-1), peptide2, mono_weight, intensity, charge, ion_name);

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

  std::vector< std::set< TheoreticalSpectrumGeneratorXLMS::LossMass > > TheoreticalSpectrumGeneratorXLMS::getForwardLossesForLinearIons_(AASequence & peptide) const
  {
    std::vector< std::set< TheoreticalSpectrumGeneratorXLMS::LossMass > > losses(peptide.size());
    for (Size i = 0; i < peptide.size(); ++i)
    {
      if (peptide[i].hasNeutralLoss())
      {
        vector<EmpiricalFormula> loss_formulas = peptide[i].getLossFormulas();
        for (Size k = 0; k != loss_formulas.size(); ++k)
        {
          TheoreticalSpectrumGeneratorXLMS::LossMass new_loss_mass;
          new_loss_mass.name = loss_formulas[k].toString();
          new_loss_mass.mass = loss_formulas[k].getMonoWeight();
          losses[i].insert(new_loss_mass);
        }
      }
    }

    // this gives us a "forward set" with incremental losses from the first to the last residue
    std::vector< std::set< TheoreticalSpectrumGeneratorXLMS::LossMass > > ion_losses(losses.size());
    ion_losses[0] = losses[0];
    for (Size i = 1; i < losses.size(); ++i)
    {
      std::set< TheoreticalSpectrumGeneratorXLMS::LossMass > new_set = ion_losses[i-1];
      new_set.insert(losses[i].begin(), losses[i].end());
      ion_losses[i] = new_set;
    }
    return ion_losses;
  }

  std::vector< std::set< TheoreticalSpectrumGeneratorXLMS::LossMass > > TheoreticalSpectrumGeneratorXLMS::getBackwardLossesForLinearIons_(AASequence & peptide) const
  {
    std::vector< std::set< TheoreticalSpectrumGeneratorXLMS::LossMass > > losses(peptide.size());
    for (Size i = 0; i < peptide.size(); ++i)
    {
      if (peptide[i].hasNeutralLoss())
      {
        vector<EmpiricalFormula> loss_formulas = peptide[i].getLossFormulas();
        for (Size k = 0; k != loss_formulas.size(); ++k)
        {
          String loss_name = loss_formulas[k].toString();
          if (loss_name == "H2O1" || loss_name == "H3N1")
          {
            TheoreticalSpectrumGeneratorXLMS::LossMass new_loss_mass;
            new_loss_mass.name = loss_formulas[k].toString();
            new_loss_mass.mass = loss_formulas[k].getMonoWeight();
            losses[i].insert(new_loss_mass);
          }
        }
      }
    }

    // this gives us a "backward set" with incremental losses from the last to the first residue
    std::vector< std::set< TheoreticalSpectrumGeneratorXLMS::LossMass > > ion_losses(losses.size());
    ion_losses[ion_losses.size()-1] = losses[losses.size()-1];
    for (Size i = ion_losses.size()-1; i > 0; --i)
    {
      std::set< TheoreticalSpectrumGeneratorXLMS::LossMass > new_set = ion_losses[i];
      new_set.insert(losses[i-1].begin(), losses[i-1].end());
      ion_losses[i-1] = new_set;
    }
    return ion_losses;
  }

}
