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
// $Maintainer: Eugen Netz $
// $Authors: Eugen Netz $
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/SimpleTSGXLMS.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/KERNEL/MSSpectrum.h>


using namespace std;

namespace OpenMS
{

  SimpleTSGXLMS::SimpleTSGXLMS() :
    DefaultParamHandler("SimpleTSGXLMS")
  {
    // TODO only partly functional (second isotopic peak if max_isotope = 2)
    defaults_.setValue("add_isotopes", "false", "If set to 1 isotope peaks of the product ion peaks are added");
    defaults_.setValidStrings("add_isotopes", ListUtils::create<String>("true,false"));

    defaults_.setValue("max_isotope", 2, "Defines the maximal isotopic peak which is added, add_isotopes must be set to 1");

    defaults_.setValue("add_charges", "true", "Adds the charges to a DataArray of the spectrum");
    defaults_.setValidStrings("add_charges", ListUtils::create<String>("true,false"));

    defaults_.setValue("add_losses", "false", "Adds common losses to those ion expect to have them, only water and ammonia loss is considered");
    defaults_.setValidStrings("add_losses", ListUtils::create<String>("true,false"));

    defaults_.setValue("add_precursor_peaks", "true", "Adds peaks of the precursor to the spectrum, which happen to occur sometimes");
    defaults_.setValidStrings("add_precursor_peaks", ListUtils::create<String>("true,false"));

    // TODO not functional yet
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

  SimpleTSGXLMS::SimpleTSGXLMS(const SimpleTSGXLMS & rhs) :
    DefaultParamHandler(rhs)
  {
  }

  SimpleTSGXLMS & SimpleTSGXLMS::operator=(const SimpleTSGXLMS & rhs)
  {
    if (this != &rhs)
    {
      DefaultParamHandler::operator=(rhs);
    }
    return *this;
  }

  SimpleTSGXLMS::~SimpleTSGXLMS()
  {
  }

  void SimpleTSGXLMS::getLinearIonSpectrum(PeakSpectrum & spectrum, AASequence & peptide, Size link_pos, int charge, Size link_pos_2) const
  {
    PeakSpectrum::IntegerDataArray charges;

    if (add_charges_)
    {
      if (spectrum.getIntegerDataArrays().size() > 0)
      {
        charges = spectrum.getIntegerDataArrays()[0];
      }
      charges.setName("Charges");
    }

    std::vector< std::set< SimpleTSGXLMS::LossMass, SimpleTSGXLMS::LossMassComparator > > forward_losses;
    std::vector< std::set< SimpleTSGXLMS::LossMass, SimpleTSGXLMS::LossMassComparator > > backward_losses;

    if (add_losses_)
    {
      forward_losses = getForwardLosses_(peptide);
      backward_losses = getBackwardLosses_(peptide);
    }

    for (Int z = 1; z <= charge; ++z)
    {
      if (add_b_ions_)
      {
        addLinearPeaks_(spectrum, charges, peptide, link_pos, Residue::BIon, forward_losses, backward_losses, z, link_pos_2);
      }
      if (add_y_ions_)
      {
        addLinearPeaks_(spectrum, charges, peptide, link_pos, Residue::YIon, forward_losses, backward_losses, z, link_pos_2);
      }
      if (add_a_ions_)
      {
        addLinearPeaks_(spectrum, charges, peptide, link_pos, Residue::AIon, forward_losses, backward_losses, z, link_pos_2);
      }
      if (add_x_ions_)
      {
        addLinearPeaks_(spectrum, charges, peptide, link_pos, Residue::XIon, forward_losses, backward_losses, z, link_pos_2);
      }
      if (add_c_ions_)
      {
        addLinearPeaks_(spectrum, charges, peptide, link_pos, Residue::CIon, forward_losses, backward_losses, z, link_pos_2);
      }
      if (add_z_ions_)
      {
        addLinearPeaks_(spectrum, charges, peptide, link_pos, Residue::ZIon, forward_losses, backward_losses, z, link_pos_2);
      }
    }

    if (add_charges_)
    {
      if (spectrum.getIntegerDataArrays().size() > 0)
      {
        spectrum.getIntegerDataArrays()[0] = charges;
      }
      else
      {
        spectrum.getIntegerDataArrays().push_back(charges);
      }
    }

    spectrum.sortByPosition();
    return;
  }

  void SimpleTSGXLMS::addLinearPeaks_(PeakSpectrum & spectrum, DataArrays::IntegerDataArray & charges, AASequence & peptide, Size link_pos, Residue::ResidueType res_type, std::vector< std::set< SimpleTSGXLMS::LossMass, SimpleTSGXLMS::LossMassComparator > > & forward_losses, std::vector< std::set< SimpleTSGXLMS::LossMass, SimpleTSGXLMS::LossMassComparator > > & backward_losses, int charge, Size link_pos_2) const
  {
    if (peptide.empty())
    {
      cout << "Warning: Attempt at creating XLink Ions Spectrum from empty string!" << endl;
      return;
    }

    // second link position, in case of a loop-link
    Size link_pos_B = link_pos_2;
    if (link_pos_2 == 0)
    {
      link_pos_B = link_pos;
    }

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

        addPeak_(spectrum, charges, pos, intensity, charge);
        if (add_losses_)
        {
          addLinearIonLosses_(spectrum, charges, mono_weight, intensity, charge, forward_losses[i]);
        }
        if (add_isotopes_ && max_isotope_ >= 2) // add second isotopic peak with fast method, if two or more peaks are asked for
        {
          pos += Constants::C13C12_MASSDIFF_U / static_cast<double>(charge);
          addPeak_(spectrum, charges, pos, intensity, charge);
        }
      }
    }
    else // if (res_type == Residue::XIon || res_type == Residue::YIon || res_type == Residue::ZIon)
    {
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

        addPeak_(spectrum, charges, pos, intensity, charge);
        if (add_losses_)
        {
          addLinearIonLosses_(spectrum, charges, pos, intensity, charge, backward_losses[i]);
        }
        if (add_isotopes_ && max_isotope_ >= 2) // add second isotopic peak with fast method, if two or more peaks are asked for
        {
          pos += Constants::C13C12_MASSDIFF_U / static_cast<double>(charge);
          addPeak_(spectrum, charges, pos, intensity, charge);
        }
      }
    }
    return;
  }

  void SimpleTSGXLMS::getXLinkIonSpectrum(PeakSpectrum & spectrum, AASequence & peptide, Size link_pos, double precursor_mass, int mincharge, int maxcharge, Size link_pos_2) const
  {
    PeakSpectrum::IntegerDataArray charges;

    if (add_charges_)
    {
      if (spectrum.getIntegerDataArrays().size() > 0)
      {
        charges = spectrum.getIntegerDataArrays()[0];
      }
      charges.setName("Charges");
    }

    std::vector< std::set< SimpleTSGXLMS::LossMass, SimpleTSGXLMS::LossMassComparator > > forward_losses;
    std::vector< std::set< SimpleTSGXLMS::LossMass, SimpleTSGXLMS::LossMassComparator > > backward_losses;

    if (add_losses_)
    {
      forward_losses = getForwardLosses_(peptide);
      backward_losses = getBackwardLosses_(peptide);
    }


    for (Int z = mincharge; z <= maxcharge; ++z)
    {
      if (add_b_ions_)
      {
        addXLinkIonPeaks_(spectrum, charges, peptide, link_pos, precursor_mass, Residue::BIon, forward_losses, backward_losses, z, link_pos_2);
      }
      if (add_y_ions_)
      {
        addXLinkIonPeaks_(spectrum, charges, peptide, link_pos, precursor_mass, Residue::YIon, forward_losses, backward_losses, z, link_pos_2);
      }
      if (add_a_ions_)
      {
        addXLinkIonPeaks_(spectrum, charges, peptide, link_pos, precursor_mass, Residue::AIon, forward_losses, backward_losses, z, link_pos_2);
      }
      if (add_x_ions_)
      {
        addXLinkIonPeaks_(spectrum, charges, peptide, link_pos, precursor_mass, Residue::XIon, forward_losses, backward_losses, z, link_pos_2);
      }
      if (add_c_ions_)
      {
        addXLinkIonPeaks_(spectrum, charges, peptide, link_pos, precursor_mass, Residue::CIon, forward_losses, backward_losses, z, link_pos_2);
      }
      if (add_z_ions_)
      {
        addXLinkIonPeaks_(spectrum, charges, peptide, link_pos, precursor_mass, Residue::ZIon, forward_losses, backward_losses, z, link_pos_2);
      }
      if (add_k_linked_ions_)
      {
        addKLinkedIonPeaks_(spectrum, charges, peptide, link_pos, precursor_mass, z);
      }
    }

    if (add_precursor_peaks_)
    {
      addPrecursorPeaks_(spectrum, charges, precursor_mass, maxcharge);
    }

    if (add_charges_)
    {
      if (spectrum.getIntegerDataArrays().size() > 0)
      {
        spectrum.getIntegerDataArrays()[0] = charges;
      }
      else
      {
        spectrum.getIntegerDataArrays().push_back(charges);
      }
    }

    spectrum.sortByPosition();
    return;
  }

  void SimpleTSGXLMS::addXLinkIonPeaks_(PeakSpectrum & spectrum, DataArrays::IntegerDataArray & charges, AASequence & peptide, Size link_pos, double precursor_mass, Residue::ResidueType res_type, std::vector< std::set< LossMass, LossMassComparator > > & forward_losses, std::vector< std::set< LossMass, LossMassComparator > > & backward_losses, int charge, Size link_pos_2) const
  {
    if (peptide.empty())
    {
      cout << "Warning: Attempt at creating XLink Ions Spectrum from empty string!" << endl;
      return;
    }

    // second link position, in case of a loop-link
    Size link_pos_B = link_pos_2;
    if (link_pos_2 == 0)
    {
      link_pos_B = link_pos;
    }

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

        addPeak_(spectrum, charges, pos, intensity, charge);
        if (add_losses_ && !forward_losses.empty() && (!forward_losses[i-1].empty()))
        {
          addXLinkIonLosses_(spectrum, charges, mono_weight, intensity, charge, forward_losses[i-1]);
        }

        if (add_isotopes_ && max_isotope_ >= 2) // add second isotopic peak with fast method, if two or more peaks are asked for
        {
          pos += Constants::C13C12_MASSDIFF_U / static_cast<double>(charge);
          addPeak_(spectrum, charges, pos, intensity, charge);
        }
      }
    }
    else // if (res_type == Residue::XIon || res_type == Residue::YIon || res_type == Residue::ZIon)
    {
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

        addPeak_(spectrum, charges, pos, intensity, charge);
        if (add_losses_ && !backward_losses.empty() && (!backward_losses[i+1].empty()))
        {
          addXLinkIonLosses_(spectrum, charges, mono_weight, intensity, charge, backward_losses[i+1]);
        }

        if (add_isotopes_ && max_isotope_ >= 2) // add second isotopic peak with fast method, if two or more peaks are asked for
        {
          pos += Constants::C13C12_MASSDIFF_U / static_cast<double>(charge);
          addPeak_(spectrum, charges, pos, intensity, charge);
        }
      }
    }
    return;
  }


  // helper to add a single peak to a spectrum (simple fragmentation)
  void SimpleTSGXLMS::addPeak_(PeakSpectrum & spectrum, DataArrays::IntegerDataArray & charges, double pos, double intensity, int charge) const
  {
    if (pos < 0)
    {
      return;
    }
    Peak1D p;
    p.setMZ(pos);
    p.setIntensity(intensity);
    spectrum.push_back(p);
    if (add_charges_)
    {
      charges.push_back(charge);
    }
  }

  void SimpleTSGXLMS::addLinearIonLosses_(PeakSpectrum & spectrum, DataArrays::IntegerDataArray& charges, double mono_weight, double intensity, int charge, std::set< SimpleTSGXLMS::LossMass, SimpleTSGXLMS::LossMassComparator > & losses) const
  {
    Peak1D p;
    p.setIntensity(intensity * rel_loss_intensity_);

    for (SimpleTSGXLMS::LossMass loss : losses)
    {
      double mass_with_loss = mono_weight - loss.mass;

      if (mass_with_loss < 0.0) { continue; }

      p.setMZ(mass_with_loss / static_cast<double>(charge));
      if (add_charges_)
      {
        charges.push_back(charge);
      }
      spectrum.push_back(p);
    }
  }

  void SimpleTSGXLMS::addPrecursorPeaks_(PeakSpectrum & spectrum, DataArrays::IntegerDataArray & charges, double precursor_mass, int charge) const
  {
    Peak1D p;

    // precursor peak
    double mono_pos = precursor_mass + (Constants::PROTON_MASS_U * static_cast<double>(charge));
    p.setMZ(mono_pos / static_cast<double>(charge));
    p.setIntensity(pre_int_);
    if (add_charges_)
    {
      charges.push_back(charge);
    }
    spectrum.push_back(p);
    if (add_isotopes_ && max_isotope_ >= 2) // add second isotopic peak with fast method, if two or more peaks are asked for
    {
      double pos = mono_pos + (Constants::C13C12_MASSDIFF_U / static_cast<double>(charge));
      p.setMZ(pos);
      p.setIntensity(pre_int_);
      if (add_charges_)
      {
        charges.push_back(charge);
      }
      spectrum.push_back(p);
    }

    // loss peaks of the precursor
    // loss of water
    mono_pos = precursor_mass + (Constants::PROTON_MASS_U * static_cast<double>(charge)) - EmpiricalFormula("H2O").getMonoWeight();
    p.setMZ(mono_pos / static_cast<double>(charge));
    p.setIntensity(pre_int_H2O_);
    if (add_charges_)
    {
      charges.push_back(charge);
    }
    spectrum.push_back(p);
    if (add_isotopes_ && max_isotope_ >= 2) // add second isotopic peak with fast method, if two or more peaks are asked for
    {
      double pos = mono_pos + (Constants::C13C12_MASSDIFF_U / static_cast<double>(charge));
      p.setMZ(pos);
      p.setIntensity(pre_int_H2O_);
      if (add_charges_)
      {
        charges.push_back(charge);
      }
      spectrum.push_back(p);
    }

    //loss of ammonia
    mono_pos = precursor_mass + (Constants::PROTON_MASS_U * static_cast<double>(charge)) - EmpiricalFormula("NH3").getMonoWeight();
    p.setMZ(mono_pos / static_cast<double>(charge));
    p.setIntensity(pre_int_NH3_);
    if (add_charges_)
    {
      charges.push_back(charge);
    }
    spectrum.push_back(p);
    if (add_isotopes_ && max_isotope_ >= 2) // add second isotopic peak with fast method, if two or more peaks are asked for
    {
      double pos = mono_pos + (Constants::C13C12_MASSDIFF_U / static_cast<double>(charge));
      p.setMZ(pos);
      p.setIntensity(pre_int_NH3_);
      if (add_charges_)
      {
        charges.push_back(charge);
      }
      spectrum.push_back(p);
    }
  }

  void SimpleTSGXLMS::addKLinkedIonPeaks_(PeakSpectrum & spectrum, DataArrays::IntegerDataArray & charges, AASequence & peptide, Size link_pos, double precursor_mass, int charge) const
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
    if (mono_weight < 0)
    {
      return;
    }

    double pos(mono_weight / static_cast<double>(charge));

    Peak1D p;
    p.setMZ(pos);
    p.setIntensity(1.0);
    spectrum.push_back(p);

    if (add_charges_)
    {
      charges.push_back(charge);
    }

    if (add_isotopes_ && max_isotope_ >= 2) // add second isotopic peak with fast method, if two or more peaks are asked for
    {
      pos += Constants::C13C12_MASSDIFF_U / static_cast<double>(charge);
      p.setMZ(pos);
      spectrum.push_back(p);
      if (add_charges_)
      {
        charges.push_back(charge);
      }
    }
  }

  void SimpleTSGXLMS::addXLinkIonLosses_(PeakSpectrum & spectrum, DataArrays::IntegerDataArray& charges, double mono_weight, double intensity, int charge, set< SimpleTSGXLMS::LossMass, SimpleTSGXLMS::LossMassComparator > & losses) const
  {
    Peak1D p;
    p.setIntensity(intensity * rel_loss_intensity_);

    for (SimpleTSGXLMS::LossMass loss : losses)
    {
      double mass_with_loss = mono_weight - loss.mass;

      if (mass_with_loss < 0.0) { continue; }

      p.setMZ(mass_with_loss / static_cast<double>(charge));
      if (add_charges_)
      {
        charges.push_back(charge);
      }
      spectrum.push_back(p);
    }
  }

  void SimpleTSGXLMS::getXLinkIonSpectrum(PeakSpectrum & spectrum, OPXLDataStructs::ProteinProteinCrossLink & crosslink, bool frag_alpha, int mincharge, int maxcharge) const
  {
    PeakSpectrum::IntegerDataArray charges;

    if (add_charges_)
    {
      if (spectrum.getIntegerDataArrays().size() > 0)
      {
        charges = spectrum.getIntegerDataArrays()[0];
      }
      charges.setName("Charges");
    }

    std::vector< std::set< SimpleTSGXLMS::LossMass, SimpleTSGXLMS::LossMassComparator > > forward_losses;
    std::vector< std::set< SimpleTSGXLMS::LossMass, SimpleTSGXLMS::LossMassComparator > > backward_losses;
    std::set< SimpleTSGXLMS::LossMass, SimpleTSGXLMS::LossMassComparator > losses_peptide2;

    if (add_losses_)
    {
      if (frag_alpha)
      {
        losses_peptide2 = getBackwardLosses_(crosslink.beta)[0];
        forward_losses = getForwardLosses_(crosslink.alpha);
        backward_losses = getBackwardLosses_(crosslink.alpha);
      }
      else
      {
        losses_peptide2 = getBackwardLosses_(crosslink.alpha)[0];
        forward_losses = getForwardLosses_(crosslink.beta);
        backward_losses = getBackwardLosses_(crosslink.beta);
      }
    }

    for (Int z = mincharge; z <= maxcharge; ++z)
    {
      if (add_b_ions_)
      {
        addXLinkIonPeaks_(spectrum, charges, crosslink, frag_alpha, Residue::BIon, forward_losses, backward_losses, losses_peptide2, z);
      }
      if (add_y_ions_)
      {
        addXLinkIonPeaks_(spectrum, charges, crosslink, frag_alpha, Residue::YIon, forward_losses, backward_losses, losses_peptide2, z);
      }
      if (add_a_ions_)
      {
        addXLinkIonPeaks_(spectrum, charges, crosslink, frag_alpha, Residue::AIon, forward_losses, backward_losses, losses_peptide2, z);
      }
      if (add_x_ions_)
      {
        addXLinkIonPeaks_(spectrum, charges, crosslink, frag_alpha, Residue::XIon, forward_losses, backward_losses, losses_peptide2, z);
      }
      if (add_c_ions_)
      {
        addXLinkIonPeaks_(spectrum, charges, crosslink, frag_alpha, Residue::CIon, forward_losses, backward_losses, losses_peptide2, z);
      }
      if (add_z_ions_)
      {
        addXLinkIonPeaks_(spectrum, charges, crosslink, frag_alpha, Residue::ZIon, forward_losses, backward_losses, losses_peptide2, z);
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
        }
        else
        {
          peptide = crosslink.beta;
          link_pos = crosslink.cross_link_position.second;
        }
        addKLinkedIonPeaks_(spectrum, charges, peptide, link_pos, precursor_mass, z);
      }
    }

    if (add_precursor_peaks_)
    {
      double precursor_mass = crosslink.alpha.getMonoWeight() + crosslink.cross_linker_mass;
      if (!crosslink.beta.empty())
      {
        precursor_mass += crosslink.beta.getMonoWeight();
      }
      addPrecursorPeaks_(spectrum, charges, precursor_mass, maxcharge);
    }

    if (add_charges_)
    {
      if (spectrum.getIntegerDataArrays().size() > 0)
      {
        spectrum.getIntegerDataArrays()[0] = charges;
      }
      else
      {
        spectrum.getIntegerDataArrays().push_back(charges);
      }
    }

    spectrum.sortByPosition();
    return;
  }

  void SimpleTSGXLMS::addXLinkIonPeaks_(PeakSpectrum & spectrum, DataArrays::IntegerDataArray & charges, OPXLDataStructs::ProteinProteinCrossLink & crosslink, bool frag_alpha, Residue::ResidueType res_type, std::vector< std::set< SimpleTSGXLMS::LossMass, SimpleTSGXLMS::LossMassComparator > > & forward_losses, std::vector< std::set< SimpleTSGXLMS::LossMass, SimpleTSGXLMS::LossMassComparator > > & backward_losses, std::set< SimpleTSGXLMS::LossMass, SimpleTSGXLMS::LossMassComparator > & losses_peptide2, int charge) const
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

    AASequence peptide;
    AASequence peptide2;
    Size link_pos;
    if (frag_alpha)
    {
      peptide = crosslink.alpha;
      peptide2 = crosslink.beta;
      link_pos = crosslink.cross_link_position.first;
    }
    else
    {
      peptide = crosslink.beta;
      peptide2 = crosslink.alpha;
      link_pos = crosslink.cross_link_position.second;
    }

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

        addPeak_(spectrum, charges, pos, intensity, charge);
        if (add_losses_ && !forward_losses.empty() && (!forward_losses[i-1].empty() || !losses_peptide2.empty()) )
        {
          std::set< SimpleTSGXLMS::LossMass, SimpleTSGXLMS::LossMassComparator > losses = losses_peptide2;
          losses.insert(forward_losses[i-1].begin(), forward_losses[i-1].end());
          addXLinkIonLosses_(spectrum, charges, mono_weight, intensity, charge, losses);
        }
        if (add_isotopes_ && max_isotope_ >= 2) // add second isotopic peak with fast method, if two or more peaks are asked for
        {
          pos += Constants::C13C12_MASSDIFF_U / static_cast<double>(charge);
          addPeak_(spectrum, charges, pos, intensity, charge);
        }
      }
    }
    else // if (res_type == Residue::XIon || res_type == Residue::YIon || res_type == Residue::ZIon)
    {
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

        addPeak_(spectrum, charges, pos, intensity, charge);
        if (add_losses_ && !backward_losses.empty() && (!backward_losses[i+1].empty() || !losses_peptide2.empty()) )
        {
          std::set< SimpleTSGXLMS::LossMass, SimpleTSGXLMS::LossMassComparator > losses = losses_peptide2;
          losses.insert(backward_losses[i+1].begin(), backward_losses[i+1].end());

          addXLinkIonLosses_(spectrum, charges, mono_weight, intensity, charge, losses);
        }

        if (add_isotopes_ && max_isotope_ >= 2) // add second isotopic peak with fast method, if two or more peaks are asked for
        {
          pos += Constants::C13C12_MASSDIFF_U / static_cast<double>(charge);
          addPeak_(spectrum, charges, pos, intensity, charge);
        }
      }
    }
    return;
  }

  std::vector< std::set< SimpleTSGXLMS::LossMass, SimpleTSGXLMS::LossMassComparator > > SimpleTSGXLMS::getForwardLosses_(AASequence & peptide) const
  {
    std::vector< std::set< SimpleTSGXLMS::LossMass, SimpleTSGXLMS::LossMassComparator > > losses(peptide.size());
    for (Size i = 0; i < peptide.size(); ++i)
    {
      if (peptide[i].hasNeutralLoss())
      {
        vector<EmpiricalFormula> loss_formulas = peptide[i].getLossFormulas();
        for (Size k = 0; k != loss_formulas.size(); ++k)
        {
          String loss_name = loss_formulas[k].toString();
          if (loss_name == "H2O1" || loss_name == "H3N1") // for now only these most common losses are considered
          {
            SimpleTSGXLMS::LossMass new_loss_mass;
            new_loss_mass.name = loss_formulas[k].toString();
            new_loss_mass.mass = loss_formulas[k].getMonoWeight();
            losses[i].insert(new_loss_mass);
          }
        }
      }
    }

    // this gives us a "forward set" with incremental losses from the first to the last residue
    std::vector< std::set< SimpleTSGXLMS::LossMass, SimpleTSGXLMS::LossMassComparator > > ion_losses(losses.size());
    ion_losses[0] = losses[0];
    for (Size i = 1; i < losses.size(); ++i)
    {
      std::set< SimpleTSGXLMS::LossMass, SimpleTSGXLMS::LossMassComparator > new_set = ion_losses[i-1];
      new_set.insert(losses[i].begin(), losses[i].end());
      ion_losses[i] = new_set;
    }
    return ion_losses;
  }

  std::vector< std::set< SimpleTSGXLMS::LossMass, SimpleTSGXLMS::LossMassComparator > > SimpleTSGXLMS::getBackwardLosses_(AASequence & peptide) const
  {
    std::vector< std::set< SimpleTSGXLMS::LossMass, SimpleTSGXLMS::LossMassComparator > > losses(peptide.size());
    for (Size i = 0; i < peptide.size(); ++i)
    {
      if (peptide[i].hasNeutralLoss())
      {
        vector<EmpiricalFormula> loss_formulas = peptide[i].getLossFormulas();
        for (Size k = 0; k != loss_formulas.size(); ++k)
        {
          String loss_name = loss_formulas[k].toString();
          if (loss_name == "H2O1" || loss_name == "H3N1") // for now only these most common losses are considered
          {
            SimpleTSGXLMS::LossMass new_loss_mass;
            new_loss_mass.name = loss_formulas[k].toString();
            new_loss_mass.mass = loss_formulas[k].getMonoWeight();
            losses[i].insert(new_loss_mass);
          }
        }
      }
    }

    // this gives us a "backward set" with incremental losses from the last to the first residue
    std::vector< std::set< SimpleTSGXLMS::LossMass, SimpleTSGXLMS::LossMassComparator > > ion_losses(losses.size());
    ion_losses[ion_losses.size()-1] = losses[losses.size()-1];
    for (Size i = ion_losses.size()-1; i > 0; --i)
    {
      std::set< SimpleTSGXLMS::LossMass, SimpleTSGXLMS::LossMassComparator > new_set = ion_losses[i];
      new_set.insert(losses[i-1].begin(), losses[i-1].end());
      ion_losses[i-1] = new_set;
    }
    return ion_losses;
  }

  // helper for mapping residue type to letter
  char SimpleTSGXLMS::residueTypeToIonLetter_(Residue::ResidueType res_type) const
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

  void SimpleTSGXLMS::updateMembers_()
  {
    add_b_ions_ = param_.getValue("add_b_ions").toBool();
    add_y_ions_ = param_.getValue("add_y_ions").toBool();
    add_a_ions_ = param_.getValue("add_a_ions").toBool();
    add_c_ions_ = param_.getValue("add_c_ions").toBool();
    add_x_ions_ = param_.getValue("add_x_ions").toBool();
    add_z_ions_ = param_.getValue("add_z_ions").toBool();
    add_first_prefix_ion_ = param_.getValue("add_first_prefix_ion").toBool();
    add_losses_ = param_.getValue("add_losses").toBool();
    add_charges_ = param_.getValue("add_charges").toBool();
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
}
