// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Eugen Netz $
// $Authors: Eugen Netz $
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGeneratorXLMS.h>
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
    defaults_.setValidStrings("add_isotopes", {"true","false"});

    defaults_.setValue("max_isotope", 2, "Defines the maximal isotopic peak which is added, add_isotopes must be set to 1");

    defaults_.setValue("add_metainfo", "true", "Adds the type of peaks as metainfo to the peaks, like y8+, [M-H2O+2H]++");
    defaults_.setValidStrings("add_metainfo", {"true","false"});

    defaults_.setValue("add_charges", "true", "Adds the charges to a DataArray of the spectrum");
    defaults_.setValidStrings("add_charges", {"true","false"});

    defaults_.setValue("add_losses", "false", "Adds common losses to those ion expect to have them, only water and ammonia loss is considered");
    defaults_.setValidStrings("add_losses", {"true","false"});

    defaults_.setValue("add_precursor_peaks", "true", "Adds peaks of the precursor to the spectrum, which happen to occur sometimes");
    defaults_.setValidStrings("add_precursor_peaks", {"true","false"});

    // TODO not functional yet
    defaults_.setValue("add_abundant_immonium_ions", "false", "Add most abundant immonium ions");
    defaults_.setValidStrings("add_abundant_immonium_ions", {"true","false"});

    defaults_.setValue("add_k_linked_ions", "true", "Add RES-Linked ions, which are specific to XLMS");
    defaults_.setValidStrings("add_k_linked_ions", {"true","false"});

    // TODO not functional yet
    defaults_.setValue("add_first_prefix_ion", "true", "If set to true e.g. b1 ions are added");
    defaults_.setValidStrings("add_first_prefix_ion", {"true","false"});

    defaults_.setValue("add_y_ions", "true", "Add peaks of y-ions to the spectrum");
    defaults_.setValidStrings("add_y_ions", {"true","false"});

    defaults_.setValue("add_b_ions", "true", "Add peaks of b-ions to the spectrum");
    defaults_.setValidStrings("add_b_ions", {"true","false"});

    defaults_.setValue("add_a_ions", "true", "Add peaks of a-ions to the spectrum");
    defaults_.setValidStrings("add_a_ions", {"true","false"});

    defaults_.setValue("add_c_ions", "false", "Add peaks of c-ions to the spectrum");
    defaults_.setValidStrings("add_c_ions", {"true","false"});

    defaults_.setValue("add_x_ions", "false", "Add peaks of  x-ions to the spectrum");
    defaults_.setValidStrings("add_x_ions", {"true","false"});

    defaults_.setValue("add_z_ions", "false", "Add peaks of z-ions to the spectrum");
    defaults_.setValidStrings("add_z_ions", {"true","false"});


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

    // preprocess loss_db_, a database of H2O and NH3 losses for all residues
    AASequence residues = AASequence::fromString("RHKDESTNQCUGPAVILMFYW");
    for (Size i = 0; i < residues.size(); ++i)
    {
      LossIndex residue_losses;
      loss_db_.insert(std::make_pair(residues[i].getOneLetterCode(), residue_losses));
      if (residues[i].hasNeutralLoss())
      {
        vector<EmpiricalFormula> loss_formulas = residues[i].getLossFormulas();
        for (Size k = 0; k != loss_formulas.size(); ++k)
        {
          String loss_name = loss_formulas[k].toString();
          if (loss_name == "H2O1") // for now only these most common losses are considered
          {
            if (loss_H2O_ < 1)
            {
              loss_H2O_ = loss_formulas[k].getMonoWeight();
            }
            loss_db_[residues[i].getOneLetterCode()].has_H2O_loss = true;
          }

          if (loss_name == "H3N1")
          {
            if (loss_NH3_ < 1)
            {
              loss_NH3_ = loss_formulas[k].getMonoWeight();
            }
            loss_db_[residues[i].getOneLetterCode()].has_NH3_loss = true;
          }
        }
      }
    }
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

  TheoreticalSpectrumGeneratorXLMS::~TheoreticalSpectrumGeneratorXLMS() = default;

  void TheoreticalSpectrumGeneratorXLMS::getLinearIonSpectrum(PeakSpectrum & spectrum, AASequence & peptide, Size link_pos, bool frag_alpha, int charge, Size link_pos_2) const
  {
    PeakSpectrum::IntegerDataArray charges;
    PeakSpectrum::StringDataArray ion_names;

    if (add_charges_)
    {
      if (!spectrum.getIntegerDataArrays().empty())
      {
        charges = spectrum.getIntegerDataArrays()[0];
      }
      charges.setName("charge");
    }
    if (add_metainfo_)
    {
      if (!spectrum.getStringDataArrays().empty())
      {
        ion_names = spectrum.getStringDataArrays()[0];
      }
      ion_names.setName(Constants::UserParam::IonNames);
    }

    std::vector< LossIndex > forward_losses;
    std::vector< LossIndex > backward_losses;

    if (add_losses_)
    {
      forward_losses = getForwardLosses_(peptide);
      backward_losses = getBackwardLosses_(peptide);
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

    if (add_charges_)
    {
      if (!spectrum.getIntegerDataArrays().empty())
      {
        spectrum.getIntegerDataArrays()[0] = charges;
      }
      else
      {
        spectrum.getIntegerDataArrays().push_back(charges);
      }
    }
    if (add_metainfo_)
    {
      if (!spectrum.getStringDataArrays().empty())
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

  void TheoreticalSpectrumGeneratorXLMS::addLinearPeaks_(PeakSpectrum & spectrum, DataArrays::IntegerDataArray & charges, DataArrays::StringDataArray & ion_names, AASequence & peptide, Size link_pos, bool frag_alpha, Residue::ResidueType res_type, std::vector< LossIndex > & forward_losses, std::vector< LossIndex > & backward_losses, int charge, Size link_pos_2) const
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
        int frag_index = i+1;

        addPeak_(spectrum, charges, ion_names, pos, intensity, res_type, frag_index, charge, ion_type);
        if (add_losses_)
        {
          addLinearIonLosses_(spectrum, charges, ion_names, mono_weight, res_type, frag_index, intensity, charge, ion_type, forward_losses[i]);
        }
        if (add_isotopes_ && max_isotope_ >= 2) // add second isotopic peak with fast method, if two or more peaks are asked for
        {
          pos += Constants::C13C12_MASSDIFF_U / static_cast<double>(charge);
          addPeak_(spectrum, charges, ion_names, pos, intensity, res_type, frag_index, charge, ion_type);
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
        int frag_index = peptide.size() - i;

        addPeak_(spectrum, charges, ion_names, pos, intensity, res_type, frag_index, charge, ion_type);
        if (add_losses_)
        {
          addLinearIonLosses_(spectrum, charges, ion_names, pos, res_type, frag_index, intensity, charge, ion_type, backward_losses[i]);
        }
        if (add_isotopes_ && max_isotope_ >= 2) // add second isotopic peak with fast method, if two or more peaks are asked for
        {
          pos += Constants::C13C12_MASSDIFF_U / static_cast<double>(charge);
          addPeak_(spectrum, charges, ion_names, pos, intensity, res_type, frag_index, charge, ion_type);
        }
      }
    }
    return;
  }

  void TheoreticalSpectrumGeneratorXLMS::getXLinkIonSpectrum(PeakSpectrum & spectrum, AASequence & peptide, Size link_pos, double precursor_mass, bool frag_alpha, int mincharge, int maxcharge, Size link_pos_2) const
  {
    PeakSpectrum::IntegerDataArray charges;
    PeakSpectrum::StringDataArray ion_names;

    if (add_charges_)
    {
      if (!spectrum.getIntegerDataArrays().empty())
      {
        charges = spectrum.getIntegerDataArrays()[0];
      }
      charges.setName("charge");
    }
    if (add_metainfo_)
    {
      if (!spectrum.getStringDataArrays().empty())
      {
        ion_names = spectrum.getStringDataArrays()[0];
      }
      ion_names.setName(Constants::UserParam::IonNames);
    }

    std::vector< LossIndex > forward_losses;
    std::vector< LossIndex > backward_losses;

    if (add_losses_)
    {
      forward_losses = getForwardLosses_(peptide);
      backward_losses = getBackwardLosses_(peptide);
    }


    for (Int z = mincharge; z <= maxcharge; ++z)
    {
      if (add_b_ions_)
      {
        addXLinkIonPeaks_(spectrum, charges, ion_names, peptide, link_pos, precursor_mass, frag_alpha, Residue::BIon, forward_losses, backward_losses, z, link_pos_2);
      }
      if (add_y_ions_)
      {
        addXLinkIonPeaks_(spectrum, charges, ion_names, peptide, link_pos, precursor_mass, frag_alpha, Residue::YIon, forward_losses, backward_losses, z, link_pos_2);
      }
      if (add_a_ions_)
      {
        addXLinkIonPeaks_(spectrum, charges, ion_names, peptide, link_pos, precursor_mass, frag_alpha, Residue::AIon, forward_losses, backward_losses, z, link_pos_2);
      }
      if (add_x_ions_)
      {
        addXLinkIonPeaks_(spectrum, charges, ion_names, peptide, link_pos, precursor_mass, frag_alpha, Residue::XIon, forward_losses, backward_losses, z, link_pos_2);
      }
      if (add_c_ions_)
      {
        addXLinkIonPeaks_(spectrum, charges, ion_names, peptide, link_pos, precursor_mass, frag_alpha, Residue::CIon, forward_losses, backward_losses, z, link_pos_2);
      }
      if (add_z_ions_)
      {
        addXLinkIonPeaks_(spectrum, charges, ion_names, peptide, link_pos, precursor_mass, frag_alpha, Residue::ZIon, forward_losses, backward_losses, z, link_pos_2);
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

    if (add_charges_)
    {
      if (!spectrum.getIntegerDataArrays().empty())
      {
        spectrum.getIntegerDataArrays()[0] = charges;
      }
      else
      {
        spectrum.getIntegerDataArrays().push_back(charges);
      }
    }
    if (add_metainfo_)
    {
      if (!spectrum.getStringDataArrays().empty())
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

  void TheoreticalSpectrumGeneratorXLMS::addXLinkIonPeaks_(PeakSpectrum & spectrum, DataArrays::IntegerDataArray & charges, DataArrays::StringDataArray & ion_names, AASequence & peptide, Size link_pos, double precursor_mass, bool frag_alpha, Residue::ResidueType res_type, std::vector< LossIndex > & forward_losses, std::vector< LossIndex > & backward_losses, int charge, Size link_pos_2) const
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
    }
    else
    {
      ion_type = "beta|xi";
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
        int frag_index = i;

        addPeak_(spectrum, charges, ion_names, pos, intensity, res_type, frag_index, charge, ion_type);
        if (add_losses_ && forward_losses.size() >= i)
        {
          String ion_name = "[" + ion_type + "$" + String(Residue::residueTypeToIonLetter(res_type)) + String(frag_index) + "]";
          addXLinkIonLosses_(spectrum, charges, ion_names, mono_weight, intensity, charge, ion_name, forward_losses[i-1]);
        }

        if (add_isotopes_ && max_isotope_ >= 2) // add second isotopic peak with fast method, if two or more peaks are asked for
        {
          pos += Constants::C13C12_MASSDIFF_U / static_cast<double>(charge);
          addPeak_(spectrum, charges, ion_names, pos, intensity, res_type, frag_index, charge, ion_type);
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
        int frag_index = peptide.size() - 1 - i;

        addPeak_(spectrum, charges, ion_names, pos, intensity, res_type, frag_index, charge, ion_type);
        if (add_losses_ && backward_losses.size() >= i+2)
        {
          String ion_name = "[" + ion_type + "$" + String(Residue::residueTypeToIonLetter(res_type)) + String(frag_index) + "]";
          addXLinkIonLosses_(spectrum, charges, ion_names, mono_weight, intensity, charge, ion_name, backward_losses[i+1]);
        }

        if (add_isotopes_ && max_isotope_ >= 2) // add second isotopic peak with fast method, if two or more peaks are asked for
        {
          pos += Constants::C13C12_MASSDIFF_U / static_cast<double>(charge);
          addPeak_(spectrum, charges, ion_names, pos, intensity, res_type, frag_index, charge, ion_type);
        }
      }
    }
    return;
  }


  // helper to add a single peak to a spectrum (simple fragmentation)
  void TheoreticalSpectrumGeneratorXLMS::addPeak_(PeakSpectrum & spectrum, DataArrays::IntegerDataArray & charges, DataArrays::StringDataArray & ion_names, double pos, double intensity, Residue::ResidueType res_type, Size frag_index, int charge, String ion_type) const
  {
    if (pos < 0) {return;}

    Peak1D p;
    p.setMZ(pos);
    p.setIntensity(intensity);
    spectrum.push_back(p);
    if (add_metainfo_)
    {
      ion_names.emplace_back("[" + ion_type + "$" + String(Residue::residueTypeToIonLetter(res_type)) + String(frag_index) + "]");
    }
    if (add_charges_)
    {
      charges.push_back(charge);
    }
  }

  void TheoreticalSpectrumGeneratorXLMS::addLinearIonLosses_(PeakSpectrum & spectrum, DataArrays::IntegerDataArray& charges, DataArrays::StringDataArray& ion_names, double mono_weight, Residue::ResidueType res_type, Size frag_index, double intensity, int charge, String ion_type, LossIndex & losses) const
  {
    Peak1D p;
    p.setIntensity(intensity * rel_loss_intensity_);

    if (losses.has_H2O_loss)
    {
      double mass_with_loss = mono_weight - loss_H2O_;
      if (mass_with_loss > 0.0)
      {
        p.setMZ(mass_with_loss / static_cast<double>(charge));
        if (add_metainfo_)
        {
          // remove final bracket, insert loss name and add the bracket again
          ion_names.emplace_back("[" + ion_type + "$" + String(Residue::residueTypeToIonLetter(res_type)) + String(frag_index) + "-H2O1]");
        }
        if (add_charges_)
        {
          charges.push_back(charge);
        }
        spectrum.push_back(p);
      }
    }

    if (losses.has_NH3_loss)
    {
      double mass_with_loss = mono_weight - loss_NH3_;
      if (mass_with_loss > 0.0)
      {
        p.setMZ(mass_with_loss / static_cast<double>(charge));
        if (add_metainfo_)
        {
          // remove final bracket, insert loss name and add the bracket again
          ion_names.emplace_back("[" + ion_type + "$" + String(Residue::residueTypeToIonLetter(res_type)) + String(frag_index) + "-H3N1]");
        }
        if (add_charges_)
        {
          charges.push_back(charge);
        }
        spectrum.push_back(p);
      }
    }
  }

  void TheoreticalSpectrumGeneratorXLMS::addPrecursorPeaks_(PeakSpectrum & spectrum, DataArrays::IntegerDataArray & charges, DataArrays::StringDataArray & ion_names, double precursor_mass, int charge) const
  {
    Peak1D p;

    // precursor peak
    double mono_pos = precursor_mass + (Constants::PROTON_MASS_U * static_cast<double>(charge));
    p.setMZ(mono_pos / static_cast<double>(charge));
    p.setIntensity(pre_int_);
    if (add_metainfo_)
    {
      ion_names.emplace_back("[M+H]");
    }
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
      if (add_metainfo_)
      {
        ion_names.emplace_back("[M+H]");
      }
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
    if (add_metainfo_)
    {
      ion_names.emplace_back("[M+H]-H2O");
    }
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
      if (add_metainfo_)
      {
        ion_names.emplace_back("[M+H]-H2O");
      }
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
    if (add_metainfo_)
    {
      ion_names.emplace_back("[M+H]-NH3");
    }
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
      if (add_metainfo_)
      {
        ion_names.emplace_back("[M+H]-NH3");
      }
      if (add_charges_)
      {
        charges.push_back(charge);
      }
      spectrum.push_back(p);
    }
  }

  void TheoreticalSpectrumGeneratorXLMS::addKLinkedIonPeaks_(PeakSpectrum & spectrum, DataArrays::IntegerDataArray & charges, DataArrays::StringDataArray & ion_names, AASequence & peptide, Size link_pos, double precursor_mass, bool frag_alpha, int charge) const
  {
    double mono_weight = precursor_mass;
    // link_pos can be zero, if the cross-link is N-terminal
    if (link_pos > 0)
    {
      mono_weight -= peptide.getPrefix(link_pos).getMonoWeight(Residue::BIon);
    }
    else
    {
      return; // this fragment type is not necessary for links on peptide terminal residues
    }
    // same here for C-terminal links
    if (link_pos < peptide.size())
    {
      mono_weight -= peptide.getSuffix(peptide.size() - link_pos - 1).getMonoWeight(Residue::XIon);
    }
    else
    {
      return;
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

    // here the ion type is reversed compared to other peak types,
    // because for this special ion type, it would not make sense to call it alpha$y(n)-alpha$a(n)
    // Only one residue is left of the fragmented Peptide, so we call it a RES-linked beta
    String ion_type;
    String ion_name;

    if (add_metainfo_)
    {
      if (frag_alpha)
      {
        ion_type = "beta";
      }
      else
      {
        ion_type = "alpha";
      }

      int l_pos = link_pos;
      if (l_pos < 1)
      {
        l_pos = 0;
      }
      ion_name = "[" + peptide[l_pos].getOneLetterCode() + "-linked-" + ion_type + "]";
      ion_names.push_back(ion_name);
    }
    if (add_charges_)
    {
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
      }
      if (add_charges_)
      {
        charges.push_back(charge);
      }
    }
  }

  void TheoreticalSpectrumGeneratorXLMS::addXLinkIonLosses_(PeakSpectrum & spectrum, DataArrays::IntegerDataArray& charges, DataArrays::StringDataArray& ion_names, double mono_weight, double intensity, int charge, String ion_name, LossIndex & losses) const
  {
    Peak1D p;
    p.setIntensity(intensity * rel_loss_intensity_);

    if (losses.has_H2O_loss)
    {
      double mass_with_loss = mono_weight - loss_H2O_;
      if (mass_with_loss > 0.0)
      {
        p.setMZ(mass_with_loss / static_cast<double>(charge));
        if (add_metainfo_)
        {
          // remove final bracket, insert loss name and add the bracket again
          ion_names.emplace_back(ion_name.prefix(ion_name.size()-1) + "-H2O1]");
        }
        if (add_charges_)
        {
          charges.push_back(charge);
        }
        spectrum.push_back(p);
      }
    }

    if (losses.has_NH3_loss)
    {
      double mass_with_loss = mono_weight - loss_NH3_;
      if (mass_with_loss > 0.0)
      {
        p.setMZ(mass_with_loss / static_cast<double>(charge));
        if (add_metainfo_)
        {
          // remove final bracket, insert loss name and add the bracket again
          ion_names.emplace_back(ion_name.prefix(ion_name.size()-1) + "-H3N1]");
        }
        if (add_charges_)
        {
          charges.push_back(charge);
        }
        spectrum.push_back(p);
      }
    }
  }

  void TheoreticalSpectrumGeneratorXLMS::getXLinkIonSpectrum(PeakSpectrum & spectrum, OPXLDataStructs::ProteinProteinCrossLink & crosslink, bool frag_alpha, int mincharge, int maxcharge) const
  {
    PeakSpectrum::IntegerDataArray charges;
    PeakSpectrum::StringDataArray ion_names;

    if (add_charges_)
    {
      if (!spectrum.getIntegerDataArrays().empty())
      {
        charges = spectrum.getIntegerDataArrays()[0];
      }
      charges.setName("charge");
    }
    if (add_metainfo_)
    {
      if (!spectrum.getStringDataArrays().empty())
      {
        ion_names = spectrum.getStringDataArrays()[0];
      }
      ion_names.setName(Constants::UserParam::IonNames);
    }

    std::vector< LossIndex > forward_losses;
    std::vector< LossIndex > backward_losses;
    LossIndex losses_peptide2;

    if (!crosslink.alpha)
    {
      return;
    }
    AASequence alpha = *crosslink.alpha;
    AASequence beta;
    if (crosslink.beta) { beta = *crosslink.beta; }

    if (add_losses_)
    {
      if (frag_alpha)
      {
        losses_peptide2 = getBackwardLosses_(beta)[0];
        forward_losses = getForwardLosses_(alpha);
        backward_losses = getBackwardLosses_(alpha);
      }
      else
      {
        losses_peptide2 = getBackwardLosses_(alpha)[0];
        forward_losses = getForwardLosses_(beta);
        backward_losses = getBackwardLosses_(beta);
      }
    }

    for (Int z = mincharge; z <= maxcharge; ++z)
    {
      if (add_b_ions_)
      {
        addXLinkIonPeaks_(spectrum, charges, ion_names, crosslink, frag_alpha, Residue::BIon, forward_losses, backward_losses, losses_peptide2, z);
      }
      if (add_y_ions_)
      {
        addXLinkIonPeaks_(spectrum, charges, ion_names, crosslink, frag_alpha, Residue::YIon, forward_losses, backward_losses, losses_peptide2, z);
      }
      if (add_a_ions_)
      {
        addXLinkIonPeaks_(spectrum, charges, ion_names, crosslink, frag_alpha, Residue::AIon, forward_losses, backward_losses, losses_peptide2, z);
      }
      if (add_x_ions_)
      {
        addXLinkIonPeaks_(spectrum, charges, ion_names, crosslink, frag_alpha, Residue::XIon, forward_losses, backward_losses, losses_peptide2, z);
      }
      if (add_c_ions_)
      {
        addXLinkIonPeaks_(spectrum, charges, ion_names, crosslink, frag_alpha, Residue::CIon, forward_losses, backward_losses, losses_peptide2, z);
      }
      if (add_z_ions_)
      {
        addXLinkIonPeaks_(spectrum, charges, ion_names, crosslink, frag_alpha, Residue::ZIon, forward_losses, backward_losses, losses_peptide2, z);
      }
      if (add_k_linked_ions_ && !beta.empty())
      {
        double precursor_mass = alpha.getMonoWeight() + crosslink.cross_linker_mass;
        precursor_mass += beta.getMonoWeight();
        AASequence peptide;
        Size link_pos;
        if (frag_alpha)
        {
          peptide = alpha;
          link_pos = crosslink.cross_link_position.first;
        }
        else
        {
          peptide = beta;
          link_pos = crosslink.cross_link_position.second;
        }
        addKLinkedIonPeaks_(spectrum, charges, ion_names, peptide, link_pos, precursor_mass, frag_alpha, z);
      }
    }

    if (add_precursor_peaks_)
    {
      double precursor_mass = alpha.getMonoWeight() + crosslink.cross_linker_mass;
      if (!beta.empty())
      {
        precursor_mass += beta.getMonoWeight();
      }
      addPrecursorPeaks_(spectrum, charges, ion_names, precursor_mass, maxcharge);
    }

    if (add_charges_)
    {
      if (!spectrum.getIntegerDataArrays().empty())
      {
        spectrum.getIntegerDataArrays()[0] = charges;
      }
      else
      {
        spectrum.getIntegerDataArrays().push_back(charges);
      }
    }
    if (add_metainfo_)
    {
      if (!spectrum.getStringDataArrays().empty())
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

  void TheoreticalSpectrumGeneratorXLMS::addXLinkIonPeaks_(PeakSpectrum & spectrum, DataArrays::IntegerDataArray & charges, DataArrays::StringDataArray & ion_names, OPXLDataStructs::ProteinProteinCrossLink & crosslink, bool frag_alpha, Residue::ResidueType res_type, std::vector< LossIndex > & forward_losses, std::vector< LossIndex > & backward_losses, LossIndex & losses_peptide2, int charge) const
  {
    if (!crosslink.alpha || crosslink.alpha->empty())
    {
      cout << "Warning: Attempt at creating XLink Ions Spectrum from empty string!" << endl;
      return;
    }

    AASequence alpha = *crosslink.alpha;
    AASequence beta;
    if (crosslink.beta)  { beta = *crosslink.beta; }

    double precursor_mass = alpha.getMonoWeight() + crosslink.cross_linker_mass;

    if (!beta.empty())
    {
      precursor_mass += beta.getMonoWeight();
    }

    String ion_type;
    AASequence peptide;
    AASequence peptide2;
    Size link_pos;
    if (frag_alpha)
    {
      ion_type = "alpha|xi";
      peptide = alpha;
      peptide2 = beta;
      link_pos = crosslink.cross_link_position.first;
    }
    else
    {
      ion_type = "beta|xi";
      peptide = beta;
      peptide2 = alpha;
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
        int frag_index = i;

        addPeak_(spectrum, charges, ion_names, pos, intensity, res_type, frag_index, charge, ion_type);
        if (add_losses_ && forward_losses.size() >= i)
        {
          String ion_name = "[" + ion_type + "$" + String(Residue::residueTypeToIonLetter(res_type)) + String(frag_index) + "]";
          LossIndex losses = losses_peptide2;
          losses.has_H2O_loss = losses_peptide2.has_H2O_loss || forward_losses[i-1].has_H2O_loss;
          losses.has_NH3_loss = losses_peptide2.has_NH3_loss || forward_losses[i-1].has_NH3_loss;
          addXLinkIonLosses_(spectrum, charges, ion_names, mono_weight, intensity, charge, ion_name, losses);
        }
        if (add_isotopes_ && max_isotope_ >= 2) // add second isotopic peak with fast method, if two or more peaks are asked for
        {
          pos += Constants::C13C12_MASSDIFF_U / static_cast<double>(charge);
          addPeak_(spectrum, charges, ion_names, pos, intensity, res_type, frag_index, charge, ion_type);
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
        int frag_index = peptide.size() - 1 - i;

        addPeak_(spectrum, charges, ion_names, pos, intensity, res_type, frag_index, charge, ion_type);
        if (add_losses_ && backward_losses.size() >= i+2)
        {
          String ion_name = "[" + ion_type + "$" + String(Residue::residueTypeToIonLetter(res_type)) + String(frag_index) + "]";
          LossIndex losses = losses_peptide2;
          losses.has_H2O_loss = losses_peptide2.has_H2O_loss || backward_losses[i+1].has_H2O_loss;
          losses.has_NH3_loss = losses_peptide2.has_NH3_loss || backward_losses[i+1].has_NH3_loss;
          addXLinkIonLosses_(spectrum, charges, ion_names, mono_weight, intensity, charge, ion_name, losses);
        }

        if (add_isotopes_ && max_isotope_ >= 2) // add second isotopic peak with fast method, if two or more peaks are asked for
        {
          pos += Constants::C13C12_MASSDIFF_U / static_cast<double>(charge);
          addPeak_(spectrum, charges, ion_names, pos, intensity, res_type, frag_index, charge, ion_type);
        }
      }
    }
    return;
  }

  std::vector< TheoreticalSpectrumGeneratorXLMS::LossIndex > TheoreticalSpectrumGeneratorXLMS::getForwardLosses_(AASequence & peptide) const
  {
    // this gives us a "forward set" with incremental losses from the first to the last residue
    std::vector< LossIndex > ion_losses(peptide.size());
    ion_losses[0] = loss_db_.at(peptide[0].getOneLetterCode());
    for (Size i = 1; i < peptide.size(); ++i)
    {
      ion_losses[i].has_H2O_loss = ion_losses[i-1].has_H2O_loss || loss_db_.at(peptide[i].getOneLetterCode()).has_H2O_loss;
      ion_losses[i].has_NH3_loss = ion_losses[i-1].has_NH3_loss || loss_db_.at(peptide[i].getOneLetterCode()).has_NH3_loss;
    }
    return ion_losses;
  }

  std::vector< TheoreticalSpectrumGeneratorXLMS::LossIndex > TheoreticalSpectrumGeneratorXLMS::getBackwardLosses_(AASequence & peptide) const
  {
    // this gives us a "backward set" with incremental losses from the last to the first residue
    std::vector< LossIndex > ion_losses(peptide.size());
    ion_losses[ion_losses.size()-1] = loss_db_.at(peptide[peptide.size()-1].getOneLetterCode());
    for (Size i = ion_losses.size()-1; i > 0; --i)
    {
      ion_losses[i-1].has_H2O_loss = ion_losses[i].has_H2O_loss || loss_db_.at(peptide[i-1].getOneLetterCode()).has_H2O_loss;
      ion_losses[i-1].has_NH3_loss = ion_losses[i].has_NH3_loss || loss_db_.at(peptide[i-1].getOneLetterCode()).has_NH3_loss;
    }
    return ion_losses;
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
