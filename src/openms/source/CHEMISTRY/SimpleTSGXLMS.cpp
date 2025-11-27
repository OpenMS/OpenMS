// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Eugen Netz $
// $Authors: Eugen Netz $
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/SimpleTSGXLMS.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/KERNEL/MSSpectrum.h>

#if OPENMS_BOOST_VERSION_MINOR >= 67 && OPENMS_BOOST_VERSION_MAJOR == 1
#define OPENMS_USE_PDQSORT
#include <boost/sort/pdqsort/pdqsort.hpp>
#endif

using namespace std;

namespace OpenMS
{

  SimpleTSGXLMS::SimpleTSGXLMS() :
    DefaultParamHandler("SimpleTSGXLMS")
  {
    // TODO only partly functional (second isotopic peak if max_isotope = 2)
    defaults_.setValue("add_isotopes", "false", "If set to 1 isotope peaks of the product ion peaks are added");
    defaults_.setValidStrings("add_isotopes", {"true","false"});

    defaults_.setValue("max_isotope", 2, "Defines the maximal isotopic peak which is added, add_isotopes must be set to 1");

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

  SimpleTSGXLMS::~SimpleTSGXLMS() = default;

  void SimpleTSGXLMS::getLinearIonSpectrum(std::vector< SimplePeak >& spectrum, AASequence& peptide, Size link_pos, int charge, Size link_pos_2) const
  {
    std::vector< LossIndex > forward_losses;
    std::vector< LossIndex > backward_losses;

    if (add_losses_)
    {
      forward_losses = getForwardLosses_(peptide);
      backward_losses = getBackwardLosses_(peptide);
    }

    for (Int z = charge; z >= 1; --z)
    {
      if (add_b_ions_)
      {
        addLinearPeaks_(spectrum, peptide, link_pos, Residue::BIon, forward_losses, backward_losses, z, link_pos_2);
      }
      if (add_y_ions_)
      {
        addLinearPeaks_(spectrum, peptide, link_pos, Residue::YIon, forward_losses, backward_losses, z, link_pos_2);
      }
      if (add_a_ions_)
      {
        addLinearPeaks_(spectrum, peptide, link_pos, Residue::AIon, forward_losses, backward_losses, z, link_pos_2);
      }
      if (add_x_ions_)
      {
        addLinearPeaks_(spectrum, peptide, link_pos, Residue::XIon, forward_losses, backward_losses, z, link_pos_2);
      }
      if (add_c_ions_)
      {
        addLinearPeaks_(spectrum, peptide, link_pos, Residue::CIon, forward_losses, backward_losses, z, link_pos_2);
      }
      if (add_z_ions_)
      {
        addLinearPeaks_(spectrum, peptide, link_pos, Residue::ZIon, forward_losses, backward_losses, z, link_pos_2);
      }
    }

#ifdef OPENMS_USE_PDQSORT
    boost::sort::pdqsort_branchless(spectrum.begin(), spectrum.end(), [](const SimplePeak& a, const SimplePeak& b) {return a.mz < b.mz;});
#else
    std::stable_sort(spectrum.begin(), spectrum.end(), [](const SimplePeak& a, const SimplePeak& b) {return a.mz < b.mz;});
#endif

    return;
  }

  void SimpleTSGXLMS::addLinearPeaks_(std::vector< SimplePeak >& spectrum, AASequence& peptide, Size link_pos, Residue::ResidueType res_type, std::vector< LossIndex >& forward_losses, std::vector< LossIndex >& backward_losses, int charge, Size link_pos_2) const
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

    if (res_type == Residue::AIon || res_type == Residue::BIon || res_type == Residue::CIon)
    {
      double mono_weight(Constants::PROTON_MASS_U * charge);
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
        double pos(mono_weight / charge);

        if (add_losses_)
        {
          addLosses_(spectrum, mono_weight, charge, forward_losses[i]);
        }
        spectrum.emplace_back(pos, charge);

        if (add_isotopes_ && max_isotope_ >= 2) // add second isotopic peak with fast method, if two or more peaks are asked for
        {
          spectrum.emplace_back(pos+(Constants::C13C12_MASSDIFF_U / charge), charge);
        }
      }
    }
    else // if (res_type == Residue::XIon || res_type == Residue::YIon || res_type == Residue::ZIon)
    {
      double mono_weight(Constants::PROTON_MASS_U * charge);
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
        double pos(mono_weight / charge);

        if (add_losses_)
        {
          addLosses_(spectrum, pos, charge, backward_losses[i]);
        }
        spectrum.emplace_back(pos, charge);

        if (add_isotopes_ && max_isotope_ >= 2) // add second isotopic peak with fast method, if two or more peaks are asked for
        {
          spectrum.emplace_back(pos+(Constants::C13C12_MASSDIFF_U / charge), charge);
        }
      }
    }
    return;
  }

  void SimpleTSGXLMS::getXLinkIonSpectrum(std::vector< SimplePeak >& spectrum, AASequence& peptide, Size link_pos, double precursor_mass, int mincharge, int maxcharge, Size link_pos_2) const
  {
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
        addXLinkIonPeaks_(spectrum, peptide, link_pos, precursor_mass, Residue::BIon, forward_losses, backward_losses, z, link_pos_2);
      }
      if (add_y_ions_)
      {
        addXLinkIonPeaks_(spectrum, peptide, link_pos, precursor_mass, Residue::YIon, forward_losses, backward_losses, z, link_pos_2);
      }
      if (add_a_ions_)
      {
        addXLinkIonPeaks_(spectrum, peptide, link_pos, precursor_mass, Residue::AIon, forward_losses, backward_losses, z, link_pos_2);
      }
      if (add_x_ions_)
      {
        addXLinkIonPeaks_(spectrum, peptide, link_pos, precursor_mass, Residue::XIon, forward_losses, backward_losses, z, link_pos_2);
      }
      if (add_c_ions_)
      {
        addXLinkIonPeaks_(spectrum, peptide, link_pos, precursor_mass, Residue::CIon, forward_losses, backward_losses, z, link_pos_2);
      }
      if (add_z_ions_)
      {
        addXLinkIonPeaks_(spectrum, peptide, link_pos, precursor_mass, Residue::ZIon, forward_losses, backward_losses, z, link_pos_2);
      }
      if (add_k_linked_ions_)
      {
        addKLinkedIonPeaks_(spectrum, peptide, link_pos, precursor_mass, z);
      }
    }

    if (add_precursor_peaks_)
    {
      addPrecursorPeaks_(spectrum, precursor_mass, maxcharge);
    }

#ifdef OPENMS_USE_PDQSORT
    std::reverse(spectrum.begin(), spectrum.end());
    boost::sort::pdqsort_branchless(spectrum.begin(), spectrum.end(), [](const SimplePeak& a, const SimplePeak& b) {return a.mz < b.mz;});
#else
    std::sort(spectrum.begin(), spectrum.end(), [](const SimplePeak& a, const SimplePeak& b) {return a.mz < b.mz;});
#endif

    return;
  }

  void SimpleTSGXLMS::addXLinkIonPeaks_(std::vector< SimplePeak >& spectrum, AASequence& peptide, Size link_pos, double precursor_mass, Residue::ResidueType res_type, std::vector< LossIndex >& forward_losses, std::vector< LossIndex >& backward_losses, int charge, Size link_pos_2) const
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

    if (res_type == Residue::AIon || res_type == Residue::BIon || res_type == Residue::CIon)
    {
      // whole mass of both peptides + cross-link (or peptide + mono-link), converted to an internal ion
      double mono_weight((Constants::PROTON_MASS_U * charge) + precursor_mass - Residue::getInternalToFull().getMonoWeight());


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
        double pos(mono_weight / charge);

        if (add_isotopes_ && max_isotope_ >= 2) // add second isotopic peak with fast method, if two or more peaks are asked for
        {
          spectrum.emplace_back(pos+(Constants::C13C12_MASSDIFF_U / charge), charge);
        }
        spectrum.emplace_back(pos, charge);
        if (add_losses_ && forward_losses.size() >= i)
        {
          addLosses_(spectrum, mono_weight, charge, forward_losses[i-1]);
        }
      }
    }
    else // if (res_type == Residue::XIon || res_type == Residue::YIon || res_type == Residue::ZIon)
    {
      // whole mass of both peptides + cross-link (or peptide + mono-link), converted to an internal ion
      double mono_weight((Constants::PROTON_MASS_U * charge) + precursor_mass - Residue::getInternalToFull().getMonoWeight()); // whole mass

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
        double pos(mono_weight / charge);

        if (add_isotopes_ && max_isotope_ >= 2) // add second isotopic peak with fast method, if two or more peaks are asked for
        {
          spectrum.emplace_back(pos+(Constants::C13C12_MASSDIFF_U / charge), charge);
        }
        spectrum.emplace_back(pos, charge);
        if (add_losses_ && backward_losses.size() >= i+2)
        {
          addLosses_(spectrum, mono_weight, charge, backward_losses[i+1]);
        }
      }
    }
    return;
  }

  void SimpleTSGXLMS::addPrecursorPeaks_(std::vector< SimplePeak >& spectrum, double precursor_mass, int charge) const
  {
    // precursor peak
    double mono_pos = precursor_mass + (Constants::PROTON_MASS_U * charge);
    if (add_isotopes_ && max_isotope_ >= 2) // add second isotopic peak with fast method, if two or more peaks are asked for
    {
      spectrum.emplace_back((mono_pos + Constants::C13C12_MASSDIFF_U) / charge, charge);
    }
    spectrum.emplace_back(mono_pos / charge, charge);

    // loss peaks of the precursor
    // loss of water
    mono_pos = precursor_mass + (Constants::PROTON_MASS_U * charge) - loss_H2O_;
    if (add_isotopes_ && max_isotope_ >= 2) // add second isotopic peak with fast method, if two or more peaks are asked for
    {
      spectrum.emplace_back((mono_pos + Constants::C13C12_MASSDIFF_U) / charge, charge);
    }
    spectrum.emplace_back(mono_pos / charge, charge);

    //loss of ammonia
    mono_pos = precursor_mass + (Constants::PROTON_MASS_U * charge) - loss_NH3_;

    if (add_isotopes_ && max_isotope_ >= 2) // add second isotopic peak with fast method, if two or more peaks are asked for
    {
      spectrum.emplace_back((mono_pos + Constants::C13C12_MASSDIFF_U) / charge, charge);
    }
    spectrum.emplace_back(mono_pos / charge, charge);
  }

  void SimpleTSGXLMS::addKLinkedIonPeaks_(std::vector< SimplePeak >& spectrum, AASequence& peptide, Size link_pos, double precursor_mass, int charge) const
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

    mono_weight += Constants::PROTON_MASS_U * charge;
    if (mono_weight < 0)
    {
      return;
    }

    if (add_isotopes_ && max_isotope_ >= 2) // add second isotopic peak with fast method, if two or more peaks are asked for
    {
      spectrum.emplace_back((mono_weight + Constants::C13C12_MASSDIFF_U) / charge, charge);
    }
    spectrum.emplace_back(mono_weight / charge, charge);
  }

  void SimpleTSGXLMS::addLosses_(std::vector< SimplePeak >& spectrum, double mono_weight, int charge, LossIndex& losses) const
  {
    if (losses.has_H2O_loss)
    {
      spectrum.emplace_back((mono_weight - loss_H2O_) / charge, charge);
    }

    if (losses.has_NH3_loss)
    {
      spectrum.emplace_back((mono_weight - loss_NH3_) / charge, charge);
    }
  }

  void SimpleTSGXLMS::getXLinkIonSpectrum(std::vector< SimplePeak >& spectrum, OPXLDataStructs::ProteinProteinCrossLink& crosslink, bool frag_alpha, int mincharge, int maxcharge) const
  {
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
        addXLinkIonPeaks_(spectrum, crosslink, frag_alpha, Residue::BIon, forward_losses, backward_losses, losses_peptide2, z);
      }
      if (add_y_ions_)
      {
        addXLinkIonPeaks_(spectrum, crosslink, frag_alpha, Residue::YIon, forward_losses, backward_losses, losses_peptide2, z);
      }
      if (add_a_ions_)
      {
        addXLinkIonPeaks_(spectrum, crosslink, frag_alpha, Residue::AIon, forward_losses, backward_losses, losses_peptide2, z);
      }
      if (add_x_ions_)
      {
        addXLinkIonPeaks_(spectrum, crosslink, frag_alpha, Residue::XIon, forward_losses, backward_losses, losses_peptide2, z);
      }
      if (add_c_ions_)
      {
        addXLinkIonPeaks_(spectrum, crosslink, frag_alpha, Residue::CIon, forward_losses, backward_losses, losses_peptide2, z);
      }
      if (add_z_ions_)
      {
        addXLinkIonPeaks_(spectrum, crosslink, frag_alpha, Residue::ZIon, forward_losses, backward_losses, losses_peptide2, z);
      }
      if (add_k_linked_ions_ && !beta.empty())
      {
        double precursor_mass = alpha.getMonoWeight() + beta.getMonoWeight() + crosslink.cross_linker_mass;
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
        addKLinkedIonPeaks_(spectrum, peptide, link_pos, precursor_mass, z);
      }
    }

    if (add_precursor_peaks_)
    {
      double precursor_mass = alpha.getMonoWeight() + crosslink.cross_linker_mass;
      if (!beta.empty())
      {
        precursor_mass += beta.getMonoWeight();
      }
      addPrecursorPeaks_(spectrum, precursor_mass, maxcharge);
    }

#ifdef OPENMS_USE_PDQSORT
    std::reverse(spectrum.begin(), spectrum.end());
    boost::sort::pdqsort_branchless(spectrum.begin(), spectrum.end(), [](const SimplePeak& a, const SimplePeak& b) {return a.mz < b.mz;});
#else
    std::sort(spectrum.begin(), spectrum.end(), [](const SimplePeak& a, const SimplePeak& b) {return a.mz < b.mz;});
#endif

    return;
  }

  void SimpleTSGXLMS::addXLinkIonPeaks_(std::vector< SimplePeak >& spectrum, OPXLDataStructs::ProteinProteinCrossLink& crosslink, bool frag_alpha, Residue::ResidueType res_type, std::vector< LossIndex >& forward_losses, std::vector< LossIndex >& backward_losses, LossIndex& losses_peptide2, int charge) const
  {
    if (!crosslink.alpha || crosslink.alpha->empty())
    {
      cout << "Warning: Attempt at creating XLink Ions Spectrum from empty string!" << endl;
      return;
    }

    AASequence alpha = *crosslink.alpha;
    AASequence beta;
    if (crosslink.beta) { beta = *crosslink.beta; }

    double precursor_mass = alpha.getMonoWeight() + crosslink.cross_linker_mass;

    if (!beta.empty())
    {
      precursor_mass += beta.getMonoWeight();
    }

    AASequence peptide;
    AASequence peptide2;
    Size link_pos;
    if (frag_alpha)
    {
      peptide = alpha;
      peptide2 = beta;
      link_pos = crosslink.cross_link_position.first;
    }
    else
    {
      peptide = beta;
      peptide2 = alpha;
      link_pos = crosslink.cross_link_position.second;
    }

    if (res_type == Residue::AIon || res_type == Residue::BIon || res_type == Residue::CIon)
    {
      double mono_weight((Constants::PROTON_MASS_U * charge) + precursor_mass - Residue::getInternalToFull().getMonoWeight());

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

        double pos(mono_weight / charge);

        if (add_isotopes_ && max_isotope_ >= 2) // add second isotopic peak with fast method, if two or more peaks are asked for
        {
          spectrum.emplace_back(pos+(Constants::C13C12_MASSDIFF_U / charge), charge);
        }
        spectrum.emplace_back(pos, charge);
        if (add_losses_ && forward_losses.size() >= i)
        {
          SimpleTSGXLMS::LossIndex losses;
          losses.has_H2O_loss = losses_peptide2.has_H2O_loss || forward_losses[i-1].has_H2O_loss;
          losses.has_NH3_loss = losses_peptide2.has_NH3_loss || forward_losses[i-1].has_NH3_loss;
          addLosses_(spectrum, mono_weight, charge, losses);
        }
      }
    }
    else // if (res_type == Residue::XIon || res_type == Residue::YIon || res_type == Residue::ZIon)
    {
      // whole mass of both peptides + cross-link (or peptide + mono-link), converted to an internal ion
      double mono_weight((Constants::PROTON_MASS_U * charge) + precursor_mass - Residue::getInternalToFull().getMonoWeight()); // whole mass

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

        double pos(mono_weight / charge);

        if (add_isotopes_ && max_isotope_ >= 2) // add second isotopic peak with fast method, if two or more peaks are asked for
        {
          spectrum.emplace_back(pos+(Constants::C13C12_MASSDIFF_U / charge), charge);
        }
        spectrum.emplace_back(pos, charge);
        if (add_losses_ && backward_losses.size() >= i+2)
        {
          SimpleTSGXLMS::LossIndex losses;
          losses.has_H2O_loss = losses_peptide2.has_H2O_loss || backward_losses[i+1].has_H2O_loss;
          losses.has_NH3_loss = losses_peptide2.has_NH3_loss || backward_losses[i+1].has_NH3_loss;
          addLosses_(spectrum, mono_weight, charge, losses);
        }
      }
    }
    return;
  }

  std::vector< SimpleTSGXLMS::LossIndex > SimpleTSGXLMS::getForwardLosses_(AASequence& peptide) const
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

  std::vector< SimpleTSGXLMS::LossIndex > SimpleTSGXLMS::getBackwardLosses_(AASequence& peptide) const
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
    add_isotopes_ = param_.getValue("add_isotopes").toBool();
    add_precursor_peaks_ = param_.getValue("add_precursor_peaks").toBool();
    add_abundant_immonium_ions_ = param_.getValue("add_abundant_immonium_ions").toBool();
    max_isotope_ = static_cast<Int>(param_.getValue("max_isotope"));
    add_k_linked_ions_ = param_.getValue("add_k_linked_ions").toBool();
  }
}
