// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/FineIsotopePatternGenerator.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/KERNEL/MSSpectrum.h>

#include <unordered_set>

using namespace std;

namespace OpenMS
{

  TheoreticalSpectrumGenerator::TheoreticalSpectrumGenerator() :
    DefaultParamHandler("TheoreticalSpectrumGenerator")
  {
    defaults_.setValue("isotope_model", "none", "Model to use for isotopic peaks ('none' means no isotopic peaks are added, 'coarse' adds isotopic peaks in unit mass distance, 'fine' uses the hyperfine isotopic generator to add accurate isotopic peaks. Note that adding isotopic peaks is very slow.");
    defaults_.setValidStrings("isotope_model", ListUtils::create<String>("none,coarse,fine"));

    defaults_.setValue("max_isotope", 2, "Defines the maximal isotopic peak which is added if 'isotope_model' is 'coarse'");
    defaults_.setValue("max_isotope_probability", 0.05, "Defines the maximal isotopic probability to cover if 'isotope_model' is 'fine'");

    defaults_.setValue("add_metainfo", "false", "Adds the type of peaks as metainfo to the peaks, like y8+, [M-H2O+2H]++");
    defaults_.setValidStrings("add_metainfo", ListUtils::create<String>("true,false"));

    defaults_.setValue("add_losses", "false", "Adds common losses to those ion expect to have them, only water and ammonia loss is considered");
    defaults_.setValidStrings("add_losses", ListUtils::create<String>("true,false"));

    defaults_.setValue("sort_by_position", "true", "Sort output by position");
    defaults_.setValidStrings("sort_by_position", ListUtils::create<String>("true,false"));

    defaults_.setValue("add_precursor_peaks", "false", "Adds peaks of the unfragmented precursor ion to the spectrum");
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


  TheoreticalSpectrumGenerator::TheoreticalSpectrumGenerator(const TheoreticalSpectrumGenerator& rhs) :
    DefaultParamHandler(rhs)
  {
  }


  TheoreticalSpectrumGenerator& TheoreticalSpectrumGenerator::operator=(const TheoreticalSpectrumGenerator& rhs)
  {
    DefaultParamHandler::operator=(rhs);
    return *this;
  }


  TheoreticalSpectrumGenerator::~TheoreticalSpectrumGenerator()
  {
  }

  void TheoreticalSpectrumGenerator::getSpectrum(PeakSpectrum& spectrum, const AASequence& peptide, Int min_charge, Int max_charge) const
  {
    if (peptide.empty())
    {
      return;
    }

    MSSpectrum::Chunks chunks(spectrum);
    PeakSpectrum::StringDataArray* ion_names;
    PeakSpectrum::IntegerDataArray* charges;

    bool charges_dynamic = false, ion_names_dynamic = false;

    if (spectrum.getIntegerDataArrays().empty())
    {
      charges = new PeakSpectrum::IntegerDataArray();
      charges_dynamic = true;
    }
    else
    {
      charges = &(spectrum.getIntegerDataArrays()[0]);
    }
    if (spectrum.getStringDataArrays().empty())
    {
      ion_names = new PeakSpectrum::StringDataArray();
      ion_names_dynamic = true;
    }
    else
    {
      ion_names = &(spectrum.getStringDataArrays()[0]);
    }
    ion_names->setName("IonNames");
    charges->setName("Charges");

    for (Int z = min_charge; z <= max_charge; ++z)
    {
      if (add_b_ions_) addPeaks_(spectrum, peptide, *ion_names, *charges, chunks, Residue::BIon, z);
      if (add_y_ions_) addPeaks_(spectrum, peptide, *ion_names, *charges, chunks, Residue::YIon, z);
      if (add_a_ions_) addPeaks_(spectrum, peptide, *ion_names, *charges, chunks, Residue::AIon, z);
      if (add_c_ions_) addPeaks_(spectrum, peptide, *ion_names, *charges, chunks, Residue::CIon, z);
      if (add_x_ions_) addPeaks_(spectrum, peptide, *ion_names, *charges, chunks, Residue::XIon, z);
      if (add_z_ions_) addPeaks_(spectrum, peptide, *ion_names, *charges, chunks, Residue::ZIon, z);
    }

    if (add_precursor_peaks_)
    {
      if (add_all_precursor_charges_)
      {
        for (Int z = min_charge; z <= max_charge; ++z)
        {
          addPrecursorPeaks_(spectrum, peptide, *ion_names, *charges, z);
          chunks.add(false);
        }
      }
      else // add_all_precursor_charges_ = false, only add precursor with highest charge
      {
        addPrecursorPeaks_(spectrum, peptide, *ion_names, *charges, max_charge);
        chunks.add(false);
      }
    }

    if (add_abundant_immonium_ions_)
    {
      addAbundantImmoniumIons_(spectrum, peptide, *ion_names, *charges);
      chunks.add(true); // this chunk is ordered, as the if-statements in addAbundantImmoniumIons_() are in ascending order (by MZ)
    }

    if (add_metainfo_)
    {
      if (spectrum.getIntegerDataArrays().empty())
      {
        spectrum.getIntegerDataArrays().push_back(std::move(*charges));
      }
      if (spectrum.getStringDataArrays().empty())
      {
        spectrum.getStringDataArrays().push_back(std::move(*ion_names));
      }
    }

    if (charges_dynamic) delete charges;
    if (ion_names_dynamic) delete ion_names;

    if (sort_by_position_) spectrum.sortByPositionPresorted(chunks.getChunks());
  }


  void TheoreticalSpectrumGenerator::addAbundantImmoniumIons_(PeakSpectrum& spectrum, const AASequence& peptide, DataArrays::StringDataArray& ion_names, DataArrays::IntegerDataArray& charges) const
  {
    // Proline immonium ion (C4H8N)
    if (peptide.has(*ResidueDB::getInstance()->getResidue('P')))
    {
      if (add_metainfo_)
      {
        ion_names.emplace_back("iP");
        charges.push_back(1);
      }
      spectrum.emplace_back(70.0656, 1.0); // emplace_back(MZ, intensity)
    }

    // Cysteine (C2H6NS)
    if (peptide.has(*ResidueDB::getInstance()->getResidue('C')))
    {
      if (add_metainfo_)
      {
        ion_names.emplace_back("iC");
        charges.push_back(1);
      }
      spectrum.emplace_back(76.0221, 1.0);
    }

    // Iso/Leucin immonium ion (same mass for immonium ion)
    if (peptide.has(*ResidueDB::getInstance()->getResidue('L')))
    {
      if (add_metainfo_)
      {
        ion_names.emplace_back("iL/I");
        charges.push_back(1);
      }
      spectrum.emplace_back(86.09698, 1.0);
    }

    // Histidin immonium ion (C5H8N3)
    if (peptide.has(*ResidueDB::getInstance()->getResidue('H')))
    {
      if (add_metainfo_)
      {
        ion_names.emplace_back("iH");
        charges.push_back(1);
      }
      spectrum.emplace_back(110.0718, 1.0);
    }

    // Phenylalanin immonium ion (C8H10N)
    if (peptide.has(*ResidueDB::getInstance()->getResidue('F')))
    {
      if (add_metainfo_)
      {
        ion_names.emplace_back("iF");
        charges.push_back(1);
      }
      spectrum.emplace_back(120.0813, 1.0);
    }

    // Tyrosine immonium ion (C8H10NO)
    if (peptide.has(*ResidueDB::getInstance()->getResidue('Y')))
    {
      if (add_metainfo_)
      {
        ion_names.emplace_back("iY");
        charges.push_back(1);
      }
      spectrum.emplace_back(136.0762, 1.0);
    }

    // Tryptophan immonium ion
    if (peptide.has(*ResidueDB::getInstance()->getResidue('W')))
    {
      if (add_metainfo_)
      {
        ion_names.emplace_back("iW");
        charges.push_back(1);
      }
      spectrum.emplace_back(159.0922, 1.0);
    }
  }


  char TheoreticalSpectrumGenerator::residueTypeToIonLetter_(const Residue::ResidueType res_type)
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
       OPENMS_LOG_ERROR << "Unknown residue type encountered. Can't map to ion letter." << endl;
    }
    return ' ';
  }


  void TheoreticalSpectrumGenerator::addIsotopeCluster_(PeakSpectrum& spectrum,
                                                        const AASequence& ion,
                                                        DataArrays::StringDataArray& ion_names,
                                                        DataArrays::IntegerDataArray& charges,
                                                        const Residue::ResidueType res_type,
                                                        Int charge,
                                                        double intensity) const
  {
    const String ion_name = String(Residue::residueTypeToIonLetter(res_type)) + String(ion.size()) + String((Size)abs(charge), '+');

    // manually compute correct sum formula (instead of using built-in assumption of hydrogen adduct)
    EmpiricalFormula f = ion.getFormula(res_type, charge) + EmpiricalFormula("H") * charge;
    f.setCharge(0);

    IsotopeDistribution dist;
    if (isotope_model_ == 1)
    {
      dist = f.getIsotopeDistribution(CoarseIsotopePatternGenerator(max_isotope_));
    }
    else if (isotope_model_ == 2)
    {
      dist = f.getIsotopeDistribution(FineIsotopePatternGenerator(max_isotope_probability_));
    }

    for (const auto& it : dist)
    {
      if (add_metainfo_) // one entry per peak
      {
        ion_names.push_back(ion_name);
        charges.push_back(charge);
      }
      spectrum.emplace_back(it.getMZ() / charge, intensity * it.getIntensity());
    }
  }

  void TheoreticalSpectrumGenerator::addLossesFaster_(PeakSpectrum& spectrum,
                         double mz,
                         const std::set<EmpiricalFormula>& f_losses,
                         int ion_ordinal,
                         DataArrays::StringDataArray& ion_names,
                         DataArrays::IntegerDataArray& charges,
                         const std::map<EmpiricalFormula, String>& formula_str_cache,
                         double intensity,
                         const Residue::ResidueType res_type,
                         bool add_metainfo,
                         int charge) const
  {
    const String charge_str((Size)abs(charge), '+');
    const String residue_str(Residue::residueTypeToIonLetter(res_type));
    const String ion_ordinal_str(String(ion_ordinal) + "-");

    for (auto& formula : f_losses)
    {
      spectrum.emplace_back((mz - formula.getMonoWeight()) / (double)charge, intensity);

      if (add_metainfo)
      {
        const String& loss_name = formula_str_cache.at(formula);
        // note: important to construct a string from char. If omitted it will perform pointer arithmetics on the "-" string literal
        ion_names.emplace_back(residue_str);
        //note: size of Residue::residueTypeToIonLetter(res_type) : 1;
        ion_names.back().reserve(1 + ion_ordinal_str.size() + loss_name.size() + charge_str.size());
        ((ion_names.back() += ion_ordinal_str) += loss_name) += charge_str;
        charges.push_back(charge);
      }
    }
  }

  void TheoreticalSpectrumGenerator::addLosses_(PeakSpectrum& spectrum,
                                                const AASequence& ion,
                                                DataArrays::StringDataArray& ion_names,
                                                DataArrays::IntegerDataArray& charges,
                                                double intensity,
                                                const Residue::ResidueType res_type,
                                                int charge) const
  {
    const String charge_str((Size)abs(charge), '+');
    const String residue_str(Residue::residueTypeToIonLetter(res_type));
    const String ion_str(String(ion.size()) + "-");

    std::set<String> losses;
    for (const auto& it : ion)
    {
      if (it.hasNeutralLoss())
      {
        for (const auto& formula : it.getLossFormulas())
        {
          losses.insert(formula.toString());
        }
      }
    }

    spectrum.reserve(spectrum.size() + losses.size());
    String ion_name;
    for (const auto& it : losses)
    {
      EmpiricalFormula loss_ion = ion.getFormula(res_type, charge) - EmpiricalFormula(it);
      // see 74e2ce6761e4a273164b29b8be487
      // thanks to Chris and Sandro
      // check for negative element frequencies (might happen if losses are not allowed for specific ions)
      bool negative_elements(false);
      for (const auto& eit : loss_ion)
      {
        if (eit.second < 0)
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
      const String& loss_name = it;

      ion_name = residue_str + ion_str + loss_name + charge_str;

      if (add_isotopes_)
      {
        // manually compute correct sum formula (instead of using built-in assumption of hydrogen adduct)
        loss_ion += EmpiricalFormula("H") * charge;
        loss_ion.setCharge(0);

        IsotopeDistribution dist;
        if (isotope_model_ == 1)
        {
          dist = loss_ion.getIsotopeDistribution(CoarseIsotopePatternGenerator(max_isotope_));
        }
        else if (isotope_model_ == 2)
        {
          dist = loss_ion.getIsotopeDistribution(FineIsotopePatternGenerator(max_isotope_probability_));
        }

        for (const auto& iso : dist)
        {
          if (add_metainfo_)
          {
            ion_names.push_back(ion_name);
            charges.push_back(charge);
          }
          spectrum.emplace_back(iso.getMZ() / (double)charge, intensity * rel_loss_intensity_ * iso.getIntensity());
        }
      }
      else
      {
        if (add_metainfo_)
        {
          ion_names.push_back(ion_name);
          charges.push_back(charge);
        }
        spectrum.emplace_back(loss_pos / (double)charge, intensity * rel_loss_intensity_);
      }
    }
  }


  void TheoreticalSpectrumGenerator::addPeaks_(PeakSpectrum& spectrum,
                                               const AASequence& peptide,
                                               DataArrays::StringDataArray& ion_names,
                                               DataArrays::IntegerDataArray& charges,
                                               MSSpectrum::Chunks& chunks,
                                               const Residue::ResidueType res_type,
                                               Int charge) const
  {
    const String charge_str((Size)abs(charge), '+');
    const String residue_str(Residue::residueTypeToIonLetter(res_type));

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

    std::set<EmpiricalFormula> fx_losses;
    std::map<EmpiricalFormula, String> formula_str_cache;

    // precompute formula_str_cache
    if (add_losses_)
    {
      for (auto& p : peptide)
      {
        for (auto& formula : p.getLossFormulas())
        {
          String& loss_name = formula_str_cache[formula];
          if (loss_name.empty())
          {
            loss_name = formula.toString();
          }
        }
      }
    }

    if (res_type == Residue::AIon || res_type == Residue::BIon || res_type == Residue::CIon)
    {
      if (peptide.hasNTerminalModification())
      {
        mono_weight += peptide.getNTerminalModification()->getDiffMonoMass();
      }
      double initial_mono_weight(mono_weight);

      static double stat_a = Residue::getInternalToAIon().getMonoWeight();
      static double stat_b = Residue::getInternalToBIon().getMonoWeight();
      static double stat_c = Residue::getInternalToCIon().getMonoWeight();

      if (!add_isotopes_) // add single peak
      {
        Size i = Size(!add_first_prefix_ion_);
        if (i == 1)
        {
          mono_weight += peptide[0].getMonoWeight(Residue::Internal);
          if (peptide[0].hasNeutralLoss())
          {
            for (const auto& formula : peptide[0].getLossFormulas()) fx_losses.insert(formula);
          }
        }
        for (; i < peptide.size() - 1; ++i)
        {
          mono_weight += peptide[i].getMonoWeight(Residue::Internal); // standard internal residue including named modifications: c
          double pos(mono_weight);

          double ion_offset = 0;
          switch (res_type)
          {
            case Residue::AIon: ion_offset = stat_a; break;
            case Residue::BIon: ion_offset = stat_b; break;
            case Residue::CIon: ion_offset = stat_c; break;
            default: break;
          }
          pos = (pos + ion_offset) / charge;

          spectrum.emplace_back(pos, intensity);
          if (add_metainfo_)
          {
            ion_names.emplace_back(residue_str);
            //note: size of Residue::residueTypeToIonLetter(res_type) : 1. size of String(i + 1) : 2;
            ion_names.back().reserve(1 + 2 + charge_str.size());
            (ion_names.back() += (i + 1)) += charge_str;
            charges.push_back(charge);
          }
        }
        chunks.add(true);

        mono_weight = initial_mono_weight;
        if (add_losses_)
        {
          for (i = Size(!add_first_prefix_ion_); i < peptide.size() - 1; ++i)
          {
            mono_weight += peptide[i].getMonoWeight(Residue::Internal); // standard internal residue including named modifications: c

            double ion_offset = 0;
            switch (res_type)
            {
              case Residue::AIon: ion_offset = stat_a; break;
              case Residue::BIon: ion_offset = stat_b; break;
              case Residue::CIon: ion_offset = stat_c; break;
              default: break;
            }
            if (peptide[i].hasNeutralLoss())
            {
              for (const auto& formula : peptide[i].getLossFormulas()) fx_losses.insert(formula);
            }
            addLossesFaster_(spectrum, mono_weight + ion_offset, fx_losses,
                              i + 1, ion_names, charges, formula_str_cache, intensity * rel_loss_intensity_,
                              res_type, add_metainfo_, charge);
            chunks.add(false); // unfortunately, the losses are not always inserted in sorted order
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
        chunks.add(true);

        if (add_losses_)
        {
          // add loss peaks (slow)
          i = add_first_prefix_ion_ ? 1 : 2;
          for (; i < peptide.size(); ++i)
          {
            const AASequence ion = peptide.getPrefix(i);
            addLosses_(spectrum, ion, ion_names, charges, intensity, res_type, charge);
          }
          chunks.add(true);
        }
      }
    }
    else // if (res_type == Residue::XIon || res_type == Residue::YIon || res_type == Residue::ZIon)
    {
      if (peptide.hasCTerminalModification())
      {
        mono_weight += peptide.getCTerminalModification()->getDiffMonoMass();
      }
      double initial_mono_weight(mono_weight);

      static double stat_x = Residue::getInternalToXIon().getMonoWeight();
      static double stat_y = Residue::getInternalToYIon().getMonoWeight();
      static double stat_z = Residue::getInternalToZIon().getMonoWeight();

      if (!add_isotopes_) // add single peak
      {
        for (Size i = peptide.size() - 1; i > 0; --i)
        {
          mono_weight += peptide[i].getMonoWeight(Residue::Internal); // standard internal residue including named modifications: c

          double pos(mono_weight);
          double ion_offset = 0;
          switch (res_type)
          {
            case Residue::XIon: ion_offset = stat_x; break;
            case Residue::YIon: ion_offset = stat_y; break;
            case Residue::ZIon: ion_offset = stat_z; break;
            default: break;
          }
          pos = (pos + ion_offset) / charge;

          spectrum.emplace_back(pos, intensity);
          if (add_metainfo_)
          {
            ion_names.emplace_back(residue_str);
            //note: size of Residue::residueTypeToIonLetter(res_type) => 1, size of String(peptide.size() - i) => 3;
            ion_names.back().reserve(1 + 3 + charge_str.size());
            (ion_names.back() += Size(peptide.size() - i)) += charge_str;
            charges.push_back(charge);
          }
        }
        chunks.add(true);

        if (add_losses_)
        {
          mono_weight = initial_mono_weight;
          for (Size i = peptide.size() - 1; i > 0; --i)
          {
            mono_weight += peptide[i].getMonoWeight(Residue::Internal); // standard internal residue including named modifications: c
            double ion_offset = 0;
            switch (res_type)
            {
              case Residue::XIon: ion_offset = stat_x; break;
              case Residue::YIon: ion_offset = stat_y; break;
              case Residue::ZIon: ion_offset = stat_z; break;
              default: break;
            }

            if (peptide[i].hasNeutralLoss())
            {
              for (const auto& formula : peptide[i].getLossFormulas()) fx_losses.insert(formula);
            }
            addLossesFaster_(spectrum, mono_weight + ion_offset, fx_losses,
                              peptide.size() - i, ion_names, charges, formula_str_cache, intensity * rel_loss_intensity_,
                              res_type, add_metainfo_, charge);
            chunks.add(false); // losses are not always added in sorted order
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
        chunks.add(true);

        if (add_losses_)
        {
          // add loss peaks (slow)
          for (Size i = 1; i < peptide.size(); ++i)
          {
            const AASequence ion = peptide.getSuffix(i);
            addLosses_(spectrum, ion, ion_names, charges, intensity, res_type, charge);
          }
          chunks.add(true);
        }
      }
    }
  }


  void TheoreticalSpectrumGenerator::addPrecursorPeaks_(PeakSpectrum& spectrum,
                                                        const AASequence& peptide,
                                                        DataArrays::StringDataArray& ion_names,
                                                        DataArrays::IntegerDataArray& charges,
                                                        Int charge) const
  {
    const String charge_str((Size)abs(charge), '+');
    const String ion_name("[M+H]" + charge_str);

    // precursor peak
    double mono_pos = peptide.getMonoWeight(Residue::Full, charge);

    if (add_isotopes_)
    {
      // manually compute correct sum formula (instead of using built-in assumption of hydrogen adduct)
      auto formula = peptide.getFormula(Residue::Full, charge) + EmpiricalFormula("H") * charge;
      formula.setCharge(0);

      IsotopeDistribution dist;
      if (isotope_model_ == 1)
      {
        dist = formula.getIsotopeDistribution(CoarseIsotopePatternGenerator(max_isotope_));
      }
      else if (isotope_model_ == 2)
      {
        dist = formula.getIsotopeDistribution(FineIsotopePatternGenerator(max_isotope_probability_));
      }

      for (IsotopeDistribution::ConstIterator it = dist.begin(); it != dist.end(); ++it)
      {
        if (add_metainfo_)
        {
          ion_names.push_back(ion_name);
          charges.push_back(charge);
        }
        spectrum.emplace_back(it->getMZ() / (double)charge, pre_int_ * it->getIntensity());
      }
    }
    else
    {
      if (add_metainfo_)
      {
        ion_names.push_back(ion_name);
        charges.push_back(charge);
      }
      spectrum.emplace_back(mono_pos / (double)charge, pre_int_);
    }
    // loss peaks of the precursor

    //loss of water
    EmpiricalFormula ion = peptide.getFormula(Residue::Full, charge) - EmpiricalFormula("H2O");
    mono_pos = ion.getMonoWeight();
    const String ion_name_h2o("[M+H]-H2O" + charge_str);
    if (add_isotopes_)
    {
      ion += EmpiricalFormula("H") * charge;
      ion.setCharge(0);

      IsotopeDistribution dist;
      if (isotope_model_ == 1)
      {
        dist = ion.getIsotopeDistribution(CoarseIsotopePatternGenerator(max_isotope_));
      }
      else if (isotope_model_ == 2)
      {
        dist = ion.getIsotopeDistribution(FineIsotopePatternGenerator(max_isotope_probability_));
      }

      for (IsotopeDistribution::ConstIterator it = dist.begin(); it != dist.end(); ++it)
      {
        if (add_metainfo_)
        {
          ion_names.push_back(ion_name_h2o);
          charges.push_back(charge);
        }
        spectrum.emplace_back(it->getMZ() / charge, pre_int_H2O_ * it->getIntensity());
      }
    }
    else
    {
      if (add_metainfo_)
      {
        ion_names.push_back(ion_name_h2o);
        charges.push_back(charge);
      }
      spectrum.emplace_back(mono_pos / (double)charge, pre_int_H2O_);
    }

    //loss of ammonia
    ion = peptide.getFormula(Residue::Full, charge) - EmpiricalFormula("NH3");
    mono_pos = ion.getMonoWeight();
    const String ion_name_nh3("[M+H]-NH3" + charge_str);
    if (add_isotopes_)
    {
      // manually compute correct sum formula (instead of using built-in assumption of hydrogen adduct)
      ion += EmpiricalFormula("H") * charge;
      ion.setCharge(0);

      IsotopeDistribution dist; 
      if (isotope_model_ == 1)
      {
        dist = ion.getIsotopeDistribution(CoarseIsotopePatternGenerator(max_isotope_));
      }
      else if (isotope_model_ == 2)
      {
        dist = ion.getIsotopeDistribution(FineIsotopePatternGenerator(max_isotope_probability_));
      }

      for (IsotopeDistribution::ConstIterator it = dist.begin(); it != dist.end(); ++it)
      {
        if (add_metainfo_)
        {
          ion_names.push_back(ion_name_nh3);
          charges.push_back(charge);
        }
        spectrum.emplace_back(it->getMZ() / (double)charge, pre_int_NH3_ * it->getIntensity());
      }
    }
    else
    {
      if (add_metainfo_)
      {
        ion_names.push_back(ion_name_nh3);
        charges.push_back(charge);
      }
      spectrum.emplace_back(mono_pos / (double)charge, pre_int_NH3_);
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
    add_first_prefix_ion_ = param_.getValue("add_first_prefix_ion").toBool();
    add_losses_ = param_.getValue("add_losses").toBool();
    add_metainfo_ = param_.getValue("add_metainfo").toBool();
    add_isotopes_ = param_.getValue("isotope_model") != "none";
    if (param_.getValue("isotope_model") == "coarse") isotope_model_ = 1;
    else if (param_.getValue("isotope_model") == "fine") isotope_model_ = 2;
    sort_by_position_ = param_.getValue("sort_by_position").toBool();
    add_precursor_peaks_ = param_.getValue("add_precursor_peaks").toBool();
    add_all_precursor_charges_ = param_.getValue("add_all_precursor_charges").toBool();
    add_abundant_immonium_ions_ = param_.getValue("add_abundant_immonium_ions").toBool();
    a_intensity_ = (double)param_.getValue("a_intensity");
    b_intensity_ = (double)param_.getValue("b_intensity");
    c_intensity_ = (double)param_.getValue("c_intensity");
    x_intensity_ = (double)param_.getValue("x_intensity");
    y_intensity_ = (double)param_.getValue("y_intensity");
    z_intensity_ = (double)param_.getValue("z_intensity");
    max_isotope_ = (Int)param_.getValue("max_isotope");
    max_isotope_probability_ = param_.getValue("max_isotope_probability");
    rel_loss_intensity_ = (double)param_.getValue("relative_loss_intensity");
    pre_int_ = (double)param_.getValue("precursor_intensity");
    pre_int_H2O_ = (double)param_.getValue("precursor_H2O_intensity");
    pre_int_NH3_ = (double)param_.getValue("precursor_NH3_intensity");
  }

} // end namespace OpenMS
