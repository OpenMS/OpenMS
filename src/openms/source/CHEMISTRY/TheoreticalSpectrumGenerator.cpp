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
#include <OpenMS/CONCEPT/LogStream.h>
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


  void TheoreticalSpectrumGenerator::addFragmentPeaks_(
    PeakSpectrum& spectrum, const vector<EmpiricalFormula>& fragment_forms,
    const String& ion_type, const EmpiricalFormula& offset, double intensity,
    Size start) const
  {
    for (Size i = start; i < fragment_forms.size(); ++i)
    {
      EmpiricalFormula fragment = fragment_forms[i] + offset;
      Peak1D peak(fragment.getMonoWeight(), intensity);
      spectrum.push_back(peak);
    }
    if (add_metainfo_)
    {
      for (Size i = start; i < fragment_forms.size(); ++i)
      {
        String ion_name = ion_type + String(i + 1);
        spectrum.getStringDataArrays()[0].push_back(ion_name);
      }
    }
  }


  void TheoreticalSpectrumGenerator::addAMinusBPeaks_(
    PeakSpectrum& spectrum, const vector<EmpiricalFormula>& fragment_forms,
    const NASequence& oligo, Size start) const
  {
    // offset: phosphate (from bond) minus 3 water (from various reactions)
    static const EmpiricalFormula offset = EmpiricalFormula("H-5P");
    // offset for first ("a1-B") ion: loss of 2 water
    static const EmpiricalFormula initial_offset = EmpiricalFormula("H-4O-2");
    // methyl group may be retained on ribose for "ambiguous" mods:
    static const EmpiricalFormula methyl_form = EmpiricalFormula("CH2");

    for (Size i = start; i < fragment_forms.size(); ++i)
    {
      EmpiricalFormula fragment = oligo[i]->getBaselossFormula();
      if (i > 0)
      {
        // base at position "i" is lost, so use fragment up to pos. "i - 1":
        fragment += fragment_forms[i - 1] + offset;
      }
      else // first ribonucleotide
      {
        fragment += initial_offset;
      }
      Peak1D peak(fragment.getMonoWeight(), aB_intensity_);
      if (oligo[i]->isAmbiguous())
      {
        // special treatment for a-B ions of "ambiguous" modifications:
        // create two peaks with half intensity, representing methyl group
        // lost/retained on backbone:
        peak.setIntensity(aB_intensity_ * 0.5);
        spectrum.push_back(peak);
        fragment += methyl_form;
        peak.setMZ(fragment.getMonoWeight());
      }
      spectrum.push_back(peak);
    }
    if (add_metainfo_)
    {
      for (Size i = start; i < fragment_forms.size(); ++i)
      {
        String ion_name = "a" + String(i + 1) + "-B";
        spectrum.getStringDataArrays()[0].push_back(ion_name);
        if (oligo[i]->isAmbiguous()) // two peaks were added
        {
          spectrum.getStringDataArrays()[0].push_back(ion_name);
        }
      }
    }
  }


  PeakSpectrum TheoreticalSpectrumGenerator::getUnchargedSpectrum_(
    const NASequence& oligo) const
  {
    // lots of code copied from "NASequence::getFormula" - can we avoid this?
    // @TODO: perform calculations with formulas or with masses?
    static const EmpiricalFormula H_form = EmpiricalFormula("H");
    // phosphate minus water:
    static const EmpiricalFormula backbone_form = EmpiricalFormula("H-1PO2");
    static const EmpiricalFormula a_ion_offset = EmpiricalFormula("H-2O-1");
    static const EmpiricalFormula b_ion_offset = EmpiricalFormula("");
    static const EmpiricalFormula c_ion_offset = backbone_form;
    static const EmpiricalFormula d_ion_offset = EmpiricalFormula("HPO3");
    static const EmpiricalFormula w_ion_offset = d_ion_offset;
    static const EmpiricalFormula x_ion_offset = c_ion_offset;
    static const EmpiricalFormula y_ion_offset = b_ion_offset;
    static const EmpiricalFormula z_ion_offset = a_ion_offset;

    PeakSpectrum spectrum;
    if (oligo.empty()) return spectrum;

    EmpiricalFormula three_prime_form, five_prime_form;
    if (oligo.getThreePrimeMod() != nullptr)
    {
      three_prime_form = oligo.getThreePrimeMod()->getFormula() - H_form;
    }
    if (oligo.getFivePrimeMod() != nullptr)
    {
      five_prime_form = oligo.getFivePrimeMod()->getFormula() - H_form;
    }

    vector<EmpiricalFormula> ribo_forms(oligo.size());
    Size index = 0;
    for (auto ribo : oligo)
    {
      ribo_forms[index] = ribo.getFormula();
      ++index;
    }

    spectrum.getStringDataArrays().resize(1);
    spectrum.getStringDataArrays()[0].setName("IonNames");

    vector<EmpiricalFormula> fragments_left, fragments_right;
    Size start = add_first_prefix_ion_ ? 0 : 1;
    if ((add_a_ions_ || add_b_ions_ || add_c_ions_ || add_d_ions_ ||
         add_aB_ions_) && (oligo.size() > start + 1))
    {
      fragments_left.resize(oligo.size() - 1);
      fragments_left[0] = ribo_forms[0] + five_prime_form;
      for (Size i = 1; i < oligo.size() - 1; ++i)
      {
        fragments_left[i] = (fragments_left[i - 1] + ribo_forms[i] +
                             backbone_form);
      }
      if (add_a_ions_)
      {
        addFragmentPeaks_(spectrum, fragments_left, "a", a_ion_offset,
                          a_intensity_, start);
      }
      if (add_b_ions_)
      {
        addFragmentPeaks_(spectrum, fragments_left, "b", b_ion_offset,
                          b_intensity_, start);
      }
      if (add_c_ions_)
      {
        addFragmentPeaks_(spectrum, fragments_left, "c", c_ion_offset,
                          c_intensity_, start);
      }
      if (add_d_ions_)
      {
        addFragmentPeaks_(spectrum, fragments_left, "d", d_ion_offset,
                          d_intensity_, start);
      }
      if (add_aB_ions_) // special case
      {
        addAMinusBPeaks_(spectrum, fragments_left, oligo, start);
      }
    }

    if ((add_w_ions_ || add_x_ions_ || add_y_ions_ || add_z_ions_) &&
        (oligo.size() > 1))
    {
      fragments_right.resize(oligo.size() - 1);
      fragments_right[0] = ribo_forms.back() + three_prime_form;
      for (Size i = 1; i < oligo.size() - 1; ++i)
      {
        Size ribo_index = oligo.size() - i - 1;
        fragments_right[i] = (fragments_right[i - 1] + ribo_forms[ribo_index] +
                              backbone_form);
      }
      if (add_w_ions_)
      {
        addFragmentPeaks_(spectrum, fragments_right, "w", w_ion_offset,
                          w_intensity_);
      }
      if (add_x_ions_)
      {
        addFragmentPeaks_(spectrum, fragments_right, "x", x_ion_offset,
                          x_intensity_);
      }
      if (add_y_ions_)
      {
        addFragmentPeaks_(spectrum, fragments_right, "y", y_ion_offset,
                          y_intensity_);
      }
      if (add_z_ions_)
      {
        addFragmentPeaks_(spectrum, fragments_right, "z", z_ion_offset,
                          z_intensity_);
      }
    }

    if (add_precursor_peaks_) // re-use what we've already calculated
    {
      Peak1D peak(0.0, pre_int_);
      bool have_left = !fragments_left.empty();
      bool have_right = !fragments_right.empty();
      if (have_left && have_right)
      {
        fragments_left[0] += fragments_right.back() + backbone_form;
        peak.setMZ(fragments_left[0].getMonoWeight());
      }
      else if (have_left)
      {
        fragments_left.back() += ribo_forms.back() + backbone_form +
          three_prime_form;
        peak.setMZ(fragments_left.back().getMonoWeight());
      }
      else if (have_right)
      {
        fragments_right.back() += ribo_forms[0] + backbone_form +
          five_prime_form;
        peak.setMZ(fragments_right.back().getMonoWeight());
      }
      else // really, no fragment ions?
      {
        peak.setMZ(oligo.getMonoWeight(NASequence::Full, 0));
      }
      spectrum.push_back(peak);
      if (add_metainfo_)
      {
        spectrum.getStringDataArrays()[0].push_back("M");
      }
    }

    return spectrum;
  }


  void TheoreticalSpectrumGenerator::addChargedSpectrum_(
    PeakSpectrum& spectrum, const PeakSpectrum& uncharged_spectrum, Int charge,
    bool add_precursor) const
  {
    // @TODO: use "Constants::PROTON_MASS_U" here instead?
    static const double H_mass = EmpiricalFormula("H").getMonoWeight();
    if (uncharged_spectrum.empty()) return;
    Size size = uncharged_spectrum.size();
    if (add_precursor_peaks_ && !add_precursor)
    {
      --size; // uncharged spectrum contains precursor peak - exclude it
    }
    for (Size i = 0; i < size; ++i)
    {
      spectrum.push_back(uncharged_spectrum[i]);
      double mass = spectrum.back().getMZ() + charge * H_mass;
      spectrum.back().setMZ(abs(mass / charge));
    }
    if (add_metainfo_)
    {
      auto& ions = spectrum.getStringDataArrays()[0];
      auto source_it = uncharged_spectrum.getStringDataArrays()[0].begin();
      ions.insert(ions.end(), source_it, source_it + size);
      auto& charges = spectrum.getIntegerDataArrays()[0];
      charges.resize(charges.size() + size, charge);
    }
  }


  void TheoreticalSpectrumGenerator::addSimpleSpectrum_(
    PeakSpectrum& spectrum, const NASequence& oligo, Int charge,
    bool add_precursor, bool sort) const
  {
    if (add_metainfo_)
    {
      if (spectrum.getIntegerDataArrays().empty())
      {
        spectrum.getIntegerDataArrays().resize(1);
        spectrum.getIntegerDataArrays()[0].setName("Charges");
      }
      if (spectrum.getStringDataArrays().empty())
      {
        spectrum.getStringDataArrays().resize(1);
        spectrum.getStringDataArrays()[0].setName("IonNames");
      }
    }

    if (add_b_ions_) addPeaks_(spectrum, oligo, NASequence::BIon, charge);
    if (add_y_ions_) addPeaks_(spectrum, oligo, NASequence::YIon, charge);
    if (add_a_ions_) addPeaks_(spectrum, oligo, NASequence::AIon, charge);
    if (add_c_ions_) addPeaks_(spectrum, oligo, NASequence::CIon, charge);
    if (add_x_ions_) addPeaks_(spectrum, oligo, NASequence::XIon, charge);
    if (add_z_ions_) addPeaks_(spectrum, oligo, NASequence::ZIon, charge);
    if (add_d_ions_) addPeaks_(spectrum, oligo, NASequence::DIon, charge);
    if (add_w_ions_) addPeaks_(spectrum, oligo, NASequence::WIon, charge);
    if (add_aB_ions_) addPeaks_(spectrum, oligo, NASequence::AminusB, charge);

    if (add_precursor) addPrecursorPeaks_(spectrum, oligo, charge);

    if (sort) spectrum.sortByPosition();
  }


  void TheoreticalSpectrumGenerator::getSpectrum(PeakSpectrum& spectrum, const NASequence& oligo, Int min_charge, Int max_charge) const
  {
    Int sign = 1;
    if (max_charge < 0 && min_charge < 0) // negative mode
    {
      sign = -1;
    }
    else if (max_charge * min_charge < 0)
    {
      // Signs don't match - we need to quit and thow error here to avoid messing up for loops below
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "min. and max. charge must both be either positive or negative");
    }
    if (abs(max_charge) < abs(min_charge))
    {
      swap(max_charge, min_charge);
    }

    for (uint z = (uint)abs(min_charge); z <= (uint)abs(max_charge) && z < (uint)oligo.size(); ++z)
    {
      bool add_precursor =
        ((add_precursor_peaks_ && add_all_precursor_charges_) ||
         (add_precursor_peaks_ && (z == abs(max_charge))));
      addSimpleSpectrum_(spectrum, oligo, z * sign, add_precursor, false);
    }

    spectrum.sortByPosition();
  }


  void TheoreticalSpectrumGenerator::getMultipleSpectra(map<Int, PeakSpectrum>& spectra, const NASequence& oligo, const set<Int>& charges, Int base_charge) const
  {
    spectra.clear();
    if (charges.empty()) return;
    bool negative_mode = *charges.begin() < 0;
    bool add_all_precursors = (add_precursor_peaks_ &&
                               add_all_precursor_charges_);
    bool add_final_precursor = (add_precursor_peaks_ &&
                                !add_all_precursor_charges_);

    if (add_metainfo_)
    {
      for (Int charge : charges)
      {
        PeakSpectrum& spectrum = spectra[charge];
        spectrum.getIntegerDataArrays().resize(1);
        spectrum.getIntegerDataArrays()[0].setName("Charges");
        spectrum.getStringDataArrays().resize(1);
        spectrum.getStringDataArrays()[0].setName("IonNames");
      }
    }

    PeakSpectrum uncharged_spectrum = getUnchargedSpectrum_(oligo);

    if (negative_mode)
    {
      if (base_charge > 0) base_charge = -base_charge;
      // in negative mode, charges are ordered high to low - iterate in reverse:
      set<Int>::const_reverse_iterator charge_it = charges.rbegin();
      // skip requested charges that are lower than "base_charge":
      while (*charge_it > base_charge) // ">" because of negative mode
      {
        ++charge_it;
        if (charge_it == charges.rend()) return;
      }
      Int charge = base_charge;
      while (charge_it != charges.rend())
      {
        PeakSpectrum& spectrum = spectra[*charge_it];
        for (; charge >= *charge_it; --charge)
        {
          addChargedSpectrum_(spectrum, uncharged_spectrum, charge,
                              add_all_precursors);
        }
        ++charge_it;
        if (charge_it != charges.rend())
        {
          spectra[*charge_it] = spectrum; // initialize next spectrum
        }
        // if we want precursor peaks only for selected charge states, add them
        // after the next spectrum has been initialized:
        if (add_final_precursor)
        {
          spectrum.push_back(uncharged_spectrum.back());
          double mass = spectrum.back().getMZ() + charge *
            Constants::PROTON_MASS_U;
          spectrum.back().setMZ(abs(mass / charge));
          if (add_metainfo_)
          {
            spectrum.getStringDataArrays()[0].push_back("M");
            spectrum.getIntegerDataArrays()[0].push_back(charge);
          }
        }
        spectrum.sortByPosition();
      }
    }
    else // positive mode
    {
      set<Int>::const_iterator charge_it = charges.begin();
      // skip requested charges that are lower than "base_charge":
      while (*charge_it < base_charge)
      {
        ++charge_it;
        if (charge_it == charges.end()) return;
      }
      Int charge = base_charge;
      while (charge_it != charges.end())
      {
        PeakSpectrum& spectrum = spectra[*charge_it];
        for (; charge <= *charge_it; ++charge)
        {
          addChargedSpectrum_(spectrum, uncharged_spectrum, charge,
                              add_all_precursors);
        }
        ++charge_it;
        if (charge_it != charges.end())
        {
          spectra[*charge_it] = spectrum; // initialize next spectrum
        }
        // if we want precursor peaks only for selected charge states, add them
        // after the next spectrum has been initialized:
        if (add_final_precursor)
        {
          spectrum.push_back(uncharged_spectrum.back());
          double mass = spectrum.back().getMZ() + charge *
            Constants::PROTON_MASS_U;
          spectrum.back().setMZ(mass / charge);
          if (add_metainfo_)
          {
            spectrum.getStringDataArrays()[0].push_back("M");
            spectrum.getIntegerDataArrays()[0].push_back(charge);
          }
        }
        spectrum.sortByPosition();
      }
    }
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
       LOG_ERROR << "Unknown residue type encountered. Can't map to ion letter." << endl;
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
      LOG_ERROR << "Unknown ribonucleotide type encountered. Can't map to ion code." << endl;
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


  void TheoreticalSpectrumGenerator::addIsotopeCluster_(PeakSpectrum& spectrum, const NASequence& ion, NASequence::NASFragmentType res_type, Int charge, double intensity) const
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
        spectrum.getStringDataArrays()[0].push_back(ion_name);
        spectrum.getIntegerDataArrays()[0].push_back(charge);
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


  void TheoreticalSpectrumGenerator::addPeaks_(PeakSpectrum& spectrum, const NASequence& oligo, NASequence::NASFragmentType res_type, Int charge) const
  {
    spectrum.reserve(oligo.size());

    // Generate the ion peaks:
    // Does not generate peaks of full sequence (therefore "<").
    // They are added via precursor mass (and neutral losses).
    // Could be changed in the future.

    double intensity(1);
    static const double methyl_weight = EmpiricalFormula("CH2").getMonoWeight();
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

    if (res_type == NASequence::AminusB || res_type == NASequence::AIon || res_type == NASequence::BIon || res_type == NASequence::CIon || res_type == NASequence::DIon)
    {
      // @TODO: a-B ions ("NASequence::AminusB") may not be relevant for
      // fragments of length 1 (unless modified?)

      if (!add_isotopes_) // add single peak
      {
        Size length = add_first_prefix_ion_ ? 1 : 2;
        for (; length < oligo.size(); ++length)
        {
          NASequence ion = oligo.getPrefix(length);
          double mass = ion.getMonoWeight(res_type, charge);
          Peak1D p;
          p.setMZ(mass / abs(charge));
          p.setIntensity(intensity);
          spectrum.push_back(p);
          String ion_name;
          if (add_metainfo_)
          {
            ion_name = ribonucleotideTypeToIonCode_(res_type, length); // + String((Size)abs(charge), charge_sign);
            spectrum.getStringDataArrays()[0].push_back(ion_name);
            spectrum.getIntegerDataArrays()[0].push_back(charge);
          }
          // special treatment for a-B ions of "ambiguous" modifications:
          // create two peaks with half intensity, representing methyl group
          // lost/retained on backbone:
          if ((res_type == NASequence::AminusB) &&
              ion[length - 1]->isAmbiguous())
          {
            spectrum.back().setIntensity(intensity / 2.0);
            mass += methyl_weight;
            p.setMZ(mass / abs(charge));
            p.setIntensity(intensity / 2.0);
            spectrum.push_back(p);
            if (add_metainfo_)
            {
              spectrum.getStringDataArrays()[0].push_back(ion_name);
              spectrum.getIntegerDataArrays()[0].push_back(charge);
            }
          }
        }
      }
      else // add isotope clusters (slow)
      {
        Size length = add_first_prefix_ion_ ? 1 : 2;
        for (; length < oligo.size(); ++length)
        {
          const NASequence ion = oligo.getPrefix(length);
          addIsotopeCluster_(spectrum, ion, res_type, charge, intensity);
        }
      }

//      if (add_losses_) // add loss peaks (slow)
//      {
//        Size i = add_first_prefix_ion_ ? 1 : 2;
//        for (; i < oligo.size(); ++i)
//        {
//          const NASequence ion = oligo.getPrefix(i);
//          addLosses_(spectrum, ion, ion_names, charges, intensity, res_type, charge);
//        }
//      }
    }

    else // WIon, XIon, YIon, ZIon
    {
      if (!add_isotopes_) // add single peak
      {
        for (Size length = 1; length < oligo.size(); ++length)
        {
          NASequence ion = oligo.getSuffix(length);
          double mass = ion.getMonoWeight(res_type, charge);
          Peak1D p;
          p.setMZ(mass / abs(charge));
          p.setIntensity(intensity);
          spectrum.push_back(p);
          if (add_metainfo_)
          {
            String ion_name = ribonucleotideTypeToIonCode_(res_type, length); // + String((Size)abs(charge), charge_sign);
            spectrum.getStringDataArrays()[0].push_back(ion_name);
            spectrum.getIntegerDataArrays()[0].push_back(charge);
          }
        }
      }
      else // add isotope clusters
      {
        for (Size length = 1; length < oligo.size(); ++length)
        {
          const NASequence ion = oligo.getSuffix(length);
          addIsotopeCluster_(spectrum, ion, res_type, charge, intensity); //TODO IMPLEMENT
        }
      }

      //if (add_losses_) // add loss peaks (slow)
      //{
      //  for (Size i = 1; i < oligo.size(); ++i)
       // {
        //  const NASequence ion = oligo.getSuffix(i);
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
          mono_weight += peptide[i].getMonoWeight(Residue::Internal); // standard internal residue including named modifications: c
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
          mono_weight += peptide[i].getMonoWeight(Residue::Internal); // standard internal residue including named modifications: c
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


  void TheoreticalSpectrumGenerator::addPrecursorPeaks_(PeakSpectrum& spectrum, const NASequence& oligo, Int charge) const
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
    double mono_pos = oligo.getMonoWeight(NASequence::Full, charge);

    if (add_isotopes_)
    {
      IsotopeDistribution dist = oligo.getFormula(NASequence::Full, charge).getIsotopeDistribution(CoarseIsotopePatternGenerator(max_isotope_));
      double j(0.0);
      for (IsotopeDistribution::ConstIterator it = dist.begin(); it != dist.end(); ++it, ++j)
      {
        p.setMZ((double)(mono_pos + j * Constants::C13C12_MASSDIFF_U * charge_unit) / (double)(abs(charge)));
        p.setIntensity(pre_int_ * it->getIntensity());
        if (add_metainfo_)
        {
          spectrum.getStringDataArrays()[0].push_back(ion_name);
          spectrum.getIntegerDataArrays()[0].push_back(charge);
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
        spectrum.getStringDataArrays()[0].push_back(ion_name);
        spectrum.getIntegerDataArrays()[0].push_back(charge);
      }
      spectrum.push_back(p);
    }
    // loss peaks of the precursor

    //loss of water
    EmpiricalFormula ion = oligo.getFormula(NASequence::Full, charge) - EmpiricalFormula("H2O");
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
          spectrum.getStringDataArrays()[0].push_back(ion_name);
          spectrum.getIntegerDataArrays()[0].push_back(charge);
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
        spectrum.getStringDataArrays()[0].push_back(ion_name);
        spectrum.getIntegerDataArrays()[0].push_back(charge);
      }
      spectrum.push_back(p);
    }

    //loss of ammonia: not applicable for nucleic acid sequences
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
