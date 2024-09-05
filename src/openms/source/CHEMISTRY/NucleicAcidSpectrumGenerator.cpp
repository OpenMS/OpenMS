// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/NucleicAcidSpectrumGenerator.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/LogStream.h>

using namespace std;

namespace OpenMS
{

  NucleicAcidSpectrumGenerator::NucleicAcidSpectrumGenerator() : DefaultParamHandler("NucleicAcidSpectrumGenerator")
  {
    defaults_.setValue("add_metainfo", "false", "Adds the type of peaks as meta information to the peaks, e.g. c1, y2, a3-B");
    defaults_.setValidStrings("add_metainfo", {"true", "false"});

    defaults_.setValue("add_precursor_peaks", "false", "Adds peaks of the unfragmented precursor ion to the spectrum");
    defaults_.setValidStrings("add_precursor_peaks", {"true", "false"});

    defaults_.setValue("add_all_precursor_charges", "false", "Adds precursor peaks with all charges in the given range");
    defaults_.setValidStrings("add_all_precursor_charges", {"true", "false"});

    defaults_.setValue("add_first_prefix_ion", "false", "If set to true a1, b1, ..., z1 ions are added");
    defaults_.setValidStrings("add_first_prefix_ion", {"true", "false"});

    defaults_.setValue("add_a_ions", "false", "Add peaks of a-ions to the spectrum");
    defaults_.setValidStrings("add_a_ions", {"true", "false"});

    defaults_.setValue("add_b_ions", "true", "Add peaks of b-ions to the spectrum");
    defaults_.setValidStrings("add_b_ions", {"true", "false"});

    defaults_.setValue("add_c_ions", "false", "Add peaks of c-ions to the spectrum");
    defaults_.setValidStrings("add_c_ions", {"true", "false"});

    defaults_.setValue("add_d_ions", "false", "Add peaks of d-ions to the spectrum"); // only for nucleotide sequences
    defaults_.setValidStrings("add_d_ions", {"true", "false"});

    defaults_.setValue("add_w_ions", "false", "Add peaks of w-ions to the spectrum"); // only for nucleotide sequences
    defaults_.setValidStrings("add_w_ions", {"true", "false"});

    defaults_.setValue("add_x_ions", "false", "Add peaks of  x-ions to the spectrum");
    defaults_.setValidStrings("add_x_ions", {"true", "false"});

    defaults_.setValue("add_y_ions", "true", "Add peaks of y-ions to the spectrum");
    defaults_.setValidStrings("add_y_ions", {"true", "false"});

    defaults_.setValue("add_z_ions", "false", "Add peaks of z-ions to the spectrum");
    defaults_.setValidStrings("add_z_ions", {"true", "false"});

    defaults_.setValue("add_a-B_ions", "false", "Add peaks of a-B-ions to the spectrum"); // only for nucleotide sequences
    defaults_.setValidStrings("add_a-B_ions", {"true", "false"});

    // intensity options of the ions
    defaults_.setValue("a_intensity", 1.0, "Intensity of the a-ions");
    defaults_.setValue("b_intensity", 1.0, "Intensity of the b-ions");
    defaults_.setValue("c_intensity", 1.0, "Intensity of the c-ions");
    defaults_.setValue("d_intensity", 1.0, "Intensity of the d-ions");
    defaults_.setValue("w_intensity", 1.0, "Intensity of the w-ions");
    defaults_.setValue("x_intensity", 1.0, "Intensity of the x-ions");
    defaults_.setValue("y_intensity", 1.0, "Intensity of the y-ions");
    defaults_.setValue("z_intensity", 1.0, "Intensity of the z-ions");
    defaults_.setValue("a-B_intensity", 1.0, "Intensity of the a-B-ions");

    // precursor intensity
    defaults_.setValue("precursor_intensity", 1.0, "Intensity of the precursor peak");

    defaultsToParam_();
  }


  NucleicAcidSpectrumGenerator::NucleicAcidSpectrumGenerator(const NucleicAcidSpectrumGenerator& source) : DefaultParamHandler(source)
  {
  }


  NucleicAcidSpectrumGenerator& NucleicAcidSpectrumGenerator::operator=(const NucleicAcidSpectrumGenerator& source)
  {
    DefaultParamHandler::operator=(source);
    return *this;
  }


  NucleicAcidSpectrumGenerator::~NucleicAcidSpectrumGenerator() = default;


  void NucleicAcidSpectrumGenerator::addFragmentPeaks_(MSSpectrum& spectrum, const vector<double>& fragment_masses, const String& ion_type, double offset, double intensity, Size start) const
  {
    for (Size i = start; i < fragment_masses.size(); ++i)
    {
      Peak1D peak(fragment_masses[i] + offset, intensity);
      spectrum.push_back(peak);
    }
    if (add_metainfo_)
    {
      for (Size i = start; i < fragment_masses.size(); ++i)
      {
        String ion_name = ion_type + String(i + 1);
        spectrum.getStringDataArrays()[0].push_back(ion_name);
      }
    }
  }


  void NucleicAcidSpectrumGenerator::addAMinusBPeaks_(MSSpectrum& spectrum, const vector<double>& fragment_masses, const NASequence& oligo, Size start) const
  {
    // offset: phosphate (from bond) minus 3 water (from various reactions)
    static const double offset = EmpiricalFormula("H-5P").getMonoWeight();
    // offset for first ("a1-B") ion: loss of 2 water
    static const double initial_offset = -EmpiricalFormula("H4O2").getMonoWeight();
    // methyl group may be retained on ribose for "ambiguous" mods:
    static const double methyl_mass = EmpiricalFormula("CH2").getMonoWeight();

    for (Size i = start; i < fragment_masses.size(); ++i)
    {
      double mass = oligo[i]->getBaselossFormula().getMonoWeight();
      if (i > 0)
      {
        // base at position "i" is lost, so use fragment up to pos. "i - 1":
        mass += fragment_masses[i - 1] + offset;
        // check if the offset should be thiol or not
        if (oligo[i-1]->getCode().back() == '*')
        {
          mass += EmpiricalFormula("SO-1").getMonoWeight();
        }
      }
      else // first ribonucleotide
      {
        mass += initial_offset;
      }
      Peak1D peak(mass, aB_intensity_);
      if (oligo[i]->isAmbiguous())
      {
        // special treatment for a-B ions of "ambiguous" modifications:
        // create two peaks with half intensity, representing methyl group
        // lost/retained on backbone:
        peak.setIntensity(aB_intensity_ * 0.5);
        spectrum.push_back(peak);
        mass += methyl_mass;
        peak.setMZ(mass);
      }
      spectrum.push_back(peak);
    }
    if (add_metainfo_)
    {
      for (Size i = start; i < fragment_masses.size(); ++i)
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


  MSSpectrum NucleicAcidSpectrumGenerator::getUnchargedSpectrum_(const NASequence& oligo) const
  {
    static const double H_mass = EmpiricalFormula("H").getMonoWeight();
    // phosphate minus water:
    static const double backbone_mass = EmpiricalFormula("H-1PO2").getMonoWeight();

    static const double a_ion_offset = -EmpiricalFormula("H2O").getMonoWeight();
    static const double b_ion_offset = 0.0;
    static const double c_ion_offset = backbone_mass;
    static const double d_ion_offset = EmpiricalFormula("HPO3").getMonoWeight();
    static const double w_ion_offset = d_ion_offset;
    static const double x_ion_offset = c_ion_offset;
    static const double y_ion_offset = b_ion_offset;
    static const double z_ion_offset = a_ion_offset;

    // a a-B w x ions have different offsets if we have phosphorothioate linkages,

    MSSpectrum spectrum;
    if (oligo.empty())
    {
      return spectrum;
    }
    double three_prime_mass = 0.0, five_prime_mass = 0.0;
    if (oligo.getThreePrimeMod() != nullptr)
    {
      three_prime_mass = oligo.getThreePrimeMod()->getMonoMass() - H_mass;
    }
    if (oligo.getFivePrimeMod() != nullptr)
    {
      five_prime_mass = oligo.getFivePrimeMod()->getMonoMass() - H_mass;
    }

    vector<double> ribo_masses(oligo.size());
    // Create a vector of doubles to represent the phosorothioate linkage mass shift
    vector<double> thiols(oligo.size(), 0.0);
    Size index = 0;
    for (const auto& ribo : oligo)
    {
      ribo_masses[index] = ribo.getMonoMass();
      // * at the end means phosphorothioate
      if (ribo.getCode().back() == '*')
      {
        thiols[index] = EmpiricalFormula("SO-1").getMonoWeight();
      }
      ++index;
    }

    spectrum.getStringDataArrays().resize(1);
    spectrum.getStringDataArrays()[0].setName(Constants::UserParam::IonNames);

    vector<double> fragments_left, fragments_right, thiol_offsets;
    Size start = add_first_prefix_ion_ ? 0 : 1;
    // Drop the final thiol, 'cause its not linking anything
    thiols.resize(oligo.size() - 1);
    if ((add_a_ions_ || add_b_ions_ || add_c_ions_ || add_d_ions_ || add_aB_ions_) && (oligo.size() > start + 1))
    {
      fragments_left.resize(oligo.size() - 1);

      fragments_left[0] = ribo_masses[0] + five_prime_mass;
      for (Size i = 1; i < oligo.size() - 1; ++i)
      {
        fragments_left[i] = (fragments_left[i - 1] + ribo_masses[i] + backbone_mass + thiols[i - 1]);
      }
      // with thiols c and d ions have a 15.99 mass shift, we calculated that above now we add it
      vector<double> frag_l_thiol(fragments_left.size());
      std::transform(fragments_left.begin(), fragments_left.end(), thiols.begin(), frag_l_thiol.begin(), std::plus<double>());

      if (add_a_ions_)
      {
        addFragmentPeaks_(spectrum, fragments_left, "a", a_ion_offset, a_intensity_, start);
      }
      if (add_b_ions_)
      {
        addFragmentPeaks_(spectrum, fragments_left, "b", b_ion_offset, b_intensity_, start);
      }
      if (add_c_ions_)
      {
        addFragmentPeaks_(spectrum, frag_l_thiol, "c", c_ion_offset, c_intensity_, start);
      }
      if (add_d_ions_)
      {
        addFragmentPeaks_(spectrum, frag_l_thiol, "d", d_ion_offset, d_intensity_, start);
      }
      if (add_aB_ions_) // special case
      {
        addAMinusBPeaks_(spectrum, fragments_left, oligo, start);
      }
    }

    if ((add_w_ions_ || add_x_ions_ || add_y_ions_ || add_z_ions_) && (oligo.size() > 1))
    {
      fragments_right.resize(oligo.size() - 1);
      fragments_right[0] = ribo_masses.back() + three_prime_mass;
      for (Size i = 1; i < oligo.size() - 1; ++i)
      {
        Size ribo_index = oligo.size() - i - 1;
        fragments_right[i] = (fragments_right[i - 1] + ribo_masses[ribo_index] + backbone_mass + thiols[ribo_index]);
      }
      // with thiols a and b ions have a 15.99 mass shift, we calculated that above now we add it
      vector<double> frag_r_thiol(fragments_right.size());
      std::reverse(thiols.begin(), thiols.end()); // Reverse, since we go from the other side
      std::transform(fragments_right.begin(), fragments_right.end(), thiols.begin(), frag_r_thiol.begin(), std::plus<double>());

      if (add_w_ions_)
      {
        addFragmentPeaks_(spectrum, frag_r_thiol, "w", w_ion_offset, w_intensity_);
      }
      if (add_x_ions_)
      {
        addFragmentPeaks_(spectrum, frag_r_thiol, "x", x_ion_offset, x_intensity_);
      }
      if (add_y_ions_)
      {
        addFragmentPeaks_(spectrum, fragments_right, "y", y_ion_offset, y_intensity_);
      }
      if (add_z_ions_)
      {
        addFragmentPeaks_(spectrum, fragments_right, "z", z_ion_offset, z_intensity_);
      }
    }

    if (add_precursor_peaks_) // re-use what we've already calculated
    {
      Peak1D peak(0.0, precursor_intensity_);
      bool have_left = !fragments_left.empty();
      bool have_right = !fragments_right.empty();
      if (have_left && have_right)
      {
        peak.setMZ(fragments_left[0] + fragments_right.back() + backbone_mass);
      }
      else if (have_left)
      {
        peak.setMZ(fragments_left.back() + ribo_masses.back() + backbone_mass + three_prime_mass);
      }
      else if (have_right)
      {
        peak.setMZ(fragments_right.back() + ribo_masses[0] + backbone_mass + five_prime_mass);
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


  void NucleicAcidSpectrumGenerator::addChargedSpectrum_(MSSpectrum& spectrum, const MSSpectrum& uncharged_spectrum, Int charge, bool add_precursor) const
  {
    if (uncharged_spectrum.empty())
      return;
    Size size = uncharged_spectrum.size();
    if (add_precursor_peaks_ && !add_precursor)
    {
      --size; // uncharged spectrum contains precursor peak - exclude it
    }
    for (Size i = 0; i < size; ++i)
    {
      spectrum.push_back(uncharged_spectrum[i]);
      spectrum.back().setMZ(std::fabs(spectrum.back().getMZ() / charge + Constants::PROTON_MASS_U));
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


  void NucleicAcidSpectrumGenerator::getSpectrum(MSSpectrum& spectrum, const NASequence& oligo, Int min_charge, Int max_charge) const
  {
    Int sign = 1;
    if (max_charge < 0 && min_charge < 0) // negative mode
    {
      sign = -1;
    }
    else if (max_charge * min_charge < 0)
    {
      // Signs don't match - we need to quit and throw error here to avoid messing up for loops below
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "min. and max. charge must both be either positive or negative");
    }
    if (abs(max_charge) < abs(min_charge))
    {
      swap(max_charge, min_charge);
    }

    if (add_metainfo_)
    {
      // @TODO: what if arrays already exist, but contain different data?
      if (spectrum.getIntegerDataArrays().empty())
      {
        spectrum.getIntegerDataArrays().resize(1);
        spectrum.getIntegerDataArrays()[0].setName("Charges");
      }
      if (spectrum.getStringDataArrays().empty())
      {
        spectrum.getStringDataArrays().resize(1);
        spectrum.getStringDataArrays()[0].setName(Constants::UserParam::IonNames);
      }
    }

    MSSpectrum uncharged_spectrum = getUnchargedSpectrum_(oligo);

    for (UInt z = (UInt)abs(min_charge); z <= (UInt)abs(max_charge) && z < (UInt)oligo.size(); ++z)
    {
      bool add_precursor = ((add_precursor_peaks_ && add_all_precursor_charges_) || (add_precursor_peaks_ && (z == (UInt)abs(max_charge))));
      addChargedSpectrum_(spectrum, uncharged_spectrum, z * sign, add_precursor);
    }

    spectrum.sortByPosition();
  }


  void NucleicAcidSpectrumGenerator::getMultipleSpectra(map<Int, MSSpectrum>& spectra, const NASequence& oligo, const set<Int>& charges, Int base_charge) const
  {
    spectra.clear();
    if (charges.empty())
      return;
    bool negative_mode = *charges.begin() < 0;
    bool add_all_precursors = (add_precursor_peaks_ && add_all_precursor_charges_);
    bool add_final_precursor = (add_precursor_peaks_ && !add_all_precursor_charges_);

    if (add_metainfo_)
    {
      for (Int charge : charges)
      {
        MSSpectrum& spectrum = spectra[charge];
        spectrum.getIntegerDataArrays().resize(1);
        spectrum.getIntegerDataArrays()[0].setName("Charges");
        spectrum.getStringDataArrays().resize(1);
        spectrum.getStringDataArrays()[0].setName(Constants::UserParam::IonNames);
      }
    }

    MSSpectrum uncharged_spectrum = getUnchargedSpectrum_(oligo);

    if (negative_mode)
    {
      if (base_charge > 0)
        base_charge = -base_charge;
      // in negative mode, charges are ordered high to low - iterate in reverse:
      set<Int>::const_reverse_iterator charge_it = charges.rbegin();
      // skip requested charges that are lower than "base_charge":
      while (*charge_it > base_charge) // ">" because of negative mode
      {
        ++charge_it;
        if (charge_it == charges.rend())
          return;
      }
      Int charge = base_charge;
      while (charge_it != charges.rend())
      {
        MSSpectrum& spectrum = spectra[*charge_it];
        for (; charge >= *charge_it; --charge)
        {
          addChargedSpectrum_(spectrum, uncharged_spectrum, charge, add_all_precursors);
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
          spectrum.back().setMZ(std::fabs(spectrum.back().getMZ() / charge + Constants::PROTON_MASS_U));
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
        if (charge_it == charges.end())
          return;
      }
      Int charge = base_charge;
      while (charge_it != charges.end())
      {
        MSSpectrum& spectrum = spectra[*charge_it];
        for (; charge <= *charge_it; ++charge)
        {
          addChargedSpectrum_(spectrum, uncharged_spectrum, charge, add_all_precursors);
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
          spectrum.back().setMZ(spectrum.back().getMZ() / charge + Constants::PROTON_MASS_U);
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


  void NucleicAcidSpectrumGenerator::updateMembers_()
  {
    add_a_ions_ = param_.getValue("add_a_ions").toBool();
    add_b_ions_ = param_.getValue("add_b_ions").toBool();
    add_c_ions_ = param_.getValue("add_c_ions").toBool();
    add_d_ions_ = param_.getValue("add_d_ions").toBool();
    add_w_ions_ = param_.getValue("add_w_ions").toBool();
    add_x_ions_ = param_.getValue("add_x_ions").toBool();
    add_y_ions_ = param_.getValue("add_y_ions").toBool();
    add_z_ions_ = param_.getValue("add_z_ions").toBool();
    add_aB_ions_ = param_.getValue("add_a-B_ions").toBool();
    add_first_prefix_ion_ = param_.getValue("add_first_prefix_ion").toBool();
    add_metainfo_ = param_.getValue("add_metainfo").toBool();
    add_precursor_peaks_ = param_.getValue("add_precursor_peaks").toBool();
    add_all_precursor_charges_ = param_.getValue("add_all_precursor_charges").toBool();
    a_intensity_ = (double)param_.getValue("a_intensity");
    b_intensity_ = (double)param_.getValue("b_intensity");
    c_intensity_ = (double)param_.getValue("c_intensity");
    d_intensity_ = (double)param_.getValue("d_intensity");
    w_intensity_ = (double)param_.getValue("w_intensity");
    x_intensity_ = (double)param_.getValue("x_intensity");
    y_intensity_ = (double)param_.getValue("y_intensity");
    z_intensity_ = (double)param_.getValue("z_intensity");
    aB_intensity_ = (double)param_.getValue("a-B_intensity");
    precursor_intensity_ = (double)param_.getValue("precursor_intensity");
  }

} // end namespace OpenMS
