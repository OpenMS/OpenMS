// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/METADATA/DataArrays.h>
#include <OpenMS/CHEMISTRY/NASequence.h>

namespace OpenMS
{
  class NASequence;

  /**
      @brief Generates theoretical spectra for nucleic acid sequences

      @htmlinclude OpenMS_NucleicAcidSpectrumGenerator.parameters

      @ingroup Chemistry
  */
  class OPENMS_DLLAPI NucleicAcidSpectrumGenerator :
    public DefaultParamHandler
  {
    public:

    /** @name Constructors and Destructors
    */
    //@{
    /// default constructor
    NucleicAcidSpectrumGenerator();

    /// copy constructor
    NucleicAcidSpectrumGenerator(const NucleicAcidSpectrumGenerator& source);

    /// destructor
    ~NucleicAcidSpectrumGenerator() override;
    //@}

    /// assignment operator
    NucleicAcidSpectrumGenerator& operator=(const NucleicAcidSpectrumGenerator& source);

    /** @name Acessors
     */
    //@{
    /// Generates a spectrum for an oligonucleotide sequence, with the ion types that are set in the tool parameters
    void getSpectrum(MSSpectrum& spectrum, const NASequence& oligo, Int min_charge, Int max_charge) const;

    /**
       @brief Generates spectra in multiple charge states for an oligonucleotide sequence

       @param spectra Output spectra
       @param oligo Target oligonucleotide sequence
       @param charges Set of charge states to generate
       @param base_charge Minimum charge for peaks in each spectrum

       One spectrum per element in @p charges is generated in @p spectra.

       All values in @p charges must be either positive or negative.

       This function is more efficient than calling getSpectrum() multiple times, because spectra of lower charge states are reused.
    */
    void getMultipleSpectra(std::map<Int, MSSpectrum>& spectra, const NASequence& oligo, const std::set<Int>& charges, Int base_charge = 1) const;

    /// overwrite
    void updateMembers_() override;
    //@}

    protected:

    /// Helper function to add (uncharged) fragment peaks to a spectrum
    void addFragmentPeaks_(MSSpectrum& spectrum, const std::vector<double>& fragment_masses, const String& ion_type, double offset, double intensity, Size start = 0) const;

    /// Special version of addFragmentPeaks_() for a-B ions
    void addAMinusBPeaks_(MSSpectrum& spectrum, const std::vector<double>& fragment_masses, const NASequence& oligo, Size start = 0) const;

    /// Generates a spectrum containing peaks for uncharged fragment masses
    MSSpectrum getUnchargedSpectrum_(const NASequence& oligo) const;

    /// Adds a charged version of an uncharged spectrum to another spectrum
    void addChargedSpectrum_(MSSpectrum& spectrum, const MSSpectrum& uncharged_spectrum, Int charge, bool add_precursor) const;

    bool add_a_ions_;
    bool add_b_ions_;
    bool add_c_ions_;
    bool add_d_ions_;
    bool add_w_ions_;
    bool add_x_ions_;
    bool add_y_ions_;
    bool add_z_ions_;
    bool add_aB_ions_;
    bool add_first_prefix_ion_;
    bool add_metainfo_;
    bool add_precursor_peaks_;
    bool add_all_precursor_charges_ ;
    double a_intensity_;
    double b_intensity_;
    double c_intensity_;
    double d_intensity_;
    double w_intensity_;
    double x_intensity_;
    double y_intensity_;
    double z_intensity_;
    double aB_intensity_;
    double precursor_intensity_;
  };
}
