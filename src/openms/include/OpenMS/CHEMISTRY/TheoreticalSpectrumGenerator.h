// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg, Eugen Netz $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/METADATA/DataArrays.h>
#include <OpenMS/KERNEL/MSSpectrum.h>

namespace OpenMS
{
  class AASequence;

  /**
      @brief Generates theoretical spectra for peptides with various options

      If the tool parameter add_metainfo is set to true,
      ion names like y8+ or [M-H2O+2H]++ are written as strings in a StringDataArray with the name 
      corresponding to the constant Constants::UserParam::IonNames
      and charges are written as integers in an IntegerDataArray with the name "Charges"
      in the returned PeakSpectrum.

      The getSpectrum function can be called with the same PeakSpectrum multiple times to add additional peaks.
      If the PeakSpectrum already has DataArrays, then the first StringDataArray and the first IntegerDataArray
      are extended. Therefore it is not recommended to add to or change the PeakSpectrum or these DataArrays
      between calls of the getSpectrum function with the same PeakSpectrum.

      @note The generation of neutral loss peaks is very slow in this class.
      Something similar to the neutral loss precalculation used in TheoreticalSpectrumGeneratorXLMS
      should be implemented here as well.

      @htmlinclude OpenMS_TheoreticalSpectrumGenerator.parameters

      @ingroup Chemistry
  */
  class OPENMS_DLLAPI TheoreticalSpectrumGenerator :
    public DefaultParamHandler
  {
    public:

    /** @name Constructors and Destructors
    */
    //@{
    /// default constructor
    TheoreticalSpectrumGenerator();

    /// copy constructor
    TheoreticalSpectrumGenerator(const TheoreticalSpectrumGenerator& source);

    /// destructor
    ~TheoreticalSpectrumGenerator() override;
    //@}

    /// assignment operator
    TheoreticalSpectrumGenerator& operator=(const TheoreticalSpectrumGenerator& tsg);

    /** @name Acessors
     */
    //@{
    /// Generates a spectrum for a peptide sequence, with the ion types that are set in the tool parameters
    /// If precursor_charge is set to 0 max_charge + 1 will be used.
    /// @throw Exception::InvalidParameter   if precursor_charge < max_charge
    virtual void getSpectrum(PeakSpectrum& spec, const AASequence& peptide, Int min_charge, Int max_charge, Int precursor_charge = 0) const;

    /// Generates a spectrum for a peptide sequence based on activation method and precursor charge.
    /// Activation method 'CID' or 'HCID' will generate only b- and y-ions.
    /// Activation method 'ECD' or 'ETD' will generate only c- and z-ions.
    /// If precursor charge is greater than 2 ions with charge 1 & 2 will be generated, else only ions of charge 1 will appear in the spectrum.
    /// @throw Exception::InvalidParameter   If fragmentation method is anything else than 'CID', 'HCID', 'ECD' or 'ETD'.
    static MSSpectrum generateSpectrum(const Precursor::ActivationMethod& fm, const AASequence& seq, int precursor_charge);

    /// overwrite
    void updateMembers_() override;
    //@}

    protected:

    /// adds peaks to a spectrum of the given ion-type, peptide, charge, and intensity, also adds charges and ion names to the DataArrays, if the add_metainfo parameter is set to true
    virtual void addPeaks_(PeakSpectrum& spectrum, const AASequence& peptide, DataArrays::StringDataArray& ion_names, DataArrays::IntegerDataArray& charges, MSSpectrum::Chunks& chunks, const Residue::ResidueType res_type, Int charge = 1) const;

    /// adds the precursor peaks to the spectrum, also adds charges and ion names to the DataArrays, if the add_metainfo parameter is set to true
    virtual void addPrecursorPeaks_(PeakSpectrum& spec, const AASequence& peptide, DataArrays::StringDataArray& ion_names, DataArrays::IntegerDataArray& charges, Int charge = 1) const;

    /// Adds the common, most abundant immonium ions to the theoretical spectra if the residue is contained in the peptide sequence, also adds charges and ion names to the DataArrays, if the add_metainfo parameter is set to true
    void addAbundantImmoniumIons_(PeakSpectrum& spec, const AASequence& peptide, DataArrays::StringDataArray& ion_names, DataArrays::IntegerDataArray& charges) const;

    /// helper to add an isotope cluster to a spectrum, also adds charges and ion names to the DataArrays, if the add_metainfo parameter is set to true
    void addIsotopeCluster_(PeakSpectrum& spectrum, const AASequence& ion, DataArrays::StringDataArray& ion_names, DataArrays::IntegerDataArray& charges, const Residue::ResidueType res_type, Int charge, double intensity) const;

    /// helper to add full neutral loss ladders (for isotope clusters), also adds charges and ion names to the DataArrays, if the add_metainfo parameter is set to true
    void addLosses_(PeakSpectrum& spectrum, const AASequence& ion, DataArrays::StringDataArray& ion_names, DataArrays::IntegerDataArray& charges, double intensity, const Residue::ResidueType res_type, int charge) const;

    /// helper to add full neutral loss ladders (for single peaks), also adds charges and ion names to the DataArrays, if the add_metainfo parameter is set to true
    void addLossesFaster_(PeakSpectrum& spectrum, double mz, const std::set<EmpiricalFormula>& f_losses, int ion_ordinal, DataArrays::StringDataArray& ion_names, DataArrays::IntegerDataArray& charges, const std::map<EmpiricalFormula, String>& formula_str_cache, double intensity, const Residue::ResidueType res_type, bool add_metainfo, int charge) const;

    bool add_b_ions_;
    bool add_y_ions_;
    bool add_a_ions_;
    bool add_c_ions_;
    bool add_x_ions_;
    bool add_z_ions_;
    bool add_zp1_ions_;
    bool add_zp2_ions_;
    bool add_first_prefix_ion_;
    bool add_losses_;
    bool add_metainfo_;
    bool add_isotopes_;
    int isotope_model_;
    bool add_precursor_peaks_;
    bool add_all_precursor_charges_ ;
    bool add_abundant_immonium_ions_;
    bool sort_by_position_;
    double a_intensity_;
    double b_intensity_;
    double c_intensity_;
    double x_intensity_;
    double y_intensity_;
    double z_intensity_;

    Int max_isotope_;
    double rel_loss_intensity_;
    double max_isotope_probability_;
    double pre_int_;
    double pre_int_H2O_;
    double pre_int_NH3_;
  };
}
