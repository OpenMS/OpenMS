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
// $Maintainer: Timo Sachsenberg, Eugen Netz $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#ifndef OPENMS_CHEMISTRY_THEORETICALSPECTRUMGENERATOR_H
#define OPENMS_CHEMISTRY_THEORETICALSPECTRUMGENERATOR_H

#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/METADATA/DataArrays.h>

namespace OpenMS
{
  class AASequence;

  /**
      @brief Generates theoretical spectra with various options

      If the tool parameter add_metainfo is set to true,
      ion names like y8+ or [M-H2O+2H]++ are written as strings in a StringDataArray with the name "IonNames"
      and charges are written as integers in an IntegerDataArray with the name "Charges"
      in the returned PeakSpectrum.

      The getSpectrum function can be called with the same PeakSpectrum multiple times to add additional peaks.
      If the PeakSpectrum already has DataArrays, then the first StringDataArray and the first IntegerDataArray
      are extended. Therefore it is not recommended to add to or change the PeakSpectrum or these DataArrays
      between calls of the getSpectrum function with the same PeakSpectrum.

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
    TheoreticalSpectrumGenerator(const TheoreticalSpectrumGenerator & source);

    /// destructor
    ~TheoreticalSpectrumGenerator() override;
    //@}

    /// assignment operator
    TheoreticalSpectrumGenerator & operator=(const TheoreticalSpectrumGenerator & tsg);

    /** @name Acessors
     */
    //@{
    /// returns a spectrum with the ion types, that are set in the tool parameters
    virtual void getSpectrum(PeakSpectrum & spec, const AASequence & peptide, Int min_charge, Int max_charge) const;

    /// overwrite
    void updateMembers_() override;

    //@}

    protected:
      /// adds peaks to a spectrum of the given ion-type, peptide, charge, and intensity, also adds charges and ion names to the DataArrays, if the add_metainfo parameter is set to true
      virtual void addPeaks_(PeakSpectrum & spectrum, const AASequence & peptide, DataArrays::StringDataArray& ion_names, DataArrays::IntegerDataArray& charges, Residue::ResidueType res_type, Int charge = 1) const;

      /// adds the precursor peaks to the spectrum, also adds charges and ion names to the DataArrays, if the add_metainfo parameter is set to true
      virtual void addPrecursorPeaks_(PeakSpectrum & spec, const AASequence & peptide, DataArrays::StringDataArray& ion_names, DataArrays::IntegerDataArray& charges, Int charge = 1) const;

      /// Adds the common, most abundant immonium ions to the theoretical spectra if the residue is contained in the peptide sequence, also adds charges and ion names to the DataArrays, if the add_metainfo parameter is set to true
      void addAbundantImmoniumIons_(PeakSpectrum & spec, const AASequence& peptide, DataArrays::StringDataArray& ion_names, DataArrays::IntegerDataArray& charges) const;

      /// helper to add an isotope cluster to a spectrum, also adds charges and ion names to the DataArrays, if the add_metainfo parameter is set to true
      void addIsotopeCluster_(PeakSpectrum & spectrum, const AASequence & ion, DataArrays::StringDataArray& ion_names, DataArrays::IntegerDataArray& charges, Residue::ResidueType res_type, Int charge, double intensity) const;

      /// helper for mapping residue type to letter
      char residueTypeToIonLetter_(Residue::ResidueType res_type) const;

      /// helper to add full neutral loss ladders, also adds charges and ion names to the DataArrays, if the add_metainfo parameter is set to true
      void addLosses_(PeakSpectrum & spectrum, const AASequence & ion, DataArrays::StringDataArray& ion_names, DataArrays::IntegerDataArray& charges, double intensity, Residue::ResidueType res_type, int charge) const;

      bool add_b_ions_;
      bool add_y_ions_;
      bool add_a_ions_;
      bool add_c_ions_;
      bool add_x_ions_;
      bool add_z_ions_;
      bool add_first_prefix_ion_;
      bool add_losses_;
      bool add_metainfo_;
      bool add_isotopes_;
      bool add_precursor_peaks_;
      bool add_all_precursor_charges_ ;
      bool add_abundant_immonium_ions_;
      double a_intensity_;
      double b_intensity_;
      double c_intensity_;
      double x_intensity_;
      double y_intensity_;
      double z_intensity_;
      Int max_isotope_;
      double rel_loss_intensity_;
      double pre_int_;
      double pre_int_H2O_;
      double pre_int_NH3_;
  };
}
#endif // THEORETICALSPECTRUMGENERATORRPLESS_H
