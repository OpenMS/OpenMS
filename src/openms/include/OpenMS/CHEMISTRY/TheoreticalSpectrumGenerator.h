// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2016.
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
// $Maintainer: Sandro Andreotti $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#ifndef OPENMS_CHEMISTRY_THEORETICALSPECTRUMGENERATOR_H
#define OPENMS_CHEMISTRY_THEORETICALSPECTRUMGENERATOR_H

#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/KERNEL/StandardTypes.h>

namespace OpenMS
{
  class AASequence;

  /**
      @brief Generates theoretical spectra with various options

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
    virtual ~TheoreticalSpectrumGenerator();
    //@}

    /// assignment operator
    TheoreticalSpectrumGenerator & operator=(const TheoreticalSpectrumGenerator & tsg);

    /** @name Acessors
     */
    //@{
    /// returns a spectrum with b and y peaks
    virtual void getSpectrum(RichPeakSpectrum & spec, const AASequence & peptide, Int charge = 1) const;

    /// adds peaks to a spectrum of the given ion-type, peptide, charge, and intensity
    virtual void addPeaks(RichPeakSpectrum & spectrum, const AASequence & peptide, Residue::ResidueType res_type, Int charge = 1) const;

    /// adds the precursor peaks to the spectrum
    virtual void addPrecursorPeaks(RichPeakSpectrum & spec, const AASequence & peptide, Int charge = 1) const;

    /// Adds the common, most abundant immonium ions to the theoretical spectra
    void addAbundantImmoniumIons(RichPeakSpectrum & spec) const;

    /// overwrite
    void updateMembers_();

    //@}

    protected:
      /// helper to add an isotope cluster to a spectrum
      void addIsotopeCluster_(RichPeakSpectrum & spectrum, const AASequence & ion, Residue::ResidueType res_type, Int charge, double intensity) const;

      /// helper to add a single peak to a spectrum
      void addPeak_(RichPeakSpectrum & spectrum, double pos, double intensity, Residue::ResidueType res_type, Size ion_index, int charge) const;
   
      /// helper for mapping residue type to letter
      char residueTypeToIonLetter_(Residue::ResidueType res_type) const;

      /// helper to add full neutral loss ladders
      void addLosses_(RichPeakSpectrum & spectrum, const AASequence & ion, double intensity, Residue::ResidueType res_type, int charge) const;

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
      bool add_precursor_peaks;
      bool add_abundant_immonium_ions;
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
#endif

