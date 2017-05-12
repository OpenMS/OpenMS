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
// $Maintainer: Eugen Netz $
// $Authors: Eugen Netz $
// --------------------------------------------------------------------------


#ifndef THEORETICALSPECTRUMGENERATORXLMS_H
#define THEORETICALSPECTRUMGENERATORXLMS_H

#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/METADATA/DataArrays.h>


namespace OpenMS
{
  class AASequence;

  /**
      @brief Generates theoretical spectra for cross-linked peptides with various options

  @htmlinclude OpenMS_TheoreticalSpectrumGeneratorXLMS.parameters

      @ingroup Chemistry
  */
  class OPENMS_DLLAPI TheoreticalSpectrumGeneratorXLMS :
    public DefaultParamHandler
  {
public:

    /** @name Constructors and Destructors
    */
    //@{
    /// default constructor
    TheoreticalSpectrumGeneratorXLMS();

    /// copy constructor
    TheoreticalSpectrumGeneratorXLMS(const TheoreticalSpectrumGeneratorXLMS & source);

    /// destructor
    virtual ~TheoreticalSpectrumGeneratorXLMS();
    //@}

    /// assignment operator
    TheoreticalSpectrumGeneratorXLMS & operator=(const TheoreticalSpectrumGeneratorXLMS & tsg);

    /** @name Acessors
     */
    //@{
    /// returns a spectrum with b and y peaks
    /// TODO crashes, if called without generating metainfo
    virtual void getCommonIonSpectrum(PeakSpectrum & spectrum, AASequence peptide, Size link_pos, bool frag_alpha, int charge = 1, Size link_pos_2 = 0) const;

    /// adds peaks to a spectrum of the given ion-type, peptide, charge, and intensity
    virtual void addCommonPeaks(PeakSpectrum & spectrum, DataArrays::IntegerDataArray & integer_array, DataArrays::StringDataArray & string_array, AASequence peptide, Size link_pos, bool frag_alpha, Residue::ResidueType res_type, int charge = 1, Size link_pos_2 = 0) const;

    virtual void getXLinkIonSpectrum(PeakSpectrum & spectrum, AASequence peptide, Size link_pos, double precursor_mass, bool frag_alpha, int mincharge, int maxcharge, Size link_pos_2 = 0) const;

    // TODO
    virtual void addXLinkIonPeaks(PeakSpectrum& spectrum, DataArrays::IntegerDataArray & integer_array, DataArrays::StringDataArray & string_array, AASequence peptide, Size link_pos, double precursor_mass, bool frag_alpha, Residue::ResidueType res_type, int charge, Size link_pos_2 = 0) const;

    /// overwrite
    void updateMembers_();

    //@}

    protected:
//      /// helper to add an isotope cluster to a spectrum
//      void addIsotopeCluster_(PeakSpectrum & spectrum, const AASequence ion, const AASequence other_peptide, double cross_linker_mass, Residue::ResidueType res_type, Int charge, double intensity, String ion_type) const;

      /// helper to add a single peak to a spectrum
      void addPeak_(PeakSpectrum & spectrum, DataArrays::IntegerDataArray & integer_array, DataArrays::StringDataArray & string_array, double pos, double intensity, Residue::ResidueType res_type, Size ion_index, int charge, String ion_type) const;

      // TODO copied from normal TSG, but it is protected over there. Move it to Residue class maybe?
      /// helper for mapping residue type to letter
      char residueTypeToIonLetter_(Residue::ResidueType res_type) const;

//      /// helper to add full neutral loss ladders
//      //void addLosses_(RichPeakSpectrum & spectrum, const AASequence & ion, double intensity, Residue::ResidueType res_type, int charge, bool fragment_alpha_chain) const;
//      void addXLinkLosses_(PeakSpectrum & spectrum, const AASequence & ion, const AASequence & second_peptide, double cross_linker_mass, Residue::ResidueType res_type, int charge, double intensity, String ion_type) const;

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
      bool add_abundant_immonium_ions_;
      bool multiple_fragmentation_mode_;
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

#endif // THEORETICALSPECTRUMGENERATORXLMS_H
