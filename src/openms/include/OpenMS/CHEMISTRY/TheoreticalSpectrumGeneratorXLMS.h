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


#ifndef OPENMS_CHEMISTRY_THEORETICALSPECTRUMGENERATORXLMS_H
#define OPENMS_CHEMISTRY_THEORETICALSPECTRUMGENERATORXLMS_H

#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/METADATA/DataArrays.h>


namespace OpenMS
{
  class AASequence;

  /**
      @brief Generates theoretical spectra for cross-linked peptides

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
      ~TheoreticalSpectrumGeneratorXLMS() override;
      //@}

      /// assignment operator
      TheoreticalSpectrumGeneratorXLMS & operator=(const TheoreticalSpectrumGeneratorXLMS & tsg);

      /**
       * @brief Generates fragment ions not containing the cross-linker for one peptide.

        B-ions are generated from the beginning of the peptide up to the first linked position,
        y-ions are generated from the second linked position up the end of the peptide.
        If link_pos_2 is 0, a mono-link or cross-link is assumed and the second position is the same as the first position.
        For a loop-link two different positions can be set and link_pos_2 must be larger than link_pos.
        The generated ion types and other additional settings are determined by the tool parameters.

       * @param spectrum The spectrum to which the new peaks are added. Does not have to be empty, the generated peaks will be pushed onto it.
       * @param peptide The peptide to fragment
       * @param link_pos The position of the cross-linker on the given peptide
       * @param frag_alpha True, if the fragmented peptide is the Alpha peptide. Used for ion-name annotation.
       * @param charge The maximal charge of the ions
       * @param link_pos_2 A second position for the linker, in case it is a loop link
       */
      virtual void getCommonIonSpectrum(PeakSpectrum & spectrum, AASequence peptide, Size link_pos, bool frag_alpha, int charge = 1, Size link_pos_2 = 0) const;

      /**
       * @brief Generates fragment ions containing the cross-linker for one peptide.
       *
          B-ions are generated from the first linked position up to the end of the peptide,
          y-ions are generated from the beginning of the peptide up to the second linked position.
          If link_pos_2 is 0, a mono-link or cross-link is assumed and the second position is the same as the first position.
          For a loop-link two different positions can be set and link_pos_2 must be larger than link_pos.
          Since in the case of a cross-link a whole second peptide is attached to the other side of the cross-link,
          a precursor mass for the two peptides and the linker is needed.
          In the case of a loop link the precursor mass is the mass of the only peptide and the linker.
          The generated ion types and other additional settings are determined by the tool parameters.

       * @param spectrum The spectrum to which the new peaks are added. Does not have to be empty, the generated peaks will be pushed onto it.
       * @param peptide The peptide to fragment
       * @param link_pos The position of the cross-linker on the given peptide
       * @param precursor_mass The mass of the whole cross-link candidate or the precursor mass of the experimental MS2 spectrum.
       * @param frag_alpha True, if the fragmented peptide is the Alpha peptide. Used for ion-name annotation.
       * @param mincharge The minimal charge of the ions
       * @param maxcharge The maximal charge of the ions
       * @param link_pos_2 A second position for the linker, in case it is a loop link
       */
      virtual void getXLinkIonSpectrum(PeakSpectrum & spectrum, AASequence peptide, Size link_pos, double precursor_mass, bool frag_alpha, int mincharge, int maxcharge, Size link_pos_2 = 0) const;

      /// overwrite
      void updateMembers_() override;

    protected:

      /**
       * @brief Adds cross-link-less ions of a specific ion type and charge to a spectrum and adds ion name and charge annotations to the DataArrays
       * @param spectrum The spectrum to which the new peaks are added. Does not have to be empty, the generated peaks will be pushed onto it.
       * @param charges A DataArray collecting the charges of the added peaks
       * @param ion_names A DataArray collecting the ion names of the added peaks
       * @param peptide The peptide to fragment
       * @param link_pos The position of the cross-linker on the given peptide
       * @param frag_alpha True, if the fragmented peptide is the Alpha peptide. Used for ion-name annotation.
       * @param res_type The ion type of the added peaks
       * @param charge The charge of the added peaks
       * @param link_pos_2 A second position for the linker, in case it is a loop link
       */
      virtual void addCommonPeaks_(PeakSpectrum & spectrum, DataArrays::IntegerDataArray & charges, DataArrays::StringDataArray & ion_names, AASequence peptide, Size link_pos, bool frag_alpha, Residue::ResidueType res_type, int charge = 1, Size link_pos_2 = 0) const;

      /**
       * @brief Adds cross-linked ions of a specific ion type and charge to a spectrum and adds ion name and charge annotations to the DataArrays
       * @param spectrum The spectrum to which the new peaks are added. Does not have to be empty, the generated peaks will be pushed onto it.
       * @param charges A DataArray collecting the charges of the added peaks
       * @param ion_names A DataArray collecting the ion names of the added peaks
       * @param peptide The peptide to fragment
       * @param link_pos The position of the cross-linker on the given peptide
       * @param precursor_mass The mass of the whole cross-link candidate or the precursor mass of the experimental MS2 spectrum.
       * @param frag_alpha True, if the fragmented peptide is the Alpha peptide. Used for ion-name annotation.
       * @param res_type The ion type of the added peaks
       * @param charge The charge of the added peaks
       * @param link_pos_2 A second position for the linker, in case it is a loop link
       */
      virtual void addXLinkIonPeaks_(PeakSpectrum& spectrum, DataArrays::IntegerDataArray & charges, DataArrays::StringDataArray & ion_names, AASequence peptide, Size link_pos, double precursor_mass, bool frag_alpha, Residue::ResidueType res_type, int charge, Size link_pos_2 = 0) const;


      /**
       * @brief Adds a single peak to a spectrum and its charge and ion name to the given DataArrays

          The ion_type is a string in this form: "alpha|xi",
          the first word can be either "alpha" or "beta" and indicates the fragmented peptide,
          the two letters at the end are either "ci" or "xi" for common ion or cross-linked ion.

       * @param spectrum The spectrum to which the new peak is added
       * @param charges A DataArray collecting the charges of the added peaks
       * @param ion_names A DataArray collecting the ion names of the added peaks
       * @param pos
       * @param intensity
       * @param res_type The ion type of the added peak
       * @param ion_index The index of the ion (fragmentation position)
       * @param charge The charge of the ion
       * @param ion_type Another cross-linking specific ion-type
       */
      void addPeak_(PeakSpectrum & spectrum, DataArrays::IntegerDataArray & charges, DataArrays::StringDataArray & ion_names, double pos, double intensity, Residue::ResidueType res_type, Size ion_index, int charge, String ion_type) const;

      // TODO copied from normal TSG, but it is protected over there. Move it to Residue class maybe?
      /// helper for mapping residue type to letter
      char residueTypeToIonLetter_(Residue::ResidueType res_type) const;

//      /// TODO helper to add full neutral loss ladders
//      //void addLosses_(RichPeakSpectrum & spectrum, const AASequence & ion, double intensity, Residue::ResidueType res_type, int charge, bool fragment_alpha_chain) const;

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
