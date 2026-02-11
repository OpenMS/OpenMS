// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Eugen Netz $
// $Authors: Eugen Netz $
// --------------------------------------------------------------------------


#pragma once

#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/METADATA/DataArrays.h>
#include <OpenMS/ANALYSIS/XLMS/OPXLDataStructs.h>


namespace OpenMS
{
  class AASequence;

  /**
      @brief Generates theoretical spectra for cross-linked peptides

      The spectra this class generates are instances of the PeakSpectrum class.
      This class generates the same peak types as SimpleTSGXLMS
      and the interface is very similar, but it is more complex and slower.
      If the parameters add_metainfo and add_charges are set to true, it will
      generate a StringDataArray for String annotations of ion types
      and an IntegerDataArray for peak charges and add them to the DataArrays
      of the produced PeakSpectrum. The spectra from this class are mainly used
      for annotation of matched experimental spectra.

  @htmlinclude OpenMS_TheoreticalSpectrumGeneratorXLMS.parameters

      @ingroup Chemistry
  */
  class OPENMS_DLLAPI TheoreticalSpectrumGeneratorXLMS :
    public DefaultParamHandler
  {
    public:

      struct LossIndex
      {
        bool has_H2O_loss = false;
        bool has_NH3_loss = false;
      };

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

        @param[in] spectrum The spectrum to which the new peaks are added. Does not have to be empty, the generated peaks will be pushed onto it.
        @param[in] peptide The peptide to fragment
        @param[in] link_pos The position of the cross-linker on the given peptide
        @param[in] frag_alpha True, if the fragmented peptide is the Alpha peptide. Used for ion-name annotation.
        @param[in] charge The maximal charge of the ions
        @param[in] link_pos_2 A second position for the linker, in case it is a loop link
       */
      virtual void getLinearIonSpectrum(PeakSpectrum & spectrum, AASequence & peptide, Size link_pos, bool frag_alpha, int charge = 1, Size link_pos_2 = 0) const;

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
          Although this function is more general, currently it is mainly used for loop-links and mono-links,
          because residues in the second, unknown peptide cannot be considered for possible neutral losses.
          The generated ion types and other additional settings are determined by the tool parameters.

        @param[in] spectrum The spectrum to which the new peaks are added. Does not have to be empty, the generated peaks will be pushed onto it.
        @param[in] peptide The peptide to fragment
        @param[in] link_pos The position of the cross-linker on the given peptide
        @param[in] precursor_mass The mass of the whole cross-link candidate or the precursor mass of the experimental MS2 spectrum.
        @param[in] frag_alpha True, if the fragmented peptide is the Alpha peptide. Used for ion-name annotation.
        @param[in] mincharge The minimal charge of the ions
        @param[in] maxcharge The maximal charge of the ions, it should be the precursor charge and is used to generate precursor ion peaks
        @param[in] link_pos_2 A second position for the linker, in case it is a loop link
       */
      virtual void getXLinkIonSpectrum(PeakSpectrum & spectrum, AASequence & peptide, Size link_pos, double precursor_mass, bool frag_alpha, int mincharge, int maxcharge, Size link_pos_2 = 0) const;

      /**
       * @brief Generates fragment ions containing the cross-linker for a pair of peptides.
       *
          B-ions are generated from the first linked position up to the end of the peptide,
          y-ions are generated from the beginning of the peptide up to the second linked position.
          This function generates neutral loss ions by considering both linked peptides.
          Only one of the peptides, decided by @p frag_alpha, is fragmented.
          This function is not suitable to generate fragments for mono-links or loop-links.
          This simplifies the function, but it has to be called twice to get all fragments of a peptide pair.
          The generated ion types and other additional settings are determined by the tool parameters.

        @param[in] spectrum The spectrum to which the new peaks are added. Does not have to be empty, the generated peaks will be pushed onto it.
        @param[in] crosslink ProteinProteinCrossLink to be fragmented
        @param[in] frag_alpha True, if the fragmented peptide is the Alpha peptide.
        @param[in] mincharge The minimal charge of the ions
        @param[in] maxcharge The maximal charge of the ions, it should be the precursor charge and is used to generate precursor ion peaks
       */
      virtual void getXLinkIonSpectrum(PeakSpectrum & spectrum, OPXLDataStructs::ProteinProteinCrossLink & crosslink, bool frag_alpha, int mincharge, int maxcharge) const;

      /// overwrite
      void updateMembers_() override;

    protected:

      /**
       * @brief Adds cross-link-less ions of a specific ion type and charge to a spectrum and adds ion name and charge annotations to the DataArrays
       * @param[in] spectrum The spectrum to which the new peaks are added. Does not have to be empty, the generated peaks will be pushed onto it.
       * @param[in] charges A DataArray collecting the charges of the added peaks
       * @param[in] ion_names A DataArray collecting the ion names of the added peaks
       * @param[in] peptide The peptide to fragment
       * @param[in] link_pos The position of the cross-linker on the given peptide
       * @param[in] frag_alpha True, if the fragmented peptide is the Alpha peptide. Used for ion-name annotation.
       * @param[in] res_type The ion type of the added peaks
       * @param[in] forward_losses vector of sets of losses generated by getForwardLosses_
       * @param[in] backward_losses vector of sets of losses generated by getBackwardLosses_
       * @param[in] charge The charge of the added peaks
       * @param[in] link_pos_2 A second position for the linker, in case it is a loop link
       */
      virtual void addLinearPeaks_(PeakSpectrum & spectrum, DataArrays::IntegerDataArray & charges, DataArrays::StringDataArray & ion_names, AASequence & peptide, Size link_pos, bool frag_alpha, Residue::ResidueType res_type, std::vector< LossIndex > & forward_losses, std::vector< LossIndex > & backward_losses, int charge = 1, Size link_pos_2 = 0) const;

      /**
       * @brief Adds a single peak to a spectrum and its charge and ion name to the given DataArrays

          The ion_type is a string in this form: "alpha|xi",
          the first word can be either "alpha" or "beta" and indicates the fragmented peptide,
          the two letters at the end are either "ci" or "xi" for linear ion or cross-linked ion.

       * @param[in] spectrum The spectrum to which the new peak is added
       * @param[in] charges A DataArray collecting the charges of the added peaks
       * @param[in] ion_names A DataArray collecting the ion names of the added peaks
       * @param[in] pos
       * @param[in] intensity
       * @param[in] res_type The ion type of the added peak
       * @param[in] ion_index The index of the ion (fragmentation position)
       * @param[in] charge The charge of the ion
       * @param[in] ion_type Another cross-linking specific ion-type
       */
      virtual void addPeak_(PeakSpectrum & spectrum, DataArrays::IntegerDataArray & charges, DataArrays::StringDataArray & ion_names, double pos, double intensity, Residue::ResidueType res_type, Size ion_index, int charge, String ion_type) const;

      /**
       * @brief Adds precursor masses including neutral losses for the given charge and adds charge and ion name to the given DataArrays

       * @param[in] spectrum The spectrum to which the peaks are added
       * @param[in] charges A DataArray collecting the charges of the added peaks
       * @param[in] ion_names A DataArray collecting the ion names of the added peaks
       * @param[in] precursor_mass The mass of the uncharged precursor
       * @param[in] charge The charge of the precursor
       */
      virtual void addPrecursorPeaks_(PeakSpectrum & spectrum, DataArrays::IntegerDataArray & charges, DataArrays::StringDataArray & ion_names, double precursor_mass, int charge) const;

      /**
       * @brief Adds losses for a linear ion

       * @param[in] spectrum The spectrum to which the new peak is added
       * @param[in] charges A DataArray collecting the charges of the added peaks
       * @param[in] ion_names A DataArray collecting the ion names of the added peaks
       * @param[in] mono_weight monoisotopic mass of the current ion
       * @param[in] res_type The ion type of the current ion
       * @param[in] frag_index The index of the ion (fragmentation position)
       * @param[in] intensity
       * @param[in] charge The charge of the ion
       * @param[in] ion_type Another cross-linking specific ion-type
       * @param[in] losses a set of LossMasses with which to modify the current ion
       */
      virtual void addLinearIonLosses_(PeakSpectrum & spectrum, DataArrays::IntegerDataArray& charges, DataArrays::StringDataArray& ion_names, double mono_weight, Residue::ResidueType res_type, Size frag_index, double intensity, int charge, String ion_type, LossIndex & losses) const;

      /**
       * @brief Adds losses for a cross-linked ion

       * @param[in] spectrum The spectrum to which the new peak is added
       * @param[in] charges A DataArray collecting the charges of the added peaks
       * @param[in] ion_names A DataArray collecting the ion names of the added peaks
       * @param[in] mono_weight monoisotopic mass of the current ion
       * @param[in] intensity
       * @param[in] charge The charge of the ion
       * @param[in] ion_type Another cross-linking specific ion-type
       * @param[in] losses a set of LossMasses with which to modify the current ion
       */
      virtual void addXLinkIonLosses_(PeakSpectrum& spectrum, DataArrays::IntegerDataArray& charges, DataArrays::StringDataArray& ion_names, double mono_weight, double intensity, int charge, String ion_type, LossIndex & losses) const;

      /**
       * @brief Adds one-residue-linked ion peaks, that are specific to XLMS

          These fragments consist of one whole peptide, the cross-linker and a part of the linked residue from the second peptide.
          The residue fragment on the linker is an internal ion from a y- and an a-fragmentation with the length of one residue.
          The function is called KLinked for now, but instead of K it is whatever the linker is attached to.

       * @param[in] spectrum The spectrum to which the peaks are added
       * @param[in] charges A DataArray collecting the charges of the added peaks
       * @param[in] ion_names A DataArray collecting the ion names of the added peaks
       * @param[in] peptide The fragmented peptide
       * @param[in] link_pos position of the linker on the fragmented peptide
       * @param[in] precursor_mass The mass of the whole cross-link candidate or the precursor mass of the experimental MS2 spectrum.
       * @param[in] frag_alpha True, if the fragmented peptide is the Alpha peptide. Used for ion-name annotation.
       * @param[in] charge The charge of the ion
       */
      virtual void addKLinkedIonPeaks_(PeakSpectrum & spectrum, DataArrays::IntegerDataArray & charges, DataArrays::StringDataArray & ion_names, AASequence & peptide, Size link_pos, double precursor_mass, bool frag_alpha, int charge) const;

      /**
       * @brief Adds cross-linked ions of a specific ion type and charge to a spectrum and adds ion name and charge annotations to the DataArrays

        This version of the function is for mono-links and loop-links.

       * @param[in] spectrum The spectrum to which the new peaks are added. Does not have to be empty, the generated peaks will be pushed onto it.
       * @param[in] charges A DataArray collecting the charges of the added peaks
       * @param[in] ion_names A DataArray collecting the ion names of the added peaks
       * @param[in] peptide The peptide to fragment
       * @param[in] link_pos The position of the cross-linker on the given peptide
       * @param[in] precursor_mass The mass of the whole cross-link candidate or the precursor mass of the experimental MS2 spectrum.
       * @param[in] frag_alpha True, if the fragmented peptide is the Alpha peptide. Used for ion-name annotation.
       * @param[in] res_type The ion type of the added peaks
       * @param[in] forward_losses  vector of sets of losses generated by getForwardLosses_
       * @param[in] backward_losses vector of sets of losses generated by getBackwardLosses_
       * @param[in] charge The charge of the added peaks
       * @param[in] link_pos_2 A second position for the linker, in case it is a loop link
       */
      virtual void addXLinkIonPeaks_(PeakSpectrum& spectrum, DataArrays::IntegerDataArray & charges, DataArrays::StringDataArray & ion_names, AASequence & peptide, Size link_pos, double precursor_mass, bool frag_alpha, Residue::ResidueType res_type, std::vector< LossIndex > & forward_losses, std::vector< LossIndex > & backward_losses, int charge, Size link_pos_2 = 0) const;

      /**
       * @brief Adds cross-linked ions of a specific ion type and charge to a spectrum and adds ion name and charge annotations to the DataArrays

        This version of the function is for cross-linked peptide pairs.

       * @param[in] spectrum The spectrum to which the new peaks are added. Does not have to be empty, the generated peaks will be pushed onto it.
       * @param[in] charges A DataArray collecting the charges of the added peaks
       * @param[in] ion_names A DataArray collecting the ion names of the added peaks
       * @param[in] crosslink The ProteinProteinCrossLink to be fragmented
       * @param[in] frag_alpha True, if the fragmented peptide is the Alpha peptide. Used for ion-name annotation.
       * @param[in] res_type The ion type of the added peaks
       * @param[in] forward_losses  vector of sets of losses generated by getForwardLosses_ for the fragmented peptide
       * @param[in] backward_losses vector of sets of losses generated by getBackwardLosses_ for the fragmented peptide
       * @param[in] losses_peptide2 set of losses for the second, not fragmented peptide, e.g. last set from getForwardLosses_ for the second peptide
       * @param[in] charge The charge of the added peaks
       */
      virtual void addXLinkIonPeaks_(PeakSpectrum & spectrum, DataArrays::IntegerDataArray & charges, DataArrays::StringDataArray & ion_names, OPXLDataStructs::ProteinProteinCrossLink & crosslink, bool frag_alpha, Residue::ResidueType res_type, std::vector< LossIndex > & forward_losses, std::vector< LossIndex > & backward_losses, LossIndex & losses_peptide2, int charge) const;

      /**
       * @brief Calculates sets of possible neutral losses for each position in the given peptide

        This function generates a vector of sets. Each set contains the possible neutral losses for a specific prefix of the peptide.

       * @param[in] peptide The peptide or ion for which to collect possible losses
       */
      std::vector< LossIndex > getForwardLosses_(AASequence & peptide) const;

      /**
       * @brief Calculates sets of possible neutral losses for each position in the given peptide

        This function generates a vector of sets. Each set contains the possible neutral losses for a specific suffix of the peptide.

       * @param[in] peptide The peptide or ion for which to collect possible losses
       */
      std::vector< LossIndex > getBackwardLosses_(AASequence & peptide) const;

      bool add_b_ions_;
      bool add_y_ions_;
      bool add_a_ions_;
      bool add_c_ions_;
      bool add_x_ions_;
      bool add_z_ions_;
      bool add_first_prefix_ion_;
      bool add_losses_;
      bool add_metainfo_;
      bool add_charges_;
      bool add_isotopes_;
      bool add_precursor_peaks_;
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
      bool add_k_linked_ions_;

      std::map< String, LossIndex > loss_db_;
      double loss_H2O_ = 0;
      double loss_NH3_ = 0;
  };
}
