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

  The spectra this class generates are vectors of SimplePeaks.
  This class generates the same peak types as TheoreticalSpectrumGeneratorXLMS
  and the interface is very similar, but it is simpler and faster.
  SimplePeak only contains an mz value and a charge. No intensity values
  or String annotations or other additional DataArrays are generated.

  @htmlinclude OpenMS_SimpleTSGXLMS.parameters

  @ingroup Chemistry
  */
  class OPENMS_DLLAPI SimpleTSGXLMS :
    public DefaultParamHandler
  {
    public:

      /**
        * @brief A simple struct to represent peaks with mz and charge and sort them easily
       */
      struct SimplePeak
      {
        double mz;
        int charge;

        SimplePeak(double mz, int charge)
        : mz(mz), charge(charge)
        {}

        SimplePeak()
        : mz(0.0), charge(0)
        {}
      };

      /**
        * @brief Comparator to sort SimplePeaks by mz
       */
      struct SimplePeakComparator
      {
        bool operator() (const SimplePeak& a, const SimplePeak& b)
        {
          return a.mz < b.mz;
        }
      };

      struct LossIndex
      {
        bool has_H2O_loss = false;
        bool has_NH3_loss = false;
      };

      /** @name Constructors and Destructors
      */
      //@{
      /// default constructor
      SimpleTSGXLMS();

      /// copy constructor
      SimpleTSGXLMS(const SimpleTSGXLMS & source);

      /// destructor
      ~SimpleTSGXLMS() override;
      //@}

      /// assignment operator
      SimpleTSGXLMS & operator=(const SimpleTSGXLMS & tsg);

      /**
       * @brief Generates fragment ions not containing the cross-linker for one peptide.

        B-ions are generated from the beginning of the peptide up to the first linked position,
        y-ions are generated from the second linked position up the end of the peptide.
        If link_pos_2 is 0, a mono-link or cross-link is assumed and the second position is the same as the first position.
        For a loop-link two different positions can be set and link_pos_2 must be larger than link_pos.
        The generated ion types and other additional settings are determined by the tool parameters.

        @param spectrum The spectrum to which the new peaks are added. Does not have to be empty, the generated peaks will be pushed onto it.
        @param peptide The peptide to fragment
        @param link_pos The position of the cross-linker on the given peptide
        @param charge The maximal charge of the ions
        @param link_pos_2 A second position for the linker, in case it is a loop link
       */
      virtual void getLinearIonSpectrum(std::vector< SimplePeak >& spectrum, AASequence& peptide, Size link_pos, int charge = 1, Size link_pos_2 = 0) const;

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

        @param spectrum The spectrum to which the new peaks are added. Does not have to be empty, the generated peaks will be pushed onto it.
        @param peptide The peptide to fragment
        @param link_pos The position of the cross-linker on the given peptide
        @param precursor_mass The mass of the whole cross-link candidate or the precursor mass of the experimental MS2 spectrum.
        @param min_charge The minimal charge of the ions
        @param max_charge The maximal charge of the ions, it should be the precursor charge and is used to generate precursor ion peaks
        @param link_pos_2 A second position for the linker, in case it is a loop link
       */
      virtual void getXLinkIonSpectrum(std::vector< SimplePeak >& spectrum, AASequence& peptide, Size link_pos, double precursor_mass, int min_charge, int max_charge, Size link_pos_2 = 0) const;

      /**
       * @brief Generates fragment ions containing the cross-linker for a pair of peptides.
       *
          B-ions are generated from the first linked position up to the end of the peptide,
          y-ions are generated from the beginning of the peptide up to the second linked position.
          This function generates neutral loss ions by considering both linked peptides.
          Only one of the peptides, decided by @p frag_alpha, is fragmented.
          This simplifies the function, but it has to be called twice to get all fragments of a peptide pair.
          The generated ion types and other additional settings are determined by the tool parameters.
          This function is not suitable to generate fragments for mono-links or loop-links.

        @param spectrum The spectrum to which the new peaks are added. Does not have to be empty, the generated peaks will be pushed onto it.
        @param crosslink ProteinProteinCrossLink to be fragmented
        @param frag_alpha True, if the fragmented peptide is the Alpha peptide.
        @param min_charge The minimal charge of the ions
        @param max_charge The maximal charge of the ions, it should be the precursor charge and is used to generate precursor ion peaks
       */
      virtual void getXLinkIonSpectrum(std::vector< SimplePeak >& spectrum, OPXLDataStructs::ProteinProteinCrossLink& crosslink, bool frag_alpha, int min_charge, int max_charge) const;

      /// overwrite
      void updateMembers_() override;

    protected:

      /**
       * @brief Adds cross-link-less ions of a specific ion type and charge to a spectrum
       * @param spectrum The spectrum to which the new peaks are added. Does not have to be empty, the generated peaks will be pushed onto it.
       * @param peptide The peptide to fragment
       * @param link_pos The position of the cross-linker on the given peptide
       * @param res_type The ion type of the added peaks
       * @param forward_losses vector of LossIndex generated by getForwardLosses_
       * @param backward_losses vector of LossIndex generated by getBackwardLosses_
       * @param charge The charge of the added peaks
       * @param link_pos_2 A second position for the linker, in case it is a loop link
       */
      virtual void addLinearPeaks_(std::vector< SimplePeak >& spectrum, AASequence& peptide, Size link_pos, Residue::ResidueType res_type, std::vector< LossIndex >& forward_losses, std::vector< LossIndex >& backward_losses, int charge = 1, Size link_pos_2 = 0) const;

      /**
       * @brief Adds precursor masses including neutral losses for the given charge

       * @param spectrum The spectrum to which the peaks are added
       * @param precursor_mass The mass of the uncharged precursor
       * @param charge The charge of the precursor
       */
      virtual void addPrecursorPeaks_(std::vector< SimplePeak >& spectrum, double precursor_mass, int charge) const;

      /**
       * @brief Adds neutral losses for an ion to a spectrum

       * @param spectrum The spectrum to which the new peak is added
       * @param mono_weight monoisotopic mass of the current ion
       * @param charge The charge of the ion
       * @param losses a LossIndex with which to modify the current ion
       */
      virtual void addLosses_(std::vector< SimplePeak >& spectrum, double mono_weight, int charge, LossIndex & losses) const;

      /**
       * @brief Adds one-residue-linked ion peaks, that are specific to XLMS

          These fragments consist of one whole peptide, the cross-linker and a part of the linked residue from the second peptide.
          The residue fragment on the linker is an internal ion from a y- and an a-fragmentation with the length of one residue.
          The function is called KLinked for now, but instead of K it is whatever residue the linker is bound to.

       * @param spectrum The spectrum to which the peaks are added
       * @param peptide The fragmented peptide
       * @param link_pos position of the linker on the fragmented peptide
       * @param precursor_mass The mass of the whole cross-link candidate or the precursor mass of the experimental MS2 spectrum.
       * @param charge The charge of the ion
       */
      virtual void addKLinkedIonPeaks_(std::vector< SimplePeak >& spectrum, AASequence & peptide, Size link_pos, double precursor_mass, int charge) const;

      /**
       * @brief Adds cross-linked ions of a specific ion type and charge to a spectrum

        This version of the function is for mono-links and loop-links.

       * @param spectrum The spectrum to which the new peaks are added. Does not have to be empty, the generated peaks will be pushed onto it.
       * @param peptide The peptide to fragment
       * @param link_pos The position of the cross-linker on the given peptide
       * @param precursor_mass The mass of the whole cross-link candidate or the precursor mass of the experimental MS2 spectrum.
       * @param res_type The ion type of the added peaks
       * @param forward_losses  vector of LossIndex generated by getForwardLosses_
       * @param backward_losses vector of LossIndex generated by getBackwardLosses_
       * @param charge The charge of the added peaks
       * @param link_pos_2 A second position for the linker, in case it is a loop link
       */
      virtual void addXLinkIonPeaks_(std::vector<SimplePeak>& spectrum, AASequence& peptide, Size link_pos, double precursor_mass, Residue::ResidueType res_type, std::vector< LossIndex > & forward_losses, std::vector< LossIndex > & backward_losses, int charge, Size link_pos_2 = 0) const;

      /**
       * @brief Adds cross-linked ions of a specific ion type and charge to a spectrum

        This version of the function is for cross-linked peptide pairs.

       * @param spectrum The spectrum to which the new peaks are added. Does not have to be empty, the generated peaks will be pushed onto it.
       * @param crosslink The ProteinProteinCrossLink to be fragmented
       * @param frag_alpha True, if the fragmented peptide is the Alpha peptide. Used for ion-name annotation.
       * @param res_type The ion type of the added peaks
       * @param forward_losses  vector of LossIndex generated by getForwardLosses_ for the fragmented peptide
       * @param backward_losses vector of LossIndex generated by getBackwardLosses_ for the fragmented peptide
       * @param losses_peptide2 set of losses for the second, not fragmented peptide, e.g. last set from getForwardLosses_ for the second peptide
       * @param charge The charge of the added peaks
       */
      virtual void addXLinkIonPeaks_(std::vector< SimplePeak >& spectrum, OPXLDataStructs::ProteinProteinCrossLink & crosslink, bool frag_alpha, Residue::ResidueType res_type, std::vector< LossIndex > & forward_losses, std::vector< LossIndex > & backward_losses, LossIndex & losses_peptide2, int charge) const;

      /**
       * @brief Calculates sets of possible neutral losses for each position in the given peptide

        This function generates a vector of sets. Each set contains the possible neutral losses for a specific prefix of the peptide.

       * @param peptide The peptide or ion for which to collect possible losses
       */
      std::vector< LossIndex > getForwardLosses_(AASequence & peptide) const;

      /**
       * @brief Calculates sets of possible neutral losses for each position in the given peptide

        This function generates a vector of sets. Each set contains the possible neutral losses for a specific suffix of the peptide.

       * @param peptide The peptide or ion for which to collect possible losses
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
      bool add_charges_;
      bool add_isotopes_;
      bool add_precursor_peaks_;
      bool add_abundant_immonium_ions_;
      Int max_isotope_;
      double pre_int_;
      double pre_int_H2O_;
      double pre_int_NH3_;
      bool add_k_linked_ions_;

      std::map< String, LossIndex > loss_db_;
      double loss_H2O_ = 0;
      double loss_NH3_ = 0;
  };
}
