// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Oliver Alka $
// $Authors: Oliver Alka $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/FORMAT/CsvFile.h>
#include <OpenMS/KERNEL/MSSpectrum.h>

namespace OpenMS
{
  class OPENMS_DLLAPI SiriusFragmentAnnotation
      {
      public:

          /**
          @brief SiriusTargetDecoySpectra holds the target and/or decoy information for one entry (subdirectory from SIRIUS)
          */
          class SiriusTargetDecoySpectra
          {
          public:

            MSSpectrum target;
            MSSpectrum decoy;

            SiriusTargetDecoySpectra() = default;
            SiriusTargetDecoySpectra(MSSpectrum target_spectrum, MSSpectrum decoy_spectrum) : target(std::move(target_spectrum)), decoy(std::move(decoy_spectrum)) {}
          };

          /**
          @brief extractAndResolveSiriusAnnotations
          Extract and resolves SIRIUS target and/or decoy annotation for mapping native_id to MSSpectrum.

          @return map native_id to annotated MSSpectrum (target or decoy)

          If there are multiple identifications for a feature with the same MS2 spectras (concatenated native ids)
          the identification with the higher SIRIUS score is chosen (currently based on the explained peak intensities).

          @param sirius_workspace_subdirs Vector of paths to SIRIUS subdirectories.
          @param score_threshold Only use spectra over a certain score threshold (0-1)
          @param use_exact_mass Option to use exact mass instead of peak mz in MSSpectrum.
          @param decoy_generation Extract decoy spectra from SIRIUS subdirectories. 
          */
          static std::vector<SiriusTargetDecoySpectra> extractAndResolveSiriusAnnotations(const std::vector<String>& sirius_workspace_subdirs, double score_threshold,
                                                                                          bool use_exact_mass, bool decoy_generation);

          /**
          @brief extractSiriusFragmentAnnotationMapping  
          Extract concatenated native ids and concatenated m_ids (unique identifier) from (./spectrum.ms) and annotations from spectra/decoy subfolder


          If @p decoy is true, uses fragment annotation (./spectra/1_sumformula.tsv) from SIRIUS output (per compound)
          else uses fragment annotation (./decoy/1_sumformula.tsv) from SIRIUS/PASSATUTTO output (per compound).

          @return annotated decoy MSSpectrum with associated native id

          MetaValues:
          peak_mz
          annotated_sumformula
          annotated_adduct

          The data is stored in a MSSpectrum, which contains a Peak1D (mz or exact mass [depending on @p use_exact_mass], int),
          a FloatDataArray for targets only (exact mass or mz [depending on @p use_exact_mass]),
          a StringDataArray (explanation), and a StringDataArray (ionization).

          <table>
          <caption> MSSpectrum </caption>
          <tr><th> Peak1D <th> <th> [FloatDataArray] <th> StringDataArray <th> StringDataArray
          <tr><td> mz <td> intensity <td> [exact_mass] <td> explanation <td> ionization
          <tr><td> 56.050855 <td> 20794.85 <td> [56.049476] <td> C3H5N <td> [M + H]+
          </table>

          @param path_to_sirius_workspace Path to SIRIUS workspace.
          @param max_rank Up to which rank to extract annotations maximally. Auto-stops at last candidate.
          @param decoy Extract annotations for decoys? Or else targets. Run twice if you want both
          @param use_exact_mass Option to use exact mass instead of peak mz in MSSpectrum.
          */
          static std::vector<MSSpectrum> extractAnnotationsFromSiriusFile(const String& path_to_sirius_workspace, Size max_rank = 1, bool decoy = false, bool use_exact_mass = false);

          /**
          @brief Extract columnname and index based in SIRIUS entries
          */
          static std::map< std::string, Size > extract_columnname_to_columnindex(const CsvFile& csvfile);

      protected:
          /**
          @brief extractConcatNativeIDsFromSiriusMS
          Extract concatenated native id from SIRIUS output (./spectrum.ms) and concatenates them.

          @return String native id of current SIRIUS compound
          
          @param path_to_sirius_workspace Path to SIRIUS workspace.
          */
          static OpenMS::String extractConcatNativeIDsFromSiriusMS_(const OpenMS::String& path_to_sirius_workspace);

          /**
          @brief extractConcatMIDsFromSiriusMS
          Extract m_ids from SIRIUS output (./spectrum.ms) and concatenates them.
          M_id is the native id + an index, which is incremented based
          on the number of possible identifications (accurate mass search).

          @return String m_id of current SIRIUS compound

          @param path_to_sirius_workspace Path to SIRIUS workspace.
          */
          static OpenMS::String extractConcatMIDsFromSiriusMS_(const String& path_to_sirius_workspace);

          /**
          @brief extractConcatMIDsFromSiriusMS
          Extract fid (i.e. original OpenMS feature ID) from SIRIUS output (./spectrum.ms).

          @return String fid of current SIRIUS workspace

          @param path_to_sirius_workspace Path to SIRIUS workspace.
          */
          static OpenMS::String extractFeatureIDFromSiriusMS_(const String& path_to_sirius_workspace);


          /**
          @brief extractCompoundRankingAndFilename
          Extract compound ranking and filename   (./formula_candidates.tsv).

          @return a map with specified rank and filename (formula_adduct.tsv) (based on the annotation)

          @param path_to_sirius_workspace Path to SIRIUS workspace.
          */
          static std::map< Size, String > extractCompoundRankingAndFilename_(const String& path_to_sirius_workspace);

          /**
          @brief extractCompoundRankingAndFilename
          Extract compound ranking and score   (./formula_candidates.tsv).

          @return a map with specified rank and score (explainedIntensity) (based on the annotation)

          @param path_to_sirius_workspace Path to SIRIUS workspace.
          */
          static std::map< Size, double > extractCompoundRankingAndScore_(const String& path_to_sirius_workspace);

  };
} // namespace OpenMS
