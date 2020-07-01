// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
// $Maintainer: Oliver Alka $
// $Authors: Oliver Alka $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/MSSpectrum.h>

namespace OpenMS
{
  class OPENMS_DLLAPI SiriusFragmentAnnotation
      {
      public:

          /**
          @brief extractAndResolveSiriusAnnotations
          Extract SIRIUS target or decoy annotation for mapping native_id to MSSpectrum.

          @return map native_id to annotated MSSpectrum (target or decoy)

          The data is stored in a MSSpectrum, which contains a Peak1D (mz, int), a FloatDataArray (exact mass), and a StringDataArray (explanation).

          <table>
          <caption id="SiriusFragmentAnnotation"> MSSpectrum </caption>
          <tr><th> Peak1D <th> <th> FloatDataArray <th> StringDataArray
          <tr><td> mz <td> intensity <td> exact_mass <td> explanation
          <tr><td> 56.050855 <td> 20794.85 <td> 56.049476 <td> C3H5N
          </table>

          @param sirius_workspace_subdirs: Vector of paths to SIRIUS subdirectories.
          @param use_exact_mass: Option to use exact mass instead of peak mz in MSSpectrum.
          @param decoy: Use for target == false, use for decoy == true.
          */
          static std::map<String, MSSpectrum> extractAndResolveSiriusAnnotations(std::vector<String> const &sirius_workspace_subdirs, bool use_exact_mass, bool decoy);

          /**
          @brief extractSiriusFragmentAnnotationMapping  
          Extract native id (./spectrum.ms) and fragment annotation (./spectra/1_sumformula.tsv) from SIRIUS output (per compound).

          @return annotated (consensus) MSSpectrum with associated native id

          MetaValues:
            peak_mz: indicates which mass input was used in Peak1D (mass or exact_mass).
            annotated_sumformula
            annotated_adduct

          The data is stored in a MSSpectrum, which contains a Peak1D (mz, int), a FloatDataArray (exact mass), and a StringDataArray (explanation).

          <table>
          <caption id="SiriusFragmentAnnotation"> MSSpectrum </caption>
          <tr><th> Peak1D <th> <th> FloatDataArray <th> StringDataArray
          <tr><td> mz <td> intensity <td> exact_mass <td> explanation
          <tr><td> 56.050855 <td> 20794.85 <td> 56.049476 <td> C3H5N
          </table>

          @param path_to_sirius_workspace: Path to SIRIUS workspace.
          @param use_exact_mass: Option to use exact mass instead of peak mz in MSSpectrum.
          */
          static void extractSiriusFragmentAnnotationMapping(const String& path_to_sirius_workspace, MSSpectrum& msspectrum_to_fill, bool use_exact_mass = false);

          /**
          @brief extractSiriusDecoyAnnotationMapping
          Extract native id (./spectrum.ms) and fragment annotation (./decoy/1_sumformula.tsv) from SIRIUS/PASSATUTTO output (per compound).

          @return annotated decoy MSSpectrum with associated native id

          MetaValues:
          peak_mz
          annotated_sumformula
          annotated_adduct

          The data is stored in a MSSpectrum, which contains a Peak1D (mz, int), a FloatDataArray (exact mass),a StringDataArray (explanation), and a StringDataArray (ionization).

          <table>
          <caption id="SiriusFragmentAnnotation"> MSSpectrum </caption>
          <tr><th> Peak1D <th> <th> FloatDataArray <th> StringDataArray
          <tr><td> mz <td> intensity <td> exact_mass <td> explanation
          <tr><td> 56.050855 <td> 20794.85 <td> 56.049476 <td> C3H5N
          </table>

          @param path_to_sirius_workspace: Path to SIRIUS workspace.
          @param use_exact_mass: Option to use exact mass instead of peak mz in MSSpectrum.
          */
          static void extractSiriusDecoyAnnotationMapping(const String& path_to_sirius_workspace, MSSpectrum& msspectrum_to_fill);

      protected:

          /**
          @brief extractNativeIDFromSiriusMS  
          Extract native id from SIRIUS output (./spectrum.ms).

          @return String native id of current SIRIUS compound
          
          @param path_to_sirius_workspace: Path to SIRIUS workspace.
          */
          static OpenMS::String extractNativeIDFromSiriusMS_(const OpenMS::String& path_to_sirius_workspace);

          /**
          @brief extractMIDFromSiriusMS
          Extract mid from SIRIUS output (./spectrum.ms).
          Mid is the native id + an index, which is incremented based
          on the number of possible identifications (accurate mass search).

          @return String mid of current SIRIUS compound

          @param path_to_sirius_workspace: Path to SIRIUS workspace.
          */
          static OpenMS::String extractMIDFromSiriusMS_(const String& path_to_sirius_workspace);

          /**
          @brief extractAnnotationFromSiriusFile  
          Extract fragment annotation from SIRIUS  (./spectra/1_sumformula.tsv).

          @return annotated (consensus) MSSpectrum (mz, int, exact mass, fragment explanation).

          MetaValues:
            peak_mz: indicates which mass input was used in Peak1D (mz or exact_mass).
            annotated_sumformula
            annotated_adduct
          
          @param path_to_sirius_workspace: Path to SIRIUS workspace.
          @param use_exact_mass: Option to use exact mass instead of peak mz in MSSpectrum.
          */
          static void extractAnnotationFromSiriusFile_(const String& path_to_sirius_workspace, MSSpectrum& msspectrum_to_fill, bool use_exact_mass = false);

          /**
          @brief extractCompoundRankingAndFilename
          Extract compound ranking and filename   (./formula_candidates.tsv).

          @return a map with specified rank and filename (formula_adduct.tsv) (based on the annotation)

          @param path_to_sirius_workspace: Path to SIRIUS workspace.
          */
          static std::map< Size, String > extractCompoundRankingAndFilename_(const String& path_to_sirius_workspace);

          /**
          @brief extractCompoundRankingAndFilename
          Extract compound ranking and score   (./formula_candidates.tsv).

          @return a map with specified rank and score (TreeIsotope_Score) (based on the annotation)

          @param path_to_sirius_workspace: Path to SIRIUS workspace.
          */
          static std::map< Size, double > extractCompoundRankingAndScore_(const String& path_to_sirius_workspace);

          /**
          @brief extractAnnotationFromDecoyFile
          Extract decoy annotation from SIRIUS/PASSATUTTO  (./decoy/1_sumformula.tsv).

          @return annotated decoy MSSpectrum (mz, rel. int, fragment explanation, ionization).

          MetaValues:
            peak_mz
            annotated_sumformula
            annotated_adduct

          @param path_to_sirius_workspace: Path to SIRIUS workspace.
          @param use_exact_mass: Option to use exact mass instead of peak mz in MSSpectrum.
          */
          static void extractAnnotationFromDecoyFile_(const String& path_to_sirius_workspace, MSSpectrum& msspectrum_to_fill);

  };
} // namespace OpenMS
