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

#ifndef OPENMS_ANALYSIS_XLMS_XQUESTXML_H
#define OPENMS_ANALYSIS_XLMS_XQUESTXML_H

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/ANALYSIS/XLMS/OPXLDataStructs.h>

namespace OpenMS
{

  /**
   * @brief The XQuestXML class contains functions to write xquest.xml and spec.xml files compatible with xQuest

      The xquest.xml file format is an XL-MS ID format used by the xQuest algorithm (O. Rinner et al., 2008, "Identification of cross-linked peptides from large sequence databases")
      and compatible with the xQuest results manager, as well as xProphet and some other post processing tools of the xQuest pipeline.
      The spec.xml format is mostly used to visualize spectra in the xQuest results manager.
   */

  class OPENMS_DLLAPI XQuestXML
  {
    public:

     /**
      * @brief Writes xquest.xml output files containing XL-MS identification results
      * @param Path and filename for the output file
      * @param The base_name should be the name of the input spectra file without the file ending. Used as part of an identifier string for the spectra.
      * @param The PeptideIdentifications from a XL-MS search
      * @param CrossLinkSpectrumMatches, from which the IDs were generated.
      * @param The spectra, that were searched as a PeakMap
      * @param The unit of the precursor mass tolerance ("Da" or "ppm")
      * @param The unit of the fragment mass tolerance ("Da" or "ppm")
      * @param The precursor mass tolerance
      * @param The fragment mass tolerance for common ions
      * @param The fragment mass tolerance for cross-linked ions
      * @param The name of the cross-link reagent
      * @param The mass of the cross-link reagent when linking two peptides
      * @param List of possible masses for the cross-link, if it is attached to a peptide on one side
      * @param Path and filename of the target protein database
      * @param Path and filename of the decoy protein database
      * @param Residues, that the first side of the cross-linker can react with
      * @param Residues, that the second side of the cross-linker can react with
      * @param The difference of mass between the two types of cross-linkers, if labeled linkers were used. Should be 0 for label-free linkers.
      * @param Name of the Enzyme used for digestion
      * @param Number of allowed missed cleavages
      */
      static void writeXQuestXML(String out_file, String base_name, const std::vector< PeptideIdentification >& peptide_ids, const std::vector< std::vector< OPXLDataStructs::CrossLinkSpectrumMatch > >& all_top_csms, const PeakMap& spectra,
                                                  String precursor_mass_tolerance_unit, String fragment_mass_tolerance_unit, double precursor_mass_tolerance, double fragment_mass_tolerance, double fragment_mass_tolerance_xlinks, String cross_link_name,
                                                  double cross_link_mass_light, DoubleList cross_link_mass_mono_link, String in_fasta, String in_decoy_fasta, StringList cross_link_residue1, StringList cross_link_residue2, double cross_link_mass_iso_shift, String enzyme_name, Size missed_cleavages);

     /**
      * @brief Writes spec.xml output containing matching peaks between heavy and light spectra after comparing and filtering
      * @param Path and filename for the output file
      * @param The base_name should be the name of the input spectra file without the file ending. Used as part of an identifier string for the spectra.
      * @param The preprocessed spectra after comparing and filtering
      * @param Indices of spectrum pairs in the input map
      * @param CrossLinkSpectrumMatches, from which the IDs were generated. Only spectra with matches are written out.
      * @param The spectra, that were searched as a PeakMap. The indices in spectrum_pairs correspond to spectra in this map.
      */
      static void writeXQuestXMLSpec(String out_file, String base_name, const OPXLDataStructs::PreprocessedPairSpectra& preprocessed_pair_spectra, const std::vector< std::pair<Size, Size> >& spectrum_pairs, const std::vector< std::vector< OPXLDataStructs::CrossLinkSpectrumMatch > >& all_top_csms, const PeakMap& spectra);

     /**
      * @brief Writes spec.xml output containing spectra for visualization. This version of the function is meant to be used for label-free linkers.
      * @param Path and filename for the output file
      * @param The base_name should be the name of the input spectra file without the file ending. Used as part of an identifier string for the spectra.
      * @param CrossLinkSpectrumMatches, from which the IDs were generated. Only spectra with matches are written out.
      * @param The spectra, that were searched as a PeakMap.
      */
      static void writeXQuestXMLSpec(String out_file, String base_name, const std::vector< std::vector< OPXLDataStructs::CrossLinkSpectrumMatch > >& all_top_csms, const PeakMap& spectra);

    protected:
     /**
      * @brief Transforms a PeakSpectrum into a base 64 encoded string, which is the format used in spec.xml for xQuest.
      * @param The spectrum
      * @param A header for the spectrum, build using the base_name parameter for writeXQuestXMLSpec and the index of the spectrum.
      */
      static String getxQuestBase64EncodedSpectrum_(const PeakSpectrum& spec, String header);

     /**
      * @brief A helper function, that takes one string containing one line and wraps it into several lines of a given width
      * @param The string as one line
      * @param The preferred line width
      * @param String in which the output is written
      */
      static void wrap_(const String& input, Size width, String & output);
  };
}

#endif // XQUESTXML_H
