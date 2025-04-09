// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Eugen Netz $
// $Authors: Lukas Zimmermann, Eugen Netz $
// --------------------------------------------------------------------------
#pragma once

#include <OpenMS/FORMAT/XMLFile.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/ANALYSIS/XLMS/OPXLDataStructs.h>

namespace OpenMS
{

  /**
    @brief Used to load and store xQuest result files

    These files are used to store results derived from chemical cross-linking
    coupled to MS experiments.

    A documented schema for this format can be found at https://github.com/OpenMS/OpenMS/tree/develop/share/OpenMS/SCHEMAS

    @ingroup FileIO
  */
  class OPENMS_DLLAPI XQuestResultXMLFile :
    public Internal::XMLFile,
    public ProgressLogger
  {
public:
    XQuestResultXMLFile();
    ~XQuestResultXMLFile() override;

    /**
      @brief Load the content of the xquest.xml file into the provided data structures.
      @param filename Filename of the file which is to be loaded.
      @param pep_ids Where the spectra with identifications of the input file will be loaded to.
      @param prot_ids Where the protein identification of the input file will be loaded to.
     */
    void load(const String & filename,
              std::vector< PeptideIdentification > & pep_ids,
              std::vector< ProteinIdentification > & prot_ids
            );

    /**
        @brief Stores the identifications in a xQuest XML file.

        @exception Exception::UnableToCreateFile is thrown if the file could not be created
    */
    void store(const String& filename,
               const std::vector<ProteinIdentification>& poid,
               const std::vector<PeptideIdentification>& peid) const;

    /**
     * @brief Returns the total number of hits in the file
     * @return total number of hits in the file
     */
    int getNumberOfHits() const;

    /**
     * @brief Returns minimum score among the hits in the file.
     * @return Minimum score among the hits in the file.
     */
    double getMinScore() const;

    /**
     * @brief Returns maximum score among the hits in the file.
     * @return Maximum score among the hits in the file.
     */
    double getMaxScore() const;

     /**
      * @brief Writes spec.xml output containing matching peaks between heavy and light spectra after comparing and filtering
       @param out_file Path and filename for the output file
       @param base_name The base_name should be the name of the input spectra file without the file ending. Used as part of an identifier string for the spectra.
       @param preprocessed_pair_spectra The preprocessed spectra after comparing and filtering
       @param spectrum_pairs Indices of spectrum pairs in the input map
       @param all_top_csms CrossLinkSpectrumMatches, from which the IDs were generated. Only spectra with matches are written out.
       @param spectra The spectra, that were searched as a PeakMap. The indices in spectrum_pairs correspond to spectra in this map.
       @param test_mode Skip base64 encoding in test mode
      */
    static void writeXQuestXMLSpec(const String& out_file, const String& base_name,
                                   const OPXLDataStructs::PreprocessedPairSpectra& preprocessed_pair_spectra,
                                   const std::vector< std::pair<Size, Size> >& spectrum_pairs,
                                   const std::vector< std::vector< OPXLDataStructs::CrossLinkSpectrumMatch > >& all_top_csms,
                                   const PeakMap& spectra, const bool& test_mode = false);

     /**
       @brief Writes spec.xml output containing spectra for visualization. This version of the function is meant to be used for label-free linkers.
       @param out_file Path and filename for the output file
       @param base_name The base_name should be the name of the input spectra file without the file ending. Used as part of an identifier string for the spectra.
       @param all_top_csms CrossLinkSpectrumMatches, from which the IDs were generated. Only spectra with matches are written out.
       @param spectra The spectra, that were searched as a PeakMap.
       @param test_mode Skip base64 encoding in test mode
      */
    static void writeXQuestXMLSpec(const String& out_file, const String& base_name,
                                   const std::vector< std::vector< OPXLDataStructs::CrossLinkSpectrumMatch > >& all_top_csms,
                                   const PeakMap& spectra, const bool& test_mode = false);



private:

     /**
      * @brief Transforms a PeakSpectrum into a base64 encoded string, which is the format used in spec.xml for xQuest.
      * @param spec The spectrum
      * @param header A header for the spectrum, build using the base_name parameter for writeXQuestXMLSpec and the index of the spectrum.
      * @param test_mode Skip base64 encoding in test mode
      */
      static String getxQuestBase64EncodedSpectrum_(const PeakSpectrum& spec, const String& header, const bool& test_mode = false);

     /**
      * @brief A helper function, that takes one string containing one line and wraps it into several lines of a given width
      * @param input The string as one line
      * @param width The preferred line width
      * @param output String in which the output is written
      */
      static void wrap_(const String& input, Size width, String& output);

    int n_hits_; ///< Total number of hits within the result file
    double min_score_; ///< Minimum score encountered in file
    double max_score_; ///< Maximum score encountered in file
  };
} // namespace OpenMS
