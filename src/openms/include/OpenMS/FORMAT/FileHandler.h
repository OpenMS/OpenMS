// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationDescription.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/FORMAT/OPTIONS/PeakFileOptions.h>
#include <OpenMS/FORMAT/OPTIONS/FeatureFileOptions.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/KERNEL/ConsensusMap.h>


namespace OpenMS
{
  class PeakFileOptions;
  class MSSpectrum;
  class MSExperiment;
  class FeatureMap;
  class ConsensusMap;
  class TargetedExperiment;

  /**
    @brief Facilitates file handling by file type recognition.

    This class provides file type recognition from the file name and
    from the file content.

    It also offer a common interface to load MSExperiment data
    and allows querying for supported file types.

    @see FileTypes

    @ingroup FileIO
  */
  class OPENMS_DLLAPI FileHandler
  {
public:
    /**
      @brief Tries to determine the file type (by name or content)

      First tries to determine the type from the file name.
      If this fails, the type is determined from the file content.

      @param filename the name of the file to check

      @return A FileTypes::Type corresponding to the extension, or FileTypes::UNKNOWN if not determinable

      @exception Exception::FileNotFound is thrown if the file is not present
    */
    static FileTypes::Type getType(const String& filename);


    /**
      @brief Try to get the file type from the filename

      @param filename the name of the file to check

      @return A FileTypes::Type corresponding to the extension, or FileTypes::UNKNOWN if not determinable

      @exception Exception::FileNotFound is thrown if the file is not present
    */
    static FileTypes::Type getTypeByFileName(const String& filename);


    /**
       @brief Check if @p filename has the extension @p type 
               
       If the extension is not known (e.g. '.tmp') this is also allowed.
       However, if the extension is another one (neither @p type nor unknown), false is returned.
    */
    static bool hasValidExtension(const String& filename, const FileTypes::Type type);


    /**
      @brief If filename contains an extension, it will be removed (including the '.'). Special extensions, known to OpenMS, e.g. '.mzML.gz' will be recognized as well.

      E.g. 'experiment.featureXML' becomes 'experiment' and 'c:\\files\\data.mzML.gz' becomes 'c:\\files\\data'
      If the extension is unknown, the everything in the basename of the file after the last '.' is removed. E.g. 'future.newEnding' becomes 'future'
      If the filename does not contain '.', but the path (if any) does, nothing is removed, e.g. '/my.dotted.dir/filename' is returned unchanged.

      @param filename the name to strip

      @return the stripped filename
      
    */
    static String stripExtension(const String& filename);


    /**
      @brief Tries to find and remove a known file extension, and append the new one.

      Internally calls 'stripExtension()' and adds the new suffix to the result.
      E.g. 'experiment.featureXML'+ FileTypes::TRANSFORMATIONXML  becomes 'experiment.trafoXML' and 'c:\\files\\data.mzML.gz' + FileTypes::FEATUREXML becomes 'c:\\files\\data.featureXML'
      If the existing extension is unknown, the everything after the last '.' is removed, e.g. 'exp.tmp' + FileTypes::IDXML becomes 'exp.idXML'

      @param filename the original @p filename
      @param new_type the @p FileTypes::Types to use to set the new extension

      @return the updated string

    */
    static String swapExtension(const String& filename, const FileTypes::Type new_type);

    
    /**
      @brief Useful function for TOPP tools which have an 'out_type' parameter and want to know what
             output format to write.
             This function makes sure that the type derived from @p output_filename and @p requested_type are consistent, i.e.
             are either identical or one of them is UNKNOWN. Upon conflict, an error message is printed and the UNKNOWN type is returned.

      @param output_filename A full filename (with none, absolute or relative paths) whose type is 
                             determined using FileHandler::getTypeByFileName() internally
      @param requested_type A type as string, usually obtained from '-out_type', e.g. "FASTA" (case insensitive).
                            The string can be empty (yields UNKNOWN for this type)
      @return A consistent file type or UNKNOWN upon conflict
    */
    static FileTypes::Type getConsistentOutputfileType(const String& output_filename, const String& requested_type);


    /**
      @brief Determines the file type of a file by parsing the first few lines

      @exception Exception::FileNotFound is thrown if the file is not present
    */
    static FileTypes::Type getTypeByContent(const String& filename);

    /// @brief Returns if the file type is supported in this build of the library
    static bool isSupported(FileTypes::Type type);

    /// @brief Mutable access to the options for loading/storing
    PeakFileOptions& getOptions();

    /// @brief Non-mutable access to the options for loading/storing
    const PeakFileOptions& getOptions() const;

    /// @brief Mutable access to the feature file options for loading/storing
    FeatureFileOptions& getFeatOptions();

    /// @brief Non-mutable access to the feature file options for loading/storing
    const FeatureFileOptions& getFeatOptions() const;

    /// @brief set options for loading/storing
    void setOptions(const PeakFileOptions&);

    /// @brief set feature file options for loading/storing
    void setFeatOptions(const FeatureFileOptions&);

    /**
      @brief Loads a file into an MSExperiment

      @param filename The file name of the file to load.
      @param exp The experiment to load the data into.
      @param allowed_types A vector of supported filetypes. If the vector is empty, load from any type that we have a handler for. Otherwise @p getType() is called internally to check the type
      @param log Progress logging mode
      @param rewrite_source_file Set's the SourceFile name and path to the current file. Note that this looses the link to the primary MS run the file originated from.
      @param compute_hash If source files are rewritten, this flag triggers a recomputation of hash values. A SHA1 string gets stored in the checksum member of SourceFile.

      @exception Exception::FileNotFound is thrown if the file could not be opened
      @exception Exception::ParseError is thrown if an error occurs during parsing
    */
    void loadExperiment(const String& filename, PeakMap& exp, const std::vector<FileTypes::Type> allowed_types = std::vector<FileTypes::Type>(),
                        ProgressLogger::LogType log = ProgressLogger::NONE, const bool rewrite_source_file = false,
                        const bool compute_hash = false);

    /**
      @brief Stores an MSExperiment to a file

      The file type to store the data in is determined by the file name. Supported formats for storing are mzML, mzXML, mzData and DTA2D. If the file format cannot be determined from the file name, the mzML format is used.

      @param filename The name of the file to store the data in.
      @param exp The experiment to store.
      @param allowed_types A vector of supported filetypes. If empty we try to guess based on the filename. If that fails we throw UnableToCreateFile. If there is only one allowed type, check whether it agrees with the filename, and throw UnableToCreateFile if they disagree.
      @param log Progress logging mode

      @exception Exception::UnableToCreateFile is thrown if the file could not be written
    */
    void storeExperiment(const String& filename, const PeakMap& exp, const std::vector<FileTypes::Type> allowed_types = {}, ProgressLogger::LogType log = ProgressLogger::NONE);

        /**
      @brief Loads a single MSSpectrum from a file

      @param filename The file name of the file to load.
      @param spec The spectrum to load the data into.
      @param allowed_types A vector of supported filetypes. If the vector is empty, load from any type that we have a handler for. Otherwise @p getType() is called internally to check the type

      @exception Exception::FileNotFound is thrown if the file could not be opened
      @exception Exception::ParseError is thrown if an error occurs during parsing
    */
    void loadSpectrum(const String& filename, MSSpectrum& spec, const std::vector<FileTypes::Type> allowed_types = {});

        /**
      @brief Stores a single MSSpectrum to a file

      @param filename The file name of the file to store.
      @param spec The spectrum to store the data from.
      @param allowed_types A vector of supported filetypes. If the vector is empty, load from any type that we have a handler for. Otherwise @p getType() is called internally to check the type

      @exception Exception::UnableToCreateFile is thrown if the file could not be written
    */
    void storeSpectrum(const String& filename, MSSpectrum& spec, const std::vector<FileTypes::Type> allowed_types = {});

    /**
      @brief Loads a file into a FeatureMap

      @param filename the file name of the file to load.
      @param map The FeatureMap to load the data into.
      @param allowed_types A vector of supported filetypes. If the vector is empty, load from any type that we have a handler for. Otherwise @p getType() is called internally to check the type
      @param log Progress logging mode
      
      @exception Exception::FileNotFound is thrown if the file could not be opened
      @exception Exception::ParseError is thrown if an error occurs during parsing
    */
    void loadFeatures(const String& filename, FeatureMap& map, const std::vector<FileTypes::Type> allowed_types = {}, ProgressLogger::LogType log = ProgressLogger::NONE);

    /**
      @brief Store a FeatureMap

      @param filename the file name of the file to write.
      @param map The FeatureMap to store.
      @param allowed_types A vector of supported filetypes. If empty we try to guess based on the filename. If that fails we throw UnableToCreateFile. If there is only one allowed type, check whether it agrees with the filename, and throw UnableToCreateFile if they disagree.
      @param log Progress logging mode

      @exception Exception::UnableToCreateFile is thrown if the file could not be written
    */
    void storeFeatures(const String& filename, const FeatureMap& map, const std::vector<FileTypes::Type> allowed_types = {}, ProgressLogger::LogType log = ProgressLogger::NONE);

    /**
      @brief Loads a file into a ConsensusMap

      @param filename the file name of the file to load.
      @param map The ConsensusMap to load the data into.
      @param allowed_types A vector of supported filetypes. If the vector is empty, load from any type that we have a handler for. Otherwise @p getType() is called internally to check the type
      @param log Progress logging mode

      @exception Exception::FileNotFound is thrown if the file could not be opened
      @exception Exception::ParseError is thrown if an error occurs during parsing
    */
    void loadConsensusFeatures(const String& filename, ConsensusMap& map, const std::vector<FileTypes::Type> allowed_types = {}, ProgressLogger::LogType log = ProgressLogger::NONE);

    /**
      @brief Store a ConsensusFeatureMap

      @param filename the file name of the file to write.
      @param map The ConsensusMap to store.
      @param allowed_types A vector of supported filetypes. If empty we try to guess based on the filename. If that fails we throw UnableToCreateFile. If there is only one allowed type, check whether it agrees with the filename, and throw UnableToCreateFile if they disagree.
      @param log Progress logging mode

      @exception Exception::UnableToCreateFile is thrown if the file could not be written
    */
    void storeConsensusFeatures(const String& filename, const ConsensusMap& map, const std::vector<FileTypes::Type> allowed_types = {}, ProgressLogger::LogType log = ProgressLogger::NONE);

    /**
      @brief Loads an identification file into a proteinIdentifications and peptideIdentifications

      @param filename the file name of the file to load.
      @param additional_proteins The proteinIdentification vector to load the data into.
      @param additional_peptides The peptideIdentification vector to load the data into.
      @param allowed_types A vector of supported filetypes. If the vector is empty, load from any type that we have a handler for. Otherwise @p getType() is called internally to check the type
      @param log Progress logging mode      

      @exception Exception::FileNotFound is thrown if the file could not be opened
      @exception Exception::ParseError is thrown if an error occurs during parsing
    */
    void loadIdentifications(const String& filename, std::vector<ProteinIdentification>& additional_proteins, std::vector<PeptideIdentification>& additional_peptides, const std::vector<FileTypes::Type> allowed_types = {}, ProgressLogger::LogType log = ProgressLogger::NONE);

    /**
      @brief Stores proteins and peptides into an Identification File

      @param filename the file name of the file to write to.
      @param additional_proteins The proteinIdentification vector to load the data from.
      @param additional_peptides The peptideIdentification vector to load the data from.
      @param allowed_types A vector of supported filetypes. If empty we try to guess based on the filename. If that fails we throw UnableToCreateFile. If there is only one allowed type, check whether it agrees with the filename, and throw UnableToCreateFile if they disagree.
      @param log Progress logging mode

      @exception Exception::UnableToCreateFile is thrown if the file could not be written
    */
    void storeIdentifications(const String& filename, const std::vector<ProteinIdentification>& additional_proteins, const std::vector<PeptideIdentification>& additional_peptides, const std::vector<FileTypes::Type> allowed_types = {}, ProgressLogger::LogType log = ProgressLogger::NONE);

    /**
      @brief Load transitions of a spectral library

      @param filename the file name of the file to read.
      @param library The TargetedExperiment to load.
      @param allowed_types A vector of supported filetypes. If the vector is empty, load from any type that we have a handler for. Otherwise @p getType() is called internally to check the type
      @param log Progress logging mode
      
      @exception Exception::FileNotFound is thrown if the file could not be opened
      @exception Exception::ParseError is thrown if an error occurs during parsing
    */
    void loadTransitions(const String& filename, TargetedExperiment& library, const std::vector<FileTypes::Type> allowed_types = {}, ProgressLogger::LogType log = ProgressLogger::NONE);

    /**
      @brief Store transitions of a spectral library

      @param filename the file name of the file to write.
      @param library The TargetedExperiment to store.
      @param allowed_types A vector of supported filetypes. If empty we try to guess based on the filename. If that fails we throw UnableToCreateFile. If there is only one allowed type, check whether it agrees with the filename, and throw UnableToCreateFile if they disagree.
      @param log Progress logging mode

      @exception Exception::UnableToCreateFile is thrown if the file could not be written
    */
    void storeTransitions(const String& filename, const TargetedExperiment& library, const std::vector<FileTypes::Type> allowed_types = {}, ProgressLogger::LogType log = ProgressLogger::NONE);

        /**
      @brief Loads a file into Transformations

      @param filename the file name of the file to load.
      @param[out] map The Transformations to load the data into.
      @param fit_model Call fitModel() on the @p map before returning?
      @param allowed_types A vector of supported filetypes. If the vector is empty, load from any type that we have a handler for. Otherwise @p getType() is called internally to check the type

      @exception Exception::FileNotFound is thrown if the file could not be opened
      @exception Exception::ParseError is thrown if an error occurs during parsing
    */

    void loadTransformations(const String& filename, TransformationDescription& map, bool fit_model=true, const std::vector<FileTypes::Type> allowed_types = {});

    /**
      @brief Store Transformations

      @param filename the file name of the file to write.
      @param map The Transformations to store.
      @param allowed_types A vector of supported filetypes. If empty we try to guess based on the filename. If that fails we throw UnableToCreateFile. If there is only one allowed type, check whether it agrees with the filename, and throw UnableToCreateFile if they disagree.

      @exception Exception::UnableToCreateFile is thrown if the file could not be written
    */
    void storeTransformations(const String& filename, const TransformationDescription& map, const std::vector<FileTypes::Type> allowed_types = {});

    /**
      @brief Store QC info

      @brief Stores QC data in mzQC file with JSON format
      @param input_file mzML input file name
      @param filename mzQC output file name
      @param exp MSExperiment to extract QC data from, prior sortSpectra() and updateRanges() required
      @param feature_map FeatureMap from feature file (featureXML)
      @param prot_ids protein identifications from ID file (idXML)
      @param pep_ids protein identifications from ID file (idXML)
      @param consensus_map an optional consensus map to store.
      @param contact_name name of the person creating the mzQC file
      @param contact_address contact address (mail/e-mail or phone) of the person creating the mzQC file
      @param description description and comments about the mzQC file contents
      @param label unique and informative label for the run
      @param remove_duplicate_features whether to remove duplicate features only for QCML for now
      @param allowed_types A vector of supported filetypes. If empty we try to guess based on the filename. If that fails we throw UnableToCreateFile. If there is only one allowed type, check whether it agrees with the filename, and throw UnableToCreateFile if they disagree.

      @exception Exception::UnableToCreateFile is thrown if the file could not be written
    */
    void storeQC(const String& input_file,
               const String& filename,
               const MSExperiment& exp,
               const FeatureMap& feature_map,
               std::vector<ProteinIdentification>& prot_ids,
               std::vector<PeptideIdentification>& pep_ids,
               const ConsensusMap& consensus_map = ConsensusMap(),
               const String& contact_name = "",
               const String& contact_address = "",
               const String& description = "",
               const String& label = "label",
               const bool remove_duplicate_features = false,
               const std::vector<FileTypes::Type> allowed_types = {});

    /**
      @brief Computes a SHA-1 hash value for the content of the given file.

      @return The SHA-1 hash of the given file.
    */
    static String computeFileHash(const String& filename);

private:
    PeakFileOptions options_;
    FeatureFileOptions f_options_;
  };

} //namespace

