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

// Forward declarations instead of full includes — all these types appear only
// as references/pointers in method signatures, so full definitions are not needed.
namespace OpenMS
{
  class PeakFileOptions;
  class MSSpectrum;
  class MSExperiment;
  class FeatureMap;
  class ConsensusMap;
  class TargetedExperiment;
  class ProteinIdentification;
  class PeptideIdentificationList;

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

      @param[in] filename the name of the file to check

      @return A FileTypes::Type corresponding to the extension, or FileTypes::UNKNOWN if not determinable

      @exception Exception::FileNotFound is thrown if the file is not present
    */
    static FileTypes::Type getType(const String& filename);


    /**
      @brief Try to get the file type from the filename

      @param[in] filename the name of the file to check

      @return A FileTypes::Type corresponding to the extension, or FileTypes::UNKNOWN if not determinable

      @exception Exception::FileNotFound is thrown if the file is not present
    */
    static FileTypes::Type getTypeByFileName(const String& filename);


    /**
       @brief Check if @p filename has the extension @p type

       If the extension is not known (e.g. '.tmp') this is also allowed.
       However, if the extension is another one (neither @p type nor unknown), false is returned.

       @param[in] filename The filename to check
       @param[in] type The expected file type
    */
    static bool hasValidExtension(const String& filename, const FileTypes::Type type);


    /**
      @brief If filename contains an extension, it will be removed (including the '.'). Special extensions, known to OpenMS, e.g. '.mzML.gz' will be recognized as well.

      E.g. 'experiment.featureXML' becomes 'experiment' and 'c:\\files\\data.mzML.gz' becomes 'c:\\files\\data'
      If the extension is unknown, the everything in the basename of the file after the last '.' is removed. E.g. 'future.newEnding' becomes 'future'
      If the filename does not contain '.', but the path (if any) does, nothing is removed, e.g. '/my.dotted.dir/filename' is returned unchanged.

      @param[in] filename the name to strip

      @return the stripped filename

    */
    static String stripExtension(const String& filename);


    /**
      @brief Tries to find and remove a known file extension, and append the new one.

      Internally calls 'stripExtension()' and adds the new suffix to the result.
      E.g. 'experiment.featureXML'+ FileTypes::TRANSFORMATIONXML  becomes 'experiment.trafoXML' and 'c:\\files\\data.mzML.gz' + FileTypes::FEATUREXML becomes 'c:\\files\\data.featureXML'
      If the existing extension is unknown, the everything after the last '.' is removed, e.g. 'exp.tmp' + FileTypes::IDXML becomes 'exp.idXML'

      @param[in] filename the original @p filename
      @param[in] new_type the @p FileTypes::Types to use to set the new extension

      @return the updated string

    */
    static String swapExtension(const String& filename, const FileTypes::Type new_type);

    
    /**
      @brief Useful function for TOPP tools which have an 'out_type' parameter and want to know what
             output format to write.
             This function makes sure that the type derived from @p output_filename and @p requested_type are consistent, i.e.
             are either identical or one of them is UNKNOWN. Upon conflict, an error message is printed and the UNKNOWN type is returned.

      @param[in] output_filename A full filename (with none, absolute or relative paths) whose type is
                             determined using FileHandler::getTypeByFileName() internally
      @param[in] requested_type A type as string, usually obtained from '-out_type', e.g. "FASTA" (case insensitive).
                            The string can be empty (yields UNKNOWN for this type)
      @return A consistent file type or UNKNOWN upon conflict
    */
    static FileTypes::Type getConsistentOutputfileType(const String& output_filename, const String& requested_type);


    /**
      @brief Determines the file type of a file by parsing the first few lines

      @param[in] filename The file to check

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

      @param[out] filename The file name of the file to load.
      @param[in] exp The experiment to load the data into.
      @param[in] allowed_types A vector of supported filetypes. If the vector is empty, load from any type that we have a handler for. Otherwise @p getType() is called internally to check the type
      @param[in] log Progress logging mode
      @param[in] rewrite_source_file Set's the SourceFile name and path to the current file. Note that this looses the link to the primary MS run the file originated from.
      @param[out] compute_hash If source files are rewritten, this flag triggers a recomputation of hash values. A SHA1 string gets stored in the checksum member of SourceFile.

      @exception Exception::FileNotFound is thrown if the file could not be opened
      @exception Exception::ParseError is thrown if an error occurs during parsing
    */
    void loadExperiment(const String& filename, PeakMap& exp, const std::vector<FileTypes::Type> allowed_types = std::vector<FileTypes::Type>(),
                        ProgressLogger::LogType log = ProgressLogger::NONE, const bool rewrite_source_file = false,
                        const bool compute_hash = false);

    /**
      @brief Stores an MSExperiment to a file

      The file type to store the data in is determined by the file name. Supported formats for storing are mzML, mzXML, mzData and DTA2D. If the file format cannot be determined from the file name, the mzML format is used.

      @param[in] filename The name of the file to store the data in.
      @param[out] exp The experiment to store.
      @param[in] allowed_types A vector of supported filetypes. If empty we try to guess based on the filename. If that fails we throw UnableToCreateFile. If there is only one allowed type, check whether it agrees with the filename, and throw UnableToCreateFile if they disagree.
      @param[in] log Progress logging mode

      @exception Exception::UnableToCreateFile is thrown if the file could not be written
    */
    void storeExperiment(const String& filename, const PeakMap& exp, const std::vector<FileTypes::Type> allowed_types = {}, ProgressLogger::LogType log = ProgressLogger::NONE);

        /**
      @brief Loads a single MSSpectrum from a file

      @param[in] filename The file name of the file to load.
      @param[out] spec The spectrum to load the data into.
      @param[in] allowed_types A vector of supported filetypes. If the vector is empty, load from any type that we have a handler for. Otherwise @p getType() is called internally to check the type

      @exception Exception::FileNotFound is thrown if the file could not be opened
      @exception Exception::ParseError is thrown if an error occurs during parsing
    */
    void loadSpectrum(const String& filename, MSSpectrum& spec, const std::vector<FileTypes::Type> allowed_types = {});

        /**
      @brief Stores a single MSSpectrum to a file

      @param[in] filename The file name of the file to store.
      @param[in] spec The spectrum to store the data from.
      @param[in] allowed_types A vector of supported filetypes. If the vector is empty, load from any type that we have a handler for. Otherwise @p getType() is called internally to check the type

      @exception Exception::UnableToCreateFile is thrown if the file could not be written
    */
    void storeSpectrum(const String& filename, MSSpectrum& spec, const std::vector<FileTypes::Type> allowed_types = {});

    /**
      @brief Loads a file into a FeatureMap

      @param[in] filename the file name of the file to load.
      @param[out] map The FeatureMap to load the data into.
      @param[in] allowed_types A vector of supported filetypes. If the vector is empty, load from any type that we have a handler for. Otherwise @p getType() is called internally to check the type
      @param[in] log Progress logging mode
      
      @exception Exception::FileNotFound is thrown if the file could not be opened
      @exception Exception::ParseError is thrown if an error occurs during parsing
    */
    void loadFeatures(const String& filename, FeatureMap& map, const std::vector<FileTypes::Type> allowed_types = {}, ProgressLogger::LogType log = ProgressLogger::NONE);

    /**
      @brief Store a FeatureMap

      @param[in] filename the file name of the file to write.
      @param[out] map The FeatureMap to store.
      @param[in] allowed_types Restricts the set of file formats this call may write.
        - Empty: infer the format from the filename extension; if that fails, throw InvalidFileType.
        - Size 1: if the filename's type is UNKNOWN (e.g. TOPPAS intermediate ".unknown" files),
          fall back to this single type. Otherwise the filename type must equal it, or InvalidFileType
          is thrown.
        - Size >1: the filename type must be one of these; UNKNOWN filenames CANNOT be disambiguated
          and will throw InvalidFileType. TOPP tools writing through TOPPAS pipelines must therefore
          pass a single concrete type here (typically derived from the input file type or an
          out_type parameter), not the full list of formats they support.
      @param[in] log Progress logging mode

      @exception Exception::UnableToCreateFile is thrown if the file could not be written
      @exception Exception::InvalidFileType is thrown if the file type cannot be determined or is not allowed
    */
    void storeFeatures(const String& filename, const FeatureMap& map, const std::vector<FileTypes::Type> allowed_types = {}, ProgressLogger::LogType log = ProgressLogger::NONE);

    /**
      @brief Loads a file into a ConsensusMap

      @param[in] filename the file name of the file to load.
      @param[in] map The ConsensusMap to load the data into.
      @param[in] allowed_types A vector of supported filetypes. If the vector is empty, load from any type that we have a handler for. Otherwise @p getType() is called internally to check the type
      @param[in] log Progress logging mode

      @exception Exception::FileNotFound is thrown if the file could not be opened
      @exception Exception::ParseError is thrown if an error occurs during parsing
    */
    void loadConsensusFeatures(const String& filename, ConsensusMap& map, const std::vector<FileTypes::Type> allowed_types = {}, ProgressLogger::LogType log = ProgressLogger::NONE);

    /**
      @brief Store a ConsensusFeatureMap

      @param[in] filename the file name of the file to write.
      @param[out] map The ConsensusMap to store.
      @param[in] allowed_types Restricts the set of file formats this call may write.
        - Empty: infer the format from the filename extension; if that fails, throw InvalidFileType.
        - Size 1: if the filename's type is UNKNOWN (e.g. TOPPAS intermediate ".unknown" files),
          fall back to this single type. Otherwise the filename type must equal it, or InvalidFileType
          is thrown.
        - Size >1: the filename type must be one of these; UNKNOWN filenames CANNOT be disambiguated
          and will throw InvalidFileType. TOPP tools writing through TOPPAS pipelines must therefore
          pass a single concrete type here (typically derived from the input file type or an
          out_type parameter), not the full list of formats they support.
      @param[in] log Progress logging mode

      @exception Exception::UnableToCreateFile is thrown if the file could not be written
      @exception Exception::InvalidFileType is thrown if the file type cannot be determined or is not allowed
    */
    void storeConsensusFeatures(const String& filename, const ConsensusMap& map, const std::vector<FileTypes::Type> allowed_types = {}, ProgressLogger::LogType log = ProgressLogger::NONE);

    /**
      @brief Loads an identification file into a proteinIdentifications and peptideIdentifications

      @param[in] filename the file name of the file to load.
      @param[in] additional_proteins The proteinIdentification vector to load the data into.
      @param[in] additional_peptides The peptideIdentification vector to load the data into.
      @param[in] allowed_types A vector of supported filetypes. If the vector is empty, load from any type that we have a handler for. Otherwise @p getType() is called internally to check the type
      @param[in] log Progress logging mode      

      @exception Exception::FileNotFound is thrown if the file could not be opened
      @exception Exception::ParseError is thrown if an error occurs during parsing
    */
    void loadIdentifications(const String& filename, std::vector<ProteinIdentification>& additional_proteins, PeptideIdentificationList& additional_peptides, const std::vector<FileTypes::Type> allowed_types = {}, ProgressLogger::LogType log = ProgressLogger::NONE);

    /**
      @brief Stores proteins and peptides into an Identification File

      @param[in] filename the file name of the file to write to.
      @param[in] additional_proteins The proteinIdentification vector to load the data from.
      @param[in] additional_peptides The peptideIdentification vector to load the data from.
      @param[in] allowed_types Restricts the set of file formats this call may write.
        - Empty: infer the format from the filename extension; if that fails, throw InvalidFileType.
        - Size 1: if the filename's type is UNKNOWN (e.g. TOPPAS intermediate ".unknown" files),
          fall back to this single type. Otherwise the filename type must equal it, or InvalidFileType
          is thrown.
        - Size >1: the filename type must be one of these; UNKNOWN filenames CANNOT be disambiguated
          and will throw InvalidFileType. TOPP tools writing through TOPPAS pipelines must therefore
          pass a single concrete type here (typically derived from the input file type or an
          out_type parameter), not the full list of formats they support.
      @param[in] log Progress logging mode

      @exception Exception::UnableToCreateFile is thrown if the file could not be written
      @exception Exception::InvalidFileType is thrown if the file type cannot be determined or is not allowed
    */
    void storeIdentifications(const String& filename, const std::vector<ProteinIdentification>& additional_proteins, const PeptideIdentificationList& additional_peptides, const std::vector<FileTypes::Type> allowed_types = {}, ProgressLogger::LogType log = ProgressLogger::NONE);

    /**
      @brief Load transitions of a spectral library

      @param[in] filename the file name of the file to read.
      @param[out] library The TargetedExperiment to load.
      @param[in] allowed_types A vector of supported filetypes. If the vector is empty, load from any type that we have a handler for. Otherwise @p getType() is called internally to check the type
      @param[in] log Progress logging mode
      
      @exception Exception::FileNotFound is thrown if the file could not be opened
      @exception Exception::ParseError is thrown if an error occurs during parsing
    */
    void loadTransitions(const String& filename, TargetedExperiment& library, const std::vector<FileTypes::Type> allowed_types = {}, ProgressLogger::LogType log = ProgressLogger::NONE);

    /**
      @brief Store transitions of a spectral library

      @param[out] filename the file name of the file to write.
      @param[out] library The TargetedExperiment to store.
      @param[in] allowed_types A vector of supported filetypes. If empty we try to guess based on the filename. If that fails we throw UnableToCreateFile. If there is only one allowed type, check whether it agrees with the filename, and throw UnableToCreateFile if they disagree.
      @param[in] log Progress logging mode

      @exception Exception::UnableToCreateFile is thrown if the file could not be written
    */
    void storeTransitions(const String& filename, const TargetedExperiment& library, const std::vector<FileTypes::Type> allowed_types = {}, ProgressLogger::LogType log = ProgressLogger::NONE);

        /**
      @brief Loads a file into Transformations

      @param[in] filename the file name of the file to load.
      @param[out] map The Transformations to load the data into.
      @param[in] fit_model Call fitModel() on the @p map before returning?
      @param[in] allowed_types A vector of supported filetypes. If the vector is empty, load from any type that we have a handler for. Otherwise @p getType() is called internally to check the type

      @exception Exception::FileNotFound is thrown if the file could not be opened
      @exception Exception::ParseError is thrown if an error occurs during parsing
    */

    void loadTransformations(const String& filename, TransformationDescription& map, bool fit_model=true, const std::vector<FileTypes::Type> allowed_types = {});

    /**
      @brief Store Transformations

      @param[out] filename the file name of the file to write.
      @param[in] map The Transformations to store.
      @param[in] allowed_types A vector of supported filetypes. If empty we try to guess based on the filename. If that fails we throw UnableToCreateFile. If there is only one allowed type, check whether it agrees with the filename, and throw UnableToCreateFile if they disagree.

      @exception Exception::UnableToCreateFile is thrown if the file could not be written
    */
    void storeTransformations(const String& filename, const TransformationDescription& map, const std::vector<FileTypes::Type> allowed_types = {});

    /**
      @brief Store QC info

      @brief Stores QC data in mzQC file with JSON format
      @param[in] input_file mzML input file name
      @param[in] filename mzQC output file name
      @param[in] exp MSExperiment to extract QC data from, prior sortSpectra() and updateRanges() required
      @param[in] feature_map FeatureMap from feature file (featureXML)
      @param[in] prot_ids protein identifications from ID file (idXML)
      @param[in] pep_ids protein identifications from ID file (idXML)
      @param[out] consensus_map an optional consensus map to store.
      @param[in] contact_name name of the person creating the mzQC file
      @param[in] contact_address contact address (mail/e-mail or phone) of the person creating the mzQC file
      @param[in] description description and comments about the mzQC file contents
      @param[in] label unique and informative label for the run
      @param[in] remove_duplicate_features whether to remove duplicate features only for QCML for now
      @param[in] allowed_types A vector of supported filetypes. If empty we try to guess based on the filename. If that fails we throw UnableToCreateFile. If there is only one allowed type, check whether it agrees with the filename, and throw UnableToCreateFile if they disagree.

      @exception Exception::UnableToCreateFile is thrown if the file could not be written
    */
    void storeQC(const String& input_file,
               const String& filename,
               const MSExperiment& exp,
               const FeatureMap& feature_map,
               std::vector<ProteinIdentification>& prot_ids,
               PeptideIdentificationList& pep_ids,
               const ConsensusMap& consensus_map,
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

