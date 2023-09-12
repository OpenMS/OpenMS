// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/FORMAT/OPTIONS/PeakFileOptions.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

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

      First the type is determined from the file name.
      If this fails, the type is determined from the file content.

      @exception Exception::FileNotFound is thrown if the file is not present
    */
    static FileTypes::Type getType(const String& filename);


    /// Determines the file type from a file name
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
    */
    static String stripExtension(const String& filename);

    /**
      @brief Tries to find and remove a known file extension, and append the new one.

      Internally calls 'stripExtension()' and adds the new suffix to the result.
      E.g. 'experiment.featureXML'+ FileTypes::TRAFOXML becomes 'experiment.trafoXML' and 'c:\\files\\data.mzML.gz' + FileTypes::FEATUREXML becomes 'c:\\files\\data.featureXML'
      If the existing extension is unknown, the everything after the last '.' is removed, e.g. 'exp.tmp'+FileTypes::IDXML becomes 'exp.idXML'
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

    /// Returns if the file type is supported in this build of the library
    static bool isSupported(FileTypes::Type type);

    /// Mutable access to the options for loading/storing
    PeakFileOptions& getOptions();

    /// Non-mutable access to the options for loading/storing
    const PeakFileOptions& getOptions() const;

    /// set options for loading/storing
    void setOptions(const PeakFileOptions&);

    /**
      @brief Loads a file into an MSExperiment

      @param filename The file name of the file to load.
      @param exp The experiment to load the data into.
      @param force_type Forces to load the file with that file type. If no type is forced, it is determined from the extension (or from the content if that fails).
      @param log Progress logging mode
      @param rewrite_source_file Set's the SourceFile name and path to the current file. Note that this looses the link to the primary MS run the file originated from.
      @param compute_hash If source files are rewritten, this flag triggers a recomputation of hash values. A SHA1 string gets stored in the checksum member of SourceFile.

      @return true if the file could be loaded, false otherwise

      @exception Exception::FileNotFound is thrown if the file could not be opened
      @exception Exception::ParseError is thrown if an error occurs during parsing
    */
    bool loadExperiment(const String& filename, MSExperiment& exp, FileTypes::Type force_type = FileTypes::UNKNOWN, ProgressLogger::LogType log = ProgressLogger::NONE, const bool rewrite_source_file = true, const bool compute_hash = true);

    /**
      @brief Stores an MSExperiment to a file

      The file type to store the data in is determined by the file name. Supported formats for storing are mzML, mzXML, mzData and DTA2D. If the file format cannot be determined from the file name, the mzML format is used.

      @param filename The name of the file to store the data in.
      @param exp The experiment to store.
      @param log Progress logging mode

      @exception Exception::UnableToCreateFile is thrown if the file could not be written
    */
    void storeExperiment(const String& filename, const MSExperiment& exp, ProgressLogger::LogType log = ProgressLogger::NONE);

    /**
      @brief Loads a file into a FeatureMap

      @param filename the file name of the file to load.
      @param map The FeatureMap to load the data into.
      @param force_type Forces to load the file with that file type. If no type is forced, it is determined from the extension (or from the content if that fails).

      @return true if the file could be loaded, false otherwise

      @exception Exception::FileNotFound is thrown if the file could not be opened
      @exception Exception::ParseError is thrown if an error occurs during parsing
    */
    bool loadFeatures(const String& filename, FeatureMap& map, FileTypes::Type force_type = FileTypes::UNKNOWN);

    /**
      @brief Store a FeatureMap

      @param filename the file name of the file to write.
      @param map The FeatureMap to store.

      @return true if the file could be stored, false otherwise

      @exception Exception::UnableToCreateFile is thrown if the file could not be written
    */
    bool storeFeatures(const String& filename, const FeatureMap& map);

    /**
      @brief Store a ConsensFeatureMap

      @param filename the file name of the file to write.
      @param map The ConsensusMap to store.

      @return true if the file could be stored, false otherwise

      @exception Exception::UnableToCreateFile is thrown if the file could not be written
    */
    bool storeConsensusFeatures(const String& filename, const ConsensusMap& map);

    /**
      @brief Loads a file into a ConsensMap

      @param filename the file name of the file to load.
      @param map The ConsensMap to load the data into.

      @return true if the file could be loaded, false otherwise

      @exception Exception::FileNotFound is thrown if the file could not be opened
      @exception Exception::ParseError is thrown if an error occurs during parsing
    */
    bool loadConsensusFeatures(const String& filename, ConsensusMap& map);

    bool loadIdentifications(const String& filename, std::vector<ProteinIdentification> additional_proteins, std::vector<PeptideIdentification> additional_peptides);

    /**
      @brief Store transitions of a spectral library

      @param filename the file name of the file to write.
      @param map The TargetedExperiment to store.

      @return true if the file could be stored, false otherwise

      @exception Exception::UnableToCreateFile is thrown if the file could not be written
    */
    bool storeTransitions(const String& filename, const TargetedExperiment& library);

    /**
      @brief Computes a SHA-1 hash value for the content of the given file.

      @return The SHA-1 hash of the given file.
    */
    static String computeFileHash(const String& filename);

private:
    PeakFileOptions options_;

  };

} //namespace

