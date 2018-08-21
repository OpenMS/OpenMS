// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/FORMAT/OPTIONS/PeakFileOptions.h>

namespace OpenMS
{
  class PeakFileOptions;
  class MSSpectrum;
  class MSExperiment;
  class FeatureMap;

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
       @brief Check if the file extension of @param type matches no other known FileType
               
       This means that may match the extension of @type or any other unknown extension (e.g., '.tmp')
    */
    static bool hasValidExtension(const String& filename, const FileTypes::Type type);

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
      @brief Computes a SHA-1 hash value for the content of the given file.

      @return The SHA-1 hash of the given file.
    */
    static String computeFileHash(const String& filename);

private:
    PeakFileOptions options_;

  };

} //namespace

