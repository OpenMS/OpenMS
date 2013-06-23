// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
// $Maintainer: Andreas Bertsch $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_FILEHANDLER_H
#define OPENMS_FORMAT_FILEHANDLER_H

#include <OpenMS/config.h>
#include <OpenMS/FORMAT/FileTypes.h>

#include <OpenMS/FORMAT/DTAFile.h>
#include <OpenMS/FORMAT/DTA2DFile.h>
#include <OpenMS/FORMAT/MzXMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/MzDataFile.h>
#include <OpenMS/FORMAT/MascotGenericFile.h>
#include <OpenMS/FORMAT/MS2File.h>
#include <OpenMS/FORMAT/XMassFile.h>

#include <OpenMS/FORMAT/MsInspectFile.h>
#include <OpenMS/FORMAT/SpecArrayFile.h>
#include <OpenMS/FORMAT/KroenikFile.h>

#include <OpenMS/KERNEL/ChromatogramTools.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>

namespace OpenMS
{
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
    static FileTypes::Type getType(const String & filename);


    /// Determines the file type from a file name
    static FileTypes::Type getTypeByFileName(const String & filename);

    /**
      @brief Determines the file type of a file by parsing the first few lines

      @exception Exception::FileNotFound is thrown if the file is not present
    */
    static FileTypes::Type getTypeByContent(const String & filename);

    /// Returns if the file type is supported in this build of the library
    static bool isSupported(FileTypes::Type type);

    /// Mutable access to the options for loading/storing
    PeakFileOptions & getOptions();

    /// Non-mutable access to the options for loading/storing
    const PeakFileOptions & getOptions() const;

    /**
        @brief Loads a file into an MSExperiment

        @param filename The file name of the file to load.
        @param exp The experiment to load the data into.
        @param force_type Forces to load the file with that file type. If no type is forced, it is determined from the extention ( or from the content if that fails).
        @param log Progress logging mode

        @return true if the file could be loaded, false otherwise

        @exception Exception::FileNotFound is thrown if the file could not be opened
        @exception Exception::ParseError is thrown if an error occurs during parsing
    */
    template <class PeakType>
    bool loadExperiment(const String & filename, MSExperiment<PeakType> & exp, FileTypes::Type force_type = FileTypes::UNKNOWN, ProgressLogger::LogType log = ProgressLogger::NONE)
    {
      //determine file type
      FileTypes::Type type;
      if (force_type != FileTypes::UNKNOWN)
      {
        type = force_type;
      }
      else
      {
        try
        {
          type = getType(filename);
        }
        catch (Exception::FileNotFound)
        {
          return false;
        }
      }

      //load right file
      switch (type)
      {
      case FileTypes::DTA:
        exp.reset();
        exp.resize(1);
        DTAFile().load(filename, exp[0]);
        break;

      case FileTypes::DTA2D:
        {
          DTA2DFile f;
          f.getOptions() = options_;
          f.setLogType(log);
          f.load(filename, exp); 
        }

        break;

      case FileTypes::MZXML:
        {
          MzXMLFile f;
          f.getOptions() = options_;
          f.setLogType(log);
          f.load(filename, exp); 
        }

        break;

      case FileTypes::MZDATA:
        {
          MzDataFile f;
          f.getOptions() = options_;
          f.setLogType(log);
          f.load(filename, exp);
        }
        break;

      case FileTypes::MZML:
        {
          MzMLFile f;
          f.getOptions() = options_;
          f.setLogType(log);
          f.load(filename, exp);
          ChromatogramTools().convertSpectraToChromatograms<MSExperiment<PeakType> >(exp, true); 
        }
        break;

      case FileTypes::MGF:
        {
          MascotGenericFile f;
          f.setLogType(log);
          f.load(filename, exp); 
        }

        break;

      case FileTypes::MS2:
        {
          MS2File f;
          f.setLogType(log);
          f.load(filename, exp); 
        }

        break;

      case FileTypes::XMASS:
          exp.reset();
          exp.resize(1);
          XMassFile().load(filename, exp[0]);
          XMassFile().importExperimentalSettings(filename, exp); 

        break;

      default:
        return false;
        break;
      }

      SourceFile src_file;
      src_file.setNameOfFile(File::basename(filename));
      src_file.setPathToFile(String("file:///") + File::path(filename));
      // this is more complicated since the data formats allowed by mzML are very verbose.
      // this is prone to changing CV's... our writer will fall back to a default if the name given here is invalid.
      src_file.setFileType(FileTypes::typeToMZML(type));
      
      exp.getSourceFiles().clear();
      exp.getSourceFiles().push_back(src_file);
      
      return true;
    }

    /**
        @brief Stores an MSExperiment to a file

        The file type to store the data in is determined by the file name. Supported formats for storing are mzML, mzXML, mzData and DTA2D. If the file format cannot be determined from the file name, the mzML format is used.

        @param filename The name of the file to store the data in.
        @param exp The experiment to store.
        @param log Progress logging mode

        @exception Exception::UnableToCreateFile is thrown if the file could not be written
    */
    template <class PeakType>
    void storeExperiment(const String & filename, const MSExperiment<PeakType> & exp, ProgressLogger::LogType log = ProgressLogger::NONE)
    {
      //load right file
      switch (getTypeByFileName(filename))
      {
      case FileTypes::DTA2D:
      {
        DTA2DFile f;
        f.getOptions() = options_;
        f.setLogType(log);
        f.store(filename, exp);
      }
      break;

      case FileTypes::MZXML:
      {
        MzXMLFile f;
        f.getOptions() = options_;
        f.setLogType(log);
        if (!exp.getChromatograms().empty())
        {
          MSExperiment<PeakType> exp2 = exp;
          ChromatogramTools().convertChromatogramsToSpectra<MSExperiment<PeakType> >(exp2);
          f.store(filename, exp2);
        }
        else
        {
          f.store(filename, exp);
        }
      }
      break;

      case FileTypes::MZDATA:
      {
        MzDataFile f;
        f.getOptions() = options_;
        f.setLogType(log);
        if (!exp.getChromatograms().empty())
        {
          MSExperiment<PeakType> exp2 = exp;
          ChromatogramTools().convertChromatogramsToSpectra<MSExperiment<PeakType> >(exp2);
          f.store(filename, exp2);
        }
        else
        {
          f.store(filename, exp);
        }
      }
      break;

      default:
      {
        MzMLFile f;
        f.getOptions() = options_;
        f.setLogType(log);
        f.store(filename, exp);
      }
      break;
      }
    }

    /**
        @brief Loads a file into a FeatureMap

        @param filename the file name of the file to load.
        @param map The FeatureMap to load the data into.
        @param force_type Forces to load the file with that file type. If no type is forced, it is determined from the extention ( or from the content if that fails).

        @return true if the file could be loaded, false otherwise

        @exception Exception::FileNotFound is thrown if the file could not be opened
        @exception Exception::ParseError is thrown if an error occurs during parsing
    */
    template <class FeatureType>
    bool loadFeatures(const String & filename, FeatureMap<FeatureType> & map, FileTypes::Type force_type = FileTypes::UNKNOWN)
    {
      //determine file type
      FileTypes::Type type;
      if (force_type != FileTypes::UNKNOWN)
      {
        type = force_type;
      }
      else
      {
        try
        {
          type = getType(filename);
        }
        catch (Exception::FileNotFound)
        {
          return false;
        }
      }

      //load right file
      if (type == FileTypes::FEATUREXML)
      {
        FeatureXMLFile().load(filename, map);
      }
      else if (type == FileTypes::TSV)
      {
        MsInspectFile().load(filename, map);
      }
      else if (type == FileTypes::PEPLIST)
      {
        SpecArrayFile().load(filename, map);
      }
      else if (type == FileTypes::KROENIK)
      {
        KroenikFile().load(filename, map);
      }
      else
      {
        return false;
      }

      return true;
    }

private:
    PeakFileOptions options_;
  };

} //namespace

#endif //OPENMS_FORMAT_FILEHANDLER_H
