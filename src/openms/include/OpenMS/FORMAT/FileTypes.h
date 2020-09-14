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
// $Maintainer: Timo Sachsenberg $
// $Authors: Stephan Aiche, Andreas Bertsch, Marc Sturm, Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <array>

namespace OpenMS
{
  /**
    @brief Centralizes the file types recognized by FileHandler.

    FileType separate from FileHandler to avoid circular inclusions by DocumentIdentifier, ExperimentalSettings
    and FileHandler and respective fileclasses (e.g. DTA2DFile).

    @ingroup FileIO
  */
  struct OPENMS_DLLAPI FileTypes
  {
    //NOTE: if you change/add something here, do not forget to change FileTypes::initializeMap_

    ///Actual file types enum.
    enum Type
    {
      UNKNOWN,            ///< Unknown file extension
      DTA,                ///< DTA file (.dta)
      DTA2D,              ///< DTA2D file (.dta2d)
      MZDATA,             ///< MzData file (.mzData)
      MZXML,              ///< MzXML file (.mzXML)
      FEATUREXML,         ///< %OpenMS feature file (.featureXML)
      IDXML,              ///< %OpenMS identification format (.idXML)
      CONSENSUSXML,       ///< %OpenMS consensus map format (.consensusXML)
      MGF,                ///< Mascot Generic Format (.mgf)
      INI,                ///< %OpenMS parameters file (.ini)
      TOPPAS,             ///< %OpenMS parameters file with workflow information (.toppas)
      TRANSFORMATIONXML,  ///< Transformation description file (.trafoXML)
      MZML,               ///< MzML file (.mzML)
      CACHEDMZML,         ///< CachedMzML file (.cachedmzML)
      MS2,                ///< MS2 file (.ms2)
      PEPXML,             ///< TPP pepXML file (.pepXML)
      PROTXML,            ///< TPP protXML file (.protXML)
      MZIDENTML,          ///< mzIdentML (HUPO PSI AnalysisXML followup format) (.mzid)
      MZQUANTML,          ///< mzQuantML (HUPO PSI AnalysisXML followup format) (.mzq)
      QCML,               ///< qcML (will undergo standardisation maybe) (.qcml)
      GELML,              ///< GelML (HUPO PSI format) (.gelML)
      TRAML,              ///< TraML (HUPO PSI format) for transitions (.traML)
      MSP,                ///< NIST spectra library file format (.msp)
      OMSSAXML,           ///< OMSSA XML file format for peptide identifications (.xml)
      MASCOTXML,          ///< Mascot XML file format for peptide identifications (.xml)
      PNG,                ///< Portable Network Graphics (.png)
      XMASS,              ///< XMass Analysis file (fid)
      TSV,                ///< any TSV file, for example msInspect file or OpenSWATH transition file (see TransitionTSVFile)
      MZTAB,              ///< mzTab file (.mzTab)
      PEPLIST,            ///< specArray file (.peplist)
      HARDKLOER,          ///< hardkloer file (.hardkloer)
      KROENIK,            ///< kroenik file (.kroenik)
      FASTA,              ///< FASTA file (.fasta)
      EDTA,               ///< enhanced comma separated files (RT, m/z, Intensity, [meta])
      CSV,                ///< general comma separated files format (might also be tab or space separated!!!), data should be regular, i.e. matrix form
      TXT,                ///< any text format, which has only loose definition of what it actually contains -- thus it is usually hard to say where the file actually came from (e.g. PepNovo).
      OBO,                ///< Controlled Vocabulary format
      HTML,               ///< any HTML format
      ANALYSISXML,        ///< analysisXML format
      XSD,                ///< XSD schema format
      PSQ,                ///< NCBI binary blast db
      MRM,                ///< SpectraST MRM List
      SQMASS,             ///< SqLite format for mass and chromatograms, see SqMassFile
      PQP,                ///< OpenSWATH Peptide Query Parameter (PQP) SQLite DB, see TransitionPQPFile
      MS,                 ///< SIRIUS file format (.ms)
      OSW,                ///< OpenSWATH OpenSWATH report (OSW) SQLite DB
      PSMS,               ///< Percolator tab-delimited output (PSM level)
      PIN,                ///< Percolator tab-delimited input (PSM level)
      PARAMXML,           ///< internal format for writing and reading parameters (also used as part of CTD)
      SPLIB,              ///< SpectraST binary spectral library file (sptxt is the equivalent text-based format, similar to the MSP format)
      NOVOR,              ///< Novor custom parameter file
      XQUESTXML,          ///< xQuest XML file format for protein-protein cross-link identifications (.xquest.xml)
      SPECXML,            ///< xQuest XML file format for matched spectra for spectra visualization in the xQuest results manager (.spec.xml)
      JSON,               ///< JavaScript Object Notation file (.json)
      RAW,                ///< Thermo Raw File (.raw)
      EXE,                ///< Executable (.exe)
      XML,                ///< any XML format
      BZ2,                ///< any BZ2 compressed file
      GZ,                 ///< any Gzipped file
      SIZE_OF_TYPE        ///< No file type. Simply stores the number of types
    };


    enum class Filter
    {
      COMPACT,    ///< make a single item, e.g. 'all readable files (*.mzML *.mzXML);;'
      ONE_BY_ONE, ///< list all types individually, e.g. 'mzML files (*.mzML);;mzXML files (*.mzXML);;'
      BOTH        ///< combine COMPACT and ONE_BY_ONE
    };
    /**
      @brief holds a known array of file types, e.g. as a way to specify supported input formats

      The array can be exported in Qt's file dialog format.

    */
    template<size_t N>
    struct FileTypeList
    {
      FileTypeList(const std::array<Type, N> types)
        : type_list(types)
      {
      }

      /// check if @p type is contained in this array
      bool contains(const Type& type)
      {
        for (const auto& t : type_list)
        {
          if (t == type) return true;
        }
        return false;
      }

      /// converts the array into a Qt-compatible filter for selecting files in a user dialog.
      /// e.g. "all readable files (*.mzML *.mzXML);;". See Filter enum.
      /// @param style Create a combined filter, or single filters, or both
      /// @param add_all_filter Add 'all files (*)' as a single filter at the end?
      String toFileDialogFilter(const Filter style, bool add_all_filter)
      {
        String out;
        if (style == Filter::COMPACT || style == Filter::BOTH)
        {
          StringList items;
          for (const auto& t : type_list)
          {
            items.push_back("*." + FileTypes::typeToName(t));
          }
          out += "all readable files (" + ListUtils::concatenate(items, " ") + ");;";
        }
        if (style == Filter::ONE_BY_ONE || style == Filter::BOTH)
        {
          StringList items;
          for (const auto& t : type_list)
          {
            items.push_back(FileTypes::typeToDescription(t) + " (*." + FileTypes::typeToName(t) + ");;");
          }
          out += ListUtils::concatenate(items, "");
        }
        if (add_all_filter) out += "all files (*);;";

        // remove the last ";;", since this will be interpreted as ' (*)' by Qt
        out = out.chop(2);
        
        return out;
      }

      std::array<Type, N> type_list;
    };

    /// Returns the name/extension of the type.
    static String typeToName(Type type);
    
    /// Returns the human-readable explanation of the type.
    /// This may or may not add information, e.g.
    /// MZML becomes "mzML raw data file", but FEATUREXML becomes "OpenMS feature map"
    static String typeToDescription(Type type);
    
    /// Converts a file type name into a Type
    static Type nameToType(const String& name);

    /// Returns the mzML name (TODO: switch to accession since they are more stable!)
    static String typeToMZML(Type type);
  };

} //namespace OpenMS

