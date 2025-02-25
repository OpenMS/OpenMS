// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Stephan Aiche, Andreas Bertsch, Marc Sturm, Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <vector>

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
      QCML,               ///< qcML (will undergo standardisation maybe) (.qcml)
      MZQC,               ///< mzQC (HUPO PSI format) (.mzQC)
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
      OMS,                ///< OpenMS database file
      EXE,                ///< Executable (.exe)
      XML,                ///< any XML format
      BZ2,                ///< any BZ2 compressed file
      GZ,                 ///< any Gzipped file
      SIZE_OF_TYPE        ///< No file type. Simply stores the number of types
    };

    enum class FileProperties
    {
      READABLE,                     // SOMETHING in OpenMS can read this (it doesn't have to be in FileHandler though)
      WRITEABLE,                    // SOMETHING in OpenMS can write this (it doesn't have to be in FileHandler though) 
      PROVIDES_SPECTRUM,            // All of the PROVIDES_x properties correspond to which FileHandlers are implemented for a file type.
      PROVIDES_EXPERIMENT,          // 
      PROVIDES_FEATURES,            //
      PROVIDES_CONSENSUSFEATURES,   //
      PROVIDES_IDENTIFICATIONS,     //
      PROVIDES_TRANSITIONS,         //
      PROVIDES_QUANTIFICATIONS,     //
      PROVIDES_TRANSFORMATIONS,     //
      PROVIDES_QC,                  //
      SIZE_OF_FILEPROPERTIES        // Not a property, just the number of 'em
    };

    /// Returns the name/extension of the type.
    static String typeToName(Type type);
    
    /// Returns the human-readable explanation of the type.
    /// This may or may not add information, e.g.
    /// MZML becomes "mzML raw data file", but FEATUREXML becomes "OpenMS feature map"
    static String typeToDescription(Type type);
    
    /// Converts a file type name into a Type 
    /// @param name A case-insensitive name (e.g. FASTA or Fasta, etc.)
    static Type nameToType(const String& name);

    /// Returns the mzML name (TODO: switch to accession since they are more stable!)
    static String typeToMZML(Type type);
  };

 enum class FilterLayout
  {
    COMPACT,    ///< make a single item, e.g. 'all readable files (*.mzML *.mzXML);;'
    ONE_BY_ONE, ///< list all types individually, e.g. 'mzML files (*.mzML);;mzXML files (*.mzXML);;'
    BOTH        ///< combine COMPACT and ONE_BY_ONE
  };
  /**
    @brief holds a vector of known file types, e.g. as a way to specify supported input formats

    The vector can be exported in Qt's file dialog format.
  */
  class OPENMS_DLLAPI FileTypeList
  {
  public:
    FileTypeList(const std::vector<FileTypes::Type>& types);

    /// check if @p type is contained in this array
    bool contains(const FileTypes::Type& type) const;

    const std::vector<FileTypes::Type>& getTypes() const
    {
      return type_list_;
    }

    /// converts the array into a Qt-compatible filter for selecting files in a user dialog.
    /// e.g. "all readable files (*.mzML *.mzXML);;". See Filter enum.
    /// @param style Create a combined filter, or single filters, or both
    /// @param add_all_filter Add 'all files (*)' as a single filter at the end?
    String toFileDialogFilter(const FilterLayout style, bool add_all_filter) const;

    /**
      @brief Convert a Qt filter back to a Type if possible.

      E.g. from a full filter such as '"mzML files (*.mzML);;mzData files (*.mzData);;mzXML files (*.mzXML);;all files (*)"',
      as created by toFileDialogFilter(), the selected @p filter could be "mzML files (*.mzML)", in which case the type is Type::MZML .
      However, for the filter "all files (*)", Type::UNKNOWN will be returned.

      If the type is UNKNOWN, then the fallback is returned (by default also UNKNOWN). This is useful if you want a default type to fall back to.

      @param filter The filter returned by 'QFileDialog::getSaveFileName' and others, i.e. an item from the result of 'toFileDialogFilter'.
      @param fallback If the filter is ambiguous, return this type instead
      @return The type associated to the filter or the fallback
      @throw Exception::ElementNotFound if the given @p filter is not a filter produced by toFileDialogFilter()
    **/
    FileTypes::Type fromFileDialogFilter(const String& filter, const FileTypes::Type fallback = FileTypes::Type::UNKNOWN) const;


    /**
      @brief Get a std::vector<FileTypes::Type> with all fileTypes that support a set of features.

      
      @param features A set of features that must be supported
      @return A std::vector<FileTypes::Type> with the files that support features
    **/
    static std::vector<FileTypes::Type> typesWithProperties(const std::vector<FileTypes::FileProperties> features);

  private:
    /// hold filter items (for Qt dialogs) along with their OpenMS type
    struct FilterElements_
    {
      std::vector<String> items;
      std::vector<FileTypes::Type> types;
    };
    /// creates Qt filters and the corresponding elements from type_list_
    /// @param style Create a combined filter, or single filters, or both
    /// @param add_all_filter Add 'all files (*)' as a single filter at the end?
    FilterElements_ asFilterElements_(const FilterLayout style, bool add_all_filter) const;

    std::vector<FileTypes::Type> type_list_;
  };

} // namespace OpenMS

