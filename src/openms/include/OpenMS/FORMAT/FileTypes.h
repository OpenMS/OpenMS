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
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch, Marc Sturm, Stephan Aiche $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_FILETYPES_H
#define OPENMS_FORMAT_FILETYPES_H

#include <OpenMS/config.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <string>
#include <map>

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
      TSV,                ///< msInspect file (.tsv)
      PEPLIST,            ///< specArray file (.peplist)
      HARDKLOER,          ///< hardkloer file (.hardkloer)
      KROENIK,            ///< kroenik file (.kroenik)
      FASTA,              ///< FASTA file (.fasta)
      EDTA,               ///< enhanced comma separated files (RT, m/z, Intensity, [meta])
      CSV,                ///< general comma separated files format (might also be tab or space separated!!!), data should be regular, i.e. matrix form
      TXT,                ///< any text format, which has only loose definition of what it actually contains -- thus it is usually hard to say where the file actually came from (e.g. PepNovo).
      OBO,                ///< Controlled Vocabulary format
      HTML,               ///< any HTML format
      XML,                ///< any XML format
      ANALYSISXML,        ///< analysisXML format
      XSD,                ///< XSD schema format
      PSQ,                ///< NCBI binary blast db
      MRM,                ///< SpectraST MRM List
      SQMASS,             ///< SqLite format for mass and chromatograms
      PQP,                ///< OpenSWATH Peptide Query Parameter (PQP) SQLite DB
      OSW,                ///< OpenSWATH OpenSWATH report (OSW) SQLite DB
      PSMS,               ///< Percolator tab-delimited output (PSM level)
      PIN,                ///< Percolator tab-delimited input (PSM level)
      PARAMXML,           ///< internal format for writing and reading parameters (also used as part of CTD)
      SPLIB,              ///< SpectraST binary spectral library file (sptxt is the equivalent text-based format, similar to the MSP format)
      NOVOR,               ///< Novor custom parameter file
      SIZE_OF_TYPE        ///< No file type. Simply stores the number of types
    };

    /// Returns the name/extension of the type.
    static String typeToName(Type type);

    /// Returns the mzML name (TODO: switch to accession since they are more stable!)
    static String typeToMZML(Type type);

    /// Converts a file type name into a Type
    static Type nameToType(const String& name);

private:
    /// Maps the FileType::Type to the preferred extension.
    static const std::map<Type, String> name_of_types_;
    
    /// Maps the FileType::Type to the preferred mzML CV name.
    static const std::map<Type, String> name_of_MZMLtypes_;

    /// Initializer for the file extension map.
    static std::map<Type, String> initializeMap_();

    /// Initializer for the file extension map.
    static std::map<Type, String> initializeMZMLMap_();

  };

} //namespace OpenMS

#endif //OPENMS_FORMAT_FILETYPES_H
