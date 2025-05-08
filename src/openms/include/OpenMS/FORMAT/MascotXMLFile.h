// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Nico Pfeifer $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/FORMAT/XMLFile.h>
#include <OpenMS/FORMAT/HANDLERS/MascotXMLHandler.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/METADATA/PeptideIdentification.h>


namespace OpenMS
{
  class ProteinIdentification;

  /**
    @brief Used to load Mascot XML files

    This class is used to load documents that implement
    the schema of Mascot XML files.

    @ingroup FileIO
  */
  class OPENMS_DLLAPI MascotXMLFile :
    public Internal::XMLFile
  {
public:

    /// Constructor
    MascotXMLFile();

    /**
      @brief Loads data from a Mascot XML file

      @param filename the file to be loaded
      @param protein_identification protein identifications belonging to the whole experiment
      @param id_data the identifications with m/z and RT
      @param lookup helper object for looking up spectrum meta data

      @exception Exception::FileNotFound is thrown if the file does not exists.
      @exception Exception::ParseError is thrown if the file does not suit to the standard.
    */
    void load(const String& filename,
              ProteinIdentification& protein_identification,
              std::vector<PeptideIdentification>& id_data,
              const SpectrumMetaDataLookup& lookup);

    /**
      @brief Loads data from a Mascot XML file

      @param filename the file to be loaded
      @param protein_identification protein identifications belonging to the whole experiment
      @param id_data the identifications with m/z and RT
      @param peptides a map of modified peptides identified by the String title
      @param lookup helper object for looking up spectrum meta data

      @exception Exception::FileNotFound is thrown if the file does not exists.
      @exception Exception::ParseError is thrown if the file does not suit to the standard.
    */
    void load(const String& filename,
              ProteinIdentification& protein_identification,
              std::vector<PeptideIdentification>& id_data, 
              std::map<String, std::vector<AASequence> >& peptides, 
              const SpectrumMetaDataLookup& lookup);

    /**
      @brief Initializes a helper object for looking up spectrum meta data (RT, m/z)

      @param lookup Helper object to initialize
      @param experiment Experiment containing the spectra
      @param scan_regex Optional regular expression for extracting information from references to spectra
    */  
    static void initializeLookup(SpectrumMetaDataLookup& lookup, const PeakMap& experiment, const String& scan_regex = "");

  };

} // namespace OpenMS

