// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/METADATA/AnnotatedMSRun.h>

#include <vector>

/**
 * @brief Stores spectral data and annotations from an AnnotatedMSRun object into an MSP file.
 *
 * @param filename Name of the output MSP file.
 * @param exp AnnotatedMSRun object containing spectra and associated annotations to be stored.
 *
 * @throw UnableToCreateFile if the file cannot be created or opened for writing.
 */
namespace OpenMS
{
  class AnnotatedMSRun;
  /**
      @brief File adapter for MSP files (NIST spectra library)


      @htmlinclude OpenMS_MSPFile.parameters

      @ingroup FileIO
  */
  class OPENMS_DLLAPI MSPFile :
    public DefaultParamHandler
  {
public:

    /** Constructors and destructors
    */
    //@{
    ///Default constructor
    MSPFile();

    /// Copy constructor
    MSPFile(const MSPFile & rhs);

    ///Destructor
    ~MSPFile() override;
    //@}

    /// assignment operator
    MSPFile & operator=(const MSPFile & rhs);

    /**
        @brief Loads a map from a MSPFile file.

        @param exp PeakMap which contains the spectra after reading
        @param filename the filename of the experiment
        @param ids output parameter which contains the peptide identifications from the spectra annotations

        @throw FileNotFound is thrown if the file could not be found
        @throw ParseError is thrown if the given file could not be parsed
        @throw ElementNotFound is thrown if a annotated modification cannot be found in ModificationsDB (PSI-MOD definitions)
    */
    void load(const String & filename, std::vector<PeptideIdentification> & ids, PeakMap & exp);

    /**
        @brief Loads a map from a MSPFile file.

        @param filename the filename of the experiment
        @param annot_exp annotated experiment with spectra and ids

        @throw FileNotFound is thrown if the file could not be found
        @throw ParseError is thrown if the given file could not be parsed
        @throw ElementNotFound is thrown if a annotated modification cannot be found in ModificationsDB (PSI-MOD definitions)
    */
    void load(const String & filename, AnnotatedMSRun & annot_exp);

    /**
        @brief Stores a map in a MSPFile file.

        @throw UnableToCreateFile is thrown if the given file could not be created
    */
    void store(const String & filename, const AnnotatedMSRun & exp) const;

protected:

    /// reads the header information and stores it as metainfo in the spectrum
    void parseHeader_(const String & header, PeakSpectrum & spec);
  };

} // namespace OpenMS

