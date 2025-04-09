// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Douglas McCloskey, Pasquale Domenico Colaianni $
// $Authors: Douglas McCloskey, Pasquale Domenico Colaianni $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h> // OPENMS_DLLAPI
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/KERNEL/MSExperiment.h>

namespace OpenMS
{
  /**
    @brief Load MSP text file and save it into an `MSExperiment`

    This class is specialized for metabolites data.
    The required fields are: Name, Num Peaks, and the peaks data

    Points (meaning x and y info) may be separated by a space or a colon.
    Peaks may be separated by a space or a semicolon.

    An example of the expected format:
    > Name: foo
    > Num Peaks: 11
    > 35 310; 36 1230; 37 27; 38 303; 47 5240;
    > 66 203; 67 68; 68 77; 82 63; 83 240;
    > 136 350;

    Another supported format:
    > Name: bar
    > Num Peaks: 11
    > 35:310 36:1230 37:27 38:303 47:5240
    > 66:203 67:68 68:77 82:63 83:240
    > 136:350
  */
  class OPENMS_DLLAPI MSPGenericFile :
    public DefaultParamHandler
  {
public:
    /// Default constructor
    MSPGenericFile();

    /// Constructor with filename and output library
    MSPGenericFile(const String& filename, MSExperiment& library);

    /// Destructor
    ~MSPGenericFile() override = default;

    /// Get the class' default parameters
    void getDefaultParameters(Param& params);

    /// To test private and protected methods
    friend class MSPGenericFile_friend;

    /**
      @brief Load the file's data and metadata, and save it into an `MSExperiment`.

      @param[in] filename Path to the MSP input file
      @param[out] library The variable into which the extracted information will be saved

      @throw FileNotFound If the file could not be found
    */
    void load(const String& filename, MSExperiment& library);

    /**
      @brief Save data and metadata into a file.

      @param[in] filename Path to the MSP input file
      @param[out] library The variable from which extracted information will be saved

      @throw FileNotWritable If the file is not writable
    */
    void store(const String& filename, const MSExperiment& library) const;
  
  private:
    /// Overrides `DefaultParamHandler`'s method
    void updateMembers_() override;

    /**
      Validate and add a spectrum to a spectral library

      The spectrum is added to the library if all following criteria are met:
      - Name field is present and not empty
      - The number of peaks parsed matches the value of Num Peaks
      - A spectrum of the same name has not already been added

      @throw MissingInformation If the spectrum doesn't have a name or Num Peaks info is missing
      @throw ParseError If Num Peaks' value doesn't match with the number of raw peaks parsed

      @param[in,out] spectrum The spectrum to be added
      @param[out] library The spectrum is added into this `MSExperiment` library
    */
    void addSpectrumToLibrary(
      MSSpectrum& spectrum,
      MSExperiment& library
    );

    /// To keep track of which spectra have already been loaded and avoid duplicates
    std::set<String> loaded_spectra_names_;

    /*
      The synonyms of a spectrum are collected into this variable and,
      when `addSpectrumtoLibrary()` is called, the elements are concatenated
      and the result is saved as a "Synon" metaValue.
      The synonyms are separated by `synonyms_separator_`.
    */
    std::vector<String> synonyms_;

    /// The separator to be used in "Synon" metaValue
    String synonyms_separator_;
  };

  class MSPGenericFile_friend
  {
public:
    MSPGenericFile_friend() = default;
    ~MSPGenericFile_friend() = default;

    void addSpectrumToLibrary(
      MSSpectrum& spectrum,
      MSExperiment& library
    )
    {
      return msp_.addSpectrumToLibrary(spectrum, library);
    }

    MSPGenericFile msp_;
  };
}
