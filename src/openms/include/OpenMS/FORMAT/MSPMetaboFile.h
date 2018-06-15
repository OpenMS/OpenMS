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
// $Maintainer: Douglas McCloskey, Pasquale Domenico Colaianni $
// $Authors: Douglas McCloskey, Pasquale Domenico Colaianni $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h> // OPENMS_DLLAPI
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
  class OPENMS_DLLAPI MSPMetaboFile
  {
public:
    /// Default constructor
    MSPMetaboFile() = default;

    /// Constructor with filename and output library
    MSPMetaboFile(const String& filename, MSExperiment& library);

    /// Destructor
    ~MSPMetaboFile() = default;

    /// To test private and protected methods
    friend class MSPMetaboFile_friend;

    /**
      @brief Load the file's data and metadata, and save it into an `MSExperiment`.

      @param[in] filename Path to the MSP input file
      @param[out] library The variable into which the extracted information will be saved

      @throw FileNotFound is thrown if the file could not be found
    */
    void load(const String& filename, MSExperiment& library);

private:
    /**
      Push a field of the MSP structure into a named `StringDataArray`

      @param[in/out] spectrum The metadata will be added/updated in this `MSSpectrum`
      @param[in] name The name of the field to add or update
      @param[in] info The value to insert
    */
    void pushParsedInfoToNamedDataArray(
      MSSpectrum& spectrum,
      const String& name,
      const String& info
    ) const;

    /**
      Validate and add a spectrum to a spectral library

      The spectrum is added to the library if all following criteria are met:
      - Name field is present and not empty
      - The number of peaks parsed matches the value of Num Peaks
      - A spectrum of the same name has not already been added

      @param[in/out] spectrum The spectrum to be added
      @param[out] library The spectrum is added into this `MSExperiment` library
    */
    void addSpectrumToLibrary(
      MSSpectrum& spectrum,
      MSExperiment& library
    );

    /// To keep track of which spectra have already been loaded and avoid duplicates
    std::set<String> loaded_spectra_names_;
  };

  class MSPMetaboFile_friend
  {
public:
    MSPMetaboFile_friend() = default;
    ~MSPMetaboFile_friend() = default;

    void pushParsedInfoToNamedDataArray(
      MSSpectrum& spectrum,
      const String& name,
      const String& info
    ) const
    {
      return msp_.pushParsedInfoToNamedDataArray(spectrum, name, info);
    }

    void addSpectrumToLibrary(
      MSSpectrum& spectrum,
      MSExperiment& library
    )
    {
      return msp_.addSpectrumToLibrary(spectrum, library);
    }

    MSPMetaboFile msp_;
  };
}
