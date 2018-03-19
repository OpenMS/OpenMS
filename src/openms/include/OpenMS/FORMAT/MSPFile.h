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
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_MSPFILE_H
#define OPENMS_FORMAT_MSPFILE_H

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/KERNEL/StandardTypes.h>

#include <vector>

namespace OpenMS
{
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
        @brief Stores a map in a MSPFile file.

        @throw UnableToCreateFile is thrown if the given file could not be created
    */
    void store(const String & filename, const PeakMap & exp) const;

protected:

    /// reads the header information and stores it as metainfo in the spectrum
    void parseHeader_(const String & header, PeakSpectrum & spec);

  };

} // namespace OpenMS

#endif // OPENMS_FORMAT_MSPFILE_H
