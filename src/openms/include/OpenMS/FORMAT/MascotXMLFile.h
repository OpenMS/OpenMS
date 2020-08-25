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

