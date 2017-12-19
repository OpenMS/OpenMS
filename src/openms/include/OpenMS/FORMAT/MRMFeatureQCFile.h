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
// $Maintainer: Douglas McCloskey $
// $Authors: Douglas McCloskey $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_MRMFEATUREQCFILE_H
#define OPENMS_FORMAT_MRMFEATUREQCFILE_H

#include <OpenMS/FORMAT/CsvFile.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/DATASTRUCTURES/StringListUtils.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMFeatureQC.h>

namespace OpenMS
{
  /**
      @brief File adapter for MRMFeatureQC files.

      loads and stores .csv or .tsv files describing an MRMFeatureQC.

      @ingroup FileIO
  */
  class OPENMS_DLLAPI MRMFeatureQCFile :
    private CsvFile,
    public ProgressLogger
  {
public:
  ///Default constructor
  MRMFeatureQCFile();
  ///Destructor
  ~MRMFeatureQCFile() override;

  /**
    @brief Loads an MRMFeatureQC file.

    @exception Exception::FileNotFound is thrown if the file could not be opened
    @exception Exception::ParseError is thrown if an error occurs during parsing
  */
  void load(const String & filename, MRMFeatureQC & mrmfqc);

  /**
    @brief Stores an MRMFeatureQC file.

    @exception Exception::UnableToCreateFile is thrown if the file could not be created
  */
  void store(const String & filename, const MRMFeatureQC & mrmfqc);


protected:
  /**
    @brief Checks if a file is valid with respect to the mapping file and the controlled vocabulary.

    @param line Header line of the .csv file.
    @param headers A map of header strings to column positions.
    @param params_headers A map of transformation model parameter header strings to column positions.
  */
  void parseHeader_(StringList & line, std::map<String, int> & headers,
    std::map<String, int> & params_headers);

  /**
    @brief parses a line into the members of MRMFeatureQC.

    @param line line of the .csv file.
    @param aqm MRMFeatureQC.
  */
  void parseLine_(StringList & line, std::map<String, int> & headers, 
    std::map<String, int> & params_headers, MRMFeatureQC & mrmfqc);

  };

} // namespace OpenMS

#endif // OPENMS_FORMAT_MRMFEATUREQCFILE_H
