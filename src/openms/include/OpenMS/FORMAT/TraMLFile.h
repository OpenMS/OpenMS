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
// $Authors: Andreas Bertsch$
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/FORMAT/XMLFile.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>

namespace OpenMS
{
  class Identification;
  class TargetedExperiment;

  /**
      @brief File adapter for HUPO PSI TraML files

      TraML files contain information about transitions used for targeted
      proteomics and metabolomics experiments:

      Deutsch et al. "TraML--a standard format for exchange of selected reaction monitoring transition lists."
      Mol Cell Proteomics. 2012 Apr;11(4):R111.015040. doi: 10.1074/mcp.R111.015040. 

      In OpenMS, TraML files can be generated from TSV or CSV files using the
      @ref OpenMS::TransitionTSVFile "TransitionTSVFile class" or the @ref
      UTILS_TargetedFileConverter "TargetedFileConverter TOPP Tool". For more information on the TSV format required by the TOPP tool, see
      see also the documentation of @ref
      OpenMS::TransitionTSVFile "TransitionTSVFile".

      @ingroup FileIO
  */
  class OPENMS_DLLAPI TraMLFile :
    public Internal::XMLFile,
    public ProgressLogger
  {
public:
    ///Default constructor
    TraMLFile();
    ///Destructor
    ~TraMLFile() override;

    /**
        @brief Loads a map from a TraML file.

        @exception Exception::FileNotFound is thrown if the file could not be opened
        @exception Exception::ParseError is thrown if an error occurs during parsing
    */
    void load(const String & filename, TargetedExperiment & id);

    /**
        @brief Stores a map in a TraML file.

        @exception Exception::UnableToCreateFile is thrown if the file could not be created
    */
    void store(const String & filename, const TargetedExperiment & id) const;

    /**
        @brief Checks if a file is valid with respect to the mapping file and the controlled vocabulary.

        @param filename File name of the file to be checked.
        @param errors Errors during the validation are returned in this output parameter.
        @param warnings Warnings during the validation are returned in this output parameter.

        @exception Exception::FileNotFound is thrown if the file could not be opened
    */
    bool isSemanticallyValid(const String & filename, StringList & errors, StringList & warnings);

private:

  };

} // namespace OpenMS


