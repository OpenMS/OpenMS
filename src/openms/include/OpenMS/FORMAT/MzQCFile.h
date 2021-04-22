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
// $Maintainer: Axel Walter $
// $Authors: Axel Walter $
// --------------------------------------------------------------------------

#pragma once
#include <OpenMS/KERNEL/MSExperiment.h>

namespace OpenMS
{
  /**
      @brief File adapter for mzQC files used to load and store mzQC files

      This Class is supposed to internally collect the data for the mzQC File

      @ingroup FileIO
  */
  class OPENMS_DLLAPI MzQCFile
  {
  public:
    // Default constructor
    MzQCFile() = default;

    // Store the mzQC file
    /**
      @brief Stores QC data in mzQC file with JSON format
      @param inputFileName mzML input file name
      @param outputFileName mzQC output file name
      @param exp MSExperiment to extract QC data from, prior sortSpectra() and updateRanges() required
      @param contactName name of the person creating the mzQC file
      @param contactAddress contact address (mail/e-mail or phone) of the person creating the mzQC file
      @param description description and comments about the mzQC file contents
      @param label unique and informative label for the run
    */
    void store(const String & inputFileName,
               const String & outputFileName,
               const MSExperiment & exp,
               const String & contactName,
               const String & contactAddress,
               const String & description,
               const String & label) const;

  };
}