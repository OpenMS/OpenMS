// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2021.
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
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/FORMAT/TextFile.h>

#include <vector>

namespace OpenMS
{

  /**
     @brief Class for storing Percolator tab-delimited input files.

  */
  class OPENMS_DLLAPI PercolatorInfile
  {
    public:
      static void store(const String& pin_file, 
        const std::vector<PeptideIdentification>& peptide_ids, 
        const StringList& feature_set, 
        const std::string& enz, 
        int min_charge, 
        int max_charge);

      // uses spectrum_reference, if empty uses spectrum_id, if also empty fall back to using index
      static String getScanIdentifier(const PeptideIdentification& pid, size_t index);
      
    protected:

      //id <tab> label <tab> scannr <tab> calcmass <tab> expmass <tab> feature1 <tab> ... <tab> featureN <tab> peptide <tab> proteinId1 <tab> .. <tab> proteinIdM
      static TextFile preparePin_(
        const std::vector<PeptideIdentification>& peptide_ids, 
        const StringList& feature_set, 
        const std::string& enz, 
        int min_charge, 
        int max_charge);


      static Int getScanNumber_(const String& scan_identifier);

      static bool isEnz_(const char& n, const char& c, const std::string& enz);

      static Size countEnzymatic_(const String& peptide, const std::string& enz);

  };
} // namespace OpenMS
