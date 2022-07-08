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
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/ANALYSIS/NUXL/NuXLReport.h>
#include <vector>

namespace OpenMS
{

/// @brief adapted FDR calculation for NA cross-links
class OPENMS_DLLAPI NuXLFDR
{
  public:
    explicit NuXLFDR(size_t report_top_hits);

    // split by meta value "NuXL:isXL" == 0
    void splitIntoPeptidesAndXLs(const std::vector<PeptideIdentification>& peptide_ids, 
      std::vector<PeptideIdentification>& pep_pi, 
      std::vector<PeptideIdentification>& xl_pi) const;

    void mergePeptidesAndXLs(const std::vector<PeptideIdentification>& pep_pi, 
      const std::vector<PeptideIdentification>& xl_pi, 
      std::vector<PeptideIdentification>& peptide_ids) const;

    // calculate PSM-level q-values (irrespective of XL/non-XL class)
    void QValueAtPSMLevel(std::vector<PeptideIdentification>& peptide_ids) const;

    // calculate PSM-level q-values for XL and non-XL class separately.
    void calculatePeptideAndXLQValueAtPSMLevel(const std::vector<PeptideIdentification>& peptide_ids, 
      std::vector<PeptideIdentification>& pep_pi, 
      std::vector<PeptideIdentification>& xl_pi) const;

    // calculate separate FDRs, filter decoys, write PSM and protein reports
    void calculatePeptideAndXLQValueAndFilterAtPSMLevel(
      const std::vector<ProteinIdentification>& protein_ids,
      const std::vector<PeptideIdentification>& peptide_ids, 
      std::vector<PeptideIdentification>& pep,
      double peptide_PSM_qvalue_threshold,
      double peptide_peptide_qvalue_threshold,
      std::vector<PeptideIdentification>& xl_pi,
      std::vector<double> xl_PSM_qvalue_thresholds,
      std::vector<double> xl_peptidelevel_qvalue_thresholds,
      const String& out_idxml,
      int decoy_factor) const;

  private:
    size_t report_top_hits_;
};

}


