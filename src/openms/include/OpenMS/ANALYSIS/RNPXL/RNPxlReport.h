// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2015.
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

#ifndef OPENMS_ANALYSIS_RNPXL_RNPXLREPORT
#define OPENMS_ANALYSIS_RNPXL_RNPXLREPORT

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/ANALYSIS/RNPXL/RNPxlMarkerIonExtractor.h>
#include <OpenMS/DATASTRUCTURES/ListUtilsIO.h>

namespace OpenMS
{

struct OPENMS_DLLAPI RNPxlReportRow
{
  bool no_id;
  double rt;
  double original_mz;
  String accessions;
  String RNA;
  String peptide;
  double best_localization_score;
  String localization_scores;
  String best_localization;
  Int charge;
  double score;
  double peptide_weight;
  double RNA_weight;
  double xl_weight;
  double abs_prec_error;
  double rel_prec_error;
  RNPxlMarkerIonExtractor::MarkerIonsType marker_ions;
  double m_H;
  double m_2H;
  double m_3H;
  double m_4H;

  String getString(String separator) const;

};

struct OPENMS_DLLAPI RNPxlReportRowHeader
{
  String getString(String separator)
  {
    StringList sl;
    sl << "#RT" << "original m/z" << "proteins" << "RNA" << "peptide" << "charge" << "score"
       << "best localization score" << "localization scores" << "best localization(s)"
       << "peptide weight" << "RNA weight" << "cross-link weight";

    // marker ion fields
    RNPxlMarkerIonExtractor::MarkerIonsType marker_ions = RNPxlMarkerIonExtractor::extractMarkerIons(PeakSpectrum(), 0.0); // call only to generate header entries
    for (RNPxlMarkerIonExtractor::MarkerIonsType::const_iterator it = marker_ions.begin(); it != marker_ions.end(); ++it)
    {
      for (Size i = 0; i != it->second.size(); ++i)
      {
        sl << String(it->first + "_" + it->second[i].first);
      }
    }
    sl << "abs prec. error Da" << "rel. prec. error ppm" << "M+H" << "M+2H" << "M+3H" << "M+4H";
    return ListUtils::concatenate(sl, separator);
  }
};

// create report
struct OPENMS_DLLAPI RNPxlReport
{
  static std::vector<RNPxlReportRow> annotate(const PeakMap& spectra, std::vector<PeptideIdentification>& peptide_ids, double marker_ions_tolerance);
};

}

#endif

