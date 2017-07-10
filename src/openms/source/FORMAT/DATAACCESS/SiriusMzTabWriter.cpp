// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2016.
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
// $Maintainer: Oliver Alka $
// $Authors: Oliver Alka, Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/DATASTRUCTURES/ListUtilsIO.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/MzTabFile.h>
#include <OpenMS/FORMAT/MzTab.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/CsvFile.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/HANDLERS/MzMLHandler.h>

#include <OpenMS/FORMAT/DATAACCESS/SiriusMzTabWriter.h>

using namespace OpenMS;
using namespace std;

MzTab SiriusMzTabWriter::store(const std::vector<String> paths, Size number)
{

  MzTab mztab;

  SiriusMzTabWriter::SiriusAdapterRun sirius_result;

  for(std::vector<String>::const_iterator it = paths.begin(); it != paths.end(); ++it)
  {

    std::string pathtosiriuscsv = *it + "/summary_sirius.csv";

    ifstream file(pathtosiriuscsv);

    //check if file is available and read input & write output
    if (file)
    {
      // read results from sirius output files
      CsvFile compounds(pathtosiriuscsv, '\t');

      // fill indentification structure containing all candidate hits for a single spectrum
      SiriusMzTabWriter::SiriusAdapterIdentification sirius_id;

      //Extract scan_index from path
      OpenMS::String path = File::path(pathtosiriuscsv);
      vector<String> substrings;
      OpenMS::String SringString;
      vector<String> newsubstrings;

      // get the digits from the end of the path (independet of their length)
      // /var/folders/mz/8dz2f2mx0fl159nrypn6p4h40000gn/T/0_out/9_2017-07-10_131258_unicorn.Informatik.Uni-Tuebingen.De_86791_1_unknown0_unknown9
      path.split('_', substrings);
      SringString = substrings[substrings.size() - 1];
      SringString.split('_', newsubstrings);
      vector<String> digit;
      OpenMS::String SringStringString = newsubstrings[newsubstrings.size() - 1];
      SringStringString.split('n', digit);
      String scan_index = digit[digit.size() - 1];

      //If there are less rows than output -> rowCount - 1 is used since the last row of the file is empty
      if (number > compounds.rowCount())
      {
        number = compounds.rowCount() - 1;
      }

      for (Size j = 1; j < number; ++j)
      {
        StringList sl;
        compounds.getRow(j, sl);
        SiriusMzTabWriter::SiriusAdapterHit sirius_hit;

        // parse single candidate hit
        sirius_hit.formula = sl[0];
        sirius_hit.adduct = sl[1];
        sirius_hit.rank = sl[2].toInt();
        sirius_hit.score = sl[3].toDouble();
        sirius_hit.treescore = sl[4].toDouble();
        sirius_hit.isoscore = sl[5].toDouble();
        sirius_hit.explainedpeaks = sl[6].toInt();
        sirius_hit.explainedintensity = sl[7].toDouble();

        sirius_id.hits.push_back(sirius_hit);
      }

      sirius_id.scan_index = scan_index;
      sirius_result.identifications.push_back(sirius_id);

      // write out results to mzTab file
      MzTabMetaData md;
      MzTabMSRunMetaData md_run;
      md_run.location = MzTabString(pathtosiriuscsv);
      md.ms_run[1] = md_run;

      MzTabSmallMoleculeSectionRows smsd;
      for (Size i = 0; i != sirius_result.identifications.size(); ++i)
      {
        const SiriusMzTabWriter::SiriusAdapterIdentification &id = sirius_result.identifications[i];
        for (Size j = 0; j != id.hits.size(); ++j)
        {
          const SiriusMzTabWriter::SiriusAdapterHit &hit = id.hits[j];
          MzTabSmallMoleculeSectionRow smsr;

          smsr.best_search_engine_score[0] = MzTabDouble(hit.score);
          smsr.best_search_engine_score[1] = MzTabDouble(hit.treescore);
          smsr.best_search_engine_score[2] = MzTabDouble(hit.isoscore);
          smsr.chemical_formula = MzTabString(hit.formula);

          MzTabOptionalColumnEntry adduct;
          adduct.first = "adduct";
          adduct.second = MzTabString(hit.adduct);

          MzTabOptionalColumnEntry rank;
          rank.first = "rank";
          rank.second = MzTabString(hit.rank);

          MzTabOptionalColumnEntry explainedPeaks;
          explainedPeaks.first = "explainedPeaks";
          explainedPeaks.second = MzTabString(hit.explainedpeaks);

          MzTabOptionalColumnEntry explainedIntensity;
          explainedIntensity.first = "explainedIntensity";
          explainedIntensity.second = MzTabString(hit.explainedintensity);

//          MzTabOptionalColumnEntry compound;
//          compound.first = "compound";
//          compound.second = MzTabString(id.scan_index);

          smsr.opt_.push_back(adduct);
          smsr.opt_.push_back(rank);
          smsr.opt_.push_back(explainedPeaks);
          smsr.opt_.push_back(explainedIntensity);
          smsd.push_back(smsr);
        }
      }

      mztab.setSmallMoleculeSectionRows(smsd);

    }
  }
  return mztab;

}

/// @endcond
