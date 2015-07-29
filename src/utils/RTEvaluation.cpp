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
// $Maintainer: Nico Pfeifer $
// $Authors: Nico Pfeifer $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FILTERING/ID/IDFilter.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CHEMISTRY/EnzymaticDigestion.h>

#include <map>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
    @page UTILS_RTEvaluation RTEvaluation

    @brief Application that evaluates TPs (true positives), TNs, FPs, and FNs for an idXML file with predicted RTs.

    The method needs an idXML file with IDs and predicted RTs. The second input file is a file containing
    the protein sequences which are considered as positive hits. This tool then evaluates the tps, fps, tns,
    and fns for the unfiltered IDs, for the IDs filtered in first RT dimension, for the IDs filtered in
    the second RT dimension as well as for the IDs filtered in both dimensions. The output is a table with
    either csv format (can be imported by Excel) or latex format (to include in tables in your latex manuscripts).

    @note Currently mzIdentML (mzid) is not directly supported as an input/output format of this tool. Convert mzid files to/from idXML using @ref TOPP_IDFileConverter if necessary.

    <B>The command line parameters of this tool are:</B>
    @verbinclude UTILS_RTEvaluation.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude UTILS_RTEvaluation.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPIdXMLEvaluation :
  public TOPPBase
{
public:
  TOPPIdXMLEvaluation() :
    TOPPBase("RTEvaluation", "Application that evaluates TPs (true positives), TNs, FPs, and FNs for an idXML file with predicted RTs.", false)
  {

  }

  enum State {TP, FP, TN, FN, NE};

protected:
  void registerOptionsAndFlags_()
  {
    registerInputFile_("in", "<file>", "", "input file");
    setValidFormats_("in", ListUtils::create<String>("idXML"));
    registerOutputFile_("out", "<file>", "", "output file ");
    setValidFormats_("out", ListUtils::create<String>("csv"));
    registerInputFile_("sequences_file", "<file>", "", "Filename of a FASTA file containing protein sequences.\n"
                                                       "All peptides that are not a substring of a sequence in this file are considered as false", false);
    setValidFormats_("sequences_file", ListUtils::create<String>("fasta"));
    registerFlag_("latex", "indicates whether the output file format of the table should be latex or csv");
    registerDoubleOption_("p_value_dim_1", "<float>", 0.01, "Significance level of first dimension RT filter", false);
    setMinFloat_("p_value_dim_1", 0);
    setMaxFloat_("p_value_dim_1", 1);
    registerDoubleOption_("p_value_dim_2", "<float>", 0.05, "Significance level of second dimension RT filter", false);
    setMinFloat_("p_value_dim_2", 0);
    setMaxFloat_("p_value_dim_2", 1);
  }

  ExitCodes main_(int, const char**)
  {
    IdXMLFile idXML_file;
    vector<ProteinIdentification> protein_identifications;
    vector<PeptideIdentification> identifications;
    vector<FASTAFile::FASTAEntry> sequences;

    bool latex = getFlag_("latex");
    bool strict = true;
    PeptideIdentification filtered_identification;
    PeptideIdentification filtered_identification_rt1;
    PeptideIdentification filtered_identification_rt2;
    PeptideIdentification filtered_identification_both;
    IDFilter filter;
    bool no_protein_identifiers = true;
    double p_value_dim_1 = getDoubleOption_("p_value_dim_1");
    double p_value_dim_2 = getDoubleOption_("p_value_dim_2");
    State state = TP;
    State state_rt1 = TP;
    State state_rt2 = TP;
    vector<double> fdrs;
    fdrs.push_back(0.);
    fdrs.push_back(0.01);
    fdrs.push_back(0.02);
    fdrs.push_back(0.03);
    fdrs.push_back(0.04);
    fdrs.push_back(0.05);
    fdrs.push_back(0.1);
    fdrs.push_back(0.15);
    fdrs.push_back(0.2);
    fdrs.push_back(0.25);
    fdrs.push_back(0.3);
    fdrs.push_back(0.35);
    fdrs.push_back(0.4);
    fdrs.push_back(0.45);
    fdrs.push_back(0.5);
    vector<Size> temp_performances;
    // ofstream tempfile("test.txt");
    vector<vector<Size> > performances;


    protein_identifications.push_back(ProteinIdentification());
    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------
    String inputfile_name = getStringOption_("in");
    String sequences_file_name = getStringOption_("sequences_file");
    String outputfile_name = getStringOption_("out");

    //-------------------------------------------------------------
    // reading input
    //-------------------------------------------------------------
    String document_id;
    idXML_file.load(inputfile_name, protein_identifications, identifications, document_id);
    if (sequences_file_name != "")
    {
      FASTAFile().load(sequences_file_name, sequences);
    }
    //-------------------------------------------------------------
    // calculations
    //-------------------------------------------------------------
    for (SignedSize j = fdrs.size() - 1; j >= 0; --j)
    {
      Size tps = 0;
      Size fps = 0;
      Size nes = 0;
      Size tps_rt1 = 0;
      Size fps_rt1 = 0;
      Size tns_rt1 = 0;
      Size fns_rt1 = 0;
      Size nes_rt1 = 0;
      Size tps_rt2 = 0;
      Size fps_rt2 = 0;
      Size tns_rt2 = 0;
      Size fns_rt2 = 0;
      Size nes_rt2 = 0;
      Size tps_both = 0;
      Size fps_both = 0;
      Size tns_both = 0;
      Size fns_both = 0;
      Size nes_both = 0;

      temp_performances.clear();
      temp_performances.resize(20, 0);
      for (Size i = 0; i < identifications.size(); i++)
      {
        PeptideIdentification temp_identification_2 = identifications[i];
        PeptideIdentification temp_identification_3;
        filter.filterIdentificationsByScore(temp_identification_2, fdrs[j], temp_identification_3);

        filtered_identification = temp_identification_3;
        filtered_identification_rt1 = temp_identification_3;
        filtered_identification_rt2 = temp_identification_3;

        if (p_value_dim_1 > 0)
        {
          PeptideIdentification temp_identification = filtered_identification_rt1;
          filter.filterIdentificationsByRTFirstDimPValues(temp_identification, filtered_identification_rt1, p_value_dim_1);
        }

        if (p_value_dim_2 > 0)
        {
          PeptideIdentification temp_identification = filtered_identification_rt2;
          filter.filterIdentificationsByRTPValues(temp_identification, filtered_identification_rt2, p_value_dim_2);
        }
        if (p_value_dim_1 > 0 && p_value_dim_2 > 0)
        {
          PeptideIdentification temp_identification = filtered_identification_rt1;
          filter.filterIdentificationsByRTPValues(temp_identification, filtered_identification_both, p_value_dim_2);
        }

        if (strict)
        {
          PeptideIdentification temp_identification = filtered_identification;
          filter.filterIdentificationsByBestHits(temp_identification, filtered_identification, strict);
          temp_identification = filtered_identification_rt1;
          filter.filterIdentificationsByBestHits(temp_identification, filtered_identification_rt1, strict);
          temp_identification = filtered_identification_rt2;
          filter.filterIdentificationsByBestHits(temp_identification, filtered_identification_rt2, strict);
          temp_identification = filtered_identification_both;
          filter.filterIdentificationsByBestHits(temp_identification, filtered_identification_both, strict);
        }
        if (!filtered_identification.getHits().empty())
        {
          if (sequences_file_name != "")
          {
            PeptideIdentification temp_identification = filtered_identification;
            filter.filterIdentificationsByProteins(temp_identification, sequences, filtered_identification, no_protein_identifiers);
          }
          if (filtered_identification.getHits().empty())
          {
            ++fps;
            state = FP;
          }
          else
          {
            ++tps;
            state = TP;
          }

        }
        else
        {
          ++nes;
          state = NE;
        }

        if (!filtered_identification_rt1.getHits().empty())
        {
          if (sequences_file_name != "")
          {
            PeptideIdentification temp_identification = filtered_identification_rt1;
            filter.filterIdentificationsByProteins(temp_identification, sequences, filtered_identification_rt1, no_protein_identifiers);
          }
          if (filtered_identification_rt1.getHits().empty())
          {
            ++fps_rt1;
            state_rt1 = FP;
          }
          else
          {
            ++tps_rt1;
            state_rt1 = TP;
            // tempfile <<     filtered_identification.getHits()[0].getSequence() << " " << filtered_identification.getMetaValue("MZ") << " " << filtered_identification.getMetaValue("RT") << endl;
          }
        }
        else if (state == FP)
        {
          ++tns_rt1;
          state_rt1 = TN;
        }
        else if (state == TP)
        {
          ++fns_rt1;
          state_rt1 = FN;
        }
        else
        {
          ++nes_rt1;
          state_rt1 = NE;
        }

        if (!filtered_identification_rt2.getHits().empty())
        {
          if (sequences_file_name != "")
          {
            PeptideIdentification temp_identification = filtered_identification_rt2;
            filter.filterIdentificationsByProteins(temp_identification, sequences, filtered_identification_rt2, no_protein_identifiers);
          }
          if (filtered_identification_rt2.getHits().empty())
          {
            ++fps_rt2;
            state_rt2 = FP;
          }
          else
          {
            ++tps_rt2;
            state_rt2 = TP;
          }
        }
        else if (state == FP)
        {
          ++tns_rt2;
          state_rt2 = TN;
        }
        else if (state == TP)
        {
          ++fns_rt2;
          state_rt2 = FN;
        }
        else
        {
          ++nes_rt2;
          state_rt2 = NE;
        }

        if (!filtered_identification_both.getHits().empty())
        {
          if (sequences_file_name != "")
          {
            PeptideIdentification temp_identification = filtered_identification_both;
            filter.filterIdentificationsByProteins(temp_identification, sequences, filtered_identification_both, no_protein_identifiers);
          }
        }
        if (state_rt1 == TP && state_rt2 == TP)
        {
          ++tps_both;
        }
        else if ((state_rt1 == TP || state_rt2 == TP) && (state_rt1 == NE || state_rt2 == NE))
        {
          ++tps_both;
        }
        else if (state_rt1 == FP && state_rt2 == FP)
        {
          ++fps_both;
        }
        else if ((state_rt1 == TN || state_rt2 == TN || state_rt1 == NE || state_rt2 == NE)
                && state == FP)
        {
          ++tns_both;
        }
        else if ((state_rt1 == FN || state_rt2 == FN || state_rt1 == NE || state_rt2 == NE) && state == TP)
        {
          ++fns_both;
        }
        else if ((state_rt1 == NE || state_rt2 == NE) && state == NE)
        {
          ++nes_both;
        }
        else if (((state_rt1 == TP && state_rt2 == FP) || (state_rt1 == FP && state_rt2 == TP))
                && (!filtered_identification_both.getHits().empty()))
        {
          ++tps_both;
        }
        else
        {
          cout << "RT1 is in state: " << state_rt1 << " and RT2 is in state: "
               << state_rt2 << endl;
        }
      }
      cout << "q-value threshold: " << fdrs[j] << " ***************" << endl;
      cout << "Unfiltered:: True positives: " << tps << " false positives: "
           << fps << " not evaluated: " << nes
           << " total: " << (tps + fps + nes) << endl;
      cout << "Filtered RT1:: TPss: " << tps_rt1 << " FPs: "
           << fps_rt1 << " TNs: " << tns_rt1 << " FNs: " << fns_rt1
           << " not evaluated: "   << nes_rt1
           << " total: " << (tps_rt1 + fps_rt1 + tns_rt1 + fns_rt1 + nes_rt1)
           << endl;
      cout << "Filtered RT2:: TPss: " << tps_rt2 << " FPs: "
           << fps_rt2 << " TNs: " << tns_rt2 << " FNs: " << fns_rt2
           << " not evaluated: " << nes_rt2
           << " total: " << (tps_rt2 + fps_rt2 + tns_rt2 + fns_rt2 + nes_rt2)
           << endl;
      cout << "Filtered both dimensions:: TPss: " << tps_both << " FPs: "
           << fps_both << " TNs: " << tns_both << " FNs: " << fns_both
           << " not evaluated: " << nes_both
           << " total: " << (tps_both + fps_both + tns_both + fns_both + nes_both)
           << endl;

      temp_performances[0] = tps;
      temp_performances[1] = fps;
      temp_performances[2] = 0;
      temp_performances[3] = 0;
      temp_performances[4] = nes;
      temp_performances[5] = tps_rt1;
      temp_performances[6] = fps_rt1;
      temp_performances[7] = tns_rt1;
      temp_performances[8] = fns_rt1;
      temp_performances[9] = nes_rt1;
      temp_performances[10] = tps_rt2;
      temp_performances[11] = fps_rt2;
      temp_performances[12] = tns_rt2;
      temp_performances[13] = fns_rt2;
      temp_performances[14] = nes_rt2;
      temp_performances[15] = tps_both;
      temp_performances[16] = fps_both;
      temp_performances[17] = tns_both;
      temp_performances[18] = fns_both;
      temp_performances[19] = nes_both;
      performances.push_back(temp_performances);
    }

    ofstream output_file(outputfile_name.c_str());
    if (latex)
    {
      output_file << "q-value_threshold & tp & fp & tn & fn & precision & tp & fp & tn & fn & precision & tp & fp & tn & fn & precision & tp & fp & tn & fn & precision" << endl;
    }
    else
    {
      output_file << "q-value_threshold ; tp ; fp ; tn ; fn ; precision ; tp ; fp ; tn ; fn ; precision ; tp ; fp ; tn ; fn ; precision ; tp ; fp ; tn ; fn ; precision" << endl;
    }

    for (SignedSize i = performances.size() - 1; i >= 0; --i)
    {
      output_file << fdrs[performances.size() - i - 1];
      for (Size j = 0; j < performances[i].size(); ++j)
      {
        if (latex)
        {
          output_file << " &";
        }
        else
        {
          output_file << " ;";
        }


        if (j % 5 == 4)
        {
          output_file << " " << (performances[i][j - 4] / ((double) performances[i][j - 4] + performances[i][j - 3]));
        }
        else
        {
          output_file << " " << performances[i][j];
        }
      }
      output_file << endl;
    }
    output_file << flush;
    output_file.close();

    return EXECUTION_OK;
  }

};


int main(int argc, const char** argv)
{
  TOPPIdXMLEvaluation tool;
  return tool.main(argc, argv);
}

/// @endcond
