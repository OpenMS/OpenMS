// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
// $Maintainer: Stephan Aiche $
// $Authors: Stephan Aiche $
// --------------------------------------------------------------------------
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/ANALYSIS/QUANTITATION/ItraqChannelExtractor.h>
#include <OpenMS/ANALYSIS/QUANTITATION/ItraqQuantifier.h>
#include <OpenMS/ANALYSIS/QUANTITATION/ItraqConstants.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/FORMAT/MzQuantMLFile.h>
#include <OpenMS/METADATA/MSQuantifications.h>

using namespace OpenMS;
using namespace std;

// NOTE: This tool is a copy if ITRAQAnalyzer with the difference that it sets, the mode
//       automatically to TMT_SIXPLEX and ignores all iTRAQ-modes.

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
  @page TOPP_TMTAnalyzer TMTAnalyzer

  @brief Extracts and normalizes TMT information from an MS experiment.

<CENTER>
    <table>
        <tr>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
      <td VALIGN="middle" ROWSPAN=3> \f$ \longrightarrow \f$ TMTAnalyzer \f$ \longrightarrow \f$</td>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_FileConverter </td>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=2> @ref TOPP_IDMapper</td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_FileFilter </td>
        </tr>
    </table>
</CENTER>

  Extract the TMT reporter ion intensities (6plex) from raw MS2 data, does isotope corrections and
  stores the resulting quantitation as consensusXML, where each consensus centroid corresponds to
  one TMT MS2 scan (e.g., CID). The position of the centroid is the precursor position, its
  sub-elements are the channels (thus having m/z's of 126-131).

  Isotope correction is done using non-negative least squares (NNLS). See @ref TOPP_ITRAQAnalyzer
  for details.

  The correction matrices can be found (and changed) in the INI file. However, these matrices for
  TMT are now stable, and every kit delivered should have the same isotope
  correction values. Thus, there should be no need to change them, but feel free to compare the
  values in the INI file with your kit's Certificate.

  After this quantitation step, you might want to annotate the consensus elements with the
  respective identifications, obtained from an identification pipeline.
  Note that quantification is solely on peptide level at this stage. In order to obtain protein
  quantifications, you can try @ref TOPP_TextExporter to obtain a simple text format which you can
  feed to other software tools (e.g., R), or you can try @ref TOPP_ProteinQuantifier.

    <B>The command line parameters of this tool are:</B>
  @verbinclude TOPP_TMTAnalyzer.cli
    <B>INI file documentation of this tool:</B>
  @htmlinclude TOPP_TMTAnalyzer.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPTMTAnalyzer :
  public TOPPBase
{
public:
  TOPPTMTAnalyzer() :
    TOPPBase("TMTAnalyzer", "Calculates TMT quantitative values for peptides", true, true)
  {
  }

protected:
  void registerOptionsAndFlags_()
  {
    registerInputFile_("in", "<file>", "", "input raw/picked data file ");
    setValidFormats_("in", StringList::create("mzML"));
    registerOutputFile_("out", "<file>", "", "output consensusXML file with quantitative information");
    setValidFormats_("out", StringList::create("consensusXML"));

    registerOutputFile_("out_mzq", "<file>", "", "Optional output file of MzQuantML.", false, true);
    setValidFormats_("out_mzq", StringList::create("mzq"));

    registerOutputFile_("out_stats", "<file>", "", "output statistics as tab-separated file (readable by R or Excel or ...)", false);
    setValidFormats_("out_stats", StringList::create("tsv"));

    addEmptyLine_();

    registerSubsection_("algorithm", "Algorithm parameters section");
  }

  Param getSubsectionDefaults_(const String & /*section*/) const
  {
    Param tmp;
    tmp.insert("Extraction:", ItraqChannelExtractor(ItraqQuantifier::TMT_SIXPLEX).getParameters());
    tmp.insert("Quantification:", ItraqQuantifier(ItraqQuantifier::TMT_SIXPLEX).getParameters());
    tmp.setValue("MetaInformation:Program", "OpenMS::TMTAnalyzer", "", StringList::create("advanced"));
    return tmp;
  }

  ExitCodes main_(int, const char **)
  {
    //-------------------------------------------------------------
    // parameter handling
    //-------------------------------------------------------------
    String in = getStringOption_("in");
    String out = getStringOption_("out");
    String out_stats = getStringOption_("out_stats");
    String out_mzq = getStringOption_("out_mzq");

    Int itraq_type = ItraqQuantifier::TMT_SIXPLEX;
    //-------------------------------------------------------------
    // loading input
    //-------------------------------------------------------------

    MzMLFile mz_data_file;
    MSExperiment<Peak1D> exp;
    mz_data_file.setLogType(log_type_);
    mz_data_file.load(in, exp);

    //-------------------------------------------------------------
    // calculations
    //-------------------------------------------------------------
    Param extract_param(getParam_().copy("algorithm:Extraction:", true));
    ItraqChannelExtractor itraq_ce(itraq_type, extract_param);

    ConsensusMap consensus_map_raw, consensus_map_quant;
    // extract raw signals
    itraq_ce.run(exp, consensus_map_raw);

    // do normalization
    Param quant_param(getParam_().copy("algorithm:Quantification:", true));
    ItraqQuantifier itraq_quant(itraq_type, quant_param);

    itraq_quant.run(consensus_map_raw, consensus_map_quant);

    // assign unique ID to output file (this might throw an exception.. but thats ok, as we want the program to quit then)
    if (getStringOption_("id_pool").trim().length() > 0) getDocumentIDTagger_().tag(consensus_map_quant);

    // annotate output file with MetaInformation
    Param metainfo_param(getParam_().copy("algorithm:MetaInformation:", true));
    for (Param::ParamIterator it = metainfo_param.begin(); it != metainfo_param.end(); ++it)
    {
      consensus_map_quant.setMetaValue(it->name, it->value);
    }


    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------

    //annotate output with data processing info
    addDataProcessing_(consensus_map_quant, getProcessingInfo_(DataProcessing::QUANTITATION));

    // add filename references
    for (ConsensusMap::FileDescriptions::iterator it = consensus_map_quant.getFileDescriptions().begin();
         it != consensus_map_quant.getFileDescriptions().end();
         ++it)
    {
      it->second.filename = in;
    }

    ConsensusXMLFile cm_file;
    cm_file.store(out, consensus_map_quant);

    if (!out_mzq.trim().empty())
    {
      MSQuantifications msq;
      std::vector<std::vector<std::pair<String, DoubleReal> > > labels;
      for (Size i = 0; i < ItraqConstants::CHANNEL_COUNT[ItraqConstants::TMT_SIXPLEX]; ++i)
      {
        std::vector<std::pair<String, DoubleReal> > one_label;
        one_label.push_back(std::make_pair<String, DoubleReal>(String("Channel ") + String(ItraqConstants::CHANNELS_TMT_SIXPLEX[i][0]), DoubleReal(ItraqConstants::CHANNELS_FOURPLEX[i][0])));
        labels.push_back(one_label);
      }
      msq.registerExperiment(exp, labels);       //add assays
      msq.assignUIDs();
      MSQuantifications::QUANT_TYPES quant_type = MSQuantifications::MS2LABEL;
      msq.setAnalysisSummaryQuantType(quant_type);      //add analysis_summary_

      msq.addConsensusMap(consensus_map_quant);      //add TMTAnalyzer result
      //~ add AuditCollection - no such concept in TOPPTools yet
      MzQuantMLFile file;
      file.store(out_mzq, msq);
    }

    std::cout << itraq_quant.getStats();
    if (!out_stats.trim().empty())
    {
      ofstream f;
      f.open(out_stats.c_str(), ios_base::out);
      f << itraq_quant.getStats();
      f.close();
    }

    return EXECUTION_OK;
  }

};

int main(int argc, const char ** argv)
{
  TOPPTMTAnalyzer tool;
  return tool.main(argc, argv);
}

/// @endcond
