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
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <OpenMS/KERNEL/MSExperiment.h>

#include <OpenMS/ANALYSIS/QUANTITATION/ItraqChannelExtractor.h>
#include <OpenMS/ANALYSIS/QUANTITATION/ItraqQuantifier.h>
#include <OpenMS/ANALYSIS/QUANTITATION/ItraqConstants.h>

// the available quantitation methods
#include <OpenMS/ANALYSIS/QUANTITATION/IsobaricQuantitationMethod.h>
#include <OpenMS/ANALYSIS/QUANTITATION/ItraqFourPlexQuantitationMethod.h>
#include <OpenMS/ANALYSIS/QUANTITATION/ItraqEightPlexQuantitationMethod.h>
#include <OpenMS/ANALYSIS/QUANTITATION/TMTSixPlexQuantitationMethod.h>
#include <OpenMS/ANALYSIS/QUANTITATION/TMTTenPlexQuantitationMethod.h>

#include <OpenMS/ANALYSIS/QUANTITATION/IsobaricChannelExtractor.h>
#include <OpenMS/ANALYSIS/QUANTITATION/IsobaricQuantifier.h>

#include <OpenMS/SYSTEM/File.h>

#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/MzQuantMLFile.h>

#include <OpenMS/METADATA/MSQuantifications.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
    @page TOPP_IsobaricAnalyzer IsobaricAnalyzer

    @brief Extracts and normalizes isobaric labeling information from an LC-MS/MS experiment.

<CENTER>
    <table>
        <tr>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
            <td VALIGN="middle" ROWSPAN=3> \f$ \longrightarrow \f$ IsobaricAnalyzer \f$ \longrightarrow \f$</td>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_PeakPickerHiRes </td>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=2> @ref TOPP_IDMapper</td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_FileFilter </td>
        </tr>
    </table>
</CENTER>

  This tool currently supports iTRAQ 4-plex and 8-plex, and TMT 6-plex and 10-plex as labeling methods.
  It extracts the isobaric reporter ion intensities from centroided MS2 data, performs isotope correction and stores the resulting quantitation in a consensus map, in which each consensus feature represents one relevant MS2 scan (e.g. HCD; see parameters @p select_activation and @p min_precursor_intensity).
  The position (RT, m/z) of the consensus centroid is the precursor position; the sub-elements correspond to the channels (with m/z values of 113-121 for iTRAQ and 126-131 for TMT, respectively).
  
  @note If none of the reporter ions can be detected in an MS2 scan, a consensus feature will still be generated, but the intensities of the overall feature and of all its sub-elements will be zero. (If desired, such features can be removed by applying an intensity filter in @ref TOPP_FileFilter.)

  The input MS2 spectra have to be in centroid mode for the tool to work properly. Use e.g. @ref TOPP_PeakPickerHiRes to perform centroiding of profile data, if necessary.
  
  Isotope correction is done using non-negative least squares (NNLS), i.e.:@n
  Minimize ||Ax - b||, subject to x >= 0, where b is the vector of observed reporter intensities (with "contaminating" isotope species), A is a correction matrix (as supplied by the manufacturer of the labeling kit) and x is the desired vector of corrected (real) reporter intensities.
  Other software tools solve this problem by using an inverse matrix multiplication, but this can yield entries in x which are negative. In a real sample, this solution cannot possibly be true, so usually negative values (= negative reporter intensities) are set to zero.
  However, a negative result usually means that noise was not properly accounted for in the calculation. We thus use NNLS to get a non-negative solution, without the need to truncate negative values. In the (usual) case that inverse matrix multiplication yields only positive values, our NNLS will give the exact same optimal solution.

  The correction matrices can be found (and changed) in the INI file (parameter @p correction_matrix of the corresponding labeling method). However, these matrices for both 4-plex and 8-plex iTRAQ are now stable, and every kit delivered should have the same isotope correction values. Thus, there should be no need to change them, but feel free to compare the values in the INI file with your kit's certificate. For TMT (6-plex and 10-plex) the values have to be adapted for each kit.

  After the quantitation, you may want to annotate the consensus features with corresponding peptide identifications, obtained from an identification pipeline. Use @ref TOPP_IDMapper to perform the annotation, but make sure to set suitably small RT and m/z tolerances for the mapping, since the identifications will come from the very same MS2 scans that are now represented by consensus features. In general it should be possible to achieve a perfect one-to-one matching of every identification to a single consensus feature.@n
  Note that quantification will be solely on peptide level after this stage. In order to obtain protein quantities, you can use @ref TOPP_TextExporter to obtain a simple text format which you can feed to other software tools (e.g., R), or you can apply @ref TOPP_ProteinQuantifier.

    <B>The command line parameters of this tool are:</B>
    @verbinclude TOPP_IsobaricAnalyzer.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude TOPP_IsobaricAnalyzer.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPIsobaricAnalyzer :
  public TOPPBase
{
private:
  std::map<String, IsobaricQuantitationMethod*> quant_methods_;
  std::map<String, String> quant_method_names_;

public:
  TOPPIsobaricAnalyzer() :
    TOPPBase("IsobaricAnalyzer", "Calculates isobaric quantitative values for peptides", true, true)
  {
    ItraqFourPlexQuantitationMethod* itraq4plex = new ItraqFourPlexQuantitationMethod();
    ItraqEightPlexQuantitationMethod* itraq8plex = new ItraqEightPlexQuantitationMethod();
    TMTSixPlexQuantitationMethod* tmt6plex = new TMTSixPlexQuantitationMethod();
    TMTTenPlexQuantitationMethod* tmt10plex = new TMTTenPlexQuantitationMethod();
    quant_methods_[itraq4plex->getName()] = itraq4plex;
    quant_methods_[itraq8plex->getName()] = itraq8plex;
    quant_methods_[tmt6plex->getName()] = tmt6plex;
    quant_methods_[tmt10plex->getName()] = tmt10plex;
    quant_method_names_[itraq4plex->getName()] = "iTRAQ 4-plex";
    quant_method_names_[itraq8plex->getName()] = "iTRAQ 8-plex";
    quant_method_names_[tmt6plex->getName()] = "TMT 6-plex";
    quant_method_names_[tmt10plex->getName()] = "TMT 10-plex";
  }

  ~TOPPIsobaricAnalyzer()
  {
    // free allocated labelers
    for (std::map<String, IsobaricQuantitationMethod*>::iterator it = quant_methods_.begin();
         it != quant_methods_.end();
         ++it)
    {
      delete it->second;
    }
  }

protected:
  void registerOptionsAndFlags_()
  {
    // initialize with the first available type
    registerStringOption_("type", "<mode>", quant_methods_.begin()->first, "Isobaric Quantitation method used in the experiment.", false);
    StringList valid_types;
    for (std::map<String, IsobaricQuantitationMethod*>::iterator it = quant_methods_.begin();
         it != quant_methods_.end();
         ++it)
    {
      valid_types.push_back(it->first);
    }
    setValidStrings_("type", valid_types);

    registerInputFile_("in", "<file>", "", "input raw/picked data file ");
    setValidFormats_("in", ListUtils::create<String>("mzML"));
    registerOutputFile_("out", "<file>", "", "output consensusXML file with quantitative information");
    setValidFormats_("out", ListUtils::create<String>("consensusXML"));

    registerSubsection_("extraction", "Parameters for the channel extraction.");
    registerSubsection_("quantification", "Parameters for the peptide quantification.");
    for (std::map<String, IsobaricQuantitationMethod*>::iterator it = quant_methods_.begin();
         it != quant_methods_.end();
         ++it)
    {
      registerSubsection_(it->second->getName(), String("Algorithm parameters for ") + quant_method_names_[it->second->getName()]);
    }
  }

  Param getSubsectionDefaults_(const String& section) const
  {
    ItraqFourPlexQuantitationMethod temp_quant;
    if (section == "extraction")
    {
      return IsobaricChannelExtractor(&temp_quant).getParameters();
    }
    else if (section == "quantification")
    {
      return IsobaricQuantifier(&temp_quant).getParameters();
    }
    else
    {
      std::map<String, IsobaricQuantitationMethod*>::const_iterator it = quant_methods_.find(section);
      if (it != quant_methods_.end())
      {
        return it->second->getParameters();
      }
      else
      {
        // this should not happen
        Param empty;
        return empty;
      }
    }
  }

  ExitCodes main_(int, const char**)
  {
    //-------------------------------------------------------------
    // parameter handling
    //-------------------------------------------------------------
    String in = getStringOption_("in");
    String out = getStringOption_("out");

    //-------------------------------------------------------------
    // loading input
    //-------------------------------------------------------------

    MzMLFile mz_data_file;
    MSExperiment<Peak1D> exp;
    mz_data_file.setLogType(log_type_);
    mz_data_file.load(in, exp);

    //-------------------------------------------------------------
    // init quant method
    //-------------------------------------------------------------
    IsobaricQuantitationMethod* quant_method = quant_methods_[getStringOption_("type")];

    // set the parameters for this method
    quant_method->setParameters(getParam_().copy(quant_method->getName() + ":", true));

    //-------------------------------------------------------------
    // calculations
    //-------------------------------------------------------------
    Param extract_param(getParam_().copy("extraction:", true));
    IsobaricChannelExtractor channel_extractor(quant_method);
    channel_extractor.setParameters(extract_param);

    ConsensusMap consensus_map_raw, consensus_map_quant;

    // extract channel information
    channel_extractor.extractChannels(exp, consensus_map_raw);

    IsobaricQuantifier quantifier(quant_method);
    Param quant_param(getParam_().copy("quantification:", true));
    quantifier.setParameters(quant_param);

    quantifier.quantify(consensus_map_raw, consensus_map_quant);

    // assign unique ID to output file (this might throw an exception... but that's ok, as we want the program to quit then)
    if (getStringOption_("id_pool").trim().length() > 0) getDocumentIDTagger_().tag(consensus_map_quant);

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

    consensus_map_quant.ensureUniqueId();
    consensus_map_quant.setPrimaryMSRunPath(exp.getPrimaryMSRunPath());
    ConsensusXMLFile().store(out, consensus_map_quant);

    return EXECUTION_OK;
  }

};

int main(int argc, const char** argv)
{
  TOPPIsobaricAnalyzer tool;
  return tool.main(argc, argv);
}

/// @endcond
