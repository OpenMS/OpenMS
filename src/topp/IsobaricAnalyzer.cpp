// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>

// the available quantitation methods
#include <OpenMS/ANALYSIS/QUANTITATION/IsobaricQuantitationMethod.h>
#include <OpenMS/ANALYSIS/QUANTITATION/ItraqFourPlexQuantitationMethod.h>
#include <OpenMS/ANALYSIS/QUANTITATION/ItraqEightPlexQuantitationMethod.h>
#include <OpenMS/ANALYSIS/QUANTITATION/TMTSixPlexQuantitationMethod.h>
#include <OpenMS/ANALYSIS/QUANTITATION/TMTTenPlexQuantitationMethod.h>
#include <OpenMS/ANALYSIS/QUANTITATION/TMTElevenPlexQuantitationMethod.h>
#include <OpenMS/ANALYSIS/QUANTITATION/TMTSixteenPlexQuantitationMethod.h>
#include <OpenMS/ANALYSIS/QUANTITATION/TMTEighteenPlexQuantitationMethod.h>
#include <OpenMS/ANALYSIS/QUANTITATION/IsobaricChannelExtractor.h>
#include <OpenMS/ANALYSIS/QUANTITATION/IsobaricQuantifier.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <memory> // for std::unique_ptr

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
            <th ALIGN = "center"> pot. predecessor tools </td>
            <td VALIGN="middle" ROWSPAN=3> &rarr; IsobaricAnalyzer &rarr;</td>
            <th ALIGN = "center"> pot. successor tools </td>
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

The input MSn spectra have to be in centroid mode for the tool to work properly. Use e.g. @ref TOPP_PeakPickerHiRes to perform centroiding of profile data, if necessary.

This tool currently supports iTRAQ 4-plex and 8-plex, and TMT 6-plex, 10-plex, 11-plex, 16-plex, and 18-plex as labeling methods.
It extracts the isobaric reporter ion intensities from centroided MS2 or MS3 data (MSn), then performs isotope correction and stores the resulting quantitation in a consensus map,
in which each consensus feature represents one relevant MSn scan (e.g. HCD; see parameters @p select_activation and @p min_precursor_intensity).
The MS level for quantification is chosen automatically, i.e. if MS3 is present, MS2 will be ignored.
For intensity, the closest non-zero m/z signal to the theoretical position is taken as reporter ion abundance.
The position (RT, m/z) of the consensus centroid is the precursor position in MS1 (from the MS2 spectrum);
the consensus sub-elements correspond to the theoretical channel m/z (with m/z values of 113-121 Th for iTRAQ and 126-131 Th for TMT, respectively).

For all labeling techniques, the search radius (@p reporter_mass_shift) should be set as small as possible, to avoid picking up false-positive ions as reporters.
Usually, Orbitraps deliver precision of about 0.0001 Th at this low mass range. Low intensity reporters might have a slightly higher deviation.
By default, the mass range is set to ~0.002 Th, which should be sufficient for all instruments (~15 ppm).
The tool will throw an Exception if you set it below 0.0001 Th (~0.7ppm).
The tool will also throw an Exception if you set @p reporter_mass_shift > 0.003 Th for TMT-10plex and TMT-11plex, since this could
lead to ambiguities with neighbouring channels (which are ~0.006 Th apart in most cases).

For quality control purposes, the tool reports the median distance between the theoretical vs. observed reporter ion peaks in each channel.
The search radius is fixed to 0.5 Th (regardless of the user defined search radius). This allows to track calibration issues.
For TMT-10plex, these results are automatically omitted if they could be confused with a neighbouring channel, i.e.
exceed the tolerance to a neighbouring channel with the same nominal mass (C/N channels).
If the distance is too large, you might have a m/z calibration problem (see @ref TOPP_InternalCalibration).

@note If none of the reporter ions can be detected in an MSn scan, a consensus feature will still be generated, 
but the intensities of the overall feature and of all its sub-elements will be zero.
(If desired, such features can be removed by applying an intensity filter in @ref TOPP_FileFilter.)
However, if the spectrum is completely empty (no ions whatsoever), no consensus feature will be generated.

Isotope correction is done using non-negative least squares (NNLS), i.e.:@n
Minimize ||Ax - b||, subject to x >= 0, where b is the vector of observed reporter intensities (with "contaminating" isotope species), 
A is a correction matrix (as supplied by the manufacturer of the labeling kit) and x is the desired vector of corrected (real) reporter intensities.
Other software tools solve this problem by using an inverse matrix multiplication, but this can yield entries in x which are negative. 
In a real sample, this solution cannot possibly be true, so usually negative values (= negative reporter intensities) are set to zero.
However, a negative result usually means that noise was not properly accounted for in the calculation.
We thus use NNLS to get a non-negative solution, without the need to truncate negative values. 
In the (usual) case that inverse matrix multiplication yields only positive values, our NNLS will give the exact same optimal solution.

The correction matrices can be found (and changed) in the INI file (parameter @p correction_matrix of the corresponding labeling method).
However, these matrices for both 4-plex and 8-plex iTRAQ are now stable, and every kit delivered should have the same isotope correction values.
Thus, there should be no need to change them, but feel free to compare the values in the INI file with your kit's certificate.
For TMT (6-plex and 10-plex) the values have to be adapted for each kit: Modify the correction matrix according to the data in the product data sheet of your charge:
<pre>
Data sheet:
Mass Tag  Repoter Ion -2      -1      Monoisotopic    +1     +2
126       126.12776   0.0%    0.0%        100%        5.0%   0.0%
127N      127.124761  0.0%    0.2%        100%        4.6%   0.0%
...
</pre>
Corresponding correction matrix:
<pre>
[0.0/0.0/5.0/0.0,
0.0/0.2/4.6/0.0,
...
</pre>

After the quantitation, you may want to annotate the consensus features with corresponding peptide identifications,
obtained from an identification pipeline. Use @ref TOPP_IDMapper to perform the annotation, but make sure to set
suitably small RT and m/z tolerances for the mapping. Since the positions of the consensus features reported here 
are taken from the precursor of the MS2 (also if quant was done in MS3), it should be possible to achieve a 
perfect one-to-one matching of every identification (from MS2) to a single consensus feature.

Note that quantification will be solely on peptide level after this stage. In order to obtain protein quantities,
you can use @ref TOPP_TextExporter to obtain a simple text format which you can feed to other software tools (e.g., R),
or you can apply @ref TOPP_ProteinQuantifier.


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
  std::map<String, std::unique_ptr<IsobaricQuantitationMethod>> quant_methods_;
  std::map<String, String> quant_method_names_;

  void addMethod_(std::unique_ptr<IsobaricQuantitationMethod> ptr, std::string name)
  {
    std::string internal_name = ptr->getMethodName();
    quant_methods_[internal_name] = std::move(ptr);
    quant_method_names_[internal_name] = name;
  }

public:
  TOPPIsobaricAnalyzer() :
    TOPPBase("IsobaricAnalyzer", "Calculates isobaric quantitative values for peptides")
  {
    addMethod_(make_unique<ItraqFourPlexQuantitationMethod>(), "iTRAQ 4-plex");
    addMethod_(make_unique<ItraqEightPlexQuantitationMethod>(), "iTRAQ 8-plex");
    addMethod_(make_unique<TMTSixPlexQuantitationMethod>(), "TMT 6-plex");
    addMethod_(make_unique<TMTTenPlexQuantitationMethod>(), "TMT 10-plex");
    addMethod_(make_unique<TMTElevenPlexQuantitationMethod>(), "TMT 11-plex");
    addMethod_(make_unique<TMTSixteenPlexQuantitationMethod>(), "TMT 16-plex");
    addMethod_(make_unique<TMTEighteenPlexQuantitationMethod>(), "TMT 18-plex");
  }

protected:
  void registerOptionsAndFlags_() override
  {
    // initialize with the first available type
    registerStringOption_("type", "<mode>", quant_methods_.begin()->first, "Isobaric Quantitation method used in the experiment.", false);
    StringList valid_types;
    for (const auto& qm : quant_methods_)
    {
      valid_types.push_back(qm.first);
    }
    setValidStrings_("type", valid_types);

    registerInputFile_("in", "<file>", "", "input raw/picked data file ");
    setValidFormats_("in", {"mzML"});
    registerOutputFile_("out", "<file>", "", "output consensusXML file with quantitative information");
    setValidFormats_("out", {"consensusXML"});

    registerSubsection_("extraction", "Parameters for the channel extraction.");
    registerSubsection_("quantification", "Parameters for the peptide quantification.");
    for (const auto& qm : quant_methods_)
    {
      registerSubsection_(qm.second->getMethodName(), String("Algorithm parameters for ") + quant_method_names_[qm.second->getMethodName()]);
    }
  }

  Param getSubsectionDefaults_(const String& section) const override
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
      const auto it = quant_methods_.find(section);
      if (it == quant_methods_.end())
      { // should not happen
        throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Invalid subsection " + section);
      }
      return it->second->getParameters();
    }
  }

  ExitCodes main_(int, const char**) override
  {
    //-------------------------------------------------------------
    // parameter handling
    //-------------------------------------------------------------
    String in = getStringOption_("in");
    String out = getStringOption_("out");

    //-------------------------------------------------------------
    // loading input
    //-------------------------------------------------------------

    PeakMap exp;
    FileHandler().loadExperiment(in, exp, {FileTypes::MZML}, log_type_);

    //-------------------------------------------------------------
    // init quant method
    //-------------------------------------------------------------
    const auto& quant_method = quant_methods_[getStringOption_("type")];

    // set the parameters for this method
    quant_method->setParameters(getParam_().copy(quant_method->getMethodName() + ":", true));

    //-------------------------------------------------------------
    // calculations
    //-------------------------------------------------------------
    Param extract_param(getParam_().copy("extraction:", true));
    IsobaricChannelExtractor channel_extractor(quant_method.get());
    channel_extractor.setParameters(extract_param);

    ConsensusMap consensus_map_raw, consensus_map_quant;

    // extract channel information
    channel_extractor.extractChannels(exp, consensus_map_raw);

    IsobaricQuantifier quantifier(quant_method.get());
    Param quant_param(getParam_().copy("quantification:", true));
    quantifier.setParameters(quant_param);

    quantifier.quantify(consensus_map_raw, consensus_map_quant);

    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------

    //annotate output with data processing info
    addDataProcessing_(consensus_map_quant, getProcessingInfo_(DataProcessing::QUANTITATION));

    // add filename references
    for (auto& column : consensus_map_quant.getColumnHeaders())
    {
      column.second.filename = in;
    }

    const auto empty_feat = [](const ConsensusFeature& c){return c.getPeptideIdentifications().empty() && c.metaValueExists("all_empty") && c.getMetaValue("all_empty") == "true";};
    consensus_map_quant.erase(remove_if(consensus_map_quant.begin(), consensus_map_quant.end(), empty_feat), consensus_map_quant.end());
    consensus_map_quant.ensureUniqueId();
    FileHandler().storeConsensusFeatures(out, consensus_map_quant);

    return EXECUTION_OK;
  }

};

int main(int argc, const char** argv)
{
  TOPPIsobaricAnalyzer tool;
  return tool.main(argc, argv);
}

/// @endcond
