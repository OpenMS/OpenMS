// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm, Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/IDMapper.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/MzIdentMLFile.h>
#include <OpenMS/config.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_IDMapper IDMapper

@brief Assigns protein/peptide identifications to features or consensus features.

<CENTER>
<table>
    <tr>
        <th ALIGN = "center"> potential predecessor tools </td>
        <td VALIGN="middle" ROWSPAN=3> &rarr; IDMapper &rarr;</td>
        <th ALIGN = "center"> potential successor tools </td>
    </tr>
    <tr>
        <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_MascotAdapter (or other ID engines) </td>
        <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_ConsensusID </td>
    </tr>
    <tr>
      <td VALIGN="middle" ALIGN="center" ROWSPAN=1> @ref TOPP_IDFilter </td>
      <td VALIGN="middle" ALIGN="center" ROWSPAN=1> @ref TOPP_MapAlignerIdentification </td>
    </tr>
</table>
</CENTER>

The mapping is based on retention times and mass-to-charge values. Roughly, a peptide identification is assigned to a (consensus) feature if its position lies within the boundaries of the feature
or close enough to the feature centroid. Peptide identifications that don't match anywhere are still recorded in the resulting map, as "unassigned peptides". Protein identifications are annotated
to the whole map, i.e. not to any particular (consensus) feature.

In all cases, tolerance in RT and m/z dimension is applied according to the parameters @p rt_tolerance and @p mz_tolerance. Tolerance is understood as "plus or minus x", so the matching range is
actually increased by twice the tolerance value.

If several features or consensus features overlap the position of a peptide identification (taking the allowed tolerances into account), the identification is annotated to all of them.

<B>Annotation of feature maps (featureXML input):</B>\n
If @em all features have at least one convex hull, peptide positions are matched against the bounding boxes of the convex hulls (of individual mass traces, if available) by default.
If not, the positions of the feature centroids are used. The respective coordinates of the centroids are also used for matching (in place of the corresponding ranges from the bounding boxes)
if @p feature:use_centroid_rt or @p feature:use_centroid_mz are true.

<B>Annotation of consensus maps (consensusXML input):</B>\n
Peptide positions are always matched against centroid positions. By default, the consensus centroids are used. However, if @p consensus:use_subelements is set, the centroids of sub-features are
considered instead. In this case, a peptide identification is mapped to a consensus feature if any of its sub-features matches.

@note Currently mzIdentML (mzid) is not directly supported as an input/output format of this tool. Convert mzid files to/from idXML using @ref TOPP_IDFileConverter if necessary.

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_IDMapper.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_IDMapper.html

On the peptide side, two sources for m/z values are possible (see parameter @p mz_reference): 1. m/z of the precursor of the MS2 spectrum that gave rise to the peptide identification;
2. theoretical masses computed from the amino acid sequences of peptide hits.
(When using theoretical masses, make sure that peptide modifications were identified correctly. OpenMS currently "forgets" mass shifts that it can't assign to modifications - if that
happens, masses computed from peptide sequences will be off.)

@deprecated The parameter handling of this tool has been reworked. For greater consistency with other tools, the parameters @p rt_delta and @p mz_delta have been renamed to @p rt_tolerance
and @p mz_tolerance. The possible values of the @p mz_reference parameter have also been renamed. The default value of @p mz_tolerance has been increased from 1 ppm to a more realistic 20 ppm.\n
Most importantly, the @p use_centroids parameter from previous versions has been split into two parameters, @p feature:use_centroid_rt and @p feature:use_centroid_mz. In OpenMS 1.6, peptide
identifications would be matched only against monoisotopic mass traces of features if @p mz_reference was @p PeptideMass; otherwise, all mass traces would be used. This implicit behaviour has
been abandoned, you can now explicitly control it with the @p feature:use_centroid_mz parameter. @p feature:use_centroid_mz does not take into account m/z deviations in the monoisotopic mass
trace, but this can be compensated by increasing @p mz_tolerance. The new implementation should work correctly even if the monoisotopic mass trace itself was not detected.

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPIDMapper : public TOPPBase
{
public:
  TOPPIDMapper() : TOPPBase("IDMapper", "Assigns protein/peptide identifications to features or consensus features.")
  {
  }

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("id", "<file>", "", "Protein/peptide identifications file");
    setValidFormats_("id", ListUtils::create<String>("mzid,idXML"));
    registerInputFile_("in", "<file>", "", "Feature map/consensus map file");
    setValidFormats_("in", ListUtils::create<String>("featureXML,consensusXML"));
    registerOutputFile_("out", "<file>", "", "Output file (the format depends on the input file format).");
    setValidFormats_("out", ListUtils::create<String>("featureXML,consensusXML"));

    addEmptyLine_();
    IDMapper mapper;
    Param p = mapper.getParameters();
    registerDoubleOption_("rt_tolerance", "<value>", p.getValue("rt_tolerance"),
                          "RT tolerance (in seconds) for the matching of peptide identifications and (consensus) features.\nTolerance is understood as 'plus or minus x', so the matching range "
                          "increases by twice the given value.",
                          false);
    setMinFloat_("rt_tolerance", 0.0);
    registerDoubleOption_("mz_tolerance", "<value>", p.getValue("mz_tolerance"),
                          "m/z tolerance (in ppm or Da) for the matching of peptide identifications and (consensus) features.\nTolerance is understood as 'plus or minus x', so the matching range "
                          "increases by twice the given value.",
                          false);
    setMinFloat_("mz_tolerance", 0.0);
    registerStringOption_("mz_measure", "<choice>", p.getEntry("mz_measure").valid_strings[0], "Unit of 'mz_tolerance'.", false);
    setValidStrings_("mz_measure", ListUtils::toStringList<std::string>(p.getEntry("mz_measure").valid_strings));
    registerStringOption_(
      "mz_reference", "<choice>", p.getEntry("mz_reference").valid_strings[1],
      "Source of m/z values for peptide identifications. If 'precursor', the precursor-m/z from the idXML is used. If 'peptide',\nmasses are computed from the sequences of peptide hits; in this "
      "case, an identification matches if any of its hits matches.\n('peptide' should be used together with 'feature:use_centroid_mz' to avoid false-positive matches.)",
      false);
    setValidStrings_("mz_reference", ListUtils::toStringList<std::string>(p.getEntry("mz_reference").valid_strings));
    registerFlag_("ignore_charge", "For feature/consensus maps: Assign an ID independently of whether its charge state matches that of the (consensus) feature.", true);

    addEmptyLine_();
    registerTOPPSubsection_("feature", "Additional options for featureXML input");
    registerStringOption_("feature:use_centroid_rt", "<choice>", "false", "Use the RT coordinates of the feature centroids for matching, instead of the RT ranges of the features/mass traces.", false);
    setValidStrings_("feature:use_centroid_rt", ListUtils::create<String>("true,false"));
    registerStringOption_("feature:use_centroid_mz", "<choice>", "true",
                          "Use the m/z coordinates of the feature centroids for matching, instead of the m/z ranges of the features/mass traces.\n(If you choose 'peptide' as 'mz_reference', you "
                          "should usually set this flag to avoid false-positive matches.)",
                          false);
    setValidStrings_("feature:use_centroid_mz", ListUtils::create<String>("true,false"));

    addEmptyLine_();
    registerTOPPSubsection_("consensus", "Additional options for consensusXML input");
    registerFlag_("consensus:use_subelements", "Match using RT and m/z of sub-features instead of consensus RT and m/z. A consensus feature matches if any of its sub-features matches.");
    registerFlag_("consensus:annotate_ids_with_subelements", "Store the map index of the sub-feature in the peptide ID.", true);

    registerTOPPSubsection_("spectra", "Additional options for mzML input");
    registerInputFile_("spectra:in", "<file>", "", "MS run used to annotated unidentified spectra to features or consensus features.", false);
    setValidFormats_("spectra:in", ListUtils::create<String>("mzML"));
  }

  ExitCodes main_(int, const char**) override
  {
    // OPENMS_LOG_DEBUG << "Starting..." << endl;

    //----------------------------------------------------------------
    // load ids
    //----------------------------------------------------------------
    // OPENMS_LOG_DEBUG << "Loading idXML..." << endl;
    String id = getStringOption_("id");
    vector<ProteinIdentification> protein_ids;
    vector<PeptideIdentification> peptide_ids;
    FileTypes::Type in_type = FileHandler::getType(id);
    FileHandler().loadIdentifications(id, protein_ids, peptide_ids, {FileTypes::IDXML, FileTypes::MZIDENTML});

    String in = getStringOption_("in");
    String spectra = getStringOption_("spectra:in");
    String out = getStringOption_("out");
    in_type = FileHandler::getType(in);
    //----------------------------------------------------------------
    // create mapper
    //----------------------------------------------------------------
    // OPENMS_LOG_DEBUG << "Creating mapper..." << endl;
    IDMapper mapper;
    Param p = mapper.getParameters();
    p.setValue("rt_tolerance", getDoubleOption_("rt_tolerance"));
    p.setValue("mz_tolerance", getDoubleOption_("mz_tolerance"));
    p.setValue("mz_measure", getStringOption_("mz_measure"));
    p.setValue("mz_reference", getStringOption_("mz_reference"));
    p.setValue("ignore_charge", getFlag_("ignore_charge") ? "true" : "false");
    mapper.setParameters(p);

    //----------------------------------------------------------------
    // consensusXML
    //----------------------------------------------------------------
    if (in_type == FileTypes::CONSENSUSXML)
    {
      // OPENMS_LOG_DEBUG << "Processing consensus map..." << endl;
      FileHandler consensusFile;
      ConsensusMap map;
      consensusFile.loadConsensusFeatures(in, map, {FileTypes::CONSENSUSXML});

      PeakMap exp;
      if (!spectra.empty())
      {
        FileHandler().loadExperiment(spectra, exp, {FileTypes::MZML});
      }

      bool measure_from_subelements = getFlag_("consensus:use_subelements");
      bool annotate_ids_with_subelements = getFlag_("consensus:annotate_ids_with_subelements");

      mapper.annotate(map, peptide_ids, protein_ids, measure_from_subelements, annotate_ids_with_subelements, exp);

      // annotate output with data processing info
      addDataProcessing_(map, getProcessingInfo_(DataProcessing::IDENTIFICATION_MAPPING));

      // sort list of peptide identifications in each consensus feature by map index
      map.sortPeptideIdentificationsByMapIndex();

      consensusFile.storeConsensusFeatures(out, map, {FileTypes::CONSENSUSXML});
    }

    //----------------------------------------------------------------
    // featureXML
    //----------------------------------------------------------------
    if (in_type == FileTypes::FEATUREXML)
    {
      // OPENMS_LOG_DEBUG << "Processing feature map..." << endl;
      FeatureMap map;
      FileHandler featureFile;
      featureFile.loadFeatures(in, map, {FileTypes::FEATUREXML});

      PeakMap exp;

      if (!spectra.empty())
      {
        FileHandler().loadExperiment(spectra, exp, {FileTypes::MZML});
      }

      mapper.annotate(map, peptide_ids, protein_ids, (getStringOption_("feature:use_centroid_rt") == "true"), (getStringOption_("feature:use_centroid_mz") == "true"), exp);

      // annotate output with data processing info
      addDataProcessing_(map, getProcessingInfo_(DataProcessing::IDENTIFICATION_MAPPING));

      featureFile.storeFeatures(out, map, {FileTypes::FEATUREXML});
    }

    return EXECUTION_OK;
  }
};


int main(int argc, const char** argv)
{
  TOPPIDMapper tool;
  return tool.main(argc, argv);
}

/// @endcond
