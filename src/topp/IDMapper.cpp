// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Authors: Marc Sturm, Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/config.h>

#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/MzIdentMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/MzQuantMLFile.h>
#include <OpenMS/METADATA/MSQuantifications.h>

#include <OpenMS/ANALYSIS/ID/IDMapper.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/LogStream.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
    @page TOPP_IDMapper IDMapper

    @brief Assigns protein/peptide identifications to features or consensus features.

    <CENTER>
    <table>
        <tr>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> potential predecessor tools </td>
            <td VALIGN="middle" ROWSPAN=3> \f$ \longrightarrow \f$ IDMapper \f$ \longrightarrow \f$</td>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> potential successor tools </td>
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

    The mapping is based on retention times and mass-to-charge values. Roughly, a peptide identification is assigned to a (consensus) feature if its position lies within the boundaries of the feature or close enough to the feature centroid.
    Peptide identifications that don't match anywhere are still recorded in the resulting map, as "unassigned peptides". Protein identifications are annotated to the whole map, i.e. not to any particular (consensus) feature.

    In all cases, tolerance in RT and m/z dimension is applied according to the parameters @p rt_tolerance and @p mz_tolerance. Tolerance is understood as "plus or minus x", so the matching range is actually increased by twice the tolerance value.

    If several features or consensus features overlap the position of a peptide identification (taking the allowed tolerances into account), the identification is annotated to all of them.

    <B>Annotation of feature maps (featureXML input):</B>\n
    If @em all features have at least one convex hull, peptide positions are matched against the bounding boxes of the convex hulls (of individual mass traces, if available) by default.
    If not, the positions of the feature centroids are used. The respective coordinates of the centroids are also used for matching (in place of the corresponding ranges from the bounding boxes)
    if @p feature:use_centroid_rt or @p feature:use_centroid_mz are true.

    <B>Annotation of consensus maps (consensusXML input):</B>\n
    Peptide positions are always matched against centroid positions. By default, the consensus centroids are used. However, if @p consensus:use_subelements is set, the centroids of sub-features are considered instead.
    In this case, a peptide identification is mapped to a consensus feature if any of its sub-features matches.

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
    been abandoned, you can now explicitly control it with the @p feature:use_centroid_mz parameter. @p feature:use_centroid_mz does not take into account m/z deviations in the monoisotopic mass trace, but this can be compensated by increasing @p mz_tolerance. The new implementation should work correctly even if the monoisotopic mass trace itself was not detected.

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPIDMapper :
  public TOPPBase
{

public:

  TOPPIDMapper() :
    TOPPBase("IDMapper", "Assigns protein/peptide identifications to features or consensus features.")
  {
  }

protected:

  void registerOptionsAndFlags_() override
  {
    registerInputFile_("id", "<file>", "", "Protein/peptide identifications file");
    setValidFormats_("id", ListUtils::create<String>("mzid,idXML"));
    registerInputFile_("in", "<file>", "", "Feature map/consensus map file");
    setValidFormats_("in", ListUtils::create<String>("featureXML,consensusXML,mzq"));
    registerOutputFile_("out", "<file>", "", "Output file (the format depends on the input file format).");
    setValidFormats_("out", ListUtils::create<String>("featureXML,consensusXML,mzq"));

    addEmptyLine_();
    IDMapper mapper;
    Param p = mapper.getParameters();
    registerDoubleOption_("rt_tolerance", "<value>", p.getValue("rt_tolerance"), "RT tolerance (in seconds) for the matching of peptide identifications and (consensus) features.\nTolerance is understood as 'plus or minus x', so the matching range increases by twice the given value.", false);
    setMinFloat_("rt_tolerance", 0.0);
    registerDoubleOption_("mz_tolerance", "<value>", p.getValue("mz_tolerance"), "m/z tolerance (in ppm or Da) for the matching of peptide identifications and (consensus) features.\nTolerance is understood as 'plus or minus x', so the matching range increases by twice the given value.", false);
    setMinFloat_("mz_tolerance", 0.0);
    registerStringOption_("mz_measure", "<choice>", p.getEntry("mz_measure").valid_strings[0], "Unit of 'mz_tolerance'.", false);
    setValidStrings_("mz_measure", p.getEntry("mz_measure").valid_strings);
    registerStringOption_("mz_reference", "<choice>", p.getEntry("mz_reference").valid_strings[0], "Source of m/z values for peptide identifications. If 'precursor', the precursor-m/z from the idXML is used. If 'peptide',\nmasses are computed from the sequences of peptide hits; in this case, an identification matches if any of its hits matches.\n('peptide' should be used together with 'feature:use_centroid_mz' to avoid false-positive matches.)", false);
    setValidStrings_("mz_reference", p.getEntry("mz_reference").valid_strings);
    registerFlag_("ignore_charge", "For feature/consensus maps: Assign an ID independently of whether its charge state matches that of the (consensus) feature.");

    addEmptyLine_();
    registerTOPPSubsection_("feature", "Additional options for featureXML input");
    registerFlag_("feature:use_centroid_rt", "Use the RT coordinates of the feature centroids for matching, instead of the RT ranges of the features/mass traces.");
    registerFlag_("feature:use_centroid_mz", "Use the m/z coordinates of the feature centroids for matching, instead of the m/z ranges of the features/mass traces.\n(If you choose 'peptide' as 'mz_reference', you should usually set this flag to avoid false-positive matches.)");

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
    // LOG_DEBUG << "Starting..." << endl;

    //----------------------------------------------------------------
    // load ids
    //----------------------------------------------------------------
    // LOG_DEBUG << "Loading idXML..." << endl;
    String id = getStringOption_("id");
    vector<ProteinIdentification> protein_ids;
    vector<PeptideIdentification> peptide_ids;
    FileTypes::Type in_type = FileHandler::getType(id);
    if (in_type == FileTypes::IDXML)
    {
      IdXMLFile().load(id, protein_ids, peptide_ids);
    }
    else if (in_type == FileTypes::MZIDENTML)
    {
      MzIdentMLFile().load(id, protein_ids, peptide_ids);
    }
    else
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__,
                                       OPENMS_PRETTY_FUNCTION,
                                       "wrong id fileformat");
    }

    String in = getStringOption_("in");
    String spectra = getStringOption_("spectra:in");
    String out = getStringOption_("out");
    in_type = FileHandler::getType(in);
    //----------------------------------------------------------------
    //create mapper
    //----------------------------------------------------------------
    // LOG_DEBUG << "Creating mapper..." << endl;
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
      // LOG_DEBUG << "Processing consensus map..." << endl;
      ConsensusXMLFile file;
      ConsensusMap map;
      file.load(in, map);

      PeakMap exp;
      if (!spectra.empty())
      {
        MzMLFile().load(spectra, exp);
      }

      bool measure_from_subelements = getFlag_("consensus:use_subelements");
      bool annotate_ids_with_subelements = getFlag_("consensus:annotate_ids_with_subelements");

      mapper.annotate(map, peptide_ids, protein_ids, 
                      measure_from_subelements, annotate_ids_with_subelements,
                      exp);

      // annotate output with data processing info
      addDataProcessing_(map, getProcessingInfo_(DataProcessing::IDENTIFICATION_MAPPING));

      // sort list of peptide identifications in each consensus feature by map index
      map.sortPeptideIdentificationsByMapIndex();

      file.store(out, map);
    }

    //----------------------------------------------------------------
    // featureXML
    //----------------------------------------------------------------
    if (in_type == FileTypes::FEATUREXML)
    {
      // LOG_DEBUG << "Processing feature map..." << endl;
      FeatureMap map;
      FeatureXMLFile file;
      file.load(in, map);

      PeakMap exp;

      if (!spectra.empty())
      {
        MzMLFile().load(spectra, exp);
      }

      mapper.annotate(map, peptide_ids, protein_ids,
                      getFlag_("feature:use_centroid_rt"),
                      getFlag_("feature:use_centroid_mz"),
                      exp);

      //annotate output with data processing info
      addDataProcessing_(map, getProcessingInfo_(DataProcessing::IDENTIFICATION_MAPPING));

      file.store(out, map);
    }

    //----------------------------------------------------------------
    // MzQuantML
    //----------------------------------------------------------------
    if (in_type == FileTypes::MZQUANTML)
    {
      // LOG_DEBUG << "Processing mzq ..." << endl;
      MSQuantifications msq;
      MzQuantMLFile file;
      file.load(in, msq);

      bool measure_from_subelements = getFlag_("consensus:use_subelements");
      for (std::vector<ConsensusMap>::iterator it = msq.getConsensusMaps().begin(); it != msq.getConsensusMaps().end(); ++it)
      {
        mapper.annotate(*it, peptide_ids, protein_ids, measure_from_subelements);
        //annotate output with data processing info
        addDataProcessing_(*it, getProcessingInfo_(DataProcessing::IDENTIFICATION_MAPPING));
      }

      //~ writeDebug_(msq.getConsensusMaps().size(),3);
      //~ writeDebug_(msq.getConsensusMaps().back().size(),3);
      //~ writeDebug_(msq.getAnalysisSummary().quant_type_,3);
      file.store(out, msq);
    }

    // LOG_DEBUG << "Done." << endl;
    return EXECUTION_OK;
  }

};


int main(int argc, const char** argv)
{
  TOPPIDMapper tool;
  return tool.main(argc, argv);
}

/// @endcond
