// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Fabian Aicheler $
// $Authors: Fabian Aicheler $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <OpenMS/SYSTEM/StopWatch.h>

#include <OpenMS/ANALYSIS/DECHARGING/MetaboliteFeatureDeconvolution.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/KERNEL/ConsensusMap.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_MetaboliteAdductDecharger MetaboliteAdductDecharger

@brief Decharges a feature map by clustering charge variants of metabolites to zero-charge entities.
<CENTER>
    <table>
        <tr>
            <th ALIGN = "center"> pot. predecessor tools </td>
            <td VALIGN="middle" ROWSPAN=2> &rarr; MetaboliteAdductDecharger &rarr;</td>
            <th ALIGN = "center"> pot. successor tools </td>
        </tr>
        <tr>
      <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_FeatureFinderMetabo </td>
        </tr>
    </table>
</CENTER>

The Decharger uses an ILP approach to group charge variants of the same metabolite, which
usually occur in ESI ionization mode. The resulting zero-charge metabolites, which are defined by RT and mass,
are written to consensusXML. Intensities of charge variants are summed up. The position of the zero charge
variant is the average of all clustered metabolites in each dimension (m and RT). For clustered metabolites,
the reported m/z is thus their neutral mass. For unclusted features with known charge, a default adduct
(protonation for positive mode, deprotonation for negative mode) is assumed to compute the neutral mass.
For unclustered features without known charge, m/z zero is reported.
It is also possible to include adducted species to the charge ladders (see 'potential_adducts' parameter).
Via this mechanism it is also possible to use this tool to find pairs/triples/quadruples/... in labeled data (by specifing the mass
tag weight as an adduct). If mass tags induce an RT shift (e.g. deuterium labeled data) you can also specify this also in the adduct list.
This will allow to tighten the RT search window, thus reducing false positive results.

This tool is derived from the method described in the following publication:

Bielow C, Ruzek S, Huber CG, Reinert K. Optimal decharging and clustering of charge ladders generated in ESI-MS. J Proteome Res 2010; 9: 2688.<br>
DOI: 10.1021/pr100177k

 <B>The command line parameters of this tool are:</B>
@verbinclude TOPP_MetaboliteAdductDecharger.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_MetaboliteAdductDecharger.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class UTILMetaboliteAdductDecharger :
  virtual public TOPPBase
{
public:
  UTILMetaboliteAdductDecharger() :
    TOPPBase("MetaboliteAdductDecharger", "Decharges and merges different feature charge variants of the same metabolite.")
  {
  }

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "input file ");
    setValidFormats_("in", ListUtils::create<String>("featureXML"));
    registerOutputFile_("out_cm", "<file>", "", "output consensus map", false);
    registerOutputFile_("out_fm", "<file>", "", "output feature map", false);
    registerOutputFile_("outpairs", "<file>", "", "output file", false);
    setValidFormats_("out_fm", ListUtils::create<String>("featureXML"));
    setValidFormats_("out_cm", ListUtils::create<String>("consensusXML"));
    setValidFormats_("outpairs", ListUtils::create<String>("consensusXML"));
    addEmptyLine_();
    registerSubsection_("algorithm", "Feature decharging algorithm section");
  }

  Param getSubsectionDefaults_(const String & /*section*/) const override
  {
    // there is only one subsection: 'algorithm' (s.a) .. and in it belongs the FeatureDecharger param
    MetaboliteFeatureDeconvolution fdc;
    Param tmp;
    tmp.insert("MetaboliteFeatureDeconvolution:", fdc.getParameters());
    return tmp;
  }

  ExitCodes main_(int, const char **) override
  {
    //-------------------------------------------------------------
    // parameter handling
    //-------------------------------------------------------------
    String infile = getStringOption_("in");
    String outfile_fm = getStringOption_("out_fm");
    String outfile_cm = getStringOption_("out_cm");
    String outfile_p = getStringOption_("outpairs");

    MetaboliteFeatureDeconvolution fdc;
    Param const & dc_param = getParam_().copy("algorithm:MetaboliteFeatureDeconvolution:", true);

    writeDebug_("Parameters passed to MetaboliteAdductDecharger", dc_param, 3);

    fdc.setParameters(dc_param);

    //-------------------------------------------------------------
    // loading input
    //-------------------------------------------------------------

    writeDebug_("Loading input file", 1);

    typedef FeatureMap FeatureMapType;
    FeatureMapType map_in, map_out;
    FileHandler().loadFeatures(infile, map_in, {FileTypes::FEATUREXML});

    //-------------------------------------------------------------
    // calculations
    //-------------------------------------------------------------
    ConsensusMap cm, cm2;
    StopWatch a;
    a.start();
    fdc.compute(map_in, map_out, cm, cm2);
    a.stop();
    //std::cerr << "took: " << a.getClockTime() << " seconds\n\n\n";

    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------

    writeDebug_("Saving output files", 1);

    cm.getColumnHeaders()[0].filename = infile;
    cm2.getColumnHeaders()[0].filename = infile;

    //annotate output with data processing info
    addDataProcessing_(map_out, getProcessingInfo_(DataProcessing::CHARGE_DECONVOLUTION));
    addDataProcessing_(cm, getProcessingInfo_(DataProcessing::CHARGE_DECONVOLUTION));
    addDataProcessing_(cm2, getProcessingInfo_(DataProcessing::CHARGE_DECONVOLUTION));


    FileHandler f;
    if (!outfile_cm.empty()) f.storeConsensusFeatures(outfile_cm, cm, {FileTypes::CONSENSUSXML});

    if (!outfile_p.empty()) f.storeConsensusFeatures(outfile_p, cm2, {FileTypes::CONSENSUSXML});
    if (!outfile_fm.empty()) FileHandler().storeFeatures(outfile_fm, map_out, {FileTypes::FEATUREXML});

    return EXECUTION_OK;
  }

};


int main(int argc, const char ** argv)
{
  UTILMetaboliteAdductDecharger tool;
  return tool.main(argc, argv);
}

/// @endcond
