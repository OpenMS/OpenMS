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
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------


#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/ANALYSIS/DENOVO/CompNovoIdentification.h>
#include <OpenMS/ANALYSIS/DENOVO/CompNovoIdentificationCID.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/ProteaseDB.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>

using namespace OpenMS;
using namespace std;

/**
 @page TOPP_CompNovo CompNovo

 @brief Performs a peptide/protein identification with the CompNovo engine.

<CENTER>
 <table>
  <tr>
   <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
   <td VALIGN="middle" ROWSPAN=2> \f$ \longrightarrow \f$ CompNovo \f$ \longrightarrow \f$</td>
   <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
  </tr>
  <tr>
   <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> any signal-/preprocessing tool @n (in mzML format)</td>
   <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_IDFilter or @n any protein/peptide processing tool</td>
  </tr>
 </table>
</CENTER>

 @experimental This TOPP-tool is not well tested and not all features might be properly implemented and tested!

 The CompNovo engine can operate in CID/ETD mode where the ETD spectra of the
  same precursor are used to support the de novo identification. In CID mode
  all spectra are assumed to be CID spectra only.

 The details are described in the publication:

 Andreas Bertsch, Andreas Leinenbach, Anton Pervukhin, Markus Lubeck,
 Ralf Hartmer,	Carsten Baessmann, Yasser A Elnakady, Rolf M&uuml;ller,
 Sebastian B&ouml;cker, Christian G Huber and Oliver Kohlbacher (2009)
 "De novo peptide sequencing by tandem MS using complementary CID and
 electron transfer dissociation"
 Electrophoresis, 30(21):3736-3747. (PubMed ID: 19862751)

 @experimental This implementation may contain bugs!

 @note Currently mzIdentML (mzid) is not directly supported as an input/output format of this tool. Convert mzid files to/from idXML using @ref TOPP_IDFileConverter if necessary.

 <B>The command line parameters of this tool are:</B>
 @verbinclude TOPP_CompNovo.cli
 <B>INI file documentation of this tool:</B>
 @htmlinclude TOPP_CompNovo.html
*/


// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPCompNovo :
  public TOPPBase
{
public:
  TOPPCompNovo() :
    TOPPBase("CompNovo", "Performs a de novo peptide identification using the CompNovo engine.")
  {
  }

protected:

  Param getSubsectionDefaults_(const String & /*section*/) const override
  {
    return CompNovoIdentification().getDefaults();
  }

  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "input file in mzML format", true);
    setValidFormats_("in", ListUtils::create<String>("mzML"));

    registerOutputFile_("out", "<file>", "", "output file in idXML format", true);
    setValidFormats_("out", ListUtils::create<String>("idXML"));

    registerSubsection_("algorithm", "Algorithm section");
    addEmptyLine_();
  }

  ExitCodes main_(int, const char **) override
  {
    //-------------------------------------------------------------
    // parameter handling
    //-------------------------------------------------------------

    //input/output files
    String in(getStringOption_("in"));
    String out(getStringOption_("out"));
    //-------------------------------------------------------------
    // loading input
    //-------------------------------------------------------------

    PeakMap exp;
    MzMLFile f;
    f.setLogType(log_type_);

    PeakFileOptions options;
    options.clearMSLevels();
    options.addMSLevel(2);
    f.getOptions() = options;
    f.load(in, exp);

    writeDebug_("Data set contains " + String(exp.size()) + " spectra", 1);

    //-------------------------------------------------------------
    // calculations
    //-------------------------------------------------------------

    vector<PeptideIdentification> pep_ids;
    CompNovoIdentification comp_novo_id;

    // set the options
    Param algorithm_param = getParam_().copy("algorithm:", true);
    comp_novo_id.setParameters(algorithm_param);
    comp_novo_id.getIdentifications(pep_ids, exp);
    algorithm_param = comp_novo_id.getParameters();

    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------

    DateTime now = DateTime::now();
    String date_string = now.get();
    String identifier("CompNovo_" + date_string);

    for (vector<PeptideIdentification>::iterator it = pep_ids.begin(); it != pep_ids.end(); ++it)
    {
      it->assignRanks();
      it->setIdentifier(identifier);
    }

    vector<ProteinIdentification> prot_ids;
    ProteinIdentification prot_id;
    prot_id.setIdentifier(identifier);
    prot_id.setDateTime(now);
    StringList ms_runs;
    exp.getPrimaryMSRunPath(ms_runs);
    prot_id.setPrimaryMSRunPath(ms_runs);

    ProteinIdentification::SearchParameters search_parameters;
    search_parameters.charges = "+2-+3";
    if (algorithm_param.getValue("tryptic_only").toBool())
    {
      search_parameters.digestion_enzyme = *(ProteaseDB::getInstance()->getEnzyme("Trypsin"));
    }
    else
    {
      search_parameters.digestion_enzyme = *(ProteaseDB::getInstance()->getEnzyme("no cleavage"));
    }
    search_parameters.mass_type = ProteinIdentification::MONOISOTOPIC;
    search_parameters.fixed_modifications = algorithm_param.getValue("fixed_modifications");
    search_parameters.variable_modifications = algorithm_param.getValue("variable_modifications");

    search_parameters.missed_cleavages = (UInt)algorithm_param.getValue("missed_cleavages");
    search_parameters.fragment_mass_tolerance = (double)algorithm_param.getValue("fragment_mass_tolerance");
    search_parameters.precursor_mass_tolerance = (double)algorithm_param.getValue("precursor_mass_tolerance");
    prot_id.setSearchParameters(search_parameters);
    prot_id.setSearchEngineVersion("0.9beta");
    prot_id.setSearchEngine("CompNovo");
    prot_ids.push_back(prot_id);

    IdXMLFile().store(out, prot_ids, pep_ids);

    return EXECUTION_OK;
  }

};

int main(int argc, const char ** argv)
{
  TOPPCompNovo tool;
  return tool.main(argc, argv);
}

/// @endcond
