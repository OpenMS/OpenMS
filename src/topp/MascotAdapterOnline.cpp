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
// $Authors: Andreas Bertsch, Daniel Jameson, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/MascotXMLFile.h>
#include <OpenMS/FORMAT/MascotRemoteQuery.h>
#include <OpenMS/FORMAT/MascotGenericFile.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/SYSTEM/File.h>

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <sstream>

#include <QtCore/QDir>
#include <QtCore/QFile>
#include <QtCore/QCoreApplication>
#include <QtCore/QTimer>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**

    @page TOPP_MascotAdapterOnline MascotAdapterOnline

    @brief Identifies peptides in MS/MS spectra via Mascot.

<CENTER>
    <table>
        <tr>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
            <td VALIGN="middle" ROWSPAN=2> \f$ \longrightarrow \f$ MascotAdapterOnline \f$ \longrightarrow \f$</td>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> any signal-/preprocessing tool @n (that writes mzML format)</td>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_IDFilter or @n any protein/peptide processing tool</td>
        </tr>
    </table>
</CENTER>

    This wrapper application generates peptide identifications for MS/MS
    spectra using the search engine Mascot. It communicates with the Mascot
    server over the network (i.e. it does not have to run on the server
    itself).

    The adapter supports Mascot security features as well as proxy connections.
    Mascot versions 2.2.x up to 2.4.1 are supported and have been successfully
    tested (to varying degrees).

    @bug Running the adapter on Mascot 2.4 (possibly also other versions) produces the following error messages, which should be ignored:\n
    MascotRemoteQuery: An error occurred (requestId=11): Request aborted (QT Error Code: 7)\n
    MascotRemoteQuery: An error occurred (requestId=12): Request aborted (QT Error Code: 7)

    @note Some Mascot server instances seem to fail without reporting back an
    error message. In such cases, try to run the search on another Mascot
    server or change/validate the search parameters (e.g. using modifications
    that are known to Mascot and can thus be set in the INI file, but which are
    unknown to Mascot, might pose a problem).

    @note Mascot returns incomplete/incorrect protein assignments for most
    identified peptides (due to protein-level grouping/filtering). By default
    the protein associations are therefore not included in the output of this
    adapter, only the peptide sequences. @ref TOPP_PeptideIndexer should be run
    after this tool to get correct assignments. The flag @p keep_protein_links
    can be used to override this behavior.

    @note Currently mzIdentML (mzid) is not directly supported as an
    input/output format of this tool. Convert mzid files to/from idXML using
    @ref TOPP_IDFileConverter if necessary.

    <B>The command line parameters of this tool are:</B>
    @verbinclude TOPP_MascotAdapterOnline.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude TOPP_MascotAdapterOnline.html

    For the parameters of the algorithm section see the algorithms documentation: @n
    @ref OpenMS::MascotRemoteQuery "Mascot_server" @n
    @ref OpenMS::MascotGenericFile "Mascot_parameters" @n

*/

// We do not want this class to show up in the doc:
/// @cond TOPPCLASSES


class TOPPMascotAdapterOnline :
  public TOPPBase
{
public:
  TOPPMascotAdapterOnline() :
    TOPPBase("MascotAdapterOnline", "Annotates MS/MS spectra using Mascot.")
  {
  }

protected:

  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "input file in mzML format.\n");
    setValidFormats_("in", ListUtils::create<String>("mzML"));
    registerOutputFile_("out", "<file>", "", "output file in idXML format.\n");
    setValidFormats_("out", ListUtils::create<String>("idXML"));

    registerSubsection_("Mascot_server", "Mascot server details");
    registerSubsection_("Mascot_parameters", "Mascot parameters used for searching");
    registerFlag_("keep_protein_links", "The Mascot response file usually returns incomplete/wrong protein hits, so re-indexing the peptide hits is required. To avoid confusion why there are so few protein hits and force re-indexing, no proteins should be reported. To see the original (wrong) list, enable this flag.", true);
  }

  Param getSubsectionDefaults_(const String& section) const override
  {
    if (section == "Mascot_server")
    {
      MascotRemoteQuery mascot_query;
      return mascot_query.getParameters();
    }

    if (section == "Mascot_parameters")
    {
      MascotGenericFile mgf_file;
      Param p = mgf_file.getParameters();
      p.remove("internal:");
      return p;
    }

    return Param();
  }

  ExitCodes main_(int argc, const char** argv) override
  {
    //-------------------------------------------------------------
    // parameter handling
    //-------------------------------------------------------------

    // input/output files
    String in(getStringOption_("in")), out(getStringOption_("out"));
    FileHandler fh;
    FileTypes::Type in_type = fh.getType(in);

    //-------------------------------------------------------------
    // loading input
    //-------------------------------------------------------------

    PeakMap exp;
    // keep only MS2 spectra
    fh.getOptions().addMSLevel(2);
    fh.loadExperiment(in, exp, in_type, log_type_, false, false);
    writeLog_("Number of spectra loaded: " + String(exp.size()));

    if (exp.getSpectra().empty())
    {
      throw OpenMS::Exception::FileEmpty(__FILE__, __LINE__, __FUNCTION__, "Error: No MS2 spectra in input file.");
    }

    // determine type of spectral data (profile or centroided)
    SpectrumSettings::SpectrumType spectrum_type = exp[0].getType();

    if (spectrum_type == SpectrumSettings::PROFILE)
    {
      if (!getFlag_("force"))
      {
        throw OpenMS::Exception::IllegalArgument(__FILE__, __LINE__, __FUNCTION__, "Error: Profile data provided but centroided MS2 spectra expected. To enforce processing of the data set the -force flag.");
      }
    }

    //-------------------------------------------------------------
    // calculations
    //-------------------------------------------------------------

    Param mascot_param = getParam_().copy("Mascot_parameters:", true);

    // overwrite default search title with filename
    if (mascot_param.getValue("search_title") == "OpenMS_search")
    {
      mascot_param.setValue("search_title", File::removeExtension(File::basename(in)));
    }

    mascot_param.setValue("internal:HTTP_format", "true");
    MascotGenericFile mgf_file;
    mgf_file.setParameters(mascot_param);

    // get the spectra into string stream
    writeDebug_("Writing MGF file to stream", 1);
    stringstream ss;
    mgf_file.store(ss, in, exp, true); // write in compact format

    // Usage of a QCoreApplication is overkill here (and ugly too), but we just use the
    // QEventLoop to process the signals and slots and grab the results afterwards from
    // the MascotRemotQuery instance
    char** argv2 = const_cast<char**>(argv);
    QCoreApplication event_loop(argc, argv2);
    MascotRemoteQuery* mascot_query = new MascotRemoteQuery(&event_loop);
    Param mascot_query_param = getParam_().copy("Mascot_server:", true);
    writeDebug_("Setting parameters for Mascot query", 1);
    mascot_query->setParameters(mascot_query_param);
    writeDebug_("Setting spectra for Mascot query", 1);
    mascot_query->setQuerySpectra(ss.str());

    // remove unnecessary spectra
    ss.clear();

    QObject::connect(mascot_query, SIGNAL(done()), &event_loop, SLOT(quit()));
    QTimer::singleShot(1000, mascot_query, SLOT(run()));
    writeLog_("Submitting Mascot query (now: " + DateTime::now().get() + ")...");
    event_loop.exec();
    writeLog_("Mascot query finished");

    if (mascot_query->hasError())
    {
      writeLog_("An error occurred during the query: " + mascot_query->getErrorMessage());
      delete mascot_query;
      return EXTERNAL_PROGRAM_ERROR;
    }

    vector<PeptideIdentification> pep_ids;
    ProteinIdentification prot_id;

    if (!mascot_query_param.exists("skip_export") ||
        !mascot_query_param.getValue("skip_export").toBool())
    {
      // write Mascot response to file
      String mascot_tmp_file_name(File::getTempDirectory() + "/" + File::getUniqueName() + "_Mascot_response");
      QFile mascot_tmp_file(mascot_tmp_file_name.c_str());
      mascot_tmp_file.open(QIODevice::WriteOnly);
      mascot_tmp_file.write(mascot_query->getMascotXMLResponse());
      mascot_tmp_file.close();

      // set up helper object for looking up spectrum meta data:
      SpectrumMetaDataLookup lookup;
      MascotXMLFile::initializeLookup(lookup, exp);

      // read the response
      MascotXMLFile().load(mascot_tmp_file_name, prot_id, pep_ids, lookup);
      writeDebug_("Read " + String(pep_ids.size()) + " peptide ids and " + String(prot_id.getHits().size()) + " protein identifications from Mascot", 5);

      // for debugging errors relating to unexpected response files
      if (this->debug_level_ >= 100)
      {
        writeDebug_(String("\nMascot Server Response file saved to: '") + mascot_tmp_file_name + "'. If an error occurs, send this file to the OpenMS team.\n", 100);
      }
      else
      {
        mascot_tmp_file.remove(); // delete file
      }

      // keep or delete protein identifications?!
      if (!getFlag_("keep_protein_links"))
      {
        // remove protein links from peptides
        for (vector<PeptideIdentification>::iterator pep_it = pep_ids.begin();
             pep_it != pep_ids.end(); ++pep_it)
        {
          for (vector<PeptideHit>::iterator hit_it = pep_it->getHits().begin();
               hit_it != pep_it->getHits().end(); ++hit_it)
          {
            hit_it->setPeptideEvidences(vector<PeptideEvidence>());
          }
        }
        // remove proteins
        prot_id.getHits().clear();
      }
    }

    String search_number = mascot_query->getSearchIdentifier();
    if (search_number.empty())
    {
      writeLog_("Error: Failed to extract the Mascot search identifier (search number).");
      if (mascot_query_param.exists("skip_export") &&
          mascot_query_param.getValue("skip_export").toBool())
      {
        return PARSE_ERROR;
      }
    }
    else prot_id.setMetaValue("SearchNumber", search_number);

    // clean up
    delete mascot_query;

    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------

    vector<ProteinIdentification> prot_ids;
    prot_ids.push_back(prot_id);
    StringList ms_runs;
    exp.getPrimaryMSRunPath(ms_runs);
    prot_id.setPrimaryMSRunPath(ms_runs);
    IdXMLFile().store(out, prot_ids, pep_ids);
    
    return EXECUTION_OK;
  }

};


int main(int argc, const char** argv)
{
  TOPPMascotAdapterOnline tool;

  return tool.main(argc, argv);
}

/// @endcond
