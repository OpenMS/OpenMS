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
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/CsvFile.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/MzXMLFile.h>
#include <OpenMS/FORMAT/MzDataFile.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/FORMAT/MSPFile.h>
#include <iostream>

#include <vector>
#include <cmath>
using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
  @page UTILS_SpecLibCreator SpecLibCreator

    @brief creates with given data a .MSP format spectral library.

    Information file should have the following information: peptide, retention time, measured weight, charge state.
    Extra information is allowed.

    @experimental This Utility is not well tested and some features might not work as expected.

    <B>The command line parameters of this tool are:</B>
    @verbinclude UTILS_SpecLibCreator.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude UTILS_SpecLibCreator.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPSpecLibCreator :
  public TOPPBase
{
public:
  TOPPSpecLibCreator() :
    TOPPBase("SpecLibCreator", "Creates an MSP formatted spectral library.", false)
  {
  }

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("info", "<file>", "", "Holds id, peptide, retention time etc.");
    setValidFormats_("info", ListUtils::create<String>("csv"));

    registerStringOption_("itemseperator", "<char>", ",", " Separator between items. e.g. ,", false);
    registerStringOption_("itemenclosed", "<bool>", "false", "'true' or 'false' if true every item is enclosed e.g. '$peptide$,$run$...", false);
    setValidStrings_("itemenclosed", ListUtils::create<String>("true,false"));

    registerInputFile_("spec", "<file>", "", "spectra");
    setValidFormats_("spec", ListUtils::create<String>("mzData,mzXML"));

    registerOutputFile_("out", "<file>", "", "output MSP formatted spectra library");
    setValidFormats_("out", ListUtils::create<String>("msp"));
  }

  ExitCodes main_(int, const char**) override
  {
    //-------------------------------------------------------------
    // parameter handling
    //-------------------------------------------------------------

    String info = getStringOption_("info");
    String itemseperator = getStringOption_("itemseperator");
    String out = getStringOption_("out");
    bool itemenclosed;
    if (getStringOption_("itemenclosed") == "true")
    {
      itemenclosed  = true;
    }
    else
    {
      itemenclosed = false;
    }

    String spec = getStringOption_("spec");
    if (info == String::EMPTY)
    {
      throw Exception::RequiredParameterNotGiven(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "info");
    }
    if (spec == String::EMPTY)
    {
      throw Exception::RequiredParameterNotGiven(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "spec");
    }


    //-------------------------------------------------------------
    // loading input
    //-------------------------------------------------------------
    Int retention_time = -1;
    Int peptide = -1;
    Int measured_weight = -1;
    //UInt first_scan;
    UInt charge_state(0), Experimental_id(0); //,found_by, track, comment, vaccination_peptid,epitope, confident, hlaallele;
    const char* sepi = itemseperator.c_str();
    char sepo = *sepi;
    CsvFile csv_file(info, sepo, itemenclosed);
    vector<StringList>  list;

    list.resize(csv_file.rowCount());

    for (UInt i = 0; i < csv_file.rowCount(); ++i)
    {
      csv_file.getRow(i, list[i]);
    }
    for (UInt i = 0; i < list[0].size(); ++i)
    {

      if (list[0][i].toLower().removeWhitespaces().compare("retentiontime") == 0)
      {
        retention_time = i;
      }
      else if (list[0][i].toLower().hasSubstring("_id"))
      {
        Experimental_id = i;
      }
      else if (list[0][i].toLower() == "last scan")
      {
        // last_scan = i;
      }
      else if (list[0][i].toLower() == "modification")
      {
        // modification = i;
      }
      else if (list[0][i].toLower().removeWhitespaces().compare("chargestate") == 0 || list[0][i].toLower().removeWhitespaces().hasSubstring("charge"))
      {
        charge_state = i;
      }
      else if (list[0][i].toLower().trim().compare("peptide") == 0)
      {
        peptide = i;
      }
      else if (list[0][i].toLower().removeWhitespaces().hasSubstring("measuredweight")  || list[0][i].removeWhitespaces().compare("measuredweight[M+nH]n+") == 0)
      {
        measured_weight = i;
      }
    }
    if (retention_time  == -1)
    {
      throw Exception::RequiredParameterNotGiven(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "unclear which parameter is retention time");
    }
    if (peptide  == -1)
    {
      throw Exception::RequiredParameterNotGiven(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "unclear which parameter is peptide");
    }
    if (measured_weight  == -1)
    {
      throw Exception::RequiredParameterNotGiven(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "unclear which parameter is measured weight");
    }
    FileHandler fh;
    FileTypes::Type in_type = fh.getType(spec);
    PeakMap msexperiment;

    if (in_type == FileTypes::UNKNOWN)
    {
      writeLog_("Warning: Could not determine input file type!");
    }
    else if (in_type == FileTypes::MZDATA)
    {
      MzDataFile mzData;
      mzData.load(spec, msexperiment);
    }
    else if (in_type == FileTypes::MZXML)
    {
      MzXMLFile mzXML;
      mzXML.load(spec, msexperiment);
    }
    if (msexperiment.getMinRT() == 0)
    {
      throw Exception::RequiredParameterNotGiven(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "EMPTY??");
    }
    PeakMap library;

    //-------------------------------------------------------------
    // creating library
    //-------------------------------------------------------------
    UInt found_counter = 0;

    for (UInt i = 1; i < list.size(); ++i)
    {
      bool no_peptide = true;
      double rt =  (60 * (list[i][retention_time].toFloat())); // from minutes to seconds
      double mz = list[i][measured_weight].toFloat();
      for (PeakMap::Iterator it = msexperiment.begin(); it < msexperiment.end(); ++it)
      {
        if ((abs(rt - it->getRT()) < 5) && (abs(mz - it->getPrecursors()[0].getMZ()) < 0.1))
        {
          //if ( ceil(rt) == ceil(it->getRT()) || ceil(rt) == floor(it->getRT()) || floor(rt) == ceil(it->getRT()) || floor(rt) == floor(it->getRT()))
          ++found_counter;
          no_peptide = false;
          cout << "Found Peptide " << list[i][peptide] << " with id: " << list[i][Experimental_id] << "\n";
          cout << "rt: " << it->getRT() << " and mz: " << it->getPrecursors()[0].getMZ() << "\n";

          MSSpectrum speci;
          speci.setRT(it->getRT());
          speci.setMSLevel(2);
          speci.setPrecursors(it->getPrecursors());
          for (UInt j = 0; j < it->size(); ++j)
          {

            Peak1D richy;
            richy.setIntensity(it->operator[](j).getIntensity());
            richy.setPosition(it->operator[](j).getPosition());
            richy.setMZ(it->operator[](j).getMZ());
            richy.setPos(it->operator[](j).getPos()); //ALIAS for setMZ???

            speci.push_back(richy);
          }
          PeptideHit hit; // = *it->getPeptideIdentifications().begin()->getHits().begin();
          AASequence aa = AASequence::fromString(list[i][peptide]);
          hit.setSequence(aa);
          hit.setCharge(list[i][charge_state].toInt());
          vector<PeptideHit> hits;
          hits.push_back(hit);
          vector<PeptideIdentification> pepi;
          PeptideIdentification pep;
          pep.setHits(hits);
          pepi.push_back(pep);
          speci.setPeptideIdentifications(pepi);
          //it->getPeptideIdentifications().begin()->setHits(hits);
          library.addSpectrum(speci);
        }
      }
      if (no_peptide)
      {
        cout << "Peptide: " << list[i][peptide] << " not found\n";
      }
    }
    cout << "Found " << found_counter << " peptides\n";

    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------
    in_type = fh.getType(out);
    if (in_type == FileTypes::MZDATA)
    {
      MzDataFile f;
      f.store(out, library);
    }
    else if (in_type == FileTypes::MZXML)
    {
      MzXMLFile f;
      f.store(out, library);
    }
    else
    {
      MSPFile msp;
      msp.store(out, library);
    }

    return EXECUTION_OK;
  }

};




int main(int argc, const char** argv)
{
  TOPPSpecLibCreator tool;
  return tool.main(argc, argv);
}

/// @endcond
