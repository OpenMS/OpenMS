// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
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
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <iostream>

#include <vector>
#include <cmath>
using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_SpecLibCreator SpecLibCreator

@brief creates with given data a .MSP format spectral library.

Information file should have the following information: peptide, retention time, measured weight, charge state.
Extra information is allowed.

@experimental This Utility is not well tested and some features might not work as expected.

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_SpecLibCreator.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_SpecLibCreator.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPSpecLibCreator :
  public TOPPBase
{
public:
  TOPPSpecLibCreator() :
    TOPPBase("SpecLibCreator", "Creates an MSP formatted spectral library.")
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
      writeLogWarn_("Warning: Could not determine input file type!");
    }
    else if (in_type == FileTypes::MZDATA || in_type == FileTypes::MZXML)
    {
      FileHandler().loadExperiment(spec, msexperiment, {FileTypes::MZDATA, FileTypes::MZXML});
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
    FileHandler().storeExperiment(out, library, {FileTypes::MZDATA, FileTypes::MZXML, FileTypes::MSP});
    return EXECUTION_OK;
  }

};




int main(int argc, const char** argv)
{
  TOPPSpecLibCreator tool;
  return tool.main(argc, argv);
}

/// @endcond
