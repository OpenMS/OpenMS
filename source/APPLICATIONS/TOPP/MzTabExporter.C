// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <OpenMS/DATASTRUCTURES/StringList.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/SVOutStream.h>

#include <boost/math/special_functions/fpclassify.hpp>

#include <vector>
#include <algorithm>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
   @page TOPP_MzTabExporter MzTabExporter

   @brief This application converts several %OpenMS XML formats (featureXML, consensusXML, and idXML) to mzTab.

  <CENTER>
    <table>
     <tr>
      <td ALIGN = "center" BGCOLOR="#EBEBEB"> potential predecessor tools </td>
         <td VALIGN="middle" ROWSPAN=2> \f$ \longrightarrow \f$ MzTabExporter \f$ \longrightarrow \f$</td>
     <td ALIGN = "center" BGCOLOR="#EBEBEB"> potential successor tools </td>
    </tr>
    <tr>
      <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> almost any TOPP tool </td>
      <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> external tools (MS Excel, OpenOffice, Notepad)</td>
    </tr>
   </table>
  </CENTER>

  See the mzTab specification for details on the format.

 <B>The command line parameters of this tool are:</B>
  @verbinclude TOPP_MzTabExporter.cli
 */

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

namespace OpenMS
{
class TOPPMzTabExporter : public TOPPBase
{
public:
  TOPPMzTabExporter() :
    TOPPBase("MzTabExporter", "Exports various XML formats to a mzTab file.")
  {
  }

protected:

  void registerOptionsAndFlags_()
  {
    registerInputFile_("in", "<file>", "", "Input file ");
    //setValidFormats_("in", StringList::create("featureXML,consensusXML,idXML,mzML"));
    setValidFormats_("in", StringList::create("idXML"));
    registerOutputFile_("out", "<file>", "", "Output file (mzTab)", true);
    registerFlag_("report_all_psm", "Report all peptide spectrum matches of a single spectrum", false);
  }

  /// Extracts, if possible a unique protein accession for a peptide hit in mzTab format. Otherwise NA is returned
  static String extractProteinAccession_(const PeptideHit& peptide_hit)
  {
    // if unique protein is present peptide can be assigned
    String accession;
    if (peptide_hit.getProteinAccessions().size() == 1)
    {
      accession = peptide_hit.getProteinAccessions()[0];
    } else
    {
      accession = "NA"; // no unique accession
    }
    return accession;
  }

  /// Extracts, modifications and positions of a peptide hit in mzTab format
  static String extractPeptideModifications_(const PeptideHit& peptide_hit)
  {
    String mods;

    const AASequence& aa_seq = peptide_hit.getSequence();
    bool first = true;
    for (Size i = 0; i != aa_seq.size(); ++i)
    {
      if ( aa_seq[i].isModified() )
      {
        String position = String( i + 1 );
        String reliability = "[1.0]";
        String unimod_name = aa_seq[i].getModification();
        if ( !first )
        {
          mods +=  ", ";
        } else
        {
          first = false;
        }
        mods += position + reliability + "UNIMOD:" + unimod_name;
      }
    }

    if (mods.length() == 0)
    {
      return "--";
    }

    return mods;
  }


  static String mapSearchEngineToCvParam_(const String& openms_search_engine_name)
  {
    if (openms_search_engine_name == "OMSSA")
    {
      return "[MS,MS:1001475,OMSSA,]";
    } else if (openms_search_engine_name == "Mascot")
    {
      return "[MS,MS:1001207,MASCOT,]";
    } else if (openms_search_engine_name == "XTandem")
    {
      return "[MS,MS:1001476,xtandem,]";
    } else if (openms_search_engine_name == "SEQUEST")
    {
      return "[MS,MS:1001208,Sequest,]";
    } else if (openms_search_engine_name == "CompNovo")
    {
      return "[MS,MS:UNKOWN,CompNovo,]";
    } else
    {
      return "--";
    }
    /*
    TODO:
    additional search engine strings in OpenMS:
    OpenMS/ConsensusID
        InsPecT
        PILIS
        In-silico digestion
        PepNovo
        TurboSEQUEST SEQUEST ???
     */
  }

  static String mapSearchEngineScoreToCvParam_(const String& openms_search_engine_name, DoubleReal score, String score_type)
  {
    String s;

    if (score_type.hasSubstring("Consensus"))
    {
      s = "[MS,MS:UNKNOWN,Consensus:score,";
    } else if (score_type == "q-value")
    {
      s = "[MS,MS:UNKOWN,qvalue,";
    } else if (score_type == "FDR")
    {
      s = "[MS,MS:UNKOWN,FDR,";
    } else if (score_type == "Posterior Error Probability")
    {
      s = "[MS,MS:UNKOWN,PEP,";
    } else if (score_type == "PhosphoScore")
    {
      s = "[MS,MS:UNKOWN,PhosphoScore,";
    } else if (openms_search_engine_name == "OMSSA")
    {
      s = "[MS,MS:1001328,OMSSA:evalue,";
    } else if (openms_search_engine_name == "Mascot")
    {
      s = "[MS,MS:1001171,MASCOT:score,";
    } else if (openms_search_engine_name == "XTandem")
    {
      s = "[MS,MS:1001330,X!Tandem:expect,";
    } else if (openms_search_engine_name == "SEQUEST")
    {
      s = "[MS,MS:1001155,Sequest:xcorr,";
    } else if (openms_search_engine_name == "CompNovo")
    {
      s = "[MS,MS:UNKOWN,CompNovo,";
    } else
    {
      return "--";
    }

    s += String::number(score, 8) + "]";
    return s;
    /*
    TODO:
    additional search engine strings in OpenMS:
    OpenMS/ConsensusID
        InsPecT
        PILIS und PILIS-E-value
        In-silico digestion
        PepNovo
        TurboSEQUEST SEQUEST ???
        InterProphet probability
        ProteinProphet probability
     */
  }


  ExitCodes main_( int, const char** )
  {
    // parameter handling
    String in = getStringOption_("in");
    String out = getStringOption_("out");

    // input file type
    FileTypes::Type in_type = FileHandler::getType( in );
    writeDebug_(String("Input file type: ") + FileHandler::typeToName(in_type), 2);

    if ( in_type == FileTypes::UNKNOWN )
    {
      writeLog_("Error: Could not determine input file type!");
      return PARSE_ERROR;
    }

    if (in_type == FileTypes::IDXML)
    {
      vector<ProteinIdentification> prot_ids;
      vector<PeptideIdentification> pep_ids;
      String document_id;
      IdXMLFile().load( in, prot_ids, pep_ids, document_id );
      bool has_coverage = true;
      try
      { // might throw Exception::MissingInformation() if no protein sequence information is added
        for (Size i = 0; i < prot_ids.size(); ++i)
        {
          prot_ids[i].computeCoverage(pep_ids);
        }
      }
      catch (Exception::MissingInformation& e)
      {
        LOG_WARN << e.what() << "\n";
        has_coverage = false;
      }

      ofstream txt_out( out.c_str() );
      SVOutStream output( txt_out, "\t", "_", String::NONE ); // according to MzTab specification

      // every ProteinIdentification corresponds to a search engine run and contains protein hits

      // write meta data of all runs
      Size run_count = 0;
      bool meta_info_printed = false;
      for (vector<ProteinIdentification>::const_iterator it = prot_ids.begin(); it != prot_ids.end(); ++it)
      {
        String UNIT_ID = File::basename(in) + "-" + String(run_count);
        String title = document_id;
        if ( title != "" )
        {
          output << "MOD" << UNIT_ID + "-title" << title;
          meta_info_printed = true;
        }
        run_count++;
      }

      if (meta_info_printed)
      {
        output << endl;
      }

      // write protein table header
      if (meta_info_printed)
      {
        output << endl;
      }
      output << "PRH" << "accession" << "unit_id" << "description" << "taxid"
             << "species" << "database" << "database_version" << "search_engine"
             << "search_engine_score" << "reliability" << "num_peptides" << "num_peptides_distinct"
             << "num_peptides_unambiguous" << "ambiguity_members" << "modifications" << "uri"
             << "go_terms" << "protein_coverage" << endl;

      // write protein table data
      run_count = 0;
      for (vector<ProteinIdentification>::const_iterator prot_id_it = prot_ids.begin();
           prot_id_it != prot_ids.end(); ++prot_id_it, ++run_count)
      {
        // TODO: maybe save these ProteinIdentification run properties in meta data
        // it->getScoreType()
        // it->isHigherScoreBetter())
        // it->getDateTime().toString(Qt::ISODate).toStdString()
        // it->getSearchEngineVersion();

        // search parameters
        const ProteinIdentification::SearchParameters& sp = prot_id_it->getSearchParameters();
        // TODO: maybe save these SearchParameters properties in a user param
        // String charges; ///< The allowed charges for the search
        // PeakMassType mass_type; ///< Mass type of the peaks
        // std::vector<String> fixed_modifications; ///< Used fixed modifications
        // std::vector<String> variable_modifications; ///< Allowed variable modifications
        // ProteinIdentification::NamesOfDigestionEnzyme[sp.enzyme]
        // UInt missed_cleavages; ///< The number of allowed missed cleavages
        // DoubleReal peak_mass_tolerance; ///< Mass tolerance of fragment ions (Dalton)
        // DoubleReal precursor_tolerance; ///< Mass tolerance of precursor ions (Dalton)

        // in OpenMS global to a ProteinIdentification
        String UNIT_ID_String = File::basename(in) + "-" + String(run_count);
        String database_String = (sp.db != "" ? sp.db : "--");
        String database_version_String = (sp.db_version != "" ? sp.db_version : "--");
        String species_String =  (sp.taxonomy != "0" ? sp.taxonomy : "--");
        String search_engine_cvParams = mapSearchEngineToCvParam_(prot_id_it->getSearchEngine());        
        String openms_search_engine_name = prot_id_it->getSearchEngine();
        //
        for (vector<ProteinHit>::const_iterator protein_hit_it = prot_id_it->getHits().begin();
             protein_hit_it != prot_id_it->getHits().end(); ++protein_hit_it)
        {
          String accession = protein_hit_it->getAccession();
          String unit_id = UNIT_ID_String; // run specific in OpenMS
          String description = "--";  // TODO: support description in protein hit
          String taxid = "--"; // TODO: mapping to NCBI taxid needed
          String species = species_String; // run specific in OpenMS
          String database = database_String; // run specific in OpenMS
          String database_version = database_version_String; // run specific in OpenMS
          String search_engine = search_engine_cvParams;
          String search_engine_score = mapSearchEngineScoreToCvParam_(openms_search_engine_name,
                                                                      protein_hit_it->getScore(),
                                                                      prot_id_it->getScoreType());

          String reliability = "--";
          String num_peptides = "TODO";
          String num_peptides_distinct = "TODO";
          String num_peptides_unambiguous = "TODO";
          String ambiguity_members = "TODO";
          String modifications = "TODO";
          String uri = in;
          String go_terms = "--";
          String protein_coverage;

          if (has_coverage)
          {
            protein_coverage = String(protein_hit_it->getCoverage() / 100.0);
          } else
          {
            protein_coverage = "NA";
          }

          output << "PRT" << accession << unit_id << description << taxid
                 << species << database << database_version << search_engine
                 << search_engine_score << reliability << num_peptides << num_peptides_distinct
                 << num_peptides_unambiguous << ambiguity_members << modifications << uri
                 << go_terms << protein_coverage << endl;
        }
      }

      // write peptide header
      output << endl;
      output << "PEH" << "sequence" << "accession" << "unit_id" << "unique" << "database"
             << "database_version" << "search_engine" << "search_engine_score"
             << "modifications" << "retention_time" << "charge"
             << "mass_to_charge" << "uri" << endl;

      // iterate over runs of peptide identifications
      for (vector<PeptideIdentification>::iterator pep_id_it = pep_ids.begin();
           pep_id_it != pep_ids.end(); ++pep_id_it)
      {        
        // find the ProteinIdentification corresponding to the current PeptideIdentification
        String common_identifier = pep_id_it->getIdentifier();
        for (Size i = 0; i != prot_ids.size(); ++i)
        {
          if (prot_ids[i].getIdentifier() == common_identifier)
          {
            run_count = i;
          }
        }

        // TODO: bad design of Protein/PeptideIdentification as search engine parameters are stored in prot.
        String openms_search_engine_name = prot_ids[run_count].getSearchEngine();
        String search_engine_cvParams = mapSearchEngineToCvParam_(openms_search_engine_name);

        const ProteinIdentification::SearchParameters& sp = prot_ids[run_count].getSearchParameters();
        String UNIT_ID_String = File::basename(in) + "-" + String(run_count);
        String database_String = (sp.db != "" ? sp.db : "--");
        String database_version_String = (sp.db_version != "" ? sp.db_version : "--");

        // sort according to score
        pep_id_it->assignRanks();

        bool report_all_psm = getFlag_("report_all_psm");
        for (vector<PeptideHit>::const_iterator peptide_hit_it = pep_id_it->getHits().begin();
             peptide_hit_it != pep_id_it->getHits().end(); ++peptide_hit_it)
        {
          String sequence = peptide_hit_it->getSequence().toString();
          String accession = extractProteinAccession_(*peptide_hit_it);
          String unit_id = UNIT_ID_String;
          String unique;
          String database = database_String;
          String database_version = database_version_String;
          String search_engine = search_engine_cvParams;
          String search_engine_score = mapSearchEngineScoreToCvParam_(openms_search_engine_name,
                                                                      peptide_hit_it->getScore(),
                                                                      pep_id_it->getScoreType());
          String modifications = extractPeptideModifications_(*peptide_hit_it); //TODO: check if terminal mods work

          // if unique protein is present peptide can be assigned
          if (peptide_hit_it->getProteinAccessions().size() == 1)
          {
            unique = "1";
          } else
          {
            unique = "0";
          }

          String retention_time;
          if (pep_id_it->metaValueExists("RT")) // Note: RT stored on pep_id_it not on hit
          {
            retention_time = String::number(String(pep_id_it->getMetaValue("RT")).toDouble(), 2);
          } else
          {
            retention_time = "--";
          }

          String mass_to_charge;
          if (pep_id_it->metaValueExists("MZ")) // Note: MZ stored on pep_id_it not on hit
          {
            mass_to_charge = String::number(String(pep_id_it->getMetaValue("MZ")).toDouble(), 10);
          } else
          {
            mass_to_charge = "--";
          }

          String charge = peptide_hit_it->getCharge();
          String uri = in;

          output << "PEP" << sequence << accession << unit_id << unique << database
                 << database_version << search_engine << search_engine_score
                 << modifications << retention_time << charge
                 << mass_to_charge << uri << endl;

          // usually only the best peptide spectrum match is reported
          if (!report_all_psm)
          {
            break;
          }

        }
      }
      txt_out.close();
    }

    return EXECUTION_OK;
  }
};
}

int main( int argc, const char** argv )
{
  TOPPMzTabExporter t;
  return t.main(argc, argv);
}

/// @endcond
