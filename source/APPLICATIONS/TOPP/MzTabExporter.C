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
#include <OpenMS/FILTERING/ID/IDFilter.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>

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
      <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
         <td VALIGN="middle" ROWSPAN=2> \f$ \longrightarrow \f$ MzTabExporter \f$ \longrightarrow \f$</td>
     <td ALIGN = "center" BGCOLOR="#EBEBEB"> potential successor tools </td>
    </tr>
    <tr>
      <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> ProteinQuantifier </td>
      <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> external tools (MS Excel, OpenOffice, Notepad)</td>
    </tr>
   </table>
  </CENTER>

  See the mzTab specification for details on the format.

  @experimental This algorithm and underlying format is work in progress and might change.

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
    TOPPBase("MzTabExporter", "Exports various XML formats to an mzTab file.")
  {
  }

protected:

  typedef map< pair< String, String>, vector<PeptideHit> > MapAccPepType;

  void registerOptionsAndFlags_()
  {
    registerInputFile_("in", "<file>", "", "Input file ");
    //setValidFormats_("in", StringList::create("featureXML,consensusXML,idXML"));
    setValidFormats_("in", StringList::create("idXML,featureXML"));
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
    String mods_string;

    const AASequence& aa_seq = peptide_hit.getSequence();
    bool first = true;

    // check terminal modifications
    if (aa_seq.hasNTerminalModification())
    {
      if ( !first )
      {
        mods_string +=  ",";
      } else
      {
        first = false;
      }
      String position = "0";
      String unimod_name = aa_seq.getNTerminalModification();
      String unimod_accession =  ModificationsDB::getInstance()->getModification(unimod_name).getUniModAccession();
      mods_string += position  + "-" + unimod_accession;
    }

    if (aa_seq.hasCTerminalModification())
    {
      if ( !first )
      {
        mods_string +=  ",";
      } else
      {
        first = false;
      }
      String position = String(aa_seq.size() + 1);
      String unimod_name = aa_seq.getCTerminalModification();
      String unimod_accession =  ModificationsDB::getInstance()->getModification(unimod_name).getUniModAccession();
      mods_string += position  + "-" + unimod_accession;
    }

    // check internal modifications
    for (Size i = 0; i != aa_seq.size(); ++i)
    {
      if ( aa_seq[i].isModified() )
      {
        if ( !first )
        {
          mods_string +=  ",";
        } else
        {
          first = false;
        }
        String position = String(i + 1);
        // find all modifications with the given name (but different residue/term specifity)
        std::set< const ResidueModification * > modis;
        ModificationsDB::getInstance()->searchModifications(modis, aa_seq[i].getModification(), ResidueModification::ANYWHERE);
        if (!modis.empty())
        {
          // all have the same unimod accession (=record_id) so just take the first one
          set<const ResidueModification*>::const_iterator mit = modis.begin();
          String unimod_accession = (*mit)->getUniModAccession();
          mods_string += position  + "-" + unimod_accession.c_str();
        }
      }
    }

    if (mods_string.length() == 0)
    {
      return "--";
    }

    return mods_string;
  }

  static String mapSearchEngineToCvParam_(const String& openms_search_engine_name)
  {
    String s = openms_search_engine_name;
    s.toUpper();

    if (s == "OMSSA")
    {
      return "[MS,MS:1001475,OMSSA,]";
    } else if (s == "MASCOT")
    {
      return "[MS,MS:1001207,MASCOT,]";
    } else if (s == "XTANDEM")
    {
      return "[MS,MS:1001476,xtandem,]";
    } else if (s == "SEQUEST")
    {
      return "[MS,MS:1001208,Sequest,]";
    } else if (s == "COMPNOVO")
    {
      return "[,,CompNovo,]";
    } else if (s == "PROTEINPROPHET")
    {
      return "[,,ProteinProphet,]";
    } else
    {
      return "NA";
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
      s = "[,,Consensus:score,";
    } else if (score_type == "q-value")
    {
      s = "[MS,MS:1001364,pep:global FDR,";
    } else if (score_type == "FDR")
    {
      s = "[MS,MS:1001364,pep:global FDR,";
    } else if (score_type == "Posterior Error Probability")
    {
      s = "[,,PEP,";
    } else if (score_type == "PhosphoScore")
    {
      s = "[,,PhosphoScore,";
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
      s = "[,,CompNovo,";
    } else if (score_type == "ProteinProphet probability")
    {
      s = "[,,ProteinProphet,";
    }
    else
    {
      return "NA";
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

  /// Extract protein and peptide identifications for each run. maps are assumed empty.
  static void partitionIntoRuns(const vector<PeptideIdentification>& pep_ids,
                                const vector<ProteinIdentification>& pro_ids,
                                map<String, vector<PeptideIdentification> >& map_run_to_pepids,
                                map<String, vector<ProteinIdentification> >& map_run_to_proids
                                )
  {
    // IdXML can be merged so we have to deal with several runs
    // Extract protein and peptide identifications for each run
    for (Size i = 0; i != pro_ids.size(); ++i)
    {
      map_run_to_proids[pro_ids[i].getIdentifier()].push_back(pro_ids[i]);
    }

    for (Size i = 0; i != pep_ids.size(); ++i)
    {
      map_run_to_pepids[pep_ids[i].getIdentifier()].push_back(pep_ids[i]);
    }

    // perform some sanity check (run ids should be identical)
    assert(map_run_to_proids.size() == map_run_to_pepids.size());
    map<String, vector<PeptideIdentification> >::const_iterator mpep_it = map_run_to_pepids.begin();
    map<String, vector<ProteinIdentification> >::const_iterator mprot_it = map_run_to_proids.begin();
    for (; mpep_it != map_run_to_pepids.end(); ++mpep_it, ++mprot_it)
    {
      assert(mpep_it->first == mprot_it->first);
    }
  }

  /// create links from protein to peptides
  static void createProteinToPeptideLinks(const map<String, vector<PeptideIdentification> >& map_run_to_pepids, MapAccPepType& map_run_accession_to_pephits)
  {
    // create links for each run
    map< String, vector<PeptideIdentification> >::const_iterator mpep_it = map_run_to_pepids.begin();

    for (; mpep_it != map_run_to_pepids.end(); ++mpep_it)
    {
      const String& run = mpep_it->first;
      const vector<PeptideIdentification>& pids = mpep_it->second;
      for (Size i = 0; i != pids.size(); ++i)
      {
        const std::vector<PeptideHit>& phits = pids[i].getHits();
        if (phits.size() == 1)
        {
          const std::vector<String>& accessions = phits[0].getProteinAccessions();
          for (Size k = 0; k != accessions.size(); ++k)
          {
            const String& accession = accessions[k];
            std::pair<String, String> key = make_pair(run, accession);
            map_run_accession_to_pephits[key].push_back(phits[0]);
          }
        }
      }
    }
  }

  static String extractNumPeptides(const String& common_identifier, const String& protein_accession,
                                   const MapAccPepType& map_run_accesion_to_peptides)
  {
    std::pair<String, String> key = make_pair(common_identifier, protein_accession);
    String ret = "0";
    MapAccPepType::const_iterator it = map_run_accesion_to_peptides.find(key);
    if (it != map_run_accesion_to_peptides.end())
    {
      const std::vector<PeptideHit>& peptide_hits = it->second;
      ret = String(peptide_hits.size());
    }
    return ret;
  }
/*
  static String extractModifications(const String& common_identifier, const String& protein_accession,
                                   const MapAccPepType& map_run_accesion_to_peptides)
  {
    std::pair<String, String> key = make_pair(common_identifier, protein_accession);
    String ret = "0";
    MapAccPepType::const_iterator it = map_run_accesion_to_peptides.find(key);
    if (it != map_run_accesion_to_peptides.end())
    {
      const std::vector<PeptideHit>& peptide_hits = it->second;
      // extract all AASequences that contain modifications for the current protein
      std::set<String> mod_sequences;
      for (vector<PeptideHit>::const_iterator pet = peptide_hits.begin(); pet != peptide_hits.end(); ++pet)
      {
        if (pet->getSequence().isModified())
        {
          mod_sequences.insert(pet->getSequence().toString()); // AASequence including Modifications
        }
      }
    }
    return ret;
  }
*/
  // mzTab definition of distinct
  static String extractNumPeptidesDistinct(String common_identifier, String protein_accession,
                                           const MapAccPepType& map_run_accesion_to_peptides)
  {
    std::pair<String, String> key = make_pair(common_identifier, protein_accession);
    String ret = "0";
    MapAccPepType::const_iterator it = map_run_accesion_to_peptides.find(key);
    if (it != map_run_accesion_to_peptides.end())
    {
      const std::vector<PeptideHit>& peptide_hits = it->second;

      // mzTab unambigous peptides are all peptides with different AA sequence OR Modifications
      std::set<String> sequences;
      for (vector<PeptideHit>::const_iterator pet = peptide_hits.begin(); pet != peptide_hits.end(); ++pet)
      {
        sequences.insert(pet->getSequence().toString()); // AASequence including Modifications
      }

      ret = String(sequences.size());
    }

    return ret;
  }

  // same as distinct but additional constraint of uniquenes (=maps to exactly one Protein)
  static String extractNumPeptidesUnambiguous(String common_identifier, String protein_accession,
                                              const MapAccPepType& map_run_accesion_to_peptides)
  {
    std::pair<String, String> key = make_pair(common_identifier, protein_accession);
    String ret = "0";
    MapAccPepType::const_iterator it = map_run_accesion_to_peptides.find(key);
    if (it != map_run_accesion_to_peptides.end())
    {
      const std::vector<PeptideHit>& peptide_hits = it->second;

      // mzTab unambigous peptides are all peptides with different AA sequence OR Modifications
      std::set<String> sequences;
      for (vector<PeptideHit>::const_iterator pet = peptide_hits.begin(); pet != peptide_hits.end(); ++pet)
      {
        // only add sequences of unique peptides
        if (pet->getProteinAccessions().size() == 1)
        {
          sequences.insert(pet->getSequence().toString()); // AASequence with Modifications
        }
      }
      ret = String(sequences.size());
    }
    return ret;
  }

  static map<String, Size> extractNumberOfSubSamples_(const map<String, vector<ProteinIdentification> >& map_run_to_proids)
  {
    map<String, set<Size> > map_run_to_subsamples_id;

    // for each run...
    for (map<String, vector<ProteinIdentification> >::const_iterator run_it = map_run_to_proids.begin();
         run_it != map_run_to_proids.end(); ++run_it)
    {
      String run = run_it->first;
      const vector<ProteinIdentification>& protein_ids = run_it->second;
      // note: per run there should only exist one protein identification
      for ( vector<ProteinIdentification>::const_iterator prot_it = protein_ids.begin();
            prot_it != protein_ids.end(); ++prot_it)
      {               
        const ProteinIdentification& protein_id = *prot_it;
        const vector< ProteinHit >& protein_hits = protein_id.getHits();
        // for each ProteinHit...
        for ( vector< ProteinHit >::const_iterator pit = protein_hits.begin(); pit != protein_hits.end(); ++pit)
        {
          vector< String > metainfo_keys;
          pit->getKeys(metainfo_keys);
          // find meta values starting with mzTab:protein_abundance_sub
          for ( vector<String>::const_iterator s_it = metainfo_keys.begin(); s_it != metainfo_keys.end(); ++s_it)
          {
            //cout << *s_it << endl;
            if (s_it->hasPrefix("mzTab:protein_abundance_sub"))
            {
              String s = *s_it;
              s = s.substitute("mzTab:protein_abundance_sub", "").remove('[').remove(']');
              Size subsample_number = (Size) s.toInt();
              map_run_to_subsamples_id[run].insert(subsample_number);
            }
          }
        }
      }
    }

    // count and return subsample set sizes
    map<String, Size> map_run_to_nsubsamples;
    for ( map<String, set<Size> >::const_iterator sub_it = map_run_to_subsamples_id.begin();
          sub_it != map_run_to_subsamples_id.end(); ++sub_it)
    {
      map_run_to_nsubsamples[sub_it->first] = sub_it->second.size();
    }

    return map_run_to_nsubsamples;
  }

  static void writePeptideHeader_( SVOutStream& output, map<String, Size> n_sub_samples)
  {
    output << "PEH" << "sequence" << "accession" << "unit_id" << "unique" << "database"
           << "database_version" << "search_engine" << "search_engine_score"
           << "modifications" << "retention_time" << "charge"
           << "mass_to_charge" << "uri" << "spectra_ref";

    // to generate sufficient number of columns the maximum of sub samples in all runs is used
    Size max_subsamples = 0;
    for (map<String, Size>::const_iterator run_it = n_sub_samples.begin();
         run_it != n_sub_samples.end(); ++run_it)
    {
      if (run_it->second > max_subsamples)
      {
        max_subsamples = run_it->second;
      }
    }

    // print column headers
    for (Size i = 1; i <= max_subsamples; ++i)
    {
      output << String("peptide_abundance_sub[") + String(i) + String("]") <<
                String("peptide_abundance_stdev_sub[") + String(i) + String("]") <<
                String("peptide_abundance_std_error_sub[") + String(i) + String("]");
    }

    output << endl;
  }

  static void writeProteinHeader_( SVOutStream& output, map<String, Size> n_sub_samples)
  {
    output << "PRH" << "accession" << "unit_id" << "description" << "taxid"
           << "species" << "database" << "database_version" << "search_engine"
           << "search_engine_score" << "reliability" << "num_peptides" << "num_peptides_distinct"
           << "num_peptides_unambiguous" << "ambiguity_members" << "modifications" << "uri"
           << "go_terms" << "protein_coverage";

    // to generate sufficient number of columns the maximum of sub samples in all runs is used
    Size max_subsamples = 0;
    for (map<String, Size>::const_iterator run_it = n_sub_samples.begin();
         run_it != n_sub_samples.end(); ++run_it)
    {
      if (run_it->second > max_subsamples)
      {
        max_subsamples = run_it->second;
      }
    }

    // print column headers
    for (Size i = 1; i <= max_subsamples; ++i)
    {
      output << String("protein_abundance_sub[") + String(i) + String("]") <<
                String("protein_abundance_stdev_sub[") + String(i) + String("]") <<
                String("protein_abundance_std_error_sub[") + String(i) + String("]");
    }

    output << endl;
  }

  static void writeProteinData_(SVOutStream& output,
                                const ProteinIdentification& prot_id,
                                Size run_count,
                                String input_filename,
                                bool has_coverage,
                                const MapAccPepType& map_run_accesion_to_peptides,
                                const map<String, Size>& map_run_to_num_sub
                                )
  {
    // TODO: maybe save these ProteinIdentification run properties in meta data
    // it->getScoreType()
    // it->isHigherScoreBetter())
    // it->getDateTime().toString(Qt::ISODate).toStdString()
    // it->getSearchEngineVersion();

    // search parameters
    const ProteinIdentification::SearchParameters& sp = prot_id.getSearchParameters();
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
    String UNIT_ID_String = File::basename(input_filename) + "-" + String(run_count);
    String database_String = (sp.db != "" ? sp.db : "--");
    database_String = "file://" +  database_String;
    String database_version_String = (sp.db_version != "" ? sp.db_version : "--");
    String species_String =  (sp.taxonomy == "0" || sp.taxonomy == ""  ? "--" : sp.taxonomy);
    String search_engine_cvParams = mapSearchEngineToCvParam_(prot_id.getSearchEngine());
    String openms_search_engine_name = prot_id.getSearchEngine();
    //
    for (vector<ProteinHit>::const_iterator protein_hit_it = prot_id.getHits().begin();
         protein_hit_it != prot_id.getHits().end(); ++protein_hit_it)
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
                                                                  prot_id.getScoreType());

      String reliability = "--";
      String num_peptides;
      String num_peptides_distinct = extractNumPeptidesDistinct(prot_id.getIdentifier(), accession, map_run_accesion_to_peptides);
      String num_peptides_unambiguous = extractNumPeptidesUnambiguous(prot_id.getIdentifier(), accession, map_run_accesion_to_peptides);
      String ambiguity_members = "NA";  //TODO
      String modifications = "NA"; // TODO
      String uri = "file://" + input_filename;
      String go_terms = "--";
      String protein_coverage;

      if (has_coverage)
      {
        protein_coverage = String(protein_hit_it->getCoverage() / 100.0);
      } else
      {
        protein_coverage = "NA";
      }

      if ( protein_hit_it->metaValueExists("num_peptides") )
      {
        num_peptides = protein_hit_it->getMetaValue("num_peptides");
      } else
      {
        num_peptides = extractNumPeptides(prot_id.getIdentifier(), accession, map_run_accesion_to_peptides);
      }

      output << "PRT" << accession << unit_id << description << taxid
             << species << database << database_version << search_engine
             << search_engine_score << reliability << num_peptides << num_peptides_distinct
             << num_peptides_unambiguous << ambiguity_members << modifications << uri
             << go_terms << protein_coverage;

      // get number of sub samples for this run
      map<String, Size>::const_iterator sub_it = map_run_to_num_sub.find(prot_id.getIdentifier());
      Size n_subsamples = sub_it->second;
      for (Size n = 1; n <= n_subsamples; ++n)
      {
        {
          String key = "mzTab:protein_abundance_sub[" + String(n) + "]";
          String abundancy_value = "--";
          if (protein_hit_it->metaValueExists(key))
          {
            abundancy_value = protein_hit_it->getMetaValue(key);
          }
          output << abundancy_value;
        }
        {
          String key = "mzTab:protein_abundance_stdev_sub[" + String(n) + "]";
          String abundancy_value = "--";
          if (protein_hit_it->metaValueExists(key))
          {
            abundancy_value = protein_hit_it->getMetaValue(key);
          }
          output << abundancy_value;
        }
        {
          String key = "mzTab:protein_abundance_std_error_sub[" + String(n) + "]";
          String abundancy_value = "--";
          if (protein_hit_it->metaValueExists(key))
          {
            abundancy_value = protein_hit_it->getMetaValue(key);
          }
          output << abundancy_value;
        }
      }
      output << endl;
    }
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

    // create output stream for MzTab file
    ofstream txt_out( out.c_str() );
    SVOutStream output( txt_out, "\t", "_", String::NONE ); // separator according to MzTab specification

    map<String, vector<PeptideIdentification> > map_run_to_pepids;
    map<String, vector<ProteinIdentification> > map_run_to_proids;

    bool has_coverage = true;

    if (in_type == FileTypes::IDXML)
    {
      String document_id;
      Size num_runs = 0;
      {
        vector<ProteinIdentification> prot_ids;
        vector<PeptideIdentification> pep_ids;

        IdXMLFile().load( in, prot_ids, pep_ids, document_id );

        // pre-filter for best PSM
        bool report_all_psm = getFlag_("report_all_psm");
        IDFilter id_filter;
        for (vector<PeptideIdentification>::iterator pep_id_it = pep_ids.begin();
             pep_id_it != pep_ids.end(); ++pep_id_it)
        {
          pep_id_it->assignRanks();
          if (!report_all_psm)
          {
            PeptideIdentification new_pep_id;
            id_filter.filterIdentificationsByBestHits(*pep_id_it, new_pep_id, false);
            pep_id_it->setHits(new_pep_id.getHits());
          }
        }

        // warn if no coverage information is available
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

        // partition into runs (as this could be a merged IdXML)
        partitionIntoRuns(pep_ids, prot_ids, map_run_to_pepids, map_run_to_proids);
        num_runs = map_run_to_pepids.size();
        cout << "IdXML contains: " << num_runs << " runs." << endl;
      }

      MapAccPepType map_run_accesion_to_peptides;
      createProteinToPeptideLinks(map_run_to_pepids, map_run_accesion_to_peptides);

      // every ProteinIdentification corresponds to a search engine run and contains protein hits

      // write meta data for each run
      bool meta_info_printed = false;
      Size run_count = 0;
      map<String, vector<ProteinIdentification> >::const_iterator mprot_it = map_run_to_proids.begin();
      map<String, vector<PeptideIdentification> >::const_iterator mpep_it = map_run_to_pepids.begin();

      for (; mprot_it != map_run_to_proids.end(); ++mprot_it) // iterate over runs
      {
        // extract ProteinIdentifications of this run (should be only 1)
        const vector<ProteinIdentification>& prot_ids = mprot_it->second;
        for (vector<ProteinIdentification>::const_iterator it = prot_ids.begin(); it != prot_ids.end(); ++it)
        {
          String UNIT_ID = File::basename(in) + "-" + String(run_count);
          String title = document_id;
          if ( title != "" )
          {
            output << "MOD" << UNIT_ID + "-title" << title;
            meta_info_printed = true;
          }
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

      // determine the number of sub samples in each run from protein ids (it is assumed that peptide ids don't introduce new sub sample categories)
      map<String, Size> map_run_to_n_subsamples = extractNumberOfSubSamples_(map_run_to_proids);
      writeProteinHeader_(output, map_run_to_n_subsamples);

      // write protein table data
      run_count = 0;     
      mprot_it = map_run_to_proids.begin();
      mpep_it = map_run_to_pepids.begin();
      for (; mprot_it != map_run_to_proids.end(); ++mprot_it, ++mpep_it) // iterate over runs
      {
        // extract ProteinIdentifications of this run (should be only 1)
        const vector<ProteinIdentification>& prot_ids = mprot_it->second;
        for (vector<ProteinIdentification>::const_iterator prot_id_it = prot_ids.begin();
             prot_id_it != prot_ids.end(); ++prot_id_it, ++run_count)
        {
          writeProteinData_(output, *prot_id_it, run_count, in, has_coverage, map_run_accesion_to_peptides, map_run_to_n_subsamples);
        }
      }

      // write peptide header
      output << endl;

      writePeptideHeader_(output, map_run_to_n_subsamples);

      run_count = 0;
      mprot_it = map_run_to_proids.begin();
      mpep_it = map_run_to_pepids.begin();
      for (; mpep_it != map_run_to_pepids.end(); ++mprot_it, ++mpep_it, ++run_count) // iterate over runs
      {
        // extract ProteinIdentifications of this run (should be only 1)
        const vector<ProteinIdentification>& prot_ids = mprot_it->second;
        // extract PeptideIdentifications of this run
        const vector<PeptideIdentification>& pep_ids = mpep_it->second;
        // iterate over runs of peptide identifications
        for (vector<PeptideIdentification>::const_iterator pep_id_it = pep_ids.begin();
             pep_id_it != pep_ids.end(); ++pep_id_it)
        {
          // TODO: check if bad design of Protein/PeptideIdentification as search engine parameters are stored in prot.
          String openms_search_engine_name = prot_ids[0].getSearchEngine();
          String search_engine_cvParams = mapSearchEngineToCvParam_(openms_search_engine_name);

          const ProteinIdentification::SearchParameters& sp = prot_ids[0].getSearchParameters();
          String UNIT_ID_String = File::basename(in) + "-" + String(run_count);
          String database_String = (sp.db != "" ? sp.db : "--");
          database_String = "file://" + database_String;
          String database_version_String = (sp.db_version != "" ? sp.db_version : "--");

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

            String spectra_ref = "--";

            if ( peptide_hit_it->metaValueExists("mzTab:unique") )
            {
              bool is_unique = peptide_hit_it->getMetaValue("mzTab:unique").toBool();
              if ( is_unique )
              {
                unique = "1";
              } else
              {
                unique = "0";
              }
            } else // no uniqueness annotation
            {
              // if unique protein is present peptide can be assigned
              if (peptide_hit_it->getProteinAccessions().size() == 1)
              {
                unique = "1";
              } else
              {
                unique = "0";
              }
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

            String charge = "NA";
            if (peptide_hit_it->getCharge() != 0)
            {
              charge = peptide_hit_it->getCharge();
            }

            String uri = "file://" + in;

            output << "PEP" << sequence << accession << unit_id << unique << database
                   << database_version << search_engine << search_engine_score
                   << modifications << retention_time << charge
                   << mass_to_charge << uri << spectra_ref;

            // get number of sub samples for this run
            map<String, Size>::const_iterator sub_it = map_run_to_n_subsamples.find(mpep_it->first);
            Size n_subsamples = sub_it->second;

            for (Size n = 1; n <= n_subsamples; ++n)
            {
              {
                String key = "mzTab:peptide_abundance_sub[" + String(n) + "]";
                String abundancy_value = "--";
                if (peptide_hit_it->metaValueExists(key))
                {
                  abundancy_value = peptide_hit_it->getMetaValue(key);
                }
                output << abundancy_value;
              }
              {
                String key = "mzTab:peptide_abundance_stdev_sub[" + String(n) + "]";
                String abundancy_value = "--";
                if (peptide_hit_it->metaValueExists(key))
                {
                  abundancy_value = peptide_hit_it->getMetaValue(key);
                }
                output << abundancy_value;
              }
              {
                String key = "mzTab:peptide_abundance_std_error_sub[" + String(n) + "]";
                String abundancy_value = "--";
                if (peptide_hit_it->metaValueExists(key))
                {
                  abundancy_value = peptide_hit_it->getMetaValue(key);
                }
                output << abundancy_value;
              }
            }

            output << endl;

          }
        }
      }
      txt_out.close();
    } else if (in_type == FileTypes::FEATUREXML)
    {
      FeatureMap<> feature_map;
      FeatureXMLFile f_file;
      f_file.load(in, feature_map);
      vector<ProteinIdentification> protein_ids = feature_map.getProteinIdentifications();
      vector<PeptideIdentification> unassigned_peptides = feature_map.getUnassignedPeptideIdentifications();
      vector<PeptideIdentification> assigned_peptides;

      // 1. sort PeptideHits of peptide identification (assigned to a feature)
      // 2. retain best PeptideHit of each peptide identification
      for (FeatureMap<>::iterator fit = feature_map.begin(); fit != feature_map.end(); ++fit)
      {
        vector<PeptideIdentification>& f_pids = fit->getPeptideIdentifications();
        for (vector<PeptideIdentification>::iterator pid_it = f_pids.begin(); pid_it != f_pids.end(); ++pid_it)
        {
          pid_it->assignRanks();
          // only retain best hit of each peptide identification
          if (pid_it->getHits().size() > 1)
          {
            vector<PeptideHit> phts(pid_it->getHits().begin(), pid_it->getHits().begin() + 1);
            pid_it->setHits(phts);
            cerr << "Warning: Feature contains peptide identification with "
                 << pid_it->getHits().size() << " hits. Only best hit is retained." << endl;
          }
        }
      }

      // extract assigned PeptideIdentifications
      for (FeatureMap<>::iterator fit = feature_map.begin(); fit != feature_map.end(); ++fit)
      {
        vector<PeptideIdentification>& f_pids = fit->getPeptideIdentifications();
        for (vector<PeptideIdentification>::iterator pid_it = f_pids.begin(); pid_it != f_pids.end(); ++pid_it)
        {
          if (pid_it->getHits().size() == 1)
          {
            assigned_peptides.push_back(*pid_it);
          }
        }
      }

      // warn if no coverage information can be calculated
      // Note: For Proteins coverage is calculated based on the assigned peptides
      try
      { // might throw Exception::MissingInformation() if no protein sequence information is added
        for (Size i = 0; i < protein_ids.size(); ++i)
        {
          protein_ids[i].computeCoverage(assigned_peptides);
        }
      }
      catch (Exception::MissingInformation& e)
      {
        LOG_WARN << e.what() << "\n";
        has_coverage = false;
      }

      // partition into runs (for featureXML this usually only contain one run)
      partitionIntoRuns(assigned_peptides, protein_ids, map_run_to_pepids, map_run_to_proids);
      Size num_runs = map_run_to_pepids.size();
      cout << "FeatureXML contains: " << num_runs << " runs." << endl;

      // write meta data for each run
      bool meta_info_printed = false;
      Size run_count = 0;

      // TODO: add meta-info

      if (meta_info_printed)
      {
        output << endl;
      }

      MapAccPepType map_run_accesion_to_peptides;
      createProteinToPeptideLinks(map_run_to_pepids, map_run_accesion_to_peptides);

      // write protein table
      map<String, Size> map_run_to_nsubsamples = extractNumberOfSubSamples_(map_run_to_proids);
      writeProteinHeader_(output, map_run_to_nsubsamples);

      for (vector<ProteinIdentification>::const_iterator prot_id_it = protein_ids.begin();
           prot_id_it != protein_ids.end(); ++prot_id_it, ++run_count)
      {
        writeProteinData_(output, *prot_id_it, run_count, in, has_coverage, map_run_accesion_to_peptides, map_run_to_nsubsamples);
      }

      // writePeptideHeader_(output, n_sub_samples);
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
