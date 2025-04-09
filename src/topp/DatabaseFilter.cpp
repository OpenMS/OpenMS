// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Oliver Alka $
// $Authors: Oliver Alka $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_DatabaseFilter DatabaseFilter

@brief The DatabaseFilter tool filters a protein database in fasta format according to one or multiple filtering criteria.

The resulting database is written as output. Depending on the reporting method (method="whitelist" or "blacklist") only entries are retained that passed all filters ("whitelist) or failed at least one filter ("blacklist").

Implemented filter criteria:

    accession: Filter database according to the set of protein accessions contained in an identification file (idXML, mzIdentML)

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_DatabaseFilter.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_DatabaseFilter.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPDatabaseFilter :
  public TOPPBase
{
public:
  TOPPDatabaseFilter() :
    TOPPBase("DatabaseFilter", "Filters a protein database (FASTA format) based on identified proteins")
  {
  }

protected:
  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "Input FASTA file, containing a protein database.");
    setValidFormats_("in", {"fasta"});
    registerInputFile_("id", "<file>", "", "Input file containing identified peptides and proteins.");
    setValidFormats_("id", {"idXML", "mzid"});
    registerStringOption_("method", "<choice>", "whitelist", "Switch between white-/blacklisting of protein IDs", false);
    setValidStrings_("method", {"whitelist", "blacklist"});
    registerOutputFile_("out", "<file>", "", "Output FASTA file where the reduced database will be written to.");
    setValidFormats_("out", {"fasta"});
  }

  void filterByProteinAccessions_(const vector<FASTAFile::FASTAEntry>& db,
                                  const vector<PeptideIdentification>& peptide_identifications,
                                  bool whitelist,
                                  vector<FASTAFile::FASTAEntry>& db_new)
  {
    set<String> id_accessions;
    for (const auto& pep_id : peptide_identifications)
    {
      for (const auto& hit : pep_id.getHits())
      {
        for (const auto& ev : hit.getPeptideEvidences())
        {
          const String& id_accession = ev.getProteinAccession();
          id_accessions.insert(id_accession);
        }
      }
    }

    OPENMS_LOG_INFO << "Number of Protein IDs: " << id_accessions.size() << endl;

    for (const auto& entry : db)
    {
      const String& fasta_accession = entry.identifier;
      const bool found = id_accessions.find(fasta_accession) != id_accessions.end();
      if ((found && whitelist) || (! found && ! whitelist)) // either found in the whitelist or not found in the blacklist
      {
        db_new.push_back(entry);
      }
    }

    //! [doxygen_snippet_Functionality_2]
  }

  ExitCodes main_(int, const char **) override
  {
    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------
    String in(getStringOption_("in"));
    String ids(getStringOption_("id"));
    String method(getStringOption_("method"));
    bool whitelist = (method == "whitelist");
    String out(getStringOption_("out"));

    //-------------------------------------------------------------
    // reading input
    //-------------------------------------------------------------

    vector<FASTAFile::FASTAEntry> db;
    FASTAFile().load(in, db);

    // Check if no filter criteria was given
    // If you add a new filter please check if it was set here as well
    if (ids.empty())
    {
      FASTAFile().store(out, db);
    }

    vector<FASTAFile::FASTAEntry> db_new;

    if (!ids.empty()) // filter by protein accessions in id files
    {
      FileHandler fh;
      FileTypes::Type ids_type = fh.getType(ids);
      vector<ProteinIdentification> protein_identifications;
      vector<PeptideIdentification> peptide_identifications;

      if (ids_type == FileTypes::IDXML || ids_type == FileTypes::MZIDENTML)
      {
        FileHandler().loadIdentifications(ids, protein_identifications, peptide_identifications, {FileTypes::IDXML, FileTypes::MZIDENTML});
      }
      else
      {
        writeLogError_("Error: Unknown input file type given. Aborting!");
        printUsage_();
        return ILLEGAL_PARAMETERS;
      }

      OPENMS_LOG_INFO << "Identifications: " << ids.size() << endl;

      // run filter
      filterByProteinAccessions_(db, peptide_identifications, whitelist, db_new);
    }

    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------

    OPENMS_LOG_INFO << "Database entries (before / after): " << db.size() << " / " << db_new.size() << endl;
    FASTAFile().store(out, db_new);

    return EXECUTION_OK;
  }

};

int main(int argc, const char ** argv)
{
  TOPPDatabaseFilter tool;
  return tool.main(argc, argv);
}

/// @endcond
