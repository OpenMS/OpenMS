//! [doxygen_snippet_TOPPexample]
// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Oliver Alka $
// $Authors: Oliver Alka $
// This file is ONLY used for code snippets in the developer tutorial
// --------------------------------------------------------------------------

//! [doxygen_snippet_Includes]

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

//! [doxygen_snippet_Includes]

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------

/**
   @page TOPP_DatabaseFilter DatabaseFilter

   @brief The DatabaseFilter tool filters a protein database in fasta format according to one or multiple filtering criteria.

   The resulting database is written as output. Depending on the reporting method (method="whitelist" or "blacklist") only entries are retained that
   passed all filters ("whitelist) or failed at least one filter ("blacklist").

   Implemented filter criteria:

       ID: Filter database according to the set of proteinIDs contained in an identification file (idXML, mzIdentML)

   <B>The command line parameters of this tool are:</B>
   @verbinclude TOPP_DatabaseFilter.cli
   <B>INI file documentation of this tool:</B>
   @htmlinclude TOPP_DatabaseFilter.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPDatabaseFilter : public TOPPBase
{
public:
  TOPPDatabaseFilter():
      TOPPBase("DatabaseFilter", "Filters a protein database (FASTA format) based on identified proteins", false) // false: mark as unofficial tool
  {
  }

protected:
  //! [doxygen_snippet_Register]

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

  //! [doxygen_snippet_Register]

  //! [doxygen_snippet_Functionality_1]

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

    //! [doxygen_snippet_Functionality_1]

    OPENMS_LOG_INFO << "Number of Protein IDs: " << id_accessions.size() << endl;

    //! [doxygen_snippet_Functionality_2]

    for (const auto entry : db)
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

  ExitCodes main_(int, const char**) override
  {

    //! [doxygen_snippet_InputParam]

    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------
    String in(getStringOption_("in"));
    String ids(getStringOption_("id"));
    String method(getStringOption_("method"));
    bool whitelist = (method == "whitelist");
    String out(getStringOption_("out"));

    //! [doxygen_snippet_InputParam]

    //-------------------------------------------------------------
    // reading input
    //-------------------------------------------------------------

    //! [doxygen_snippet_InputRead]

    vector<FASTAFile::FASTAEntry> db;
    FASTAFile().load(in, db);

    //! [doxygen_snippet_InputRead]

    vector<ProteinIdentification> protein_identifications;
    vector<PeptideIdentification> peptide_identifications;

    FileHandler().loadIdentifications(ids, protein_identifications, peptide_identifications);

    OPENMS_LOG_INFO << "Identifications: " << ids.size() << endl;

    // run filter
    vector<FASTAFile::FASTAEntry> db_new;
    filterByProteinAccessions_(db, peptide_identifications, whitelist, db_new);

    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------

    OPENMS_LOG_INFO << "Database entries (before / after): " << db.size() << " / " << db_new.size() << endl;
    //! [doxygen_snippet_output]

    FASTAFile().store(out, db_new);

    //! [doxygen_snippet_output]

    return EXECUTION_OK;
  }
};

int main(int argc, const char** argv)
{
  TOPPDatabaseFilter tool;
  OPENMS_LOG_FATAL_ERROR << "THIS IS TEST CODE AND SHOULD NEVER BE RUN OUTSIDE OF TESTING" << endl;
  tool.main(argc, argv);
  return 0;
}

/// @endcond

//! [doxygen_snippet_TOPPexample]
