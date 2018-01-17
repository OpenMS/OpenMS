//! [final]

// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2016.
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
// $Maintainer: Oliver Alka $
// $Authors: Oliver Alka $
// --------------------------------------------------------------------------

//! [Includes]

#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/FORMAT/MzIdentMLFile.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>

//! [Includes]

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
   @page UTILS_DatabaseFilter DatabaseFilter

   @brief The DatabaseFilter tool filters a protein database in fasta format according to one or multiple filtering criteria.

   The resulting database is written as output. Depending on the reporting method (method="whitelist" or "blacklist") only entries are retained that passed all filters ("whitelist) or failed at least one filter ("blacklist").

   Implemented filter criteria:

       ID: Filter database according to the set of proteinIDs contained in an identification file (idXML, mzIdentML)

   <B>The command line parameters of this tool are:</B>
   @verbinclude UTILS_DatabaseFilter.cli
   <B>INI file documentation of this tool:</B>
   @htmlinclude UTILS_DatabaseFilter.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPDatabaseFilter :
 public TOPPBase
{
public:
 TOPPDatabaseFilter() :
   TOPPBase("DatabaseFilter", "Filters a protein database (FASTA format) based on identified proteins", false)
 {
 }

protected:

//! [Register]

 void registerOptionsAndFlags_()
 {
   registerInputFile_("in", "<file>", "","Input FASTA file, containing a database.");
   setValidFormats_("in", ListUtils::create<String>("fasta"));
   registerInputFile_("id", "<file>", "", "Input file containing identified peptides and proteins.");
   setValidFormats_("id", ListUtils::create<String>("idXML,mzid"));
   registerStringOption_("method", "<choice>", "whitelist", "Switch between white-/blacklisting of protein IDs", false);
   setValidStrings_("method", ListUtils::create<String>("whitelist,blacklist"));
   registerOutputFile_("out", "<file>", "", "Output FASTA file where the reduced database will be written to.");
   setValidFormats_("out", ListUtils::create<String>("fasta"));
 }

//! [Register]

//! [Functionality_1]

 void filterByProteinIDs_(const vector<FASTAFile::FASTAEntry>& db, const vector<PeptideIdentification>& peptide_identifications, bool whitelist, vector<FASTAFile::FASTAEntry>& db_new)
 {
   set<String> id_accessions;
   for (Size i = 0; i != peptide_identifications.size(); ++i)
   {
     const PeptideIdentification& id = peptide_identifications[i];
     const vector<PeptideHit>& hits = id.getHits();
     for (Size k = 0; k != hits.size(); ++k)
     {
       const vector<PeptideEvidence>& evidences = hits[k].getPeptideEvidences();
       for (Size m = 0; m != evidences.size(); ++m)
       {
         const String& id_accession = evidences[m].getProteinAccession();
         id_accessions.insert(id_accession);
       }
     }
   }

 //! [Functionality_1]

   LOG_INFO << "Protein IDs: " << id_accessions.size() << endl;

 //! [Functionality_2]

   for (Size i = 0; i != db.size() ; ++i)
   {
     const String& fasta_accession = db[i].identifier;
     const bool found = id_accessions.find(fasta_accession) != id_accessions.end();
     if ((found && whitelist) || (!found && !whitelist)) //either found in the whitelist or not found in the blacklist
     {
       db_new.push_back(db[i]);
     }
   }

 //! [Functionality_2]
 }

 ExitCodes main_(int, const char **)
 {

   //! [InputParam]

   //-------------------------------------------------------------
   // parsing parameters
   //-------------------------------------------------------------
   String in(getStringOption_("in"));
   String ids(getStringOption_("id"));
   String method(getStringOption_("method"));
   bool whitelist = (method == "whitelist");
   String out(getStringOption_("out"));

   //! [InputParam]

   //-------------------------------------------------------------
   // reading input
   //-------------------------------------------------------------

   //! [InputRead]

   vector<FASTAFile::FASTAEntry> db;
   FASTAFile().load(in, db);

   //! [InputRead]

   // Check if no filter criteria was given
   // If you add a new filter please check if it was set here as well
   if (ids.empty())
   {
     FASTAFile().store(out, db);
   }

   vector<FASTAFile::FASTAEntry> db_new;

   if (!ids.empty()) // filter by protein IDs
   {
     FileHandler fh;
     FileTypes::Type ids_type = fh.getType(ids);
     vector<ProteinIdentification> protein_identifications;
     vector<PeptideIdentification> peptide_identifications;

     if (ids_type == FileTypes::IDXML)
     {
       IdXMLFile().load(ids, protein_identifications, peptide_identifications);
     }
     else if (ids_type == FileTypes::MZIDENTML)
     {
       MzIdentMLFile().load(ids, protein_identifications, peptide_identifications);
     }
     else
     {
       writeLog_("Error: Unknown input file type given. Aborting!");
       printUsage_();
       return ILLEGAL_PARAMETERS;
     }

     LOG_INFO << "Identifications: " << ids.size() << endl;

     // run filter
     filterByProteinIDs_(db, peptide_identifications, whitelist, db_new);
   }

   //-------------------------------------------------------------
   // writing output
   //-------------------------------------------------------------

   LOG_INFO << "Database entries (before / after): " << db.size() << " / " << db_new.size() << endl;
   //! [output]  

   FASTAFile().store(out, db_new);

   //! [output]

   return EXECUTION_OK;
 }

};

int main(int argc, const char ** argv)
{
 TOPPDatabaseFilter tool;
 return tool.main(argc, argv);
}

/// @endcond

//! [final]

