// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Andreas Bertsch $
// $Authors: Andreas Bertsch, Knut Reinert $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <algorithm>

#include <OpenMS/DATASTRUCTURES/SeqanIncludeWrapper.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page TOPP_PeptideIndexer PeptideIndexer

	@brief Refreshes the protein references for all peptide hits from a idXML file.

<CENTER>
	<table>
		<tr>
			<td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
			<td VALIGN="middle" ROWSPAN=2> \f$ \longrightarrow \f$ PeptideIndexer \f$ \longrightarrow \f$</td>
			<td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
		</tr>
		<tr>
			<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_IDFilter or @n any protein/peptide processing tool </td>
			<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_FalseDiscoveryRate </td>
		</tr>
	</table>
</CENTER>

	Each peptide hit is annotated by a target_decoy string,
	indicating if the peptide sequence is found in a 'target', a 'decoy' or in both 'target+decoy' protein. This information is
	crucial for the @ref TOPP_FalseDiscoveryRate tool.

	<B>The command line parameters of this tool are:</B>
	@verbinclude TOPP_PeptideIndexer.cli
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPPeptideIndexer
	: public TOPPBase
{
	public:
		TOPPPeptideIndexer()
			: TOPPBase("PeptideIndexer","Refreshes the protein references for all peptide hits.", false)
		{

		}

	protected:
		void registerOptionsAndFlags_()
		{
			registerInputFile_("in","<file>","","Input idXML file containing the identifications.");
			setValidFormats_("in", StringList::create("IdXML"));
			registerInputFile_("fasta", "<file>", "", "Input sequence database in fasta format.");
			registerOutputFile_("out","<file>","","Output idXML file.");
			setValidFormats_("in", StringList::create("IdXML"));
			registerStringOption_("decoy_string", "<string>", "_rev", "String that was appended to the accession of the protein database to indicate a decoy protein.", false);
			registerFlag_("write_protein_sequence", "If set, the protein sequences are added to the protein hits.");
			registerFlag_("keep_unreferenced_proteins", "If set, protein hits which are not referenced by any peptide are kept.");
		}

		ExitCodes main_(int , const char**)
		{
			//-------------------------------------------------------------
			// parsing parameters
			//-------------------------------------------------------------
			String in(getStringOption_("in"));
			String fasta(getStringOption_("fasta"));
			String out(getStringOption_("out"));
			bool write_protein_sequence(getFlag_("write_protein_sequence"));
			bool keep_unreferenced_proteins(getFlag_("keep_unreferenced_proteins"));
			String decoy_string(getStringOption_("decoy_string"));

			//-------------------------------------------------------------
			// reading input
			//-------------------------------------------------------------

			// we stream the Fasta file
			vector<FASTAFile::FASTAEntry> proteins;
			FASTAFile().load(fasta, proteins);

			vector<ProteinIdentification> prot_ids;
			vector<PeptideIdentification> pep_ids;
			IdXMLFile().load(in, prot_ids, pep_ids);

			//-------------------------------------------------------------
			// calculations
			//-------------------------------------------------------------


			writeDebug_("Collecting peptides...", 1);
			// collect the peptides in a Seqan StringSet
			seqan::StringSet<seqan::String<char> > needle;

			// store for each run the protein idx and number of peptides that hit this protein
			Map<String, Map<Size, Size> > prot_idx_hits;

			// map the number of the peptide to the corresponding iterator in vector<PeptideHits>
			Size needle_count = (numeric_limits<Size>::max)();
			Map<String, Size> peptide_to_idx;

			for (vector<PeptideIdentification>::const_iterator it1 = pep_ids.begin(); it1 != pep_ids.end(); ++it1)
			{
				String run_id = it1->getIdentifier();
				vector<PeptideHit> hits = it1->getHits();
				for (vector<PeptideHit>::iterator it2 = hits.begin(); it2 != hits.end(); ++it2)
				{
					it2->setProteinAccessions(vector<String>());
					String seq = it2->getSequence().toUnmodifiedString();
					seqan::appendValue(needle, seq.c_str());

					peptide_to_idx[seq] = ++needle_count;
				}
			}


			// read and concatenate all proteins
			seqan::String<char> all_protein_sequences;
			// build map accessions to proteins
			Map<String, vector<Size> > acc_to_prot;
			Size pos(0);
			Map<Size, Size> idx_to_protein; // stores the begin indices of the 'all_protein_sequences' string and the corresponding protein indices (proteins vector)
			vector<Size> protein_idx_vector; // contains all begin indices of the proteins in the 'all_protein_sequences' string
			for (Size i = 0; i != proteins.size(); ++i)
			{
				protein_idx_vector.push_back(pos);
				idx_to_protein[pos] = i;
				pos += proteins[i].sequence.size() + 1; // consider the terminating '$'
				all_protein_sequences += (proteins[i].sequence + "$").c_str();

				String acc = proteins[i].identifier;
				if (acc_to_prot.has(acc))
				{
					writeLog_(String("PeptideIndexer: error, identifiers of proteins should by unique to a database, identifier '") + acc + String("' found multiply."));
				}
				acc_to_prot[acc].push_back(i);
			}

			// Aho Corasick Call
			seqan::Finder<seqan::String<char>  > finder(all_protein_sequences);
			seqan::Pattern<seqan::StringSet<seqan::String<char> >, seqan::AhoCorasick > pattern(needle);

			seqan::String<seqan::Pair<Size, Size> > pat_hits;
			Map<Size, vector<Size> > peptide_to_indices;
			writeDebug_("Finding peptide/protein matches...", 1);
			while (find(finder, pattern))
			{
				seqan::appendValue(pat_hits, seqan::Pair<Size, Size>(position(pattern), position(finder)));
				peptide_to_indices[position(pattern)].push_back(position(finder));
			}
			writeDebug_("Ended finding", 1);

			writeDebug_("Reindexing peptide/protein matches...", 1);
			for (vector<PeptideIdentification>::iterator it1 = pep_ids.begin(); it1 != pep_ids.end(); ++it1)
			{
				String run_id = it1->getIdentifier();
				vector<PeptideHit> hits = it1->getHits();
				for (vector<PeptideHit>::iterator it2 = hits.begin(); it2 != hits.end(); ++it2)
				{
					it2->setProteinAccessions(vector<String>());
					String seq = it2->getSequence().toUnmodifiedString();


					//seqan::Pair<Size> hits;
					//hits = seqan::equalRangeSA(all_protein_sequences, suffix_array, seq.c_str());


					for (vector<Size>::const_iterator it = peptide_to_indices[peptide_to_idx[seq]].begin(); it != peptide_to_indices[peptide_to_idx[seq]].end(); ++it)
					{
						vector<Size>::const_iterator lower_bound_iter = lower_bound(protein_idx_vector.begin(), protein_idx_vector.end(), *it) - 1;
						Size prot_idx = idx_to_protein[*lower_bound_iter];
						it2->addProteinAccession(proteins[prot_idx].identifier);

						if (prot_idx_hits.has(run_id))
						{
							if (prot_idx_hits[run_id].has(prot_idx))
							{
								++prot_idx_hits[run_id][prot_idx];
							}
							else
							{
								prot_idx_hits[run_id][prot_idx] = 1;
							}
						}
						else
						{
							prot_idx_hits[run_id][prot_idx] = 1;
						}
					}

					// add information whether this is a decoy hit
					bool matches_target(false);
					bool matches_decoy(false);
					for (vector<String>::const_iterator it = it2->getProteinAccessions().begin(); it != it2->getProteinAccessions().end(); ++it)
					{
						if (it->hasSuffix(decoy_string))
						{
							matches_decoy = true;
						}
						else
						{
							matches_target = true;
						}
					}
					String target_decoy = "";
					if (matches_decoy && !matches_target)
					{
						target_decoy = "decoy";
					}
					if (!matches_decoy && matches_target)
					{
						target_decoy = "target";
					}
					if (matches_decoy && matches_target)
					{
						target_decoy = "target+decoy";
					}
					it2->setMetaValue("target_decoy", target_decoy);
					if (it2->getProteinAccessions().size() == 1)
					{
						it2->setMetaValue("protein_references", "unique");
					}
					else
					{
						if (it2->getProteinAccessions().size() > 1)
						{
							it2->setMetaValue("protein_references", "non-unique");
						}
						else
						{
							it2->setMetaValue("protein_references", "unmatched");
						}
					}
				}
				it1->setHits(hits);
			}

			// all peptides contain the correct protein hit references, now update the protein hits
			vector<ProteinIdentification> new_prot_ids;
			for (vector<ProteinIdentification>::iterator it1 = prot_ids.begin(); it1 != prot_ids.end(); ++it1)
			{
				String run_id = it1->getIdentifier();
				ProteinIdentification new_prot_id = *it1;
				new_prot_id.setHits(vector<ProteinHit>());
				vector<ProteinHit> protein_hits;
				set<String> acc_done;
				// walk through already existing protein hits and update them
				for (vector<ProteinHit>::iterator it2 = it1->getHits().begin(); it2 != it1->getHits().end(); ++it2)
				{
					String acc = it2->getAccession();
					acc_done.insert(acc);
					if (acc_to_prot.has(acc))
					{
						if (write_protein_sequence)
						{
							for (Size i = 0; i != acc_to_prot[acc].size(); ++i)
							{
								it2->setSequence(proteins[acc_to_prot[acc][i]].sequence);
							}
						}
						protein_hits.push_back(*it2);
					}
					else
					{
						if (keep_unreferenced_proteins)
						{
							protein_hits.push_back(*it2);
						}
					}
				}

				// go through new referenced proteins
				if (prot_idx_hits.has(run_id))
				{
					for (Map<Size, Size>::const_iterator it = prot_idx_hits[run_id].begin(); it != prot_idx_hits[run_id].end(); ++it)
					{
						if (it->second == 0)
						{
							continue; // should not happen
						}
						String acc = proteins[it->first].identifier;
						if (acc_done.find(acc) == acc_done.end())
						{
							ProteinHit hit;
							hit.setAccession(acc);
							if(write_protein_sequence)
							{
								hit.setSequence(proteins[it->first].sequence);
							}
							protein_hits.push_back(hit);
						}
					}
				}

				new_prot_id.setHits(protein_hits);
				new_prot_ids.push_back(new_prot_id);
			}
			writeDebug_("Ended reindexing", 1);

			//-------------------------------------------------------------
      // writing output
      //-------------------------------------------------------------

			IdXMLFile().store(out, new_prot_ids, pep_ids);

			return EXECUTION_OK;
		}
};


int main( int argc, const char** argv )
{
	TOPPPeptideIndexer tool;
	return tool.main(argc,argv);
}

/// @endcond





