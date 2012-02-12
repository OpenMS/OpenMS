// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <OpenMS/ANALYSIS/QUANTITATION/PeptideAndProteinQuant.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/SVOutStream.h>

using namespace OpenMS;
using namespace std;


//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page TOPP_ProteinQuantifier ProteinQuantifier
	
	@brief Compute peptide and protein abundances from annotated feature/consensus maps.

<CENTER>
	<table>
		<tr>
			<td ALIGN = "center" BGCOLOR="#EBEBEB"> potential predecessor tools </td>
			<td VALIGN="middle" ROWSPAN=3> \f$ \longrightarrow \f$ ProteinQuantifier \f$ \longrightarrow \f$</td>
			<td ALIGN = "center" BGCOLOR="#EBEBEB"> potential successor tools </td>
		</tr>
		<tr>
			<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_IDMapper </td>
			<td VALIGN="middle" ALIGN = "center" ROWSPAN=2> external tools @n e.g. for statistical analysis</td>
		</tr>
		<tr>
			<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_FeatureLinkerUnlabeled @n (or another feature grouping tool) </td>
		</tr>
	</table>
</CENTER>

	Quantification is based on the intensity values of the features in the input. Feature intensities are first accumulated to peptide abundances, according to the peptide identifications annotated to the features/feature groups. Then, abundances of the peptides of a protein are averaged to compute the protein abundance.

	The peptide-to-protein step uses the (e.g. 3) most abundant proteotypic peptides per protein to compute the protein abundances. This is a general version of the "top 3 approach" (but only for relative quantification) described in:\n
	Silva <em>et al.</em>: Absolute quantification of proteins by LCMS<sup>E</sup>: a virtue of parallel MS acquisition (Mol. Cell. Proteomics, 2006).

	Only features/feature groups with unambiguous peptide annotation are used for peptide quantification, and generally only proteotypic peptides (i.e. those matching to exactly one protein) are used for protein quantification. As an exception to this rule, if ProteinProphet results for the whole sample set are provided with the @p protxml option, or are already included in a featureXML input, also groups of indistinguishable proteins will be quantified. The reported quantity then refers to the total for the whole group.

	Peptide/protein IDs from multiple identification runs can be handled, but will not be differentiated (i.e. protein accessions for a peptide will be accumulated over all identification runs).

	Peptides with the same sequence, but with different modifications are quantified separately on the peptide level, but treated as one peptide for the protein quantification (i.e. the contributions of differently-modified variants of the same peptide are accumulated).

	More information below the parameter specification.

	<B>The command line parameters of this tool are:</B>
	@verbinclude TOPP_ProteinQuantifier.cli

	<B>Output format</B>

	The output files produced by this tool have a table format, with columns as described below:

	<b>Protein output</b> (one protein/set of indistinguishable proteins per line):
	- @b protein: Protein accession(s) (as in the annotations in the input file; separated by "/" if more than one).
	- @b n_proteins: Number of indistinguishable proteins quantified (usually "1").
	- @b protein_score: Protein score, e.g. ProteinProphet probability (if available).
	- @b n_peptides: Number of proteotypic peptides observed for this protein (or group of indistinguishable proteins) across all samples. Note that not necessarily all of these peptides contribute to the protein abundance (depending on parameter @p top).
	- @b abundance: Computed protein abundance. For consensusXML input, there will be one column  per sample ("abundance_1", "abundance_2", etc.).

	<b>Peptide output</b> (one peptide or - if @p filter_charge is set - one charge state of a peptide per line):
	- @b peptide: Peptide sequence. Only peptides that occur in unambiguous annotations of features are reported.
	- @b protein: Protein accession(s) for the peptide (separated by "/" if more than one).
	- @b n_proteins: Number of proteins this peptide maps to. (Same as the number of accessions in the previous column.)
	- @b charge: Charge state quantified in this line. "0" (for "all charges") unless @p filter_charge was set.
	- @b abundance: Computed abundance for this peptide. If the charge in the preceding column is 0, this is the total abundance of the peptide over all charge states; otherwise, it is only the abundance observed for the indicated charge (in this case, there may be more than one line for the peptide sequence). Again, for consensusXML input, there will be one column  per sample ("abundance_1", "abundance_2", etc.). Also for consensusXML, the reported values are already normalized if @p consensus:normalize was set.

	<B>Protein quantification examples</B>

	While quantification on the peptide level is fairly straight-forward, a number of options influence quantification on the protein level - especially for consensusXML input. The three parameters @p top, @p include_all and @p consensus:fix_peptides determine which peptides are used to quantify proteins in different samples.

	As an example, consider a protein with four proteotypic peptides. Each peptide is detected in a subset of three samples, as indicated in the table below. The peptides are ranked by abundance (1: highest, 4: lowest; assuming for simplicity that the order is the same in all samples).

<CENTER>
	<table>
		<tr>
		  <td></td>
			<td ALIGN="center" BGCOLOR="#EBEBEB"> sample 1 </td>
			<td ALIGN="center" BGCOLOR="#EBEBEB"> sample 2 </td>
			<td ALIGN="center" BGCOLOR="#EBEBEB"> sample 3 </td>
		</tr>
		<tr>
			<td ALIGN="center" BGCOLOR="#EBEBEB"> peptide 1 </td>
			<td ALIGN="center"> X </td>
			<td></td>
			<td ALIGN="center"> X </td>
		</tr>
		<tr>
			<td ALIGN="center" BGCOLOR="#EBEBEB"> peptide 2 </td>
			<td ALIGN="center"> X </td>
			<td ALIGN="center"> X </td>
			<td></td>
		</tr>
		<tr>
			<td ALIGN="center" BGCOLOR="#EBEBEB"> peptide 3 </td>
			<td ALIGN="center"> X </td>
			<td ALIGN="center"> X </td>
			<td ALIGN="center"> X </td>
		</tr>
		<tr>
			<td ALIGN="center" BGCOLOR="#EBEBEB"> peptide 4 </td>
			<td ALIGN="center"> X </td>
			<td ALIGN="center"> X </td>
			<td></td>
		</tr>
	</table>
</CENTER>

  Different parameter combinations lead to different quantification scenarios, as shown here:

<CENTER>
	<table>
		<tr>
		  <td ALIGN="center" BGCOLOR="#EBEBEB" COLSPAN=3> @b parameters \n "*": no effect in this case </td>
		  <td ALIGN="center" BGCOLOR="#EBEBEB" COLSPAN=3> <b>peptides used for quantification</b> \n "(...)": not quantified here because ... </td>
		  <td ALIGN="center" VALIGN="middle" BGCOLOR="#EBEBEB" ROWSPAN=2> explanation </td>
		</tr>
		<tr>
			<td ALIGN="center" BGCOLOR="#EBEBEB"> @p top </td>
			<td ALIGN="center" BGCOLOR="#EBEBEB"> @p include_all </td>
			<td ALIGN="center" BGCOLOR="#EBEBEB"> @p c.:fix_peptides </td>
			<td ALIGN="center" BGCOLOR="#EBEBEB"> sample 1 </td>
			<td ALIGN="center" BGCOLOR="#EBEBEB"> sample 2 </td>
			<td ALIGN="center" BGCOLOR="#EBEBEB"> sample 3 </td>
		</tr>
		<tr>
			<td ALIGN="center"> 0 </td>
			<td ALIGN="center"> * </td>
			<td ALIGN="center"> no </td>
			<td ALIGN="center"> 1, 2, 3, 4 </td>
			<td ALIGN="center"> 2, 3, 4 </td>
			<td ALIGN="center"> 1, 3 </td>
			<td> all peptides </td>
		</tr>
		<tr>
			<td ALIGN="center"> 1 </td>
			<td ALIGN="center"> * </td>
			<td ALIGN="center"> no </td>
			<td ALIGN="center"> 1 </td>
			<td ALIGN="center"> 2 </td>
			<td ALIGN="center"> 1 </td>
			<td> single most abundant peptide </td>
		</tr>
		<tr>
			<td ALIGN="center"> 2 </td>
			<td ALIGN="center"> * </td>
			<td ALIGN="center"> no </td>
			<td ALIGN="center"> 1, 2 </td>
			<td ALIGN="center"> 2, 3 </td>
			<td ALIGN="center"> 1, 3 </td>
			<td> two most abundant peptides </td>
		</tr>
		<tr>
			<td ALIGN="center"> 3 </td>
			<td ALIGN="center"> no </td>
			<td ALIGN="center"> no </td>
			<td ALIGN="center"> 1, 2, 3 </td>
			<td ALIGN="center"> 2, 3, 4 </td>
			<td ALIGN="center"> (too few peptides) </td>
			<td> three most abundant peptides </td>
		</tr>
		<tr>
			<td ALIGN="center"> 3 </td>
			<td ALIGN="center"> yes </td>
			<td ALIGN="center"> no </td>
			<td ALIGN="center"> 1, 2, 3 </td>
			<td ALIGN="center"> 2, 3, 4 </td>
			<td ALIGN="center"> 1, 3 </td>
			<td> three or fewer most abundant peptides </td>
		</tr>
		<tr>
			<td ALIGN="center"> 4 </td>
			<td ALIGN="center"> no </td>
			<td ALIGN="center"> * </td>
			<td ALIGN="center"> 1, 2, 3, 4 </td>
			<td ALIGN="center"> (too few peptides) </td>
			<td ALIGN="center"> (too few peptides) </td>
			<td> four most abundant peptides </td>
		</tr>
		<tr>
			<td ALIGN="center"> 4 </td>
			<td ALIGN="center"> yes </td>
			<td ALIGN="center"> * </td>
			<td ALIGN="center"> 1, 2, 3, 4 </td>
			<td ALIGN="center"> 2, 3, 4 </td>
			<td ALIGN="center"> 1, 3 </td>
			<td> four or fewer most abundant peptides </td>
		</tr>
		<tr>
			<td ALIGN="center"> 0 </td>
			<td ALIGN="center"> * </td>
			<td ALIGN="center"> yes </td>
			<td ALIGN="center"> 3 </td>
			<td ALIGN="center"> 3 </td>
			<td ALIGN="center"> 3 </td>
			<td> all peptides present in every sample </td>
		</tr>
		<tr>
			<td ALIGN="center"> 1 </td>
			<td ALIGN="center"> * </td>
			<td ALIGN="center"> yes </td>
			<td ALIGN="center"> 3 </td>
			<td ALIGN="center"> 3 </td>
			<td ALIGN="center"> 3 </td>
			<td> single peptide present in most samples </td>
		</tr>
		<tr>
			<td ALIGN="center"> 2 </td>
			<td ALIGN="center"> no </td>
			<td ALIGN="center"> yes </td>
			<td ALIGN="center"> 1, 3 </td>
			<td ALIGN="center"> (peptide 1 missing) </td>
			<td ALIGN="center"> 1, 3 </td>
			<td> two peptides present in most samples </td>
		</tr>
		<tr>
			<td ALIGN="center"> 2 </td>
			<td ALIGN="center"> yes </td>
			<td ALIGN="center"> yes </td>
			<td ALIGN="center"> 1, 3 </td>
			<td ALIGN="center"> 3 </td>
			<td ALIGN="center"> 1, 3 </td>
			<td> two or fewer peptides present in most samples </td>
		</tr>
		<tr>
			<td ALIGN="center"> 3 </td>
			<td ALIGN="center"> no </td>
			<td ALIGN="center"> yes </td>
			<td ALIGN="center"> 1, 2, 3 </td>
			<td ALIGN="center"> (peptide 1 missing) </td>
			<td ALIGN="center"> (peptide 2 missing) </td>
			<td> three peptides present in most samples </td>
		</tr>
		<tr>
			<td ALIGN="center"> 3 </td>
			<td ALIGN="center"> yes </td>
			<td ALIGN="center"> yes </td>
			<td ALIGN="center"> 1, 2, 3 </td>
			<td ALIGN="center"> 2, 3 </td>
			<td ALIGN="center"> 1, 3 </td>
			<td> three or fewer peptides present in most samples </td>
		</tr>
	</table>
</CENTER>

	<B>Further considerations for parameter selection</B>

	With @p filter_charge and @p average, there is a trade-off between comparability of protein abundances within a sample and of abundances for the same protein across different samples.\n
	Setting @p filter_charge may increase reproducibility between samples, but will distort the proportions of protein abundances within a sample. The reason is that ionization properties vary between peptides, but should remain constant across samples. Filtering by charge state can help to reduce the impact of feature detection differences between samples.\n
	For @p average, there is a qualitative difference between @p mean/median and @p sum in the effect that missing peptide abundances have (only if @p include_all is set or @p top is 0): @p mean and @p median ignore missing cases, averaging only present values. If low-abundant peptides are not detected in some samples, the computed protein abundances for those samples may thus be too optimistic. @p sum implicitly treats missing values as zero, so this problem does not occur and comparability across samples is ensured. However, with @p sum the total number of peptides ("summands") available for a protein may affect the abundances computed for it (depending on @p top), so results within a sample may become unproportional.

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES


class TOPPProteinQuantifier: public TOPPBase
{
public:
		
	TOPPProteinQuantifier(): 
		TOPPBase("ProteinQuantifier", "Compute peptide and protein abundances"),
		algo_params_(), proteins_(), peptides_(), files_() {}

protected:

	typedef PeptideAndProteinQuant::PeptideQuant PeptideQuant;
	typedef PeptideAndProteinQuant::ProteinQuant ProteinQuant;
	typedef PeptideAndProteinQuant::SampleAbundances SampleAbundances;
	typedef PeptideAndProteinQuant::Statistics Statistics;

	Param algo_params_; // parameters for PeptideAndProteinQuant algorithm
	ProteinIdentification proteins_; // ProteinProphet results (proteins)
	PeptideIdentification peptides_; // ProteinProphet results (peptides)
	ConsensusMap::FileDescriptions files_; // Information about files involved

	void registerOptionsAndFlags_()
  {
		registerInputFile_("in", "<file>", "", "Input file");
		setValidFormats_("in", StringList::create("featureXML,consensusXML"));
		registerInputFile_("protxml", "<file>", "", "ProteinProphet results (protXML converted to idXML) for the identification runs that were used to annotate the input.\nInformation about indistinguishable proteins will be used for protein quantification.", false);
		setValidFormats_("protxml", StringList::create("idXML"));
		registerOutputFile_("out", "<file>", "", "Output file for protein abundances", false);
		registerOutputFile_("peptide_out", "<file>", "", "Output file for peptide abundances", false);
		registerOutputFile_("id_out", "<file>", "", "Output file for peptide and protein abundances (annotated idXML) - suitable for export to mzTab.\nEither 'out', 'peptide_out', or 'id_out' are required. They can be used together.", false);
		setValidFormats_("id_out", StringList::create("idXML"));

		// algorithm parameters:
		addEmptyLine_();
		Param temp = PeptideAndProteinQuant().getParameters();
		registerFullParam_(temp);

		registerTOPPSubsection_("format", "Output formatting options");
		registerStringOption_("format:separator", "<sep>", "", "Character(s) used to separate fields; by default, the 'tab' character is used", false);
		registerStringOption_("format:quoting", "<method>", "double", "Method for quoting of strings: 'none' for no quoting, 'double' for quoting with doubling of embedded quotes,\n'escape' for quoting with backslash-escaping of embedded quotes", false);
		setValidStrings_("format:quoting", StringList::create("none,double,escape"));
		registerStringOption_("format:replacement", "<x>", "_", "If 'quoting' is 'none', used to replace occurrences of the separator in strings before writing", false);

		// registerSubsection_("algorithm", "Algorithm parameters section");
	}


	// Param getSubsectionDefaults_(const String& /* section */) const
	// {
	// 	 return PeptideAndProteinQuant().getParameters();
	// }


	/// Write a table of peptide results.
	void writePeptideTable_(SVOutStream& out, const PeptideQuant& quant)
	{
		// write header:
		out << "peptide" << "protein" << "n_proteins" << "charge";
		if (files_.size() <= 1)
		{
			out << "abundance";
		}
		else
		{
			for (Size i = 1; i <= files_.size(); ++i)
			{
				out << "abundance_" + String(i);
			}
		}
		out << endl;

		bool filter_charge = algo_params_.getValue("filter_charge") == "true";
		for (PeptideQuant::const_iterator q_it = quant.begin(); q_it != quant.end();
				 ++q_it)
		{
			if (q_it->second.total_abundances.empty()) continue; // not quantified

			StringList accessions;
			for (set<String>::const_iterator acc_it = 
						 q_it->second.accessions.begin(); acc_it != 
						 q_it->second.accessions.end(); ++acc_it)
			{
				String acc = *acc_it;
				accessions << acc.substitute('/', '_');
			}
			String protein = accessions.concatenate("/");
			if (filter_charge)
			{
				// write individual abundances (one line for each charge state):
				for (map<Int, SampleAbundances>::const_iterator ab_it = 
							 q_it->second.abundances.begin(); ab_it != 
							 q_it->second.abundances.end(); ++ab_it)
				{
					out << q_it->first.toString() << protein << accessions.size() 
							<< ab_it->first;
					for (ConsensusMap::FileDescriptions::iterator file_it = 
								 files_.begin(); file_it != files_.end(); ++file_it)
					{
						// write abundance for the sample if it exists, 0 otherwise:
						SampleAbundances::const_iterator pos = 
							ab_it->second.find(file_it->first);
						out << (pos != ab_it->second.end() ? pos->second : 0.0);
					}
					out << endl;
				}
			}
			else
			{
				// write total abundances (accumulated over all charge states):
				out << q_it->first.toString() << protein << accessions.size() << 0;
				for (ConsensusMap::FileDescriptions::iterator file_it = 
								 files_.begin(); file_it != files_.end(); ++file_it)
				{
					// write abundance for the sample if it exists, 0 otherwise:
					SampleAbundances::const_iterator pos = 
						q_it->second.total_abundances.find(file_it->first);
					out << (pos != q_it->second.total_abundances.end() ? 
									pos->second : 0.0);
				}
				out << endl;					
			}
		}
	}


	/// Write a table of protein results.
	void writeProteinTable_(SVOutStream& out, const ProteinQuant& quant)
	{
		// write header:
		out << "protein" << "n_proteins" << "protein_score" << "n_peptides";
		if (files_.size() <= 1)
		{
			out << "abundance";
		}
		else
		{
			for (Size i = 1; i <= files_.size(); ++i) 
			{
				out << "abundance_" + String(i);
			}
		}
		out << endl;

		map<String, StringList> leader_to_accessions;
		if (!proteins_.getIndistinguishableProteins().empty())
		{
			for (vector<ProteinIdentification::ProteinGroup>::iterator group_it = 
						 proteins_.getIndistinguishableProteins().begin(); group_it !=
						 proteins_.getIndistinguishableProteins().end(); ++group_it)
			{
				StringList& accessions = leader_to_accessions[group_it->
																											accessions[0]];
				accessions = group_it->accessions;
				for (StringList::Iterator acc_it = accessions.begin(); 
						 acc_it != accessions.end(); ++acc_it)
				{
					acc_it->substitute('/', '_'); // to allow concatenation later
				}
			}
		}

		for (ProteinQuant::const_iterator q_it = quant.begin(); q_it != quant.end();
				 ++q_it)
		{
			if (q_it->second.total_abundances.empty()) continue; // not quantified

			if (leader_to_accessions.empty())
			{
				out << q_it->first << 1;
			}
			else
			{
				out << leader_to_accessions[q_it->first].concatenate('/')
						<< leader_to_accessions[q_it->first].size();
			}
			if (proteins_.getHits().empty())
			{
				out << 0;
			}
			else
			{
				vector<ProteinHit>::iterator pos = proteins_.findHit(q_it->first);
				out << pos->getScore();
			}
			Size n_peptide = q_it->second.abundances.size();
			out << n_peptide;
			// make a copy to allow using "operator[]" below:
			SampleAbundances total_abundances = q_it->second.total_abundances;
			for (ConsensusMap::FileDescriptions::iterator file_it = files_.begin(); 
					 file_it != files_.end(); ++file_it)
			{
				out << total_abundances[file_it->first];
			}
			out << endl;
		}
	}

		
	/// Write comment lines before a peptide/protein table.
	void writeComments_(SVOutStream& out, const bool proteins=true)
	{
		String what = (proteins ? "Protein" : "Peptide"); 
		bool old = out.modifyStrings(false);
		out << "# " + what + " abundances computed from file '" + 
			getStringOption_("in") + "'" << endl;
		StringList relevant_params;
		if (proteins) // parameters relevant only for protein output
		{
			relevant_params << "top" << "average" << "include_all";
		}
		relevant_params << "filter_charge"; // also relevant for peptide output
		if (files_.size() > 1) // flags only for consensusXML input
		{
			relevant_params << "consensus:normalize";
			if (proteins) relevant_params << "consensus:fix_peptides";
		}
		String params;
		for (StringList::Iterator it = relevant_params.begin(); 
				 it != relevant_params.end(); ++it)
		{
			String value = algo_params_.getValue(*it);
			if (value != "false") params += *it + "=" + value + ", ";
		}
		if (params.empty()) params = "(none)";
		else params.resize(params.size() - 2); // remove trailing ", "
		out << "# Parameters (relevant only): " + params << endl;

		if (files_.size() > 1)
		{
			String desc = "# Files/samples associated with abundance values below: ";
			Size counter = 1;
			for (ConsensusMap::FileDescriptions::iterator it = files_.begin(); 
					 it != files_.end(); ++it, ++counter)
			{
				if (counter > 1) desc += ", ";
				desc += String(counter) + ": '" + it->second.filename + "'";
				String label = it->second.label;
				if (!label.empty()) desc += " ('" + label + "')";
			}
			out << desc << endl;
		}
		out.modifyStrings(old);
	}


	/// Write processing statistics.
	void writeStatistics_(const Statistics& stats)
	{
		LOG_INFO << "\nProcessing summary - number of..." 
						 << "\n...features: " << stats.quant_features
						 << " used for quantification, " << stats.total_features 
						 << " total (" << stats.blank_features << " no annotation, "
						 << stats.ambig_features << " ambiguous annotation)"
						 << "\n...peptides: "  << stats.quant_peptides 
						 << " quantified, " << stats.total_peptides 
						 << " identified (considering best hits only)";
		if (!getStringOption_("out").empty() || !getStringOption_("id_out").empty())
		{
			bool include_all = algo_params_.getValue("include_all") == "true";
			Size top = algo_params_.getValue("top");
			LOG_INFO << "\n...proteins/protein groups: " << stats.quant_proteins
							 << " quantified";
			if (top > 1)
			{
				if (include_all) LOG_INFO << " (incl. ";
				else LOG_INFO << ", ";
				LOG_INFO << stats.too_few_peptides << " with fewer than " << top 
								 << " peptides";
				if (stats.n_samples > 1) LOG_INFO << " in every sample";
				if (include_all) LOG_INFO << ")";
			}
		}
		LOG_INFO << endl;
	}

	
	/// Annotate a Protein-/PeptideHit with abundance values (for mzTab export).
	template <typename HitType>
	void storeAbundances_(HitType& hit, SampleAbundances& total_abundances, 
												const String& what = "protein")
	{
		Size counter = 1;
		for (ConsensusMap::FileDescriptions::iterator file_it = files_.begin(); 
				 file_it != files_.end(); ++file_it, ++counter)
		{
			String field[2] = {"mzTab:" + what + "_abundance_", 
												 "sub[" + String(counter) + "]"};
			DoubleReal value = total_abundances[file_it->first];
			if (value > 0) hit.setMetaValue(field[0] + field[1], value);
			else hit.setMetaValue(field[0] + field[1], "--"); // missing value
			// TODO: compute std. deviations/std. errors (how?)
			// hit.setMetaValue(field[0] + "stdev_" + field[1], "--");
			// hit.setMetaValue(field[0] + "std_error_" + field[1], "--");
		}
	}


	void prepareMzTab_(const ProteinQuant& prot_quant, 
										 const PeptideQuant& pep_quant, 
										 vector<DataProcessing>& processing)
	{
		// proteins:
		// mapping: protein accession -> position in list of protein hits
		typedef map<String, vector<ProteinHit>::iterator> AccessionMap;
		AccessionMap accession_map;
		for (vector<ProteinHit>::iterator ph_it = proteins_.getHits().begin();
				 ph_it != proteins_.getHits().end(); ++ph_it)
		{
			accession_map[ph_it->getAccession()] = ph_it;
			}
			// indistinguishable proteins:
			map<String, StringList> leader_to_accessions;
			for (vector<ProteinIdentification::ProteinGroup>::iterator group_it = 
						 proteins_.getIndistinguishableProteins().begin(); group_it !=
						 proteins_.getIndistinguishableProteins().end(); ++group_it)
			{
				if (group_it->accessions.size() > 1)
				{
					StringList& acc = leader_to_accessions[group_it->accessions[0]];
					acc.insert(acc.end(), ++group_it->accessions.begin(), 
										 group_it->accessions.end()); // insert all but first
				}
			}
		// annotate protein hits with abundances:
		vector<ProteinHit> quantified_prot;
		for (ProteinQuant::const_iterator q_it = prot_quant.begin(); 
				 q_it != prot_quant.end(); ++q_it)
			{
			if (accession_map.empty()) // generate a new hit
				{
				ProteinHit hit;
				hit.setAccession(q_it->first); // no further data
				quantified_prot.push_back(hit);
				}
			else // copy existing hit
			{
				AccessionMap::iterator pos = accession_map.find(q_it->first);
				if (pos == accession_map.end()) continue; // not in list, skip
				quantified_prot.push_back(*(pos->second));
				// annotate with indistinguishable proteins:
				map<String, StringList>::iterator la_it = 
					leader_to_accessions.find(q_it->first);
				if (la_it != leader_to_accessions.end()) // protein is group member
				{
					quantified_prot.back().setMetaValue("mzTab:ambiguity_members", 
																				 la_it->second.concatenate(","));
				}
			}
			// annotate with abundance values:
			SampleAbundances total_abundances = q_it->second.total_abundances;
			storeAbundances_(quantified_prot.back(), total_abundances);
			quantified_prot.back().setMetaValue("num_peptides", 
																					q_it->second.id_count);
		}
		proteins_.getHits().swap(quantified_prot);
		// set meta values:
		UInt64 id = UniqueIdGenerator::getUniqueId();
		proteins_.setMetaValue("mzTab:unit_id", "OpenMS_" + String(id));
		proteins_.setMetaValue("mzTab:title", 
													 "Quantification by OpenMS/ProteinQuantifier");
		processing.push_back(getProcessingInfo_(DataProcessing::QUANTITATION));
		for (Size i = 0; i < processing.size(); ++i)
		{
			Software sw = processing[i].getSoftware();
			String param = "[" + sw.getName() + "," + sw.getVersion() + "]";
			proteins_.setMetaValue("mzTab:software[" + String(i+1) + "]", param);
		}
		// mzTab:ms_file[1]-format and mzTab:ms_file[1]-location: set here to
		// input featureXML/consensusXML, or in MzTabExporter to input idXML?
		for (Size i = 0; i < files_.size(); ++i)
		{
			String loc = "mzTab:ms_file[" + String(i+1) + "]-location";
			proteins_.setMetaValue(loc, files_[i].filename);
			if (!files_[i].label.empty())
			{
				String desc =  "mzTab:sub[" + String(i+1) + "]-description";
				proteins_.setMetaValue(desc, "label: " + files_[i].label);
			}
		}
		// peptides:
		// mapping: unmodified peptide seq. -> position in list of peptide hits
		typedef map<String, vector<PeptideHit>::const_iterator> SequenceMap;
		SequenceMap sequence_map;
		for (vector<PeptideHit>::const_iterator ph_it = peptides_.getHits().begin();
				 ph_it != peptides_.getHits().end(); ++ph_it)
		{
			// ProteinProphet results list unmodified sequences anyway...
			sequence_map[ph_it->getSequence().toUnmodifiedString()] = ph_it;
		}
		map<String, String> pep2prot; // unmod. peptides used for protein quant.
		for (ProteinQuant::const_iterator q_it = prot_quant.begin();
				 q_it != prot_quant.end(); ++q_it)
		{
			for (map<String, SampleAbundances>::const_iterator ab_it = 
						 q_it->second.abundances.begin(); ab_it != 
						 q_it->second.abundances.end(); ++ ab_it)
			{		
				pep2prot[ab_it->first] = q_it->first;
			}
		}
		// annotate peptide hits with abundances:
		vector<PeptideHit> quantified_pep;
			for (PeptideQuant::const_iterator q_it = pep_quant.begin(); 
					 q_it != pep_quant.end(); ++q_it)
			{
			if (sequence_map.empty()) // generate a new hit
			{
				PeptideHit hit;
				quantified_pep.push_back(hit);
			}
			else // copy existing hit
			{
				SequenceMap::iterator pos = 
					sequence_map.find(q_it->first.toUnmodifiedString());
				if (pos == sequence_map.end()) continue; // not in list, skip
				quantified_pep.push_back(*(pos->second));
		}
			quantified_pep.back().setSequence(q_it->first);
			// set protein accession only for proteotypic peptides:
			map<String, String>::iterator pos = 
				pep2prot.find(q_it->first.toUnmodifiedString());
			if (pos == pep2prot.end()) // not proteotypic
		{
				quantified_pep.back().setProteinAccessions(vector<String>());
				quantified_pep.back().setMetaValue("mzTab:unique", "false");
			}
			else // proteotypic (therefore used for quantification)
			{
				vector<String> accessions(1, pos->second);
				quantified_pep.back().setProteinAccessions(accessions);
				quantified_pep.back().setMetaValue("mzTab:unique", "true");
			}
			if (algo_params_.getValue("filter_charge") != "true") // all charges
			{
				SampleAbundances total_abundances = q_it->second.total_abundances;
				storeAbundances_(quantified_pep.back(), total_abundances, "peptide");
			}
			else // generate hits for individual charge states
			{
				for (map<Int, SampleAbundances>::const_iterator ab_it = 
							 q_it->second.abundances.begin(); ab_it != 
							 q_it->second.abundances.end(); ++ab_it)
				{
					if (ab_it != q_it->second.abundances.begin())
					{
						quantified_pep.push_back(quantified_pep.back()); // copy last entry
					}
					quantified_pep.back().setCharge(ab_it->first);
					SampleAbundances charge_abundances = ab_it->second;
					storeAbundances_(quantified_pep.back(), charge_abundances, "peptide");
				}
			}
		}
		peptides_.setHits(quantified_pep);
			
		// remove possibly outdated meta data:
		proteins_.getProteinGroups().clear();
		proteins_.getIndistinguishableProteins().clear();
		// make sure identifiers match:
		if (proteins_.getIdentifier().empty())
		{
			proteins_.setIdentifier(String(UniqueIdGenerator::getUniqueId()));
	}
		peptides_.setIdentifier(proteins_.getIdentifier());
	}


	ExitCodes main_(int, const char**)
  {
		String in = getStringOption_("in");
		String out = getStringOption_("out");
		String peptide_out = getStringOption_("peptide_out");
		String id_out = getStringOption_("id_out");

		if (out.empty() && peptide_out.empty() && id_out.empty())
		{
			throw Exception::RequiredParameterNotGiven(__FILE__, __LINE__,
																								 __PRETTY_FUNCTION__, 
																								 "out/peptide_out/id_out");
		}

		String protxml = getStringOption_("protxml");

		PeptideAndProteinQuant quantifier;
		// algo_params_ = getParam_().copy("algorithm:", true);
		algo_params_ = quantifier.getParameters();
		Logger::LogStream nirvana; // avoid parameter update messages
		algo_params_.update(getParam_(), false, false, nirvana);
		// algo_params_.update(getParam_());
		quantifier.setParameters(algo_params_);

		vector<DataProcessing> processing;
		FileTypes::Type in_type = FileHandler::getType(in);

		if (in_type == FileTypes::FEATUREXML)
		{
			FeatureMap<> features;
			FeatureXMLFile().load(in, features);
			if (!id_out.empty()) 
			{
				processing = features.getDataProcessing();
			}
			files_[0].filename = in;
			// ProteinProphet results in the featureXML?
			if (protxml.empty() && 
					(features.getProteinIdentifications().size() == 1) &&
					(!features.getProteinIdentifications()[0].getHits().empty()))
			{
				proteins_ = features.getProteinIdentifications()[0];
			}
			quantifier.quantifyPeptides(features);
		}
		else // consensusXML
		{
			ConsensusMap consensus;
			ConsensusXMLFile().load(in, consensus);
			files_ = consensus.getFileDescriptions(); 
			if (!id_out.empty()) processing = consensus.getDataProcessing();
			// ProteinProphet results in the consensusXML?
			if (protxml.empty() && 
					(consensus.getProteinIdentifications().size() == 1) &&
					(!consensus.getProteinIdentifications()[0].getHits().empty()))
			{
				proteins_ = consensus.getProteinIdentifications()[0];
			}
			quantifier.quantifyPeptides(consensus);
		}

		if (!out.empty() || !id_out.empty()) // quantify on protein level
		{
			if (!protxml.empty()) // read ProteinProphet data
			{
				vector<ProteinIdentification> proteins;
				vector<PeptideIdentification> peptides;
				IdXMLFile().load(protxml, proteins, peptides);
				if ((proteins.size() == 1) && (peptides.size() == 1))
				{
					proteins_ = proteins[0];
					peptides_ = peptides[0];
				}
				else
				{
					throw Exception::InvalidParameter(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Expected a converted protXML file (with only one 'ProteinIdentification' and one 'PeptideIdentification' instance) in file '" + protxml + "'");
				}
			}
			quantifier.quantifyProteins(proteins_);
		}

		// output:
		String separator = getStringOption_("format:separator");
		String replacement = getStringOption_("format:replacement");
		String quoting = getStringOption_("format:quoting");
		if (separator == "") separator = "\t";
		String::QuotingMethod quoting_method;
		if (quoting == "none") quoting_method = String::NONE;
		else if (quoting == "double") quoting_method = String::DOUBLE;
		else quoting_method = String::ESCAPE;

		if (!peptide_out.empty())
		{
			ofstream outstr(peptide_out.c_str());
			SVOutStream output(outstr, separator, replacement, quoting_method);
			writeComments_(output, false);
			writePeptideTable_(output, quantifier.getPeptideResults());
			outstr.close();
		}
		if (!out.empty())
		{
			ofstream outstr(out.c_str());
			SVOutStream output(outstr, separator, replacement, quoting_method);
			writeComments_(output);
			writeProteinTable_(output, quantifier.getProteinResults());
			outstr.close();
		}
		if (!id_out.empty())
		{
			prepareMzTab_(quantifier.getProteinResults(), 
										quantifier.getPeptideResults(), processing);
			vector<ProteinIdentification> proteins(1, proteins_);
			// create one peptide identification for each peptide hit:
			PeptideIdentification temp = peptides_;
			temp.setHits(vector<PeptideHit>());
			vector<PeptideIdentification> peptides(peptides_.getHits().size(), temp);
			for (Size i = 0; i < peptides.size(); ++i)
			{
				peptides[i].insertHit(peptides_.getHits()[i]);
			}
			IdXMLFile().store(id_out, proteins, peptides);
		}
		
		writeStatistics_(quantifier.getStatistics());
		
		return EXECUTION_OK;				
	}
};
	

int main(int argc, const char** argv)
{
  TOPPProteinQuantifier t;
  return t.main(argc, argv);
}

/// @endcond
