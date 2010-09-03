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

#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/SVOutStream.h>

#include <numeric> // for "accumulate"
#include <algorithm> // for "equal"

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
			<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_FeatureLinker </td>
		</tr>
	</table>
</CENTER>

	Quantification is based on the intensity values of the features in the input. Feature intensities are first accumulated to peptide abundances, according to the peptide identifications annotated to the features/feature groups. Then, abundances of the peptides of a protein are averaged to compute the protein abundance.

	The peptide-to-protein step implements a general version of the "top 3 approach" (but only for relative quantification) described in:\n
	Silva <em>et al.</em>: Absolute quantification of proteins by LCMS<sup>E</sup>: a virtue of parallel MS acquisition (Mol. Cell. Proteomics, 2006).

	Only features/feature groups with unambiguous peptide annotation are used for peptide quantification, and generally only proteotypic peptides (i.e. those matching to exactly one protein) are used for protein quantification. As an exception to this rule, if ProteinProphet results for the whole sample set are provided with the @p protxml option, or are already included in a featureXML input, also groups of indistinguishable proteins will be quantified. The reported quantity then refers to the total for the whole group.

	Peptide/protein IDs from multiple identification runs can be handled, but will not be differentiated (i.e. protein accessions for a peptide will be accumulated over all identification runs).

	More information below the parameter specification.

	<B>The command line parameters of this tool are:</B>
	@verbinclude TOPP_ProteinQuantifier.cli


	The output files produced by this tool have a table format, with columns as described below:

	<b>Protein output</b> (one protein/set of indistinguishable proteins per line):
	- @b protein: Protein accession(s) (as in the annotations in the input file; separated by "/" if more than one).
	- @b n_proteins: Number of indistinguishable proteins quantified (usually "1").
	- @b protein_score: Protein score, e.g. ProteinProphet probability (if available).
	- @b n_peptides: Number of proteotypic peptides observed for this protein (or group of indistinguishable proteins) across all samples. Note that not necessarily all of these peptides contribute to the protein abundance (depending on parameter @p top).
	- @b abundance: Computed protein abundance. For consensusXML input, there will be one column  per sample ("abundance_0", "abundance_1", etc.).

	<b>Peptide output</b> (one peptide or - if @p filter_charge is set - one charge state of a peptide per line):
	- @b peptide: Peptide sequence. Only peptides that occur in unambiguous annotations of features are reported.
	- @b protein: Protein accession(s) for the peptide (separated by "/" if more than one).
	- @b n_proteins: Number of proteins this peptide maps to. (Same as the number of accessions in the previous column.)
	- @b charge: Charge state quantified in this line. "0" (for "all charges") unless @p filter_charge was set.
	- @b abundance: Computed abundance for this peptide. If the charge in the preceding column is 0, this is the total abundance of the peptide over all charge states; otherwise, it is only the abundance observed for the indicated charge (in this case, there may be more than one line for the peptide sequence). Again, for consensusXML input, there will be one column  per sample ("abundance_0", "abundance_1", etc.). Also for consensusXML, the reported values are already normalized if @p consensus:normalize was set.


	In addition to the information above, consider the following for parameter selection: With @p filter_charge and @p average, there is a trade-off between comparability of protein abundances within a sample and of abundances for the same protein across different samples.\n
	Setting @p filter_charge may increase reproducibility between samples, but will distort the proportions of protein abundances within a sample. The reason is that ionization properties vary between peptides, but should remain constant across samples. Filtering by charge state can help to reduce the impact of feature detection differences between samples.\n
	For @p average, there is a qualitative difference between @p mean/median and @p sum in the effect that missing peptide abundances have (only if @p include_all is set): @p mean and @p median ignore missing cases, averaging only present values. If low-abundant peptides are not detected in some samples, the computed protein abundances for those samples may thus be too optimistic. @p sum implicitly treats missing values as zero, so this problem does not occur and comparability across samples is ensured. However, with @p sum the total number of peptides ("summands") available for a protein may affect the abundances computed for it (depending on @p top), so results within a sample may become unproportional.

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

namespace OpenMS
{

  class TOPPProteinQuantifier 
  	: public TOPPBase
  {
	public:
		
		TOPPProteinQuantifier() :
			TOPPBase("ProteinQuantifier", "Compute peptide and protein abundances"),
			proteins_()
      {
      }

	protected:

		typedef map<UInt64, DoubleReal> sample_abundances; // abundance by sample

		/// Quantitative and associated data for a peptide
		struct peptide_data
		{
			map<Int, sample_abundances> abundances; // charge -> sample -> abundance
			sample_abundances total_abundances; // sample -> total abundance
			set<String> accessions; // protein accessions
		};
		typedef map<AASequence, peptide_data> peptide_quant; // by peptide sequence

		/// Quantitative and associated data for a protein
		struct protein_data
		{
			// peptide -> sample -> abundance:
			map<AASequence, sample_abundances> abundances;
			sample_abundances total_abundances; // sample -> total abundance
		};
		typedef map<String, protein_data> protein_quant; // by protein accession

		ProteinIdentification proteins_; // protein information from protXML

		/// Statistics for processing summary
		struct statistics 
		{
			Size quant_proteins, too_few_peptides; // protein statistics
			Size quant_peptides, total_peptides; // peptide statistics
			// feature statistics:
			Size quant_features, total_features, blank_features, ambig_features;
			statistics(): quant_proteins(0), too_few_peptides(0), quant_peptides(0),
										total_peptides(0), quant_features(0), total_features(0), 
										blank_features(0), ambig_features(0) {}
		} stats_; // for output in the end


	  void registerOptionsAndFlags_()
    {
		  registerInputFile_("in", "<file>", "", "Input file");
		  setValidFormats_("in", StringList::create("featureXML,consensusXML"));
		  registerInputFile_("protxml", "<file>", "", "ProteinProphet results (protXML converted to idXML) for the identification runs that were used to annotate the input.\nInformation about indistinguishable proteins will be used for protein quantification.", false);
		  setValidFormats_("protxml", StringList::create("idXML"));
      registerOutputFile_("out", "<file>", "", "Output file for protein abundances", false);
		  registerOutputFile_("peptide_out", "<file>", "", "Output file for peptide abundances\nEither 'out' or 'peptide_out' are required. They can be used together.", false);

		  addEmptyLine_();
		  registerIntOption_("top", "<number>", 3, "Calculate protein abundance from this number of proteotypic peptides (best first; '0' for all)", false);
		  setMinInt_("top", 0);
		  registerStringOption_("average", "<method>", "median", "Averaging method used to compute protein abundances from peptide abundances", false);
		  setValidStrings_("average", StringList::create("median,mean,sum"));
		  registerFlag_("include_all", "Include results for proteins with fewer than 'top' proteotypic peptides");
		  registerFlag_("filter_charge", "Distinguish between charge states of a peptide. For peptides, abundances will be reported separately for each charge;\nfor proteins, abundances will be computed based only on the most prevalent charge of each peptide.\nBy default, abundances are summed over all charge states.");

		  registerTOPPSubsection_("consensus", "Additional options for consensusXML input");
		  registerFlag_("consensus:normalize", "Scale peptide abundances so that medians of all samples are equal");
		  registerFlag_("consensus:fix_peptides", "Use the same peptides for protein quantification across all samples.\nThe 'top' peptides that occur each in the highest number of samples are selected (breaking ties by total abundance),\nbut there is no guarantee that these will be the best co-ocurring peptides.");

		  registerTOPPSubsection_("format", "Output formatting options");
		  registerStringOption_("format:separator", "<sep>", "", "Character(s) used to separate fields; by default, the 'tab' character is used", false);
		  registerStringOption_("format:quoting", "<method>", "double", "Method for quoting of strings: 'none' for no quoting, 'double' for quoting with doubling of embedded quotes,\n'escape' for quoting with backslash-escaping of embedded quotes", false);
		  setValidStrings_("format:quoting", StringList::create("none,double,escape"));
		  registerStringOption_("format:replacement", "<x>", "_", "If 'quoting' is 'none', used to replace occurrences of the separator in strings before writing", false);
    }

		
		/**
			 @brief Compute the median of a list of values (possibly already sorted)

			 Note that the list @p values must not be empty!
		*/
		DoubleReal median_(DoubleList values, bool sorted=FALSE)
			{
				if (!sorted) sort(values.begin(), values.end());
				Size size = values.size();
				if (size % 2 == 0)  // even size => average two middle values
				{
					size /= 2;
					return (values[size - 1] + values[size]) / 2.0;
				}
				else return values[(size - 1)/2];
			}


		/**
			 @brief Get the "canonical" annotation (a single peptide hit) of a feature/consensus feature from the associated list of peptide identifications.

			 Only the best-scoring peptide hit of each ID in @p peptides is taken into account. If there's more than one ID and the best hits are not identical by sequence, or if there's no peptide ID, an empty peptide hit (for "ambiguous/no annotation") is returned.
			 Protein accessions from identical peptide hits are accumulated.
		*/
		PeptideHit getAnnotation_(vector<PeptideIdentification>& peptides)
			{
				if (peptides.empty()) return PeptideHit();
				peptides.front().sort();
				PeptideHit hit = peptides.front().getHits().front();
				for (vector<PeptideIdentification>::iterator pep_it = 
							 ++peptides.begin(); pep_it != peptides.end(); ++pep_it)
				{
					pep_it->sort();
					const PeptideHit& current = pep_it->getHits().front();
					if (hit.getSequence() != current.getSequence())
					{
						return PeptideHit();
					}
					else
					{
						// add protein accessions:
						for (vector<String>::const_iterator acc_it = 
									 current.getProteinAccessions().begin(); acc_it !=
									 current.getProteinAccessions().end(); ++acc_it)
						{
							hit.addProteinAccession(*acc_it);
						}
					}
				}
				return hit;
			}


		/**
			 @brief Gather quantitative information from a feature.

			 Store quantitative information from @p feature in @ quant, based on the peptide annotation in @p hit. If @p hit is empty ("ambiguous/no annotation"), nothing is stored.
		*/
		void quantifyFeature_(const FeatureHandle& feature, const PeptideHit& hit, 
													peptide_quant& quant)
			{
				if (hit == PeptideHit()) 
				{
					return; // annotation for the feature ambiguous or missing
				}
				stats_.quant_features++;
				const AASequence& seq = hit.getSequence();
				quant[seq].abundances[feature.getCharge()][feature.getMapIndex()] += 
					feature.getIntensity(); // new map element is initialized with 0
				// maybe TODO: check charge of peptide hit against charge of feature
				quant[seq].accessions.insert(hit.getProteinAccessions().begin(),
																		 hit.getProteinAccessions().end());
			}
		

		/**
			 @brief Order keys (charges/peptides for peptide/protein quantification) according to how many samples they allow to quantify, breaking ties by total abundance.

			 The keys of @p abundances are stored ordered in @p result, best first.
		*/
		template <typename T>
		void orderBest_(const map<T, sample_abundances> abundances, 
										vector<T>& result)
			{
				typedef pair<Size, DoubleReal> pair_type;
				multimap<pair_type, T, greater<pair_type> > order;
				for (typename map<T, sample_abundances>::const_iterator ab_it = 
							 abundances.begin(); ab_it != abundances.end(); ++ab_it)
				{
					DoubleReal total = 0.0;
					for (sample_abundances::const_iterator samp_it = 
								 ab_it->second.begin(); samp_it != 
								 ab_it->second.end(); ++samp_it)
					{
						total += samp_it->second;
					}
					pair_type key = make_pair(ab_it->second.size(), total);
					order.insert(make_pair(key, ab_it->first));
				}
				result.clear();
				for (typename multimap<pair_type, T, greater<pair_type> >::iterator 
							 ord_it = order.begin(); ord_it != order.end(); ++ord_it)
			  {
					result.push_back(ord_it->second);
				}
			}


		/**
			 @brief Compute overall peptide quantities.

			 Based on quantitative data for individual charge states (derived from annotated features) in @p quant, compute overall abundances for all peptides and store them also in @p quant.
		*/
		void quantifyPeptides_(peptide_quant& quant)
			{
				for (peptide_quant::iterator q_it = quant.begin(); q_it != quant.end();
						 ++q_it)
				{
					if (getFlag_("filter_charge"))
					{
						// find charge state with abundances for highest number of samples
						// (break ties by total abundance):
						IntList charges; // sorted charge states (best first)
						orderBest_(q_it->second.abundances, charges);
						Int best_charge = charges.front();

						// quantify according to the best charge state only:
						for (sample_abundances::iterator samp_it = 
									 q_it->second.abundances[best_charge].begin(); samp_it != 
									 q_it->second.abundances[best_charge].end(); ++samp_it)
						{
							q_it->second.total_abundances[samp_it->first] = samp_it->second;
						}
					}
					else
					{
						// sum up abundances over all charge states:
						for (map<Int, sample_abundances>::iterator ab_it = 
									 q_it->second.abundances.begin(); ab_it != 
									 q_it->second.abundances.end(); ++ab_it)
						{
							for (sample_abundances::iterator samp_it = ab_it->second.begin();
									 samp_it != ab_it->second.end(); ++samp_it)
							{
								q_it->second.total_abundances[samp_it->first] += 
									samp_it->second;
							}
						}
					}
				}
			}
					
		
		/**
			 @brief Normalize peptide abundances across samples by (multiplicative) scaling to equal medians.
		*/
		void normalizePeptides_(peptide_quant& quant)
			{
				// gather data:
				map<UInt64, DoubleList> abundances; // all peptide abundances by sample
				for (peptide_quant::iterator q_it = quant.begin();
						 q_it != quant.end(); ++q_it)
				{
					// maybe TODO: treat missing abundance values as zero
					for (sample_abundances::iterator samp_it = 
								 q_it->second.total_abundances.begin(); samp_it !=
								 q_it->second.total_abundances.end(); ++samp_it)
					{
						abundances[samp_it->first] << samp_it->second;
					}
				}
				if (abundances.size() <= 1) return;
				
				// compute scale factors for all samples:
				sample_abundances medians; // median abundances by sample
				for (map<UInt64, DoubleList>::iterator ab_it = abundances.begin();
						 ab_it != abundances.end(); ++ab_it)
				{
					medians[ab_it->first] = median_(ab_it->second);
				}
				DoubleList all_medians;
				for (sample_abundances::iterator med_it = medians.begin();
						 med_it != medians.end(); ++med_it)
				{
					all_medians << med_it->second;
				}
				DoubleReal overall_median = median_(all_medians);
				sample_abundances scale_factors;
				for (sample_abundances::iterator med_it = medians.begin();
						 med_it != medians.end(); ++med_it)
				{
					scale_factors[med_it->first] = overall_median / med_it->second;
				}

				// scale all abundance values:
				for (peptide_quant::iterator q_it = quant.begin();
						 q_it != quant.end(); ++q_it)
				{
					for (sample_abundances::iterator tot_it = 
								 q_it->second.total_abundances.begin(); tot_it !=
								 q_it->second.total_abundances.end(); ++tot_it)
					{
						tot_it->second *= scale_factors[tot_it->first];
					}
					for (map<Int, sample_abundances>::iterator ab_it =
								 q_it->second.abundances.begin(); ab_it !=
								 q_it->second.abundances.end(); ++ab_it)
					{
						for (sample_abundances::iterator samp_it = ab_it->second.begin();
								 samp_it != ab_it->second.end(); ++samp_it)
						{
							samp_it->second *= scale_factors[samp_it->first];
						}
					}
				}				
			}


		/**
			 @brief Get the "canonical" protein accession from the list of protein accessions of a peptide.

			 @param pep_accessions Protein accessions of a peptide
			 @param accession_to_leader Captures information about indistinguishable proteins (maps accession to accession of group leader)

			 If there is no information about indistinguishable proteins (from protXML) available, a canonical accession exists only for proteotypic peptides - it's the single accession for this peptide.

			 If the information is available, a peptide has a canonical accession if it maps only to proteins of one indistinguishable group. In this case, the canonical accession is that of the group leader.

			 If there is no canonical accession, the empty string is returned.
		*/
		String getAccession_(const set<String>& pep_accessions, 
												 map<String, String>& accession_to_leader)
			{
				if (accession_to_leader.empty())
				{ 
					// no info about indistinguishable proteins available
					if (pep_accessions.size() == 1)	return *pep_accessions.begin();
				}
				else
				{ 
					// if all accessions belong to the same group of indistinguishable
					// proteins, return accession of the group leader
					StringList leaders;
					for (set<String>::const_iterator it = pep_accessions.begin(); 
							 it != pep_accessions.end(); ++it)
					{
						map<String, String>::const_iterator pos = 
							accession_to_leader.find(*it);
						if (pos != accession_to_leader.end()) leaders << pos->second;
						// if the protein accession was not found, this is not an error:
						// if there's not enough evidence for a protein, it won't occur in
						// the protXML - so we also won't quantify it
					}
					if (leaders.empty()) return "";
					bool all_equal = equal(leaders.begin(), --leaders.end(), 
																 ++leaders.begin());
					if (all_equal) return leaders[0];
				}
				return "";
			}

		
		/**
			 @brief Compute protein quantities.

			 Based on quantitative data for peptides in @p pep_quant, compute protein abundances and store them in @p prot_quant.
		*/
		void quantifyProteins_(const peptide_quant& pep_quant, 
													 protein_quant& prot_quant)
			{
				// if information about indistinguishable proteins is available, map
				// each accession to the accession of the leader of its group of
				// indistinguishable proteins:
				map<String, String> accession_to_leader;
				if (!proteins_.getIndistinguishableProteins().empty())
				{
					for (vector<ProteinIdentification::ProteinGroup>::iterator group_it = 
								 proteins_.getIndistinguishableProteins().begin(); group_it !=
								 proteins_.getIndistinguishableProteins().end(); ++group_it)
					{
						for (StringList::Iterator acc_it = group_it->accessions.begin();
								 acc_it != group_it->accessions.end(); ++acc_it)
						{
							// each accession should only occur once, but we don't check...
							accession_to_leader[*acc_it] = group_it->accessions[0];
						}
					}
				}

				for (peptide_quant::const_iterator pep_it = pep_quant.begin(); 
						 pep_it != pep_quant.end(); ++pep_it)
				{
					String accession = getAccession_(pep_it->second.accessions, 
																					 accession_to_leader);
					if (!accession.empty()) // proteotypic peptide
					{
						for (sample_abundances::const_iterator tot_it = 
									 pep_it->second.total_abundances.begin(); tot_it !=
									 pep_it->second.total_abundances.end(); ++tot_it)
						{
							prot_quant[accession].abundances[pep_it->first][tot_it->first] =
								tot_it->second;
						}
					}
				}

				Size top = getIntOption_("top");
				String average = getStringOption_("average");
				bool include_all = getFlag_("include_all"), 
					fix_peptides = getFlag_("consensus:fix_peptides");

				for (protein_quant::iterator prot_it = prot_quant.begin();
						 prot_it != prot_quant.end(); ++prot_it)
				{
					if (prot_it->second.abundances.size() < top)
					{
						if (include_all) stats_.too_few_peptides++;
						else continue; // not enough proteotypic peptides
					}

					vector<AASequence> peptides; // peptides selected for quantification
					if (fix_peptides && (prot_it->second.abundances.size() > top))
					{
						// consider only "top" best peptides
						orderBest_(prot_it->second.abundances, peptides);
						peptides.resize(top);
					}
					else // consider all peptides
					{
						for (map<AASequence, sample_abundances>::iterator ab_it =
									 prot_it->second.abundances.begin(); ab_it !=
									 prot_it->second.abundances.end(); ++ab_it)
						{
							peptides.push_back(ab_it->first);							
						}
					}

					map<UInt64, DoubleList> abundances; // all pept. abundances by sample
					// consider only the peptides selected above for quantification:
					for (vector<AASequence>::iterator pep_it = peptides.begin();
							 pep_it != peptides.end(); ++pep_it)
					{
						sample_abundances& current_ab = prot_it->second.abundances[*pep_it];
						for (sample_abundances::iterator samp_it = current_ab.begin();
								 samp_it != current_ab.end(); ++samp_it)
						{
							abundances[samp_it->first] << samp_it->second;
						}
					}

					for (map<UInt64, DoubleList>::iterator ab_it = abundances.begin();
							 ab_it != abundances.end(); ++ab_it)
					{
						if (!include_all && (ab_it->second.size() < top))
						{
							continue; // not enough peptide abundances for this sample
						}
						if ((top > 0) && (ab_it->second.size() > top))
						{
							sort(ab_it->second.begin(), ab_it->second.end(), 
									 greater<double>()); // sort descending
							ab_it->second.resize(top); // remove all but best "top" values
						}

						DoubleReal result;
						if (average == "median")
						{
							result = median_(ab_it->second);
						}
						else // "mean" or "sum"
						{
							result = accumulate(ab_it->second.begin(), ab_it->second.end(), 
																	0.0);
							if (average == "mean") result /= ab_it->second.size();
						}
						prot_it->second.total_abundances[ab_it->first] = result;
					}
				}
			}
		

		/// Write a table of peptide results.
		void writePeptideTable_(SVOutStream& out, const peptide_quant& quant, 
														const vector<UInt64> samples)
			{
				// write header:
				out << "peptide" << "protein" << "n_proteins" << "charge";
				if (samples.size() <= 1)
				{
					out << "abundance";
				}
				else
				{
					for (Size i = 0; i < samples.size(); ++i) {
						out << "abundance_" + String(i);
					}
				}
				out << endl;

				bool filter_charge = getFlag_("filter_charge");
				for (peptide_quant::const_iterator q_it = quant.begin(); 
						 q_it != quant.end(); ++q_it)
				{
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
						for (map<Int, sample_abundances>::const_iterator ab_it = 
									 q_it->second.abundances.begin(); ab_it != 
									 q_it->second.abundances.end(); ++ab_it)
						{
							out << q_it->first.toString() << protein << accessions.size() 
									<< ab_it->first;
							for (vector<UInt64>::const_iterator samp_it = samples.begin(); 
									 samp_it != samples.end(); ++samp_it)
							{
								// write abundance for the sample if it exists, 0 otherwise:
								sample_abundances::const_iterator pos = 
									ab_it->second.find(*samp_it);
								out << (pos != ab_it->second.end() ? pos->second : 0.0);
							}
							out << endl;
						}
					}
					else
					{
						// write total abundances (accumulated over all charge states):
						out << q_it->first.toString() << protein << accessions.size() << 0;
						for (vector<UInt64>::const_iterator samp_it = samples.begin(); 
								 samp_it != samples.end(); ++samp_it)
						{
							// write abundance for the sample if it exists, 0 otherwise:
							sample_abundances::const_iterator pos = 
								q_it->second.total_abundances.find(*samp_it);
							out << (pos != q_it->second.total_abundances.end() ? 
											pos->second : 0.0);
						}
						out << endl;					
					}
				}
			}


		/// Write a table of protein results.
		void writeProteinTable_(SVOutStream& out, protein_quant& quant,
														const vector<UInt64> samples)
			{
				// write header:
				out << "protein" << "n_proteins" << "protein_score" << "n_peptides";
				if (samples.size() <= 1)
				{
					out << "abundance";
				}
				else
				{
					for (Size i = 0; i < samples.size(); ++i) 
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

				for (protein_quant::iterator q_it = quant.begin(); q_it != quant.end();
						 ++q_it)
				{
					if (q_it->second.total_abundances.empty())
					{
						stats_.too_few_peptides++;
						continue;
					}

					stats_.quant_proteins++;
					Size n_peptide = q_it->second.abundances.size();
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
					out << n_peptide;
					for (vector<UInt64>::const_iterator samp_it = samples.begin(); 
							 samp_it != samples.end(); ++samp_it)
					{
						out << q_it->second.total_abundances[*samp_it];
					}
					out << endl;
				}
			}

		
		/// Write comment lines before a peptide/protein table.
		void writeComments_(SVOutStream& out, 
												const ConsensusMap::FileDescriptions files,
												const bool proteins=true)
			{
				String what = (proteins ? "Protein" : "Peptide"); 
				bool old = out.modifyStrings(false);
				out << "# " + what + " abundances computed from file '" + 
					getStringOption_("in") + "'" << endl;
				String params;
				StringList flags;
				if (proteins) // parameters relevant only for protein output
				{
					params = "top=" + String(getIntOption_("top")) + ", average=" + 
						getStringOption_("average") + ", ";
					flags << "include_all";
				}
				flags << "filter_charge"; // also relevant for peptide output
				if (files.size() > 1) // flags only for consensusXML input
				{
					flags << "consensus:normalize";
					if (proteins) flags << "consensus:fix_peptides";
				}
				for (StringList::Iterator it = flags.begin(); it != flags.end(); ++it)
				{
					if (getFlag_(*it)) params += *it + ", ";
				}
				if (params.empty()) params = "(none)";
				else params.resize(params.size() - 2); // remove trailing ", "
				out << "# Parameters (relevant only): " + params << endl;

				if (files.size() > 1)
				{
					String desc = 
						"# Files/samples associated with abundance values below: ";
					Size counter = 0;
					for (ConsensusMap::FileDescriptions::const_iterator it = 
								 files.begin(); it != files.end(); ++it, ++counter)
					{
						if (counter > 0) desc += ", ";
						desc += String(counter) + ": '" + it->second.filename + "'";
						String label = it->second.label;
						if (!label.empty()) desc += " ('" + label + "')";
					}
					out << desc << endl;
				}
				out.modifyStrings(old);
			}


		/// Collect peptide sequences of best hits from a list of peptide identifications.
		void collectPeptideSequences_(vector<PeptideIdentification>& peptides,
																	set<AASequence>& sequences)
			{
				for (vector<PeptideIdentification>::iterator pep_it = 
							 peptides.begin(); pep_it != peptides.end(); ++pep_it)
				{
					if (!pep_it->getHits().empty())
					{
						pep_it->sort();
						sequences.insert(pep_it->getHits()[0].getSequence());
					}
				}
			}


		/// Get the number of distinct identified peptides in a feature map.
		Size countPeptides_(FeatureMap<>& features)
			{
				set<AASequence> sequences;
				for (FeatureMap<>::Iterator feat_it = features.begin(); 
						 feat_it != features.end(); ++feat_it)
				{
					collectPeptideSequences_(feat_it->getPeptideIdentifications(),
																	 sequences);
				}
				collectPeptideSequences_(features.getUnassignedPeptideIdentifications(),
																 sequences);

				return sequences.size();
			}


		/// Get the number of distinct identified peptides in a consensus map.
		Size countPeptides_(ConsensusMap& consensus)
			{
				set<AASequence> sequences;
				for (ConsensusMap::Iterator cons_it = consensus.begin(); 
						 cons_it != consensus.end(); ++cons_it)
				{
					collectPeptideSequences_(cons_it->getPeptideIdentifications(),
																	 sequences);
				}
				collectPeptideSequences_(
					consensus.getUnassignedPeptideIdentifications(), sequences);

				return sequences.size();
			}


    /// Write processing statistics.
		void writeStatistics_(Size n_samples)
			{
				stats_.ambig_features = stats_.total_features - stats_.blank_features - 
					stats_.quant_features;
				LOG_INFO << "\nProcessing summary - number of..." 
								 << "\n...features: " << stats_.quant_features
								 << " used for quantification, " << stats_.total_features 
								 << " total (" << stats_.blank_features << " no annotation, "
								 << stats_.ambig_features << " ambiguous annotation)"
								 << "\n...peptides: "  << stats_.quant_peptides 
								 << " quantified, " << stats_.total_peptides 
								 << " identified (considering best hits only)";
				if (!getStringOption_("out").empty())
				{
					bool include_all = getFlag_("include_all");
					LOG_INFO << "\n...proteins/protein groups: " << stats_.quant_proteins
									 << " quantified";
					if (include_all) LOG_INFO << " (incl. ";
					else LOG_INFO << ", ";
					LOG_INFO << stats_.too_few_peptides << " with fewer than " 
									 << getIntOption_("top") << " peptides";
					if (include_all) LOG_INFO << ")";
					else if (n_samples > 1) LOG_INFO << " in every sample";
				}
				LOG_INFO << endl;
			}

	
		ExitCodes main_(int, const char**)
      {
				String in = getStringOption_("in"), out = getStringOption_("out"), 
					peptide_out = getStringOption_("peptide_out");

				if (out.empty() && peptide_out.empty())
				{
					throw Exception::RequiredParameterNotGiven(__FILE__, __LINE__,
																										 __PRETTY_FUNCTION__, 
																										 "out/peptide_out");				
				}

				FileTypes::Type in_type = FileHandler::getType(in);

				String separator = getStringOption_("format:separator"), 
					replacement = getStringOption_("format:replacement"), 
					quoting = getStringOption_("format:quoting");
        if (separator == "") separator = "\t";
				String::QuotingMethod quoting_method;
				if (quoting == "none") quoting_method = String::NONE;
				else if (quoting == "double") quoting_method = String::DOUBLE;
				else quoting_method = String::ESCAPE;

				String protxml = getStringOption_("protxml");

				bool normalize = false;
				vector<UInt64> samples;
				ConsensusMap::FileDescriptions files;
				peptide_quant pep_quant;

				if (in_type == FileTypes::FEATUREXML)
				{
					FeatureMap<> features;
          FeatureXMLFile().load(in, features);
					samples.push_back(0);

					if (protxml.empty() && 
							(features.getProteinIdentifications().size() == 1) &&
							(!features.getProteinIdentifications()[0].getHits().empty()))
					{
						proteins_ = features.getProteinIdentifications()[0];
					}

					stats_.total_features = features.size();
					for (FeatureMap<>::Iterator feat_it = features.begin(); 
							 feat_it != features.end(); ++feat_it)
					{
						PeptideHit hit = 
							getAnnotation_(feat_it->getPeptideIdentifications());
						FeatureHandle handle(0, *feat_it);
						quantifyFeature_(handle, hit, pep_quant);
						if (feat_it->getPeptideIdentifications().empty())
						{
							stats_.blank_features++;
						}
					}

					stats_.total_peptides = countPeptides_(features);
				}

				else // consensusXML
				{
					normalize = getFlag_("consensus:normalize");
					ConsensusMap consensus;
          ConsensusXMLFile().load(in, consensus);

					files = consensus.getFileDescriptions(); 
					for (ConsensusMap::FileDescriptions::Iterator file_it = files.begin();
							 file_it != files.end(); ++file_it)
					{
						samples.push_back(file_it->first);
					}

					if (protxml.empty() && 
							(consensus.getProteinIdentifications().size() == 1) &&
							(!consensus.getProteinIdentifications()[0].getHits().empty()))
					{
						proteins_ = consensus.getProteinIdentifications()[0];
					}

					for (ConsensusMap::Iterator cons_it = consensus.begin(); 
							 cons_it != consensus.end(); ++cons_it)
					{
						PeptideHit hit = 
							getAnnotation_(cons_it->getPeptideIdentifications());
            for (ConsensusFeature::HandleSetType::const_iterator feat_it = 
									 cons_it->getFeatures().begin(); feat_it !=
									 cons_it->getFeatures().end(); ++feat_it)
						{
							stats_.total_features++;
							quantifyFeature_(*feat_it, hit, pep_quant);
							if (cons_it->getPeptideIdentifications().empty())
							{
								stats_.blank_features++;
							}
						}
					}

					stats_.total_peptides = countPeptides_(consensus);
				}

				if (!protxml.empty())
				{
					vector<ProteinIdentification> proteins;
					vector<PeptideIdentification> peptides;
					IdXMLFile().load(protxml, proteins, peptides);
					if (proteins.size() == 1) 
					{
						proteins_ = proteins[0];
					}
					else
					{
						throw Exception::InvalidParameter(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Expected a converted protXML file (with only one 'ProteinIdentification' instance) in file '" + protxml + "'");
					}
				}

				quantifyPeptides_(pep_quant);
				if (normalize) normalizePeptides_(pep_quant);
				stats_.quant_peptides = pep_quant.size();

				// output:
				if (!peptide_out.empty())
				{
					ofstream outstr(peptide_out.c_str());
					SVOutStream output(outstr, separator, replacement, quoting_method);
					writeComments_(output, files, false);
					writePeptideTable_(output, pep_quant, samples);
					outstr.close();
				}
				if (!out.empty())
				{
					protein_quant prot_quant;
					quantifyProteins_(pep_quant, prot_quant);
					ofstream outstr(out.c_str());
					SVOutStream output(outstr, separator, replacement, quoting_method);
					writeComments_(output, files);
					writeProteinTable_(output, prot_quant, samples);
					outstr.close();
				}

				writeStatistics_(samples.size());

				return EXECUTION_OK;				
			}
	};
	
} 


int main(int argc, const char** argv)
{
  TOPPProteinQuantifier t;
  return t.main(argc, argv);
}

/// @endcond
