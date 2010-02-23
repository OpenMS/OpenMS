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
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/SVOutStream.h>

#include <numeric> // for 'accumulate'

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page TOPP_ProteinQuantifier ProteinQuantifier
	
	@brief Application to compute protein abundances from annotated feature/consensus maps

	Peptide/protein IDs from multiple identification runs can be handled, but will not be differentiated (i.e. protein accessions for a peptide will be accumulated over all identification runs).


	<B>The command line parameters of this tool are:</B>
	@verbinclude TOPP_ProteinQuantifier.cli
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
			TOPPBase("ProteinQuantifier", "Compute protein abundances")
      {
      }

	protected:

		typedef map<UInt64, DoubleReal> sample_abundances; // abundance per sample

		struct peptide_data
		{
			map<Int, sample_abundances> abundances; // charge -> sample -> abundance
			sample_abundances total_abundances; // sample -> total abundance
			set<String> accessions; // protein accessions
		};
		typedef map<AASequence, peptide_data> peptide_quant; // by peptide sequence

		struct protein_data
		{
			// peptide -> sample -> abundance:
			map<AASequence, sample_abundances> abundances;
			sample_abundances total_abundances; // sample -> total abundance
		};
		typedef map<String, protein_data> protein_quant; // by protein accession


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


		PeptideHit getAnnotation_(vector<PeptideIdentification>& peptides)
			{
				// only the best peptide hit in each peptide ID is taken into account!
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

		
		void quantifyFeature_(const FeatureHandle& feature, const PeptideHit& hit, 
													peptide_quant& quant)
			{
				if (hit == PeptideHit()) 
				{
					return; // annotation for the feature ambiguous or missing
				}
				const AASequence& seq = hit.getSequence();
				// note: sample number stored in "feature::map_index_"
				quant[seq].abundances[feature.getCharge()][feature.getMapIndex()] += 
					feature.getIntensity(); // new map element is initialized with 0
				// maybe TODO: check charge of peptide hit against charge of feature
				quant[seq].accessions.insert(hit.getProteinAccessions().begin(),
																		 hit.getProteinAccessions().end());
			}
		

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
/*
						Size best_size = 0;
						Int best_charge;
						DoubleReal best_abundance;
						for (map<Int, sample_abundances>::iterator ab_it = 
									 q_it->second.abundances.begin(); ab_it != 
									 q_it->second.abundances.end(); ++ab_it)
						{
							Size size = ab_it->second.size();
							if (size < best_size) continue;
							DoubleReal current_abundance = 0.0;
							for (sample_abundances::iterator samp_it = ab_it->second.begin();
									 samp_it != ab_it->second.end(); ++samp_it)
							{
								current_abundance += samp_it->second;
							}
							if ((size > best_size) || (current_abundance > best_abundance))
							{
								best_size = size;
								best_charge = ab_it->first;
								best_abundance = current_abundance;
							}
						}
*/
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
		

		void quantifyProteins_(const peptide_quant& pep_quant, 
													 protein_quant& prot_quant)
			{
				for (peptide_quant::const_iterator pep_it = pep_quant.begin(); 
						 pep_it != pep_quant.end(); ++pep_it)
				{
					if (pep_it->second.accessions.size() == 1) // proteotypic peptide
					{
						String accession = *(pep_it->second.accessions.begin());
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
				bool include_fewer = getFlag_("include_fewer"), 
					fix_peptides = getFlag_("fix_peptides");

				for (protein_quant::iterator prot_it = prot_quant.begin();
						 prot_it != prot_quant.end(); ++prot_it)
				{
					vector<AASequence> peptides; // peptides selected for quantification
					if (fix_peptides) // consider only "top" best peptides
					{
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
						if (!include_fewer && (ab_it->second.size() < top))
						{
							continue; // not enough proteotypic peptides
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


		void writeProteinTable_(SVOutStream& out, protein_quant& quant,
														const vector<UInt64> samples)
			{
				// write header:
				out << "protein" << "n_peptides";
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

				bool include_fewer = getFlag_("include_fewer");
				Size top = getIntOption_("top");
				for (protein_quant::iterator q_it = quant.begin(); q_it != quant.end();
						 ++q_it)
				{
					Size n_peptide = q_it->second.abundances.size();
					if (include_fewer || (n_peptide >= top))
					{
						out << q_it->first << n_peptide;
						for (vector<UInt64>::const_iterator samp_it = samples.begin(); 
								 samp_it != samples.end(); ++samp_it)
						{
							out << q_it->second.total_abundances[*samp_it];
						}
						out << endl;
					}
				}
			}

		
		void writeComments_(SVOutStream& out, 
												const ConsensusMap::FileDescriptions files,
												const bool proteins=true)
			{
				String what = (proteins ? "Protein" : "Peptide"); 
				bool old = out.modifyStrings(false);
				out << "#" + what + " abundances computed from file '" + 
					getStringOption_("in") + "'" << endl;
				String params = "#Parameters: top=" + String(getIntOption_("top")) + 
					", average=" + getStringOption_("average");
				StringList flags = StringList::create("include_fewer,filter_charge");
				if (files.size() > 1)
				{
					flags << "fix_peptides" << "normalize";
				}
				for (StringList::Iterator it = flags.begin(); it != flags.end(); ++it)
				{
					if (getFlag_(*it)) params += ", " + *it;
				}
				out << params << endl;
				if (files.size() > 1)
				{
					String desc = 
						"#Files/samples associated with abundance values below: ";
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

	
		void registerOptionsAndFlags_()
      {
				registerInputFile_("in", "<file>", "", "Input file");
				setValidFormats_("in", StringList::create("featureXML,consensusXML"));
        registerOutputFile_("out", "<file>", "", "Output file for protein abundances");
				registerOutputFile_("peptide_out", "<file>", "", "Output file for peptide abundances", false);
				registerIntOption_("top", "<number>", 3, "Calculate protein abundance from this number of proteotypic peptides (best first; '0' for all)", false);
				setMinInt_("top", 0);
				registerFlag_("include_fewer", "Include results for proteins with fewer than 'top' proteotypic peptides");
				registerStringOption_("average", "<method>", "median", "Averaging method for computing protein abundances from peptide abundances", false);
				setValidStrings_("average", StringList::create("median,mean,sum"));
				registerFlag_("filter_charge", "Set this flag to distinguish between charge states of a peptide. For peptides, abundances will be reported separately for each charge;\nfor proteins, abundances will be computed based only on the most prevalent charge of each peptide (this may increase reproducibility between samples).\nBy default, abundances are summed over all charge states.");
				addEmptyLine_();
        addText_("Options for consensusXML files:");
				registerFlag_("fix_peptides", "Use the same peptides for protein quantification across all samples");
				registerFlag_("normalize", "Scale peptide abundances so that medians of all samples are equal");
				addEmptyLine_();
				addText_("Output formatting options:");
				registerStringOption_("separator", "<sep>", "", "The used separator character(s); if not set the 'tab' character is used", false);
				registerStringOption_("quoting", "<method>", "double", "Method for quoting of strings: 'none' for no quoting, 'double' for quoting with doubling of embedded quotes,\n'escape' for quoting with backslash-escaping of embedded quotes", false);
				setValidStrings_("quoting", StringList::create("none,double,escape"));
				registerStringOption_("replacement", "<string>", "_", "Used to replace occurrences of the separator in strings before writing, if 'quoting' is 'none'", false);
      }

		
		ExitCodes main_(int, const char**)
      {
				String in = getStringOption_("in"), out = getStringOption_("out"), 
					peptide_out = getStringOption_("peptide_out");
				FileTypes::Type in_type = FileHandler::getType(in);

				String separator = getStringOption_("separator"), 
					replacement = getStringOption_("replacement"), 
					quoting = getStringOption_("quoting");
        if (separator == "") separator = "\t";
				String::QuotingMethod quoting_method;
				if (quoting == "none") quoting_method = String::NONE;
				else if (quoting == "double") quoting_method = String::DOUBLE;
				else quoting_method = String::ESCAPE;

				bool normalize = false;
				vector<UInt64> samples;
				ConsensusMap::FileDescriptions files;
				peptide_quant pep_quant;
				
				if (in_type == FileTypes::FEATUREXML)
				{
					FeatureMap<> features;
          FeatureXMLFile().load(in, features);
					samples.push_back(0);

					for (FeatureMap<>::Iterator feat_it = features.begin(); 
							 feat_it != features.end(); ++feat_it)
					{
						PeptideHit hit = 
							getAnnotation_(feat_it->getPeptideIdentifications());
						FeatureHandle handle(0, *feat_it);
						quantifyFeature_(handle, hit, pep_quant);
					}
				}

				else // consensusXML
				{
					ConsensusMap consensus;
          ConsensusXMLFile().load(in, consensus);

					normalize = getFlag_("normalize");
					files = consensus.getFileDescriptions(); 
					for (ConsensusMap::FileDescriptions::Iterator file_it = files.begin();
							 file_it != files.end(); ++file_it)
					{
						samples.push_back(file_it->first);
					}

					for (ConsensusMap::Iterator cons_it = consensus.begin(); 
							 cons_it != consensus.end(); ++cons_it)
					{
						PeptideHit hit = 
							getAnnotation_(cons_it->getPeptideIdentifications());
						for (set<FeatureHandle>::iterator feat_it = 
									 cons_it->getFeatures().begin(); feat_it !=
									 cons_it->getFeatures().end(); ++feat_it)
						{
							quantifyFeature_(*feat_it, hit, pep_quant);
						}
					}
				}

				quantifyPeptides_(pep_quant);
				if (normalize) normalizePeptides_(pep_quant);

				// output:
				if (!peptide_out.empty())
				{
					ofstream outstr(peptide_out.c_str());
					SVOutStream output(outstr, separator, replacement, quoting_method);
					writeComments_(output, files, false);
					writePeptideTable_(output, pep_quant, samples);
					outstr.close();
				}

				protein_quant prot_quant;
				quantifyProteins_(pep_quant, prot_quant);
				ofstream outstr(out.c_str());
				SVOutStream output(outstr, separator, replacement, quoting_method);
				writeComments_(output, files);
				writeProteinTable_(output, prot_quant, samples);
				outstr.close();

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
