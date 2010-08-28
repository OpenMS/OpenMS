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
// $Maintainer: David Wojnar $
// $Authors: David Wojnar $
// --------------------------------------------------------------------------
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FORMAT/MzDataFile.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/CONCEPT/Factory.h>
#include <OpenMS/FORMAT/MSPFile.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/COMPARISON/SPECTRA/BinnedSpectrum.h>
#include <OpenMS/COMPARISON/SPECTRA/SpectraSTSimilarityScore.h>
#include <OpenMS/COMPARISON/SPECTRA/CompareFouriertransform.h>
#include <OpenMS/COMPARISON/SPECTRA/ZhangSimilarityScore.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>

#include <ctime>
#include <vector>
#include <map>
#include <cmath>
using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
  @page TOPP_SpecLibSearcher SpecLibSearcher

  @brief Identifies peptide MS/MS spectra by spectral matching with a searchable spectral library.

<CENTER>
	<table>
		<tr>
			<td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
			<td VALIGN="middle" ROWSPAN=2> \f$ \longrightarrow \f$ SpecLibSearcher \f$ \longrightarrow \f$</td>
			<td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
		</tr>
		<tr>
			<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref UTILS_SpecLibCreator </td>
			<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_IDFilter or @n any protein/peptide processing tool</td>
		</tr>
	</table>
</CENTER>

	@experimental This TOPP-tool is not well tested and not all features might be properly implemented and tested.

	<B>The command line parameters of this tool are:</B>
	@verbinclude TOPP_SpecLibSearcher.cli
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPSpecLibSearcher
	: public TOPPBase
	{
	public:
		TOPPSpecLibSearcher()
		: TOPPBase("SpecLibSearcher","Identifies peptide MS/MS spectra by spectral matching with a searchable spectral library.")
		{
		}

	protected:
		void registerOptionsAndFlags_()
		{
			registerInputFileList_("in","<files>",StringList::create(""),"Input files");
			setValidFormats_("in",StringList::create("mzData"));
			registerInputFile_("lib","<file>","","searchable spectral library(MSP format)");
			registerOutputFileList_("out","<files>",StringList::create(""),"Output files. Have to be as many as input files");
			setValidFormats_("out",StringList::create("idXML"));
			registerDoubleOption_("precursor_mass_tolerance","<tolerance>",3,"Precursor mass tolerance, (Th)",false);
			registerIntOption_("round_precursor_to_integer","<number>",10,"many precursor m/z multipling number lead to the same number; are packed in the same vector for faster search.Should be higher for high-resolution data",false,true);
		//	registerDoubleOption_("fragment_mass_tolerance","<tolerance>",0.3,"Fragment mass error",false);

   //   registerStringOption_("precursor_error_units", "<unit>", "Da", "parent monoisotopic mass error units", false);
     // registerStringOption_("fragment_error_units", "<unit>", "Da", "fragment monoisotopic mass error units", false);
     // vector<String> valid_strings;
     // valid_strings.push_back("Da");
     // setValidStrings_("precursor_error_units", valid_strings);
     // setValidStrings_("fragment_error_units", valid_strings);
		//	registerIntOption_("min_precursor_charge", "<charge>", 1, "minimum precursor ion charge", false);
     // registerIntOption_("max_precursor_charge", "<charge>", 3, "maximum precursor ion charge", false);
     registerStringOption_("compare_function","<string>","ZhangSimilarityScore","function for similarity comparisson",false);
     PeakSpectrumCompareFunctor::registerChildren();
     setValidStrings_("compare_function",Factory<PeakSpectrumCompareFunctor>::registeredProducts());
		 registerIntOption_("top_hits","<number>",10,"save the first <number> top hits. For all type -1",false);
     addEmptyLine_();
     addText_("Filtering options. Most are especially useful when the query spectra are raw.");
     registerIntOption_("min_peaks","<number>",5, "required mininum number of peaks for a query spectrum",false);
     registerDoubleOption_("remove_peaks_below_threshold","<threshold>",2.01,"All peaks of a query spectrum with intensities below <threshold> will be zeroed.",false);
     registerIntOption_("max_peaks","<number>",150,"Use only the top <number> of peaks.",false);
     registerIntOption_("cut_peaks_below","<number>",1000,"Remove all peaks which are lower than 1/<number> of the highest peaks. Default equals all peaks which are lower than 0.001 of the maximum intensity peak",false);

			vector<String> all_mods;
			ModificationsDB::getInstance()->getAllSearchModifications(all_mods);
     	registerStringList_("fixed_modifications", "<mods>", StringList::create(""), "fixed modifications, specified using UniMod (www.unimod.org) terms, e.g. 'Carbamidomethyl (C)' or 'Oxidation (M)'", false);
			setValidStrings_("fixed_modifications", all_mods);

			registerStringList_("variable_modifications", "<mods>", StringList::create(""), "variable modifications, specified using UniMod (www.unimod.org) terms, e.g. 'Carbamidomethyl (C)' or 'Oxidation (M)'", false);
			setValidStrings_("variable_modifications", all_mods);
			addEmptyLine_();
			addText_("");
		}

		ExitCodes main_(int , const char**)
		{
			//-------------------------------------------------------------
			// parameter handling
			//-------------------------------------------------------------

			StringList in_spec = getStringList_("in");
			StringList out = getStringList_("out");
			String in_lib = getStringOption_("lib");
			String compare_function = getStringOption_("compare_function");
			Int precursor_mass_multiplier = getIntOption_("round_precursor_to_integer");
			Real precursor_mass_tolerance = getDoubleOption_("precursor_mass_tolerance");
			//Int min_precursor_charge = getIntOption_("min_precursor_charge");
			//Int max_precursor_charge = getIntOption_("max_precursor_charge");
			Real remove_peaks_below_threshold = getDoubleOption_("remove_peaks_below_threshold");
			UInt min_peaks = getIntOption_("min_peaks");
			UInt max_peaks= getIntOption_("max_peaks");
			Int cut_peaks_below = getIntOption_("cut_peaks_below");
			StringList fixed_modifications = getStringList_("fixed_modifications");
			StringList variable_modifications = getStringList_("variable_modifications");
			Int top_hits  = getIntOption_("top_hits");
			if(top_hits < -1 )
			{
				writeLog_("top_hits (should be  >= -1 )");
				return ILLEGAL_PARAMETERS;
			}
			//-------------------------------------------------------------
			// loading input
			//-------------------------------------------------------------
			if(out.size() != in_spec.size())
			{
				writeLog_("out (should be as many as input files)");
				return ILLEGAL_PARAMETERS;
			}

			time_t prog_time = time(NULL);
			MSPFile spectral_library;
			RichPeakMap query, library;
			//spectrum which will be identified
			MzDataFile spectra;
			spectra.setLogType(log_type_);

			time_t start_build_time = time(NULL);
			//-------------------------------------------------------------
			//building map for faster search
			//-------------------------------------------------------------

			//library containing already identified peptide spectra
			vector< PeptideIdentification > ids;
			spectral_library.load(in_lib,ids,library);

			map<Size, vector<PeakSpectrum> > MSLibrary;
			{
			RichPeakMap::iterator s;
			vector<PeptideIdentification>::iterator i;
			ModificationsDB* mdb = ModificationsDB::getInstance();
			for(s = library.begin(), i = ids.begin() ; s < library.end();++s,++i)
			{
				DoubleReal precursor_MZ = (*s).getPrecursors()[0].getMZ();
				Size MZ_multi = (Size)precursor_MZ*precursor_mass_multiplier;
				map<Size,vector<PeakSpectrum> >::iterator found;
				found = MSLibrary.find(MZ_multi);

				PeakSpectrum librar;
				bool variable_modifications_ok = true;
				bool fixed_modifications_ok = true;
				const AASequence& aaseq= i->getHits()[0].getSequence();
				//variable fixed modifications
				if(!fixed_modifications.empty())
				{
					for(Size i = 0; i< aaseq.size(); ++i)
					{
						const	Residue& mod  = aaseq.getResidue(i);
						for(Size s = 0; s < fixed_modifications.size();++s)
						{
						 if(mod.getOneLetterCode() ==mdb->getModification(fixed_modifications[s]).getOrigin() && fixed_modifications[s] != mod.getModification())
						 {
							 fixed_modifications_ok = false;
								break;
							}
						}
					}
				}
				//variable modifications
				if(aaseq.isModified() && (!variable_modifications.empty()))
				{
					for(Size i = 0; i< aaseq.size(); ++i)
					{
						if(aaseq.isModified(i))
						{
							const	Residue& mod  = aaseq.getResidue(i);
							for(Size s = 0; s < variable_modifications.size();++s)
							{
								if(mod.getOneLetterCode() ==mdb->getModification(variable_modifications[s]).getOrigin() && variable_modifications[s] != mod.getModification())
								{
									variable_modifications_ok = false;
									break;
								}
							}
						}
					}
				}
				if(variable_modifications_ok && fixed_modifications_ok)
				{
					PeptideIdentification& translocate_pid = *i;
					librar.getPeptideIdentifications().push_back(translocate_pid);
					librar.setPrecursors(s->getPrecursors());
					//library entry transformation
					for(UInt l = 0; l< s->size(); ++l)
					{
						Peak1D peak;
						if((*s)[l].getIntensity() >  remove_peaks_below_threshold)
						{
							const String& info = (*s)[l].getMetaValue("MSPPeakInfo");
							if(info[0] == '?')
							{
								peak.setIntensity(sqrt(0.2*(*s)[l].getIntensity()));
							}
							else
							{
								peak.setIntensity(sqrt((*s)[l].getIntensity()));
							}

							peak.setMZ((*s)[l].getMZ());
							peak.setPosition((*s)[l].getPosition());
							librar.push_back(peak);
						}
					}
					if(found != MSLibrary.end())
					{
						found->second.push_back(librar);
					}
					else
					{
						vector<PeakSpectrum> tmp;
						tmp.push_back(librar);
						MSLibrary.insert(make_pair(MZ_multi,tmp));
					}
				}
			}
			}
			time_t end_build_time = time(NULL);
		  cout<<"Time needed for preprocessing data: "<<(end_build_time - start_build_time)<<"\n";
			//compare function
			PeakSpectrumCompareFunctor* comparor = Factory<PeakSpectrumCompareFunctor>::create(compare_function);
			//-------------------------------------------------------------
			// calculations
			//-------------------------------------------------------------
			DoubleReal score;
			StringList::iterator in, out_file;
			for(in  = in_spec.begin() , out_file  = out.begin(); in < in_spec.end(); ++in,++out_file)
			{
				time_t start_time = time(NULL);
				spectra.load(*in,query);
				//Will hold valuable hits
				vector<PeptideIdentification> peptide_ids;
				vector<ProteinIdentification> protein_ids;
				// Write parameters to ProteinIdentifcation
				ProteinIdentification prot_id;
				//Parameters of identificaion
				prot_id.setIdentifier("test");
				prot_id.setSearchEngineVersion("SpecLibSearcher");
				prot_id.setDateTime(DateTime::now());
				prot_id.setScoreType(compare_function);
				ProteinIdentification::SearchParameters searchparam;
				searchparam.precursor_tolerance = precursor_mass_tolerance;
				prot_id.setSearchParameters(searchparam);
				/***********SEARCH**********/
				for(UInt j = 0; j < query.size(); ++j)
				{
					//Set identifier for each identifications
					PeptideIdentification pid;
					pid.setIdentifier("test");
					pid.setScoreType(compare_function);
					ProteinHit pr_hit;
					pr_hit.setAccession(j);
					prot_id.insertHit(pr_hit);
					//RichPeak1D to Peak1D transformation for the compare function query
					PeakSpectrum quer;
					bool peak_ok = true;
					query[j].sortByIntensity(true);
					DoubleReal min_high_intensity = 0;
					if(!query[j].empty())
					{
						min_high_intensity = (1/cut_peaks_below)*query[j][0].getIntensity();
					}
					query[j].sortByPosition();
					for(UInt k = 0; k < query[j].size() && k < max_peaks; ++k)
					{
						if(query[j][k].getIntensity() >  remove_peaks_below_threshold && query[j][k].getIntensity() >= min_high_intensity )
						{
							Peak1D peak;
							peak.setIntensity(sqrt(query[j][k].getIntensity()));
							peak.setMZ(query[j][k].getMZ());
							peak.setPosition(query[j][k].getPosition());
							quer.push_back(peak);
						}
					}
					if(quer.size() >= min_peaks)
					{
						peak_ok = true;
					}
					else
					{
						peak_ok = false;
					}
					DoubleReal query_MZ = query[j].getPrecursors()[0].getMZ();
					if(peak_ok)
					{
						bool charge_one = false;
						Int percent = (Int)Math::round((query[j].size()/100.0) * 3.0);
						Int margin  = (Int)Math::round((query[j].size()/100.0)* 1.0);
						for(vector<RichPeak1D>::iterator peak = query[j].end()-1; percent >= 0; --peak , --percent)
						{
							if(peak->getMZ() < query_MZ)
							{
								break;
							}
						}
						if(percent > margin)
						{
							charge_one = true;
						}
						Real min_MZ = (query_MZ - precursor_mass_tolerance) *precursor_mass_multiplier;
						Real max_MZ = (query_MZ + precursor_mass_tolerance) *precursor_mass_multiplier;
						for(Size mz = (Size)min_MZ; mz <= ((Size)max_MZ)+1; ++mz)
						{
							map<Size, vector<PeakSpectrum> >::iterator found;
								found = MSLibrary.find(mz);
								if(found != MSLibrary.end())
								{
									vector<PeakSpectrum>& library = found->second;
									for(Size i = 0; i < library.size(); ++i)
									{
										Real this_MZ  = library[i].getPrecursors()[0].getMZ()*precursor_mass_multiplier;
										if(this_MZ >= min_MZ && max_MZ >= this_MZ && ((charge_one == true && library[i].getPeptideIdentifications()[0].getHits()[0].getCharge() == 1) ||charge_one == false ))
										{
											PeptideHit hit = library[i].getPeptideIdentifications()[0].getHits()[0];
											PeakSpectrum& librar = library[i];
											//Special treatment for SpectraST score as it computes a score based on the whole library
											if(compare_function == "SpectraSTSimilarityScore")
											{
												SpectraSTSimilarityScore* sp= static_cast<SpectraSTSimilarityScore*>(comparor);
												BinnedSpectrum quer_bin = sp->transform(quer);
												BinnedSpectrum librar_bin = sp->transform(librar);
												score = (*sp)(quer,librar);//(*sp)(quer_bin,librar_bin);
												double dot_bias = sp->dot_bias(quer_bin,librar_bin,score);
												hit.setMetaValue("DOTBIAS",dot_bias);
											}
											else
											{
												if(compare_function =="CompareFouriertransform")
												{
													CompareFouriertransform* ft = static_cast<CompareFouriertransform*>(comparor);
													ft->transform(quer);
													ft->transform(librar);
												}
												score = (*comparor)(quer,librar);
											}

											DataValue RT(library[i].getRT());
											DataValue MZ(library[i].getPrecursors()[0].getMZ());
											hit.setMetaValue("RT",RT);
											hit.setMetaValue("MZ",MZ);
											hit.setScore(score);
											hit.addProteinAccession(pr_hit.getAccession());
											pid.insertHit(hit);
										}
									}
								}
							}
						}
						pid.setHigherScoreBetter(true);
						pid.sort();
						if(compare_function =="SpectraSTSimilarityScore")
						{
							if(!pid.empty() && !pid.getHits().empty())
							{
								vector<PeptideHit> final_hits;
								final_hits.resize(pid.getHits().size());
								SpectraSTSimilarityScore* sp= static_cast<SpectraSTSimilarityScore*>(comparor);
								Size runner_up = 1;
								for( ; runner_up < pid.getHits().size() ; ++runner_up)
								{
									if(pid.getHits()[0].getSequence().toUnmodifiedString() != pid.getHits()[runner_up].getSequence().toUnmodifiedString()|| runner_up > 5)
									{
										break;
									}
								}
								double delta_D = sp->delta_D(pid.getHits()[0].getScore(), pid.getHits()[runner_up].getScore());
								for(Size s = 0; s < pid.getHits().size();++s)
								{
									final_hits[s] = pid.getHits()[s];
									final_hits[s].setMetaValue("delta D",delta_D);
									final_hits[s].setMetaValue("dot product",pid.getHits()[s].getScore());
									final_hits[s].setScore(sp->compute_F(pid.getHits()[s].getScore(),delta_D,pid.getHits()[s].getMetaValue("DOTBIAS")));

									//final_hits[s].removeMetaValue("DOTBIAS");
							}
							pid.setHits(final_hits);
							pid.sort();
							pid.setMetaValue("MZ",query[j].getPrecursors()[0].getMZ());
							pid.setMetaValue("RT",query_MZ);
						}
					}
					if(top_hits != -1 && (UInt)top_hits < pid.getHits().size())
					{
						vector<PeptideHit> hits;
						hits.resize(top_hits);
						for(Size i = 0; i < (UInt)top_hits; ++i)
						{
							hits[i] = pid.getHits()[i];
						}
						pid.setHits(hits);
					}
					peptide_ids.push_back(pid);
				}
				protein_ids.push_back(prot_id);
				//-------------------------------------------------------------
				// writing output
				//-------------------------------------------------------------
				IdXMLFile id_xml_file;
				id_xml_file.store(*out_file,protein_ids,peptide_ids);
				time_t end_time = time(NULL);
				cout<<"Search time: "<<difftime(end_time,start_time)<<" seconds for "<<*in<<"\n";
			}
			time_t end_time = time(NULL);
			cout<<"Total time: "<<difftime(end_time,prog_time)<<" secconds\n";
			return EXECUTION_OK;
		}

	};




int main( int argc, const char** argv )
{
    TOPPSpecLibSearcher tool;
    return tool.main(argc,argv);
}

/// @endcond
