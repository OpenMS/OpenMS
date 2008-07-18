// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_MSPFILE_H
#define OPENMS_FORMAT_MSPFILE_H

#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/FORMAT/HANDLERS/MzDataHandler.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <cctype>
#include <fstream>

namespace OpenMS
{
	/**
		@brief File adapter for MzData files
	
		@todo add ProgressLogger to File? (Andreas)
		@todo is TextFile a good idea? (Andreas)
		@ingroup FileIO
	*/
	class MSPFile /*
			public ProgressLogger*/
	{
		public:

			///Default constructor
			MSPFile();

			///Destructor
			virtual ~MSPFile();
			
			/**
				@brief Loads a map from a MSPFile file.

				@p exp has to be a MSExperiment or have the same interface.
				@param filename the filename of the experiment
				@param read_headers if set to true the header information is also read from the file and stored into metadata
				@param ids output parameter which contains the peptide identifications from the spectra anntations
				@param exp output parameter which contains the spectra 
				@throw FileNotFound is thrown if the file could not be found
				@throw ParseError is thrown if the given file could not be parsed
			*/
			template <typename MapType>
			void load(const String& filename, std::vector<PeptideIdentification>& ids, MapType& exp, bool read_headers = false)
			{
				if (!File::exists(filename))
				{
					throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
				}
				if (!File::readable(filename))
				{
					throw Exception::FileNotReadable(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
				}
				String line;
				std::ifstream is(filename.c_str());

				Map<String, double> mod_to_mass;
				mod_to_mass["Oxidation"] = 15.994915;
				mod_to_mass["Carbamidomethyl "] = 57.02146;
				mod_to_mass["ICAT_light"] = 227.12;
				mod_to_mass["ICAT_heavy"] = 236.12;
				mod_to_mass["AB_old_ICATd0"] = 442.20;
				mod_to_mass["AB_old_ICATd8"] = 450.20;
				mod_to_mass["Acetyl"] = 42.0106;
				mod_to_mass["Deamidation"] = 0.9840;
				mod_to_mass["Pyro-cmC"] = -17.026549;
				mod_to_mass["Pyro-glu"] = -17.026549;
				mod_to_mass["Pyro_glu"] = -18.010565;
				mod_to_mass["Amide"] = -0.984016;
				mod_to_mass["Phospho"] = 79.9663;
				mod_to_mass["Methyl"] = 14.0157;
				mod_to_mass["Carbamyl"] = 43.00581;
				
				typename MapType::SpectrumType spec;

				while (getline(is, line))
				{
					if (line.hasPrefix("Name:"))
					{
						std::vector<String> split, split2;
						line.split(' ', split);
						split[1].split('/', split2);
						String peptide = split2[0];
						// remove damn (O), also defined in 'Mods=' comment
						peptide.substitute("(O)", "");
						PeptideIdentification id;
						id.insertHit(PeptideHit(0, 0, split2[1].toInt(), peptide));
						ids.push_back(id);
					}
					if (line.hasPrefix("MW:"))
					{
						// skip that as it is not necessary and might not be available at all
					}
					if (line.hasPrefix("Comment:"))
					{
						// slow, but we need the modifications from the header
						std::vector<String> split;
						line.split(' ', split);
						for (std::vector<String>::const_iterator it = split.begin(); it != split.end(); ++it)
						{
							if (*it == "Mods=0")
							{
								break;
							}
							if (it->hasPrefix("Mods="))
							{
								String mods = it->suffix('=');
								// e.g. Mods=2/7,K,Carbamyl/22,K,Carbamyl
								std::vector<String> mod_split;
								mods.split('/', mod_split);
								AASequence peptide = ids.back().getHits().begin()->getSequence();
								for (UInt i = 1; i <= (UInt)mod_split[0].toInt(); ++i)
								{
									std::vector<String> single_mod;
									mod_split[i].split(',', single_mod);
									
									// fix errornous names
									String mod_name = single_mod[2];
									if (mod_name == "Pyro-glu")
									{
										mod_name += " from " + single_mod[1];
									}
									String psi_mod = ModificationsDB::getInstance()->getModification(single_mod[1], mod_name).getId();
									peptide.setModification(single_mod[0].toInt(), psi_mod);
								}
								std::vector<PeptideHit> hits(ids.back().getHits());
								hits.begin()->setSequence(peptide);
								ids.back().setHits(hits);
							}
						}
										
						if (read_headers)
						{
							parseHeader_(line, spec);
						}
					}
					if (line.hasPrefix("Num peaks:"))
					{
						while (getline(is, line) && line.size() > 0 && std::isdigit(line[0]))
						{
							std::vector<String> split;
							line.split('\t', split);
							typename MapType::SpectrumType::PeakType peak;
							if (split.size() != 3)
							{
								throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, line, "not <mz> <intensity> <comment>");
							}
							peak.setMZ(split[0].toFloat());
							peak.setIntensity(split[1].toFloat());
							peak.setMetaValue("MSPPeakInfo", split[2]);
							spec.push_back(peak);
						}
						
						exp.push_back(spec);
						spec.clear();
					}
				}
			}

			/**
				@brief Stores a map in a MSPFile file.
				
				@param filename the filename of the MSPFile which should be written
				@param exp has to be a MSExperiment or have the same interface.
				@throw UnableToCreateFile is thrown if the given file could not be created
			*/
			template <typename MapType> void store(const String& filename, const MapType& exp) const
			{
				throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
			}

			protected:
				
				template <typename SpectrumType> void parseHeader_(const String& header, SpectrumType& spec)
				{
					// first header from std_protein of NIST spectra DB
					// Spec=Consensus Pep=Tryptic Fullname=R.AAANFFSASCVPCADQSSFPK.L/2 Mods=0 Parent=1074.480 Inst=it Mz_diff=0.500 Mz_exact=1074.4805 Mz_av=1075.204 Protein="TRFE_BOVIN" Organism="Protein Standard" Se=2^X23:ex=3.1e-008/1.934e-005,td=5.14e+007/2.552e+019,sd=0/0,hs=45.8/5.661,bs=6.3e-021,b2=1.2e-015,bd=5.87e+020^O22:ex=3.24e-005/0.0001075,td=304500/5.909e+297,pr=3.87e-007/1.42e-006,bs=1.65e-301,b2=1.25e-008,bd=1.3e+299 Sample=1/bovine-serotransferrin_cam,23,26 Nreps=23/34 Missing=0.3308/0.0425 Parent_med=1074.88/0.23 Max2med_orig=22.1/9.5 Dotfull=0.618/0.029 Dot_cons=0.728/0.040 Unassign_all=0.161 Unassigned=0.000 Dotbest=0.70 Naa=21 DUScorr=2.3/3.8/0.61 Dottheory=0.86 Pfin=4.3e+010 Probcorr=1 Tfratio=8e+008 Specqual=0.0
					
					std::vector<String> split;
					header.split(' ', split);
				
					
					for (std::vector<String>::const_iterator it = split.begin(); it != split.end(); ++it)
					{
						std::vector<String> split2;
						String tmp = *it;
						tmp.trim();
						tmp.split('=', split2);
						if (split2.size() == 2)
						{
							spec.setMetaValue(split2[0], split2[1]);
						}
					}
				}
		
	};

} // namespace OpenMS

#endif // OPENMS_FORMAT_MSPFILE_H
