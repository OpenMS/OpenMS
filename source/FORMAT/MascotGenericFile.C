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
// $Maintainer: Andreas Bertsch $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/MascotGenericFile.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>

using namespace std;

namespace OpenMS 
{

	MascotGenericFile::MascotGenericFile()
		: ProgressLogger(),
			DefaultParamHandler("MascotInfile")
	{
		defaults_.setValue("database", "MSDB", "Name of the sequence database");
		defaults_.setValue("search_type", "MIS", "Name of the search type for the query", StringList::create("advanced"));
		defaults_.setValidStrings("search_type", StringList::create("MIS,SQ,PMF"));
		defaults_.setValue("enzyme", "Trypsin", "Name of enzyme used for the digestion");
		defaults_.setValue("instrument", "Default", "Instrument definition which specifies the fragmentation rules");
		defaults_.setValue("missed_cleavages", 1, "Number of missed cleavages allowed for the enzyme");
		defaults_.setMinInt("missed_cleavages", 0);
		defaults_.setValue("precursor_mass_tolerance", 3.0, "Tolerance of the precursor peaks");
		defaults_.setMinFloat("precursor_mass_tolerance", 0.0);
		defaults_.setValue("precursor_error_units", "Da", "Units of the precursor mass tolerance");
		defaults_.setValidStrings("precursor_error_units", StringList::create("%,ppm,mmu,Da"));
		defaults_.setValue("fragment_mass_tolerance", 0.3, "Tolerance of the peaks in the fragment spectrum");
		defaults_.setMinFloat("fragment_mass_tolerance", 0.0);
		defaults_.setValue("fragment_error_units", "Da", "Units of the fragment peaks tolerance");
		defaults_.setValidStrings("fragment_error_units", StringList::create("mmu,Da"));
		defaults_.setValue("charges", "1,2,3", "Allowed charge states, given as a comma separated list of integers");
		defaults_.setValue("taxonomy", "All entries", "Taxonomy specification of the sequences");
		defaults_.setValue("fixed_modifications", StringList::create(""), "List of fixed modifications, according to UniMod definitions.");
		vector<String> all_mods;
		ModificationsDB::getInstance()->getAllSearchModifications(all_mods);
		defaults_.setValidStrings("fixed_modifications", all_mods);
		defaults_.setValue("variable_modifications", StringList::create(""), "Variable modifications given as UniMod definitions.");
		defaults_.setValidStrings("variable_modifications", all_mods);
	
		defaults_.setValue("mass_type", "monoisotopic", "Defines the mass type, either monoisotopic or average");
		defaults_.setValidStrings("mass_type", StringList::create("monoisotopic,average"));
		defaults_.setValue("number_of_hits", 10, "Number of hits which should be returned, if 0 AUTO mode is enabled.");
		defaults_.setMinInt("number_of_hits", 0);
		defaults_.setValue("skip_spectrum_charges", "false", "Sometimes precursor charges are given for each spectrum but are wrong, setting this to 'true' does not write any charge information to the spectrum, the general charge information is however kept.");
		defaults_.setValidStrings("skip_spectrum_charges", StringList::create("true,false"));
		
		defaults_.setValue("search_title", "OpenMS_search", "Sets the title of the search.", StringList::create("advanced"));
		defaults_.setValue("username", "OpenMS", "Sets the username which is mentioned in the results file.", StringList::create("advanced"));
		defaults_.setValue("format", "Mascot generic", "Sets the format type of the peak list, this should not be changed.", StringList::create("advanced"));
		defaults_.setValue("form_version", "1.01", "Sets the version of the peak list format, this should not be changed.", StringList::create("advanced"));
		defaults_.setValue("peaklists_only", "false", "Skip the parameter header and the MIME parts; just write the peak lists", StringList::create("advanced"));
		defaults_.setValue("boundary", "GZWgAaYKjHFeUaLOLEIOMq", "MIME boundary", StringList::create("advanced"));
		defaults_.setValidStrings("peaklists_only", StringList::create("true,false"));

		defaultsToParam_();
	}
	
	MascotGenericFile::~MascotGenericFile()
	{
		
	}
	
	void MascotGenericFile::store(const String& filename, const PeakMap& experiment)
	{
		if (!File::writable(filename))
		{
			throw Exception::FileNotWritable(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
		}
		ofstream os(filename.c_str());
		store(os, filename, experiment);
		os.close();
	}

	void MascotGenericFile::store(ostream& os, const String& filename, const PeakMap& experiment)
	{
		if (!param_.getValue("peaklists_only").toBool())
		{
			writeHeader_(os);
		}
		writeMSExperiment_(os, filename, experiment);

		//close file
		if (!param_.getValue("peaklists_only").toBool())
		{
			os << "\n\n" << "--" << param_.getValue("boundary") << "--\n";
		}
	}
	
	void MascotGenericFile::writeParameterHeader_(const String& name, ostream& os)
	{
		os << "--" << param_.getValue("boundary") << "\n" << "Content-Disposition: form-data; name=\"" << name << "\"" << "\n\n";
	}
		
	void MascotGenericFile::writeHeader_(ostream& os)
	{
		// search title
		if (param_.getValue("search_title") != "")
		{
			writeParameterHeader_("COM", os);
			os << param_.getValue("search_title") << "\n";
		}
		
		// user name
		writeParameterHeader_("USERNAME", os);
		os  << param_.getValue("username") << "\n";

		// format
		writeParameterHeader_("FORMAT", os);
		os << param_.getValue("format") << "\n";

		// precursor mass tolerance unit : Da
		writeParameterHeader_("TOLU", os);
		os << param_.getValue("precursor_error_units") << "\n";

		// ion mass tolerance unit : Da
		writeParameterHeader_("ITOLU", os);
		os << param_.getValue("fragment_error_units") << "\n";

		// format version
		writeParameterHeader_("FORMVER", os);
		os << param_.getValue("form_version") << "\n";
		
		//db name
		writeParameterHeader_("DB", os);
		os << param_.getValue("database") << "\n";
		
		// search type
		writeParameterHeader_("SEARCH", os);
		os << param_.getValue("search_type") << "\n";

		// number of peptide candidates in the list
		writeParameterHeader_("REPORT", os);
		UInt num_hits((UInt)param_.getValue("number_of_hits"));
		if (num_hits != 0)
		{
			os << param_.getValue("number_of_hits") << "\n";
		}
		else
		{
			os << "AUTO" << "\n";
		}
		
		//cleavage enzyme
		writeParameterHeader_("CLE", os);
		os << param_.getValue("enzyme") << "\n";

		//average/monoisotopic
		writeParameterHeader_("MASS", os);
		os << param_.getValue("mass_type") << "\n";
		
		//fixed modifications	
		StringList fixed_mods((StringList)param_.getValue("fixed_modifications"));
		for(StringList::const_iterator it = fixed_mods.begin(); it != fixed_mods.end(); ++it)
		{
			writeParameterHeader_("MODS", os);
			os << *it << "\n";
		}

		//variable modifications
		StringList var_mods((StringList)param_.getValue("variable_modifications"));
		for(StringList::const_iterator it = var_mods.begin(); it != var_mods.end(); ++it)
		{
			writeParameterHeader_("IT_MODS", os);
			os << *it << "\n";
		}

		//instrument
		writeParameterHeader_("INSTRUMENT", os);
		os << param_.getValue("instrument") << "\n";
		
		//missed cleavages
		writeParameterHeader_("PFA", os);
		os << param_.getValue("missed_cleavages") << "\n";

		//precursor mass tolerance
		writeParameterHeader_("TOL", os);
		os << param_.getValue("precursor_mass_tolerance") << "\n";

		//ion mass tolerance_
		writeParameterHeader_("ITOL", os);
		os << param_.getValue("fragment_mass_tolerance") << "\n";

		//taxonomy
		writeParameterHeader_("TAXONOMY", os);
		os << param_.getValue("taxonomy") << "\n";

		//charge
		writeParameterHeader_("CHARGE", os);
		os << param_.getValue("charges") << "\n";
	}
	
	void MascotGenericFile::writeSpectrum_(ostream& os,	const PeakSpectrum& spec)
	{
		Precursor precursor;
		if (spec.getPrecursors().size()>0)
		{
			precursor = spec.getPrecursors()[0];
		}
		if (spec.getPrecursors().size()>1)
		{
			cerr << "Warning: The spectrum written to Mascot file has more than one precursor. The first precursor is used!\n";
		}
		DoubleReal mz(precursor.getMZ()), rt(spec.getRT());
		int charge(precursor.getCharge());

		if (mz == 0)
		{
			//retention time
			cout << "No precursor m/z information for spectrum with rt: " << rt << " present, skipping spectrum!\n";
		}
		else
		{
			os << "\n";
			os << "BEGIN IONS\n";
			os << "TITLE=" << precisionWrapper(mz) << "_" << precisionWrapper(rt) << "\n";
			os << "PEPMASS=" << precisionWrapper(mz) <<  "\n";
			os << "RTINSECONDS=" << precisionWrapper(rt) << "\n";

			bool skip_spectrum_charges(param_.getValue("skip_spectrum_charges").toBool());
			if (charge != 0)
			{
				if (!skip_spectrum_charges)
				{
					os << "CHARGE=" << charge << "\n";
				}
			}

			os << "\n";

			for (PeakSpectrum::const_iterator it = spec.begin() ; it != spec.end();++it)
			{
				os << precisionWrapper(it->getMZ()) << " " << precisionWrapper(it->getIntensity()) << "\n";
			}
			os << "END IONS\n";
		}
	}

	void MascotGenericFile::writeMSExperiment_(ostream& os, const String& filename, const PeakMap& experiment)
	{
		if (!param_.getValue("peaklists_only").toBool())
		{
			os << "--" << param_.getValue("boundary") << "\n" << "Content-Disposition: form-data; name=\"FILE\"; filename=\"" << filename << "\"" << "\n" << "\n";
		}

		for (Size i = 0; i < experiment.size(); i++)
		{
			if (experiment[i].getMSLevel() == 0)
			{
				cout << "MascotGenericFile: MSLevel is set to 0, ignoring this spectrum!" << "\n";
			}
			
			if (experiment[i].getMSLevel() == 2)
			{
				writeSpectrum_(os, experiment[i]);
			}
		}
	}

	bool MascotGenericFile::getNextSpectrum_(istream& is, vector<pair<DoubleReal, DoubleReal> >& spectrum, UInt& charge, DoubleReal& precursor_mz, DoubleReal& precursor_int, DoubleReal& rt, String& title)
	{
		bool ok(false);
		spectrum.clear();
		charge = 0;
		precursor_mz = 0;
		precursor_int = 0;
		
		String line;
		// seek to next peak list block
		while (getline(is, line, '\n'))
		{
			// found peak list block?
			if (line.trim() == "BEGIN IONS")
			{
				ok = false;
				while (getline(is, line, '\n'))
				{
					// parse precursor position
					if (line.trim().hasPrefix("PEPMASS"))
					{
						String tmp = line.substr(8);
						tmp.substitute('\t', ' ');
						vector<String> split;
						tmp.split(' ', split);
						if (split.size() == 1)
						{
							precursor_mz = split[0].trim().toDouble();
						}
						else
						{
							if (split.size() == 2)
							{
								precursor_mz = split[0].trim().toDouble();
								precursor_int = split[1].trim().toDouble();
							}
							else
							{
								throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "cannot parse PEPMASS: " + line, "");
							}
						}
					}
					if (line.trim().hasPrefix("CHARGE"))
					{
						String tmp = line.substr(7);
						tmp.remove('+');
						charge = tmp.toInt();
					}
					if (line.trim().hasPrefix("RTINSECONDS"))
					{
						String tmp = line.substr(12);
						rt = tmp.toDouble();
					}
					if (line.trim().hasPrefix("TITLE"))
					{
						// test if we have a line like "TITLE= Cmpd 1, +MSn(595.3), 10.9 min"
						if (line.hasSubstring("min"))
						{
							try
							{
								vector<String> split;
								line.split(',', split);
								if (split.size() > 0)
								{
									for (Size i = 0; i != split.size(); ++i)
									{
										if (split[i].hasSubstring("min"))
										{
											vector<String> split2;
											split[i].trim().split(' ', split2);
											if (split2.size() > 0)
											{
												rt = split2[0].trim().toDouble() * 60.0;
											}
										}
									}
								}
							}
              catch (Exception::BaseException& /*e*/)
              {
                // just do nothing and write the whole title to spec
                vector<String> split;
                line.split('=', split);
                if (split.size() >= 2)
                {
                  title = split[1];
                }
              }
						}
						else // just write the title as metainfo to the spectrum
						{
							vector<String> split;
							line.split('=', split);
							if (split.size() == 2)
							{
								title = split[1];
							}
							// TODO concatenate the other parts if the title contains additional '=' chars
						}
					}
					if (line.trim().size() > 0 && isdigit(line[0]))
					{
						do
						{
							line.substitute('\t', ' ');
							vector<String> split;
							line.split(' ', split);
							if (split.size() == 2)
							{
								spectrum.push_back(make_pair(split[0].toDouble(), split[1].toDouble()));
							}
							else
							{
								if (split.size() == 3)
								{
									spectrum.push_back(make_pair(split[0].toDouble(), split[1].toDouble()));
									// @improvement add meta info e.g. charge, name... (Andreas)
								}
								else
								{
									throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "the line (" + line + ") should contain m/z and intensity value separated by whitespace!", "");
								}
							}
						}while(getline(is, line, '\n') && line.trim() != "END IONS");
						if (line.trim() == "END IONS")
						{
							// found spectrum
							return true;
						}
						else
						{
							throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Found \"BEGIN IONS\" but not the corresponding \"END IONS\"!", "");
						}
					}
				}
			}
		}
		
		return ok;
	}
	
} // namespace OpenMS
