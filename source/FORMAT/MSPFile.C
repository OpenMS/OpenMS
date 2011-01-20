// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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

#include <OpenMS/FORMAT/MSPFile.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/KERNEL/RichPeak1D.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/SYSTEM/File.h>
#include <fstream>

using namespace std;

namespace OpenMS
{
	MSPFile::MSPFile()
		: DefaultParamHandler("MSPFile")
	{
		defaults_.setValue("parse_headers", "false", "Flag whether header information should be parsed an stored for each spectrum");
		vector<String> parse_strings;
		parse_strings.push_back("true");
		parse_strings.push_back("false");
		defaults_.setValidStrings("parse_headers", parse_strings);
		defaults_.setValue("parse_peakinfo", "true", "Flag whether the peak annotation information should be parsed and stored for each peak");
		defaults_.setValidStrings("parse_peakinfo", parse_strings);
		defaults_.setValue("instrument", "", "If instrument given, only spectra of these type of instrument (Inst= in header) are parsed");
		defaults_.setValidStrings("instrument", StringList::create(",it,qtof,toftof"));

		defaultsToParam_();
	}


	MSPFile::MSPFile(const MSPFile& rhs)
		: DefaultParamHandler(rhs)
	{
	}

	MSPFile& MSPFile::operator = (const MSPFile& rhs)
	{
		if (this != &rhs)
		{
			DefaultParamHandler::operator = (rhs);
		}
		return *this;
	}

	MSPFile::~MSPFile()
	{
	}


  void MSPFile::load(const String& filename, vector<PeptideIdentification>& ids, RichPeakMap& exp)
  {
    if (!File::exists(filename))
    {
    	throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
    }
    if (!File::readable(filename))
    {
      throw Exception::FileNotReadable(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
    }

		exp.reset();

		//set DocumentIdentifier
		exp.setLoadedFileType(filename);
		exp.setLoadedFilePath(filename);

    String line;
    ifstream is(filename.c_str());

    Map<String, DoubleReal> mod_to_mass;
    mod_to_mass["Oxidation"] = 15.994915;
    mod_to_mass["Carbamidomethyl"] = 57.02146;
    mod_to_mass["ICAT_light"] = 227.126991;
    mod_to_mass["ICAT_heavy"] = 236.157185;
    mod_to_mass["AB_old_ICATd0"] = 442.224991;
    mod_to_mass["AB_old_ICATd8"] = 450.275205;
    mod_to_mass["Acetyl"] = 42.0106;
    mod_to_mass["Deamidation"] = 0.9840;
    mod_to_mass["Pyro-cmC"] = -17.026549;
    mod_to_mass["Pyro-glu"] = -18.010565;
		mod_to_mass["Gln->pyro-Glu"] = -18.010565;
    mod_to_mass["Amide"] = -0.984016;
    mod_to_mass["Phospho"] = 79.9663;
    mod_to_mass["Methyl"] = 14.0157;
    mod_to_mass["Carbamyl"] = 43.00581;
		mod_to_mass["di-Methylation"] = 28.031300;

		Map<String, String> modname_to_unimod;
		modname_to_unimod["Pyro-glu"] = "Gln->pyro-Glu";
		modname_to_unimod["AB_old_ICATd8"] = "ICAT-D:2H(8)";
		modname_to_unimod["AB_old_ICATd0"] = "ICAT-D";

    RichPeakSpectrum spec;

		bool parse_headers(param_.getValue("parse_headers").toBool());
		bool parse_peakinfo(param_.getValue("parse_peakinfo").toBool());
		String instrument((String)param_.getValue("instrument"));
		bool inst_type_correct(true);
		Size spectrum_number = 0;

		// line number counter
		Size line_number = 0;
				
		while (getline(is, line))
    {
			++line_number;
			
      if (line.hasPrefix("Name:"))
      {
        vector<String> split, split2;
        line.split(' ', split);
        split[1].split('/', split2);
        String peptide = split2[0];
        // remove damn (O), also defined in 'Mods=' comment
        peptide.substitute("(O)", "");
        PeptideIdentification id;
        id.insertHit(PeptideHit(0, 0, split2[1].toInt(), peptide));
        ids.push_back(id);
				inst_type_correct = true;
      }
			else if (line.hasPrefix("MW:"))
      {
				vector<String> split;
				line.split(' ', split);
				if (split.size() == 2)
				{
					UInt charge = ids.back().getHits().begin()->getCharge();
					spec.getPrecursors().resize(1);
        	spec.getPrecursors()[0].setMZ((split[1].toDouble() + (DoubleReal)charge * Constants::PROTON_MASS_U) / (DoubleReal)charge);
				}
      }
			else if (line.hasPrefix("Comment:"))
      {
        // slow, but we need the modifications from the header and the instrument type
        vector<String> split;
        line.split(' ', split);
        for (vector<String>::const_iterator it = split.begin(); it != split.end(); ++it)
        {
					if (!inst_type_correct)
					{
						break;
					}
					if (instrument != "" && it->hasPrefix("Inst="))
					{
						String inst_type = it->suffix('=');
						if (instrument != inst_type)
						{
							inst_type_correct = false;
							ids.erase(--ids.end());
						}
						break;
					}

          if (it->hasPrefix("Mods=") && *it != "Mods=0")
          {
            String mods = it->suffix('=');
            // e.g. Mods=2/7,K,Carbamyl/22,K,Carbamyl
            vector<String> mod_split;
            mods.split('/', mod_split);
            AASequence peptide = ids.back().getHits().begin()->getSequence();
            for (Size i = 1; i <= (UInt)mod_split[0].toInt(); ++i)
            {
              vector<String> single_mod;
              mod_split[i].split(',', single_mod);

              String mod_name = single_mod[2];
							if (modname_to_unimod.has(mod_name))
							{
								mod_name = modname_to_unimod[mod_name];
							}
							String origin  = single_mod[1];
							Size position = single_mod[0].toInt();

							//cerr << "MSP modification: " << origin << " " << mod_name << " " << position << "\n";

							if (position > 0 && position < peptide.size() - 1)
							{
								peptide.setModification(position, mod_name);
							}
							else if (position == 0)
							{
								// we must decide whether this can be a terminal mod
								try 
								{
              		peptide.setNTerminalModification(mod_name);
								}
								catch (Exception::ElementNotFound /*e*/)
								{
									peptide.setModification(position, mod_name);
								}
							}
							else if (position == peptide.size() - 1)
							{
								// we must decide whether this can be a terminal mod
								try 
								{
									peptide.setCTerminalModification(mod_name);
								}
								catch (Exception::ElementNotFound /*e*/)
								{
									peptide.setModification(position, mod_name);
								}
							}
							else
							{
								cerr << "MSPFile: Error: ignoring modification: '" << line << "' in line " << line_number << "\n";
							}
            }
            vector<PeptideHit> hits(ids.back().getHits());
            hits.begin()->setSequence(peptide);
            ids.back().setHits(hits);
          }
        }

        if (parse_headers && inst_type_correct)
        {
          parseHeader_(line, spec);
        }
      }
			else if (line.hasPrefix("Num peaks:"))
			{
				if (!inst_type_correct)
				{
					while (getline(is, line) && ++line_number && line.size() > 0 && isdigit(line[0])) ;
				}
				else
				{
        	while (getline(is, line) && ++line_number && line.size() > 0 && isdigit(line[0]))
        	{
          	vector<String> split;
          	line.split('\t', split);
          	RichPeak1D peak;
          	if (split.size() != 3)
          	{
            	throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, line, "not <mz><tab><intensity><tab>\"<comment>\" in line " + String(line_number));
          	}
          	peak.setMZ(split[0].toFloat());
          	peak.setIntensity(split[1].toFloat());
						if (parse_peakinfo)
						{
           		peak.setMetaValue("MSPPeakInfo", split[2]);
						}
          	spec.push_back(peak);
        	}
					spec.setNativeID(String("index=")+spectrum_number);
       		exp.push_back(spec);
        	spec.clear(true);
      	}
      	spectrum_number++;
			}
    }

		// last spectrum, if available


  }




	void MSPFile::parseHeader_(const String& header, RichPeakSpectrum& spec)
  {
    // first header from std_protein of NIST spectra DB
    // Spec=Consensus Pep=Tryptic Fullname=R.AAANFFSASCVPCADQSSFPK.L/2 Mods=0 Parent=1074.480 Inst=it Mz_diff=0.500 Mz_exact=1074.4805 Mz_av=1075.204 Protein="TRFE_BOVIN" Organism="Protein Standard" Se=2^X23:ex=3.1e-008/1.934e-005,td=5.14e+007/2.552e+019,sd=0/0,hs=45.8/5.661,bs=6.3e-021,b2=1.2e-015,bd=5.87e+020^O22:ex=3.24e-005/0.0001075,td=304500/5.909e+297,pr=3.87e-007/1.42e-006,bs=1.65e-301,b2=1.25e-008,bd=1.3e+299 Sample=1/bovine-serotransferrin_cam,23,26 Nreps=23/34 Missing=0.3308/0.0425 Parent_med=1074.88/0.23 Max2med_orig=22.1/9.5 Dotfull=0.618/0.029 Dot_cons=0.728/0.040 Unassign_all=0.161 Unassigned=0.000 Dotbest=0.70 Naa=21 DUScorr=2.3/3.8/0.61 Dottheory=0.86 Pfin=4.3e+010 Probcorr=1 Tfratio=8e+008 Specqual=0.0
		vector<String> split;
    header.split(' ', split);

    for (vector<String>::const_iterator it = split.begin(); it != split.end(); ++it)
    {
    	vector<String> split2;
      String tmp = *it;
      tmp.trim();
      tmp.split('=', split2);
      if (split2.size() == 2)
      {
        spec.setMetaValue(split2[0], split2[1]);
      }
    }
  }

	void MSPFile::store(const String& filename, const RichPeakMap& exp) const
	{
		if (!File::writable(filename))
		{
			throw Exception::FileNotWritable(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
		}

		ofstream out(filename.c_str());

		for (RichPeakMap::ConstIterator it = exp.begin(); it != exp.end(); ++it)
		{
			if (it->getPeptideIdentifications().size() > 0 && it->getPeptideIdentifications().begin()->getHits().size() > 0)
			{
				PeptideHit hit = *it->getPeptideIdentifications().begin()->getHits().begin();
				String peptide;
				for (AASequence::ConstIterator pit = hit.getSequence().begin(); pit != hit.getSequence().end(); ++pit)
				{
					if (pit->isModified() && pit->getOneLetterCode() == "M" &&
							fabs(ModificationsDB::getInstance()->getModification("M", pit->getModification(), ResidueModification::ANYWHERE).getDiffFormula().getMonoWeight() - 16.0) < 0.01)
					{
						peptide += "M(O)";
					}
					else
					{
						peptide += pit->getOneLetterCode();
					}
				}
				out << "Name: " << peptide << "/" << hit.getCharge() << "\n";
				out << "MW: " << hit.getSequence().getMonoWeight() << "\n";
				out << "Comment:";
			
				// modifications
				// e.g. 2/9,C,Carbamidomethyl/12,C,Carbamidomethyl
				Size num_mods(0);
				vector<String> modifications;
				if (hit.getSequence().hasNTerminalModification())
				{
					String mod = hit.getSequence().getNTerminalModification();
					++num_mods;
					String modification = "0," + hit.getSequence().begin()->getOneLetterCode() + ",";
					modification += ModificationsDB::getInstance()->getTerminalModification(mod, ResidueModification::N_TERM).getId();
					modifications.push_back(modification);
				}
				
				// @improvement improve writing support (Andreas)
				UInt pos(0);
				for (AASequence::ConstIterator pit = hit.getSequence().begin(); pit != hit.getSequence().end(); ++pit, ++pos)
				{
					if (!pit->isModified())
					{
						continue;	
					}

					String mod = pit->getModification();
					String res = pit->getOneLetterCode();
					++num_mods;
					String modification = String(pos) + "," + res + ",";
					modification += ModificationsDB::getInstance()->getModification(res, mod, ResidueModification::ANYWHERE).getId();
					modifications.push_back(modification);
				}
			
				String mods;
				mods.concatenate(modifications.begin(), modifications.end(), "/");
				if (mods.size() != 0)
				{
					out << " Mods=" << String(num_mods)  << "/" << mods;
				}
				else
				{
					out << " Mods=0";
				}
				out << " Inst=it\n"; // @improvement write instrument type, protein...and other information
				out << "Num peaks: " << it->size() << "\n";

				// normalize to 10,000
				RichPeakSpectrum rich_spec = *it;
				DoubleReal max_int(0);
				for (RichPeakSpectrum::ConstIterator sit = rich_spec.begin(); sit != rich_spec.end(); ++sit)
				{
					if (sit->getIntensity() > max_int)
					{
						max_int = sit->getIntensity();
					}
				}

				if (max_int != 0)
				{
					for (RichPeakSpectrum::Iterator sit = rich_spec.begin(); sit != rich_spec.end(); ++sit)
					{
						sit->setIntensity(sit->getIntensity() / max_int * 10000.0);
					}
				}
				else
				{
					cerr << "MSPFile: spectrum contains only zero intensities!" << endl;
				}
				
				for (RichPeakSpectrum::ConstIterator sit = rich_spec.begin(); sit != rich_spec.end(); ++sit)
				{
					out << sit->getPosition()[0] << "\t" << sit->getIntensity() << "\t";
					if (sit->metaValueExists("IonName"))
					{
						out << "\"" << sit->getMetaValue("IonName") << "\"";
					}
					else
					{
						out << "\"?\"";
					}
					out << "\n";
				}
				out << "\n";
			}
		}
	}

}// namespace OpenMS

