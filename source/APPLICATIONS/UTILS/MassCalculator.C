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
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/FORMAT/SVOutStream.h>
#include <OpenMS/SYSTEM/File.h>
#include <iostream>
#include <ostream>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page UTILS_MassCalculator MassCalculator
	
	@brief Calculates masses and mass-to-charge ratios of peptide sequences.

	Given a peptide sequence and a charge state, the charged mass (including H+ adducts) and the mass-to-charge ratio are computed. The peptide sequence can include modifications (for information on valid notation see the @ref OpenMS::AASequence "AASequence" class documentation). Neutral masses can be computed by using "0" as charge state.

	Input can be given directly as values of the parameters: @p in for peptide sequences and @p charge for charge states. Alternatively, it can be read from a file with the following format: A peptide sequence at the beginning of each line, optionally followed by any number of charge states. Whitespace, commas or semicolons can de used to delimit the different items. Parts of the input that cannot be understood will be skipped. If charge states are given in the input file as well as via the @p charge parameter, results are returned for the union of both sets of charge states.

	Output can be written to a file or to the screen (see parameter @p out). Results for different charge states are always ordered from lowest to highest charge.
A number of different output formats are available via the parameter @p format:
	- @p list writes a human-readable list of the form "ABCDEF: z=1 m=566.192 m/z=566.192, z=2 m=567.199 m/z=283.599";
	- @p table produces a CSV-like table (using parameter @p separator to delimit fields) with the columns "peptide", "charge", "mass", and "mass-to-charge", and with one row per peptide and charge state;
	- @p mass_only writes only mass values (one line per peptide, values for different charge states separated by spaces);
	- @p mz_only writes only mass-to-charge ratios (one line per peptide, values for different charge states separated by spaces).
	

	<B>The command line parameters of this tool are:</B>
	@verbinclude UTILS_MassCalculator.cli
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPMassCalculator
	: public TOPPBase
{
public:

	TOPPMassCalculator()
		: TOPPBase("MassCalculator", "Calculates masses and mass-to-charge ratios of peptide sequences", false), use_avg_mass_(false), output_(0), format_(), res_type_(Residue::Full)
		{
			for (Size i = 0; i < Residue::SizeOfResidueType; i++)
			{
				Residue::ResidueType res_type = Residue::ResidueType(i);
				res_type_names_[Residue::getResidueTypeName(res_type)] = res_type;
			}
		}


protected:

	bool use_avg_mass_;
	ostream* output_; // pointer to output stream (stdout or file)
	String format_, separator_;
	Residue::ResidueType res_type_;
	map<String, Residue::ResidueType> res_type_names_;

	void registerOptionsAndFlags_()
		{
			registerStringList_("in", "<peptides/file>", StringList(), "List of peptide sequences, or single input file containing peptide sequences (and potentially charge numbers)");
			registerOutputFile_("out", "<file>", "", "Output file; if empty, output is written to the screen", false);
			registerIntList_("charge", "<numbers>", IntList(), "List of charge states; required if 'in' is a list of peptide sequences", false);
			registerStringOption_("format", "<choice>", "list", "Output format ('list': human-readable list, 'table': CSV-like table, 'mass_only': mass values only, 'mz_only': m/z values only)\n", false);
			setValidStrings_("format", StringList::create("list,table,mass_only,mz_only"));
			registerFlag_("average_mass", "Compute average (instead of monoisotopic) peptide masses");
			registerStringOption_("fragment_type", "<choice>", "full", "For what type of sequence/fragment the mass should be computed\n", false);
			setValidStrings_("fragment_type", StringList::create("full,internal,N-terminal,C-terminal,a-ion,b-ion,c-ion,x-ion,y-ion,z-ion"));
			registerStringOption_("separator", "<sep>", "", "Field separator for 'table' output format; by default, the 'tab' character is used", false);
		}
	
	DoubleReal computeMass_(const AASequence& seq, Int charge) const
		{
			if (use_avg_mass_) return seq.getAverageWeight(res_type_, charge);
			else return seq.getMonoWeight(res_type_, charge);
		}

	void writeTable_(const AASequence& seq, const set<Int>& charges)
		{
			SVOutStream sv_out(*output_, separator_);
			for (set<Int>::const_iterator it = charges.begin(); it != charges.end(); 
					 ++it)
			{
				DoubleReal mass = computeMass_(seq, *it);
				sv_out << seq.toString() << *it << mass;
				sv_out.writeValueOrNan(mass / *it);
				sv_out << endl;
			}
		}

	void writeList_(const AASequence& seq, const set<Int>& charges)
		{
			*output_ << seq.toString() << ": ";
			for (set<Int>::const_iterator it = charges.begin(); it != charges.end();
					 ++it)
			{
				DoubleReal mass = computeMass_(seq, *it);
				if (it != charges.begin()) *output_ << ", ";
				*output_ << "z=" << *it << " m=" << mass << " m/z=";
				if (*it != 0) *output_ << (mass / *it);
				else *output_ << "inf";
			}
			*output_ << endl;
		}

	void writeMassOnly_(const AASequence& seq, const set<Int>& charges, 
											bool mz=false)
		{
			for (set<Int>::const_iterator it = charges.begin(); it != charges.end(); 
					 ++it)
			{
				DoubleReal mass = computeMass_(seq, *it);
				if (it != charges.begin()) *output_ << " ";
				if (!mz) *output_ << mass;
				else if (*it == 0) *output_ << "inf";
				else *output_ << mass / *it;
			}
			*output_ << endl;
		}

	void writeLine_(const AASequence& seq, const set<Int>& charges)
		{
			if (format_ == "list") writeList_(seq, charges);
			else if (format_ == "table") writeTable_(seq, charges);
			else if (format_ == "mass_only") writeMassOnly_(seq, charges);
			else writeMassOnly_(seq, charges, true); // "mz_only"
		}

	String getItem_(String& line, const String& skip=" \t,;")
		{
			Size pos = line.find_first_of(skip);
			String prefix = line.substr(0, pos);
			pos = line.find_first_not_of(skip, pos);
			if (pos == String::npos) line = "";
			else line = line.substr(pos);
			return prefix;
		}

	void readFile_(const String& filename, const set<Int>& charges)
		{
			ifstream input(filename.c_str());
			String line;
			while (getline(input, line))
			{
				String item = getItem_(line);
				if ((item[0] == '"') && (item[item.size() - 1] == '"'))
				{
					item.unquote();
				}
				AASequence seq(item);
				if (!seq.isValid())
				{
					LOG_ERROR << "Error: '" << item
										<< "' is not a valid peptide sequence - skipping";
					continue;
				}
				set<Int> local_charges(charges);
				while (!line.empty())
				{
					item = getItem_(line);
					try
					{
						local_charges.insert(item.toInt());
					}
					catch (Exception::ConversionError) {};
				}
				if (local_charges.empty())
				{
					LOG_ERROR << "Error: No charge state specified - skipping";
					continue;
				}
				writeLine_(seq, local_charges);
			}
			input.close();
		}
	

	ExitCodes main_(int, const char**)
		{
			StringList in = getStringList_("in");
			String out = getStringOption_("out");
			IntList charge_list = getIntList_("charge");
			set<Int> charges(charge_list.begin(), charge_list.end());
			use_avg_mass_ = getFlag_("average_mass");
			res_type_ = res_type_names_[getStringOption_("fragment_type")];

			ofstream outfile;
			if (out.empty())
			{
				output_ = &cout;
			}
			else
			{
				outputFileWritable_(out, "out");
				outfile.open(out.c_str());
				output_ = &outfile;
			}

			format_ = getStringOption_("format");
			if (format_ == "table")
			{
				separator_ = getStringOption_("separator");
				if (separator_.empty()) separator_ = "\t";
				// write header:
				SVOutStream sv_out(*output_, separator_);
				sv_out << "peptide" << "charge" << "mass" << "mass-to-charge" << endl;
			}

			if ((in.size() == 1) && File::exists(in[0]))
			{
				inputFileReadable_(in[0], "in");
				readFile_(in[0], charges);
			}
			else
			{
				if (charges.empty())
				{
					LOG_ERROR << "Error: No charge state specified";
					return ILLEGAL_PARAMETERS;
				}
				for (StringList::iterator it = in.begin(); it != in.end(); ++it)
				{
					AASequence seq(*it);
					if (!seq.isValid())
					{
						LOG_ERROR << "Error: '" << *it
											<< "' is not a valid peptide sequence - skipping";
						continue;
					}
					writeLine_(seq, charges);
				}
			}
			
			if (!out.empty()) outfile.close();

			return EXECUTION_OK;
		}
};


int main( int argc, const char** argv )
{
	TOPPMassCalculator tool;
	return tool.main(argc, argv);
}
  
/// @endcond
