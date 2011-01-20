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
// $Maintainer: Nico Pfeifer $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

#include <OpenMS/APPLICATIONS/TOPPBase.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page TOPP_IDRTCalibration IDRTCalibration

	@brief Can be used to calibrate the RTs of peptide hits linearly to standards.

<CENTER>
	<table>
		<tr>
			<td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
			<td VALIGN="middle" ROWSPAN=2> \f$ \longrightarrow \f$ IDRTCalibration \f$ \longrightarrow \f$</td>
			<td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
		</tr>
		<tr>
			<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_MascotAdapter (or other ID engines) </td>
			<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_PeptideIndexer (or other tools operating @n with identifications (in idXML format))</td>
		</tr>
	</table>
</CENTER>

	This tool can be used to linearly align RTs of the IdXML-File to a reference. If only calibrant_1_input and
  calibrant_2_input are given, the first calibrant will result at RT 0.1 and calibrant_2_input will be at 0.9.
	If one wants to align the RTs of this IdXML file to the IDs of a reference file one can also give the RTs
	of the same calibrant in the reference file (calibrant_1_reference, calibrant_2_reference). If these calibrants
	are given, the linear transformation (shift and scale) will be calculated such that calibrant_1_input will
	be at the same RT as calibrant_1_reference and calibrant_2_input will
	be at the same RT as calibrant_2_reference. This only applies if calibrant_1* has a smaller RT than calibrant_2*.
	Otherwise the values are swapped.

	<B>The command line parameters of this tool are:</B>
	@verbinclude TOPP_IDRTCalibration.cli
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPIDRTCalibration
	: public TOPPBase
{
 public:
	TOPPIDRTCalibration()
		: TOPPBase("IDRTCalibration","Can be used to calibrate RTs of peptide hits linearly to standards.")
	{

	}

 protected:
	void registerOptionsAndFlags_()
	{
		registerInputFile_("in","<file>","", "input file ");
		setValidFormats_("in",StringList::create("idXML"));
		registerOutputFile_("out","<file>","","output file ");
		setValidFormats_("out",StringList::create("idXML"));
		registerDoubleOption_("calibrant_1_reference","<RT>",0.1,"The RT of the first calibrant in the reference file", false);
		registerDoubleOption_("calibrant_2_reference","<RT>",0.9,"The RT of the second calibrant in the reference file", false);
		registerDoubleOption_("calibrant_1_input","<RT>",std::numeric_limits<double>::quiet_NaN(),"The RT of the first calibrant in the input file");
		registerDoubleOption_("calibrant_2_input","<RT>",std::numeric_limits<double>::quiet_NaN(),"The RT of the second calibrant in the input file");
	}

	ExitCodes main_(int , const char**)
	{
		//-------------------------------------------------------------
		// parameter handling
		//-------------------------------------------------------------

		String in_file = getStringOption_("in");
		String out_file = getStringOption_("out");

		DoubleReal rt_calibrant_1_input = getDoubleOption_("calibrant_1_input");
		DoubleReal rt_calibrant_2_input =  getDoubleOption_("calibrant_2_input");
		DoubleReal rt_calibrant_1_reference =  getDoubleOption_("calibrant_1_reference");
		DoubleReal rt_calibrant_2_reference =  getDoubleOption_("calibrant_2_reference");

		if (rt_calibrant_1_input == rt_calibrant_2_input)
		{
			cout << "rt_calibrant_1_input and rt_calibrant_2_input must not have the same value";
			return ILLEGAL_PARAMETERS;
		}
		if (rt_calibrant_1_reference == rt_calibrant_2_reference)
		{
			cout << "rt_calibrant_1_reference and rt_calibrant_2_reference must not have the same value";
			return ILLEGAL_PARAMETERS;
		}

		//-------------------------------------------------------------
		// testing whether input and output files are accessible
		//-------------------------------------------------------------

		if (rt_calibrant_1_input > rt_calibrant_2_input)
		{
			DoubleReal temp = rt_calibrant_1_input;
			rt_calibrant_1_input = rt_calibrant_2_input;
			rt_calibrant_2_input = temp;
		}
		if (rt_calibrant_1_reference > rt_calibrant_2_reference)
		{
			DoubleReal temp = rt_calibrant_1_reference;
			rt_calibrant_1_reference = rt_calibrant_2_reference;
			rt_calibrant_2_reference = temp;
		}

		//-------------------------------------------------------------
		// calculations
		//-------------------------------------------------------------
		IdXMLFile file;
		vector<ProteinIdentification> protein_identifications;
		vector<PeptideIdentification> identifications;
		String document_id;
		file.load(in_file, protein_identifications, identifications, document_id);

		for (Size i = 0; i < identifications.size(); ++i)
		{
			if (identifications[i].metaValueExists("RT"))
			{
				DoubleReal temp_rt = identifications[i].getMetaValue("RT");
				temp_rt = (temp_rt - rt_calibrant_1_input) / (rt_calibrant_2_input - rt_calibrant_1_input)
					* (rt_calibrant_2_reference - rt_calibrant_1_reference) + rt_calibrant_1_reference;
				identifications[i].setMetaValue("RT", temp_rt);
			}
		}

		//-------------------------------------------------------------
		// writing output
		//-------------------------------------------------------------

		file.store(out_file,
							protein_identifications,
							identifications);

		return EXECUTION_OK;
	}
};


int main( int argc, const char** argv )
{
	TOPPIDRTCalibration tool;
	return tool.main(argc,argv);
}

/// @endcond
