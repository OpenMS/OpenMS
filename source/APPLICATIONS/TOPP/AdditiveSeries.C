// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Chris Bielow $
// $Authors: Clemens Groepl, Chris Bielow$
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/KERNEL/FeatureMap.h>

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/MATH/STATISTICS/LinearRegression.h>

#include <OpenMS/SYSTEM/File.h>

#include <map>
#include <vector>
#include <algorithm>
#include <gsl/gsl_math.h>

using namespace OpenMS;
using namespace Math;
using namespace std;

typedef Feature::CoordinateType CoordinateType;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page TOPP_AdditiveSeries AdditiveSeries

	@brief Computes an additive series to quantify a peptide in a set of samples.
<CENTER>
	<table>
		<tr>
			<td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
			<td VALIGN="middle" ROWSPAN=3> \f$ \longrightarrow \f$ AdditiveSeries \f$ \longrightarrow \f$</td>
			<td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
		</tr>
		<tr>
      <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_FeatureFinderCentroided </td>
			<td VALIGN="middle" ALIGN = "center" ROWSPAN=2> - </td>
		</tr>
		<tr>
			<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_IDMapper </td>
		</tr>
	</table>
</CENTER>
	This module computes an additive series for an absolute
	quantification of a peptide in a set of samples. The
	output consists of a GNUplot script which can be used
	to visualize the results and some XML output for further processing.

	In this version, the application computes the additive
	series as a ratio of the intensities of two different peptides.
	One of these peptides serves as internal standard for
	calibration.

	<B>The command line parameters of this tool are:</B>
	@verbinclude TOPP_AdditiveSeries.cli
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class AdditiveSeries
	: public TOPPBase
{
 public:
	AdditiveSeries()
		: TOPPBase("AdditiveSeries","Computes an additive series to quantify a peptide in a set of samples.")
	{
	}

 protected:
	void registerOptionsAndFlags_()
	{
		registerInputFileList_("in","<files>",StringList(),"input files separated by blanks",true);
		setValidFormats_("in",StringList::create("featureXML"));
		registerOutputFile_("out","<file>","","output XML file containg regression line and confidence interval");
		registerDoubleOption_("mz_tolerance","<tol>",1.0, "Tolerance in m/z dimension",false);
		registerDoubleOption_("rt_tolerance","<tol>",1.0, "Tolerance in RT dimension",false);
    registerDoubleList_("concentrations","<concentrations>", DoubleList(), "Spiked concentrations");
		addEmptyLine_();
		addText_("  Feature/standard position:");
		registerDoubleOption_("feature_rt","<rt>",std::numeric_limits<double>::quiet_NaN(), "RT position of the feature");
		registerDoubleOption_("feature_mz","<mz>",std::numeric_limits<double>::quiet_NaN(), "m/z position of the feature");
		registerDoubleOption_("standard_rt","<rt>",std::numeric_limits<double>::quiet_NaN(), "RT position of the standard");
		registerDoubleOption_("standard_mz","<mz>",std::numeric_limits<double>::quiet_NaN(), "m/z position of the standard");

		addEmptyLine_();
		addText_("  GNUplot options:");
		registerFlag_("write_gnuplot_output","Flag that activates the GNUplot output");
		registerStringOption_("out_gp","<name>","","base file name (3 files with different extensions are created)",false);
		registerStringOption_("mz_unit","<unit>","Thomson","the m/z unit of the plot",false);
		registerStringOption_("rt_unit","<unit>","seconds","the RT unit of the plot",false);
  }


  // searches for a features with coordinates within the tolerance in this map
	// NOTE: It might happen that there are several features at similar coordinates.
	// In this case, the program cannot be sure which one is the correct. So we decided
	// to use the one with the strongest intensity.
	bool readMapFile_(String filename, vector<double>& intensities,
										CoordinateType tol_mz, CoordinateType tol_rt,
										DPosition<2> fpos1, DPosition<2> fpos2)
	{

		if (!File::exists(filename) )
		{
			cout << "File " << filename << " not found. " << endl;
			return false;
		}

		cout << "Reading from " << filename << endl;

		FeatureXMLFile map_file;
		FeatureMap<> map;
		map_file.load(filename,map);

		Feature* feat1 = 0;
		Feature* feat2 = 0;

		FeatureMap<>::iterator iter = map.begin();
		while (iter!= map.end() )
		{

// 			cout << "Read: " << *iter << endl;

			if ( (iter->getRT() <  fpos1[Feature::RT] + tol_rt) &&
					 (iter->getRT() >  fpos1[Feature::RT] - tol_rt) &&
					 (iter->getMZ() <  fpos1[Feature::MZ] + tol_mz) &&
					 (iter->getMZ() >  fpos1[Feature::MZ] - tol_mz) )
			{
// 				cout << "Found feature1 at " << endl;
// 				cout << iter->getRT() << " " << iter->getMZ()  << " " << iter->getIntensity() <<  endl;
				// feature at correct position found, save intensity
				if (!feat1)
				{
					feat1 = &(*iter);
				}
				else if (feat1->getIntensity() <  iter->getIntensity() )
				{
					feat1 = &(*iter);
				}
				// 				f1_sum += iter->getIntensity();

			}

			if ( (iter->getRT() <  fpos2[Feature::RT] + tol_rt) &&
					 (iter->getRT() >  fpos2[Feature::RT] - tol_rt) &&
					 (iter->getMZ() <  fpos2[Feature::MZ] + tol_mz) &&
					 (iter->getMZ() >  fpos2[Feature::MZ] - tol_mz) )
			{
// 				cout << "Found feature2 at " << endl;
// 				cout << iter->getRT() << " " << iter->getMZ() << " " << iter->getIntensity() <<  endl;
				// same as above
				if (!feat2)
				{
					feat2 = &(*iter);
				}
				else if (feat2->getIntensity() <  iter->getIntensity() )
				{
					feat2 = &(*iter);
				}

				// 				f2_sum += iter->getIntensity();
			}

			++iter;
		}	// end of while

		if (feat1 != 0 && feat2 != 0)  //(f1_sum != 0 && f2_sum != 0)
		{
			cout << "Feature 1: " << *feat1 << endl;
			cout << "Feature 2: " << *feat2 << endl;
			cout << "Intensity ratio : " << ( feat1->getIntensity() / feat2->getIntensity() ) << endl;
			intensities.push_back( feat1->getIntensity() / feat2->getIntensity());

			return true;
		}
		if (!feat1)
			writeDebug_(String("Feature 1 was not found. "),1);

		if (!feat2)
			writeDebug_(String("Feature 2 was not found. "),1);

		return false;
	}

	/*
		 Computes the linear regression for a series of measurements, the
		 x-axis intercept of the regression line and its confidence interval, and
		 writes a couple of files from which a nice plot of all this can be
		 generated using the gnuplot program.
	*/
	bool computeRegressionAndWriteGnuplotFiles_ ( vector<double>::const_iterator const conc_vec_begin,
																								vector<double>::const_iterator const conc_vec_end,
																								vector<double>::const_iterator const area_vec_begin,
																								double const confidence_p,
																								String const filename_prefix,
																								String const output_filename,
																								String const format = "",
																								bool const write_gnuplot = true
																							)
	{

		try
		{
			LinearRegression linreg;
			linreg.computeRegression ( confidence_p, conc_vec_begin, conc_vec_end, area_vec_begin );

			if (write_gnuplot)
			{

				// the peak data goes here
				String datafilename(filename_prefix);
				datafilename+=String(".dat");
				ofstream dataout(datafilename.c_str());

				// the gnuplot commands go here
				String commandfilename(filename_prefix);
				commandfilename+=String(".cmd");
				ofstream cmdout(commandfilename.c_str());

				// the error bar for the x-axis intercept goes here
				String errorbarfilename(filename_prefix);
				errorbarfilename+=String(".err");
				ofstream errout(errorbarfilename.c_str());

				// writing the commands
				cmdout <<
					"set ylabel \"ion count\"\n"
					"set xlabel \"concentration\"\n"
					"set key left Left reverse\n"
					;

				if ( ! format.empty() )
				{
					if ( format == "png" )
					{
						cmdout <<
							"set terminal png \n"
							"set output \"" << filename_prefix << ".png\"\n" ;
					}
					else if ( format == "eps" )
					{
						cmdout <<
							"set terminal postscript eps \n"
							"set output \"" << filename_prefix << ".eps\"\n" ;
					}

				}
				cmdout <<
					"plot \""  << datafilename <<"\"  w points ps 2 pt 1 lt 8 title \"data\" " // want data on first line of key
					",  " << linreg.getIntercept() << "+" <<  linreg.getSlope() << "*x lt 2 lw 3 title \"linear regression: "
							 << linreg.getIntercept() << " + " <<  linreg.getSlope() << " * x\" "
					", \""  << datafilename <<"\"  w points ps 2 pt 1 lt 8 notitle " // draw data a second time, on top of reg. line
					", \"" << errorbarfilename << "\"  using ($1):(0) w points pt 13 ps 2 lt 1 title \"x-intercept: " << linreg.getXIntercept() << "\" "
					", \"" << errorbarfilename << "\"  w xerrorbars lw 3 lt 1 title \"95% interval: [ " << linreg.getLower() << ", " << linreg.getUpper() << " ]\"\n"
					;
				cmdout.close();

				// writing the x-axis intercept error bar
				errout << linreg.getXIntercept() << " 0 " << linreg.getLower() << " " << linreg.getUpper() << endl;
				errout.close();

				// writing the peak data points
				vector<double>::const_iterator cit = conc_vec_begin;
				vector<double>::const_iterator ait = area_vec_begin;
				dataout.precision(writtenDigits<DoubleReal>());
				for ( ;cit != conc_vec_end; ++cit, ++ait )
				{
					dataout << *cit << ' ' << *ait << '\n';
				}
				dataout.close();

			} // end if (write_gnuplot)

			// write results to XML file
			ofstream results;
			results.open(output_filename.c_str());

			results << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>" << endl;
			results << "<results_additiveseries>" << endl;
			results << "\t<slope>" << linreg.getSlope() << "</slope>" << endl;
			results << "\t<intercept>" << linreg.getIntercept() << "</intercept>" << endl;
			results << "\t<x_intercept>" << linreg.getXIntercept() << "</x_intercept>" << endl;
			results << "\t<confidence_lowerlimit>" << linreg.getLower() << "</confidence_lowerlimit>" << endl;
			results << "\t<confidence_upperlimit>" << linreg.getUpper() << "</confidence_upperlimit>" << endl;
			results << "\t<pearson_squared>" << linreg.getRSquared() << "</pearson_squared>" << endl;
			results << "\t<std_residuals>" << linreg.getStandDevRes() << "</std_residuals>" << endl;
			results << "\t<t_statistic>" << linreg.getTValue() << "</t_statistic>" << endl;
			results << "</results_additiveseries>" << endl;

			results.close();
		}
		catch (string& s)
		{
			cout << s <<  endl;
			return 1;
		}

		return 0;
	}

	ExitCodes main_(int , const char**)
	{
		//-------------------------------------------------------------
		// parsing parameters
		//-------------------------------------------------------------
		Param const& add_param =  getParam_();
		writeDebug_("Used parameters", add_param, 3);

		CoordinateType tol_mz = getDoubleOption_("mz_tolerance");
		CoordinateType tol_rt = getDoubleOption_("rt_tolerance");

		String out_f  = getStringOption_("out");

		if (add_param.getValue("feature_mz").isEmpty() || add_param.getValue("feature_rt").isEmpty() )
		{
			writeLog_("Feature coordinates not given. Aborting.");
			return ILLEGAL_PARAMETERS;
		}
		DPosition<2> feat_pos1;
		feat_pos1[Feature::MZ] = (CoordinateType) add_param.getValue("feature_mz");
		feat_pos1[Feature::RT] = (CoordinateType) add_param.getValue("feature_rt");

		if (add_param.getValue("standard_mz").isEmpty() || add_param.getValue("standard_rt").isEmpty() )
		{
			writeLog_("Standard coordinates not given. Aborting.");
			return ILLEGAL_PARAMETERS;
		}
		DPosition<2> feat_pos2;
		feat_pos2[Feature::MZ] = (CoordinateType) add_param.getValue("standard_mz");
		feat_pos2[Feature::RT] = (CoordinateType) add_param.getValue("standard_rt");

		writeDebug_(String("Setting tolerances to ") + tol_mz + " " + tol_rt,1);

		// introduce a flag for each concentration. true => the corresponding feature was found
		vector<bool> flags;

		// fetching list of files
		StringList files = getStringList_("in");

		// collect features
		vector<DoubleReal> intensities;
		vector<String>::const_iterator cit = files.begin();
		while (cit != files.end())
		{
			if (readMapFile_(*cit,intensities,tol_mz,tol_rt,feat_pos1,feat_pos2) )
			{
				flags.push_back(true);
			}
			else
			{
				flags.push_back(false);
			}
			++cit;
		}

		// read the spiked concentrations
    DoubleList sp_concentrations = getDoubleList_("concentrations");

		vector<double> sp_concentrations2;
		for (Size i=0; i<sp_concentrations.size(); i++)
		{
			if (flags.at(i) == true )
			{
				sp_concentrations2.push_back( sp_concentrations.at(i) );
			}
		}

		cout << "Found feature pairs: " <<  intensities.size() << endl;
		cout << "Spiked concentrations: " << sp_concentrations.size() << endl;

		if (intensities.empty() || sp_concentrations.empty() )
		{

			writeLog_("Did not find any data. Aborting!");
			return ILLEGAL_PARAMETERS;
		}

		// set prefix of gnuplot output
		String filename_prefix = getStringOption_("out_gp");;
		if (getFlag_("write_gnuplot_output"))
		{
			writeDebug_(String("Writing gnuplot output"),1);
			computeRegressionAndWriteGnuplotFiles_ (sp_concentrations2.begin(), sp_concentrations2.end(),
																							intensities.begin(), 0.95, filename_prefix, out_f, "eps", true);
		}
		else
		{
			writeDebug_(" No GNUplot output is written...",1);
			computeRegressionAndWriteGnuplotFiles_ (sp_concentrations2.begin(), sp_concentrations2.end(),
																							intensities.begin(), 0.95, filename_prefix, out_f, "eps", false);
		}

		return EXECUTION_OK;
	}

};


/// @endcond


int main( int argc, const char** argv )
{
    AdditiveSeries tool;
    return tool.main(argc,argv);
}


