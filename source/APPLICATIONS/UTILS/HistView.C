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
// $Maintainer: $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

/**
	@page UTILS_HistView HistView
  
	@brief A viewer for histograms.
  	
	<B>The command line parameters of this tool are:</B>
	@verbinclude UTILS_HistView.cli
*/

//QT
#include <QApplication>
#include <QMainWindow>
#include <QStyleFactory>

//OpenMS
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/VISUAL/HistogramWidget.h>
#include <OpenMS/MATH/STATISTICS/Histogram.h>

//STL
#include <iostream>
#include <fstream>
#include <string>

using namespace OpenMS;
using namespace OpenMS::Math;
using namespace std;

void print_usage()
{
	cerr << endl
       << "HistView -- A viewer for histograms." << endl
       << endl
       << "Usage:" << endl
			 << " HistView <input> [options]" << endl
			 << endl
			 << "Options are:" << endl
			 << "  -bins <int>   The number of bins (default: 100)" << endl
			 << "  -min <float>  Start of value range (default: data minimum)" << endl
			 << "  -max <float>  End of value range (default: data maximum)" << endl
			 << "  -v            Prints verbose information to the command line" << endl
			 << "  --help        Shows this help" << endl
			 << endl
			 << "Note: <input> must contain one number per line!" << endl
			 << endl;
}

int main( int argc, const char** argv )
{
	//list of all the valid options
	Map<String,String> options, flags, option_lists;
	flags["--help"] = "help";
	flags["-v"] = "v";
	options["-bins"] = "bins";
	options["-min"] = "min";
	options["-max"] = "max";
	
	Param param;
	param.parseCommandLine(argc, argv, options, flags, option_lists);
	
	// '--help' given
	if (param.exists("help"))
	{
		print_usage();
		return 0;
	}	

	// test if unknown options were given
	if (param.exists("unknown"))
	{
		cout << "Unknown option(s) '" << param.getValue("unknown").toString() << "' given. Aborting!" << endl;
		print_usage();
		return 1;
	}

	// test if input file was given
	if (!param.exists("misc"))
	{
		cout << "No input file given. Aborting!" << endl;
		print_usage();
		return 1;
	}

  //set plastique style unless windows / mac style is available
  QApplication a( argc, const_cast<char**>(argv));
  if (QStyleFactory::keys().contains("windowsxp",Qt::CaseInsensitive))
  {
		a.setStyle("windowsxp");
  }
  else if (QStyleFactory::keys().contains("macintosh",Qt::CaseInsensitive))
  {
		a.setStyle("macintosh");
  }
  else if (QStyleFactory::keys().contains("plastique",Qt::CaseInsensitive))
  {
		a.setStyle("plastique");
  }
	
	bool verbose = false;
	if (param.exists("v")) verbose = true;
	
	//load input data
	if (verbose) cout << "Loading input data" << endl;
	vector<DoubleReal> input_data;
	StringList filenames = param.getValue("misc");
	ifstream is(filenames[0].c_str());
  if (!is)
  {
    cerr << "File '" << filenames[0].c_str() << "' not found!" << endl;
    return 1;
  }
  String str;
  while(getline(is,str,'\n'))
  {
  	try
  	{
  		DoubleReal value = str.toDouble();
  		input_data.push_back(value);
  	}
  	catch(Exception::ConversionError& /*e*/)
  	{
  		cerr << "Invalid input data line '" << str << "' is ignored!" << endl;
  	}
  }		
	
	//determine min and max
	if (verbose) cout << "Determining data minimum and maximum:" << endl;
	DoubleReal min = input_data[0];
	DoubleReal max = input_data[0];
	DoubleReal avg = 0.0;
	for (Size i=0; i<input_data.size(); ++i)
	{
		if(input_data[i]>max) max = input_data[i];
		if(input_data[i]<min) min = input_data[i];
		avg += input_data[i];
	}
	//overwrite by command line arguments
	if (param.exists("min")) min = param.getValue("min").toString().toDouble();
	if (param.exists("max")) max = param.getValue("max").toString().toDouble();
	if (verbose) cout << " - minimum: " << min << endl;
	if (verbose) cout << " - maximum: " << max << endl;
	if (verbose) cout << " - average: " << avg/input_data.size() << endl;
		
	//determine number of bin size
	if (verbose) cout << "Bins:" << endl;
	DoubleReal bins = 100.0;
	if (param.exists("bins")) bins = param.getValue("bins").toString().toDouble();
	DoubleReal bin_size = (max-min)/bins;
	if (verbose) cout << " - bins: " << bins << endl;
	if (verbose) cout << " - size: " << bin_size << endl;
	
	//create histogram	
	if (verbose) cout << "Creating histogram:" << endl;
	Histogram<> hist(min,max,bin_size);
	for (Size i=0; i<input_data.size(); ++i)
	{
		if (input_data[i]>=min && input_data[i]<=max)
		{
			hist.inc(input_data[i]);
		}
	}
	input_data.clear();
	if (verbose) cout << hist << endl;
	
  HistogramWidget* mw = new HistogramWidget(hist);
  mw->show();
  
  a.connect( &a, SIGNAL(lastWindowClosed()), &a, SLOT(quit()) );

  int result = a.exec();
  delete(mw);
  
  return result;
}

