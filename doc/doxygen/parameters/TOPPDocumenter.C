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
// $Maintainer: $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <QtCore/QProcess>

#include <fstream>

using namespace std;
using namespace OpenMS;

int main (int , char** )
{	
	//TOPP tools
	map<String,StringList> topp_tools = TOPPBase::getToolList();
	topp_tools["TOPPView"] = StringList();
	topp_tools["TOPPAS"] = StringList();
  bool errors_occured = false;
	for (map<String,StringList>::const_iterator it=topp_tools.begin(); it!=topp_tools.end(); ++it)
	{
		//start process
		QProcess process;
		process.setProcessChannelMode(QProcess::MergedChannels);
		process.start((it->first + " --help").toQString());
		process.waitForFinished();

    ofstream f((String("output/TOPP_") + it->first + ".cli").c_str());
    if(process.error() != QProcess::UnknownError)
    {
      // error while generation cli docu
      f << "Errors occured while generating the command line documentation for " << it->first << "!" << endl;
      f << "Please check your PATH variable if it contains the path to the " << it->first << " executable." << endl;
      errors_occured |= true;
    }
    else
    {
      // write output
      f << QString(process.readAllStandardOutput()).toStdString();
    }
    f.close();
	}

	//UTILS
	map<String,StringList> util_tools = TOPPBase::getUtilList();
	for (map<String,StringList>::const_iterator it=util_tools.begin(); it!=util_tools.end(); ++it)
	{
		//start process
		QProcess process;
		process.setProcessChannelMode(QProcess::MergedChannels);
		process.start((it->first + " --help").toQString());
		process.waitForFinished();
		//write output
		ofstream f((String("output/UTILS_") + it->first + ".cli").c_str());
    if(process.error() != QProcess::UnknownError)
    {
      // error while generation cli docu
      f << "Errors occured while generating the command line documentation for " << it->first << endl;
      f << "Please check your PATH variable if it contains the path to the " << it->first << " executable." << endl;
      errors_occured |= true;
    }
    else
    {
      // write output
      f << QString(process.readAllStandardOutput()).toStdString();
    }
    f.close();
  }


  if(errors_occured)
  {
    // errors occured while generating the TOPP CLI docu .. tell the user
    cerr << "Errors occured while generating the command line documentation for some of the " << endl;
    cerr << "TOPP tools/UTILS. Please check your PATH variable if it contains the TOPP tool directory." << endl;
    return EXIT_SUCCESS;
  }
  else
  {
    return EXIT_FAILURE;
  }
}

