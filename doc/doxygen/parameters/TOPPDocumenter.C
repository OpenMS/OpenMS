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
// $Maintainer: Mathias Walzer $
// $Authors: Marc Sturm, Mathias Walzer $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <QtCore/QProcess>

#include <fstream>

#include <boost/algorithm/string.hpp>

using namespace std;
using namespace OpenMS;

int main (int , char** )
{
	size_t screen_length = 110;

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
			errors_occured = true;
		}
		else
		{
			// write output
			std::string lines = QString(process.readAllStandardOutput()).toStdString();
			std::vector<std::string> splits;
			boost::split(splits, lines, boost::is_any_of("\n"));
			for(size_t t = 0; t < splits.size(); ++t)
			{
				if(splits[t].size()>screen_length)
				{
					size_t n=0;
					f << splits[t].substr(n,screen_length) << " ..." << "\n" ;
					n += screen_length;
					while(n<splits[t].size())
					{
						f << " ... " << splits[t].substr(n,screen_length) << "\n" ;
						n += screen_length;
					}
				}
				else
				{
					f << splits[t] << "\n" ;
				}
			}
			f.close();
		}
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
			errors_occured = true;
		}
		else
		{
			// write output
			std::string lines = QString(process.readAllStandardOutput()).toStdString();
			std::vector<std::string> splits;
			boost::split(splits, lines, boost::is_any_of("\n"));
			for(size_t t = 0; t < splits.size(); ++t)
			{
				if(splits[t].size()>screen_length)
				{
					size_t n=0;
					f << splits[t].substr(n,screen_length) << " ..." << "\n" ;
					n += screen_length;
					while(n<splits[t].size())
					{
						f << " ... " << splits[t].substr(n,screen_length) << "\n" ;
						n += screen_length;
					}
				}
				else
				{
					f << splits[t] << "\n" ;
				}
			}
			//~ f << QString(process.readAllStandardOutput()).toStdString();
		}
		f.close();
	}

	if(errors_occured)
	{
		// errors occured while generating the TOPP CLI docu .. tell the user
		cerr << "Errors occured while generating the command line documentation for some of the " << endl;
		cerr << "TOPP tools/UTILS. Please check your PATH variable if it contains the TOPP tool directory." << endl;
		return EXIT_FAILURE;
	}
	else
	{
		return EXIT_SUCCESS;
	}
}

