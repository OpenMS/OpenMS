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
// $Authors: Marc Sturm, Mathias Walzer, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <QtCore/QProcess>

#include <fstream>

#include <boost/algorithm/string.hpp>

using namespace std;
using namespace OpenMS;


bool generate(const map<String,StringList>& tools, const String& prefix)
{
	bool errors_occured = false;
	for (map<String,StringList>::const_iterator it=tools.begin(); it!=tools.end(); ++it)
	{
		//start process
		QProcess process;
		process.setProcessChannelMode(QProcess::MergedChannels);
    QStringList env = QProcess::systemEnvironment();
    env << String("COLUMNS=110").toQString(); // Add an environment variable (used by each TOPP tool to determine width of help text (see TOPPBase))
    process.setEnvironment(env);
 		process.start((it->first + " --help").toQString());
		process.waitForFinished();

		ofstream f((String("output/")+ prefix + it->first + ".cli").c_str());
		std::string lines = QString(process.readAll()).toStdString();
		if(process.error() != QProcess::UnknownError)
		{
			// error while generation cli docu
			f << "Errors occured while generating the command line documentation for " << it->first << "!" << endl;
			f << "Please check your PATH variable if it contains the path to the " << it->first << " executable." << endl;
      f << "Output was: \n" << lines << endl;
			errors_occured = true;
		}
		else
		{
			// write output
			f << lines;
		}
		f.close();
	}
  return errors_occured;
}

int main (int , char** )
{
	//TOPP tools
	map<String,StringList> topp_tools = TOPPBase::getToolList();
	topp_tools["TOPPView"] = StringList();
	topp_tools["TOPPAS"] = StringList();
	//UTILS
	map<String,StringList> util_tools = TOPPBase::getUtilList();

  bool errors_occured = generate(topp_tools,"TOPP_") || generate(util_tools, "UTILS_");

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

