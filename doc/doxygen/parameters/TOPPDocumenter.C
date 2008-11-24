// -*- mode: C++; tab-width: 2; -*-
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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <QtCore/QProcess>

#include <fstream>

using namespace std;
using namespace OpenMS;

int main (int , char** )
{
	bool all_calls_worked = true;
	
	ofstream f;
	f.open("TOPPParameters.doxygen");
	
	StringList tools = TOPPBase::getToolList();
	for (UInt i=0; i<tools.size(); ++i)
	{
		String tool = tools[i];
		f << "/**" << endl;
		f << " @page TOPP_" << tool << "_CLI " << tool << " command line interface" << endl;
		f << "Command line interface of the TOPP tool @ref TOPP_" << tool << " :\n";
		
		QProcess process;
		process.setProcessChannelMode(QProcess::MergedChannels);
		process.start((tool + " --help").toQString());
		process.waitForFinished();
		String output = QString(process.readAllStandardOutput());
		//abort of call did not work
		if (process.exitStatus()!=QProcess::NormalExit)
		{
			all_calls_worked = false;
		}
		f << " @verbatim " << output << " @endverbatim" << endl;
		f << "*/" << endl;
		f << endl;
	}
	
	if (!all_calls_worked)
	{
		cerr << "TOPPDocumenter: Not all tools could be called successfully!" << endl;
		return 1;	
  }
  return 0;
}
