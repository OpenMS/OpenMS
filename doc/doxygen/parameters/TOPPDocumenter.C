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
	for (map<String,StringList>::const_iterator it=topp_tools.begin(); it!=topp_tools.end(); ++it)
	{
		//start process
		QProcess process;
		process.setProcessChannelMode(QProcess::MergedChannels);
		process.start((it->first + " --help").toQString());
		process.waitForFinished();
		//write output
		ofstream f((String("output/TOPP_") + it->first + ".cli").c_str());
		f << QString(process.readAllStandardOutput()).toStdString();
	}
	
	//UTILS
	StringList utils_tools = StringList::create("IDMassAccuracy,DecoyDatabase,MapAlignmentEvaluation,CaapConvert,CVInspector,DecoyDatabase,Digestor,FFEval,FuzzyDiff,HistView,IDExtractor,LabeledEval,SemanticValidator,SequenceCoverageCalculator,XMLValidator,IdXMLEvaluation,MSSimulator,ERPairFinder,SpecLibCreator,SpectrumGeneratorNetworkTrainer,MRMPairFinder,DeMeanderize");
	for (Size i=0; i<utils_tools.size(); ++i)
	{
		//start process
		QProcess process;
		process.setProcessChannelMode(QProcess::MergedChannels);
		process.start((utils_tools[i] + " --help").toQString());
		process.waitForFinished();
		//write output
		ofstream f((String("output/UTILS_") + utils_tools[i] + ".cli").c_str());
		f << QString(process.readAllStandardOutput()).toStdString();
	}
	
  return 0;
}

