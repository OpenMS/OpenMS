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
// $Maintainer: Johannes Junker $
// $Authors: Johannes Junker $
// --------------------------------------------------------------------------

#include <OpenMS/config.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/VISUAL/TOPPASScene.h>
#include <OpenMS/SYSTEM/File.h>

#include <QtGui/QApplication>
#include <QtCore/QDir>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page TOPP_ExecutePipeline ExecutePipeline

	@brief Executes workflows created by TOPPAS.

	<B>The command line parameters of this tool are:</B>
	@verbinclude TOPP_ExecutePipeline.cli
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPExecutePipeline
	: public TOPPBase
{
 public:
	TOPPExecutePipeline()
		: TOPPBase("ExecutePipeline",
							 "Executes workflows created by TOPPAS.")
	{
	}

 protected:

	void registerOptionsAndFlags_()
	{
		registerInputFile_("in", "<file>", "", "The workflow to be executed (valid formats: \"toppas\")");
		registerStringOption_ ("out_dir", "<directory>", "", "The directory where the output files will be written", false);
	}

	ExitCodes main_(int argc, const char** argv)
	{
		QString toppas_file = getStringOption_("in").toQString();
		QString out_dir_name = getStringOption_("out_dir").toQString();
		
		QApplication a(argc, const_cast<char**>(argv), false);
		TOPPASScene ts(0, QDir::tempPath()+QDir::separator(), false);
		a.connect (&ts, SIGNAL(entirePipelineFinished()), &a, SLOT(quit()));
		a.connect (&ts, SIGNAL(pipelineExecutionFailed()), &a, SLOT(quit()));
		ts.load(toppas_file);
		
		if (out_dir_name != "")
		{
			if (QDir::isRelativePath(out_dir_name))
			{
				out_dir_name = QDir::currentPath() + QDir::separator() + out_dir_name;
			}
			
			if (File::exists(out_dir_name) && File::isDirectory(out_dir_name))
			{
				ts.setOutDir(out_dir_name);
			}
			else
			{
				cout << "The specified output directory does not exist." << endl;
				return CANNOT_WRITE_OUTPUT_FILE;
			}
		}
		else
		{
			cout << "No output directory specified. Using current directory..." << endl;
			
			if (!File::writable("test_file_in_the_current_directory"))
			{
				cout << "You do not have permission to write in the current directory." << endl;
				return CANNOT_WRITE_OUTPUT_FILE;
			}
		}
		
		ts.runPipeline();
		
		if (a.exec() == 0)
		{
			return EXECUTION_OK;
		}
		
		return UNKNOWN_ERROR;
	}

};


int main( int argc, const char** argv )
{
	TOPPExecutePipeline tool;
	return tool.main(argc, argv);
}

/// @endcond
