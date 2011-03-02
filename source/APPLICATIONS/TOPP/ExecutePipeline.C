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
// $Maintainer: Johannes Junker $
// $Authors: Johannes Junker $
// --------------------------------------------------------------------------

#include <OpenMS/config.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/VISUAL/TOPPASScene.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/VISUAL/TOPPASResources.h>

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
 
  This tool is the non-GUI, i.e. command line version for non-interactive execution of TOPPAS pipelines.
  In order to really use this tool in batch-mode, you can provide a TOPPAS resource file (.trf) which specifies the
  input files for the input nodes in your pipeline.

  <B> *.trf files </B>

 A TOPPAS resource file (<TT>*.trf</TT>) specifies the locations of input files for a pipeline. 
 It is an XML file following the normal TOPP INI file schema, i.e. it can be edited using the INIFileEditor or filled using a script (we do NOT provide one - sorry).
 It can be exported from TOPPAS (<TT>File -> Save TOPPAS resource file</TT>). For two input nodes 1 and 2 with files (<TT>dataA.mzML</TT>, <TT>dataB.mzML</TT>) and (<TT>dataC.mzML</TT>) respectively it has the following format. 
 
  \code
<?xml version="1.0" encoding="ISO-8859-1"?>
<PARAMETERS version="1.3" xsi:noNamespaceSchemaLocation="http://open-ms.sourceforge.net/schemas/Param_1_3.xsd" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
  <NODE name="1" description="">
    <ITEMLIST name="url_list" type="string" description="">
      <LISTITEM value="file:///Users/jeff/dataA.mzML"/>
      <LISTITEM value="file:///Users/jeff/dataB.mzML"/>
    </ITEMLIST>
  </NODE>
  <NODE name="2" description="">
    <ITEMLIST name="url_list" type="string" description="">
      <LISTITEM value="file:///Users/jeff/dataC.mzML"/>
    </ITEMLIST>
  </NODE>
</PARAMETERS> 
  \endcode

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
    registerStringOption_ ("out_dir", "<directory>", "", "Directory for output files (default: user's home directory)", false);
		registerStringOption_ ("resource_file", "<file>", "", "A TOPPAS resource file (*.trf) specifying the files this workflow is to be applied to", false);
	}

	ExitCodes main_(int argc, const char** argv)
	{
		QString toppas_file = getStringOption_("in").toQString();
		QString out_dir_name = getStringOption_("out_dir").toQString();
		QString resource_file = getStringOption_("resource_file").toQString();

		QApplication a(argc, const_cast<char**>(argv), false);
		TOPPASScene ts(0, QDir::tempPath()+QDir::separator(), false);
		a.connect (&ts, SIGNAL(entirePipelineFinished()), &a, SLOT(quit()));
		a.connect (&ts, SIGNAL(pipelineExecutionFailed()), &a, SLOT(quit()));
		ts.load(toppas_file);

		if (resource_file != "")
		{
			TOPPASResources resources;
			resources.load(resource_file);
			ts.loadResources(resources);
		}

		if (out_dir_name != "")
		{
			if (QDir::isRelativePath(out_dir_name))
			{
				out_dir_name = QDir::currentPath() + QDir::separator() + out_dir_name;
			}
      out_dir_name = QDir::cleanPath(out_dir_name);
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
      QFileInfo fi(ts.getSaveFileName().toQString());
      out_dir_name = QDir::cleanPath( ts.getOutDir() + QDir::separator() + String(fi.baseName()).toQString() + QDir::separator() );
			cout << "No output directory specified. Using the user's home directory (" << out_dir_name.toStdString() << ")" << endl;
      ts.setOutDir(out_dir_name);
      QDir qd;
      if (!(qd.exists(out_dir_name) || qd.mkdir(out_dir_name)) || !File::writable(out_dir_name + "test_file_in_the_current_directory"))
			{
				cerr << "You do not have permission to write to " << out_dir_name.toStdString() << endl;
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
