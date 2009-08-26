// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Andreas Bertsch $
// $Authors: Andreas Bertsch, Daniel Jameson$
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/MascotXMLFile.h>
#include <OpenMS/FORMAT/MascotRemoteQuery.h>
#include <OpenMS/FORMAT/MascotGenericFile.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <sstream>

#include <QtCore/QDir>
#include <QtCore/QFile>
#include <QtCore/QCoreApplication>
#include <QtCore/QTimer>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page TOPP_MascotAdapterOnline MascotAdapterOnline

	@brief Identifies peptides in MS/MS spectra via Mascot.

	@experimental This tool has not been tested thoroughly and might behave not as expected!

	This wrapper application serves for getting peptide identifications
	for MS/MS spectra.

	<B>The command line parameters of this tool are:</B>
	@verbinclude TOPP_MascotAdapterOnline.cli
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES


class TOPPMascotAdapterOnline
	: public TOPPBase
{
	public:
		TOPPMascotAdapterOnline()
			: TOPPBase("MascotAdapterOnline","Annotates MS/MS spectra using Mascot.")
		{
		}

	protected:

		void registerOptionsAndFlags_()
		{
			registerInputFile_("in", "<file>", "", "input file in mzML format.\n");
			setValidFormats_("in", StringList::create("mzML"));
			registerOutputFile_("out", "<file>", "", "output file in IdXML format.\n");
			setValidFormats_("out", StringList::create("idXML"));

			registerSubsection_("Mascot_server", "Mascot server details");
			registerSubsection_("Mascot_parameters", "Mascot parameters used for searching");
		}

    Param getSubsectionDefaults_(const String& section) const
    {
			if (section == "Mascot_server")
			{
				MascotRemoteQuery mascot_query;
				return mascot_query.getParameters();
			}

			if (section == "Mascot_parameters")
			{
				MascotGenericFile mascot_infile;
				return mascot_infile.getParameters();
			}

      return Param();
    }


		ExitCodes main_(int argc, const char** argv)
		{
			//-------------------------------------------------------------
      // parameter handling
      //-------------------------------------------------------------

      //input/output files
			String in(getStringOption_("in")), out(getStringOption_("out"));
			FileHandler fh;
			FileTypes::Type in_type = fh.getType(in);

      //-------------------------------------------------------------
      // loading input
      //-------------------------------------------------------------

			PeakMap exp;
			fh.loadExperiment(in, exp, in_type, log_type_);

      //-------------------------------------------------------------
      // calculations
      //-------------------------------------------------------------

			Param mascot_param = getParam_().copy("Mascot_parameters:", true);
      MascotGenericFile mascot_infile;
			mascot_infile.setParameters(mascot_param);

			// get the spectra into string stream
			writeDebug_("Writing Mascot mgf file to stringstream", 1);
			stringstream ss;
			mascot_infile.store(ss, in, exp);

			// Usage of a QCoreApplication is overkill here (and ugly too), but we just use the
			// QEventLoop to process the signals and slots and grab the results afterwards from
			// the MascotRemotQuery instance
			char** argv2 = const_cast<char**>(argv);
			QCoreApplication event_loop(argc, argv2);
			MascotRemoteQuery* mascot_query = new MascotRemoteQuery(&event_loop);
			Param mascot_query_param = getParam_().copy("Mascot_server:", true);
			writeDebug_("Setting parameters for Mascot query", 1);
			mascot_query->setParameters(mascot_query_param);
			writeDebug_("Setting spectra for Mascot query", 1);
			mascot_query->setQuerySpectra(ss.str());

			// remove unnecessary spectra
			ss.clear();

		  QObject::connect(mascot_query, SIGNAL(done()), &event_loop, SLOT(quit()));
			QTimer::singleShot(1000, mascot_query, SLOT(run()));
			writeDebug_("Fire off Mascot query", 1);
			event_loop.exec();
			writeDebug_("Mascot query finished", 1);

			if (mascot_query->hasError())
			{
				writeLog_("An error occurred during the query: " + mascot_query->getErrorMessage());
				delete mascot_query;
				return EXTERNAL_PROGRAM_ERROR;
			}

			// write Mascot response to file
			String unique_name = File::getUniqueName(); // body for the tmp files
			String mascot_tmp_file_name(unique_name + "_Mascot_response");
			QFile mascot_tmp_file(mascot_tmp_file_name.c_str());
			mascot_tmp_file.open(QIODevice::WriteOnly);
			mascot_tmp_file.write(mascot_query->getMascotXMLResponse());
			mascot_tmp_file.close();

			// clean up
			delete mascot_query;

			vector<PeptideIdentification> pep_ids;
			ProteinIdentification prot_id;

			// read the response
			MascotXMLFile().load(mascot_tmp_file_name, prot_id, pep_ids);

			// delete file
			mascot_tmp_file.remove();

			vector<ProteinIdentification> prot_ids;
			prot_ids.push_back(prot_id);

			writeDebug_("Read " + String(pep_ids.size()) + " peptide ids and " + String(prot_id.getHits().size()) + " protein identifications", 5);


      //-------------------------------------------------------------
      // writing output
      //-------------------------------------------------------------

      IdXMLFile().store(out, prot_ids, pep_ids);

			return EXECUTION_OK;
		}

};


int main( int argc, const char** argv )
{
	TOPPMascotAdapterOnline tool;

	return tool.main(argc,argv);
}

/// @endcond
