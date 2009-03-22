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
// $Maintainer: Stephan Aiche$
// $Authors: Ole Schulz-Trieglaff $
// --------------------------------------------------------------------------

#include <iostream>
#include <string>

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/FileHandler.h>

#include <OpenMS/SYSTEM/StopWatch.h>

#include <OpenMS/SIMULATION/LCMSSample.h>
#include <OpenMS/SIMULATION/LCMSSim.h>


using namespace OpenMS;
using namespace std;


class TOPPMapSimulator
	: public TOPPBase
{
	public:
		TOPPMapSimulator()
			: TOPPBase("MapSimulator","This application simulates an LC-MS run.",false)
		{ }

	protected:
	
	void registerOptionsAndFlags_()
	{
		// I/O settings
		registerStringOption_("in","<file>","","input protein sequences in FASTA format",true);
		registerStringOption_("out","<file>","","output (simulated LC-MS map) in mzData format",true);
		registerStringOption_("rt_model","<file>","none","SVM model for retention time prediction ",true);
		registerStringOption_("pd_model","<file>","none","SVM model for peptide detectability prediction ",true);
		
 		registerSubsection_("simulation","LC-MS simulation settings");
		registerSubsection_("sample","Protein sample settings");
	}
		
	Param getSubsectionDefaults_(const String& section) const
	{
		if (section == "simulation")
		{
 			return LCMSSim().getDefaults(); 
		}
		else if (section == "sample")
		{
			return LCMSSample().getDefaults();
		}
		
		// should not happen
		Param tmp;
		return tmp;
	}
	

		ExitCodes main_(int, const char**)
		{
			//-------------------------------------------------------------
			// parsing parameters
			//-------------------------------------------------------------
						
			String inputfile_name = getStringOption_("in");			
			inputFileReadable_(inputfile_name);
			String outputfile_name = getStringOption_("out");	
			outputFileWritable_(outputfile_name);
			
			String rtmodel_file = getStringOption_("rt_model");
      if(rtmodel_file != "none" && rtmodel_file != "1D")
      {
        inputFileReadable_(rtmodel_file);
      }
			
			String pdmodel_file = getStringOption_("pd_model");
			inputFileReadable_(pdmodel_file);

			writeLog_(String("Reading from file: ") + inputfile_name);
			//-------------------------------------------------------------
			// Init simulation
			//-------------------------------------------------------------
			
			LCMSSim sim;
						
			Param const& sim_param = getParam_().copy("simulation:",true);
			sim.setParameters(sim_param);
			sim.setRTModelFile(rtmodel_file);
			
			StopWatch w;
			{			
			w.start();
			LCMSSample sample;
			Param const& digest_param = getParam_().copy("sample:",true);
			sample.setParameters(digest_param);
			sample.setPdModelFile(pdmodel_file);
			sample.loadFASTA(inputfile_name);
			sample.digest();
			sample.clearProteins();
			
			writeLog_(String("Peptides: ") + String(sample.size()) );
			
			// debug output
// 			sample.printProteins();
// 			sample.printPeptides();
			sim.setSample(sample);
			w.stop();
			writeLog_(String("Pre-processing took ") + String(w.getClockTime()) + String(" seconds"));   	  	
			}
			
			//-------------------------------------------------------------
			// Run simulation
			//-------------------------------------------------------------
			w.reset();
			w.start();			
			sim.run();
			w.stop();
			writeLog_(String("Simulation took ") + String(w.getClockTime()) + String(" seconds"));   	  	
			//-------------------------------------------------------------
			// output
			//-------------------------------------------------------------
			writeLog_(String("Writing LC-MS map and feature lists.. "));
			
			sim.exportMzData(outputfile_name);
			sim.exportFeatureMap(outputfile_name);		
			
			return EXECUTION_OK;
		}
};


int main( int argc, const char** argv )
{
	TOPPMapSimulator tool;
	return tool.main(argc,argv);
}

