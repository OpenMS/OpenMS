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
// $Authors: Ole Schulz-Trieglaff, Stephan Aiche, Chris Bielow $
// --------------------------------------------------------------------------

#include <iostream>
#include <string>

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FASTAFile.h>

#include <OpenMS/SYSTEM/StopWatch.h>

#include <OpenMS/SIMULATION/LCMSSample.h>
#include <OpenMS/SIMULATION/LCMSSim.h>

#include <OpenMS/SIMULATION/MSSim.h>
#include <OpenMS/SIMULATION/SimTypes.h>

// GSL includes (random number generation)
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

// file types
#include <OpenMS/FORMAT/DTA2DFile.h>
#include <OpenMS/FORMAT/MzXMLFile.h>
#include <OpenMS/FORMAT/MzDataFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
using namespace OpenMS;
using namespace std;


class TOPPMSSimulator
	: public TOPPBase
{
  public:
		TOPPMSSimulator()
    : TOPPBase("MSSimulator","\n\nWARNING: EXPERIMENTAL\n\nHighly configurable simulator for mass spectrometry experiments.",false)
		{ }
    
  protected:
    
    void registerOptionsAndFlags_()
    {
      // I/O settings
      registerInputFile_("in","<file>","","Input protein sequences in FASTA format",true);
      registerOutputFile_("out","<file>","","output (simulated MS map) in mzData format",true);
      registerOutputFile_("out_fm","<file>","","output (simulated MS map) in featureXML format",false);
      registerOutputFile_("out_cm","<file>","","output (simulated MS map) in consensusXML format",false);

      registerSubsection_("algorithm","Algorithm parameters section");    
    }
  
    Param getSubsectionDefaults_(const String& /*section*/) const
    { 
      Param tmp;
      tmp.insert("MSSim:", MSSim().getParameters());      
      return tmp;
    }
  
  
    // Load proteins from FASTA file
    void loadFASTA(const String filename, SampleProteins & proteins)
    {
      writeLog_(String("Loading sequence data from ") + filename +  String(" ..") );
      
      FASTAFile fastafile;
      typedef std::vector< FASTAFile::FASTAEntry > FASTAdata;
      FASTAdata fastadata;
      
      // load FASTA file contents
      fastafile.load(filename, fastadata);
           
      // add data from file to protein storage
      String::size_type index;
      SimIntensityType relativeQuantity;
      for (FASTAdata::iterator it = fastadata.begin(); it != fastadata.end(); ++it)
      {
        // remove all ambiguous characters from FASTA entry
        // TODO: this is somehow problematic since we modfiy user input
        it->sequence.remove('X');
        it->sequence.remove('B');
        it->sequence.remove('Z');
        
        // Look for a relative quantity given in the first line of a FASTA entry
        index = (it->identifier).find_first_of("#");
        // if found, extract and set relative quantity accordingly
        if (index != string::npos)
        {
          stringstream strm ((it->identifier).substr(0, index));
          // we clean the identifier for later use
          it->identifier = (it->identifier).suffix('#');
          strm >> relativeQuantity;
        }
        else
        {
          relativeQuantity = 100;
        }
        AASequence aaseq(it->sequence);
        proteins.push_back(make_pair(*it, relativeQuantity ));
      }
      
      writeLog_(String("done (") + fastadata.size() + String(" protein(s) loaded)"));
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
			
      // read proteins 
      SampleProteins proteins;
      loadFASTA(inputfile_name,proteins);
      
      // initialize the random number generator
      gsl_rng_default_seed = time(0);
      gsl_rng* rnd_gen_ = gsl_rng_alloc(gsl_rng_mt19937);
      
      // read contaminants

      // select contaminants?? -> should this be done by MSSim??
      
      // start simulation
      writeLog_("Starting simulation");
      StopWatch w;
      MSSim ms_simulation;
      ms_simulation.setParameters(getParam_().copy("algorithm:MSSim:",true));

      w.start();
      ms_simulation.simulate(rnd_gen_, proteins);
      w.stop();
			writeLog_(String("Simulation took ") + String(w.getClockTime()) + String(" seconds"));   	  	
      
      writeLog_(String("Storing simulated map in: ") + outputfile_name);
      MzDataFile().store(outputfile_name, ms_simulation.getExperiment());
      
      String fxml_out = getStringOption_("out_fm");
			if (File::writable(fxml_out))
			{
				writeLog_(String("Storing simulated features in: ") + fxml_out);
				FeatureXMLFile().store(fxml_out, ms_simulation.getSimulatedFeatures());
			}

      String cxml_out = getStringOption_("out_cm");
			if (File::writable(cxml_out))
			{
				writeLog_(String("Storing simulated consensus features in: ") + cxml_out);
				ConsensusXMLFile().store(cxml_out, ms_simulation.getSimulatedConsensus());
			}
      
      // free random number generator
      gsl_rng_free(rnd_gen_);
      
			return EXECUTION_OK;
		}
};


int main( int argc, const char** argv )
{
	TOPPMSSimulator tool;
	return tool.main(argc,argv);
}

