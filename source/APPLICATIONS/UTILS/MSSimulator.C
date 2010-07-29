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
// $Maintainer: Stephan Aiche$
// $Authors: Ole Schulz-Trieglaff, Stephan Aiche, Chris Bielow $
// --------------------------------------------------------------------------

#include <iostream>
#include <string>

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FASTAFile.h>

#include <OpenMS/SYSTEM/StopWatch.h>

#include <OpenMS/SIMULATION/MSSim.h>
#include <OpenMS/SIMULATION/SimTypes.h>

// GSL includes (random number generation)
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

// file types
#include <OpenMS/FORMAT/DTA2DFile.h>
#include <OpenMS/FORMAT/MzXMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
using namespace OpenMS;
using namespace std;


/**
	@page UTILS_MSSimulator MSSimulator
		
	@brief A highly configurable simulator for mass spectrometry experiments.
	
  This implementation is a rewritten and extended version of the concepts and ideas presented in:<br>
  <p>
  Ole Schulz-Trieglaff, Nico Pfeifer, Clemens Gropl, Oliver Kohlbacher, and Knut Reinert.<br>
  LC-MSsim - A simulation software for liquid chromatography mass spectrometry data.<br>
  <em>BMC Bioinformatics</em> <b>9</b>:423, 2008.
  </p>

  The electronic version of this article is the complete one and can be found online at: <br>
  <a href="http://www.biomedcentral.com/1471-2105/9/423">http://www.biomedcentral.com/1471-2105/9/423</a>

  Added features are:
  <ul>
    <li>Simulation of MS/MS spectra with configurable precursor-selection strategy</li>
    <li>Simulation of Capillary electrophoresis as separation step</li>
    <li>Simulation of iTRAQ labels</li>
    <li>Simulation of 1D spectra</li>
  </ul>

	@note This tool is experimental!	

	<B>The command line parameters of this tool are:</B>
	@verbinclude UTILS_MSSimulator.cli
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPMSSimulator
	: public TOPPBase
{
  public:
		TOPPMSSimulator()
    : TOPPBase("MSSimulator","A highly configurable simulator for mass spectrometry experiments.",false)
		{
		}
    
  protected:
    
    void registerOptionsAndFlags_()
    {
      // I/O settings
      registerInputFileList_("in","<files>",StringList::create(""),"Input protein sequences in FASTA format",true,false);
      registerOutputFile_("out","<file>","","output (simulated MS map) in mzML format",true);
      registerOutputFile_("out_fm","<file>","","output (simulated MS map) in featureXML format",false);
      registerOutputFile_("out_cm","<file>","","output (simulated MS map) in consensusXML format (grouping charge variants from a parent peptide from ESI)",false);

      registerStringOption_("type","<name>","","Labeling type\n",true);
      setValidStrings_("type", getUtilList()[toolName_()] );

			addEmptyLine_();
  		addText_("To specify intensity values for certain proteins,\nadd an abundance tag for the corresponding protein\nin the FASTA input file:");
			addEmptyLine_();
  		addText_("- add '[# xx]' at the end of the > line to specify");
  		addText_("  xx total abundance units.");
			addEmptyLine_();
			addText_("e.g. >seq2 optional comment [#45]");
			addText_("     ASQKRPSQRHGSKYLATASTMDHARHGFLPRHRDTGILDSIGRFFGGDRGAPK");


      registerSubsection_("algorithm","Algorithm parameters section");    
    }
  
    Param getSubsectionDefaults_(const String& /*section*/) const
    { 
      Param tmp;
      String type = getStringOption_("type");
      tmp.insert("MSSim:", MSSim().getParameters(type));
      return tmp;
    }
  
  
    // Load proteins from FASTA file
    void loadFASTA_(const String& filename, SampleProteins & proteins )
    {
      writeLog_(String("Loading sequence data from ") + filename +  String(" ..") );
      
      FASTAFile fastafile;
      typedef std::vector< FASTAFile::FASTAEntry > FASTAdata;
      FASTAdata fastadata;
      
      // load FASTA file contents
      fastafile.load(filename, fastadata);
           
      // add data from file to protein storage
      String::size_type index;
            
      // re-parse fasta description to obtain quantitation info
      for (FASTAdata::iterator it = fastadata.begin(); it != fastadata.end(); ++it)
      {
        // parsed abundance
        SimIntensityType abundance = 100;


        // remove all ambiguous characters from FASTA entry
        // TODO: this is somehow problematic since we modfiy user input
        it->sequence.remove('X');
        it->sequence.remove('B');
        it->sequence.remove('Z');
        
        // Look for a relative quantity given in the comment line of a FASTA entry
				// e.g. >BSA [#120]
				index = (it->description).find("[#");
        // if found, extract and set relative quantity accordingly
        if (index != string::npos)
        {
					String::size_type index_end = (it->description).find(']', index);
					if (index_end == string::npos) throw Exception::InvalidParameter(__FILE__, __LINE__, __PRETTY_FUNCTION__,"MSSimulator: Invalid entry (" + it->identifier + ") in FASTA file; abundance section has open tag '[#' but missing close tag ']'.");
					
					//std::cout << (it->description).substr(index+2,index_end-index-2) << std::endl;
					StringList abundances = StringList::create((it->description).substr(index+2,index_end-index-2));
					if (abundances.size() == 0) throw Exception::InvalidParameter(__FILE__, __LINE__, __PRETTY_FUNCTION__,"MSSimulator: Invalid entry (" + it->identifier + ") in FASTA file; abundance section is missing abundance value.");
          abundance = abundances[0].toDouble();
					

					if (abundances.size() > 1)
          {	// additional abundances (e.g. iTRAQ) given... this is not supported (new syntax required)
            throw Exception::InvalidParameter(__FILE__, __LINE__, __PRETTY_FUNCTION__, String("MSSimulator (line ") + __LINE__ + "): Invalid entry (" + it->identifier + ") in FASTA file.");
					}
		    }
        
        proteins.push_back(make_pair(*it, abundance ));
      }
      
      writeLog_(String("done (") + fastadata.size() + String(" protein(s) loaded)"));
    }
	
		ExitCodes main_(int, const char**)
		{
			//-------------------------------------------------------------
			// parsing parameters
			//-------------------------------------------------------------
      String labeling_type = getStringOption_("type");

      StringList input_files = getStringList_("in");
			String outputfile_name = getStringOption_("out");	

      MSSim ms_simulation;
      ms_simulation.setParameters(getParam_().copy("algorithm:MSSim:",true));
			
      // read proteins 
      SampleChannels channels;
      for(Size i = 0 ; i < input_files.size() ; ++i)
      {
        SampleProteins proteins;
        loadFASTA_(input_files[i], proteins);
        channels.push_back(proteins);
      }
      
      // initialize the random number generator
      gsl_rng_default_seed = time(0);
      gsl_rng* rnd_gen_ = gsl_rng_alloc(gsl_rng_mt19937);
      
      // read contaminants

      // select contaminants?? -> should this be done by MSSim??
      
      // start simulation
      writeLog_("Starting simulation");
      StopWatch w;

      w.start();
      ms_simulation.simulate(rnd_gen_, channels, labeling_type);
      w.stop();
			writeLog_(String("Simulation took ") + String(w.getClockTime()) + String(" seconds"));   	  	
      
      writeLog_(String("Storing simulated map in: ") + outputfile_name);
      MzMLFile().store(outputfile_name, ms_simulation.getExperiment());
      
      String fxml_out = getStringOption_("out_fm");
			if (fxml_out != "" && File::writable(fxml_out))
			{
				writeLog_(String("Storing simulated features in: ") + fxml_out);
				FeatureXMLFile().store(fxml_out, ms_simulation.getSimulatedFeatures());
			}

      String cxml_out = getStringOption_("out_cm");
			if (cxml_out != "" && File::writable(cxml_out))
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

/// @endcond
