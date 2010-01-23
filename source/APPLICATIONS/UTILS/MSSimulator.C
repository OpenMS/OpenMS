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
      registerInputFile_("in","<file>","","Input protein sequences in FASTA format",true);
      registerOutputFile_("out","<file>","","output (simulated MS map) in mzML format",true);
      registerOutputFile_("out_fm","<file>","","output (simulated MS map) in featureXML format",false);
      registerOutputFile_("out_cm","<file>","","output (simulated MS map) in consensusXML format (grouping charge variants from a parent peptide from ESI)",false);

			addEmptyLine_();
  		addText_("To specify intensity values for certain proteins,\nadd an abundance tag for the corresponding protein\nin the FASTA input file:");
			addEmptyLine_();
  		addText_("- add '[# xx]' at the end of the > line to specify");
  		addText_("  xx total abundance units.");
  		addText_("- add '[# xx, itraq113:yy, itraq115:zz]' to specify");
  		addText_("  xx total abundance units and yy itraq(channel 113) units");
  		addText_("  and zz itraq(channel 115) units. All Itraq is normalized");
  		addText_("  to 1 and xx is distributed accordingly.");
			addEmptyLine_();
			addText_("e.g. >seq2 optional comment [#45, itraq119:20, itraq121:25]");
			addText_("     ASQKRPSQRHGSKYLATASTMDHARHGFLPRHRDTGILDSIGRFFGGDRGAPK");


      registerSubsection_("algorithm","Algorithm parameters section");    
    }
  
    Param getSubsectionDefaults_(const String& /*section*/) const
    { 
      Param tmp;
      tmp.insert("MSSim:", MSSim().getParameters());      
      return tmp;
    }
  
  
    // Load proteins from FASTA file
    void loadFASTA_(const String& filename, const MSSim& mssim, SampleProteins & proteins )
    {
      writeLog_(String("Loading sequence data from ") + filename +  String(" ..") );
      
      FASTAFile fastafile;
      typedef std::vector< FASTAFile::FASTAEntry > FASTAdata;
      FASTAdata fastadata;
      
      // load FASTA file contents
      fastafile.load(filename, fastadata);
           
      // add data from file to protein storage
      String::size_type index;
      FASTAEntryEnhanced quant_info;

			// test for valid iTRAQ-channels
			String iTRAQ_status = mssim.getParameters().getValue("Global:iTRAQ");
      ItraqConstants::ChannelMapType abundance_itraq;
			if (iTRAQ_status!="off")
			{
				ItraqConstants::initChannelMap(iTRAQ_status=="4plex"? ItraqConstants::FOURPLEX : ItraqConstants::EIGHTPLEX, abundance_itraq);
			}
			//
      ItraqConstants::ChannelMapType abundance_itraq_8;
      ItraqConstants::initChannelMap(ItraqConstants::EIGHTPLEX, abundance_itraq_8);

      
      // re-parse fasta description to obtain quantitation info
      for (FASTAdata::iterator it = fastadata.begin(); it != fastadata.end(); ++it)
      {
        // remove all ambiguous characters from FASTA entry
        // TODO: this is somehow problematic since we modfiy user input
        it->sequence.remove('X');
        it->sequence.remove('B');
        it->sequence.remove('Z');
        
        // Look for a relative quantity given in the first line of a FASTA entry
				// e.g. [#120,itraq117:34,itraq119:23]
				index = (it->description).find("[#");
        // if found, extract and set relative quantity accordingly
        if (index != string::npos)
        {
					String::size_type index_end = (it->description).find(']', index);
					if (index_end == string::npos) throw Exception::InvalidParameter(__FILE__, __LINE__, __PRETTY_FUNCTION__,"MSSimulator: Invalid entry (" + it->identifier + ") in FASTA file; abundance section has open tag '[#' but missing close tag ']'.");
					
					//std::cout << (it->description).substr(index+2,index_end-index-2) << std::endl;
					StringList abundances = StringList::create((it->description).substr(index+2,index_end-index-2));
					if (abundances.size() == 0) throw Exception::InvalidParameter(__FILE__, __LINE__, __PRETTY_FUNCTION__,"MSSimulator: Invalid entry (" + it->identifier + ") in FASTA file; abundance section is missing abundance value.");
					quant_info["intensity"] = abundances[0].toDouble();
					
					if (abundances.size() > 1)
					{	// additional abundances (e.g. iTRAQ) given...
						
						//TODO: normalize itraq to 1
						SimIntensityType itraq_abundance_sum = 0;
						for (Size i = 1; i < abundances.size(); ++i)
						{
							StringList parts;
							abundances[i].split(':', parts);
							if (parts.size()!=2)
							{
								throw Exception::InvalidParameter(__FILE__, __LINE__, __PRETTY_FUNCTION__,"MSSimulator: Invalid entry (" + it->identifier + ") in FASTA file; expected one semicolon ('" + abundances[i] + "')");
							}
							parts[0] = parts[0].trim();
							parts[1] = parts[1].trim();
							if (parts[0]==String::EMPTY || parts[1]==String::EMPTY)
							{
								throw Exception::InvalidParameter(__FILE__, __LINE__, __PRETTY_FUNCTION__,"MSSimulator: Invalid entry (" + it->identifier + ") in FASTA file; key or value is empty ('" + abundances[i] + "')");
							}
							
							const string itraq_prefix = "itraq";
							// iTRAQ specific stuff
							if (parts[0].hasPrefix(itraq_prefix))
							{
								Int channel = parts[0].substr(itraq_prefix.length()).toInt();
								if (abundance_itraq_8.find(channel) == abundance_itraq_8.end())
								{
									throw Exception::InvalidParameter(__FILE__, __LINE__, __PRETTY_FUNCTION__, String("MSSimulator (line ") + __LINE__ + "): Invalid entry (" + it->identifier + ") in FASTA file; channel is not valid ('" + String(channel) + "')");
								}
								
								if (abundance_itraq.find(channel) != abundance_itraq.end())
								{	// current iTRAQ supports this channel --> save it
									quant_info["intensity_itraq"+String(channel)]	= parts[1].toDouble();
									itraq_abundance_sum += quant_info["intensity_itraq"+String(channel)];
								}
							}
							// ICAT?! ...
							/* else if (parts[0].hasPrefix("heavy"))  */
							else
							{
								throw Exception::InvalidParameter(__FILE__, __LINE__, __PRETTY_FUNCTION__,"MSSimulator: Invalid entry (" + it->identifier + ") in FASTA file; quantitation method not supported (" + parts[0] + ")");
							}							
							
							
						} // ! abundance values
						
						if (itraq_abundance_sum>0) // signals found
						{
							// normalize iTRAQ abundances to 1 and distribute overall abundance onto channels
							for (FASTAEntryEnhanced::iterator it_q = quant_info.begin(); it_q!=quant_info.end(); ++it_q)
							{
								if (it_q->first.hasPrefix("intensity_itraq"))
								{
									it_q->second = quant_info["intensity"] * (it_q->second / itraq_abundance_sum);
								}
							}
						}
						
					}
		    }
        else
        {
          quant_info["intensity"] = 100;
        }
        
        proteins.push_back(make_pair(*it, quant_info ));
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

      MSSim ms_simulation;
      ms_simulation.setParameters(getParam_().copy("algorithm:MSSim:",true));
			
      // read proteins 
      SampleProteins proteins;
      loadFASTA_(inputfile_name, ms_simulation, proteins);
      
      // initialize the random number generator
      gsl_rng_default_seed = time(0);
      gsl_rng* rnd_gen_ = gsl_rng_alloc(gsl_rng_mt19937);
      
      // read contaminants

      // select contaminants?? -> should this be done by MSSim??
      
      // start simulation
      writeLog_("Starting simulation");
      StopWatch w;

      w.start();
      ms_simulation.simulate(rnd_gen_, proteins);
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
