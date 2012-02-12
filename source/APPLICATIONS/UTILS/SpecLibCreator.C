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
// $Maintainer: David Wojnar $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FORMAT/CsvFile.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/FORMAT/MSPFile.h>
#include <OpenMS/KERNEL/RichPeak1D.h>
#include <iostream>

#include <vector>
#include <cmath>
using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
  @page UTILS_SpecLibCreator SpecLibCreator
 
 	@brief creates with given data a msp format spectral library.

	@experimental This Utility is not well tested and some features might not work as expected.
	
	<B>The command line parameters of this tool are:</B>
	@verbinclude UTILS_SpecLibCreator.cli 
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPSpecLibCreator
	: public TOPPBase
	{
	public:
		TOPPSpecLibCreator()
		: TOPPBase("SpecLibCreator","Creates an MSP formated spectral library.",false)
		{
		}
		
	protected:
		void registerOptionsAndFlags_()
		{
			registerInputFile_("info","<file>","","Holds id, peptide, retention time etc.");
			//setValidFormats_("info",StringList::create("csv"));
			registerStringOption_("itemseperator","<char>",","," Seperator between items. e.g. ,",false);
			registerStringOption_("itemenclosed","<bool>","false","'true' or 'false' if true every item is enclosed e.g. '$peptide$,$run$...",false);
			//setValidFormats_("itemenclosed",StringList::create("true,false"));
			registerInputFile_("spec","<file>","","spectra");
			setValidFormats_("spec",StringList::create("mzData,mzXML"));
			registerOutputFile_("out","<file>","","output MSP formated spectra library");
		//	setValidFormats_("out",StringList::create("MSP"));

			addEmptyLine_();
			addText_("Note: information file should have the following information: peptide, retention time, measured weight, charge state");
			addText_("Extra information is allowed");
		}
		
		ExitCodes main_(int , const char**)
		{
			//-------------------------------------------------------------
			// parameter handling
			//-------------------------------------------------------------
		
		  String info = getStringOption_("info");
			String itemseperator = getStringOption_("itemseperator");
			String out = getStringOption_("out");
			bool itemenclosed;
			if(getStringOption_("itemenclosed") == "true")
			{
				itemenclosed  = true;
			}
			else
			{
				itemenclosed = false;
			}

			String spec = getStringOption_("spec");
			if(info == String::EMPTY)
			{
				throw Exception::RequiredParameterNotGiven(__FILE__,__LINE__,__PRETTY_FUNCTION__,"info");
			}
			if(spec == String::EMPTY)
			{
				throw Exception::RequiredParameterNotGiven(__FILE__,__LINE__,__PRETTY_FUNCTION__,"spec");
			}

			
			//-------------------------------------------------------------
			// loading input
			//-------------------------------------------------------------
			Int retention_time = -1 ;
			Int peptide =-1;
			Int measured_weight = -1;
			//UInt first_scan;
			UInt charge_state(0), Experimental_id(0);//,found_by, track, comment, vaccination_peptid,epitope, confident, hlaallele;
			const char* sepi = itemseperator.c_str();
			char sepo = *sepi;
			CsvFile csv_file(info,sepo,itemenclosed);
			vector<StringList>	list;

			list.resize(csv_file.size());

			for(UInt i= 0 ; i< csv_file.size(); ++i)
			{
				csv_file.getRow(i,list[i]);
			}
			for(UInt i = 0;i < list[0].size(); ++i)
			{

				if(list[0][i].toLower().removeWhitespaces().compare("retentiontime") == 0)
				{
					retention_time = i;
				}
				else if(list[0][i].toLower().hasSubstring("_id"))
				{
					Experimental_id = i;
				}
				else if(list[0][i].toLower() == "last scan")
				{
					// last_scan = i;
				}
				else if(list[0][i].toLower() == "modification")
				{
					// modification = i;
				}
				else if(list[0][i].toLower().removeWhitespaces().compare("chargestate")==0 || list[0][i].toLower().removeWhitespaces().hasSubstring("charge") )
				{
					charge_state = i;
				}
				else if(list[0][i].toLower().trim().compare("peptide") == 0)
				{
					peptide = i; 
				}
				else if(list[0][i].toLower().removeWhitespaces().hasSubstring("measuredweight")  || list[0][i].removeWhitespaces().compare("measuredweight[M+nH]n+") == 0)
				{
					measured_weight = i;
				}
			}
			if(retention_time  == -1)
			{
				throw Exception::RequiredParameterNotGiven(__FILE__,__LINE__,__PRETTY_FUNCTION__,"unclear which parameter is retention time");
			}
			if(peptide  == -1)
			{
				throw Exception::RequiredParameterNotGiven(__FILE__,__LINE__,__PRETTY_FUNCTION__,"unclear which parameter is peptide");
			}
			if(measured_weight  == -1)
			{
				throw Exception::RequiredParameterNotGiven(__FILE__,__LINE__,__PRETTY_FUNCTION__,"unclear which parameter is measured weight");
			}
			FileHandler fh;
			FileTypes::Type in_type = fh.getType(spec);
			/*MSExperiment<>*/PeakMap msexperiment;
			
			if(in_type == FileTypes::UNKNOWN)
			{
				writeLog_("Warning: Could not determine input file type!");
			}			
			else if(in_type == FileTypes::MZDATA)
			{
				MzDataFile mzData;
				mzData.load(spec,msexperiment);
			}
			else if(in_type == FileTypes::MZXML)
			{
				MzXMLFile mzXML;
				mzXML.load(spec, msexperiment);
			}
			if (msexperiment.getMinRT() == 0)
			{
				throw Exception::RequiredParameterNotGiven(__FILE__,__LINE__,__PRETTY_FUNCTION__,"EMPTY??");
			}			
			RichPeakMap library;
			
			//-------------------------------------------------------------
			// creating library
			//-------------------------------------------------------------
			UInt found_counter = 0;

				for(UInt i = 1; i < list.size(); ++i)
				{
					bool no_peptide = true;
					DoubleReal rt =  (60* (list[i][retention_time].toFloat())); // from minutes to seconds
					DoubleReal mz = list[i][measured_weight].toFloat();
					for(MSExperiment<>::Iterator it=msexperiment.begin(); it < msexperiment.end(); ++it)
					{
						//cout<<"i =" <<i<<endl; 
						//cout<<rt <<" (rt) - " << it->getRT()<<" (getRT) = "<<(rt - it->getRT())<<endl;
						if((abs(rt - it->getRT()) < 5) && (abs(mz - it->getPrecursors()[0].getMZ()) < 0.1 ))
						//if ( ceil(rt) == ceil(it->getRT()) || ceil(rt) == floor(it->getRT()) || floor(rt) == ceil(it->getRT()) || floor(rt) == floor(it->getRT()))
						{
							++found_counter;
							no_peptide = false;
							cout<<"Found Peptide " <<list[i][peptide] << " with id: " << list[i][Experimental_id]<<"\n";
							cout<<"rt: "<<it->getRT()<<" and mz: "<<it->getPrecursors()[0].getMZ()<<"\n";
							
						//		MSSpectrum<RichPeak1D> spec;
						//	for(UInt k = 0; k < it->size(); ++k)
						//		{
						//	spec.push_back(it->operator[](k));
						//			
						//		}
							MSSpectrum<RichPeak1D> speci;
							speci.setRT(it->getRT());
							speci.setMSLevel(2);
							speci.setPrecursors(it->getPrecursors());
							for(UInt j=0 ;j < it->size()	;++j)
							{
							
								RichPeak1D richy;
								richy.setIntensity(it->operator[](j).getIntensity());
								richy.setPosition(it->operator[](j).getPosition());
								richy.setMZ(it->operator[](j).getMZ());
								richy.setPos(it->operator[](j).getPos());//ALIAS for setMZ??? 
							
								speci.push_back(richy);
							}
							PeptideHit hit;// = *it->getPeptideIdentifications().begin()->getHits().begin();
							AASequence aa(list[i][peptide]);
							hit.setSequence(aa);
							hit.setCharge(list[i][charge_state].toInt());
							vector<PeptideHit> hits;
							hits.push_back(hit);
							vector<PeptideIdentification>pepi;
							PeptideIdentification pep;
							pep.setHits(hits);
							pepi.push_back(pep);
							speci.setPeptideIdentifications(pepi);
							//it->getPeptideIdentifications().begin()->setHits(hits);
							library.push_back(speci);
						}
					}
						if(no_peptide)
						{
							cout<<"Peptide: "<<list[i][peptide] <<" not found\n";
						}
				}
			cout<<"Found "<<found_counter <<" peptides\n";
			//library = static_cast<MSExperiment<MSSpectrum<RichPeak1D> > >(msexperiment);

			//-------------------------------------------------------------
			// writing output
			//-------------------------------------------------------------
			 in_type = fh.getType(out);
			if(in_type == FileTypes::MZDATA)
			{
				MzDataFile mzData;
				mzData.store(out,library);
			}
			else if(in_type == FileTypes::MZXML)
			{
				MzXMLFile mzXML;
				mzXML.store(out,library);
			}
			else
			{
				MSPFile msp;
				msp.store(out,library);
			}
			
			return EXECUTION_OK;
		}
		
	};




int main( int argc, const char** argv )
{
    TOPPSpecLibCreator tool;
    return tool.main(argc,argv);
}

/// @endcond

