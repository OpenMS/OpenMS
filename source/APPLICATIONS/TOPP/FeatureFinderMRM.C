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
// $Maintainer: Ole Schulz-Trieglaff $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <OpenMS/DATASTRUCTURES/ConvexHull2D.h>
#include <OpenMS/DATASTRUCTURES/DPosition.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/MzDataFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/Feature.h>

#include <iostream>
#include <fstream>
#include <math.h>

using namespace OpenMS;
using namespace std;

typedef MSExperiment<Peak1D> LCMSmap;
typedef LCMSmap::SpectrumType Spectrum;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page TOPP_FeatureFinderMRM FeatureFinderMRM
	
	@brief Peptide quantitation based on Multiple-Reaction-Monitoring (MRM).
	
	Multiple-Reaction-Monitoring (MRM) is a method to quantify peptides based on
	peak signal intensities in MS/MS spectra.	In short, the abundance of a peptide
	is estimated by summing the intensities of selected fragment ions in its MS/MS spectra.
	
	The advantages of this method are high sensitivity and accuracy of quantitation.
	It stems from drug testing and research, but it is increasingly applied in proteomics.
	For details, see Anderson and Hunter (2005) "Quantitative Mass Spectrometric Multiple
  Reaction Monitoring Assays for Major Plasma Proteins" Molecular and Cellular Proteomics
  5.4 pp. 573-588.

	Obviously, we need to know the precursor masses and charges of the peptide ions
	we want to monitor as well as the m/z values of the fragment ion used for quantitation.
	These values are usually determined in a first LC-MS run of the sample using targeted 
	MS/MS.
	
	The input to this program consists of a list of precursor m/zs and fragment ion m/zs.
	It performs a quantitation as explained above and writes a list of peptide features 
	with the estimate abundance.

	<B>The command line parameters of this tool are:</B>
	@verbinclude TOPP_FeatureFinderMRM.cli
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPFeatureFinderMRM
	: public TOPPBase
{
 public:
	TOPPFeatureFinderMRM()
		: TOPPBase("FeatureFinderMRM","Quantitates peptides based on multiple reaction monitoring.")
	{
	}
	
 protected:
 
	void registerOptionsAndFlags_()
	{
		registerInputFile_("in","<file>","","input file");
		setValidFormats_("in",StringList::create("mzData"));
		registerOutputFile_("out","<file>","","output file");
		setValidFormats_("out",StringList::create("featureXML"));
		
		registerDoubleOption_("p_mz_tol","<float>",50.0,"Precursor m/z tolerance (in ppm)",false);		
		registerDoubleOption_("msms_mz_tol","<float>",50.0,"Fragment ion m/z tolerance (in ppm)",false);
		
		registerFlag_("d","Write elution curves of fragment ions to file");
	
		addEmptyLine_();
		addText_("You have to define the list of precursor and fragment ion m/z values in the INI file.");
		
		registerSubsection_("precursor_mz_list","Precursor mz list");
		registerSubsection_("msms_mz_list","Fragment ion list");
			
	}
	
	Param getSubsectionDefaults_(const String& section) const
  {
  	Param tmp;
		
		// One fragment ion per precursor m/z, we do not check
		// for an equal number of entries.
		if (section == "precursor_mz_list")
		{
			tmp.setValue("1",1200.0);
			tmp.setValue("2",1500.0);
		}
		else if (section == "msms_mz_list")
		{
			tmp.setValue("1",300.0);
			tmp.setValue("2",420.0);
		}
		
		return tmp;
  }

	ExitCodes main_(int , const char**)
	{
		// input file names and types
		String in = getStringOption_("in");	
		String out = getStringOption_("out");
		
		// read input data
		MSExperiment<Peak1D> exp;
		
		MzDataFile mz_file;
		mz_file.setLogType(log_type_);
		mz_file.load(in,exp);
		
		Param const& prec_mzs = getParam_().copy("precursor_mz_list:",true);
		writeDebug_("precursor_mz_list:", prec_mzs, 2);
			
		Param const& msms_mzs = getParam_().copy("msms_mz_list:",true);
		writeDebug_("msms_mz_list:", msms_mzs, 2);

		Feature::CoordinateType prec_mz_tol    = getDoubleOption_("p_mz_tol");
		Feature::CoordinateType msms_mz_tol = getDoubleOption_("msms_mz_tol");
		
		//cout << "Tolerances : " << prec_mz_tol << " " << msms_mz_tol << endl;

		FeatureMap< > features;
		
		bool dump_profile = getFlag_("d");
		
		Param::ParamIterator pit1 = prec_mzs.begin();	
		Param::ParamIterator pit2 = msms_mzs.begin();	
		for (;pit1 != prec_mzs.end() && pit2 != msms_mzs.end(); ++pit1, ++pit2)
		{	
			Feature::CoordinateType p_mz        = pit1->value;
			Feature::CoordinateType msms_mz = pit2->value;
		
			Feature f;
			f.setMZ(p_mz);
			
			vector<DPosition<2> > points; // to estimate the convex hull
			
			Feature::CoordinateType rt_center = 0.0;
			Feature::IntensityType int_sum     = 0.0;
			UInt scan_count = 0;
						
			cout << "Searching for precursor " << p_mz << " and ms/ms ion " << msms_mz << endl;
		
			String rtfile = String("precursor_") + String(p_mz);
			ofstream db_out;
			if (dump_profile)
			{
				db_out.open(rtfile.c_str());
			}
		
			for(MSExperiment<Peak1D>::iterator sit = exp.begin();
	    	  	sit != exp.end();
						++sit)
			{
	
				Feature::CoordinateType mz_err = pow(10.0, 6.0) * (p_mz - sit->getPrecursorPeak().getPosition()[0]) / sit->getPrecursorPeak().getPosition()[0];
				
				if (sit->getMSLevel() == 2 && fabs(mz_err) <= prec_mz_tol) 
				{			
					rt_center += sit->getRT();
					++scan_count;
					
					for(Spectrum::const_iterator spit=sit->begin(); spit != sit->end();++spit)
					{			
						mz_err = pow(10.0, 6.0) * (msms_mz - spit->getMZ()) / spit->getMZ();
										
						if (fabs(mz_err) <= msms_mz_tol)
						{				
							int_sum += spit->getIntensity();			
							points.push_back( DPosition<2>(sit->getRT(),sit->getPrecursorPeak().getPosition()[0]) );
							if (dump_profile)
							{
								db_out << sit->getRT() << " " << 	spit->getIntensity() << endl;
							}
						}
					}			
				
				}	// end if (it->getMSLevel() == 2)
			} // end for each scan
			
			db_out.close();
			
			// estimate rt coordinate as medium rt of all MS/MS scans
			rt_center /= scan_count;
			f.setRT(rt_center);
			f.setIntensity(int_sum);
			
			// compute convex hull
			ConvexHull2D hull = points;
			f.getConvexHulls().push_back(hull);
			
			if (int_sum > 0)
			{
				features.push_back(f);
			}
		} // end while
		
		FeatureXMLFile().store(out,features);
		
		return EXECUTION_OK;
	}
};


int main( int argc, const char** argv )
{
	TOPPFeatureFinderMRM t;
	return t.main(argc,argv);
}

/// @endcond

