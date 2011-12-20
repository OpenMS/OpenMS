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
// $Maintainer: Alexandra Zerck $
// $Authors: Eva Lange $
// --------------------------------------------------------------------------
#include <OpenMS/config.h>

#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FILTERING/SMOOTHING/GaussFilter.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/PeakTypeEstimator.h>
#include <OpenMS/DATASTRUCTURES/StringList.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page TOPP_NoiseFilterGaussian NoiseFilterGaussian

	@brief  Executes a Gaussian filter to reduce the noise in an MS experiment.

<CENTER>
  <table>
    <tr>
      <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
      <td VALIGN="middle" ROWSPAN=4> \f$ \longrightarrow \f$ NoiseFilterGaussian \f$ \longrightarrow \f$</td>
      <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
    </tr>
    <tr>
      <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_FileConverter </td>
      <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_PeakPickerWavelet</td>
    </tr>
    <tr>
      <td VALIGN="middle" ALIGN = "center" ROWSPAN=2> @ref TOPP_Resampler </td>
      <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_PeakPickerHiRes</td>
    </tr>
    <tr>
      <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_BaselineFilter</td>
    </tr>
  </table>
</CENTER>

	The Gaussian filter is a peak area preserving low-pass filter and is characterized by narrow bandwidths,
	sharp cutoffs, and low passband ripple.

	@note The Gaussian filter works for uniform as well as for non-uniform data.

	<B>The command line parameters of this tool are:</B>
  @verbinclude TOPP_NoiseFilterGaussian.cli

	<B>The algorithm parameters for the Gaussian filter are:</B>
  @htmlinclude OpenMS_GaussFilter.parameters

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES


class TOPPNoiseFilterGaussian
      : public TOPPBase
{
  public:
    TOPPNoiseFilterGaussian()
        : TOPPBase("NoiseFilterGaussian","Removes noise from profile spectra by using Gaussian filter.")
    {
    }

    void registerOptionsAndFlags_()
    {
	  	registerInputFile_("in","<file>","","input raw data file ");
			setValidFormats_("in",StringList::create("mzML"));
			registerOutputFile_("out","<file>","","output raw data file ");
	  	setValidFormats_("out",StringList::create("mzML"));
			addEmptyLine_();
	  	addText_("Parameters for the algorithms can be given in the INI file only.");
			addEmptyLine_();
			addText_("Note: The Gaussian filter works for uniform as well as for non-uniform data.");
    	registerSubsection_("algorithm","Algorithm parameters section");
    }

    Param getSubsectionDefaults_(const String& /*section*/) const
    {
      return GaussFilter().getDefaults();
    }

    ExitCodes main_(int , const char**)
    {
      //-------------------------------------------------------------
      // parameter handling
      //-------------------------------------------------------------
      String in = getStringOption_("in");
      String out = getStringOption_("out");

      //-------------------------------------------------------------
      // loading input
      //-------------------------------------------------------------
      MzMLFile mz_data_file;
      mz_data_file.setLogType(log_type_);
      MSExperiment<Peak1D> exp;
      mz_data_file.load(in,exp);

			if (exp.empty())
			{
				LOG_WARN << "The given file does not contain any conventional peak data, but might"
					          " contain chromatograms. This tool currently cannot handle them, sorry.";
				return INCOMPATIBLE_INPUT_DATA;
			}
			//check for peak type (profile data required)
			if (PeakTypeEstimator().estimateType(exp[0].begin(),exp[0].end())==SpectrumSettings::PEAKS)
			{
				writeLog_("Warning: OpenMS peak type estimation indicates that this is not profile data!");
			}

			//check if spectra are sorted
			for (Size i=0; i< exp.size(); ++i)
			{
				if (!exp[i].isSorted())
				{
					writeLog_("Error: Not all spectra are sorted according to peak m/z positions. Use FileFilter to sort the input!");
					return INCOMPATIBLE_INPUT_DATA;
				}
			}

      //-------------------------------------------------------------
      // calculations
      //-------------------------------------------------------------
    	Param filter_param = getParam_().copy("algorithm:",true);
			writeDebug_("Parameters passed to filter", filter_param,3);

      GaussFilter gauss;
      gauss.setLogType(log_type_);
      gauss.setParameters(filter_param);
      try
      {
        gauss.filterExperiment(exp);
      }
      catch(Exception::IllegalArgument& e)
      {
        writeLog_(String("Error: ") + e.getMessage()) ;
        return INCOMPATIBLE_INPUT_DATA;
      }

      //-------------------------------------------------------------
      // writing output
      //-------------------------------------------------------------
      
      //annotate output with data processing info
      addDataProcessing_(exp, getProcessingInfo_(DataProcessing::SMOOTHING));

      mz_data_file.store(out,exp);

      return EXECUTION_OK;
    }
};


int main( int argc, const char** argv )
{
  TOPPNoiseFilterGaussian tool;
  return tool.main(argc,argv);
}

/// @endcond
