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
// $Maintainer: Andreas Bertsch$
// $Authors: $
// --------------------------------------------------------------------------


#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <OpenMS/CONCEPT/Factory.h>
#include <OpenMS/FILTERING/TRANSFORMERS/PreprocessingFunctor.h>

#include <OpenMS/FORMAT/MzMLFile.h>

#include <typeinfo>

using namespace OpenMS;
using namespace std;

/**
	@page TOPP_SpectraFilter SpectraFilter

	@brief Applies different spectrum modification filters to the data, i.e. removes certain peaks from the spectra given.

<CENTER>
	<table>
		<tr>
			<td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
			<td VALIGN="middle" ROWSPAN=2> \f$ \longrightarrow \f$ SpectraFilter \f$ \longrightarrow \f$</td>
			<td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
		</tr>
		<tr>
			<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_PeakPicker </td>
			<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> any tool operating on MS peak data @n (in mzML format)</td>
		</tr>
	</table>
</CENTER>

	The filters available will do the following:
	<UL>
		<LI> BernNorm -- does the Bern et al. normalization
		<LI> MarkerMower -- removes all peaks marked (e.g. with NeutralLossMarker)
		<LI> NLargest -- keeps the n most intensive peaks of each spectrum
		<LI> Normalizer -- normalizes the peaks in the spectrum with different modes (to_one, to_TIC)
		<LI> ParentPeakMower -- reduces the intensity of the parent peak
		<LI> Scaler -- scales the peaks according to their rank
		<LI> SqrtMower -- set each intensity to the square root of the original intensity
		<LI> ThresholdMower -- removes all peaks below a Threshold
		<LI> WindowMower -- keeps the biggest peaks in a sliding window
	</UL>

	<B>The command line parameters of this tool are:</B>
	@verbinclude TOPP_SpectraFilter.cli

	For the parameters of the algorithm section see the class documentation of each filter
	respectively: @n
		@ref OpenMS::BernNorm @n
		@ref OpenMS::MarkerMower @n
		@ref OpenMS::NLargest @n
		@ref OpenMS::Normalizer @n
		@ref OpenMS::ParentPeakMower @n
		@ref OpenMS::Scaler @n
		@ref OpenMS::SqrtMower @n
		@ref OpenMS::ThresholdMower @n
		@ref OpenMS::WindowMower @n
*/


// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPSpectraFilter
	: public TOPPBase
{
	public:
		TOPPSpectraFilter()
			: TOPPBase("SpectraFilter", "Applies a filter to peak spectra.")
		{
		}

	protected:

		void registerOptionsAndFlags_()
		{
			registerInputFile_("in", "<file>", "", "input file ");
			setValidFormats_("in",StringList::create("mzML"));
			registerOutputFile_("out", "<file>", "", "output file ");
	  	setValidFormats_("out",StringList::create("mzML"));
			registerStringOption_("type","<name>","","Filter type",true);
			setValidStrings_("type", ToolHandler::getTypes(toolName_()));

			addEmptyLine_();
			addText_("Parameters for the filter can only be given in the INI file.");

			// register one section for each algorithm
			registerSubsection_("algorithm","Algorithm parameter subsection.");
		}

		Param getSubsectionDefaults_(const String& /*section*/) const
		{
			String type = getStringOption_("type");
			return Factory<PreprocessingFunctor>::create(type)->getParameters();
		}

		ExitCodes main_(int , const char**)
		{
			//-------------------------------------------------------------
			// parameter handling
			//-------------------------------------------------------------

			//input/output files
			String in(getStringOption_("in"));
			String out(getStringOption_("out"));
			String type = getStringOption_("type");

      //-------------------------------------------------------------
      // loading input
      //-------------------------------------------------------------

      MSExperiment<> exp;
      MzMLFile f;
      f.setLogType(log_type_);
      f.load(in, exp);

      //-------------------------------------------------------------
      // if meta data arrays are present, remove them and warn
      //-------------------------------------------------------------
			if (exp.clearMetaDataArrays())
			{
				writeLog_("Warning: Spectrum meta data arrays cannot be sorted. They are deleted.");
			}

      //-------------------------------------------------------------
      // filter
      //-------------------------------------------------------------
			Param filter_param = getParam_().copy("algorithm:", true);
			writeDebug_("Used filter parameters", filter_param, 3);
			PreprocessingFunctor* filter = Factory<PreprocessingFunctor>::create(type);
			filter->setParameters(filter_param);
			filter->filterPeakMap(exp);

			//-------------------------------------------------------------
			// writing output
			//-------------------------------------------------------------

			//annotate output with data processing info
			addDataProcessing_(exp, getProcessingInfo_(DataProcessing::FILTERING));

			f.store(out, exp);

			return EXECUTION_OK;
		}
};

/// @endcond


int main( int argc, const char** argv )
{
	TOPPSpectraFilter tool;
	return tool.main(argc,argv);
}

