// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/KERNEL/Feature.h>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmSH.h>

#include <OpenMS/APPLICATIONS/TOPPBase.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page TOPP_FeatureFinderSH FeatureFinderSH

  A feature finder based on the original SuperHirn codebase.

  Proteomics. 2007 Oct;7(19):3470-80.
  SuperHirn - a novel tool for high resolution LC-MS-based peptide/protein profiling.
  Mueller LN, Rinner O, Schmidt A, Letarte S, Bodenmiller B, Brusniak MY, Vitek O, Aebersold R, Müller M.
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

typedef FeatureFinderAlgorithmSH<Peak1D,Feature> FFSH;

class TOPPFeatureFinderSH
    : public TOPPBase
{
public:
    TOPPFeatureFinderSH()
        : TOPPBase("FeatureFinderSH","Finds mass spectrometric features in mass spectra.", false)
    {
          // TODO finish, make official, add to source/APPLICATIONS/ToolHandler.C
    }

protected:

    void registerOptionsAndFlags_()
    {
      registerInputFile_("in","<file>","","input profile data file ");
      setValidFormats_("in",StringList::create("mzML"));
      registerOutputFile_("out","<file>","","output peak file ");
      setValidFormats_("out",StringList::create("featureXML"));

      addEmptyLine_();
      addText_("Parameters for the peak picker algorithm can be given in the 'algorithm' part of INI file.");
      registerSubsection_("algorithm","Algorithm parameters section");
    }

    Param getSubsectionDefaults_(const String& /*section*/) const
    {
        return FFSH().getDefaults();
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
        MzMLFile mzMLFile;
        mzMLFile.setLogType(log_type_);
        MSExperiment<Peak1D > input;
        mzMLFile.getOptions().addMSLevel(1);
        mzMLFile.load(in, input);

        if (input.empty())
        {
            LOG_WARN << "The given file does not contain any conventional peak data, but might"
                    " contain chromatograms. This tool currently cannot handle them, sorry.";
            return INCOMPATIBLE_INPUT_DATA;
        }

        //check if spectra are sorted
        for (Size i=0; i<input.size(); ++i)
        {
            if (!input[i].isSorted())
            {
                writeLog_("Error: Not all spectra are sorted according to peak m/z positions. Use FileFilter to sort the input!");
                return INCOMPATIBLE_INPUT_DATA;
            }
        }

        //-------------------------------------------------------------
        // pick
        //-------------------------------------------------------------
        FeatureMap<> output;

        FeatureFinder ff;
        Param param = getParam_().copy("algorithm:",true);

        FFSH ffsh;
        ffsh.setParameters(param);
        ffsh.setData(input, output, ff);
        ffsh.run();

        //-------------------------------------------------------------
        // writing output
        //-------------------------------------------------------------
        //annotate output with data processing info
        addDataProcessing_(output, getProcessingInfo_(DataProcessing::PEAK_PICKING));
        FeatureXMLFile().store(out,output);

        return EXECUTION_OK;
    }
};

int main( int argc, const char** argv )
{
    TOPPFeatureFinderSH tool;
    return tool.main(argc,argv);
}

/// @endcond
