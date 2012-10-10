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
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/TraMLFile.h>
#include <OpenMS/ANALYSIS/TARGETED/TargetedExperiment.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>

using namespace OpenMS;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
  @page TOPP_MRMMapper MRMMapper

    @brief MRMMapper maps measured chromatograms (mzML) and the transitions used (TraML).
  
    This tool reads an mzML containing chromatograms (presumably measured on an
    SRM instrument) and a TraML file that contains the data that was used to
    generate the instrument method to measure said data. It then maps the
    transitions in the TraML file to the chromatograms found in the mzML file
    and stores the mapping by replacing the "id" paramter in the mzML with the
    "id" of the transition in the TraML file. It removes chromatograms for
    which it cannot find a mapping and throws an error if more than one
    transitions maps to a chromatogram.
    In strict mode (default) it also throws an error it not all chromatograms
    could be found in the TraML file.

    The thus mapped file can then be used in a downstream analysis.

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES


class TOPPMRMMapper 
  : public TOPPBase
{

public:

  TOPPMRMMapper() :
    TOPPBase("MRMMapper", "MRMMapper maps measured chromatograms (mzML) and the transitions used (TraML)")
  {
  }

protected:

  void registerOptionsAndFlags_()
  {
    registerInputFile_("in", "<file>", "", "Input file containing chromatograms (converted mzXML file)");
    setValidFormats_("in", StringList::create("mzML"));

    registerInputFile_("tr", "<file>", "", "transition file");
    setValidFormats_("tr", StringList::create("TraML"));

    registerOutputFile_("out", "<file>", "", "Output file containing mapped chromatograms");
    setValidFormats_("out", StringList::create("mzML"));

    registerDoubleOption_("precursor_tolerance", "<double>", 0.1, "Precursor tolerance when mapping (in Th)", false);
    registerDoubleOption_("product_tolerance", "<double>", 0.1, "Product tolerance when mapping (in Th)", false);

    registerFlag_("no-strict", "run in non-strict mode and allow some chromatograms to not be mapped.");
  }

  ExitCodes main_(int, const char **)
  {

    String in = getStringOption_("in");
    String tr_file = getStringOption_("tr");
    String out = getStringOption_("out");
    DoubleReal map_precursor_tol_ = getDoubleOption_("precursor_tolerance");
    DoubleReal map_product_tol_ = getDoubleOption_("product_tolerance");
    bool nostrict = getFlag_("no-strict");

    OpenMS::TargetedExperiment targeted_exp;
    OpenMS::MSExperiment<ChromatogramPeak> chromatogram_map;
    OpenMS::MSExperiment<ChromatogramPeak> output;

    TraMLFile().load(tr_file, targeted_exp);
    MzMLFile().load(in, chromatogram_map);

    // copy all meta data from old chromatogram
    output = chromatogram_map;
    output.clear(false);
    std::vector<MSChromatogram<ChromatogramPeak> > empty_chromats;
    output.setChromatograms(empty_chromats);

    int notmapped = 0;
    for (Size i = 0; i < chromatogram_map.getChromatograms().size(); i++)
    {
      // try to find the best matching transition for this chromatogram
      bool mapped_already = false;
      MSChromatogram<ChromatogramPeak> chromatogram = chromatogram_map.getChromatograms()[i];
      for (Size j = 0; j < targeted_exp.getTransitions().size(); j++)
      {

        if (fabs(chromatogram.getPrecursor().getMZ() - targeted_exp.getTransitions()[j].getPrecursorMZ()) < map_precursor_tol_ &&
            fabs(chromatogram.getProduct().getMZ()   - targeted_exp.getTransitions()[j].getProductMZ())   < map_product_tol_)
        {

          // std::cout << "Mapping chromatogram " << i << " to transition " << j << " (" << targeted_exp.getTransitions()[j].getNativeID() << ")"
          //    " with precursor mz " << chromatogram.getPrecursor().getMZ() << " / " <<  targeted_exp.getTransitions()[j].getPrecursorMZ() <<
          //    " and product mz " << chromatogram.getProduct().getMZ() << " / " <<  targeted_exp.getTransitions()[j].getProductMZ() << std::endl;

          // ensure: map every chromatogram to only one transition
          if (mapped_already)
          {
            throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Already mapped chromatogram " + String(i) + \
             " with " + String(chromatogram.getPrecursor().getMZ()) + \
              " -> " + String(chromatogram.getProduct().getMZ()) +  \
                "! Maybe try to decrease your mapping tolerance.");
          }
          mapped_already = true;

          // Create precursor and set the peptide sequence
          Precursor precursor = chromatogram.getPrecursor();
          String pepref = targeted_exp.getTransitions()[j].getPeptideRef();
          for (Size pep_idx = 0; pep_idx < targeted_exp.getPeptides().size(); pep_idx++)
          {
            const OpenMS::TargetedExperiment::Peptide * pep = &targeted_exp.getPeptides()[pep_idx];
            if (pep->id == pepref)
            {
              precursor.setMetaValue("peptide_sequence", pep->sequence);
              break;
            }
          }
          // add precursor to spectrum
          chromatogram.setPrecursor(precursor);

          // Set the id of the chromatogram, using the id of the transition (this gives directly the mapping of the two)
          chromatogram.setNativeID(targeted_exp.getTransitions()[j].getNativeID());
        }
      }

      // ensure: map every chromatogram to at least one transition
      if (!mapped_already)
      {
        std::cerr << "Did not find a mapping for chromatogram " + String(i) + " with " + String(chromatogram.getPrecursor().getMZ()) + \
          " -> " + String(chromatogram.getProduct().getMZ()) +  "! Maybe try to increase your mapping tolerance." << std::endl;
        notmapped++;
        if (!nostrict)
        {
          throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Did not find a mapping for chromatogram " + String(i) + "! Maybe try to increase your mapping tolerance.");
        }
      }
      else
      {
        output.addChromatogram(chromatogram);
      }
    }

    if (notmapped > 0)
    {
      std::cerr << "Could not find mapping for " << notmapped  << " chromatogram(s) " << std::endl;
    }

    // add all data processing information to all the chromatograms
    DataProcessing dp;
    dp = getProcessingInfo_(DataProcessing::FORMAT_CONVERSION);
    std::vector<MSChromatogram<ChromatogramPeak> > chromatograms = output.getChromatograms();
    for (Size i=0; i<chromatograms.size(); ++i)
    {
      chromatograms[i].getDataProcessing().push_back(dp);
    }
    output.setChromatograms(chromatograms);

    MzMLFile().store(out, output);
    return EXECUTION_OK;
  }

};

int main(int argc, const char ** argv)
{

  TOPPMRMMapper tool;
  return tool.main(argc, argv);
}

/// @endcond
