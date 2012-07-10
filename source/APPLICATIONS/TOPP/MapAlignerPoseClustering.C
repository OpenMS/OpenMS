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
// $Maintainer: Chris Bielow $
// $Authors: Marc Sturm, Clemens Groepl, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmPoseClustering.h>
#include <OpenMS/APPLICATIONS/MapAlignerBase.h>

#ifdef _OPENMP 
#include <omp.h>
#endif

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------

/**
  @page TOPP_MapAlignerPoseClustering MapAlignerPoseClustering

  @brief Corrects retention time distortions between maps, using a pose clustering approach.

<CENTER>
  <table>
    <tr>
      <td ALIGN = "center" BGCOLOR="#EBEBEB"> potential predecessor tools </td>
      <td VALIGN="middle" ROWSPAN=2> \f$ \longrightarrow \f$ MapAlignerPoseClustering \f$ \longrightarrow \f$</td>
      <td ALIGN = "center" BGCOLOR="#EBEBEB"> potential successor tools </td>
    </tr>
    <tr>
      <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_FeatureFinderCentroided @n (or another feature finding algorithm) </td>
      <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_FeatureLinkerUnlabeled or @n @ref TOPP_FeatureLinkerUnlabeledQT </td>
    </tr>
  </table>
</CENTER>

  This tool provides an algorithm to align the retention time scales of multiple input files, correcting shifts and distortions between them.
  Retention time adjustment may be necessary to correct for chromatography differences e.g. before data from multiple LC-MS runs can be combined
  (feature grouping), or when one run should be annotated with peptide identifications obtained in a different run.

  All map alignment tools (MapAligner...) collect retention time data from the input files and - by fitting a model to this data
  - compute transformations that map all runs to a common retention time scale. They can apply the transformations right away and
  return output files with aligned time scales (parameter @p out), and/or return descriptions of the transformations in trafoXML
  format (parameter @p trafo_out). Transformations stored as trafoXML can be applied to arbitrary files with the @ref TOPP_MapRTTransformer tool.

  The map alignment tools differ in how they obtain retention time data for the modeling of transformations, and consequently what types
  of data they can be applied to. The alignment algorithm implemented here is the pose clustering algorithm as described in
  doi:10.1093/bioinformatics/btm209. It is used to find an affine transformation, which is further refined by a feature grouping step.
  This algorithm can be applied to features (featureXML) and peaks (mzML), but it has mostly been developed and tested on features.
  For more details and algorithm-specific parameters (set in the INI file) see "Detailed Description" in the
  @ref OpenMS::MapAlignmentAlgorithmPoseClustering "algorithm documentation".

  @see @ref TOPP_MapAlignerPoseClustering @ref TOPP_MapAlignerSpectrum @ref TOPP_MapRTTransformer

  This algorithm uses an affine transformation model.

  To speed up the alignment, consider reducing 'max_number_of_peaks_considered'.
  If your alignment is not good enough, consider increasing this number (the alignment will take longer though).

  <B>The command line parameters of this tool are:</B> @n
  @verbinclude TOPP_MapAlignerIdentification.cli
	<B>INI file documentation of this tool:</B>
	@htmlinclude TOPP_MapAlignerIdentification.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPMapAlignerPoseClustering
: public TOPPMapAlignerBase
{

public:
  TOPPMapAlignerPoseClustering() :
    TOPPMapAlignerBase("MapAlignerPoseClustering", "Corrects retention time distortions between maps using a pose clustering approach.")
  {}

protected:
  void registerOptionsAndFlags_()
  {
    TOPPMapAlignerBase::registerOptionsAndFlags_("mzML,featureXML", true);
    registerSubsection_("algorithm", "Algorithm parameters section");
  }

  Param getSubsectionDefaults_(const String & section) const
  {
    if (section == "algorithm")
    {
      MapAlignmentAlgorithmPoseClustering algo;
      return algo.getParameters();
    }
    return Param();     // shouldn't happen
  }

  ExitCodes main_(int, const char **)
  {
    MapAlignmentAlgorithmPoseClustering algorithm;
    ExitCodes ret = TOPPMapAlignerBase::initialize_(&algorithm, true);
    if (ret!=EXECUTION_OK) return ret;
    
    StringList in_files = getStringList_("in");
    StringList out_files = getStringList_("out");
    StringList out_trafos = getStringList_("trafo_out");

    Size reference_index = getIntOption_("reference:index");
    String reference_file = getStringOption_("reference:file");

    FileTypes::Type in_type = FileHandler::getType(in_files[0]);
    String file;
    if (!reference_file.empty())
    {
      file = reference_file;
      reference_index = in_files.size(); // points to invalid index
    }
    else if (reference_index > 0)  //  normal reference (index was checked before)
    {
      file = in_files[--reference_index]; // ref index is 1-based in parameters, but should be 0-based here
    }
    else if (reference_index == 0) // no reference given
    {
			LOG_INFO << "Picking a reference (by size) ..." << std::flush;
      // use map with highest number of features as reference:
      Size max_count(0);
      FeatureXMLFile f;
      for (Size m = 0; m < in_files.size(); ++m)
      {
        Size s(0);
        if (in_type==FileTypes::FEATUREXML) s = f.loadSize(in_files[m]);
        else if (in_type==FileTypes::MZML)
        { // this is expensive!
          MSExperiment<> exp;
          MzMLFile().load(in_files[m], exp);
          exp.updateRanges(1);
          s = exp.getSize();
        }
        if (s > max_count)
        {
          max_count = s;
          reference_index = m;
        }
      }
			LOG_INFO << " done" << std::endl;
      file = in_files[reference_index];
    }

    FeatureXMLFile f_fxml;
    if (out_files.size()==0) // no need to store featureXML, thus we can load only minimum required information
    {
      f_fxml.getOptions().setLoadConvexHull(false);
      f_fxml.getOptions().setLoadSubordinates(false);
    }
    if (in_type==FileTypes::FEATUREXML)
    {
      FeatureMap<> map_ref;
      FeatureXMLFile f_fxml_tmp; // for the reference, we never need CH or subordinates
      f_fxml_tmp.getOptions().setLoadConvexHull(false);
      f_fxml_tmp.getOptions().setLoadSubordinates(false);
      f_fxml_tmp.load(file, map_ref);
      algorithm.setReference(map_ref);
    }
    else if (in_type==FileTypes::MZML)
    {
      MSExperiment<> map_ref;
      MzMLFile().load(file, map_ref);
      algorithm.setReference(map_ref);
    }

		ProgressLogger plog;
		plog.setLogType(log_type_);
		
		plog.startProgress(0, in_files.size(), "Aligning input maps");
		Size progress(0); // thread-safe progress
    // TODO: it should all work on FeatureXML files, since we might need them for output anyway. Converting to ConsensusXML is just wasting memory!
    #ifdef _OPENMP 
    #pragma omp parallel for schedule(dynamic, 1)
    #endif
    for (Int i=0; i<in_files.size(); ++i)
    {
      TransformationDescription trafo;
      if (in_type==FileTypes::FEATUREXML)
      {
        FeatureMap<> map;
        // workaround for loading: use temporary FeatureXMLFile since it is not thread-safe
        FeatureXMLFile f_fxml_tmp; // do not use OMP-firstprivate, since FeatureXMLFile has no copy c'tor
        f_fxml_tmp.getOptions() = f_fxml.getOptions();
        f_fxml_tmp.load(in_files[i], map);
        if (i==reference_index) trafo.fitModel("identity");
        else algorithm.align(map, trafo);
        if (out_files.size())
        {
          MapAlignmentTransformer::transformSingleFeatureMap(map, trafo);
          // annotate output with data processing info
          addDataProcessing_(map, getProcessingInfo_(DataProcessing::ALIGNMENT));
          f_fxml_tmp.store(out_files[i], map);
        }
      }
      else if (in_type==FileTypes::MZML)
      {
        MSExperiment<> map;
        MzMLFile().load(in_files[i], map);
        if (i==reference_index) trafo.fitModel("identity");
        else algorithm.align(map, trafo);
        if (out_files.size())
        {
          MapAlignmentTransformer::transformSinglePeakMap(map, trafo);
          // annotate output with data processing info
          addDataProcessing_(map, getProcessingInfo_(DataProcessing::ALIGNMENT));
          MzMLFile().store(out_files[i], map);
        }
      }      
      
      if (out_trafos.size())
      {
        TransformationXMLFile().store(out_trafos[i], trafo);
      }

		  #ifdef _OPENMP 
      #pragma omp critical (MAPose_Progress)
      #endif
			{
      plog.setProgress(++progress); // thread safe progress counter
			}

    }
    
		plog.endProgress();
    return EXECUTION_OK;

  }
};

int main(int argc, const char ** argv)
{
  TOPPMapAlignerPoseClustering tool;
  return tool.main(argc, argv);
}

/// @endcond
