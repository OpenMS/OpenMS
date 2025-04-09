// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Marc Sturm, Clemens Groepl, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmPoseClustering.h>
#include <OpenMS/APPLICATIONS/MapAlignerBase.h>
//TODO remove when we get loadsize support in handler
#include <OpenMS/FORMAT/FeatureXMLFile.h>

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
      <th ALIGN = "center"> potential predecessor tools </td>
      <td VALIGN="middle" ROWSPAN=2> &rarr; MapAlignerPoseClustering &rarr;</td>
      <th ALIGN = "center"> potential successor tools </td>
    </tr>
    <tr>
      <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_FeatureFinderCentroided @n (or another feature finding algorithm) </td>
      <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_FeatureLinkerUnlabeled or @n @ref TOPP_FeatureLinkerUnlabeledQT </td>
    </tr>
  </table>
</CENTER>

This tool provides an algorithm to align the retention time scales of
multiple input files, correcting shifts and distortions between them.
Retention time adjustment may be necessary to correct for chromatography
differences e.g. before data from multiple LC-MS runs can be combined
(feature grouping), or when one run should be annotated with peptide
identifications obtained in a different run.

All map alignment tools (MapAligner...) collect retention time data from the
input files and - by fitting a model to this data - compute transformations
that map all runs to a common retention time scale. They can apply the
transformations right away and return output files with aligned time scales
(parameter @p out), and/or return descriptions of the transformations in
trafoXML format (parameter @p trafo_out). Transformations stored as trafoXML
can be applied to arbitrary files with the @ref TOPP_MapRTTransformer tool.

The map alignment tools differ in how they obtain retention time data for the
modeling of transformations, and consequently what types of data they can be
applied to. The alignment algorithm implemented here is the pose clustering
algorithm as described in doi:10.1093/bioinformatics/btm209. It is used to
find an affine transformation, which is further refined by a feature grouping
step.  This algorithm can be applied to features (featureXML) and peaks
(mzML), but it has mostly been developed and tested on features.  For more
details and algorithm-specific parameters (set in the INI file) see "Detailed
Description" in the @ref OpenMS::MapAlignmentAlgorithmPoseClustering "algorithm documentation".

@see @ref TOPP_MapAlignerPoseClustering @ref TOPP_MapRTTransformer

This algorithm uses an affine transformation model.

To speed up the alignment, consider reducing 'max_number_of_peaks_considered'.
If your alignment is not good enough, consider increasing this number (the alignment will take longer though).

<B>The command line parameters of this tool are:</B> @n
@verbinclude TOPP_MapAlignerPoseClustering.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_MapAlignerPoseClustering.html
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPMapAlignerPoseClustering :
  public TOPPMapAlignerBase
{

public:
  TOPPMapAlignerPoseClustering() :
    TOPPMapAlignerBase("MapAlignerPoseClustering", "Corrects retention time distortions between maps using a pose clustering approach.")
  {}

protected:
  void registerOptionsAndFlags_() override
  {
    TOPPMapAlignerBase::registerOptionsAndFlagsMapAligners_("featureXML,mzML",
                                                            REF_RESTRICTED);
    registerSubsection_("algorithm", "Algorithm parameters section");
  }

  Param getSubsectionDefaults_(const String& section) const override
  {
    if (section == "algorithm")
    {
      MapAlignmentAlgorithmPoseClustering algo;
      return algo.getParameters();
    }
    return Param(); // shouldn't happen
  }

  ExitCodes main_(int, const char**) override
  {
    ExitCodes ret = TOPPMapAlignerBase::checkParameters_();
    if (ret != EXECUTION_OK)
    {
      return ret;
    }
    MapAlignmentAlgorithmPoseClustering algorithm;
    Param algo_params = getParam_().copy("algorithm:", true);
    algorithm.setParameters(algo_params);
    algorithm.setLogType(log_type_);

    StringList in_files = getStringList_("in");
    if (in_files.size() == 1)
    {
      OPENMS_LOG_WARN << "Only one file provided as input to MapAlignerPoseClustering." << std::endl;
    }
    
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
    else if (reference_index > 0) // normal reference (index was checked before)
    {
      file = in_files[--reference_index]; // ref. index is 1-based in parameters, but should be 0-based here
    }
    else if (reference_index == 0) // no reference given
    {
      OPENMS_LOG_INFO << "Picking a reference (by size) ..." << std::flush;
      // use map with highest number of features as reference:
      Size max_count(0);
      FeatureXMLFile f;
      for (Size i = 0; i < in_files.size(); ++i)
      {
        Size s = 0;
        if (in_type == FileTypes::FEATUREXML) 
        {
          s = f.loadSize(in_files[i]);
        }
        else if (in_type == FileTypes::MZML) // this is expensive!
        {
          PeakMap exp;
          FileHandler().loadExperiment(in_files[i], exp, {FileTypes::MZML});
          exp.updateRanges(1);
          s = exp.getSize();
        }
        if (s > max_count)
        {
          max_count = s;
          reference_index = i;
        }
      }
      OPENMS_LOG_INFO << " done" << std::endl;
      file = in_files[reference_index];
    }

    FileHandler f_fxml;
    if (out_files.empty()) // no need to store featureXML, thus we can load only minimum required information
    {
      f_fxml.getFeatOptions().setLoadConvexHull(false);
      f_fxml.getFeatOptions().setLoadSubordinates(false);
    }
    if (in_type == FileTypes::FEATUREXML)
    {
      FeatureMap map_ref;
      FileHandler f_fxml_tmp; // for the reference, we never need CH or subordinates
      f_fxml_tmp.getFeatOptions().setLoadConvexHull(false);
      f_fxml_tmp.getFeatOptions().setLoadSubordinates(false);
      f_fxml_tmp.loadFeatures(file, map_ref, {FileTypes::FEATUREXML});
      algorithm.setReference(map_ref);
    }
    else if (in_type == FileTypes::MZML)
    {
      PeakMap map_ref;
      FileHandler().loadExperiment(file, map_ref);
      algorithm.setReference(map_ref);
    }

    ProgressLogger plog;
    plog.setLogType(log_type_);

    plog.startProgress(0, in_files.size(), "Aligning input maps");
    Size progress(0); // thread-safe progress
    // TODO: it should all work on featureXML files, since we might need them for output anyway. Converting to consensusXML is just wasting memory!
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1)
#endif
    for (int i = 0; i < static_cast<int>(in_files.size()); ++i)
    {
      TransformationDescription trafo;
      if (in_type == FileTypes::FEATUREXML)
      {
        FeatureMap map;
        // workaround for loading: use temporary FeatureXMLFile since it is not thread-safe
        FileHandler f_fxml_tmp; // do not use OMP-firstprivate, since FeatureXMLFile has no copy c'tor
        f_fxml_tmp.getFeatOptions() = f_fxml.getFeatOptions();
        f_fxml_tmp.loadFeatures(in_files[i], map);
        if (i == static_cast<int>(reference_index)) 
        {
          trafo.fitModel("identity");
        }
        else 
        {
          try
          {
            algorithm.align(map, trafo);
          }
          catch (Exception::IllegalArgument& e)
          {
            OPENMS_LOG_ERROR << "Aligning " << in_files[i] << " to reference " << in_files[reference_index]
                             << " failed. No transformation will be applied (RT not changed for this file)." << endl;
            writeLogError_("Illegal argument (" + String(e.getName()) + "): " + String(e.what()) + ".");
            trafo.fitModel("identity");
          }
        }

        if (!out_files.empty())
        {
          MapAlignmentTransformer::transformRetentionTimes(map, trafo);
          // annotate output with data processing info
          addDataProcessing_(map, getProcessingInfo_(DataProcessing::ALIGNMENT));
          f_fxml_tmp.storeFeatures(out_files[i], map, {FileTypes::FEATUREXML});
        }
      }
      else if (in_type == FileTypes::MZML)
      {
        PeakMap map;
        FileHandler().loadExperiment(in_files[i], map, {FileTypes::MZML});
        if (i == static_cast<int>(reference_index))
        {
          trafo.fitModel("identity");
        }
        else
        {
          algorithm.align(map, trafo);
        }
        if (!out_files.empty())
        {
          MapAlignmentTransformer::transformRetentionTimes(map, trafo);
          // annotate output with data processing info
          addDataProcessing_(map, getProcessingInfo_(DataProcessing::ALIGNMENT));
          FileHandler().storeExperiment(out_files[i], map, {FileTypes::MZML});
        }
      }

      if (!out_trafos.empty())
      {
        FileHandler().storeTransformations(out_trafos[i], trafo, {FileTypes::TRANSFORMATIONXML});
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

int main(int argc, const char** argv)
{
  TOPPMapAlignerPoseClustering tool;
  return tool.main(argc, argv);
}

/// @endcond
