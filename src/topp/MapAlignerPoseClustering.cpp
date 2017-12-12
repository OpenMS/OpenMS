// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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

  @see @ref TOPP_MapAlignerPoseClustering @ref TOPP_MapAlignerSpectrum @ref TOPP_MapRTTransformer

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
  void registerOptionsAndFlags_()
  {
    TOPPMapAlignerBase::registerOptionsAndFlags_("featureXML,mzML",
                                                 REF_RESTRICTED);
    registerSubsection_("algorithm", "Algorithm parameters section");
  }

  Param getSubsectionDefaults_(const String& section) const
  {
    if (section == "algorithm")
    {
      MapAlignmentAlgorithmPoseClustering algo;
      return algo.getParameters();
    }
    return Param(); // shouldn't happen
  }

  ExitCodes main_(int, const char**)
  {
    ExitCodes ret = TOPPMapAlignerBase::checkParameters_();
    if (ret != EXECUTION_OK) return ret;

    MapAlignmentAlgorithmPoseClustering algorithm;
    Param algo_params = getParam_().copy("algorithm:", true);
    algorithm.setParameters(algo_params);
    algorithm.setLogType(log_type_);

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
    else if (reference_index > 0) // normal reference (index was checked before)
    {
      file = in_files[--reference_index]; // ref. index is 1-based in parameters, but should be 0-based here
    }
    else if (reference_index == 0) // no reference given
    {
      LOG_INFO << "Picking a reference (by size) ..." << std::flush;
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
          MzMLFile().load(in_files[i], exp);
          exp.updateRanges(1);
          s = exp.getSize();
        }
        if (s > max_count)
        {
          max_count = s;
          reference_index = i;
        }
      }
      LOG_INFO << " done" << std::endl;
      file = in_files[reference_index];
    }

    FeatureXMLFile f_fxml;
    if (out_files.empty()) // no need to store featureXML, thus we can load only minimum required information
    {
      f_fxml.getOptions().setLoadConvexHull(false);
      f_fxml.getOptions().setLoadSubordinates(false);
    }
    if (in_type == FileTypes::FEATUREXML)
    {
      FeatureMap map_ref;
      FeatureXMLFile f_fxml_tmp; // for the reference, we never need CH or subordinates
      f_fxml_tmp.getOptions().setLoadConvexHull(false);
      f_fxml_tmp.getOptions().setLoadSubordinates(false);
      f_fxml_tmp.load(file, map_ref);
      algorithm.setReference(map_ref);
    }
    else if (in_type == FileTypes::MZML)
    {
      PeakMap map_ref;
      MzMLFile().load(file, map_ref);
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
        FeatureXMLFile f_fxml_tmp; // do not use OMP-firstprivate, since FeatureXMLFile has no copy c'tor
        f_fxml_tmp.getOptions() = f_fxml.getOptions();
        f_fxml_tmp.load(in_files[i], map);
        if (i == static_cast<int>(reference_index)) trafo.fitModel("identity");
        else algorithm.align(map, trafo);
        if (out_files.size())
        {
          MapAlignmentTransformer::transformRetentionTimes(map, trafo);
          // annotate output with data processing info
          addDataProcessing_(map, getProcessingInfo_(DataProcessing::ALIGNMENT));
          f_fxml_tmp.store(out_files[i], map);
        }
      }
      else if (in_type == FileTypes::MZML)
      {
        PeakMap map;
        MzMLFile().load(in_files[i], map);
        if (i == static_cast<int>(reference_index)) trafo.fitModel("identity");
        else algorithm.align(map, trafo);
        if (out_files.size())
        {
          MapAlignmentTransformer::transformRetentionTimes(map, trafo);
          // annotate output with data processing info
          addDataProcessing_(map, getProcessingInfo_(DataProcessing::ALIGNMENT));
          MzMLFile().store(out_files[i], map);
        }
      }

      if (!out_trafos.empty())
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

int main(int argc, const char** argv)
{
  TOPPMapAlignerPoseClustering tool;
  return tool.main(argc, argv);
}

/// @endcond
