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
// $Maintainer: Hendrik Weisser $
// $Authors: Marc Sturm, Clemens Groepl, Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmIdentification.h>
#include <OpenMS/APPLICATIONS/MapAlignerBase.h>
#include <OpenMS/METADATA/ExperimentalDesign.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------

/**
    @page TOPP_MapAlignerIdentification MapAlignerIdentification

        @brief Corrects retention time distortions between maps, using information from peptides identified in different maps.

<CENTER>
    <table>
        <tr>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> potential predecessor tools </td>
            <td VALIGN="middle" ROWSPAN=4> \f$ \longrightarrow \f$ MapAlignerIdentification \f$ \longrightarrow \f$</td>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> potential successor tools </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_XTandemAdapter @n (or another search engine adapter) </td>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_IDMerger </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_IDFileConverter </td>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=2> @ref TOPP_FeatureLinkerUnlabeled or @n @ref TOPP_FeatureLinkerUnlabeledQT </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_IDMapper </td>
        </tr>
    </table>
</CENTER>

    Reference:\n
		Weisser <em>et al.</em>: <a href="http://dx.doi.org/10.1021/pr300992u">An automated pipeline for high-throughput label-free quantitative proteomics</a> (J. Proteome Res., 2013, PMID: 23391308).

    This tool provides an algorithm to align the retention time scales of multiple input files, correcting shifts and distortions between them. Retention time adjustment may be necessary to correct for chromatography differences e.g. before data from multiple LC-MS runs can be combined (feature grouping), or when one run should be annotated with peptide identifications obtained in a different run.

    All map alignment tools (MapAligner...) collect retention time data from the input files and - by fitting a model to this data - compute transformations that map all runs to a common retention time scale. They can apply the transformations right away and return output files with aligned time scales (parameter @p out), and/or return descriptions of the transformations in trafoXML format (parameter @p trafo_out). Transformations stored as trafoXML can be applied to arbitrary files with the @ref TOPP_MapRTTransformer tool.

    The map alignment tools differ in how they obtain retention time data for the modeling of transformations, and consequently what types of data they can be applied to. The alignment algorithm implemented here is based on peptide identifications, and thus applicable to files containing peptide IDs (idXML, annotated featureXML/consensusXML). It finds peptide sequences that different input files have in common and uses them as points of correspondence between the inputs. For more details and algorithm-specific parameters (set in the INI file) see "Detailed Description" in the @ref OpenMS::MapAlignmentAlgorithmIdentification "algorithm documentation".

    @see @ref TOPP_MapAlignerPoseClustering @ref TOPP_MapAlignerSpectrum @ref TOPP_MapRTTransformer

		Note that alignment is based on the sequence including modifications, thus an exact match is required. I.e., a peptide with oxidised methionine will not be matched to its unmodified version. This behavior is generally desired since (some) modifications can cause retention time shifts.

    Since %OpenMS 1.8, the extraction of data for the alignment has been separate from the modeling of RT transformations based on that data. It is now possible to use different models independently of the chosen algorithm. This algorithm has been tested mostly with the "b_spline" model. The different available models are:
    - @ref OpenMS::TransformationModelLinear "linear": Linear model.
    - @ref OpenMS::TransformationModelBSpline "b_spline": Smoothing spline (non-linear).
    - @ref OpenMS::TransformationModelLowess "lowess": Local regression (non-linear).
    - @ref OpenMS::TransformationModelInterpolated "interpolated": Different types of interpolation.

    The following parameters control the modeling of RT transformations (they can be set in the "model" section of the INI file):
    @htmlinclude OpenMS_MapAlignerIdentificationModel.parameters @n

    @note Currently mzIdentML (mzid) is not directly supported as an input/output format of this tool. Convert mzid files to/from idXML using @ref TOPP_IDFileConverter if necessary.

    <B>The command line parameters of this tool are:</B> @n
    @verbinclude TOPP_MapAlignerIdentification.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude TOPP_MapAlignerIdentification.html

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPMapAlignerIdentification :
  public TOPPMapAlignerBase
{

public:
  TOPPMapAlignerIdentification() :
    TOPPMapAlignerBase("MapAlignerIdentification", "Corrects retention time distortions between maps based on common peptide identifications.")
  {
  }

private:
  template <typename MapType, typename FileType>
  void loadInitialMaps_(vector<MapType>& maps, StringList& ins, 
                        FileType& input_file)
  {
    // custom progress logger for this task:
    ProgressLogger progresslogger;
    progresslogger.setLogType(TOPPMapAlignerBase::log_type_);
    progresslogger.startProgress(0, ins.size(), "loading input files");
    for (Size i = 0; i < ins.size(); ++i)
    {
      progresslogger.setProgress(i);
      input_file.load(ins[i], maps[i]);
    }
    progresslogger.endProgress();
  }

  // helper function to avoid code duplication between consensusXML and
  // featureXML storage operations:
  template <typename MapType, typename FileType>
  void storeTransformedMaps_(vector<MapType>& maps, StringList& outs, 
                             FileType& output_file)
  {
    // custom progress logger for this task:
    ProgressLogger progresslogger;
    progresslogger.setLogType(log_type_);
    progresslogger.startProgress(0, outs.size(), "writing output files");
    for (Size i = 0; i < outs.size(); ++i)
    {
      progresslogger.setProgress(i);
      // annotate output with data processing info:
      addDataProcessing_(maps[i], 
                         getProcessingInfo_(DataProcessing::ALIGNMENT));
      output_file.store(outs[i], maps[i]);
    }
    progresslogger.endProgress();
  }

  template <typename DataType>
  void performAlignment_(MapAlignmentAlgorithmIdentification& algorithm,
                         vector<DataType>& data,
                         vector<TransformationDescription>& transformations,
                         Int reference_index)
  {
    algorithm.align(data, transformations, reference_index);

    // find model parameters:
    Param model_params = getParam_().copy("model:", true);
    String model_type = model_params.getValue("type");
    if (model_type != "none")
    {
      model_params = model_params.copy(model_type + ":", true);
      for (vector<TransformationDescription>::iterator it =
             transformations.begin(); it != transformations.end(); ++it)
      {
        it->fitModel(model_type, model_params);
      }
    }
  }

  template <typename DataType>
  void applyTransformations_(vector<DataType>& data,
    const vector<TransformationDescription>& transformations)
  {
    for (Size i = 0; i < data.size(); ++i)
    {
      MapAlignmentTransformer::transformRetentionTimes(data[i],
        transformations[i]);
    }
  }

  void storeTransformationDescriptions_(const vector<TransformationDescription>&
                                        transformations, StringList& trafos)
  {
    // custom progress logger for this task:
    ProgressLogger progresslogger;
    progresslogger.setLogType(log_type_);
    progresslogger.startProgress(0, trafos.size(), 
                                 "writing transformation files");
    for (Size i = 0; i < transformations.size(); ++i)
    {
      TransformationXMLFile().store(trafos[i], transformations[i]);
    }
    progresslogger.endProgress();
  }

  Int getReference_(MapAlignmentAlgorithmIdentification& algorithm)
  {
    // consistency of reference parameters has already been checked via
    // "TOPPMapAlignerBase::checkParameters_"

    Size reference_index = getIntOption_("reference:index");
    String reference_file = getStringOption_("reference:file");

    if (!reference_file.empty())
    {
      FileTypes::Type filetype = FileHandler::getType(reference_file);
      if (filetype == FileTypes::MZML)
      {
        PeakMap experiment;
        MzMLFile().load(reference_file, experiment);
        algorithm.setReference(experiment);
      }
      else if (filetype == FileTypes::FEATUREXML)
      {
        FeatureMap features;
        FeatureXMLFile().load(reference_file, features);
        algorithm.setReference(features);
      }
      else if (filetype == FileTypes::CONSENSUSXML)
      {
        ConsensusMap consensus;
        ConsensusXMLFile().load(reference_file, consensus);
        algorithm.setReference(consensus);
      }
      else if (filetype == FileTypes::IDXML)
      {
        vector<ProteinIdentification> proteins;
        vector<PeptideIdentification> peptides;
        IdXMLFile().load(reference_file, proteins, peptides);
        algorithm.setReference(peptides);
      }
    }

    return Int(reference_index) - 1; // internally, we count from zero
  }

  void registerOptionsAndFlags_() override
  {
    String formats = "featureXML,consensusXML,idXML";
    TOPPMapAlignerBase::registerOptionsAndFlags_(formats, REF_FLEXIBLE);
    // TODO: potentially move to base class so every aligner has to support design
    registerInputFile_("design", "<file>", "", "input file containing the experimental design", false);
    setValidFormats_("design", ListUtils::create<String>("tsv"));
    registerSubsection_("algorithm", "Algorithm parameters section");
    registerSubsection_("model", "Options to control the modeling of retention time transformations from data");
  }

  Param getSubsectionDefaults_(const String & section) const override
  {
    if (section == "algorithm")
    {
      MapAlignmentAlgorithmIdentification algo;
      return algo.getParameters();
    }
    if (section == "model")
    {
      return TOPPMapAlignerBase::getModelDefaults("b_spline");
    }

    return Param(); // this shouldn't happen
  }

  ExitCodes main_(int, const char**) override
  {
    ExitCodes return_code = TOPPMapAlignerBase::checkParameters_();
    if (return_code != EXECUTION_OK) return return_code;

    // set up alignment algorithm:
    MapAlignmentAlgorithmIdentification algorithm;
    Param algo_params = getParam_().copy("algorithm:", true);
    algorithm.setParameters(algo_params);
    algorithm.setLogType(log_type_);

    Int reference_index = getReference_(algorithm);

    // handle in- and output files:
    StringList input_files = getStringList_("in");
    StringList output_files = getStringList_("out");
    StringList trafo_files = getStringList_("trafo_out");
    FileTypes::Type in_type = FileHandler::getType(input_files[0]);

    vector<TransformationDescription> transformations;

    //-------------------------------------------------------------
    // perform feature alignment
    //-------------------------------------------------------------
    if (in_type == FileTypes::FEATUREXML)
    {
      vector<FeatureMap> feature_maps(input_files.size());
      FeatureXMLFile fxml_file;
      if (output_files.empty())
      {
        // store only transformation descriptions, not transformed data =>
        // we can load only minimum required information:
        fxml_file.getOptions().setLoadConvexHull(false);
        fxml_file.getOptions().setLoadSubordinates(false);
      }
      loadInitialMaps_(feature_maps, input_files, fxml_file);

      //-------------------------------------------------------------
      // Extract (optional) fraction identifiers and associate with featureXMLs
      //-------------------------------------------------------------
      String design_file = getStringOption_("design");

      // determine map of fractions to runs
      map<unsigned, set<unsigned> > frac2run;

      // TODO: check if can be put in common helper function
      if (!design_file.empty())
      {
        // parse design file and determine fractions
        ExperimentalDesign ed;
        ExperimentalDesign().load(design_file, ed);

        // determine if design defines more than one fraction (note: fraction and run IDs are one-based)
        frac2run = ed.getFractionToRunsMapping();

        // check if all fractions have the same number of MS runs associated
        if (!ed.sameNrOfRunsPerFraction())
        {
          writeLog_("Error: Number of runs must match for every fraction!");
          return ILLEGAL_PARAMETERS;
        }
      }
      else // no design file given
      {
        for (Size i = 0; i != input_files.size(); ++i)
        {
          frac2run[1].insert(i + 1); // associate each run with fraction 1
        }
      }

      // TODO: check and handle if featureXML order differs from run order

      // perform fraction-based alignment
      if (frac2run.size() == 1) // group one fraction
      {
        performAlignment_(algorithm, feature_maps, transformations,
          reference_index);
        applyTransformations_(feature_maps, transformations);
      }
      else // group multiple fractions
      {
        for (Size i = 1; i <= frac2run.size(); ++i)
        {
          vector<FeatureMap> fraction_maps;
          vector<TransformationDescription> fraction_transformations;
          for (set<unsigned>::const_iterator sit = frac2run[i].begin(); sit != frac2run[i].end(); ++sit)
          {
            fraction_maps.push_back(feature_maps[*sit - 1]); // TODO: *sit is currently the run identifier but we need to know the corresponding feature index in ins
          }
          performAlignment_(algorithm, fraction_maps, fraction_transformations,
            reference_index);
          applyTransformations_(fraction_maps, fraction_transformations);

          // copy into transformations and feature maps
          transformations.insert(transformations.end(), fraction_transformations.begin(), fraction_transformations.end());
          Size f = 0;
          for (set<unsigned>::const_iterator sit = frac2run[i].begin(); sit != frac2run[i].end(); ++sit, ++f)
          {
            feature_maps[*sit - 1].swap(fraction_maps[f]); // TODO: *sit is currently the run identifier but we need to know the corresponding feature index in ins
          }
        }
      }

      if (!output_files.empty())
      {
        storeTransformedMaps_(feature_maps, output_files, fxml_file);
      }
    }

    //-------------------------------------------------------------
    // perform consensus alignment
    //-------------------------------------------------------------
    else if (in_type == FileTypes::CONSENSUSXML)
    {
      std::vector<ConsensusMap> consensus_maps(input_files.size());
      ConsensusXMLFile cxml_file;
      loadInitialMaps_(consensus_maps, input_files, cxml_file);

      performAlignment_(algorithm, consensus_maps, transformations,
                        reference_index);
      applyTransformations_(consensus_maps, transformations);

      if (!output_files.empty())
      {
        storeTransformedMaps_(consensus_maps, output_files, cxml_file);
      }
    }

    //-------------------------------------------------------------
    // perform peptide alignment
    //-------------------------------------------------------------
    else if (in_type == FileTypes::IDXML)
    {
      vector<vector<ProteinIdentification> > protein_ids(input_files.size());
      vector<vector<PeptideIdentification> > peptide_ids(input_files.size());
      IdXMLFile idxml_file;
      ProgressLogger progresslogger;
      progresslogger.setLogType(log_type_);
      progresslogger.startProgress(0, input_files.size(),
                                   "loading input files");
      for (Size i = 0; i < input_files.size(); ++i)
      {
        progresslogger.setProgress(i);
        idxml_file.load(input_files[i], protein_ids[i], peptide_ids[i]);
      }
      progresslogger.endProgress();

      performAlignment_(algorithm, peptide_ids, transformations,
                        reference_index);
      applyTransformations_(peptide_ids, transformations);

      if (!output_files.empty())
      {
        progresslogger.startProgress(0, output_files.size(), 
                                     "writing output files");
        for (Size i = 0; i < output_files.size(); ++i)
        {
          progresslogger.setProgress(i);
          idxml_file.store(output_files[i], protein_ids[i], peptide_ids[i]);
        }
        progresslogger.endProgress();
      }
    }

    if (!trafo_files.empty())
    {
      storeTransformationDescriptions_(transformations, trafo_files);
    }

    return EXECUTION_OK;
  }

};

int main(int argc, const char** argv)
{
  TOPPMapAlignerIdentification tool;
  return tool.main(argc, argv);
}

/// @endcond
