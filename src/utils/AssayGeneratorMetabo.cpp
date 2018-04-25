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
// $Maintainer: Oliver Alka $
// $Authors: Oliver Alka $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/ANALYSIS/QUANTITATION/KDTreeFeatureMaps.h>
#include <vector>
#include <OpenMS/FORMAT/TextFile.h>

using namespace OpenMS;
using namespace std;


//-------------------------------------------------------------
//Doxygen docu
//----------------------------------------------------------
/**
  @page UTILS_AssayGeneratorMetabo AssayGeneratorMetabo

  @brief Assay Library generator using metabolomics DDA Data

  mzml blablabla

  featureXML blabalbla

  // For DDA data -> FeatureFinderMetabo -> AccurateMassSearch -> AssayGeneratorMetabo
  // Right now only works with DDA data (later for DIA data?)
      // take featureXML from accurate mass search get the feature information
    // get all MS2 spectra of same precursor (make a consensusspectrum(? - scoring))
    // reduce to one MS2
        // highest intenstiy
        // total ion count ms2
        // consensus_spectrum
    // iterate over peaks in spectra an search for highest one
    // then second and third highest if possible
    // extract the intensity (from consesusspectrum -> transistion)
    // every further Transision which is at least of 10% of highest intensity
    // integrate further information from .tsv for exampl
    // extract feature information (precursor(?)

  <B>The command line parameters of this tool are:</B>
  @verbinclude UTILS_SiriusAdapter.cli
  <B>INI file documentation of this tool:</B>
  @htmlinclude UTILS_SiriusAdapter.html
 */

/// @cond TOPPCLASSES


/// struct to hold assay information of one row 
struct AssayRow
{
  double precursor_mz;
  double product_mz;
  float library_int;
  double normalized_rt;
  String compound_name;
  String smiles;
  String sumformula;
  String transition_group_id;
  String transition_id;
  bool decoy; 
};

class TOPPAssayGeneratorMetabo :
  public TOPPBase
{
public:
  TOPPAssayGeneratorMetabo() :
    TOPPBase("AssayGeneratorMetabo", "Tool for assay library generation from metabolomics DDA data", false)
    {}

protected:

  void registerOptionsAndFlags_() override
  {

    registerInputFile_("in", "<file>", "", "MzML Input file");
    setValidFormats_("in", ListUtils::create<String>("mzml"));

    registerInputFile_("in_id", "<file>", "", "FeatureXML Input with id information (accurate mass search)", false);
    setValidFormats_("in_id", ListUtils::create<String>("featurexml"));

    registerOutputFile_("out", "<file>", "", "Assay library output file");
    setValidFormats_("out", ListUtils::create<String>("tsv"));
  }

  // map precursors to closest feature and retrieve annotated metadata (if possible)
  // extract meta information from featureXML (MetaboliteAdductDecharger)
  map<const BaseFeature*, std::vector<size_t>> extractMetaInformation(const PeakMap & spectra, const FeatureMap & feature_map)
  {
    map<const BaseFeature*, vector<size_t>> map_feature_to_precursors;
    KDTreeFeatureMaps fp_map_kd;
    vector<FeatureMap> fp_map;
    fp_map.push_back(feature_map);
    fp_map_kd.addMaps(fp_map);

    // map precursors to closest feature and retrieve annotated metadata (if possible)
    for (size_t index = 0; index != spectra.size(); ++index)
    {
      if (spectra[index].getMSLevel() != 2) { continue; }

      // get precursor meta data (m/z, rt)
      const vector<Precursor> & pcs = spectra[index].getPrecursors();

      if (!pcs.empty())
      {
        const double mz = pcs[0].getMZ();
        const double rt = spectra[index].getRT();

        // query features in tolerance window
        vector<Size> matches;

        // TODO: optimize query mz
        fp_map_kd.queryRegion(rt - 5.0, rt + 5.0, mz - 0.2, mz + 0.2, matches, true);

        // no precuros matches the feature information found
        if (matches.empty()) { continue; }

        // in the case of multiple features in tolerance window, select the one closest in m/z to the precursor
        Size min_distance_feature_index(0);
        double min_distance(1e11);
        for (auto const & k_idx : matches)
        {
          const double f_mz = fp_map_kd.mz(k_idx);
          const double distance = fabs(f_mz - mz);
          if (distance < min_distance)
          {
            min_distance = distance;
            min_distance_feature_index = k_idx;
          }
        }
        const BaseFeature* min_distance_feature = fp_map_kd.feature(min_distance_feature_index);

        map_feature_to_precursors[min_distance_feature].push_back(index);
      }
    }
    return map_feature_to_precursors;
  }

  // precursor correction (highest intensity)
  Int getHighestIntensityPeakInMZRange(double test_mz, const MSSpectrum& spectrum1, double left_tolerance, double right_tolerance)
  {
    MSSpectrum::ConstIterator left = spectrum1.MZBegin(test_mz - left_tolerance);
    MSSpectrum::ConstIterator right = spectrum1.MZEnd(test_mz + right_tolerance);

    // no MS1 precursor peak in +- tolerance window found
    if (left == right || left->getMZ() > test_mz + right_tolerance)
    {
      return -1;
    }

    MSSpectrum::ConstIterator max_intensity_it = std::max_element(left, right, Peak1D::IntensityLess());

    if (max_intensity_it == right)
    {
      return -1;
    }

    return max_intensity_it - spectrum1.begin();
  }

  // annotate precursor intenstiy based on precursor spectrum and highest intensity peak in tolerance window
  void annotatePrecursorIntensity(PeakMap& spectra, double left_tolerance, double right_tolerance)
  {
    for (PeakMap::Iterator s_it = spectra.begin(); s_it != spectra.end(); ++s_it)
    {
      // process only MS2 spectra
      if (s_it->getMSLevel() != 2)
      {
        continue;
      }

      MSSpectrum& spectrum = *s_it;
      vector<Precursor>& precursor = spectrum.getPrecursors();
      //TODO: Throw exception if precursor vector is empty

      double test_mz = precursor[0].getMZ();

      // find corresponding precursor specturm
      PeakMap::ConstIterator s_it2 = spectra.getPrecursorSpectrum(s_it);

      // no precursor spectrum found
      if (s_it2 == spectra.end())
      {
        LOG_WARN << "No MS1 Spectrum was found to the specific precursor" << std::endl;
        continue;
      }

      const MSSpectrum& precursor_spectrum = *s_it2;
      Int mono_index = getHighestIntensityPeakInMZRange(test_mz, precursor_spectrum, left_tolerance, right_tolerance);
      if (mono_index == -1)
      {
        LOG_WARN << "No precusor peak in MS1 spectrum found." << std::endl;
        continue;
      }
      const Peak1D& max_mono_peak = precursor_spectrum[mono_index];
      precursor[0].setMZ(max_mono_peak.getMZ());
      precursor[0].setIntensity(max_mono_peak.getIntensity());
    }
  }

  ExitCodes main_(int, const char **) override
  {
    //-------------------------------------------------------------
    // Parsing parameters
    //-------------------------------------------------------------

    String in = getStringOption_("in");
    String id = getStringOption_("in_id");
    String out = getStringOption_("out");
 
    MzMLFile f;
    PeakMap spectra;
    f.load(in, spectra);

    // reannotate precursor mz and intensity
    annotatePrecursorIntensity(spectra, 0.2, 0.2);

    // read FeatureXML in KDTree for range query
    map<const BaseFeature*, std::vector<size_t> > map_feature_to_precursors;
    std::ifstream id_file(id);
    if (id_file)
    {
      FeatureXMLFile fxml;
      FeatureMap feature_map;
      fxml.load(id, feature_map);
      map_feature_to_precursors = extractMetaInformation(spectra, feature_map);
    }

    std::vector<AssayRow> assaylib;

    int transition_group_counter = 0;
    for (std::map<const BaseFeature*, std::vector<size_t>>::iterator it = map_feature_to_precursors.begin();
         it != map_feature_to_precursors.end();
         ++it)
    {

      AssayRow row;
      String description = "";
      String sumformula = "";

      const BaseFeature* min_distance_feature = it->first;

      std::cout << "des: " << min_distance_feature->getIntensity() << std::endl;

      // extract metadata from featureXML
      if (!min_distance_feature->getPeptideIdentifications().empty() &&
          !min_distance_feature->getPeptideIdentifications()[0].getHits().empty())
      {
        description = min_distance_feature->getPeptideIdentifications()[0].getHits()[0].getMetaValue("description");
        sumformula = min_distance_feature->getPeptideIdentifications()[0].getHits()[0].getMetaValue("chemical_formula");
      }
      else
      {
        description = "UNKOWN";
        sumformula = "UNKOWN";
      }

      double spectrum_rt = 0.0;
      double highest_precursor_mz = 0.0;
      float highest_precursor_int = 0.0;
      MSSpectrum highest_precursor_int_spectrum;

      // find precursor/spectrum with highest intensity precursor
      std::vector<size_t> index = it->second;
      for (std::vector<size_t>::iterator index_it = index.begin(); index_it != index.end(); ++index_it)
      {
        const MSSpectrum &spectrum = spectra[*index_it];
        const vector<Precursor> &precursor = spectrum.getPrecursors();

        // get m/z and intensity of precursor
        double precursor_mz = precursor[0].getMZ();
        float precursor_int = precursor[0].getIntensity();

        // TODO: add other methods
        // e.g. pearson correlation and consensusspectrum

        // spectrum with highest intensity precursor
        if (precursor_int > highest_precursor_int)
        {
          highest_precursor_int = precursor_int;
          highest_precursor_mz = precursor_mz;
          highest_precursor_int_spectrum = spectra[*index_it];
          spectrum_rt = spectra[*index_it].getRT();
        }
      }
      // extract transistions by  iterating over peaks in in MS2 and get one with highest int
      // TODO: ggf. sort by intensity (?)
      if (!highest_precursor_int_spectrum.isSorted())
      {
        highest_precursor_int_spectrum.sortByIntensity();
      }

      // calculate max intensity peak and threshold
      float max_int = 0.0;
      for(MSSpectrum::iterator spec_it = highest_precursor_int_spectrum.begin(); spec_it != highest_precursor_int_spectrum.end(); ++spec_it)
      {
        //find the max intensity peak
        if(spec_it->getIntensity() > max_int)
        {
          max_int = spec_it->getIntensity();
        }
      }
      // threshold should be at 25 % of the maximum intensity
      float threshold = max_int * 0.25;

      int transition_counter = 1;
      for(MSSpectrum::iterator spec_it = highest_precursor_int_spectrum.begin(); spec_it != highest_precursor_int_spectrum.end(); ++spec_it)
      {
        float current_int = spec_it->getIntensity();
        double current_mz = spec_it->getMZ();

        //write row for each transistion
        if (current_int > threshold)
        {
          float rel_int = current_int/max_int;
          row.precursor_mz = highest_precursor_mz;
          row.product_mz = current_mz;
          row.library_int = rel_int;
          row.normalized_rt = spectrum_rt;
          row.compound_name = description;
          //TODO: add simles to AccurateMassSearchMetadata
          row.smiles = "none";
          row.sumformula = sumformula;
          row.transition_group_id = String(transition_group_counter)+"_"+sumformula;
          row.transition_id = String(transition_counter)+"_"+sumformula;
          row.decoy = 0;
          transition_counter += 1;
        }
        else
        {
          continue;
        }
      assaylib.push_back(row);
      }
       transition_group_counter += 1;
    }

    // write output
    TextFile tf;
    tf.addLine("PrecursorMz\tProductMz\tLibraryIntensity\tRetentionTime\tCompoundName\tSMILES\tSumFormula\tTransitionId\tTransitionGroupId\tDecoy\n");
    for (auto & entry : assaylib)
    {
      tf.addLine(String(entry.precursor_mz)+"\t"+String(entry.product_mz)+"\t"+String(entry.library_int)+"\t"+String(entry.normalized_rt)+
                 "\t"+entry.compound_name+"\t"+entry.smiles+"\t"+entry.sumformula+"\t"+entry.transition_group_id+"\t"+entry.transition_id+
                 "\t"+entry.decoy+"\n");
    }
    tf.store(out);

    return EXECUTION_OK;
  }
};

int main(int argc, const char ** argv)
{
  TOPPAssayGeneratorMetabo tool;
  return tool.main(argc, argv);
}

/// @endcond
