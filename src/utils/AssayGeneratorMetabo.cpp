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
#include <OpenMS/ANALYSIS/QUANTITATION/KDTreeFeatureMaps.h>

using namespace OpenMS;
using namespace std;


//-------------------------------------------------------------
//Doxygen docu
//----------------------------------------------------------
/**
  @page UTILS_AssayGeneratorMetabo AssayGeneratorMetabo

  @brief Assay Library generator using DDA Metabolomics Data

  mzml blablabla

  featureXML blabalbla

  // For DDA data -> FeatureFinderMetabo -> AccurateMassSearch -> AssayGeneratorMetabo

  // TODO: DIA data -> FFM -> Unkonwn Features -> AccurateMassSearch would be possible / but maybe better unknown

  <B>The command line parameters of this tool are:</B>
  @verbinclude UTILS_SiriusAdapter.cli
  <B>INI file documentation of this tool:</B>
  @htmlinclude UTILS_SiriusAdapter.html
 */

/// @cond TOPPCLASSES


/// struct to hold assay information of one row 
struct assayrow
{
  int precursor_mz;
  int product_mz;
  int library_int; // relative intensity
  int normalized_rt;
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

  // TODO: rename function
  // extract adduct information from featureXML (MetaboliteAdductDecharger)
  void extractMetaInformation(const PeakMap & spectra, const FeatureMap & feature_map, map< size_t, vector<StringList> > & map_precursor_to_metadata)
  {
    KDTreeFeatureMaps metadata_map_kd;
    vector<FeatureMap> metadata_map;
    metadata_map.push_back(feature_map);
    metadata_map_kd.addMaps(metadata_map);

    // map precursors to closest feature and retrieve annotated metadata
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
        metadata_map_kd.queryRegion(rt - 5.0, rt + 5.0, mz - 0.2, mz + 0.2, matches, true);

        // TODO: regrads precursors with feature but no metadata as unkowns
        // no metadata information found
        if (matches.empty()) { continue; }

        // in the case of multiple features in tolerance window, select the one closest in m/z to the precursor
        Size min_distance_feature_index(0);
        double min_distance(1e11);
        for (auto const & k_idx : matches)
        {
          const double f_mz = metadata_map_kd.mz(k_idx);
          const double distance = fabs(f_mz - mz);
          if (distance < min_distance)
          {
            min_distance = distance;
            min_distance_feature_index = k_idx;
          }
        }
        const BaseFeature * min_distance_feature = metadata_map_kd.feature(min_distance_feature_index);

        // extract metadata from featureXML (AccurateMassSearch) and associate with precursor
        if (min_distance_feature->metaValueExists("description"))
        {
          vector<StringList> metadata;
          StringList description = min_distance_feature->getMetaValue("description");
          metadata.push_back(description);
          StringList sumformula = min_distance_feature->getMetaValue("chemical_formula");
          metadata.push_back(sumformula);
          map_precursor_to_metadata[index] = metadata;
        }
      }
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

    vector<assayrow> assaylib;
 
    MzMLFile f;
    PeakMap exp;
    f.load(in, exp);

    // Read FeatureXML in KDTree for range query
    map<size_t, vector<StringList> > map_precursor_to_metadata;
    std::ifstream id_file(id);
    if (id_file)
    {
      FeatureXMLFile fxml;
      FeatureMap feature_map;
      fxml.load(id, feature_map);
      extractMetaInformation(exp, feature_map, map_precursor_to_metadata);
    }

    for(auto& kv : map_precursor_to_metadata)
    {
      std::cout << kv.first << std::endl;
//      for (auto& t : it.second)
//      {
//        std::cout << t.first << std::endl;
//        std::cout << t.second << std::endl;
//      }
    }

    //    for (MSExperiment::ConstIterator spec_iter = exp.begin(); spec_iter != exp.end(); ++spec_iter)
//    {
//      if(spec_iter->getMSLevel() != 2)
//      {
//        continue;
//      }
//
//      const MSSpectrum& spectrum = *spec_iter;
//      const vector<Precursor>& precursor = spectrum.getPrecursors();
//
//    }

    // TODO: rewrite featureXML information from AccurateMassSearch
    // woher informationen in featureXML (MetaboliteAdductDecharger/AccurateMassSearch)

    // take featureXML from accurate mass search to look at the precursor/MS2
    // extract precursor information from featureXML/tsv
    // get all MS2 spectra of same precursor (ggf. consensusspectrum)
    // Iteratre over peaks in spectra an search for highest one
    // then second and third highest if possible
    // extract the intensity (from consesusspectrum -> transistion)
    // every further Transision which is at least of 10% of highest intensity
    // integrate further information from .tsv for example
   
    //-------------------------------------------------------------
    // Calculations
    //-------------------------------------------------------------


    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------

  }
};

int main(int argc, const char ** argv)
{
  TOPPAssayGeneratorMetabo tool;
  return tool.main(argc, argv);
}

/// @endcond
