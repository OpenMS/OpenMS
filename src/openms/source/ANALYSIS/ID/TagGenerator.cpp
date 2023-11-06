//
// Created by trapho on 10/23/23.
//
#include <OpenMS/ANALYSIS/ID/TagGenerator.h>
#include <OpenMS/ANALYSIS/ID/FragmentIndexTDScorer.h>
#include <OpenMS/ANALYSIS/ID/TagGeneratorNode.h>

#include <OpenMS/CHEMISTRY/AAIndex.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/DigestionEnzyme.h>
#include <OpenMS/CHEMISTRY/EnzymaticDigestion.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/ModifiedPeptideGenerator.h>
#include <OpenMS/CHEMISTRY/ProteaseDB.h>
#include <OpenMS/CHEMISTRY/ProteaseDigestion.h>
#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/DATASTRUCTURES/StringView.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>
#include <OpenMS/QC/QCBase.h>
#include <OpenMS/DATASTRUCTURES/MultiPeak.h>

#include <functional>


using namespace std;


namespace OpenMS
{
  TagGenerator::TagGenerator(const OpenMS::MSSpectrum& spectrum)
  {
    spectrum_ = spectrum;
    selected_peaks_.resize(spectrum_.size());
    n = spectrum_.getPrecursors()[0].getMZ()/ 35;  // TODO: Does this make sense??
  }
  TagGenerator::~TagGenerator() = default;

  void TagGenerator::setMSSpectrum(const MSSpectrum& spectrum)
  {
    spectrum_.clear(true);
    selected_peaks_.clear();
    dag_.clear();

    spectrum_ = spectrum;
    selected_peaks_.resize(spectrum_.size());
    n = spectrum_.getPrecursors()[0].getMZ()/ 35;  // TODO: Does this make sense??
  }


  void TagGenerator::globalSelection()
  {
    vector<IdxAndIntensity> idx_intensity;
    uint32_t idx_count = 0;
    for(Peak1D peak: spectrum_){
      idx_intensity.emplace_back(idx_count, peak.getIntensity());
      idx_count++;
    }
    //we can use the same heap structure as in the scoring to get the top k hits
    FragmentIndexTDScorer fitds;
    fitds.buildKMinHeapforTag(idx_intensity, n);

    for(auto iai = idx_intensity.begin();iai != idx_intensity.end() && iai != idx_intensity.begin() + n; iai++){
      selected_peaks_[iai->idx_] = true;

    }

  }

  void TagGenerator::localSelection()
  {
    size_t lower_peak = 0;
    size_t upper_peak = 0;
    for(double window_da = 0; window_da < spectrum_.getPrecursors()[0].getUnchargedMass(); window_da += 35.0){ //window increment = 35 Da
      while (upper_peak < spectrum_.size() && spectrum_[upper_peak].getMZ() <= (window_da +70)) upper_peak++;   //window size      = 70 Da
      size_t sum = 0; // how many peaks are in the 70-Da window
      vector<IdxAndIntensity> temp_snippet; // collect all idx within the window and sort them afterward
      for(auto sele_iter = selected_peaks_.begin() + lower_peak; sele_iter != selected_peaks_.begin()+upper_peak; sele_iter++){
        sum += static_cast<size_t>(*sele_iter);
        size_t current_peak =sele_iter - selected_peaks_.begin();
        temp_snippet.emplace_back(current_peak, spectrum_[current_peak].getIntensity());
      }
      sort(temp_snippet.begin(), temp_snippet.end(),[](IdxAndIntensity a, IdxAndIntensity b){return (a.intensity_ > b.intensity_);});
      while((sum < 2) && (sum < temp_snippet.size())){                 // now select new peaks if needed for this window
        selected_peaks_[temp_snippet[sum].idx_] = true;
        sum++;
      }

      lower_peak = upper_peak;
    }
  }

  void TagGenerator::generateDirectedAcyclicGraph(double fragment_tolerance)
  {
    //first extract the max intensity to later normalize the intensities
    double max = 0;
    for(Peak1D p: spectrum_){
      if(p.getIntensity() > max) max = p.getIntensity();
    }

    for(size_t i = selected_peaks_.size(); i > 0; i--){
      if(selected_peaks_[i-1]){
        shared_ptr<TagGeneratorNode> newNode = make_shared<TagGeneratorNode>(spectrum_[i-1]);
        newNode->calculateConfidence(spectrum_, max);
        for(auto dag_iter = dag_.end() - 1; dag_iter != dag_.begin()-1 ; dag_iter--){
          if(!newNode->generateConnection(*dag_iter, fragment_tolerance))
            break;
        }
        dag_.push_back(newNode);

      }
    }
  }

  void TagGenerator::generateAllMultiPeaks(std::vector<MultiPeak>& quad_peaks)
  {
    for(shared_ptr<TagGeneratorNode> node: dag_){
      vector<MultiPeak> quad_peaks_per_node;
      node->generateAllMultiPeaks(quad_peaks_per_node, 3);
      quad_peaks.insert(quad_peaks.end(), quad_peaks_per_node.begin(), quad_peaks_per_node.end());
    }
  }



}
