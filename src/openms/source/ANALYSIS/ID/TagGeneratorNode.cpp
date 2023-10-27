
//
// Created by trapho on 10/20/23.
//
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

#include <functional>


using namespace std;



namespace OpenMS
{


  TagGeneratorNode::TagGeneratorNode(const OpenMS::Peak1D& peak):
    charge_(1), peak_(peak), norm_intensity_(0), confidence_(0), connected_nodes_(), distance_to_nodes(), connected_AA()
  {
  }

  TagGeneratorNode::TagGeneratorNode(const OpenMS::TagGeneratorNode& cp) = default;

  TagGeneratorNode::TagGeneratorNode(const OpenMS::TagGeneratorNode& cp, uint16_t charge) :
  charge_(charge), peak_(cp.peak_), norm_intensity_(cp.norm_intensity_), confidence_(0)
  {
  }

  TagGeneratorNode::~TagGeneratorNode() = default;

  TagGeneratorNode& TagGeneratorNode::operator=(const OpenMS::TagGeneratorNode& source)
  {
    if(this != &source)
    {
      peak_ = source.peak_;
      charge_ = source.charge_;
      connected_nodes_ = source.connected_nodes_;
    }
    return *this;
  }

  void TagGeneratorNode::generateConnection(const TagGeneratorNode& other, double fragment_tolerance)
  {
    double delta_mass = abs(peak_.getMZ() * charge_ - other.peak_.getMZ() * other.charge_);
    for(size_t mass = 0; mass < monoisotopicMasses.size(); mass++){
      if(monoisotopicMasses[mass] - fragment_tolerance <= delta_mass && delta_mass <= monoisotopicMasses[mass]+ fragment_tolerance){
        this->connected_nodes_.push_back(other);
        this->distance_to_nodes.push_back(monoisotopicMasses[mass]);
        this->connected_AA.push_back(aminoAcids[mass]);
        break;
      }
    }
  }

  void TagGeneratorNode::generateAllQuadPeaksRecursion(vector<QuadPeak>& quad_peaks, QuadPeak quad_peak, uint32_t recursion_step, double prev_mz, string prev_AA)
  {
    quad_peak.sequence+= prev_AA;
    quad_peak.quad_peak_identifier_ += prev_mz;

    quad_peak.score_ += confidence_;

    if(recursion_step == 0){
      cout << quad_peak.score_ << quad_peak.sequence << endl;

      quad_peaks.push_back(quad_peak);

      return;
    }

    for(size_t i = 0; i < connected_nodes_.size(); i++){
      connected_nodes_[i].generateAllQuadPeaksRecursion(quad_peaks, quad_peak, recursion_step - 1, distance_to_nodes[i], connected_AA[i]);
    }
  }

  void TagGeneratorNode::generateAllQuadPeaks(std::vector<TagGeneratorNode::QuadPeak>& quad_peaks, uint8_t recursion_step)
  {

    for(size_t i = 0; i < connected_nodes_.size(); i++){
      cout << peak_.getMZ() << endl;
      QuadPeak qp{confidence_, 0};
      connected_nodes_[i].generateAllQuadPeaksRecursion(quad_peaks, qp, recursion_step - 1, distance_to_nodes[i], connected_AA[i]);
    }

  }

  void TagGeneratorNode::calculateConfidence(const OpenMS::MSSpectrum& spectrum, double max_intensity)
  {
    norm_intensity_ = peak_.getIntensity() / max_intensity;
    confidence_ = norm_intensity_;
  }



  uint16_t TagGeneratorNode::getCharge() const
  {
    return charge_;
  }
  const Peak1D& TagGeneratorNode::getPeak() const
  {
    return peak_;
  }
}
