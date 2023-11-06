
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
#include <OpenMS/DATASTRUCTURES/MultiPeak.h>


#include <functional>


using namespace std;



namespace OpenMS
{

  std::vector<std::string> aminoAcids = {
    "A", "C", "D", "E", "F", "G", "H", "I", "K", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"
  };
  std::vector<double> monoisotopicMasses = {
    71.03711, 103.00919, 115.02694, 129.04259, 147.06841, 57.02146, 137.05891, 113.08406, 128.09496, 131.04049, 114.04293, 97.05276, 128.05858, 156.10111, 87.03203, 101.04768, 99.06841, 186.07932, 163.06332
  };

  TagGeneratorNode::TagGeneratorNode(const OpenMS::Peak1D& peak):
    charge_(1), peak_(peak), norm_intensity_(0), confidence_(0), connected_nodes_({}), distance_to_nodes({}), connected_AA({})
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

  bool TagGeneratorNode::generateConnection(const shared_ptr<TagGeneratorNode>& other, double fragment_tolerance)
  {

    double delta_mass = abs(peak_.getMZ() * charge_ - other->peak_.getMZ() * other->charge_);
    if(delta_mass > 190.0)   // the largest monoisotopic mass
      return false;
    for(size_t mass = 0; mass < monoisotopicMasses.size(); mass++){
      if(monoisotopicMasses[mass] - fragment_tolerance <= delta_mass && delta_mass <= monoisotopicMasses[mass]+ fragment_tolerance){
        this->connected_nodes_.push_back(other);
        this->distance_to_nodes.push_back(delta_mass);  // we want the actual mass, not the theoretical
        this->connected_AA.push_back(aminoAcids[mass]);
        break;
      }
    }
    return true;
  }

  void TagGeneratorNode::generateAllQuadPeaksRecursion(vector<MultiPeak>& quad_peaks, MultiPeak quad_peak, uint32_t recursion_step, double delta_mz, string prev_AA)
  {
    quad_peak.addFollowUpPeak(delta_mz, prev_AA);


    quad_peak.addScore(confidence_);

    if(recursion_step == 0){


      quad_peaks.push_back(quad_peak);

      return;
    }

    for(size_t i = 0; i < connected_nodes_.size(); i++){
      connected_nodes_[i]->generateAllQuadPeaksRecursion(quad_peaks, quad_peak, recursion_step - 1, distance_to_nodes.at(i), connected_AA.at(i));
    }
  }

  void TagGeneratorNode::generateAllMultiPeaks(std::vector<MultiPeak>& quad_peaks, uint8_t recursion_step)
  {

    for(size_t i = 0; i < connected_nodes_.size(); i++){


      MultiPeak qp(peak_, confidence_);
      connected_nodes_[i]->generateAllQuadPeaksRecursion(quad_peaks, qp, recursion_step - 1, distance_to_nodes.at(i), connected_AA.at(i));
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
