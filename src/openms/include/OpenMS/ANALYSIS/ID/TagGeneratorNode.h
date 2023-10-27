//
// Created by trapho on 10/20/23.
//

#pragma once

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/Peak1D.h>

#include <vector>
#include <functional>



namespace OpenMS
{


  class OPENMS_DLLAPI TagGeneratorNode
  {
  private:
    uint16_t charge_;
    Peak1D peak_;
    double norm_intensity_;
    double confidence_;
    std::vector<TagGeneratorNode> connected_nodes_; // all peaks to the right, which have a mz distance of an AA-mass
    std::vector<double> distance_to_nodes;
    std::vector<std::string> connected_AA; //for debugging only


  public:

    std::map<std::string , double> aaMap = {
      {"A", 71.03711},   // Alanine
      {"C", 103.00919},  // Cysteine
      {"D", 115.02694},  // Aspartic Acid
      {"E", 129.04259},  // Glutamic Acid
      {"F", 147.06841},  // Phenylalanine
      {"G", 57.02146},   // Glycine
      {"H", 137.05891},  // Histidine
      {"I", 113.08406},  // Isoleucine (same as Leucine)
      {"K", 128.09496},  // Lysine
      {"M", 131.04049},  // Methionine
      {"N", 114.04293},  // Asparagine
      {"P", 97.05276},   // Proline
      {"Q", 128.05858},  // Glutamine
      {"R", 156.10111},  // Arginine
      {"S", 87.03203},   // Serine
      {"T", 101.04768},  // Threonine
      {"V", 99.06841},   // Valine
      {"W", 186.07931},  // Tryptophan
      {"Y", 163.06333}   // Tyrosine

    };

    std::vector<std::string> aminoAcids = {
      "A", "C", "D", "E", "F", "G", "H", "I", "K", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"
    };
    std::vector<double> monoisotopicMasses = {
      71.03711, 103.00919, 115.02694, 129.04259, 147.06841, 57.02146, 137.05891, 113.08406, 128.09496, 131.04049, 114.04293, 97.05276, 128.05858, 156.10111, 87.03203, 101.04768, 99.06841, 186.07932, 163.06332
    };

    /// Constructor
    TagGeneratorNode(const Peak1D &peak);

    /// copy consturctor
    TagGeneratorNode(const TagGeneratorNode& cp);

    /// copy a charged version
    TagGeneratorNode(const TagGeneratorNode& cp, uint16_t charge);

    /// assignment operator
    TagGeneratorNode& operator=(const TagGeneratorNode& source);

    /// destructor
    ~TagGeneratorNode();


    /// getter
    uint16_t getCharge() const;
    const Peak1D& getPeak() const;


    struct QuadPeak{
      double score_; //sum of intensities of this peak + of supporting peaks
      Peak1D peak_;
      double quad_peak_identifier_; // distance of first two peak rounded up + 2 to 3 normal + 3 to 4 rounded down
      std::string sequence; // for debuggin only

      QuadPeak(double score, double quad_peak_identifier){
        score_ = score;
        quad_peak_identifier_ = quad_peak_identifier;
        sequence = "";
      }
    };

    /** @brief gets one of the peaks TO THE RIGHT and checks if the distance is in the AA-mass range
     *
     *
     * @param other: a peak to the right !
     */
    void generateConnection(const TagGeneratorNode& other, double fragment_tolerance);

    void generateAllQuadPeaks(std::vector<TagGeneratorNode::QuadPeak>& quad_peaks, uint8_t recursion_step);

    void generateAllQuadPeaksRecursion(std::vector<TagGeneratorNode::QuadPeak>& quad_peaks, TagGeneratorNode::QuadPeak quad_peak, uint32_t recursion_step, double prev_mz, std::string prev_AA);

    void calculateConfidence(const MSSpectrum& spectrum, double max_intensity);


  };
}
