//
// Created by trapho on 10/20/23.
//

#pragma once

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/DATASTRUCTURES/MultiPeak.h>

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
    bool include_;
    std::vector<std::shared_ptr<TagGeneratorNode>> connected_nodes_; // all peaks to the right, which have a mz distance of an AA-mass
    std::vector<double> distance_to_nodes;  // This is what we actually want
    std::vector<std::string> connected_AA; //for debugging only //TODO: can later be removed


  public:



    /// Constructor
    TagGeneratorNode(const Peak1D& peak);

    /// Constructor with charge
    TagGeneratorNode(const Peak1D& peak, uint32_t charge, bool include);

    /// copy consturctor
    TagGeneratorNode(const TagGeneratorNode& cp);

    /// copy a charged version
    TagGeneratorNode(const TagGeneratorNode& cp, uint16_t charge);

    /// assignment operator
    TagGeneratorNode& operator=(const TagGeneratorNode& source);

    /// destructor
    ~TagGeneratorNode();


    /// getter
    [[nodiscard]] uint16_t getCharge() const;
    [[nodiscard]] const Peak1D getPeak() const;
    double calculateMass();


    /** @brief gets one of the peaks TO THE RIGHT and checks if the distance is in the AA-mass range
     * @param other: a peak to the right !
     */
    bool generateConnection( const std::shared_ptr<TagGeneratorNode>& other, double fragment_tolerance);

    /**
     * @brief Starts the recursiv generation of all Multi Peaks with the origin in this node
     * @param multi_peaks output
     * @param recursion_step the number of recursion steps to follow (equals the depth)
     */
    void generateAllMultiPeaks(std::vector<MultiPeak>& multi_peaks, uint8_t recursion_step);

    /**
     * @brief recursiv chain that generates the multipeaks
     * @param multi_peaks vector in which all results are safed
     * @param multi_peak The current multipeak which is constructed
     * @param recursion_step current step we are in
     * @param delta_mz the delta mz between the last node and this node
     * @param prev_AA the AA which resembles the delta_mz  // TODO: This is for debugging and can be removed in the final version
     */
    void generateAllMultiPeaksRecursion(std::vector<MultiPeak>& multi_peaks, MultiPeak multi_peak, uint32_t recursion_step, double delta_mz, std::string prev_AA);

    void calculateConfidence(const MSSpectrum& spectrum, double max_intensity);


  };
}
