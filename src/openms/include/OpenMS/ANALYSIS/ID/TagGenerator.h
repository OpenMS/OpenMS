//
// Created by trapho on 10/20/23.
//
#pragma once

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/ANALYSIS/ID/TagGeneratorNode.h>

#include <vector>
#include <functional>


namespace OpenMS
{
  class OPENMS_DLLAPI TagGenerator
  {
  private:
    MSSpectrum spectrum_;
    std::vector<bool> selected_peaks_;
    uint32_t n; // the number of globally selected peaks;
    std::vector<TagGeneratorNode> dag_; // directed acyclic graph containing all peaks



  public:

    struct IdxAndIntensity{
      uint32_t idx_;
      double intensity_;
      IdxAndIntensity(uint32_t idx, double intensity){
        idx_ = idx;
        intensity_ = intensity;
      }
    };


    TagGenerator(const MSSpectrum& spectrum);
    ~TagGenerator();


    /**@brief top N peaks are selected according to their intensities over the entire m/z range of a spectrum where N is
     *related to a precursor ion mass
     *
     */
    void globalSelection();

    /**@brief peaks are selected by sliding a window of 70 Da (window increment, 35 Da) when fewer than two peaks are selected in any window during the global selection
     *
     */
    void localSelection();


    void generateDirectedAcyclicGraph(double fragment_tolerance);

    void generateAllQuadPeaks(std::vector<TagGeneratorNode::QuadPeak>& quad_peaks);

  };
}