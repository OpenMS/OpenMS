//
// Created by trapho on 11/3/23.
//
#pragma once

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/Peak1D.h>


#include <utility>
#include <vector>
#include <functional>


namespace OpenMS
{

  class MultiPeak
  {
  public:
    MultiPeak();

    MultiPeak(Peak1D peak, double score);

    /// Copy
    MultiPeak(const MultiPeak& other);
    /// Assignment
    MultiPeak& operator=(const MultiPeak& other);
    /// Destructor
    virtual ~MultiPeak() = default;

    [[nodiscard]] const Peak1D& getPeak() const;
    double getScore() const;
    const std::string& getFollowUpPeaksAa() const;
    const std::vector<double>& getFollowUpPeaks() const;

    void addFollowUpPeak(double distance, const std::string& AA);
    void addScore(double score);

  protected:
    Peak1D peak_;
    double score_;
    std::string follow_up_peaks_AA;
    std::vector<double> follow_up_peaks;
  };
}