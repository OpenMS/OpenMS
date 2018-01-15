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
// $Maintainer: Timo Sachsenberg $
// $Authors: Lukas Mueller, Markus Mueller $
// --------------------------------------------------------------------------
//
///////////////////////////////////////////////////////////////////////////
//
//  PEAK DETECTION OF FOURIER TRANSFORME MS INSTRUMENT DATA
//
//  written by Markus Mueller, markus.mueller@imsb.biol.ethz.ch
//  ( and Lukas Mueller, Lukas.Mueller@imsb.biol.ethz.ch)
//  October 2005
//
//  Ported to OpenMS by Florian Zeller, florian.zeller@bsse.ethz.ch
//  December 2010
//
//  Group of Prof. Ruedi Aebersold, IMSB, ETH Hoenggerberg, Zurich
//
//

#include <string>
#include <vector>
#include <map>
#include <cmath>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/MSPeak.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/BackgroundIntensityBin.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/SuperHirnParameters.h>

namespace OpenMS
{

  using namespace std;

  BackgroundIntensityBin & BackgroundIntensityBin::operator=(const BackgroundIntensityBin & bib)
  {
    if (this == &bib)
      return *this;

    this->mzCoord_ = bib.mzCoord_;
    this->trCoord_ = bib.trCoord_;
    this->zCoord_ = bib.zCoord_;
    this->mean_ = bib.mean_;
    this->intensityMap_ = bib.intensityMap_;
    this->intensityHist_ = bib.intensityHist_;

    return *this;
  }

  BackgroundIntensityBin::BackgroundIntensityBin(const BackgroundIntensityBin & bib)
  {
    this->mzCoord_ = bib.mzCoord_;
    this->trCoord_ = bib.trCoord_;
    this->zCoord_ = bib.zCoord_;
    this->mean_ = bib.mean_;
    this->intensityMap_ = bib.intensityMap_;
    this->intensityHist_ = bib.intensityHist_;
  }

  BackgroundIntensityBin::BackgroundIntensityBin(double mz, double tr)
  {
    mzCoord_ = mz;
    trCoord_ = tr;
    zCoord_ = -1;
    mean_ = 0;
  }

// check if a peak belongs to this intenity bin
  bool BackgroundIntensityBin::checkBelonging(MSPeak * peak)
  {

    // check charge state:
    if (zCoord_ != -1)
    {
      if (peak->get_charge_state() != zCoord_)
      {
        return false;
      }
    }

    double TR_BINS = SuperHirnParameters::instance()->getBackgroundIntensityBinsTR();

    // check tr:
    double tr = peak->get_retention_time();
    if ((tr < (trCoord_ - TR_BINS / 2.0)) || (tr > (trCoord_ + TR_BINS / 2.0)))
    {
      return false;
    }

    double MZ_BINS = SuperHirnParameters::instance()->getBackgroundIntensityBinsMZ();

    double mz = peak->get_MZ();

    if ((mz < (mzCoord_ - MZ_BINS / 2.0)) || (mz > (mzCoord_ + MZ_BINS / 2.0)))
    {
      return false;
    }

    addIntensity(peak->get_intensity());
    peak = nullptr;
    return true;
  }

  void BackgroundIntensityBin::addMSPeak(MSPeak * peak)
  {
    addIntensity(peak->get_intensity());
    peak = nullptr;
  }

  void BackgroundIntensityBin::addIntensity(double intens)
  {
    intensityMap_.push_back(intens);
  }

// copied from simple_math
  double simple_math_WEIGHTED_AVERAGE(map<double, double> * in)
  {


    if (in->size() > 1)
    {
      double AVERAGE = 0;
      double TOT_WEIGHT = 0;

      map<double, double>::iterator start = in->begin();
      while (start != in->end())
      {
        TOT_WEIGHT += start->second;
        AVERAGE += (start->first * start->second);
        ++start;
      }

      return AVERAGE / TOT_WEIGHT;
    }
    else
    {
      return (*in->begin()).first;
    }
  }

// process collected intensities in the map
  void BackgroundIntensityBin::processIntensities()
  {

    computeIntensityHist();

    if (!intensityHist_.empty())
    {
      mean_ = simple_math_WEIGHTED_AVERAGE(&intensityHist_);
    }
    else
    {
      mean_ = 0;
    }
  }

// compute an intensity histogram
  void BackgroundIntensityBin::computeIntensityHist()
  {

    double constraint = SuperHirnParameters::instance()->getBackgroundIntensityBinsIntens();

    // insert into the histogram map
    vector<double>::iterator P = intensityMap_.begin();
    while (P != intensityMap_.end())
    {

      // intensity to bin:
      double intens = (*P);

      // find a key:
      map<double, double>::iterator F = intensityHist_.lower_bound(intens);
      if (F != intensityHist_.end())
      {

        // check this one:
        map<double, double>::iterator check = F;
        double mainLow = fabs(check->first - intens);
        if (check != intensityHist_.begin())
        {
          --check;
          double deltaHigh = fabs(check->first - intens);
          if (mainLow > deltaHigh)
          {
            mainLow = deltaHigh;
            F = check;
          }
        }
        if (mainLow > constraint)
        {
          F = intensityHist_.end();
        }
        else
        {
          F->second += 1.0;
        }
      }

      if (F == intensityHist_.end())
      {
        intensityHist_.insert(make_pair(intens, 1.0));
      }

      ++P;
    }

    // filter out bins of only 1 counts:
    map<double, double>::iterator F = intensityHist_.begin();
    while (F != intensityHist_.end())
    {

      if (F->second == SuperHirnParameters::instance()->getBackgroundIntensityBinsMinBinCount())
      {
        intensityHist_.erase(F++);
      }
      else
      {
        ++F;
      }
    }
  }

  BackgroundIntensityBin::~BackgroundIntensityBin()
  {
    intensityMap_.clear();
  }

  std::vector<double> * BackgroundIntensityBin::getIntensityMap()
  {   return &intensityMap_; }
  std::map<double, double> * BackgroundIntensityBin::getIntensityHist()
  {   return &intensityHist_; }
  double BackgroundIntensityBin::getMean()
  {   return mean_; }
}
