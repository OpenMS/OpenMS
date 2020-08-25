// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
// $Authors: Markus Mueller $
// --------------------------------------------------------------------------
//
/*
 *  CentroidData.h
 *  PeakDetection
 *
 *  Created by Markus Mueller on 10/19/06.
 *
 *  Ported to OpenMS by Florian Zeller, florian.zeller@bsse.ethz.ch
 *  December 2010
 *
 */

#pragma once

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/SuperHirnConfig.h>

#include <OpenMS/CONCEPT/Types.h>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/CentroidPeak.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/RawData.h>

#include <list>
#include <vector>

#include <boost/shared_ptr.hpp>

namespace OpenMS
{

  class SUPERHIRN_DLLAPI CentroidData
  {
public:

    /*
     static double sfMassTolPpm;
     static double sfMassTolDa;
     static double sfMinIntensity;
     static double sfIntensityFloor;
     */
    // debugging parameters => used by other classes :( :( :(
    /*
     static bool MonoIsoDebugging;
     static double DebugMonoIsoMassMin;
     static double DebugMonoIsoMassMax;
     */

    CentroidData(int, boost::shared_ptr<RawData>, bool);
    CentroidData(int, boost::shared_ptr<RawData>, double, bool);
    virtual ~CentroidData();

    void get(std::list<CentroidPeak> &);
    void set(boost::shared_ptr<RawData>);
    void set(std::vector<double> &, std::vector<double> &);

    void setWidth(int pWidth);

    int getWidth();

    void setNoise(double);
    double getNoise();

    void removeNoise();

    bool getNextPeakGroup(std::list<CentroidPeak>::iterator &, std::list<CentroidPeak>::iterator &);
    void resetPeakGroupIter();

    bool centroidDataModus_;

protected:

    void calcCentroids(boost::shared_ptr<RawData>);

    int fWindowWidth;
    double fNoise;
    double fScanRetentionTime;
    std::list<CentroidPeak> fCentroidPeaks;
    std::list<CentroidPeak>::iterator fPeakGroupStart;
  };

  std::ostream & operator<<(std::ostream &, CentroidData &);

} // ns

