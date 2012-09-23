// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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
// $Maintainer: Florian Zeller $
// $Authors: Markus Mueller $
// --------------------------------------------------------------------------
//
/*
 *  Deisotoper.h
 *  PeakDetection
 *
 *  Created by Markus Mueller on 10/19/06.
 *
 *  Ported to OpenMS by Florian Zeller, florian.zeller@bsse.ethz.ch
 *  December 2010
 *
 */

#ifndef _DEISOTOPER_H_
#define _DEISOTOPER_H_

#include <OpenMS/CONCEPT/Types.h>

#include <list>
#include <iostream>

namespace OpenMS
{

  class CentroidData;
  class DeconvPeak;

  class OPENMS_DLLAPI Deisotoper
  {
public:
//  static int sfMinCharge;
//  static int	sfMaxCharge;

    Deisotoper();
    Deisotoper(CentroidData &);
    virtual ~Deisotoper();

    std::list<DeconvPeak> & getDeconvPeaks() {return fDeconvPeaks; }

    void go(CentroidData &);
    void cleanDeconvPeaks();

    inline int getMinPeakGroupSize() {return fMinPeakGroupSize; }
    inline double getTheta() {return fTheta; }
    inline int getScanNumber() {return fScanNumber; }
    inline bool getShortReportFlag() {return fShortReportFlag; }

    inline void setMinPeakGroupSize(int pMinPeakGroupSize) {fMinPeakGroupSize = pMinPeakGroupSize; }
    inline void setTheta(double pTheta) {fTheta = pTheta; }
    inline void setScanNumber(int pScanNumber) {fScanNumber = pScanNumber; }
    inline void setShortReportFlag(bool pShortReportFlag) {fShortReportFlag = pShortReportFlag; }

protected:

    std::list<DeconvPeak> fDeconvPeaks;

    int     fMinPeakGroupSize;
    double  fTheta;
    int     fScanNumber;
    bool    fShortReportFlag;
  };

  std::ostream & operator<<(std::ostream &, Deisotoper &);

} // ns

#endif
