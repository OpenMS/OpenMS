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
// $Authors: Lukas Mueller, Markus Mueller $
// --------------------------------------------------------------------------
//
/*
 *  RawData.h
 *  PeakDetection
 *
 *  Created by Markus Mueller on 10/19/06.
 *
 *  Ported to OpenMS by Florian Zeller, florian.zeller@bsse.ethz.ch
 *  December 2010
 *
 */

#ifndef _RAWDATA_h_
#define _RAWDATA_h_
#include <vector>
#include <ostream>

#include <OpenMS/config.h>

namespace OpenMS
{

// Class for the storage of raw MS data
  class OPENMS_DLLAPI RawData
  {

public:

    RawData() {}
    RawData(std::vector<double> &, std::vector<double> &);
    virtual ~RawData();

    friend std::ostream & operator<<(std::ostream &, RawData &);

    /*
     * @brief Retrieve raw data as mass and intensity vectors. First argument: Mass values in profile mode
     * Second argument: Intensity values in profile mode
     */
    void get(std::vector<double> &, std::vector<double> &);

    /*
     * @brief Set raw data as mass and intensity vectors. First argument: Mass values in profile mode
     * Second argument: Intensity values in profile mode
     */
    void set(std::vector<double> &, std::vector<double> &);

    // Virtual functions
    virtual void smooth()
    {
    }

protected:
    std::vector<double> profileMasses_;
    std::vector<double> profileIntensities_;
  };

  std::ostream & operator<<(std::ostream & out, RawData & data);

  inline RawData::RawData(std::vector<double> & masses, std::vector<double> & intensities)
  {
    profileMasses_ = masses;
    profileIntensities_ = intensities;
  }

  inline void RawData::get(std::vector<double> & masses, std::vector<double> & intensities)
  {
    masses = profileMasses_;
    intensities = profileIntensities_;
  }

  inline void RawData::set(std::vector<double> & masses, std::vector<double> & intensities)
  {
    profileMasses_ = masses;
    profileIntensities_ = intensities;
  }

} // ns

#endif
