// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/FILTERING/DATAREDUCTION/PeakPattern.h>

#include <vector>
#include <algorithm>
#include <iostream>

using namespace std;

namespace OpenMS
{

	PeakPattern::PeakPattern(std::vector<double> ms, int msi, int c, int ppp)
    : massShifts_(ms), massShiftIndex_(msi), charge_(c), peaksPerPeptide_(ppp)
	{				
        // generate m/z shifts
        for (unsigned i = 0; i < massShifts_.size(); ++i)
        {
            for (int j = -1; j < peaksPerPeptide_; ++j)
            {
                // j=-1 shift corresponds to the zeroth peak
                mzShifts_.push_back((massShifts_[i] + j * Constants::C13C12_MASSDIFF_U)/charge_);
            }
        }
	}
    
    int PeakPattern::getCharge()
    {
        return charge_;
    }
    
    double PeakPattern::getMzShiftAt(int i)
    {
        return mzShifts_[i];
    }

    unsigned PeakPattern::getMzShiftCount()
    {
        return mzShifts_.size();
    }

    double PeakPattern::getMassShiftAt(int i)
    {
        return massShifts_[i];
    }

    unsigned PeakPattern::getMassShiftCount()
    {
        return massShifts_.size();
    }

    int PeakPattern::getMassShiftIndex()
    {
        return massShiftIndex_;
    }

    int PeakPattern::getPeaksPerPeptide()
    {
        return peaksPerPeptide_;
    }

    std::vector<double> PeakPattern::getMassShifts()
    {
        return massShifts_;
    }

}
