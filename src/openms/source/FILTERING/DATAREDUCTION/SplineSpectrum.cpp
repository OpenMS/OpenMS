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
#include <OpenMS/FILTERING/DATAREDUCTION/SplinePackage.h>
#include <OpenMS/FILTERING/DATAREDUCTION/SplineSpectrum.h>
#include <OpenMS/MATH/MISC/Spline2d.h>

#include <vector>
#include <algorithm>
#include <iostream>

using namespace std;

namespace OpenMS
{

	SplineSpectrum::SplineSpectrum(std::vector<double> mz, std::vector<double> intensity) 
	{				
        SplineSpectrum::init(mz, intensity);
	}
    
    SplineSpectrum::SplineSpectrum(MSSpectrum<Peak1D> rawSpectrum) 
	{		
        std::vector<double> mz;
        std::vector<double> intensity;
        for (MSSpectrum<Peak1D>::Iterator it = rawSpectrum.begin(); it != rawSpectrum.end(); ++it)
        {
            mz.push_back(it->getMZ());
            intensity.push_back(it->getIntensity());
        }
		SplineSpectrum::init(mz, intensity);
    }
    
    void SplineSpectrum::init(std::vector<double> mz, std::vector<double> intensity) {
        
        if (!(mz.size() == intensity.size() && mz.size() > 2))
        {
            throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__,"m/z and intensity vectors either not of the same size or too short.");
        }
        
        const double newPackage = 2;    // start a new package if delta m/z is greater than newPackage times previous one
        
        //assert(mz.size() == intensity.size());
 		//assert(mz.size() > 2);
       
		// remove unnecessary zeros, i.e. zero intensity data points with zeros to the left and right
        std::vector<double> mzSlim;
        std::vector<double> intensitySlim;
        for (unsigned i=0; i<mz.size(); ++i)
        {
            int left = ((int)i-1 < 0)? 0 : ((int)i-1);
            int right = (i+1 > mz.size()-1)? (mz.size()-1) : (i+1);
            if (intensity[left] != 0 || intensity[i] != 0 || intensity[right] != 0) {
                mzSlim.push_back(mz[i]);
                intensitySlim.push_back(intensity[i]);
            }
         }
         
         // subdivide spectrum into packages
         std::vector<double> deltaMz;
         std::vector<bool> startPackage;
         deltaMz.push_back(0);
         for (unsigned i=1; i<mzSlim.size(); ++i) {
            deltaMz.push_back(mzSlim.at(i) - mzSlim[i-1]);
         }
         startPackage.push_back(true);
         for (unsigned i=1; i<mzSlim.size(); i++) {
            startPackage.push_back(deltaMz[i]/deltaMz[i-1] > newPackage);
         }
        
        // fill the packages
        std::vector<double> mzPackage;
        std::vector<double> intensityPackage;
        for (unsigned i=0; i<mz.size(); i++) {
            if (startPackage[i] && i > 0) {
                if (intensityPackage.size() > 2) {
                    // Three or more data points in package. At least one of them will be non-zero since unnecessary zeros removed above.
                    SplinePackage * package = new SplinePackage(mzPackage, intensityPackage);
                    packages_.push_back(*package);
                }
                mzPackage.clear();
                intensityPackage.clear();
            }
            mzPackage.push_back(mz[i]);
            intensityPackage.push_back(intensity[i]);
        }
        // add the last package
        if (intensityPackage.size() > 2) {
            SplinePackage * package = new SplinePackage(mzPackage, intensityPackage);
            packages_.push_back(*package);
        }
        mzPackage.clear();
        intensityPackage.clear();
   }

    SplineSpectrum::Navigator::Navigator(const std::vector<SplinePackage> * packages)
    {
        packages_ = packages;
        lastPackage_ = 0;
    }
    
    double SplineSpectrum::Navigator::eval(double mz) {
        
        SplinePackage start = (*packages_)[lastPackage_];
        if (mz < start.getMzMin()) {
            for (int i = lastPackage_; i >= 0; --i) {
                lastPackage_ = i;
                
                SplinePackage package = (*packages_)[i];
                if (mz > package.getMzMax()) {
                    return 0.0;
                }
                if (mz >= package.getMzMin()) {
                    return package.eval(mz);
                }
            }
        } else {
            for (int i = lastPackage_; i < (int)(*packages_).size(); ++i) {
                lastPackage_ = i;
                
                SplinePackage package = (*packages_)[i];
                if (mz < package.getMzMin()) {
                    return 0.0;
                }
                if (mz <= package.getMzMax()) {
                    return package.eval(mz);
                }
            }
        }
        return 0.0;
    }
    
    double SplineSpectrum::Navigator::getNextMz(double mz) {
        
        int minIndex = 0;
        int maxIndex = (*packages_).size() - 1;
        int i = lastPackage_;
        SplinePackage package = (*packages_)[i];
        
        // find correct package
        while (!(package.isInPackage(mz))) {
            if (mz < package.getMzMin()) {
                --i;
                package = (*packages_)[i];
                // check index limit
                if (i < minIndex) {
                    lastPackage_ = minIndex;
                    package = (*packages_)[minIndex];
                    return package.getMzMin();
                }
                // m/z in the gap?
                if (mz > package.getMzMax()) {
                    lastPackage_ = i;
                    package = (*packages_)[i + 1];
                    return package.getMzMin();
                }
            }
            else if (mz > package.getMzMax()) {
                i++;
                package = (*packages_)[i];
                // check index limit
                if (i > maxIndex) {
                    lastPackage_ = maxIndex;
                    package = (*packages_)[maxIndex];
                    return package.getMzMax();
               }
               // m/z in the gap?
               if (mz < package.getMzMin()) {
                   lastPackage_ = i;
                   package = (*packages_)[i];
                   return package.getMzMin();
               }
            }
        }
        
        // find m/z in the package
        if (mz + package.getMzStepWidth() > package.getMzMax()) {
            // The next step gets us outside the current package.
            // Let's move to the package to the right. 
            i++;
            package = (*packages_)[i];
            // check index limit
            if (i > maxIndex) {
                lastPackage_ = maxIndex;
                package = (*packages_)[maxIndex];
                return package.getMzMax();
            }
            // jump to min m/z of next package
            lastPackage_ = i;
            return package.getMzMin();
        } else {
            // make a small step within the package
            lastPackage_ = i;
            return mz + package.getMzStepWidth();
        }
    }
    
    SplineSpectrum::Navigator SplineSpectrum::getNavigator() {
        SplineSpectrum::Navigator * nav = new Navigator(&packages_);
        return *nav;
    }

}
