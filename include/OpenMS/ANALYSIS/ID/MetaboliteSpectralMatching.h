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
// $Maintainer: Erhan Kenar $
// $Authors: Erhan Kenar $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_ID_METABOLITESPECTRALMATCHING
#define OPENMS_ANALYSIS_ID_METABOLITESPECTRALMATCHING

#include <OpenMS/KERNEL/MassTrace.h>
#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/FORMAT/MzTab.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>


#include <vector>
#include <algorithm>

namespace OpenMS
{
struct OPENMS_DLLAPI PrecursorMassComparator
{
    bool operator() (const MSSpectrum<Peak1D>& a, const MSSpectrum<Peak1D>& b)
    {
        return a.getPrecursors()[0].getMZ() < b.getPrecursors()[0].getMZ();
    }

} PrecursorMZLess;

class OPENMS_DLLAPI MetaboliteSpectralMatching :
public DefaultParamHandler,
public ProgressLogger
{
public:
    /// Default constructor
    MetaboliteSpectralMatching();

    // Explicit constructor
    MetaboliteSpectralMatching(const String& map_fname);

    /// Default destructor
    virtual ~MetaboliteSpectralMatching();


    /// main method of MetaboliteSpectralMatching
    void run(MSExperiment<>&, MzTab&);


protected:
    virtual void updateMembers_();

private:
    /// private member functions

    DoubleReal precursor_mz_error;
    DoubleReal fragment_mz_error;
    String precursor_mz_error_unit;
};


}




#endif // OPENMS_ANALYSIS_ID_METABOLITESPECTRALMATCHING
