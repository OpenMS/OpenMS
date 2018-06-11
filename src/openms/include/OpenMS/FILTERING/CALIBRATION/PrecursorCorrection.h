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
// $Authors: Timo Sachsenberg, Oliver Alka $
// --------------------------------------------------------------------------
//

#pragma once

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/FeatureMap.h>

#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>

#include <OpenMS/METADATA/Precursor.h>

#include <fstream>
#include <iostream>
#include <string>
#include <algorithm>

namespace OpenMS
{

class OPENMS_DLLAPI PrecursorCorrection
{
    public:
    static const std::string csv_header;

    static void getPrecursors(
        const MSExperiment & exp, 
        std::vector<Precursor> & precursors, 
        std::vector<double> & precursors_rt, 
        std::vector<Size> precursor_scan_index);
    
    static void writeHist(const String& out_csv, const std::vector<double> & deltaMZs, const std::vector<double> & mzs, const std::vector<double> & rts);

    static std::set<Size> correctToNearestMS1Peak(
      MSExperiment & exp, 
      double mz_tolerance, 
      bool ppm, 
      std::vector<double> & deltaMZs, 
      std::vector<double> & mzs, 
      std::vector<double> & rts);

    //Selection of the peak with the highest intensity as corrected precursor mass in a given mass range (e.g. precursor mass +/- 0.2 Da)
    static std::set<Size> correctToHighestIntensityMS1Peak(
      MSExperiment & exp, 
      double mz_tolerance, 
      std::vector<double> & deltaMZs, 
      std::vector<double> & mzs, 
      std::vector<double> & rts);

    // Wrong assignment of the mono-isotopic mass for precursors are assumed:
    // - if precursor_mz matches the mz of a non-monoisotopic feature mass trace
    // - and in the case that believe_charge is true: if feature_charge matches the precursor_charge
    // In the case of wrong mono-isotopic assignment several options for correction are available:
    // keep_original will create a copy of the precursor and tandem spectrum for the new mono-isotopic mass trace and retain the original one
    // all_matching_features does this not for only the closest feature but all features in a question
    static std::set<Size> correctToNearestFeature(
      const FeatureMap& features, 
      MSExperiment & exp, 
      double rt_tolerance_s = 0.0, 
      double mz_tolerance = 0.0, 
      bool ppm = true, 
      bool believe_charge = false, 
      bool keep_original = false, 
      bool all_matching_features = false, 
      int max_trace = 2,
      int debug_level = 0);

  protected:
    static bool overlaps_(const Feature& feature, const double rt, const double pc_mz, const double rt_tolerance);
 
    static bool compatible_(const Feature& feature, double pc_mz, double mz_tolerance, Size max_trace_number = 2, int debug_level = 0);
};

const std::string PrecursorCorrection::csv_header = "RT,uncorrectedMZ,correctedMZ,deltaMZ";

}
