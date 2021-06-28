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
// $Maintainer: Axel Walter $
// $Authors: Axel Walter $
// --------------------------------------------------------------------------


#include <OpenMS/QC/DetectedCompounds.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/KERNEL/FeatureMap.h>

using namespace std;

namespace OpenMS
{ 

  DetectedCompounds::Result DetectedCompounds::compute(const String& pathToFeatureXMLFile)
  {
    DetectedCompounds::Result result;
    FeatureMap map;
    FeatureXMLFile f;
    f.load(pathToFeatureXMLFile, map);
    float sum_rt_deviations;
    UInt rt_count = 0;
    map.updateRanges();
    result.detected_compounds = map.size();
    for (const auto& f : map)
    {
      // if feature has peak_apex_position, get the meassured rt 
      if (f.getSubordinates()[0].metaValueExists("peak_apex_position"))
      {
        float rt_meassured = f.getSubordinates()[0].getMetaValue("peak_apex_position");
        // if feature has native id, get substring with theoretical rt, convert to float
        // and add absolute rt deviation for this feature to sum_rt_deviations, increment rt_count
        if (f.getSubordinates()[0].metaValueExists("native_id"))
        {
          String native_id = f.getSubordinates()[0].getMetaValue("native_id");
          UInt start = native_id.find("_rt") + 3;
          UInt end = native_id.find("_i");
          float rt_th = stof(native_id.substr(start, end-start));
          sum_rt_deviations += abs(rt_th - rt_meassured);
          rt_count += 1;
        }      
      }

    }
    // calculate mean rt shift (sec)
    result.rt_shift_mean = sum_rt_deviations/rt_count;
    return result;
  }

  bool DetectedCompounds::Result::operator==(const Result& rhs) const
  {
    return detected_compounds == rhs.detected_compounds
          && rt_shift_mean == rhs.rt_shift_mean;
  }

  /// Returns the name of the metric
  const String& DetectedCompounds::getName() const
  {
    return name_;
  }

  /// Returns required file input i.e. MzML.
  /// This is encoded as a bit in a Status object.
  QCBase::Status DetectedCompounds::requires() const
  {
    return QCBase::Status(QCBase::Requires::PREFDRFEAT);
  }
}
