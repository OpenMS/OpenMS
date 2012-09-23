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

#include <OpenMS/FILTERING/CALIBRATION/TOFCalibration.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <iostream>

using namespace OpenMS;
using namespace std;

Int main()
{
  TOFCalibration ec;
  PeakMap exp_raw, calib_exp;
  MzMLFile mzml_file;
  mzml_file.load("data/Tutorial_TOFCalibration_peak.mzML", calib_exp);
  mzml_file.load("data/Tutorial_TOFCalibration_raw.mzML", exp_raw);

  vector<DoubleReal> ref_masses;
  TextFile ref_file;
  ref_file.load("data/Tutorial_TOFCalibration_masses.txt", true);
  for (TextFile::Iterator iter = ref_file.begin(); iter != ref_file.end(); ++iter)
  {
    ref_masses.push_back(String(iter->c_str()).toDouble());
  }

  std::vector<DoubleReal> ml1;
  ml1.push_back(418327.924993827);

  std::vector<DoubleReal> ml2;
  ml2.push_back(253.645187196031);

  std::vector<DoubleReal> ml3;
  ml3.push_back(-0.0414243465397252);

  ec.setML1s(ml1);
  ec.setML2s(ml2);
  ec.setML3s(ml3);

  Param param;
  param.setValue("PeakPicker:peak_width", 0.1);
  ec.setParameters(param);
  ec.pickAndCalibrate(calib_exp, exp_raw, ref_masses);

  return 0;
} //end of main
