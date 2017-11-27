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
// $Maintainer: Douglas McCloskey, Pasquale Domenico Colaianni $
// $Authors: Douglas McCloskey, Pasquale Domenico Colaianni $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/ANALYSIS/OPENSWATH/PeakIntegrator.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

void setup_toy_chromatogram(MSChromatogram& chromatogram)
{
  // Toy chromatogram
  // data is taken from raw LC-MS/MS data points acquired for L-Glutamate in RBCs
  std::vector<double> times = {
    2.23095,2.239716667,2.248866667,2.25765,2.266416667,
    2.275566667,2.2847,2.293833333,2.304066667,2.315033333,2.325983333,2.336566667,
    2.3468,2.357016667,2.367283333,2.377183333,2.387083333,2.39735,2.40725,2.4175,
    2.4274,2.4373,2.44755,2.45745,2.4677,2.477966667,2.488216667,2.498516667,2.5084,
    2.5183,2.5282,2.538466667,2.548366667,2.558266667,2.568516667,2.578783333,
    2.588683333,2.59895,2.6092,2.619466667,2.630066667,2.64065,2.65125,2.662116667,
    2.672716667,2.6833,2.6939,2.7045,2.715083333,2.725683333,2.736266667,2.746866667,
    2.757833333,2.768416667,2.779016667,2.789616667,2.8002,2.810116667,2.820033333,
    2.830316667,2.840216667,2.849766667,2.859316667,2.868866667,2.878783333,2.888683333,
    2.898233333,2.907783333,2.916033333,2.924266667,2.93215,2.940383333,2.947933333,
    2.955816667,2.964066667,2.97195,2.979833333,2.987716667,2.995616667,3.003516667,
    3.011416667,3.01895,3.026833333,3.034366667,3.042266667,3.0498,3.05735,3.065233333,
    3.073133333,3.080666667,3.0882,3.095733333,3.103633333,3.111533333,3.119066667,
    3.126966667,3.134866667,3.14275,3.15065,3.15855,3.166433333,3.174333333,3.182233333,
    3.190133333,3.198016667,3.205916667,3.213166667
  };

  std::vector<double> intensities = {
    1447,2139,1699,755,1258,1070,944,1258,1573,1636,
    1762,1447,1133,1321,1762,1133,1447,2391,692,1636,2957,1321,1573,1196,1258,881,
    1384,2076,1133,1699,1384,692,1636,1133,1573,1825,1510,2391,4342,10382,17618,
    51093,153970,368094,632114,869730,962547,966489,845055,558746,417676,270942,
    184865,101619,59776,44863,31587,24036,20450,20324,11074,9879,10508,7928,7110,
    6733,6481,5726,6921,6670,5537,4971,4719,4782,5097,5789,4279,5411,4530,3524,
    2139,3335,3083,4342,4279,3083,3649,4216,4216,3964,2957,2202,2391,2643,3524,
    2328,2202,3649,2706,3020,3335,2580,2328,2894,3146,2769,2517
  };

  chromatogram.clear(true);

  for (size_t i=0; i<times.size(); ++i)
  {
    ChromatogramPeak peak;
    peak.setMZ(times[i]);
    peak.setIntensity(intensities[i]);
    chromatogram.push_back(peak);
  }
}

START_TEST(PeakIntegrator, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

PeakIntegrator* ptr = 0;
PeakIntegrator* null_ptr = 0;

START_SECTION(PeakIntegrator())
{
  ptr = new PeakIntegrator();
  TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~PeakIntegrator())
{
  delete ptr;
}
END_SECTION

ptr = new PeakIntegrator();

START_SECTION(getParameters())
{
  Param params = ptr->getParameters();
  TEST_EQUAL(params.getValue("integration_type"), "trapezoid")
  TEST_EQUAL(params.getValue("baseline_type"), "vertical_division")
  TEST_EQUAL(params.getValue("peak_model"), "none")
}
END_SECTION

START_SECTION(setIntegrationType())
{
  TEST_EQUAL(ptr->getIntegrationType(), "trapezoid")
  ptr->setIntegrationType("intensity_sum");
  TEST_EQUAL(ptr->getIntegrationType(), "intensity_sum")
}
END_SECTION

START_SECTION(setBaselineType())
{
  TEST_EQUAL(ptr->getBaselineType(), "vertical_division")
  ptr->setBaselineType("base_to_base");
  TEST_EQUAL(ptr->getBaselineType(), "base_to_base")
}
END_SECTION

START_SECTION(getPeakModel())
{
  TEST_EQUAL(ptr->getPeakModel(), "none")
  // TODO update this to also test the setter
  //ptr->setPeakModel("base_to_base");
  //TEST_EQUAL(ptr->getPeakModel(), "base_to_base")
}
END_SECTION

START_SECTION(integratePeak())
{
  // inputs
  MSChromatogram chromatogram;
  setup_toy_chromatogram(chromatogram);
  double left = 2.477966667;
  double right = 3.01895;

  // outputs
  double peak_area;
  double peak_height;
  double peak_apex_pos; // (coordinates of mz or rt)

  ptr->setIntegrationType("simpson");
  ptr->integratePeak(chromatogram, left, right, peak_area, peak_height, peak_apex_pos);
  cout << "simpson: " << endl;
  TEST_REAL_SIMILAR(peak_area, 71720.443144994)
  TEST_REAL_SIMILAR(peak_height, 966489.0)
  TEST_REAL_SIMILAR(peak_apex_pos, 2.7045)

  ptr->setIntegrationType("intensity_sum");
  ptr->integratePeak(chromatogram, left, right, peak_area, peak_height, peak_apex_pos);
  cout << "intensity_sum: " << endl;
  TEST_REAL_SIMILAR(peak_area, 118750.49122807)
  TEST_REAL_SIMILAR(peak_height, 966489.0)
  TEST_REAL_SIMILAR(peak_apex_pos, 2.7045)

  ptr->setIntegrationType("trapezoid");
  ptr->integratePeak(chromatogram, left, right, peak_area, peak_height, peak_apex_pos);
  cout << "trapezoid: " << endl;
  TEST_REAL_SIMILAR(peak_area, 71540.2)
  TEST_REAL_SIMILAR(peak_height, 966489.0)
  TEST_REAL_SIMILAR(peak_apex_pos, 2.7045)

  // ofstream outfile;
  // outfile.open(
  //   OPENMS_GET_TEST_DATA_PATH("PeakIntegrator_integratePeak_output.html"),
  //   ios::out | ios::trunc
  // );
  // if (outfile.is_open())
  // {
  //   const string header = ""
  //   "<!doctype html>"
  //   "<html>"
  //   "<head>"
  //   "  <script src=\"https://cdn.plot.ly/plotly-latest.min.js\"></script>"
  //   "</head>"
  //   "<body>"
  //   "  <div id=\"plot-here\" style=\"height: 800px\"></div>"
  //   "left: " + to_string(left) + "<br>"
  //   "right: " + to_string(right) + "<br>"
  //   "area: " + to_string(peak_area) + "<br>"
  //   "</body>"
  //   "<script>"
  //   "const data = [";
  //
  //   outfile << header;
  //   outfile << "  {" << endl <<  "    x: [";
  //   for (auto p : chromatogram)
  //   {
  //     outfile << p.getRT() << ", ";
  //   }
  //   outfile << "]," << endl << "    y: [";
  //   for (auto p : chromatogram) {
  //     outfile << p.getIntensity() << ", ";
  //   }
  //   outfile << "]," << endl;
  //   string points_layout = ""
  //   "    mode: 'markers',"
  //   "    name: 'points',"
  //   "    type: 'scatter',"
  //   "    marker: {"
  //   "      color: 'red',"
  //   "      size: 4,"
  //   "      symbol: 'cross'"
  //   "    }"
  //   "  },";
  //
  //   outfile << points_layout << endl << "{" << endl << "    x: [";
  //   for (auto it=chromatogram.RTBegin(left); it!=chromatogram.RTEnd(right); ++it)
  //   {
  //     outfile << it->getRT() << ", ";
  //   }
  //   outfile << "]," << endl << "    y: [";
  //   for (auto it=chromatogram.RTBegin(left); it!=chromatogram.RTEnd(right); ++it) {
  //     outfile << it->getIntensity() << ", ";
  //   }
  //   outfile << "]," << endl;
  //   string lines_layout = ""
  //   "    name: 'area',"
  //   "    type: 'scatter',"
  //   "    mode: 'none',"
  //   "    fill: 'tozeroy',"
  //   "    fillcolor: 'lightskyblue'"
  //   "  },";
  //   outfile << lines_layout << endl;
  //
  //   const string footer = ""
  //   "];"
  //   "const layout = {"
  //   "  title: 'Chromatogram',"
  //   "  xaxis: {"
  //   "    title: 'Retention time'"
  //   "  },"
  //   "  yaxis: {"
  //   "    title: 'Intensity'"
  //   "  },"
  //   "  hovermode: 'closest'"
  //   "};"
  //   "Plotly.plot(document.getElementById('plot-here'), data, layout);"
  //   "</script>"
  //   "</html>";
  //
  //   outfile << footer;
  //   outfile.close();
  // }
  // else
  // {
  //   cout << "Unable to open file to write the spectrum";
  // }
}
END_SECTION

delete ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
