// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2021.
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
// $Maintainer: Jihyung Kim$
// $Authors: Jihyung Kim$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHIda.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(FLASHIda, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

FLASHIda* ptr = 0;
START_SECTION(~FLASHIda())
{
  delete ptr;
}
END_SECTION

/// < public methods without tests >
/// - default constructors and operators are not used (copy, move, assignment)
/// - default constructor taking string input argument


FLASHIda *f_ida = new FLASHIda("");
    START_SECTION((int getPeakGroups(const double *mzs, const double *intensities, const int length, const double rt, const int ms_level, const char *name)))
        {
          double mz_array[4] = {1.1, 2.5, 3.6, 4.7};
          double inty_array[4] = {1000, 1500, 2000, 1500};
          const char *spec_name = "test_spec";

          int pg_size = f_ida->getPeakGroups(mz_array, inty_array, 4, 5., 2, spec_name);

          TEST_EQUAL(pg_size, 0);
          // TODO : add MS1 example
        }
    END_SECTION

//    START_SECTION((void getIsolationWindows(double *window_start, double *window_end, double *qscores, int *charges, int *min_charges, int *max_charges, double *mono_masses, double *charge_cos, double *charge_snrs, double *iso_cos, double *snrs, double *charge_scores, double *ppm_errors, double *precursor_intensities, double *peakgroup_intensities)))
//        {
//          // TODO
//        }
//    END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
