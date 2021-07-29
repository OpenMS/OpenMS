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
#include <OpenMS/test_config.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>

///////////////////////////
#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolutedSpectrum.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(DeconvolutedSpectrum, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

DeconvolutedSpectrum* ptr = 0;
DeconvolutedSpectrum* null_ptr = 0;
START_SECTION(DeconvolutedSpectrum())
{
  ptr = new DeconvolutedSpectrum();
  TEST_NOT_EQUAL(ptr, null_ptr){
  ptr = new DeconvolutedSpectrum();
  TEST_NOT_EQUAL(ptr, null_ptr)
}
}
END_SECTION

START_SECTION(~DeconvolutedSpectrum())
{
  delete ptr;
}
END_SECTION


// test data
PeakMap input;
MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("FLASHDeconv_sample_input.mzML"), input);
MSSpectrum test_spec = input[0];

/// detailed constructor
START_SECTION((DeconvolutedSpectrum(const MSSpectrum &spectrum, const int scan_number)))
{
  DeconvolutedSpectrum tmp_spec = DeconvolutedSpectrum(test_spec, 1);
  TEST_EQUAL(tmp_spec.getScanNumber(), 1);
  TEST_EQUAL(tmp_spec.getOriginalSpectrum().size(), 16124);
}
END_SECTION

DeconvolutedSpectrum test_deconv_spec = DeconvolutedSpectrum(test_spec, 1);

///
START_SECTION((int getScanNumber() const))
{
  int tmp_num = test_deconv_spec.getScanNumber();
  TEST_EQUAL(tmp_num, 1);
}
END_SECTION

START_SECTION((const MSSpectrum& getOriginalSpectrum() const))
{
  MSSpectrum tmp_s = test_deconv_spec.getOriginalSpectrum();
  TEST_EQUAL(tmp_s.size(), 16124);
}
END_SECTION

    START_SECTION((double getCurrentMinMass(const double min_mass) const))
        {
          // TODO
        }
    END_SECTION

///
START_SECTION((void writeDeconvolutedMasses(std::fstream &fs, const String &file_name, const FLASHDeconvHelperStructs::PrecalculatedAveragine &avg, const bool write_detail)))
    {
      // TODO
    }
END_SECTION

START_SECTION((void writeTopFD(std::fstream &fs, const FLASHDeconvHelperStructs::PrecalculatedAveragine &avg, const double snr_threshold=1.0, const double harmonic_factor=1.0, const double precursor_offset=.0)))
    {
      // TODO
    }
END_SECTION

START_SECTION((MSSpectrum toSpectrum(const int mass_charge)))
    {
      // TODO
    }
END_SECTION

START_SECTION((bool registerPrecursor(const std::vector< DeconvolutedSpectrum > &survey_scans, const std::map< int, std::vector< std::vector< double >>> &precursor_map_for_real_time_acquisition)))
    {
      // TODO
    }
END_SECTION



START_SECTION((PeakGroup getPrecursorPeakGroup() const))
    {
      // TODO
    }
END_SECTION

START_SECTION((int getPrecursorCharge() const))
    {
      // TODO
    }
END_SECTION

START_SECTION((const Precursor getPrecursor() const))
    {
      // TODO
    }
END_SECTION

START_SECTION((double getCurrentMaxMass(const double max_mass) const))
    {
      // TODO
    }
END_SECTION


START_SECTION((int getCurrentMaxAbsCharge(const int max_abs_charge) const))
    {
      // TODO
    }
END_SECTION



START_SECTION((int getPrecursorScanNumber() const))
    {
      // TODO
    }
END_SECTION

START_SECTION((static void writeDeconvolutedMassesHeader(std::fstream &fs, const int ms_level, const bool detail)))
    {
      // TODO
    }
END_SECTION


    /// copy constructor
    START_SECTION((DeconvolutedSpectrum(const DeconvolutedSpectrum &)))
        {
          // TODO
          DeconvolutedSpectrum tmp_spec(test_deconv_spec);

          TEST_EQUAL(tmp_spec.getScanNumber(), 1);
          TEST_EQUAL(tmp_spec.getOriginalSpectrum().size(), 16124);
        }
    END_SECTION

    ///// move constructor
    //START_SECTION((DeconvolutedSpectrum(DeconvolutedSpectrum &&other)))
    //    {
    //      // TODO
    //    }
    //END_SECTION

    /// assignment operator
    START_SECTION((DeconvolutedSpectrum& operator=(const DeconvolutedSpectrum &deconvoluted_spectrum)=default))
        {
          // TODO
        }
    END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
