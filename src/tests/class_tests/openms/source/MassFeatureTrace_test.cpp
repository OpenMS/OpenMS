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
#include <OpenMS/ANALYSIS/TOPDOWN/MassFeatureTrace.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(MassFeatureTrace, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MassFeatureTrace* ptr = 0;
MassFeatureTrace* null_ptr = 0;
START_SECTION(MassFeatureTrace())
{
  ptr = new MassFeatureTrace();
  TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~MassFeatureTrace())
{
  delete ptr;
}
END_SECTION

/// sample input for testing



/// copy constructor test
START_SECTION((MassFeatureTrace(const MassFeatureTrace &)=default))
{
  // TODO : copy constructor
}
END_SECTION

START_SECTION((MassFeatureTrace(MassFeatureTrace &&other)=default))
{
  // TODO : assignment operator
}
END_SECTION

START_SECTION((MassFeatureTrace& operator=(const MassFeatureTrace &fd)=default))
{
  // TODO : assignment operator
}
END_SECTION

START_SECTION((void storeInformationFromDeconvolutedSpectrum(DeconvolutedSpectrum &deconvoluted_spectrum)))
    {
      // TODO
    }
END_SECTION

START_SECTION((void findFeatures(const String &file_name, const bool promex_out, const bool topfd_feature_out, const std::unordered_map< int, PeakGroup > &precursor_peak_groups, int &feature_cntr, int &feature_index, std::fstream &fsf, std::fstream &fsp, std::vector< std::fstream > &fst, const PrecalculatedAveragine &averagine)))
    {
      // TODO
    }
END_SECTION

START_SECTION((static void writeHeader(std::fstream &fs)))
    {
      // TODO
    }
END_SECTION

START_SECTION((static void writePromexHeader(std::fstream &fs)))
    {
      // TODO
    }
END_SECTION

START_SECTION((static void writeTopFDFeatureHeader(std::vector< std::fstream > &fs)))
    {
      // TODO
    }
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST


