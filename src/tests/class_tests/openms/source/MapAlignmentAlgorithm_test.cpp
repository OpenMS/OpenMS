// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2014.
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
// $Maintainer: Clemens Groepl $
// $Authors: Marc Sturm, Clemens Groepl $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithm.h>

#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmPoseClustering.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmSpectrumAlignment.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmIdentification.h>

#include <OpenMS/CONCEPT/Factory.h>

using namespace std;
using namespace OpenMS;

/////////////////////////////////////////////////////////////

START_TEST(MapAlignmentAlgorithm, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////


MapAlignmentAlgorithm* ptr = 0;
MapAlignmentAlgorithm* nullPointer = 0;
START_SECTION((MapAlignmentAlgorithm()))
	ptr = new MapAlignmentAlgorithm();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((virtual ~MapAlignmentAlgorithm()))
	delete ptr;
END_SECTION

START_SECTION((virtual void alignPeakMaps(std::vector< MSExperiment<> > &, std::vector< TransformationDescription > &)))
  MapAlignmentAlgorithm ma;
  std::vector< MSExperiment<> > maps;
  std::vector<TransformationDescription> transformations;
  TEST_EXCEPTION(Exception::NotImplemented, ma.alignPeakMaps(maps,transformations));
END_SECTION

START_SECTION((virtual void alignFeatureMaps(std::vector< FeatureMap > &, std::vector< TransformationDescription > &)))
  MapAlignmentAlgorithm ma;
  std::vector< FeatureMap > maps;
  std::vector<TransformationDescription> transformations;
  TEST_EXCEPTION(Exception::NotImplemented, ma.alignFeatureMaps(maps,transformations));
END_SECTION

START_SECTION((virtual void alignPeptideIdentifications(std::vector< std::vector< PeptideIdentification > >&, std::vector<TransformationDescription>&)))
  MapAlignmentAlgorithm ma;
  std::vector< std::vector< PeptideIdentification > > maps;
  std::vector<TransformationDescription> transformations;
  TEST_EXCEPTION(Exception::NotImplemented, ma.alignPeptideIdentifications(maps,transformations));
END_SECTION

START_SECTION((static void registerChildren()))
{
  TEST_EQUAL(Factory<MapAlignmentAlgorithm>::registeredProducts().size(), 3)
  // I do not know why the classes show up in this particular order (sorted by name?).
	TEST_STRING_EQUAL(Factory<MapAlignmentAlgorithm>::registeredProducts()[0],MapAlignmentAlgorithmIdentification::getProductName());
	TEST_STRING_EQUAL(Factory<MapAlignmentAlgorithm>::registeredProducts()[1],MapAlignmentAlgorithmPoseClustering::getProductName());
	TEST_STRING_EQUAL(Factory<MapAlignmentAlgorithm>::registeredProducts()[2],MapAlignmentAlgorithmSpectrumAlignment::getProductName());
}
END_SECTION

START_SECTION((virtual void setReference(Size, const String&)))
{
	MapAlignmentAlgorithm ma;
	ma.setReference(); // no exception, nothing happens
	TEST_EXCEPTION(Exception::InvalidParameter, ma.setReference(1));
	TEST_EXCEPTION(Exception::InvalidParameter, ma.setReference(0, "test"));
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
