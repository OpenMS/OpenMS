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
// $Maintainer: Timo Sachsenberg$
// $Authors: Stephan Aiche, Chris Bielow$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/SIMULATION/RawMSSignalSimulation.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(RawMSSignalSimulation, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

RawMSSignalSimulation* ptr = nullptr;
RawMSSignalSimulation* nullPointer = nullptr;
SimTypes::MutableSimRandomNumberGeneratorPtr empty_rnd_gen (new SimTypes::SimRandomNumberGenerator);
//const unsigned long rnd_gen_seed = 1;

START_SECTION((RawMSSignalSimulation(SimRandomNumberGeneratorPtr rng)))
{
  ptr = new RawMSSignalSimulation(empty_rnd_gen);
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~RawMSSignalSimulation())
{
	delete ptr;
}
END_SECTION

START_SECTION((RawMSSignalSimulation(const RawMSSignalSimulation &source)))
{
  RawMSSignalSimulation source(empty_rnd_gen);
  Param p = source.getParameters();
  p.setValue("peak_fwhm",0.3);
  source.setParameters(p);

  RawMSSignalSimulation target(source);
  TEST_EQUAL(source.getParameters(), target.getParameters())
}
END_SECTION

START_SECTION((RawMSSignalSimulation& operator=(const RawMSSignalSimulation &source)))
{
  RawMSSignalSimulation source(empty_rnd_gen);
  RawMSSignalSimulation target(source);

  Param p = source.getParameters();
  p.setValue("peak_fwhm",0.3);
  source.setParameters(p);
  TEST_NOT_EQUAL(source.getParameters(), target.getParameters())

  target = source;

  TEST_EQUAL(source.getParameters(), target.getParameters())
}
END_SECTION

START_SECTION((void generateRawSignals(SimTypes::FeatureMapSim &features, SimTypes::MSSimExperiment &experiment, SimTypes::MSSimExperiment &experiment_ct, SimTypes::FeatureMapSim &contaminants)))
{
  // TODO
}
END_SECTION


START_SECTION((void loadContaminants()))
{
  // TODO
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



