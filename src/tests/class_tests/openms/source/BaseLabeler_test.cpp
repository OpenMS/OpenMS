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
// $Maintainer: Timo Sachsenberg$
// $Authors: Stephan Aiche$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/SIMULATION/LABELING/BaseLabeler.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

class DerivedLabeler
  : public BaseLabeler
{
	  public:
	void preCheck(Param &) const override
  {
    throw Exception::NotImplemented(__FILE__,__LINE__,OPENMS_PRETTY_FUNCTION);
	}

	void setUpHook(SimTypes::FeatureMapSimVector & /* features */) override
	{
		throw Exception::NotImplemented(__FILE__,__LINE__,OPENMS_PRETTY_FUNCTION);
	}

	void postDigestHook(SimTypes::FeatureMapSimVector & /* features_to_simulate */) override
	{
		throw Exception::NotImplemented(__FILE__,__LINE__,OPENMS_PRETTY_FUNCTION);
	}

	void postRTHook(SimTypes::FeatureMapSimVector & /* features_to_simulate */) override
	{
		throw Exception::NotImplemented(__FILE__,__LINE__,OPENMS_PRETTY_FUNCTION);
	}

	void postDetectabilityHook(SimTypes::FeatureMapSimVector & /* features_to_simulate */) override
	{
		throw Exception::NotImplemented(__FILE__,__LINE__,OPENMS_PRETTY_FUNCTION);
	}

	void postIonizationHook(SimTypes::FeatureMapSimVector & /* features_to_simulate */) override
	{
		throw Exception::NotImplemented(__FILE__,__LINE__,OPENMS_PRETTY_FUNCTION);
	}

	void postRawMSHook(SimTypes::FeatureMapSimVector & /* features_to_simulate */) override
	{
		throw Exception::NotImplemented(__FILE__,__LINE__,OPENMS_PRETTY_FUNCTION);
	}

	void postRawTandemMSHook(SimTypes::FeatureMapSimVector & /* features_to_simulate */, SimTypes::MSSimExperiment & /* simulated map */) override
	{
		throw Exception::NotImplemented(__FILE__,__LINE__,OPENMS_PRETTY_FUNCTION);
	}

};

START_TEST(BaseLabeler, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

BaseLabeler* ptr = nullptr;
BaseLabeler* nullPointer = nullptr;
START_SECTION(BaseLabeler())
{
	ptr = new DerivedLabeler();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~BaseLabeler())
{
	delete ptr;
}
END_SECTION

DerivedLabeler labeler;
SimTypes::FeatureMapSimVector empty_fmsv;
SimTypes::MSSimExperiment empty_experiment;

START_SECTION((virtual void setUpHook(SimTypes::FeatureMapSimVector &)))
{
  TEST_EXCEPTION(Exception::NotImplemented, labeler.setUpHook(empty_fmsv))
}
END_SECTION

START_SECTION((virtual void postDigestHook(SimTypes::FeatureMapSimVector &)))
{
  TEST_EXCEPTION(Exception::NotImplemented, labeler.postDigestHook(empty_fmsv))
}
END_SECTION

START_SECTION((virtual void postRTHook(SimTypes::FeatureMapSimVector &)))
{
  TEST_EXCEPTION(Exception::NotImplemented, labeler.postRTHook(empty_fmsv))
}
END_SECTION

START_SECTION((virtual void postDetectabilityHook(SimTypes::FeatureMapSimVector &)))
{
  TEST_EXCEPTION(Exception::NotImplemented, labeler.postDetectabilityHook(empty_fmsv))
}
END_SECTION

START_SECTION((virtual void postIonizationHook(SimTypes::FeatureMapSimVector &)))
{
  TEST_EXCEPTION(Exception::NotImplemented, labeler.postIonizationHook(empty_fmsv))
}
END_SECTION

START_SECTION((virtual void postRawMSHook(SimTypes::FeatureMapSimVector &)))
{
  TEST_EXCEPTION(Exception::NotImplemented, labeler.postRawMSHook(empty_fmsv))
}
END_SECTION

START_SECTION((virtual void postRawTandemMSHook(SimTypes::FeatureMapSimVector &, SimTypes::MSSimExperiment &)))
{
  TEST_EXCEPTION(Exception::NotImplemented, labeler.postRawTandemMSHook(empty_fmsv, empty_experiment))
}
END_SECTION

START_SECTION((virtual Param getDefaultParameters() const ))
{
  Param p; // empty parameters
  TEST_EQUAL(labeler.getDefaultParameters(), p) // BaseLabeler should not have any parameters
}
END_SECTION

START_SECTION((virtual void setRnd(const SimRandomNumberGenerator &rng)))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((virtual void preCheck(Param &param) const =0))
{
  Param p;
  TEST_EXCEPTION(Exception::NotImplemented, labeler.preCheck(p))
}
END_SECTION

START_SECTION((ConsensusMap& getConsensus() ))
{
  ConsensusMap cm;
  TEST_EQUAL(labeler.getConsensus(), cm) // Consensus should be empty
}
END_SECTION

START_SECTION((String getChannelIntensityName(const Size channel_index) const ))
{
  TEST_STRING_EQUAL(labeler.getChannelIntensityName(1), "channel_1_intensity")
  TEST_STRING_EQUAL(labeler.getChannelIntensityName(100), "channel_100_intensity")
}
END_SECTION

START_SECTION((void registerChildren()))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((const String & getDescription() const))
{
  TEST_STRING_EQUAL(labeler.getDescription(), "")
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



