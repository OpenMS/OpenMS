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
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/SYSTEM/StopWatch.h>

/////////////////////////////////////////////////////////////

using namespace OpenMS;

START_TEST(StopWatch, "$Id$")

/////////////////////////////////////////////////////////////

START_SECTION((StopWatch& operator = (const StopWatch& stop_watch)))
  // ???
END_SECTION

START_SECTION((StopWatch()))
  // ???
END_SECTION

START_SECTION((StopWatch(const StopWatch& stop_watch)))
  // ???
END_SECTION

START_SECTION((bool isRunning() const))
  // ???
END_SECTION

START_SECTION((bool operator != (const StopWatch& stop_watch) const))
  // ???
END_SECTION

START_SECTION((bool operator < (const StopWatch& stop_watch) const))
  // ???
END_SECTION

START_SECTION((bool operator <= (const StopWatch& stop_watch) const))
  // ???
END_SECTION

START_SECTION((bool operator == (const StopWatch& stop_watch) const))
  // ???
END_SECTION

START_SECTION((bool operator > (const StopWatch& stop_watch) const))
  // ???
END_SECTION

START_SECTION((bool operator >= (const StopWatch& stop_watch) const))
  // ???
END_SECTION

START_SECTION((bool start()))
  // ???
END_SECTION

START_SECTION((bool stop()))
  // ???
END_SECTION

START_SECTION((DoubleReal getCPUTime() const ))
  // ???
END_SECTION

START_SECTION((DoubleReal getClockTime() const ))
  // ???
END_SECTION

START_SECTION((DoubleReal getSystemTime() const ))
  // ???
END_SECTION

START_SECTION((DoubleReal getUserTime() const ))
  // ???
END_SECTION

START_SECTION((void clear()))
  // ???
END_SECTION

START_SECTION((void reset()))
  // ???
END_SECTION

START_SECTION((virtual ~StopWatch()))
  // ???
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
