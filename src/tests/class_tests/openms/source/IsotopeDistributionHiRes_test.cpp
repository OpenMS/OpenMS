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
// $Maintainer: Nikos Patikas $
// $Authors: Nikos Patikas $
// --------------------------------------------------------------------------
//

#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopePatternGenerator.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/MIDAsFFTID.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/MIDAsPolynomialID.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/Ecipex.h>



#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/test_config.h>
#include <OpenMS/SYSTEM/SysInfo.h>
#include <OpenMS/SYSTEM/StopWatch.h>

#include <iostream>
#include <fstream>
/////////////////////////////////////////////////////////////

START_TEST(IsotopeDistributionHires, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wshadow"


using namespace OpenMS;
using namespace std;


vector< pair<EmpiricalFormula, double> > formulas;
formulas.push_back(make_pair(EmpiricalFormula("C1000H1000N1000O1000"),0));

IsotopePatternGenerator* id;
double probability_cutoff = 0.00005;
double grid_resolution = 0.5;
SysInfo::MemUsage memory_watch;
StopWatch timer;


START_SECTION(Ecipex(double,double))
  
  for(auto& f : formulas)
  {
    LOG_INFO << f.first.toString() << endl;
    id = new Ecipex(0.00001, probability_cutoff);
    id->run(f.first);
    id->merge(grid_resolution);
    delete id;
  }
  

END_SECTION

START_SECTION(MIDAsFFTID(double,double))
  
  for(auto& f : formulas)
  {
    LOG_INFO << f.first.toString() << endl;
    id = new MIDAsFFTID(0.00001, probability_cutoff);
    id->run(f.first);
    id->merge(grid_resolution);
    delete id;
  }
  

END_SECTION

START_SECTION(MIDAsPolynomialID(double,double))
  
  for(auto& f : formulas)
  {
    LOG_INFO << f.first.toString() << endl;
    id = new MIDAsPolynomialID(0.00001, probability_cutoff);
    id->run(f.first);
    id->merge(grid_resolution);
    delete id;
  }
  
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

#pragma clang diagnostic pop
