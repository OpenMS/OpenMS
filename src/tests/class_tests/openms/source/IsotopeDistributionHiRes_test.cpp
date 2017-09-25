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
//formulas.push_back(make_pair(EmpiricalFormula("C"),0));
//formulas.push_back(make_pair(EmpiricalFormula("(13)C100(2)H100(15)N100"),0));
//formulas.push_back(make_pair(EmpiricalFormula("CHNO"),0));

//formulas.push_back(make_pair(EmpiricalFormula("C10H10"),0));
//formulas.push_back(make_pair(EmpiricalFormula("C100H100"),0));
//formulas.push_back(make_pair(EmpiricalFormula("C1000H1000"),0));
//formulas.push_back(make_pair(EmpiricalFormula("C10000H10000"),0));

//formulas.push_back(make_pair(EmpiricalFormula("C10H10N10"),0));
//formulas.push_back(make_pair(EmpiricalFormula("C100H100N100"),0));
//formulas.push_back(make_pair(EmpiricalFormula("C1000H1000N1000"),0));

//formulas.push_back(make_pair(EmpiricalFormula("C10H10N10O10"),0));
//formulas.push_back(make_pair(EmpiricalFormula("C100H100N100O100"),0));
formulas.push_back(make_pair(EmpiricalFormula("C1000H1000N1000O1000"),0));
//formulas.push_back(make_pair(EmpiricalFormula("C10000H10000N10000O10000"),0));


IsotopePatternGenerator* id;
double probability_cutoff = 0.00005;
double grid_resolution = 0.5;
SysInfo::MemUsage memory_watch;
StopWatch timer;

// ofstream out("/home/main/fgid_results/results.csv");
// out << "Algorithm,Formula,TheoreticalIsotopes,IsotopesFound,MeanSkew,Runtime,Memory"<<endl;
// START_SECTION(Ecipex(double,double))
  
//   LOG_INFO <<"--------------" <<endl;
//   LOG_INFO <<"---Ecipex-----" <<endl;
  
//   for(auto& f : formulas)
//   {
//     LOG_INFO << f.first.toString() << endl;
//     f.second = f.first.calculateTheoreticalIsotopesNumber();
//     timer.reset();
//     memory_watch.reset();
//     timer.start();
    

//     id = new Ecipex(0.00001, probability_cutoff);
//     id->run(f.first);
//     id->merge(grid_resolution);
//     timer.stop();
//     UInt64 memory_kb;
//     SysInfo::getProcessMemoryConsumption(memory_kb);
//     double time_lapsed = timer.getClockTime();
    
//     out << "Ecipex," << f.first.toString() << "," << f.second << "," << id->size() << "," <<id->averageMass() - f.first.getAverageWeight()<<","<< time_lapsed <<"," << memory_kb << endl;
//     LOG_INFO << "Size " << id->size() << endl;
//     delete id;
//   }
  

// END_SECTION

// START_SECTION(MIDAsFFTID(double, double))

//   LOG_INFO <<"--------------" <<endl;
//   LOG_INFO <<"--MIDAsFFT----" <<endl;
//   for(auto& f : formulas)
//   {
//     LOG_INFO << f.first.toString();
//     f.second = f.first.calculateTheoreticalIsotopesNumber();
//     timer.reset();
//     memory_watch.reset();
//     timer.start();
    
//     id = new MIDAsFFTID(0.00001, probability_cutoff);
//     id->run(f.first);
//     LOG_INFO << id->averageMass() - f.first.getAverageWeight()<< endl;
//     id->merge(grid_resolution);
//     timer.stop();
//     UInt64 memory_kb;
//     SysInfo::getProcessMemoryConsumption(memory_kb);
//     double time_lapsed = timer.getClockTime();
//     out << "MIDAsFFT," << f.first.toString() << "," << f.second << "," << id->size() << "," <<id->averageMass() - f.first.getAverageWeight()<<","<< time_lapsed <<"," << memory_kb<< endl;
//     LOG_INFO << "Size " << id->size() << endl;
//     delete id;
//   }
// END_SECTION

// START_SECTION(MIDAsPolynomialID(double,double))
//       LOG_INFO <<"--------------" <<endl;
//   LOG_INFO <<"--MIDAsPolynomiaID----" <<endl;
//   for(auto& f : formulas)
//   {
//     LOG_INFO << f.first.toString();
//     f.second = f.first.calculateTheoreticalIsotopesNumber();
//     timer.reset();
//     memory_watch.reset();
//     timer.start();
//     id = new MIDAsPolynomialID(0.00001, probability_cutoff);
//     id->run(f.first);
//     LOG_INFO << id->averageMass() - f.first.getAverageWeight()<< endl;
//     id->merge(grid_resolution);
//     timer.stop();
//     UInt64 memory_kb;
//     SysInfo::getProcessMemoryConsumption(memory_kb);
//     double time_lapsed = timer.getClockTime();
    
//     out << "MIDAsPolynomial," << f.first.toString() << "," << f.second << "," << id->size() << "," <<id->averageMass() - f.first.getAverageWeight()<<","<< time_lapsed << ","<< memory_kb << endl;
//     LOG_INFO << "Size " << id->size() << endl;
//     delete id;
//   }
// END_SECTION

// out.close();


START_SECTION(Ecipex(double,double))
  
  LOG_INFO <<"--------------" <<endl;
  LOG_INFO <<"---Ecipex-----" <<endl;
  
  for(auto& f : formulas)
  {
    LOG_INFO << f.first.toString() << endl;
    f.second = f.first.calculateTheoreticalIsotopesNumber();
    timer.reset();
    memory_watch.reset();
    timer.start();



    id = new Ecipex(0.00001, probability_cutoff);
    id->run(f.first);
    LOG_INFO << id->averageMass() - f.first.getAverageWeight()<< endl;
    id->merge(grid_resolution);
    ofstream out("/home/main/fgid_results/ecipex.txt");
        for(auto& sample : id->getContainer())
          {
            out << sample.getMZ() << " "<<sample.getIntensity()<<endl;
          }

        out.close();

    timer.stop();
    UInt64 memory_kb;
    SysInfo::getProcessMemoryConsumption(memory_kb);
    double time_lapsed = timer.getClockTime();
    
    //out << "Ecipex," << f.first.toString() << "," << f.second << "," << id->size() << "," <<id->averageMass() - f.first.getAverageWeight()<<","<< time_lapsed <<"," << memory_kb << endl;
    LOG_INFO << "Size " << id->size() << endl;
    delete id;
  }
  

END_SECTION

START_SECTION(MIDAsFFTID(double, double))

  LOG_INFO <<"--------------" <<endl;
  LOG_INFO <<"--MIDAsFFT----" <<endl;
  for(auto& f : formulas)
  {
    LOG_INFO << f.first.toString();
    f.second = f.first.calculateTheoreticalIsotopesNumber();
    timer.reset();
    memory_watch.reset();
    timer.start();
    
    id = new MIDAsFFTID(0.00001, probability_cutoff);
    id->run(f.first);
    LOG_INFO << id->averageMass() - f.first.getAverageWeight()<< endl;
    id->merge(grid_resolution);
        ofstream out("/home/main/fgid_results/midas_fft.txt");
        for(auto& sample : id->getContainer())
          {
            out << sample.getMZ() << " "<<sample.getIntensity()<<endl;
          }

        out.close();
    timer.stop();
    UInt64 memory_kb;
    SysInfo::getProcessMemoryConsumption(memory_kb);
    double time_lapsed = timer.getClockTime();
    //out << "MIDAsFFT," << f.first.toString() << "," << f.second << "," << id->size() << "," <<id->averageMass() - f.first.getAverageWeight()<<","<< time_lapsed <<"," << memory_kb<< endl;
    LOG_INFO << "Size " << id->size() << endl;
    delete id;
  }
END_SECTION

START_SECTION(MIDAsPolynomialID(double,double))
      LOG_INFO <<"--------------" <<endl;
  LOG_INFO <<"--MIDAsPolynomiaID----" <<endl;
  for(auto& f : formulas)
  {
    LOG_INFO << f.first.toString();
    f.second = f.first.calculateTheoreticalIsotopesNumber();
    timer.reset();
    memory_watch.reset();
    timer.start();
    id = new MIDAsPolynomialID(0.00001, probability_cutoff);
    id->run(f.first);
    LOG_INFO << id->averageMass() - f.first.getAverageWeight()<< endl;
    id->merge(grid_resolution);
    ofstream out("/home/main/fgid_results/midas_pol.txt");
        for(auto& sample : id->getContainer())
          {
            out << sample.getMZ() << " "<<sample.getIntensity()<<endl;
          }

        out.close();
    timer.stop();
    UInt64 memory_kb;
    SysInfo::getProcessMemoryConsumption(memory_kb);
    double time_lapsed = timer.getClockTime();
    
    //out << "MIDAsPolynomial," << f.first.toString() << "," << f.second << "," << id->size() << "," <<id->averageMass() - f.first.getAverageWeight()<<","<< time_lapsed << ","<< memory_kb << endl;
    LOG_INFO << "Size " << id->size() << endl;
    delete id;
  }
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

#pragma clang diagnostic pop
