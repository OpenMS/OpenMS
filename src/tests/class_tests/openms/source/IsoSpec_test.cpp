// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Marc Sturm, Michał Startek, Mateusz Łącki $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsoSpecWrapper.h>
///////////////////////////

#include <OpenMS/CHEMISTRY/Element.h>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>

using namespace OpenMS;
using namespace std;

// ----------------------------------------------------------------------------------------------------------------------
// Setup for tests
// ----------------------------------------------------------------------------------------------------------------------


typedef Peak1D isopair;
std::vector<isopair> fructose_expected_oms; // A few initial isotopologues for the fructose molecule

#define ISOSPEC_TEST_EPSILON 0.0000001
// Test with more precision than TEST::isRealSimilar, without side effects, and w/o being chatty about it.
bool my_real_similar(double a, double b)
{
  return a * (1.0-ISOSPEC_TEST_EPSILON) <= b && b <= a * (1.0+ISOSPEC_TEST_EPSILON);
}

#define ISOSPEC_TEST_ASSERTION(b) \
if(!(b)) \
{ \
  std::cout << "Failing assertion in line: " << __LINE__ << std::endl; \
  return false; \
}

bool compare_to_reference(IsotopeDistribution& ID, const std::vector<Peak1D>& reference)
{
  std::sort(ID.begin(), ID.end(),  [](isopair a, isopair b) {return a.getIntensity() > b.getIntensity();});

  for (Size i = 0; i != reference.size(); ++i)
  {
    ISOSPEC_TEST_ASSERTION(my_real_similar(ID[i].getPos(), reference[i].getPos()));
    ISOSPEC_TEST_ASSERTION(my_real_similar(ID[i].getIntensity(), reference[i].getIntensity()));
  }

  return true;
}


Size generator_length(IsoSpecGeneratorWrapper& IW)
{
  Size i = 0;
  while(IW.nextConf())
    i++;
  return i;
}

// With empty vector as reference this function will just run some sanity checks on the generator output
// confs_to_extract == -1 will test the generator until exhaustion, >0 will just test the initial n confs.
bool compare_generator_to_reference(IsoSpecGeneratorWrapper& IW, const std::vector<Peak1D>& reference, UInt32 confs_to_extract)
{
  Size matches_count = 0;
  std::vector<Peak1D> isoResult;
  while(IW.nextConf() && confs_to_extract != 0)
  {
    Peak1D p = IW.getConf();
    ISOSPEC_TEST_ASSERTION(p.getPos() == IW.getMass());
    ISOSPEC_TEST_ASSERTION(p.getIntensity() == static_cast<float>(IW.getIntensity()));
    ISOSPEC_TEST_ASSERTION(my_real_similar(IW.getIntensity(), exp(IW.getLogIntensity())))

    for(auto it = reference.begin(); it != reference.end(); it++)
      if(my_real_similar(it->getPos(), IW.getMass()) && my_real_similar(it->getIntensity(), IW.getIntensity()))
        matches_count++;

    confs_to_extract--;
  }
  ISOSPEC_TEST_ASSERTION(matches_count == reference.size());
  return true;
}

START_TEST(IsoSpecWrapper, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

// A few initial isotopologues for the fructose molecule
fructose_expected_oms.push_back(Peak1D( 180.06339038280000863778696 , 0.92263317941561073798339976    ));
fructose_expected_oms.push_back(Peak1D( 181.06674538280000774648215 , 0.059873700450778437331944559   ));
fructose_expected_oms.push_back(Peak1D( 181.06760738279999145561305 , 0.0021087279716237856790062022  ));
fructose_expected_oms.push_back(Peak1D( 181.06966713090000098418386 , 0.001273380225742650438680581   ));
fructose_expected_oms.push_back(Peak1D( 182.06764438280001172643097 , 0.011376032168236337518973933   ));
fructose_expected_oms.push_back(Peak1D( 182.07010038280000685517734 , 0.0016189442373783591386932068  ));
fructose_expected_oms.push_back(Peak1D( 182.07096238279999056430825 , 0.00013684457672024180033450158 ));
fructose_expected_oms.push_back(Peak1D( 182.07302213090000009287905 , 8.2635209633747614224076605e-05 ));
fructose_expected_oms.push_back(Peak1D( 183.07099938280001083512616 , 0.00073824045954083467209472236 ));
fructose_expected_oms.push_back(Peak1D( 183.07186138279999454425706 , 2.1667113372227351532611078e-05 ));
fructose_expected_oms.push_back(Peak1D( 183.07345538280000596387254 , 2.3346748674918047107952959e-05 ));
fructose_expected_oms.push_back(Peak1D( 183.07392113090000407282787 , 1.5700729969000005987024918e-05 ));
fructose_expected_oms.push_back(Peak1D( 184.07189838280001481507497 , 5.8444185791655326584186775e-05 ));
fructose_expected_oms.push_back(Peak1D( 184.07435438280000994382135 , 1.9961521148266482778097647e-05 ));

std::sort(fructose_expected_oms.begin(), fructose_expected_oms.end(),  [](isopair a, isopair b) {return a.getIntensity() > b.getIntensity();});

EmpiricalFormula ef_fructose("C6H12O6");

std::vector<int> fructose_isotopeNumbers;
std::vector<int> fructose_atomCounts;
std::vector<std::vector<double> > fructose_isotopeMasses;
std::vector<std::vector<double> > fructose_isotopeProbabilities;

for (const auto& elem : ef_fructose)
{
  fructose_atomCounts.push_back( elem.second );

  std::vector<double> masses;
  std::vector<double> probs;
  for (const auto& iso : elem.first->getIsotopeDistribution())
  {
    if (iso.getIntensity() <= 0.0) continue; // Note: there will be an Isospec exception if one of the intensities is zero!
    masses.push_back(iso.getMZ());
    probs.push_back(iso.getIntensity());
  }
  fructose_isotopeNumbers.push_back( masses.size() );
  fructose_isotopeMasses.push_back(masses);
  fructose_isotopeProbabilities.push_back(probs);
}


// Create an invalid molecule: where one of the isotopic intensities is defined to be zero
std::vector<int> invalid_isotopeNumbers = fructose_isotopeNumbers;
std::vector<int> invalid_atomCounts = fructose_atomCounts;
std::vector<std::vector<double> > invalid_isotopeMasses = fructose_isotopeMasses;
std::vector<std::vector<double> > invalid_isotopeProbabilities = fructose_isotopeProbabilities;

invalid_isotopeNumbers[0] += 1;
invalid_isotopeMasses[0].push_back(3.0160492699999998933435563);
invalid_isotopeProbabilities[0].push_back(0.0);


// ----------------------------------------------------------------------------------------------------------------------
// Tests: IsoSpecThresholdGeneratorWrapper
// ----------------------------------------------------------------------------------------------------------------------


{
  IsoSpecGeneratorWrapper* ptr = nullptr;
  IsoSpecGeneratorWrapper* ptr2 = nullptr;
  IsoSpecGeneratorWrapper* nullPointer = nullptr;
  START_SECTION((IsoSpecThresholdGeneratorWrapper::IsoSpecThresholdGeneratorWrapper(const EmpiricalFormula&, double, bool)))
    ptr = new IsoSpecThresholdGeneratorWrapper(EmpiricalFormula("C10"), 0.5, false);
    TEST_NOT_EQUAL(ptr, nullPointer)
  END_SECTION

  START_SECTION((IsoSpecThresholdGeneratorWrapper(std::vector<int>, std::vector<int>,
                                         std::vector<std::vector<double> >,
                                         std::vector<std::vector<double> >, 
                                         double, bool)))
    ptr2 = new IsoSpecThresholdGeneratorWrapper(fructose_isotopeNumbers, fructose_atomCounts, fructose_isotopeMasses,
                                       fructose_isotopeProbabilities, 0.5, false);
    TEST_NOT_EQUAL(ptr2, nullPointer)
    TEST_EXCEPTION(Exception::IllegalArgument, IsoSpecThresholdGeneratorWrapper(invalid_isotopeNumbers, invalid_atomCounts, invalid_isotopeMasses, invalid_isotopeProbabilities, 0.5, false));
  END_SECTION

  START_SECTION((IsoSpecThresholdGeneratorWrapper::~IsoSpecThresholdGeneratorWrapper()))
    delete ptr;
    delete ptr2;
  END_SECTION
}

START_SECTION(( bool IsoSpecThresholdGeneratorWrapper::nextConf() ))
{
  double threshold = 1e-5;
  bool absolute = false;


  IsoSpecThresholdGeneratorWrapper ITW(EmpiricalFormula("C6H12O6"), threshold, absolute);
  TEST_EQUAL(compare_generator_to_reference(ITW, fructose_expected_oms, -1), true);

  IsoSpecThresholdGeneratorWrapper ITW2(EmpiricalFormula("C6H12O6"), threshold, absolute);
  TEST_EQUAL(generator_length(ITW2), 14);

  IsoSpecThresholdGeneratorWrapper ITW3(fructose_isotopeNumbers, fructose_atomCounts, fructose_isotopeMasses,
                   fructose_isotopeProbabilities, threshold, absolute);
  TEST_EQUAL(compare_generator_to_reference(ITW3, fructose_expected_oms, -1), true);


  // human insulin
  IsoSpecThresholdGeneratorWrapper ITW4(EmpiricalFormula("C520H817N139O147S8"), threshold, absolute);
  TEST_EQUAL(generator_length(ITW4), 5513);

  IsoSpecThresholdGeneratorWrapper ITW5(EmpiricalFormula("C520H817N139O147S8"), 0.01, absolute);
  TEST_EQUAL(generator_length(ITW5), 267)
}
END_SECTION

START_SECTION(( Peak1D IsoSpecThresholdGeneratorWrapper::getConf() ))
  NOT_TESTABLE; // Tested with nextConf(), above
END_SECTION

START_SECTION(( double IsoSpecThresholdGeneratorWrapper::getMass() ))
  NOT_TESTABLE; // Tested with nextConf(), above
END_SECTION

START_SECTION(( double IsoSpecThresholdGeneratorWrapper::getIntensity() ))
  NOT_TESTABLE; // Tested with nextConf(), above
END_SECTION

START_SECTION(( double IsoSpecThresholdGeneratorWrapper::getLogIntensity() ))
  NOT_TESTABLE; // Tested with nextConf(), above
END_SECTION


// ----------------------------------------------------------------------------------------------------------------------
// Tests: IsoSpecTotalProbGeneratorWrapper
// ----------------------------------------------------------------------------------------------------------------------


{
  IsoSpecGeneratorWrapper* ptr = nullptr;
  IsoSpecGeneratorWrapper* ptr2 = nullptr;
  IsoSpecGeneratorWrapper* nullPointer = nullptr;
  START_SECTION((IsoSpecTotalProbGeneratorWrapper::IsoSpecTotalProbGeneratorWrapper(const EmpiricalFormula&, double)))
    ptr = new IsoSpecTotalProbGeneratorWrapper(EmpiricalFormula("C10"), 0.5);
    TEST_NOT_EQUAL(ptr, nullPointer);
  END_SECTION

  START_SECTION((IsoSpecTotalProbGeneratorWrapper::IsoSpecTotalProbGeneratorWrapper(std::vector<int>, std::vector<int>,
                                         std::vector<std::vector<double> >,
                                         std::vector<std::vector<double> >,
                                         double, bool)))
    ptr2 = new IsoSpecTotalProbGeneratorWrapper(fructose_isotopeNumbers, fructose_atomCounts, fructose_isotopeMasses,
                                    fructose_isotopeProbabilities, 0.5);
    TEST_NOT_EQUAL(ptr2, nullPointer);

    TEST_EXCEPTION(Exception::IllegalArgument, IsoSpecTotalProbGeneratorWrapper(invalid_isotopeNumbers, invalid_atomCounts,
                                                          invalid_isotopeMasses, invalid_isotopeProbabilities, 0.5));

  END_SECTION


  START_SECTION((IsoSpecTotalProbGeneratorWrapper::~IsoSpecTotalProbGeneratorWrapper()))
    delete ptr;
    delete ptr2;
  END_SECTION
}


START_SECTION(( bool IsoSpecTotalProbGeneratorWrapper::nextConf() ))
{
  double total_prob = 0.99999;

  IsoSpecTotalProbGeneratorWrapper ITPW(EmpiricalFormula("C6H12O6"), total_prob);
  TEST_EQUAL(compare_generator_to_reference(ITPW, fructose_expected_oms, -1), true);

  IsoSpecTotalProbGeneratorWrapper ITPW2(EmpiricalFormula("C6H12O6"), total_prob);
  TEST_EQUAL(generator_length(ITPW2), 2548); // Total prob is a hint, if we keep extracting configurations
                                             // we will go over the entire configuration space

  IsoSpecTotalProbGeneratorWrapper ITPW3(fructose_isotopeNumbers, fructose_atomCounts, fructose_isotopeMasses, 
                               fructose_isotopeProbabilities, total_prob);
  TEST_EQUAL(compare_generator_to_reference(ITPW3, fructose_expected_oms, -1), true);
}
END_SECTION


START_SECTION(( Peak1D IsoSpecTotalProbGeneratorWrapper::getConf() ))
  NOT_TESTABLE; // Tested with nextConf(), above
END_SECTION

START_SECTION(( double IsoSpecTotalProbGeneratorWrapper::getMass() ))
  NOT_TESTABLE; // Tested with nextConf(), above
END_SECTION

START_SECTION(( double IsoSpecTotalProbGeneratorWrapper::getIntensity() ))
  NOT_TESTABLE; // Tested with nextConf(), above
END_SECTION

START_SECTION(( double IsoSpecTotalProbGeneratorWrapper::getLogIntensity() ))
  NOT_TESTABLE; // Tested with nextConf(), above
END_SECTION

// ----------------------------------------------------------------------------------------------------------------------


{
  IsoSpecGeneratorWrapper* ptr = nullptr;
  IsoSpecGeneratorWrapper* ptr2 = nullptr;
  IsoSpecGeneratorWrapper* nullPointer = nullptr;
  START_SECTION((IsoSpecOrderedGeneratorWrapper::IsoSpecOrderedGeneratorWrapper(const EmpiricalFormula&)))
    ptr = new IsoSpecOrderedGeneratorWrapper(EmpiricalFormula("C10"));
    TEST_NOT_EQUAL(ptr, nullPointer)
  END_SECTION

  START_SECTION((IsoSpecOrderedGeneratorWrapper::IsoSpecOrderedGeneratorWrapper(std::vector<int>, std::vector<int>,
                                         std::vector<std::vector<double> >,
                                         std::vector<std::vector<double> >)))
    ptr2 = new IsoSpecOrderedGeneratorWrapper(fructose_isotopeNumbers, fructose_atomCounts, fructose_isotopeMasses,
                                    fructose_isotopeProbabilities);
    TEST_NOT_EQUAL(ptr2, nullPointer);

    TEST_EXCEPTION(Exception::IllegalArgument, IsoSpecOrderedGeneratorWrapper(invalid_isotopeNumbers, invalid_atomCounts,
                                                          invalid_isotopeMasses, invalid_isotopeProbabilities));

  END_SECTION

  START_SECTION((~IsoSpecOrderedGeneratorWrapper()))
    delete ptr;
    delete ptr2;
  END_SECTION
}


START_SECTION(( bool IsoSpecOrderedGeneratorWrapper::nextConf() ))
{
  IsoSpecOrderedGeneratorWrapper IOGW(EmpiricalFormula("C6H12O6"));
  TEST_EQUAL(compare_generator_to_reference(IOGW, fructose_expected_oms, fructose_expected_oms.size()), true);

  IsoSpecOrderedGeneratorWrapper IOGW2(EmpiricalFormula("C6H12O6"));
  TEST_EQUAL(generator_length(IOGW2), 2548);

  IsoSpecOrderedGeneratorWrapper IOGW3(fructose_isotopeNumbers, fructose_atomCounts, fructose_isotopeMasses, 
                                       fructose_isotopeProbabilities);
  TEST_EQUAL(compare_generator_to_reference(IOGW3, fructose_expected_oms, -1), true);

  // human insulin
  IsoSpecOrderedGeneratorWrapper IOGW4(EmpiricalFormula("C520H817N139O147S8"));
  TEST_EQUAL(compare_generator_to_reference(IOGW4, std::vector<Peak1D>(), 10000), true);
}
END_SECTION


START_SECTION(( Peak1D IsoSpecOrderedGeneratorWrapper::getConf() ))
  NOT_TESTABLE; // Tested with nextConf(), above
END_SECTION

START_SECTION(( double IsoSpecOrderedGeneratorWrapper::getMass() ))
  NOT_TESTABLE; // Tested with nextConf(), above
END_SECTION

START_SECTION(( double IsoSpecOrderedGeneratorWrapper::getIntensity() ))
  NOT_TESTABLE; // Tested with nextConf(), above
END_SECTION

START_SECTION(( double IsoSpecOrderedGeneratorWrapper::getLogIntensity() ))
  NOT_TESTABLE; // Tested with nextConf(), above
END_SECTION



// ----------------------------------------------------------------------------------------------------------------------
// Tests: IsoSpecThresholdWrapper
// ----------------------------------------------------------------------------------------------------------------------


{
  IsoSpecWrapper* ptr = nullptr;
  IsoSpecWrapper* ptr2 = nullptr;
  IsoSpecWrapper* nullPointer = nullptr;
  START_SECTION((IsoSpecThresholdWrapper::IsoSpecThresholdWrapper(const EmpiricalFormula&, double, bool)))
    ptr = new IsoSpecThresholdWrapper(EmpiricalFormula("C10"), 0.5, false);
    TEST_NOT_EQUAL(ptr, nullPointer)
  END_SECTION

  START_SECTION((IsoSpecThresholdWrapper(std::vector<int>, std::vector<int>,
                                         std::vector<std::vector<double> >,
                                         std::vector<std::vector<double> >, 
                                         double, bool)))
    ptr2 = new IsoSpecThresholdWrapper(fructose_isotopeNumbers, fructose_atomCounts, fructose_isotopeMasses,
                                       fructose_isotopeProbabilities, 0.5, false);
    TEST_NOT_EQUAL(ptr2, nullPointer)
    TEST_EXCEPTION(Exception::IllegalArgument, IsoSpecThresholdWrapper(invalid_isotopeNumbers, invalid_atomCounts, invalid_isotopeMasses, invalid_isotopeProbabilities, 0.5, false));
  END_SECTION

  START_SECTION((IsoSpecThresholdWrapper::~IsoSpecThresholdWrapper()))
    delete ptr;
    delete ptr2;
  END_SECTION
}


START_SECTION(( void IsoSpecThresholdWrapper::run() ))
{
  {
  double threshold = 1e-5;
  bool absolute = false;


  IsotopeDistribution iso_result(IsoSpecThresholdWrapper(EmpiricalFormula("C6H12O6"), threshold, absolute).run());
  TEST_EQUAL(iso_result.size(), 14);
  TEST_EQUAL(compare_to_reference(iso_result, fructose_expected_oms), true);



  IsotopeDistribution iso_expl(IsoSpecThresholdWrapper(fructose_isotopeNumbers, fructose_atomCounts, fructose_isotopeMasses,
                   fructose_isotopeProbabilities, threshold, absolute).run());
  TEST_EQUAL(iso_expl.size(), 14);
  TEST_EQUAL(compare_to_reference(iso_expl, fructose_expected_oms), true);


  // human insulin
  IsotopeDistribution iso_result2 = IsoSpecThresholdWrapper(EmpiricalFormula("C520H817N139O147S8"), threshold, absolute).run();
  TEST_EQUAL(iso_result2.size(), 5513)

  IsotopeDistribution iso_result3 = IsoSpecThresholdWrapper(EmpiricalFormula("C520H817N139O147S8"), 0.01, absolute).run();
  TEST_EQUAL(iso_result3.size(), 267)
  }

  {
  double threshold = 1e-5;
  bool absolute = true;
  IsotopeDistribution iso_result(IsoSpecThresholdWrapper(EmpiricalFormula("C6H12O6"), threshold, absolute).run());

  TEST_EQUAL(iso_result.size(), 14)

  TEST_EQUAL(compare_to_reference(iso_result, fructose_expected_oms), true);

  // human insulin
  IsotopeDistribution iso_result2(IsoSpecThresholdWrapper(EmpiricalFormula("C520H817N139O147S8"), threshold, absolute).run());
  TEST_EQUAL(iso_result2.size(), 1734)

  IsotopeDistribution iso_result3(IsoSpecThresholdWrapper(EmpiricalFormula("C520H817N139O147S8"), 0.01, absolute).run());
  TEST_EQUAL(iso_result3.size(), 21)
  }
}
END_SECTION




// ----------------------------------------------------------------------------------------------------------------------
// Tests: IsoSpecTotalProbWrapper
// ----------------------------------------------------------------------------------------------------------------------


{
  IsoSpecWrapper* ptr = nullptr;
  IsoSpecWrapper* ptr2 = nullptr;
  IsoSpecWrapper* nullPointer = nullptr;
  START_SECTION((IsoSpecTotalProbWrapper::IsoSpecTotalProbWrapper(const EmpiricalFormula&, double, bool)))
    ptr = new IsoSpecTotalProbWrapper(EmpiricalFormula("C10"), 0.5, true);
    TEST_NOT_EQUAL(ptr, nullPointer);
  END_SECTION

  START_SECTION((IsoSpecTotalProbWrapper::IsoSpecTotalProbWrapper(std::vector<int>, std::vector<int>,
                                         std::vector<std::vector<double> >,
                                         std::vector<std::vector<double> >,
                                         double, bool)))
    ptr2 = new IsoSpecTotalProbWrapper(fructose_isotopeNumbers, fructose_atomCounts, fructose_isotopeMasses,
                                    fructose_isotopeProbabilities, 0.5, false);
    TEST_NOT_EQUAL(ptr2, nullPointer);

    TEST_EXCEPTION(Exception::IllegalArgument, IsoSpecTotalProbWrapper(invalid_isotopeNumbers, invalid_atomCounts,
                                                          invalid_isotopeMasses, invalid_isotopeProbabilities, 0.5, false));

  END_SECTION


  START_SECTION((IsoSpecTotalProbWrapper::~IsoSpecTotalProbWrapper()))
    delete ptr;
    delete ptr2;
  END_SECTION
}



START_SECTION(( void IsoSpecTotalProbWrapper::run() ))
{
  double total_prob = 0.99999;
  bool do_trim = true; // With do_trim == false the size of results is actually undefined, and may change as the underlying
                       // non-trimming heuristic changes

  IsotopeDistribution iso_result(IsoSpecTotalProbWrapper(EmpiricalFormula("C6H12O6"), total_prob, do_trim).run());
  TEST_EQUAL(iso_result.size(), 17);
  TEST_EQUAL(compare_to_reference(iso_result, fructose_expected_oms), true);

  IsotopeDistribution iso_result2(IsoSpecTotalProbWrapper(fructose_isotopeNumbers, fructose_atomCounts, fructose_isotopeMasses,
                                       fructose_isotopeProbabilities, total_prob, do_trim).run());
  TEST_EQUAL(iso_result2.size(), 17);
  TEST_EQUAL(compare_to_reference(iso_result2, fructose_expected_oms), true);

  // human insulin
  IsotopeDistribution iso_result3 = IsoSpecTotalProbWrapper(EmpiricalFormula("C520H817N139O147S8"), total_prob, do_trim).run();
  TEST_EQUAL(iso_result3.size(), 19615);

  IsotopeDistribution iso_result4 = IsoSpecTotalProbWrapper(EmpiricalFormula("C520H817N139O147S8"), 0.99, do_trim).run();
  TEST_EQUAL(iso_result4.size(), 1756);
}
END_SECTION



#if 0
START_SECTION(( [STRESSTEST] void run(const std::string&) ))
{
  // Do some stress testing of the library...
  // this is close to the performance of IsoSpec by itself
  int sum = 0;
  for (Size k = 0; k < 2e5; k++)
  {
    double threshold = 1e-2;
    bool absolute = false;
    IsoSpecThresholdWrapper iso("C520H817N139O147", threshold, absolute);
    auto res = iso.run();
    sum += res.size();
  }
  TEST_EQUAL(sum, 140*2*1e5)
}
END_SECTION
#endif

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
