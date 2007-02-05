// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Clemens Groepl, Marcel Grunert $
// --------------------------------------------------------------------------


#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/ExtendedModelFitter.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeaFiTraits.h>
#include <OpenMS/CONCEPT/Exception.h>

#include <OpenMS/KERNEL/DimensionDescription.h>

///////////////////////////

START_TEST(ExtendedModelFitter, "$Id: ExtendedModelFitter_test.C")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

enum DimensionId
{
	RT = DimensionDescription < LCMS_Tag >::RT,
	MZ = DimensionDescription < LCMS_Tag >::MZ
};

// default ctor
ExtendedModelFitter* ptr = 0;
CHECK(ExtendedModelFitter())
	ptr = new ExtendedModelFitter();
  TEST_EQUAL(ptr->getName(), "ExtendedModelFitter")
	TEST_NOT_EQUAL(ptr, 0)
RESULT

// destructor
CHECK(~ExtendedModelFitter())
	delete ptr;
RESULT

CHECK(void ExtendedModelFitter::setParameters(const Param& param))
	ExtendedModelFitter* fitter = new ExtendedModelFitter();
	Param p1;

	// change default settings
	p1.setValue("quality:minimum",0.0f);
	p1.setValue("isotope_model:stdev:first",0.08f);
	p1.setValue("isotope_model:stdev:last",0.12f);
	p1.setValue("isotope_model:stdev:step",0.02f);
	fitter->setParameters(p1);
	Param p2 = fitter->getParameters();
	// check changes
	TEST_EQUAL(p2.getValue("quality:minimum"),DataValue(0.0f))
	TEST_EQUAL(p2.getValue("isotope_model:stdev:first"),DataValue(0.08f))
	TEST_EQUAL(p2.getValue("isotope_model:stdev:last"),DataValue(0.12f))
	TEST_EQUAL(p2.getValue("isotope_model:stdev:step"),DataValue(0.02f))
	// check defaults
	TEST_EQUAL(p2.getValue("intensity_cutoff_factor"),DataValue(0.05f))
	TEST_EQUAL(p2.getValue("mz:interpolation_step"),DataValue(0.03f))
	TEST_EQUAL(p2.getValue("rt:interpolation_step"),DataValue(0.2f))
	TEST_EQUAL(p2.getValue("rt:max_iteration"),DataValue(500))
	TEST_EQUAL(p2.getValue("rt:deltaAbsError"),DataValue(0.0001))
	TEST_EQUAL(p2.getValue("rt:deltaRelError"),DataValue(0.0001))
	TEST_EQUAL(p2.getValue("rt:profile"),DataValue("EMG"))
	TEST_EQUAL(p2.getValue("min_num_peaks:final"),DataValue(5))
	TEST_EQUAL(p2.getValue("min_num_peaks:extended"),DataValue(10))
	TEST_EQUAL(p2.getValue("quality:type"),DataValue("Correlation"))
	TEST_EQUAL(p2.getValue("tolerance_stdev_bounding_box"),DataValue(3.0f))

	Param p3 = fitter->getParameters();
	TEST_EQUAL(p3.getValue("quality:minimum"),DataValue(0.0f))
	TEST_EQUAL(p3.getValue("isotope_model:stdev:first"),DataValue(0.08f))
	TEST_EQUAL(p3.getValue("isotope_model:stdev:last"),DataValue(0.12f))
	TEST_EQUAL(p3.getValue("isotope_model:stdev:step"),DataValue(0.02f))
RESULT


CHECK( DFeature<2> fit(const IndexSet& set) throw (UnableToFit))

	// Test EMG Fitting (mz/rt)
	const double default_precision = 0.1;
	PRECISION(default_precision)

	FeaFiTraits traits;
	double mzs[] = {675, 675.5, 676, 676.5, 677, 677.5, 678};
	const Size mz_num = 7;
	double rts[] = { 1260, 1260.5, 1261, 1261.5, 1262, 1262.5, 1263, 1263.5, 1264, 1264.5, 1265};
	const Size rt_num = 11;

	// Samples of Gaussian distribution N(mean,stdev) with scaling factor 20000
	double mean[2];	mean[MZ] = 676.5; mean[RT] = 1262.5;
	double stdev[2]; stdev[MZ] = 0.5; stdev[RT] = 0.9;

	double intens[] = {4.95329, 9.80589, 19.4003, 36.7884, 62.005, 77.2534, 62.0497, 36.7776, 19.3924, 9.80986, 4.95027, 60.9693, 120.699, 238.795, 452.823, 763.211, 950.901, 763.761, 452.69, 238.699, 120.748, 60.9322, 274.564, 543.548, 1075.37, 2039.21, 3436.98, 4282.21, 3439.46, 2038.61, 1074.94, 543.767, 274.397, 453.538, 897.857, 1776.35, 3368.46, 5677.37, 7073.55, 5681.46, 3367.46, 1775.63, 898.219, 453.262, 274.566, 543.55, 1075.38, 2039.22, 3437, 4282.23, 3439.48, 2038.62, 1074.94, 543.77, 274.398, 60.9465, 120.654, 238.706, 452.654, 762.926, 950.545, 763.476, 452.52, 238.61, 120.703, 60.9094, 4.95376, 9.80683, 19.4021, 36.7919, 62.011, 77.2608, 62.0557, 36.7811, 19.3943, 9.81079, 4.95074};	

	DPeak<2> p;
	DPeakArray<2> peak_array;
	for (Size mz=0; mz<mz_num; mz++) for (Size rt=0; rt<rt_num; rt++)
	{
		p.getPosition()[MZ] = mzs[mz];
		p.getPosition()[RT] = rts[rt];
		p.getIntensity() = intens[mz*rt_num+rt];
		peak_array.push_back(p);
	}
	
	peak_array.sortByPosition();
	MSExperimentExtern<DPeak<1> > exp;
	exp.set2DData(peak_array);
	traits.setData(exp.begin(), exp.end(),100);
	
	ExtendedModelFitter fitter;
	fitter.setTraits(&traits);
	Param param = fitter.getParameters();
	param.setValue("intensity_cutoff_factor",0.0f);
	fitter.setParameters(param);
	FeaFiModule::IndexSet  set;
	for (UnsignedInt i=0; i<exp.size(); ++i) 
	{
		for (UnsignedInt j=0; j<exp[i].size(); ++j) 
		{
			set.insert(std::make_pair(i,j));
		}
	}
	DFeature<2> feature = fitter.fit(set);

	TEST_REAL_EQUAL(feature.getPosition()[MZ], mean[MZ]);
	TEST_REAL_EQUAL(feature.getPosition()[RT], mean[RT]);
	TEST_REAL_EQUAL(feature.getIntensity(), 78602.6);
	TEST_EQUAL(feature.getCharge(), 0);
	PRECISION(0.01)
	TEST_REAL_EQUAL(feature.getOverallQuality(), 0.99);
		
	ProductModel<2>* model = dynamic_cast< ProductModel<2>* >
		(feature.getModelDescription().createModel());

	BaseModel<1>* mz_model = model->getModel(MZ);
	TEST_REAL_EQUAL(mz_model->getParameters().getValue("statistics:mean"),mean[MZ]);
	PRECISION(stdev[MZ]*stdev[MZ]*0.05)		// Variances can differ by 5%
	TEST_REAL_EQUAL(mz_model->getParameters().getValue("statistics:variance"),stdev[MZ]*stdev[MZ]);
	PRECISION(default_precision)

	BaseModel<1>* rt_model = model->getModel(RT);
	TEST_REAL_EQUAL(rt_model->getParameters().getValue("statistics:mean"),mean[RT]);
	PRECISION(stdev[RT]*stdev[RT])
	TEST_REAL_EQUAL(rt_model->getParameters().getValue("statistics:variance"),stdev[RT]*stdev[RT]);
	PRECISION(default_precision)
	
	// test predicted intensities
	DPosition<2> pos;
	for (Size mz=0; mz<mz_num; mz++) for (Size rt=0; rt<rt_num; rt++)
	{
	 	pos[MZ] = mzs[mz];
		pos[RT] = rts[rt];
		PRECISION(intens[mz*rt_num+rt]*0.08)		// Intensities can differ by 8%
		TEST_REAL_EQUAL(model->getIntensity(pos),intens[mz*rt_num+rt])
	}

RESULT


CHECK( DFeature<2> fit(const IndexSet& set) throw (UnableToFit))

	// Test Isotope/Bigauss Fitting (mz/rt)
	const double default_precision = 0.1;
	PRECISION(default_precision)

	FeaFiTraits traits;
	double mzs[] = {338, 338.1, 338.2, 338.3, 338.4, 338.5, 338.6, 338.7, 338.8,
		338.9, 339, 339.1, 339.2, 339.3, 339.4,	339.5, 339.6, 339.7, 339.8, 339.9,
		340, 340.1, 340.2, 340.3, 340.4 };
	const Size mz_num = 25;
	double rts[] = { 1261.6, 1261.8, 1262, 1262.2, 1262.4, 1262.6, 1262.8, 1263};
	const Size rt_num = 8;

	// Samples of theoretical isotope distribution in mz (charge=2, monoMass=mean[MZ], stdev[2])
	// asymmetrical retention profile (bigaussian with stdev[0] and stdev[1]
	// scaling factor 20000
	double mean[2];	mean[MZ] = 338.5; mean[RT] = 1262.4;
	double stdev[3]; stdev[0] = 0.2; stdev[1] = 0.3; stdev[2] = 0.1;

	double intens[] = { 0.002340574, 0.210691772, 6.97715327, 84.99912758, 380.9396643, 628.0641208, 381.0115632, 87.38019912, 35.98454301, 130.2127941, 214.3397749, 130.0205003, 29.61635618, 9.799801456, 33.32034304, 54.81824895, 33.25192853, 7.534121353, 2.014721947, 6.318548333, 10.38741682, 6.300717685, 1.424225194, 0.340398214, 1.011894924, 0.01108898, 0.998198173, 33.05578366, 402.7018848, 1804.784651, 2975.590602, 1805.125288, 413.98273, 170.4846121, 616.9114803, 1015.48138, 616.0004463, 140.3139396, 46.42869438, 157.8623843, 259.7133971, 157.5382557, 35.69454129, 9.545184149, 29.93549928, 49.21265019, 29.85102271, 6.747577139, 1.6127107, 4.794072654, 0.033685347, 3.032258312, 100.4146044, 1223.300312, 5482.451686, 9039.046129, 5483.486448, 1257.568494, 517.8865237, 1874.011608, 3084.760056, 1871.244131, 426.2361132, 141.0379203, 479.5435813, 788.9396394, 478.5589655, 108.4304424, 28.99570921, 90.93601745, 149.4948313, 90.67940027, 20.4973295, 4.89898254, 14.56310685, 0.065610097, 5.906032735, 195.5809433, 2382.663661, 10678.35778, 17605.65784, 10680.37322, 2449.408965, 1008.705212, 3650.076202, 6008.29217, 3644.685893, 830.1945873, 274.7043585, 934.0233574, 1536.644592, 932.1055877, 211.1936637, 56.47592987, 177.1191767, 291.176172, 176.6193547, 39.92334641, 9.54191506, 28.36505895, 0.081937096, 7.375742301, 244.2510398, 2975.586818, 13335.65503, 21986.80589, 13338.17202, 3058.941616, 1259.720363, 4558.393536, 7503.448881, 4551.661855, 1036.787571, 343.0642274, 1166.454014, 1919.036861, 1164.059009, 263.748968, 70.52990115, 221.1950835, 363.6350331, 220.5708814, 49.85822601, 11.91640983, 35.42367178, 0.049697361, 4.473613844, 148.1457443, 1804.784636, 8088.483645, 13335.67188, 8090.010272, 1855.341876, 764.0590226, 2764.805439, 4551.0718, 2760.722468, 628.8434496, 208.0789721, 707.4901223, 1163.954693, 706.0374786, 159.9718356, 42.77854747, 134.1615999, 220.5557965, 133.7830022, 30.24054271, 7.227667916, 21.48554302, 0.01108898, 0.998198173, 33.05578366, 402.7018848, 1804.784651, 2975.590602, 1805.125288, 413.98273, 170.4846121, 616.9114803, 1015.48138, 616.0004463, 140.3139396, 46.42869438, 157.8623843, 259.7133971, 157.5382557, 35.69454129, 9.545184149, 29.93549928, 49.21265019, 29.85102271, 6.747577139, 1.6127107, 4.794072654, 0.000910239, 0.081937096, 2.713383956, 33.05578366, 148.1457456, 244.2513505, 148.1737067, 33.98177182, 13.99422915, 50.63917801, 83.35578764, 50.56439579, 11.51766954, 3.811099314, 12.9581336, 21.31857384, 12.9315275, 2.929986373, 0.783516428, 2.457255417, 4.039620323, 2.450321158, 0.55387486, 0.132379356, 0.393521447};

	DPeak<2> p;
	DPeakArray<2> peak_array;
	for (Size rt=0; rt<rt_num; rt++) for (Size mz=0; mz<mz_num; mz++)
	{
		p.getPosition()[MZ] = mzs[mz];
		p.getPosition()[RT] = rts[rt];
		p.getIntensity() = intens[rt*mz_num+mz];
		peak_array.push_back(p);
	}
	peak_array.sortByPosition();
	MSExperimentExtern<DPeak<1> > exp;
	exp.set2DData(peak_array);
	traits.setData(exp.begin(), exp.end(),100);

	ExtendedModelFitter fitter;
	fitter.setTraits(&traits);
	Param param;
	param.setValue("quality:minimum",0.0f);
	param.setValue("isotope_model:stdev:first",0.06f);
 	param.setValue("isotope_model:stdev:last",0.14f);
	param.setValue("isotope_model:stdev:step",0.02f);
	param.setValue("rt:interpolation_step",0.05f);
	param.setValue("intensity_cutoff_factor",0.0f);
	fitter.setParameters(param);
	FeaFiModule::IndexSet  set;
	for (UnsignedInt i=0; i<exp.size(); ++i) 
	{
		for (UnsignedInt j=0; j<exp[i].size(); ++j) 
		{
			set.insert(std::make_pair(i,j));
		}
	}
	DFeature<2> feature = fitter.fit(set);

	TEST_REAL_EQUAL(feature.getPosition()[MZ], mean[MZ]);
	TEST_REAL_EQUAL(feature.getPosition()[RT], mean[RT]);
	TEST_REAL_EQUAL(feature.getIntensity(), 249316.7855);
	TEST_EQUAL(feature.getCharge(), 2);
	TEST_REAL_EQUAL(feature.getOverallQuality(), 0.9);

	ProductModel<2>* model = dynamic_cast< ProductModel<2>* >
		(feature.getModelDescription().createModel());

	BaseModel<1>* rt_model = model->getModel(RT);
	PRECISION(mean[RT]*0.01)		// Mean can differ by 1%
	TEST_REAL_EQUAL(rt_model->getParameters().getValue("statistics:mean"),mean[RT]);
	PRECISION(stdev[1])		// Variances can differ by 15%
	TEST_REAL_EQUAL(sqrt((double)(rt_model->getParameters().getValue("statistics:variance"))),stdev[1]);
	PRECISION(default_precision)

	BaseModel<1>* mz_model = model->getModel(MZ);
	TEST_REAL_EQUAL(mz_model->getParameters().getValue("isotope:stdev"),stdev[2]);

	// test predicted intensities
	DPosition<2> pos;
	for (Size rt=0; rt<rt_num; rt++)
		for (Size mz=0; mz<mz_num; mz++)
			if(intens[rt*mz_num+mz]>1000.0)
	{
	 	pos[MZ] = mzs[mz];
		pos[RT] = rts[rt];
		PRECISION(intens[rt*mz_num+mz]*0.50)		// individual Intensities can differ by 50%
		TEST_REAL_EQUAL(model->getIntensity(pos),intens[rt*mz_num+mz])
	}

RESULT


// checked by other check-methods 
// It is not necessarily to test the methods again.
CHECK(AssymStatistics())
  // ???
RESULT

CHECK((template< typename ProbabilityIterator, typename CoordinateIterator > void update( ProbabilityIterator const probability_begin, ProbabilityIterator const probability_end, CoordinateIterator const coordinate_begin)))
 // ???
RESULT

CHECK(static BaseModelFitter* create())
 // ???
RESULT

CHECK(static const String getName())
 // ???
RESULT

CHECK(coordinate_type variance1() const throw())
// ???
RESULT
 
CHECK(coordinate_type variance2() const throw())
// ???
RESULT

CHECK(void setData(const IndexSet& set))
// ???
RESULT

CHECK((int evaluate(const gsl_vector* x, void* params, gsl_vector* f, gsl_matrix* J)))
// ???
RESULT

CHECK((int jacobian(const gsl_vector* x, void* /* params */, gsl_matrix* J)))
// ???
RESULT

CHECK((int residual(const gsl_vector* x, void* /* params */, gsl_vector* f)))
// ???
RESULT

CHECK(void optimize())
// ???
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
