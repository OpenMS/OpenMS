// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Clemens Groepl $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/ModelFitter.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinder_impl.h>

///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(ModelFitter, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

typedef Peak1D PeakType;
typedef Feature FeatureType;
typedef ModelFitter<PeakType,FeatureType> ModelFitterType;
typedef FeatureFinderDefs::ChargedIndexSet ChargedIndexSet;

enum
{
	RT = Peak2D::RT,
	MZ = Peak2D::MZ
};

ModelFitterType* ptr = 0;
ModelFitterType* nullPointer = 0;
START_SECTION((ModelFitter(const MSExperiment<PeakType>* map, FeatureMap<FeatureType>* features, FeatureFinder* ff)))
	MSExperiment<PeakType> input;
	FeatureMap<FeatureType> features;
	FeatureFinder ff;
	ptr = new ModelFitterType(&input,&features,&ff);
  TEST_EQUAL(ptr->getName(), "ModelFitter")
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((virtual ~ModelFitter()))
	delete ptr;
END_SECTION

START_SECTION((void setMonoIsotopicMass(CoordinateType mz)))
{
  // dummy subtest
	TEST_EQUAL(1,1)
}
END_SECTION

START_SECTION([EXTRA](static const String getName()))
	MSExperiment<PeakType> input;
	FeatureMap<FeatureType> features;
	FeatureFinder ff;
	ModelFitterType model_fitter(&input,&features,&ff);
	TEST_EQUAL(model_fitter.getName(),"ModelFitter");
END_SECTION

#if 1
START_SECTION(([EXTRA]void ModelFitter::setParameters(const Param& param)))
{
	MSExperiment<PeakType> input;
	FeatureMap<FeatureType> features;
	FeatureFinder ff;
	ModelFitterType * model_fitter = new ModelFitterType(&input,&features,&ff);
	
	Param p1;

	// change default settings
	p1.setValue("quality:minimum",0.0f);
	p1.setValue("isotope_model:stdev:first",0.08f);
	p1.setValue("isotope_model:stdev:last",0.12f);
	p1.setValue("isotope_model:stdev:step",0.02f);
	model_fitter->setParameters(p1);
	Param p2 = model_fitter->getParameters();
	// check changes
	TEST_EQUAL(p2.getValue("quality:minimum"),DataValue(0.0f))
	TEST_EQUAL(p2.getValue("isotope_model:stdev:first"),DataValue(0.08f))
	TEST_EQUAL(p2.getValue("isotope_model:stdev:last"),DataValue(0.12f))
	TEST_EQUAL(p2.getValue("isotope_model:stdev:step"),DataValue(0.02f))
	// check defaults
	TEST_EQUAL(p2.getValue("intensity_cutoff_factor"),DataValue(0.05f))
	TEST_EQUAL(p2.getValue("mz:interpolation_step"),DataValue(0.03f))
	TEST_EQUAL(p2.getValue("rt:interpolation_step"),DataValue(0.2f))

	TEST_EQUAL(p2.getValue("max_iteration"),DataValue(500))
	TEST_EQUAL(p2.getValue("deltaAbsError"),DataValue(0.0001))
	TEST_EQUAL(p2.getValue("deltaRelError"),DataValue(0.0001))
	
	TEST_EQUAL(p2.getValue("min_num_peaks:final"),DataValue(5))
	TEST_EQUAL(p2.getValue("min_num_peaks:extended"),DataValue(10))
	TEST_EQUAL(p2.getValue("quality:type"),DataValue("Correlation"))
	TEST_EQUAL(p2.getValue("tolerance_stdev_bounding_box"),DataValue(3.0f))

	Param p3 = model_fitter->getParameters();
	TEST_EQUAL(p3.getValue("quality:minimum"),DataValue(0.0f))
	TEST_EQUAL(p3.getValue("isotope_model:stdev:first"),DataValue(0.08f))
	TEST_EQUAL(p3.getValue("isotope_model:stdev:last"),DataValue(0.12f))
	TEST_EQUAL(p3.getValue("isotope_model:stdev:step"),DataValue(0.02f))
}
END_SECTION
#endif

START_SECTION(Feature fit(const ChargedIndexSet& index_set))
{
	// Test BiGauss Fitting (mz/rt)

	MSExperiment<PeakType> input;
	FeatureMap<FeatureType> features;
	FeatureFinder ff;

	const double default_precision = 0.1;
	TOLERANCE_ABSOLUTE(default_precision);

	double mzs[] = {675, 675.5, 676, 676.5, 677, 677.5, 678};
	const UInt mz_num = 7;
	double rts[] = { 1260, 1260.5, 1261, 1261.5, 1262, 1262.5,
											1263, 1263.5, 1264, 1264.5, 1265};
	const UInt rt_num = 11;

	// Samples of Gaussian distribution N(mean,stdev) with scaling factor 20000
	double mean[2];	mean[MZ] = 676.5; mean[RT] = 1262.5;
	double stdev[2]; stdev[MZ] = 0.5; stdev[RT] = 0.9;

	float intens[] = { 1.65879841f, 6.652431187f, 19.59411554f, 42.38668296f, 67.34288093f,
		78.58007608f, 67.34288093f, 42.38668296f, 19.59411554f, 6.652431187f, 1.65879841f,
		20.20830161f, 81.04320276f, 238.7051942f, 516.3755092f, 820.4042402f, 957.3013023f,
		820.4042402f, 516.3755092f, 238.7051942f, 81.04320276f, 20.20830161f, 90.56732447f,
		363.210436f, 1069.80246f, 2314.234476f, 3676.796717f, 4290.326784f, 3676.796717f,
		2314.234476f, 1069.80246f, 363.210436f, 90.56732447f, 149.3202743f, 598.8327716f,
		1763.806071f, 3815.527605f, 6062.012955f, 7073.553026f, 6062.012955f, 3815.527605f,
		1763.806071f, 598.8327716f, 149.3202743f, 90.56732447f, 363.210436f, 1069.80246f,
		2314.234476f, 3676.796717f, 4290.326784f, 3676.796717f, 2314.234476f, 1069.80246f,
		363.210436f, 90.56732447f, 20.20830161f, 81.04320276f, 238.7051942f, 516.3755092f,
		820.4042402f, 957.3013023f, 820.4042402f, 516.3755092f, 238.7051942f, 81.04320276f,
		20.20830161f, 1.65879841f, 6.652431187f, 19.59411554f, 42.38668296f, 67.34288093f,
		78.58007608f, 67.34288093f, 42.38668296f, 19.59411554f, 6.652431187f, 1.65879841f};

	
	Peak2D p;
	std::vector<Peak2D> peak_array;
	for (Size mz=0; mz<mz_num; mz++) 
	{
		for (Size rt=0; rt<rt_num; rt++)
		{
			p.setMZ(mzs[mz]);
			p.setRT(rts[rt]);
			p.setIntensity(intens[mz*rt_num+rt]);
			peak_array.push_back(p);
		}
	}
	
	std::sort(peak_array.begin(),peak_array.end(),Peak2D::PositionLess());
		
	input.set2DData(peak_array);
	input.updateRanges(-1);
	
	ModelFitterType model_fitter(&input,&features,&ff);
	
	Param param = model_fitter.getParameters();
	param.setValue("intensity_cutoff_factor",0.0f);
  param.setValue("mz:model_type:first",0);
	param.setValue("fit_algorithm","simplest");
	model_fitter.setParameters(param);
	ChargedIndexSet  set;

	for (Size mz=0; mz<mz_num; mz++) 
	{
		for (Size rt=0; rt<rt_num; rt++)
		{
			set.insert(std::make_pair(rt,mz));
		}
	}
	FeatureType feature = model_fitter.fit(set);

	TEST_REAL_SIMILAR(feature.getMZ(), mean[MZ]);
	TEST_REAL_SIMILAR(feature.getRT(), mean[RT]);
	TEST_REAL_SIMILAR(feature.getIntensity(), 79820.9);
	TEST_EQUAL(feature.getCharge(), 0);
	TOLERANCE_ABSOLUTE(0.01)
	TEST_REAL_SIMILAR(feature.getOverallQuality(), 0.99);

	ProductModel<2>* model = static_cast< ProductModel<2>* >
													(feature.getModelDescription().createModel());

	BaseModel<1>* mz_model = model->getModel(MZ);
	TEST_REAL_SIMILAR(mz_model->getParameters().getValue("statistics:mean"),mean[MZ]);
	TOLERANCE_ABSOLUTE(stdev[MZ]*stdev[MZ]*0.05)		// Variances can differ by 5%
	TEST_REAL_SIMILAR(mz_model->getParameters().getValue("statistics:variance"),stdev[MZ]*stdev[MZ]);
	TOLERANCE_ABSOLUTE(default_precision)

	BaseModel<1>* rt_model = model->getModel(RT);
	TEST_REAL_SIMILAR(rt_model->getParameters().getValue("statistics:mean"),mean[RT]);
	TOLERANCE_ABSOLUTE(stdev[RT]*stdev[RT]*0.05)		// Variances can differ by 5%
	TEST_REAL_SIMILAR(rt_model->getParameters().getValue("statistics:variance1"),stdev[RT]*stdev[RT]);
	TEST_REAL_SIMILAR(rt_model->getParameters().getValue("statistics:variance2"),stdev[RT]*stdev[RT]);
	TOLERANCE_ABSOLUTE(default_precision)

	// test predicted intensities
	DPosition<2> pos;
	for (Size mz=0; mz<mz_num; mz++) for (Size rt=0; rt<rt_num; rt++)
	{
	 	pos[MZ] = mzs[mz];
		pos[RT] = rts[rt];
		TOLERANCE_ABSOLUTE(intens[mz*rt_num+rt]*0.1);		// Intensities can differ by 10%
		TEST_REAL_SIMILAR(model->getIntensity(pos),intens[mz*rt_num+rt]);
	}
}
END_SECTION

#if 1
START_SECTION(([EXTRA]Feature fit(const ChargedIndexSet& index_set) throw (UnableToFit)))
{
	// Test Isotope/Bigauss Fitting (mz/rt)

	const double default_precision = 0.1;
	TOLERANCE_ABSOLUTE(default_precision)

	double mzs[] = {338, 338.1, 338.2, 338.3,	338.4, 338.5, 338.6, 338.7, 338.8,
		338.9, 339, 339.1, 339.2, 339.3, 339.4,	339.5, 339.6, 339.7, 339.8, 339.9,
		340, 340.1, 340.2, 340.3, 340.4 };
	const UInt mz_num = sizeof(mzs) / sizeof(*mzs); // currently 25
	
	double rts[] = { 1261.6, 1261.8, 1262, 1262.2, 1262.4, 1262.6, 1262.8, 1263};
	const UInt rt_num = sizeof(rts) / sizeof(*rts); // currently 8

	// Samples of theoretical isotope distribution in mz (charge=2, monoMass=mean[MZ], stdev[2])
	// asymmetrical retention profile (bigaussian with stdev[0] and stdev[1]
	// scaling factor 20000
	double mean[2];	mean[MZ] = 338.5; mean[RT] = 1262.4;
	double stdev[3]; stdev[0] = 0.231; stdev[1] = 0.3; stdev[2] = 0.1;

	double intens[] = { 0.002340574, 0.210691772, 6.97715327, 84.99912758, 380.9396643, 628.0641208, 381.0115632, 87.38019912, 35.98454301, 130.2127941, 214.3397749, 130.0205003, 29.61635618, 9.799801456, 33.32034304, 54.81824895, 33.25192853, 7.534121353, 2.014721947, 6.318548333, 10.38741682, 6.300717685, 1.424225194, 0.340398214, 1.011894924, 0.01108898, 0.998198173, 33.05578366, 402.7018848, 2814.6645, 4522.9635, 2717.9924, 413.98273, 170.4846121, 616.9114803, 1015.48138, 616.0004463, 140.3139396, 46.42869438, 157.8623843, 259.7133971, 157.5382557, 35.69454129, 9.545184149, 29.93549928, 49.21265019, 29.85102271, 6.747577139, 1.6127107, 4.794072654, 0.033685347, 3.032258312, 100.4146044, 1223.300312, 5482.451686, 9039.046129, 5483.486448, 1257.568494, 517.8865237, 1874.011608, 3084.760056, 1871.244131, 426.2361132, 141.0379203, 479.5435813, 788.9396394, 478.5589655, 108.4304424, 28.99570921, 90.93601745, 149.4948313, 90.67940027, 20.4973295, 4.89898254, 14.56310685, 0.065610097, 5.906032735, 195.5809433, 2382.663661, 10678.35778, 17605.65784, 10680.37322, 2449.408965, 1008.705212, 3650.076202, 6008.29217, 3644.685893, 830.1945873, 274.7043585, 934.0233574, 1536.644592, 932.1055877, 211.1936637, 56.47592987, 177.1191767, 291.176172, 176.6193547, 39.92334641, 9.54191506, 28.36505895, 0.081937096, 7.375742301, 244.2510398, 2975.586818, 13335.65503, 21986.80589, 13338.17202, 3058.941616, 1259.720363, 4558.393536, 7503.448881, 4551.661855, 1036.787571, 343.0642274, 1166.454014, 1919.036861, 1164.059009, 263.748968, 70.52990115, 221.1950835, 363.6350331, 220.5708814, 49.85822601, 11.91640983, 35.42367178, 0.049697361, 4.473613844, 148.1457443, 1804.784636, 8088.483645, 13335.67188, 8090.010272, 1855.341876, 764.0590226, 2764.805439, 4551.0718, 2760.722468, 628.8434496, 208.0789721, 707.4901223, 1163.954693, 706.0374786, 159.9718356, 42.77854747, 134.1615999, 220.5557965, 133.7830022, 30.24054271, 7.227667916, 21.48554302, 0.01108898, 0.998198173, 33.05578366, 402.7018848, 1804.784651, 2975.590602, 1805.125288, 413.98273, 170.4846121, 616.9114803, 1015.48138, 616.0004463, 140.3139396, 46.42869438, 157.8623843, 259.7133971, 157.5382557, 35.69454129, 9.545184149, 29.93549928, 49.21265019, 29.85102271, 6.747577139, 1.6127107, 4.794072654, 0.000910239, 0.081937096, 2.713383956, 33.05578366, 148.1457456, 244.2513505, 148.1737067, 33.98177182, 13.99422915, 50.63917801, 83.35578764, 50.56439579, 11.51766954, 3.811099314, 12.9581336, 21.31857384, 12.9315275, 2.929986373, 0.783516428, 2.457255417, 4.039620323, 2.450321158, 0.55387486, 0.132379356, 0.393521447};

	Peak2D p;
	std::vector<Peak2D> peak_array;
	for (Size rt=0; rt<rt_num; rt++) for (Size mz=0; mz<mz_num; mz++)
	{
		p.setMZ(mzs[mz]);
		p.setRT(rts[rt]);
		p.setIntensity(intens[rt*mz_num+mz]);
		peak_array.push_back(p);
	}
	std::sort(peak_array.begin(),peak_array.end(),Peak2D::PositionLess());

	MSExperiment<PeakType> input;
	FeatureMap<FeatureType> features;
	FeatureFinder ff;

	input.set2DData(peak_array);
	input.updateRanges(-1);
	
	ModelFitterType model_fitter(&input,&features,&ff);

	Param param = model_fitter.getParameters();
	param.setValue("quality:minimum",0.0f);
	param.setValue("isotope_model:stdev:first",0.06f);
 	param.setValue("isotope_model:stdev:last",0.14f);
	param.setValue("isotope_model:stdev:step",0.02f);
	param.setValue("rt:interpolation_step",0.05f);
	param.setValue("intensity_cutoff_factor",0.0f);
	param.setValue("fit_algorithm","simplest");
	param.setValue("mz:model_type:first", 0);
	param.setValue("mz:model_type:last", 4);

	model_fitter.setParameters(param);
	ChargedIndexSet  set;
	
	for (Size i=0; i<input.size(); ++i) 
	{
		for (Size j=0; j<input[i].size(); ++j) 
		{
			set.insert(std::make_pair(i,j));
		}
	}
	Feature feature = model_fitter.fit(set);

	TOLERANCE_ABSOLUTE(2.0);
		
	TEST_REAL_SIMILAR(feature.getMZ(), mean[MZ]);
	TEST_REAL_SIMILAR(feature.getRT(), mean[RT]);
	TEST_REAL_SIMILAR(feature.getIntensity(), 252787);
	TEST_EQUAL(feature.getCharge(), 2);
	TEST_REAL_SIMILAR(feature.getOverallQuality(), 0.9);

	ProductModel<2>* model = static_cast< ProductModel<2>* >
													(feature.getModelDescription().createModel());

	// std::cout << model->getParameters() << std::endl;
	
	BaseModel<1>* rt_model = model->getModel(RT);
	TOLERANCE_ABSOLUTE(mean[RT]*0.01)		// Mean can differ by 1%
	TEST_REAL_SIMILAR(rt_model->getParameters().getValue("statistics:mean"),mean[RT]);
	TOLERANCE_ABSOLUTE(stdev[1]*0.15)		// Variances can differ by 15%
	TEST_REAL_SIMILAR(sqrt(double(rt_model->getParameters().getValue("statistics:variance1"))),stdev[1]);
	TOLERANCE_ABSOLUTE(stdev[0]*0.15)		// Variances can differ by 15%
	TEST_REAL_SIMILAR(sqrt(double(rt_model->getParameters().getValue("statistics:variance2"))),stdev[0]);
	TOLERANCE_ABSOLUTE(default_precision)

	BaseModel<1>* mz_model = model->getModel(MZ);
	TEST_REAL_SIMILAR(mz_model->getParameters().getValue("isotope:mode:GaussianSD"),stdev[2]);

	// test predicted intensities
	DPosition<2> pos;
	for (Size rt=0; rt<rt_num; rt++)
		for (Size mz=0; mz<mz_num; mz++)
			if(intens[rt*mz_num+mz]>1000.0)
	{
	 	pos[MZ] = mzs[mz];
		pos[RT] = rts[rt];
		TOLERANCE_ABSOLUTE(intens[rt*mz_num+mz]*0.50)		// individual Intensities can differ by 50%
		TEST_REAL_SIMILAR(model->getIntensity(pos),intens[rt*mz_num+mz])
	}
}
END_SECTION
#endif

// "monster test" (?)
#if 0
START_SECTION(([EXTRA]Feature fit(const ChargedIndexSet& index_set) throw (UnableToFit)))
{
	// Test Isotope/Bigauss Fitting (mz/rt)

	const double default_precision = 0.1;
	TOLERANCE_ABSOLUTE(default_precision)

	double mzs[] = {675, 675.5, 676, 676.5, 677, 677.5, 678};
	const UInt mz_num = sizeof(mzs) / sizeof(*mzs); // currently 7 
	STATUS(mz_num);

	double rts[] = { 1260, 1260.5, 1261, 1261.5, 1262, 1262.5, 1263, 1263.5, 1264, 1264.5, 1265};
	const UInt rt_num = sizeof(rts) / sizeof(*rts); // currently 11 
	STATUS(rt_num);

	// Samples of Gaussian distribution N(mean,stdev) with scaling factor 20000
	double mean[2];	mean[MZ] = 676.5; mean[RT] = 1262.5;
	double stdev[2]; stdev[MZ] = 0.5; stdev[RT] = 0.9;

	double intens[] = {4.95329, 9.80589, 19.4003, 36.7884, 62.005, 77.2534, 62.0497, 36.7776, 19.3924, 9.80986, 4.95027, 60.9693, 120.699, 238.795, 452.823, 763.211, 950.901, 763.761, 452.69, 238.699, 120.748, 60.9322, 274.564, 543.548, 1075.37, 2039.21, 3436.98, 4282.21, 3439.46, 2038.61, 1074.94, 543.767, 274.397, 453.538, 897.857, 1776.35, 3368.46, 5677.37, 7073.55, 5681.46, 3367.46, 1775.63, 898.219, 453.262, 274.566, 543.55, 1075.38, 2039.22, 3437, 4282.23, 3439.48, 2038.62, 1074.94, 543.77, 274.398, 60.9465, 120.654, 238.706, 452.654, 762.926, 950.545, 763.476, 452.52, 238.61, 120.703, 60.9094, 4.95376, 9.80683, 19.4021, 36.7919, 62.011, 77.2608, 62.0557, 36.7811, 19.3943, 9.81079, 4.95074};	

	STATUS(sizeof(intens)/sizeof(*intens));
	
	Peak2D p;
	std::vector<Peak2D> peak_array;
	for (Size rt=0; rt<rt_num; rt++) for (Size mz=0; mz<mz_num; mz++)
	{
		p.setMZ(mzs[mz]);
		p.setRT(rts[rt]);
		p.setIntensity(intens[rt*mz_num+mz]);
		peak_array.push_back(p);
	}
	std::sort(peak_array.begin(),peak_array.end(),Peak2D::PositionLess());

	MSExperiment<PeakType> input;
	FeatureMap<FeatureType> features;
	FeatureFinder ff;

	input.set2DData(peak_array);
	input.updateRanges(-1);
	
	// need to run ff so that flags are initialized
	ff.run("none",input,features,Param());
	
	ModelFitterType model_fitter(&input,&features,&ff);

	Param param = model_fitter.getParameters();
	param.setValue("rt:max_iteration",500);
	param.setValue("rt:deltaAbsError",0.0001);
	param.setValue("rt:deltaRelError",0.0001);
	param.setValue("rt:profile","EMG");
	param.setValue("quality:minimum",0.3f);
	param.setValue("intensity_cutoff_factor",0.0f);
	model_fitter.setParameters(param);
	ChargedIndexSet  set;
	
	for (Size i=0; i<input.size(); ++i) 
	{
		for (Size j=0; j<input[i].size(); ++j) 
		{
			set.insert(std::make_pair(i,j));
		}
	}

	Feature feature = model_fitter.fit(set);

	TEST_REAL_SIMILAR(feature.getMZ(), mean[MZ]);
	TEST_REAL_SIMILAR(feature.getRT(), mean[RT]);
	TEST_REAL_SIMILAR(feature.getIntensity(), 78602.6);
	TEST_EQUAL(feature.getCharge(), 0);
	TOLERANCE_ABSOLUTE(0.01)
	TEST_REAL_SIMILAR(feature.getOverallQuality(), 0.99);
		
	ProductModel<2>* model = dynamic_cast< ProductModel<2>* > (feature.getModelDescription().createModel());

	BaseModel<1>* mz_model = model->getModel(MZ);
	TEST_REAL_SIMILAR(mz_model->getParameters().getValue("statistics:mean"),mean[MZ]);
	TOLERANCE_ABSOLUTE(stdev[MZ]*stdev[MZ]*0.05)		// Variances can differ by 5%
	TEST_REAL_SIMILAR(mz_model->getParameters().getValue("statistics:variance"),stdev[MZ]*stdev[MZ]);
	TOLERANCE_ABSOLUTE(default_precision)

	BaseModel<1>* rt_model = model->getModel(RT);
	TEST_REAL_SIMILAR(rt_model->getParameters().getValue("statistics:mean"),mean[RT]);
	TOLERANCE_ABSOLUTE(stdev[RT]*stdev[RT])
	TEST_REAL_SIMILAR(rt_model->getParameters().getValue("statistics:variance"),stdev[RT]*stdev[RT]);
	TOLERANCE_ABSOLUTE(default_precision)
	
	// test predicted intensities
	DPosition<2> pos;
	for (Size mz=0; mz<mz_num; mz++) for (Size rt=0; rt<rt_num; rt++)
	{
	 	pos[MZ] = mzs[mz];
		pos[RT] = rts[rt];
		TOLERANCE_ABSOLUTE(intens[mz*rt_num+rt]*0.08)		// Intensities can differ by 8%
		TEST_REAL_SIMILAR(model->getIntensity(pos),intens[mz*rt_num+rt])
	}

}
END_SECTION
#endif

// "monster test"
#if 0
START_SECTION(([EXTRA]void ExtendedModelFitter::optimize()))
{
	// *************************************************
	// *** check parameter optimization at EMG model ***
	// *************************************************
	EmgModel em1;
	Param tmp;

	// model parameter
	double height_ = 1000.0;
	double width_ = 2;
	double symmetry_ = 1.3;
	double retention_ = 700.0;
	double min_ = 650;
	double max_ = 750;

	// iterate symmetry from 1.3 to 8.3, i.e. fronted, symmetric and tailed peaks
	while (symmetry_ < 8.4) 
	{
		// set model parameter
		tmp.setValue("bounding_box:min", min_);
		tmp.setValue("bounding_box:max", max_);
		tmp.setValue("statistics:variance",  0.0);
		tmp.setValue("emg:height",height_);
		tmp.setValue("emg:width",width_);
		tmp.setValue("emg:symmetry",symmetry_);
		tmp.setValue("emg:retention",retention_);
		em1.setParameters(tmp);

		// get samples from model
		std::vector<Peak1D> dpa1;
		em1.getSamples(dpa1);

		// save samples
		Peak2D p;
		std::vector<Peak2D> peak_array;
		for (Size i=0; i<dpa1.size(); ++i) {
			p.setMZ(1);
			p.setRT(dpa1[i].getPosition()[0]);
			p.setIntensity(dpa1[i].getIntensity());
			peak_array.push_back(p);
		}

		// set traits
		std::sort(peak_array.begin(),peak_array.end(),Peak2D::PositionLess());
		MSExperiment<PeakType> input;
		FeatureMap<FeatureType> features;
		FeatureFinder ff;

		input.set2DData(peak_array);
		input.updateRanges(-1);

		// set traits and parameter in ExtendedModelFitter
		ModelFitterType model_fitter(&input,&features,&ff);
		Param param = model_fitter.getParameters();
		param.setValue("rt:max_iteration",10000);
		model_fitter.setParameters(param);
		ChargedIndexSet  set;

		// construct indexSet
		for (Size i=0; i<input.size(); ++i) 
			for (Size j=0; j<input[i].size(); ++j) 
				set.insert(std::make_pair(i,j));

		// compute start parameter
		model_fitter.setInitialParameters(set);
		// optimize parameter with Levenberg-Maruardt algorithm
		model_fitter.optimize();

		// test
		TEST_REAL_SIMILAR(model_fitter.getSymmetry(), symmetry_); 
		TEST_REAL_SIMILAR(model_fitter.getHeight(), height_);
		TEST_REAL_SIMILAR(model_fitter.getWidth(), width_);
		TEST_REAL_SIMILAR(model_fitter.getRT(), retention_);
		TEST_EQUAL(model_fitter.getGSLStatus(), "success");

		symmetry_ += 0.5;

		// clear vector
		dpa1.clear();
		peak_array.clear();
	}

	// ************************************************************
	// *** check parameter optimization with noise at EMG model ***
	// ************************************************************

	EmgModel em2;
	em2.setInterpolationStep(0.1);

	// set model parameter
	tmp.setValue("bounding_box:min", 650);
	tmp.setValue("bounding_box:max", 750);
	tmp.setValue("statistics:variance",  0.0);
	tmp.setValue("emg:height",1000.0);
	tmp.setValue("emg:width",2.0);
	tmp.setValue("emg:symmetry",1.3);
	tmp.setValue("emg:retention",700);
	em2.setParameters(tmp);

	//get samples from model
	std::vector<Peak1D> dpa2;
	em2.getSamples(dpa2);

	// save samples
	Peak2D p2;
	std::vector<Peak2D> peak_array2;

	// noise value		
	Real noise = 10;

	//String fname2 = "samples2.dta2d";
	//ofstream file2(fname2.c_str()); 	
	for (Size i=0; i<dpa2.size(); ++i)
	{

		//if (i%3 == 0 && (dpa2[i].getPosition()[1] > 1))
		//	file2 << dpa2[i].getPosition()[0] << "	" << noise << "\n";
		//else
		//	file2 << dpa2[i].getPosition()[0] << "	" << (dpa2[i].getPosition()[1]) << "\n";

		p2.setMZ(1);
		p2.setRT(dpa2[i].getPosition()[0]);

		if (i%5 == 0 && (dpa2[i].getIntensity() > 1))
			p2.setIntensity(noise);
		else
			p2.setIntensity(dpa2[i].getIntensity());

		peak_array2.push_back(p2);
	}
	//file2.close();

	// set traits
	std::sort(peak_array2.begin(),peak_array2.end(),Peak2D::PositionLess());
	MSExperiment<PeakType> input2;
	FeatureMap<FeatureType> features;
	FeatureFinder ff;

	input2.set2DData(peak_array2);
	input2.updateRanges(-1);

	// set traits and parameter in ExtendedModelFitter
	ModelFitterType fitter2(&input2,&features,&ff);
	Param param2 = fitter2.getParameters();
	param2.setValue("rt:max_iteration",10000);
	fitter2.setParameters(param2);
	ChargedIndexSet  set2;

	// construct indexSet
	for (Size i=0; i<input2.size(); ++i) 
		for (Size j=0; j<input2[i].size(); ++j) 
			set2.insert(std::make_pair(i,j));

	// compute start parameter
	fitter2.setInitialParameters(set2);
	// optimize parameter with Levenberg-Maruardt algorithm
	fitter2.optimize();
	// test
	TOLERANCE_ABSOLUTE(0.5)
		TEST_REAL_SIMILAR(fitter2.getSymmetry(), 1.3); 
	TEST_REAL_SIMILAR(fitter2.getWidth(), 2);
	TEST_REAL_SIMILAR(fitter2.getRT(), 700);
	TEST_EQUAL(fitter2.getGSLStatus(), "success");

	// clear vector
	dpa2.clear();
	peak_array2.clear();


	// *******************************************************
	// *** check parameter optimization at LogNormal model ***
	// *******************************************************
	LogNormalModel logm1;	
	logm1.setInterpolationStep(0.1);

	tmp.setValue("bounding_box:min", 1.0);
	tmp.setValue("bounding_box:max", 70.0);
	tmp.setValue("statistics:variance",  0.0);
	tmp.setValue("emg:height",  850.0);
	tmp.setValue("emg:width",  20.0);
	tmp.setValue("emg:symmetry",  1.5);
	tmp.setValue("emg:retention",  30.0);
	tmp.setValue("lognormal:r",  2.0);
	logm1.setParameters(tmp);

	// get samples from model
	std::vector<Peak1D> dpa3;
	logm1.getSamples(dpa3);

	Peak2D p3;
	std::vector<Peak2D> peak_array3;

	//String fname3 = "samples3.dta2d";
	//ofstream file3(fname3.c_str()); 
	// save samples
	for (Size i=0; i<dpa3.size(); ++i) {
		//file3 << dpa3[i].getPosition()[0] << "	" << (dpa3[i].getPosition()[1]) << "\n";
		p3.setMZ(1);
		p3.setRT(dpa3[i].getPosition()[0]);
		p3.setIntensity(dpa3[i].getIntensity());
		peak_array3.push_back(p3);
	}
	//file3.close();

	// set traits
	std::sort(peak_array3.begin(),peak_array3.end(),Peak2D::PositionLess());
	MSExperiment<PeakType > exp3;
	exp3.set2DData(peak_array3);

	// set traits and parameter in ExtendedModelFitter
	ModelFitterType fitter3(&exp3,&features,&ff);
	Param param3 = fitter3.getParameters();
	param3.setValue("rt:profile","LogNormal");
	fitter3.setParameters(param3);

	// construct indexSet
	ChargedIndexSet set3;
	for (Size i=0; i<exp3.size(); ++i) 
		for (Size j=0; j<exp3[i].size(); ++j) 
			set3.insert(std::make_pair(i,j));

	// compute start parameter
	fitter3.setInitialParameters(set3);
	// optimize parameter with Levenberg-Maruardt algorithm
	fitter3.optimize();

	// test
	TEST_EQUAL(fitter3.getGSLStatus(), "success");
	TOLERANCE_ABSOLUTE(50.0)
		TEST_REAL_SIMILAR(fitter3.getHeight(), 850.0);
	TOLERANCE_ABSOLUTE(1.0)
		TEST_REAL_SIMILAR(fitter3.getWidth(), 20.0);
	TEST_REAL_SIMILAR(fitter3.getSymmetry(), 1.5);
	TEST_REAL_SIMILAR(fitter3.getRT(), 30.0);

	// clear vector
	dpa3.clear();
	peak_array3.clear();


	// ******************************************************
	// *** check parameter optimization at LmaGauss model ***
	// ******************************************************
	LmaGaussModel lm1;
	lm1.setInterpolationStep(0.1);

	// model parameter
	double scale_factor = 1000.0;
	double standard_deviation = 0.5;
	double expected_value = 665.0;
	double min = 650;
	double max = 700;
	Param tmp2;

	// iterate expected_value
	while (expected_value < max)
	{
		// iterate standard_deviation from 0.5 to 9.5
		while (standard_deviation < 10) 
		{
			// set model parameter
			tmp2.setValue("bounding_box:min", min);
			tmp2.setValue("bounding_box:max", max);
			tmp2.setValue("statistics:variance", 0.0);
			tmp2.setValue("lma:scale_factor", scale_factor);
			tmp2.setValue("lma:standard_deviation", standard_deviation);
			tmp2.setValue("lma:expected_value", expected_value);	
			lm1.setParameters(tmp2);

			// get samples from model
			std::vector<Peak1D > dpa4;
			lm1.getSamples(dpa4);

			Peak2D p4;
			std::vector<Peak2D> peak_array4;

			//String fname4 = "samples4.dta2d";
			//ofstream file4(fname4.c_str()); 
			// save samples
			for (Size i=0; i<dpa4.size(); ++i)
			{

				//file4 << dpa4[i].getPosition()[0] << "	" << (dpa4[i].getPosition()[1]) << "\n";
				p4.setMZ(1);
				p4.setRT(dpa4[i].getPosition()[0]);
				p4.setIntensity(dpa4[i].getIntensity());
				peak_array4.push_back(p4);
			}
			//file4.close();

			// set traits
			std::sort(peak_array4.begin(),peak_array4.end(),Peak2D::PositionLess());
			MSExperiment<PeakType> exp4;
			exp4.set2DData(peak_array4);

			// set traits and parameter in ExtendedModelFitter
			ModelFitterType fitter4(&exp4,&features,&ff);
			Param param4 = fitter4.getParameters();
			param4.setValue("rt:profile","LmaGauss");
			fitter4.setParameters(param4);

			// construct indexSet
			ChargedIndexSet set4;
			for (Size i=0; i<exp4.size(); ++i) 
				for (Size j=0; j<exp4[i].size(); ++j) 
					set4.insert(std::make_pair(i,j));

			// compute start parameter
			fitter4.setInitialParameters(set4);
			// optimize parameter with Levenberg-Maruardt algorithm
			fitter4.optimize();

			// test
			TEST_REAL_SIMILAR(fitter4.getExpectedValue(), expected_value); 
			TEST_REAL_SIMILAR(fitter4.getStandardDeviation(), standard_deviation);
			TEST_REAL_SIMILAR(fitter4.getScaleFactor(), scale_factor);
			TEST_EQUAL(fitter4.getGSLStatus(), "success");

			standard_deviation += 2;

			// clear vector
			dpa4.clear();
			peak_array4.clear();

		}
		standard_deviation = 0.5;
		expected_value += 5;
	} 
}
END_SECTION
#endif 


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



