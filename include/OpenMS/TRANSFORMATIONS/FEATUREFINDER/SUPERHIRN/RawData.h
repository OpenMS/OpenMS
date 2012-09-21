/*
 *  RawData.h
 *  PeakDetection
 *
 *  Created by Markus Mueller on 10/19/06.
 *  Copyright 2006 ETHZ, IMSB, Switzerland. All rights reserved.
 *
 *  Ported to OpenMS by Florian Zeller, florian.zeller@bsse.ethz.ch
 *  December 2010
 *
 */

#ifndef _RAWDATA_h_
#define _RAWDATA_h_
#include <vector>
#include <ostream>

#include <OpenMS/config.h>

using namespace std;

// Class for the storage of raw MS data
class OPENMS_DLLAPI RawData{

public:
  
  double LOW_INTENSITY_MS_SIGNAL_THRESHOLD;
   
		RawData();
		RawData(vector<double>&,vector<double>&);
		virtual ~RawData();
		
		friend ostream& operator<<(ostream& pOut, RawData& pRawData);
  
		void get(vector<double>&,vector<double>&);
		void set(vector<double>&,vector<double>&);
		
		// Virtual functions
		virtual void smooth() {};
  
protected:
		vector<double> fProfileMasses;
		vector<double> fProfileIntens;
};

ostream& operator<<(ostream& pOut, RawData& pRawData);

#endif
