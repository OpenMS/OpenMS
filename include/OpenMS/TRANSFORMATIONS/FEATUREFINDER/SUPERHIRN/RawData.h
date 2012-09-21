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

namespace OpenMS
{

// Class for the storage of raw MS data
class OPENMS_DLLAPI OPENMS_DLLAPI RawData{

public:
  
  double LOW_INTENSITY_MS_SIGNAL_THRESHOLD;
   
		RawData();
    RawData(std::vector<double>&,std::vector<double>&);
		virtual ~RawData();
		
		friend std::ostream& operator<<(std::ostream& pOut, RawData& pRawData);
  
		void get(std::vector<double>&,std::vector<double>&);
		void set(std::vector<double>&,std::vector<double>&);
		
		// Virtual functions
		virtual void smooth() {};
  
protected:
		std::vector<double> fProfileMasses;
		std::vector<double> fProfileIntens;
};

std::ostream& operator<<(std::ostream& pOut, RawData& pRawData);

} // ns

#endif
