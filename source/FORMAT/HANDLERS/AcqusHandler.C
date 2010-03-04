// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                    Flex series file support
// --------------------------------------------------------------------------
//  Copyright (C) 2009 -- Guillaume Belz (guillaume.belz@chu-lyon.fr)
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
// $Maintainer: Guillaume Belz$
// $Authors: Guillaume Belz$
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/HANDLERS/AcqusHandler.h>

#include <fstream>
#include <cmath>

using namespace std;

namespace OpenMS
{
	namespace Internal
	{
	
    AcqusHandler::AcqusHandler(const String& filename)
    {  
      params_.clear();
      
      std::ifstream is(filename.c_str());
      if (!is)
      {
        throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
      }
      
      //temporary variables
      String line;
      std::vector<String> strings(2);

      //read lines
      while (getline(is, line, '\n'))
      {
        if (line.size() < 5) continue; // minimal string = "##x=x"
        if (line.prefix(2) != String("##")) continue; 
      
        if (line.split('=', strings))
        {       
          if (strings.size() == 2)
          {
            params_[strings[0].substr(2)] = strings[1].trim();
          }
        }      
      }
     
      // TOF calibration params
      dw_ = params_[String("$DW")].toDouble();
      delay_ = (Size)params_[String("$DELAY")].toInt();
      ml1_ = params_[String("$ML1")].toDouble();
      ml2_ = params_[String("$ML2")].toDouble();
      ml3_ = params_[String("$ML3")].toDouble();
      td_ = (Size) params_[String("$TD")].toInt();

      is.close();  	
    }
    
    AcqusHandler::~AcqusHandler()
    {
      params_.clear();
    }

    Size AcqusHandler::getSize()
    {
      return td_;
    }
        
    DoubleReal AcqusHandler::getPosition(const Size index)
    {
      DoubleReal sqrt_mz_;
      DoubleReal tof_ = dw_ * index + delay_;
      DoubleReal a_ = ml3_;
      DoubleReal b_ = sqrt(1000000000000.0/ml1_);
      DoubleReal c_ = ml2_ - tof_;
      
      if(ml3_== 0.0)
			{
        sqrt_mz_ = c_ / b_;
			}
      else
			{
        sqrt_mz_ = (sqrt(b_ * b_ - 4 * a_ * c_) - b_) / (2 * a_);
			}
      return sqrt_mz_* sqrt_mz_;
    }

    String AcqusHandler::getParam(const String& param)
    {
      if (param == String("mzMax"))
			{
        return String(getPosition(td_ - 1));
			}
      else if (param == String("mzMin"))
			{
        return String(getPosition(0));
			}
      return params_[param];
    }

	} // namespace Internal
} // namespace OpenMS

