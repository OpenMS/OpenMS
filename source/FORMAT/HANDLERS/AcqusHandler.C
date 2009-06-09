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
// $Maintainer: Guillaume Belz
// $Authors: Guillaume Belz
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/HANDLERS/AcqusHandler.h>

#include <fstream>
#include <math.h>

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
        try
        {        
          if( line.empty() ) continue;
          if( line.prefix(2) != String("##") ) continue; 
        
          if(line.split('=', strings))
          {
            if( strings.size() != 2 )
            {
              //throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, std::string("Bad data line: \"")+line+"\"" , filename);
            }
            params_[ strings[0].substr(2) ] = strings[1].trim();
          } 
        }
				catch (...)
				{
					//throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, std::string("Bad data line: \"")+line+"\"" ,filename);
				}        
      }
     
      // TOF calibration params
      DW_ = params_[String("$DW")].toDouble();
      DELAY_ = (unsigned int) params_[String("$DELAY")].toInt();
      ML1_ = params_[String("$ML1")].toDouble();
      ML2_ = params_[String("$ML2")].toDouble();
      ML3_ = params_[String("$ML3")].toDouble();
      TD_ = (unsigned int) params_[String("$TD")].toInt();

      is.close();  	
    }
    
    AcqusHandler::~AcqusHandler()
    {
      params_.clear();
    }

    unsigned int AcqusHandler::getSize()
    {
      return TD_;
    }
        
    double AcqusHandler::getPosition(const unsigned int index)
    {
      double SqrtMZ_;
      double TOF_ = DW_ * index + DELAY_;
      double A_ = ML3_;
      double B_ = sqrt(1000000000000.0/ML1_);
      double C_ = ML2_ - TOF_;
      
      if(ML3_== 0.0) 
        SqrtMZ_ = C_ / B_;
      else 
        SqrtMZ_ = (sqrt(B_*B_-4*A_*C_)-B_) / (2*A_);
        
      return SqrtMZ_*SqrtMZ_;
    }

    String AcqusHandler::getParam(const String& param)
    {
      if( param == String("mzMax") )
        return String( getPosition(TD_-1) );
      else if( param == String("mzMin") )
        return String( getPosition(0) );
      else
        return params_[param];
    }

	} // namespace Internal
} // namespace OpenMS
