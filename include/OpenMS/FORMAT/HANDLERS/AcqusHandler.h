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

#ifndef OPENMS_FORMAT_HANDLERS_ACQUSHANDLER_H
#define OPENMS_FORMAT_HANDLERS_ACQUSHANDLER_H

#include <OpenMS/DATASTRUCTURES/String.h>

#include <map>

using namespace std;

namespace OpenMS
{
	namespace Internal
	{
	  /**
		  @brief Base class for Acqus file handler.
	  */
    class OPENMS_DLLAPI AcqusHandler
    {
      public:
        AcqusHandler(const String& filename);
        ~AcqusHandler();
        double getPosition(const unsigned int index);
        String getParam(const String& param);
        unsigned int getSize();
        
      private:
        AcqusHandler();
        
        map<String, String> params_;
        double DW_;
        unsigned int DELAY_;
        double ML1_;
        double ML2_;
        double ML3_;
        unsigned int TD_;
	  };
	
	} // namespace Internal
} // namespace OpenMS

#endif // OPENMS_FORMAT_HANDLERS_ACQUSHANDLER_H
