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
// $Id: Types.sip,v 1.1.1.1 2005/04/04 13:03:04 marc_sturm Exp $
// $Author: marc_sturm $
// --------------------------------------------------------------------------
//

//
// We map the standard OpenMS typedefs (which we require for portability)
// back to the default Python types and vice versa. Looks ugly, but is
// quite simple.
//

%MappedType Distance
{
%TypeHeaderCode
	#include <OpenMS/CONCEPT/Types.h>
	#include "sipOpenMSDistance.h"
	using namespace OpenMS;
%End

%ConvertToTypeCode
	if (sipIsErr == NULL)
	{
		return PyInt_Check(sipPy);
	}
	if (PyInt_Check(sipPy))
	{
		*sipCppPtr = new Distance;
		**sipCppPtr = PyInt_AsLong(sipPy);
	
		return 1;
	}

	return 0;
%End
%ConvertFromTypeCode
	if (sipCpp == NULL)
	{
		return PyInt_FromLong(0);
	}
	return PyInt_FromLong(*sipCpp);
%End
};


%MappedType Handle
{
%TypeHeaderCode
	#include <OpenMS/CONCEPT/Types.h>
	#include "sipOpenMSHandle.h"
	using namespace OpenMS;
%End

%ConvertToTypeCode
	if (sipIsErr == NULL)
	{
		return PyInt_Check(sipPy);
	}
	if (PyInt_Check(sipPy))
	{
		*sipCppPtr = new Handle;
		**sipCppPtr = PyInt_AsLong(sipPy);
	
		return 1;
	}

	return 0;
%End
%ConvertFromTypeCode
	if (sipCpp == NULL)
	{
		return PyInt_FromLong(0);
	}
	return PyInt_FromLong(*sipCpp);
%End
};

%MappedType SignedInt
{
%TypeHeaderCode
	#include <OpenMS/CONCEPT/Types.h>
	#include "sipOpenMSSignedInt.h"
	using namespace OpenMS;
%End

%ConvertToTypeCode
	if (sipIsErr == NULL)
	{
		return PyInt_Check(sipPy);
	}
	if (PyInt_Check(sipPy))
	{
		*sipCppPtr = new SignedInt;
		**sipCppPtr = PyInt_AsLong(sipPy);
	
		return 1;
	}

	return 0;
%End
%ConvertFromTypeCode
	if (sipCpp == NULL)
	{
		return PyInt_FromLong(0);
	}
	return PyInt_FromLong(*sipCpp);
%End
};

%MappedType UnsignedInt
{
%TypeHeaderCode
	#include <OpenMS/CONCEPT/Types.h>
	#include "sipOpenMSUnsignedInt.h"
	using namespace OpenMS;
%End

%ConvertToTypeCode
	if (sipIsErr == NULL)
	{
		return PyInt_Check(sipPy);
	}
	if (PyInt_Check(sipPy))
	{
		*sipCppPtr = new UnsignedInt;
		**sipCppPtr = PyInt_AsLong(sipPy);
	
		return 1;
	}

	return 0;
%End
%ConvertFromTypeCode
	if (sipCpp == NULL)
	{
		return PyInt_FromLong(0);
	}
	return PyInt_FromLong(*sipCpp);
%End
};

%MappedType Index
{
%TypeHeaderCode
	#include <OpenMS/CONCEPT/Types.h>
	#include "sipOpenMSIndex.h"
	using namespace OpenMS;
%End

%ConvertToTypeCode
	if (sipIsErr == NULL)
	{
		return PyInt_Check(sipPy);
	}
	if (PyInt_Check(sipPy))
	{
		*sipCppPtr = new Index;
		**sipCppPtr = PyInt_AsLong(sipPy);
	
		return 1;
	}

	return 0;
%End
%ConvertFromTypeCode
	if (sipCpp == NULL)
	{
		return PyInt_FromLong(0);
	}
	return PyInt_FromLong(*sipCpp);
%End
};

%MappedType Size
{
%TypeHeaderCode
	#include <OpenMS/CONCEPT/Types.h>
	#include "sipOpenMSSize.h"
	using namespace OpenMS;
%End

%ConvertToTypeCode
	if (sipIsErr == NULL)
	{
		return PyInt_Check(sipPy);
	}
	if (PyInt_Check(sipPy))
	{
		*sipCppPtr = new Size;
		**sipCppPtr = PyInt_AsLong(sipPy);
	
		return 1;
	}

	return 0;
%End
%ConvertFromTypeCode
	if (sipCpp == NULL)
	{
		return PyInt_FromLong(0);
	}
	return PyInt_FromLong(*sipCpp);
%End
};

%MappedType Time
{
%TypeHeaderCode
	#include <OpenMS/CONCEPT/Types.h>
	#include "sipOpenMSTime.h"
	using namespace OpenMS;
%End

%ConvertToTypeCode
	if (sipIsErr == NULL)
	{
		return PyInt_Check(sipPy);
	}
	if (PyInt_Check(sipPy))
	{
		*sipCppPtr = new Time;
		**sipCppPtr = PyInt_AsLong(sipPy);
	
		return 1;
	}

	return 0;
%End
%ConvertFromTypeCode
	if (sipCpp == NULL)
	{
		return PyInt_FromLong(0);
	}
	return PyInt_FromLong(*sipCpp);
%End
};

%MappedType HashIndex
{
%TypeHeaderCode
	#include <OpenMS/CONCEPT/Types.h>
	#include "sipOpenMSHashIndex.h"
	using namespace OpenMS;
%End

%ConvertToTypeCode
	if (sipIsErr == NULL)
	{
		return PyInt_Check(sipPy);
	}
	if (PyInt_Check(sipPy))
	{
		*sipCppPtr = new HashIndex;
		**sipCppPtr = PyInt_AsLong(sipPy);
	
		return 1;
	}

	return 0;
%End
%ConvertFromTypeCode
	if (sipCpp == NULL)
	{
		return PyInt_FromLong(0);
	}
	return PyInt_FromLong(*sipCpp);
%End
};

%MappedType Position
{
%TypeHeaderCode
	#include <OpenMS/CONCEPT/Types.h>
	#include "sipOpenMSDistance.h"
	using namespace OpenMS;
%End

%ConvertToTypeCode
	if (sipIsErr == NULL)
	{
		return PyInt_Check(sipPy);
	}
	if (PyInt_Check(sipPy))
	{
		*sipCppPtr = new Position;
		**sipCppPtr = PyInt_AsLong(sipPy);
	
		return 1;
	}

	return 0;
%End
%ConvertFromTypeCode
	if (sipCpp == NULL)
	{
		return PyInt_FromLong(0);
	}
	return PyInt_FromLong(*sipCpp);
%End
};

%MappedType Real
{
%TypeHeaderCode
	#include <OpenMS/CONCEPT/Types.h>
	#include "sipOpenMSReal.h"
	using namespace OpenMS;
%End

%ConvertToTypeCode
	if (sipIsErr == NULL)
	{
		return PyFloat_Check(sipPy);
	}
	if (PyInt_Check(sipPy))
	{
		*sipCppPtr = new Real;
		**sipCppPtr = PyFloat_AsDouble(sipPy);
	
		return 1;
	}

	return 0;
%End
%ConvertFromTypeCode
	if (sipCpp == NULL)
	{
		return PyFloat_FromDouble(0.0);
	}
	return PyFloat_FromDouble(*sipCpp);
%End
};

%MappedType DoubleReal
{
%TypeHeaderCode
	#include <OpenMS/CONCEPT/Types.h>
	#include "sipOpenMSDoubleReal.h"
	using namespace OpenMS;
%End

%ConvertToTypeCode
	if (sipIsErr == NULL)
	{
		return PyFloat_Check(sipPy);
	}
	if (PyInt_Check(sipPy))
	{
		*sipCppPtr = new DoubleReal;
		**sipCppPtr = PyFloat_AsDouble(sipPy);
	
		return 1;
	}

	return 0;
%End
%ConvertFromTypeCode
	if (sipCpp == NULL)
	{
		return PyFloat_FromDouble(0.0);
	}
	return PyFloat_FromDouble(*sipCpp);
%End
};

%MappedType Property
{
%TypeHeaderCode
	#include <OpenMS/CONCEPT/Types.h>
	#include "sipOpenMSProperty.h"
	using namespace OpenMS;
%End

%ConvertToTypeCode
	if (sipIsErr == NULL)
	{
		return PyInt_Check(sipPy);
	}
	if (PyInt_Check(sipPy))
	{
		*sipCppPtr = new Property;
		**sipCppPtr = PyInt_AsLong(sipPy);
	
		return 1;
	}

	return 0;
%End
%ConvertFromTypeCode
	if (sipCpp == NULL)
	{
		return PyInt_FromLong(0);
	}
	return PyInt_FromLong(*sipCpp);
%End
};

%MappedType ErrorCode
{
%TypeHeaderCode
	#include <OpenMS/CONCEPT/Types.h>
	#include "sipOpenMSErrorCode.h"
	using namespace OpenMS;
%End

%ConvertToTypeCode
	if (sipIsErr == NULL)
	{
		return PyInt_Check(sipPy);
	}
	if (PyInt_Check(sipPy))
	{
		*sipCppPtr = new ErrorCode;
		**sipCppPtr = PyInt_AsLong(sipPy);
	
		return 1;
	}

	return 0;
%End
%ConvertFromTypeCode
	if (sipCpp == NULL)
	{
		return PyInt_FromLong(0);
	}
	return PyInt_FromLong(*sipCpp);
%End
};

%MappedType Byte
{
%TypeHeaderCode
	#include <OpenMS/CONCEPT/Types.h>
	#include "sipOpenMSByte.h"
	using namespace OpenMS;
%End

%ConvertToTypeCode
	if (sipIsErr == NULL)
	{
		return PyInt_Check(sipPy);
	}
	if (PyInt_Check(sipPy))
	{
		*sipCppPtr = new Byte;
		long x = PyInt_AsLong(sipPy);
		if ((x > 255) && (x < 0))
		{
			x = 0;
		}
		**sipCppPtr = (Byte)x;
	
		return 1;
	}

	return 0;
%End
%ConvertFromTypeCode
	if (sipCpp == NULL)
	{
		return PyInt_FromLong(0);
	}
	return PyInt_FromLong(*sipCpp);
%End
};

%MappedType PointerSizeUInt
{
%TypeHeaderCode
	#include <OpenMS/CONCEPT/Types.h>
	#include "sipOpenMSPointerSizeUInt.h"
	using namespace OpenMS;
%End

%ConvertToTypeCode
	if (sipIsErr == NULL)
	{
		return PyInt_Check(sipPy);
	}
	if (PyInt_Check(sipPy))
	{
		*sipCppPtr = new PointerSizeUInt;
		**sipCppPtr = PyInt_AsLong(sipPy);
	
		return 1;
	}

	return 0;
%End
%ConvertFromTypeCode
	if (sipCpp == NULL)
	{
		return PyInt_FromLong(0);
	}
	return PyInt_FromLong(*sipCpp);
%End
};

%MappedType PointerSizeInt
{
%TypeHeaderCode
	#include <OpenMS/CONCEPT/Types.h>
	#include "sipOpenMSPointerSizeInt.h"
	using namespace OpenMS;
%End

%ConvertToTypeCode
	if (sipIsErr == NULL)
	{
		return PyInt_Check(sipPy);
	}
	if (PyInt_Check(sipPy))
	{
		*sipCppPtr = new PointerSizeInt;
		**sipCppPtr = PyInt_AsLong(sipPy);
	
		return 1;
	}

	return 0;
%End
%ConvertFromTypeCode
	if (sipCpp == NULL)
	{
		return PyInt_FromLong(0);
	}
	return PyInt_FromLong(*sipCpp);
%End
};
