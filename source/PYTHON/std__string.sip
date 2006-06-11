// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// $Id: std__string.sip,v 1.1.1.1 2005/04/04 13:03:04 marc_sturm Exp $
//

namespace std
{
class string
{
%TypeHeaderCode
	#include <string>
%End
  public:

  SIP_PYOBJECT __str__();
%MethodCode
	if (sipCpp == 0)
	{
		sipRes = PyString_FromString("");
	}
	else
	{
		sipRes = PyString_FromString(sipCpp->c_str());
	}
%End

SIP_PYOBJECT __repr__();
%MethodCode
	if (sipCpp == 0)
	{
		sipRes = PyString_FromString("");
	}
	else
	{
		sipRes = PyString_FromString(sipCpp->c_str());
	}
%End

%ConvertToTypeCode
  if (sipIsErr == NULL)
    return (PyString_Check(sipPy) || sipIsSubClassInstance(sipPy,sipClass_std_string));

  if (sipPy == Py_None)
  {
    *sipCppPtr = new std::string;

    return 1;
  }
  if (PyString_Check(sipPy))
  {
    *sipCppPtr = new std::string(PyString_AS_STRING(sipPy));

    return 1;
  }

  *sipCppPtr = (std::string *)sipConvertToCpp(sipPy,sipClass_std_string,sipIsErr);

  return 0;
%End
};

}; // namespace std
