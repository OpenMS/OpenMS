// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Stephan Aiche$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_CONCEPT_MACROS_H
#define OPENMS_CONCEPT_MACROS_H


#include <OpenMS/config.h>
#include <OpenMS/CONCEPT/Exception.h>

#include <string>
#include <cstring>

/**
	@defgroup Conditions Condition macros

	@brief Macros used for to enforce preconditions and postconditions.

	These macros are enabled if debug info is enabled and optimization is disabled in configure.
	Otherwise they are replaced by an empty string, so they won't cost any performance.

	The macros throw Exception::Precondition or Exception::Postcondition respectively if the condition fails.

	@ingroup Concept

	@{
*/

#ifdef OPENMS_ASSERTIONS

/**
	@brief Precondition macro.

	@hideinitializer
*/
#define OPENMS_PRECONDITION(condition, message)\
	if (!(condition))\
	{\
		Exception::Precondition e(__FILE__, __LINE__, __PRETTY_FUNCTION__, #condition);\
		if (std::strcmp(message,"")!=0)\
		{\
      ::std::string tmp(e.getMessage());\
			tmp += " ";\
			tmp += ::std::string(message);\
			e.setMessage(tmp);\
		}\
		throw e;\
	}\

/**
	@brief Postcondition macro.

	@hideinitializer
*/
#define OPENMS_POSTCONDITION(condition, message)\
	if (!(condition))\
	{\
		Exception::Postcondition e(__FILE__, __LINE__, __PRETTY_FUNCTION__, #condition);\
		if (std::strcmp(message,"")!=0)\
		{\
      std::string tmp(e.getMessage());\
			tmp += " ";\
			tmp += std::string(message);\
			e.setMessage(tmp);\
		}\
		throw e;\
	}\

#else

/**
	@brief Precondition macro.

	@hideinitializer
*/
#define OPENMS_PRECONDITION(condition, message)

/**
	@brief Postcondition macro.

	@hideinitializer
*/
#define OPENMS_POSTCONDITION(condition, message)

#endif

/** @} */ // end of Conditions

#endif //OPENMS_CONCEPT_MACROS_H
