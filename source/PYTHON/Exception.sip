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
// $Id: Exception.sip,v 1.1.1.1 2005/04/04 13:03:04 marc_sturm Exp $
// $Author: marc_sturm $
// --------------------------------------------------------------------------
//

namespace Exception 
{
	
	class Base 
		:	std::exception
	{
%TypeHeaderCode
	#include <OpenMS/CONCEPT/Exception.h>
	using namespace OpenMS;
%End
		public:
		Base() throw();
		Base(const char* /file/, int /line/) throw();
		Base(const char* /file/, int /line/, const std::string& /name/ , const std::string& /message/) throw();
		Base(const Exception::Base& /exception/)	throw();
		virtual ~Base() throw();
		const char* getName() const throw();
		virtual const char* what() const throw();
		int getLine() const	throw();
		const char* getFile() const	throw();
	};		

	class IndexUnderflow 
		:  Exception::Base
	{
%TypeHeaderCode
	#include <OpenMS/CONCEPT/Exception.h>
	using namespace OpenMS;
%End
		public:
		IndexUnderflow(const char* /file/, int /line/, Index /index/ = 0, Size /size/ = 0) throw();
		virtual ~IndexUnderflow() throw();
	};

	class SizeUnderflow 
		: Exception::Base
	{
%TypeHeaderCode
	#include <OpenMS/CONCEPT/Exception.h>
	using namespace OpenMS;
%End
		public:
		SizeUnderflow(const char* /file/, int /line/, Size /size/ = 0) throw();
		virtual ~SizeUnderflow() throw();
	};

	class IndexOverflow 
		: Exception::Base
	{
%TypeHeaderCode
	#include <OpenMS/CONCEPT/Exception.h>
	using namespace OpenMS;
%End
		public:
		IndexOverflow(const char* /file/, int /line/, Index /index/ = 0, Size /size/ = 0) throw();
		virtual ~IndexOverflow() throw();
	};
	class InvalidRange 
		: Exception::Base
	{
%TypeHeaderCode
	#include <OpenMS/CONCEPT/Exception.h>
	using namespace OpenMS;
%End
		public:
		InvalidRange(const char* /file/, int /line/) throw();
		virtual ~InvalidRange() throw();
	};
	class InvalidSize 
		: Exception::Base
	{
%TypeHeaderCode
	#include <OpenMS/CONCEPT/Exception.h>
	using namespace OpenMS;
%End
		public:
		InvalidSize(const char* /file/, int /line/, Size /size/ = 0) throw();
		virtual ~InvalidSize() throw();
	};
	class OutOfRange 
		: Exception::Base
	{
%TypeHeaderCode
	#include <OpenMS/CONCEPT/Exception.h>
	using namespace OpenMS;
%End
		public:
		OutOfRange(const char* /file/, int /line/) throw();
		virtual ~OutOfRange() throw();
	};
	class InvalidFormat 
		: Exception::Base
	{
%TypeHeaderCode
	#include <OpenMS/CONCEPT/Exception.h>
	using namespace OpenMS;
%End
		public:
		InvalidFormat(const char* /file/, int /line/, const std::string& /s/)	throw();
		virtual ~InvalidFormat() throw();
	};
	class IllegalSelfOperation
		: Exception::Base
	{
%TypeHeaderCode
	#include <OpenMS/CONCEPT/Exception.h>
	using namespace OpenMS;
%End
		public:
		IllegalSelfOperation(const char* /file/, int /line/) throw();
		virtual ~IllegalSelfOperation() throw();
	};
	class NullPointer 
		: Exception::Base
	{
%TypeHeaderCode
	#include <OpenMS/CONCEPT/Exception.h>
	using namespace OpenMS;
%End
		public:
		NullPointer(const char* /file/, int /line/)	throw();
		virtual ~NullPointer() throw();
	};
	class InvalidIterator
		: Exception::Base
	{
%TypeHeaderCode
	#include <OpenMS/CONCEPT/Exception.h>
	using namespace OpenMS;
%End
		public:
		InvalidIterator(const char* /file/, int /line/) throw();
		virtual ~InvalidIterator() throw();
	};
	class IncompatibleIterators
		: Exception::Base
	{
%TypeHeaderCode
	#include <OpenMS/CONCEPT/Exception.h>
	using namespace OpenMS;
%End
		public:
		IncompatibleIterators(const char* /file/, int /line/)	throw();
		virtual ~IncompatibleIterators() throw();
	};
	class NotImplemented
		: Exception::Base
	{
%TypeHeaderCode
	#include <OpenMS/CONCEPT/Exception.h>
	using namespace OpenMS;
%End
		public:
		NotImplemented(const char* /file/, int /line/) throw();
		virtual ~NotImplemented() throw();
	};
	class IllegalTreeOperation
		: Exception::Base
	{
%TypeHeaderCode
	#include <OpenMS/CONCEPT/Exception.h>
	using namespace OpenMS;
%End
		public:
		IllegalTreeOperation(const char* /file/, int /line/) throw();
		virtual ~IllegalTreeOperation() throw();
	};
	class OutOfMemory
		: Exception::Base, std::bad_alloc
	{
%TypeHeaderCode
	#include <OpenMS/CONCEPT/Exception.h>
	using namespace OpenMS;
%End
		public:
		OutOfMemory(const char* /file/, int /line/, Size /size/ = 0) throw();			
		virtual ~OutOfMemory() throw();
	};
	class BufferOverflow 
		: Exception::Base
	{
%TypeHeaderCode
	#include <OpenMS/CONCEPT/Exception.h>
	using namespace OpenMS;
%End
		public:
		BufferOverflow(const char* /file/, int /line/) throw();
		virtual ~BufferOverflow() throw();
	};
	class DivisionByZero 
		: Exception::Base
	{
%TypeHeaderCode
	#include <OpenMS/CONCEPT/Exception.h>
	using namespace OpenMS;
%End
		public:
		DivisionByZero(const char* /file/, int /line/) throw();
		virtual ~DivisionByZero() throw();
	};
	class OutOfGrid 
		: Exception::Base
	{
%TypeHeaderCode
	#include <OpenMS/CONCEPT/Exception.h>
	using namespace OpenMS;
%End
		public:
		OutOfGrid(const char* /file/, int /line/)	throw();
		virtual ~OutOfGrid() throw();
	};
	class FileNotFound 
		: Exception::Base
	{
%TypeHeaderCode
	#include <OpenMS/CONCEPT/Exception.h>
	using namespace OpenMS;
%End
		public:
		FileNotFound(const char* /file/, int /line/, const std::string& /filename/)	throw();
		virtual ~FileNotFound() throw();
		std::string getFilename() const	throw();
	};
	class IllegalPosition 
		: Exception::Base
	{
%TypeHeaderCode
	#include <OpenMS/CONCEPT/Exception.h>
	using namespace OpenMS;
%End
		public:
		IllegalPosition(const char* /file/, int /line/, float /x/, float /y/, float /z/) throw();
		virtual ~IllegalPosition() throw();
	};
	class ParseError
		: Exception::Base
	{
%TypeHeaderCode
	#include <OpenMS/CONCEPT/Exception.h>
	using namespace OpenMS;
%End
		public:
		ParseError(const char* /file/, int /line/, const std::string& /expression/,	const std::string& /message/)	throw();
		virtual ~ParseError() throw();
	};
	class UnableToCreateFile
		: Exception::Base
	{
%TypeHeaderCode
	#include <OpenMS/CONCEPT/Exception.h>
	using namespace OpenMS;
%End
		public:
		UnableToCreateFile(const char* /file/, int /line/, const std::string& /filename/)	throw();
		virtual ~UnableToCreateFile() throw();
		std::string getFilename() const	throw();
	};

	class GlobalExceptionHandler
	{
%TypeHeaderCode
	#include <OpenMS/CONCEPT/Exception.h>
	using namespace OpenMS;
%End
		public:
		GlobalExceptionHandler() throw();
		static void setName(const std::string& /name/) throw();
		static void setMessage(const std::string& /message/) throw();
		static void setLine(int /line/)	throw();
		static void setFile(const std::string& /file/) throw();
		static void set(const std::string& /file/, int /line/, 
										 const std::string& /name/, const std::string& /message/)
			throw();
	};

	// GlobalExceptionHandler globalHandler;
	// std::ostream& operator << (std::ostream& /os/, const Exception::Base& /e/);

};
