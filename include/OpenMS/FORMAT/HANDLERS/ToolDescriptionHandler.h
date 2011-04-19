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
// $Maintainer:  Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_HANDLERS_TOOLDESCRIPTIONHANDLER_H
#define OPENMS_FORMAT_HANDLERS_TOOLDESCRIPTIONHANDLER_H

#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/FORMAT/ControlledVocabulary.h>
#include <OpenMS/DATASTRUCTURES/ToolDescription.h>
#include <OpenMS/FORMAT/HANDLERS/ParamXMLHandler.h>

namespace OpenMS
{
	class ProgressLogger;

	namespace Internal
	{

		/**
			@brief XML handler for ToolDescriptionFile

			@note Do not use this class. It is only needed in ToolDescriptionFile.
		*/
		class OPENMS_DLLAPI ToolDescriptionHandler
			: public ParamXMLHandler
		{
		 public:
      /**@name Constructors and destructor */
      //@{

      /// Constructor 
      ToolDescriptionHandler(const String& filename, const String& version);

      /// Destructor
      virtual ~ToolDescriptionHandler();
      //@}


			// Docu in base class
      virtual void endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname);

			// Docu in base class
      virtual void startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes);

			// Docu in base class
      virtual void characters(const XMLCh* const chars, const XMLSize_t length);

			// NOT IMPLEMENTED
			virtual void writeTo(std::ostream& os);

      // Retrieve parsed tool description
      const std::vector<ToolDescription>& getToolDescriptions() const;

      // Set tool description for writing
      void setToolDescriptions(const std::vector<ToolDescription>& td);

		 protected:

      Param p_;

      Internal::ToolExternalDetails tde_;
      Internal::ToolDescription td_;
      std::vector<Internal::ToolDescription> td_vec_;

      String tag_;

      bool in_ini_section_;

		private:

			ToolDescriptionHandler();
			ToolDescriptionHandler(const ToolDescriptionHandler& rhs);
			ToolDescriptionHandler& operator = (const ToolDescriptionHandler& rhs);

		};
	} // namespace Internal
} // namespace OpenMS

#endif // OPENMS_FORMAT_HANDLERS_TOOLDESCRIPTIONHANDLER_H
