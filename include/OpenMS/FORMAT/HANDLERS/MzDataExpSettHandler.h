// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_HANDLERS_MZDATAEXPSETTHANDLER_H
#define OPENMS_FORMAT_HANDLERS_MZDATAEXPSETTHANDLER_H

#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>

#include <xercesc/sax2/Attributes.hpp>

namespace OpenMS
{
	class ExperimentalSettings;
	class ContactPerson;
	class MassAnalyzer;
	
	namespace Internal
	{

  /**
  	@brief XML handler for experimental settings of MzDataFile

		MapType has to be a MSExperiment or have the same interface.
  	Do not use this class. It is only needed in MzDataFile.
  */
  class MzDataExpSettHandler
		: public XMLHandler
  {
    public:
      /**@name Constructors and destructor */
      //@{
      /// Constructor for a write-only handler
      MzDataExpSettHandler(ExperimentalSettings& exp, const String& filename);

      /// Constructor for a read-only handler
      MzDataExpSettHandler(const ExperimentalSettings& exp, const String& filename);

      /// Destructor
      virtual ~MzDataExpSettHandler();
      //@}

			// Docu in base class
      virtual void endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname);
			
			// Docu in base class
      virtual void startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes);
			
			// Docu in base class
      virtual void characters(const XMLCh* const chars, unsigned int length);

  		///Writes the contents to a stream
			void writeTo(std::ostream& os);

    protected:
		/// map pointer for reading
		ExperimentalSettings* exp_;
		/// map pointer for writing
		const ExperimentalSettings* cexp_;

		/** @brief read attributes of MzData's cvParamType

			Example:
			&lt;cvParam cvLabel="psi" accession="PSI:1000001" name="@p name" value="@p value"/&gt;
			@p name and sometimes @p value are defined in the MzData ontology.
		*/
		void cvParam_(const XMLCh* name, const XMLCh* value);
  };

	} // namespace Internal

} // namespace OpenMS

#endif //OPENMS_FORMAT_HANDLERS_MZDATAEXPSETTHANDLER_H
