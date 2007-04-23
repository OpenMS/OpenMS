// -*- mode: C++; tab-width: 2; -*-
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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_CONCEPT_PROGRESS_H
#define OPENMS_CONCEPT_PROGRESS_H

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/String.h>

class QProgressDialog;

namespace OpenMS
{	
	/**
		@brief Base class for all classes that want to report their progess. 
	*/
	class ProgressLogger
	{
		public:
			///Constrructor
			ProgressLogger();
			
			//Destructor
			~ProgressLogger();
			
			///Possible log types
			enum LogType
			{
				CMD,    ///< Command line progress
				GUI,      ///< Progress dialog
				NONE  ///< No progress logging
			};
			
			///Sets the progress log that should be used. The default type is NONE!
			void setLogType(LogType type) const;
			
			/**	
				@brief Initializes the progress display
				
				Sets the progress range from @p begin to @p end.
				If @p begin equals @p end, setProgress only indicates that
				the program is still running, but without showing any absolute progress value.
				
				Sets the label to @p label.

				@note Make sure to call setLogType first!
			*/
			void initProgress(UInt begin, UInt end, const String& label) const;
			
			/// Sets the current progress
			void setProgress(UInt value) const;
		
		protected:
			mutable LogType type_;
			mutable UInt begin_;
			mutable UInt end_;
			mutable UInt value_;
			mutable QProgressDialog* dlg_;
	};

}//namespace OpenMS

#endif //OPENMS_CONCEPT_PROGRESS_H

