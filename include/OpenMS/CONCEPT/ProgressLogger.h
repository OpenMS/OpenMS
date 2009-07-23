// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_CONCEPT_PROGRESSLOGGER_H
#define OPENMS_CONCEPT_PROGRESSLOGGER_H

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/SYSTEM/StopWatch.h>

class QProgressDialog;

namespace OpenMS
{
	class String;

	/**
	@brief Base class for all classes that want to report their progess.

	Per default the progress log is disabled. Use setLogType to enable it.

	Use startProgress, setProgress and endProgress for the actual logging.

	@note All methods are const, so it can be used through a const reference or in const methods as well!
	*/
	class OPENMS_DLLAPI ProgressLogger
	{
	 public:
		/// Constructor
		ProgressLogger();

		/// Destructor
		~ProgressLogger();

		///Possible log types
		enum LogType
			{
				CMD,    ///< Command line progress
				GUI,      ///< Progress dialog
				NONE  ///< No progress logging
			};

		/// Sets the progress log that should be used. The default type is NONE!
		void setLogType(LogType type) const;

		/// Returns the type of progress log being used.
		LogType getLogType() const;

		/**
		@brief Initializes the progress display

		Sets the progress range from @p begin to @p end.
		If @p begin equals @p end, setProgress only indicates that
		the program is still running, but without showing any absolute progress value.

		Sets the label to @p label.

		@note Make sure to call setLogType first!
		*/
		void startProgress(SignedSize begin, SignedSize end, const String& label) const;

		/// Sets the current progress
		void setProgress(SignedSize value) const;

		/// Ends the progress display
		void endProgress() const;

	 protected:
		mutable LogType type_;
		mutable SignedSize begin_;
		mutable SignedSize end_;
		mutable SignedSize value_;
		mutable QProgressDialog* dlg_;
		mutable StopWatch stop_watch_;
		mutable time_t last_invoke_;
		static int recursion_depth_;
	};

} // namespace OpenMS

#endif //OPENMS_CONCEPT_PROGRESSLOGGER_H

