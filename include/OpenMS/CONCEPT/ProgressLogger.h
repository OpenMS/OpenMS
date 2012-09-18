// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
// 
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution 
//    may be used to endorse or promote products derived from this software 
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS. 
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING 
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 
// --------------------------------------------------------------------------
// $Maintainer: Stephan Aiche$
// $Authors: Marc Sturm $
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

