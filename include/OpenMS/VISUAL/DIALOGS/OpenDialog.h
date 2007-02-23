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
// $Maintainer: Marc Sturm$
// --------------------------------------------------------------------------


#ifndef OPENMS_VISUAL_DIALOGS_OPENDIALOG_H
#define OPENMS_VISUAL_DIALOGS_OPENDIALOG_H

#include <OpenMS/VISUAL/DIALOGS/UIC/OpenDialogTemplate.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/FORMAT/FileHandler.h>

#include <vector>

namespace OpenMS 
{
	class Param;
	
	/**
		@brief Open spectrum dialog
		
		Allows selecting a spectrum/map from a file or the DB.
		
		@todo Implement meta data preview with DB load metadata only (Marc)
		
		@ingroup Dialogs
	*/
	class OpenDialog
		: public QDialog,
  		public Ui::OpenDialogTemplate
	{
		Q_OBJECT
		
		public:
			/// Data source
			enum DataSource 
			{
				FILE, 	///< getNames() contains file names
				DB      ///< getNames() contains database ids
			};

			/// Available mowers
			enum Mower
			{
				NO_MOWER, 	     ///< none
				NOISE_ESTIMATOR  ///< cut of according to noise estimator
			};
			
			/// constructor
			OpenDialog(Param& preferences, QWidget* parent = 0 );
			/// destructor
			virtual ~OpenDialog();    
			
			/// returns the list of file names / database ids
			const std::vector<String>& getNames() const;
			/// returns the source of data
			DataSource getSource() const;
			/// returns true, if 2D mode is to be used for maps
			bool isViewMaps2D() const;
			/// returns the mower to be used
			Mower getMower() const;
			/// returns the location where the spectra are to be opened
			bool isOpenAsNewTab() const;
			/// returns the forced file type (FileHandler::UNKNOWN if none)
			FileHandler::Type forcedFileType() const;
			
		protected slots:
			/// browses filesystem or database for spectra to open
			virtual void browse();
			/// shows the metadata of the experiment, if exactly one spectrum is selected
			virtual void showMetadata();
		
		protected:
			///list of file names or database ids
			std::vector<String> names_;
			///preferences
			Param& prefs_;
			
			///returns the preference entry @p name (for convenience only)
			const DataValue& getPref_(const String& name) const;
			///returns the preference entry @p name as SignedInt (for convenience only)
			SignedInt getPrefAsInt_(const String& name) const;

	};

}
#endif // OPENMS_VISUAL_DIALOGS_OPENDIALOG_H

