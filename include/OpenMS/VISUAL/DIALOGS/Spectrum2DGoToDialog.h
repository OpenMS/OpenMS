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
// $Maintainer: $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------


#ifndef OPENMS_VISUAL_DIALOGS_SPECTRUM2DGOTODIALOG_H
#define OPENMS_VISUAL_DIALOGS_SPECTRUM2DGOTODIALOG_H

#include <OpenMS/VISUAL/DIALOGS/UIC/ui_Spectrum2DGoToDialog.h>
#include <OpenMS/CONCEPT/Types.h>

namespace OpenMS 
{

	/**
		@brief GoTo dialog used to zoom to a m/z and retention time range or to a feature.
		
		@ingroup Dialogs
	*/
	class OPENMS_DLLAPI Spectrum2DGoToDialog
		: public QDialog,
			public Ui::Spectrum2DGoToDialogTemplate
	{
		Q_OBJECT
		
		public:
			///Constructor
			Spectrum2DGoToDialog(QWidget* parent=0);
			///Destructor
			~Spectrum2DGoToDialog();
			
			/// Returns if a range should be display (true) or if a feature should be displayed (false)
			bool showRange() const;
			
			///@name Methods for ranges
			//@{
	    ///Sets the data range to display initially
	    void setRange(Real min_rt, Real max_rt, Real min_mz, Real max_mz);
			///Returns the lower RT bound
	    Real getMinRT() const;
			///Returns the upper RT bound
	    Real getMaxRT() const;
			///Returns the lower m/z bound
	    Real getMinMZ() const;
			///Returns the upper m/z bound
	    Real getMaxMZ() const;
	    //@}
	    
			///@name Methods for feature numbers
			//@{
	    ///Returns the selected feature numbers. If a number is retuned, the feature rather than the range should be displayed.
	    UInt getFeatureNumber() const;
	    ///Disables the feature number field
	    void enableFeatureNumber(bool);
	    //@}

	};

}
#endif // OPENMS_VISUAL_DIALOGS_SPECTRUM2DGOTODIALOG_H

