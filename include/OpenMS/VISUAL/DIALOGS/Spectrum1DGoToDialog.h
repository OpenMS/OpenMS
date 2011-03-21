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
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------


#ifndef OPENMS_VISUAL_DIALOGS_SPECTRUM1DGOTODIALOG_H
#define OPENMS_VISUAL_DIALOGS_SPECTRUM1DGOTODIALOG_H

#include <QtGui/QDialog>
#include <OpenMS/VISUAL/DIALOGS/UIC/ui_Spectrum1DGoToDialog.h>
#include <OpenMS/CONCEPT/Types.h>

namespace OpenMS 
{

	/**
		@brief simple goto/set visible area dialog for exact placement of the viewing window
		
		@ingroup Dialogs
	*/
	class OPENMS_GUI_DLLAPI Spectrum1DGoToDialog
		: public QDialog,
			public Ui::Spectrum1DGoToDialogTemplate
	{
		Q_OBJECT
		
		public:
			///Constructor
			Spectrum1DGoToDialog(QWidget* parent = 0);
			///Destructor
			~Spectrum1DGoToDialog();
			
			///Sets the m/z range displayed initially
	    void setRange(Real min, Real max);
	    ///Returns the lower m/z bound
	    Real getMin() const;
	    ///Returns the upper m/z bound
	    Real getMax() const;
	};

}
#endif // OPENMS_VISUAL_DIALOGS_SPECTRUM1DGOTODIALOG_H

