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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------


#ifndef OPENMS_VISUAL_DIALOGS_SPECTRUM1DGOTODIALOG_H
#define OPENMS_VISUAL_DIALOGS_SPECTRUM1DGOTODIALOG_H

#include <QtGui/QDialog>
#include <OpenMS/VISUAL/DIALOGS/UIC/Spectrum1DGoToDialogTemplate.h>

namespace OpenMS 
{

	/**
		@brief simple goto/set visible area dialog for exact placement of the viewing window
		
		@ingroup Dialogs
	*/
	class Spectrum1DGoToDialog
		: public QDialog,
			public Ui::Spectrum1DGoToDialogTemplate
			
	{
		Q_OBJECT
		
		public:
			Spectrum1DGoToDialog( QWidget* parent = 0 );
			~Spectrum1DGoToDialog();    
	    void setMin(double value);
	    void setMax(double value);
	    float getMin();
	    float getMax();
	};

}
#endif // OPENMS_VISUAL_DIALOGS_SPECTRUM1DGOTODIALOG_H

