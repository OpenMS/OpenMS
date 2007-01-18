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


#ifndef OPENMS_VISUAL_DIALOGS_SPECTRUM2DGOTODIALOG_H
#define OPENMS_VISUAL_DIALOGS_SPECTRUM2DGOTODIALOG_H

#include <OpenMS/VISUAL/DIALOGS/UIC/Spectrum2DGoToDialogTemplate.h>
#include <OpenMS/DATASTRUCTURES/DRange.h>

namespace OpenMS 
{

	/**
		@brief GoTo dialog used to zoom to a m/z and retention time range.
		
		@ingroup Dialogs
	*/
	class Spectrum2DGoToDialog: public Spectrum2DGoToDialogTemplate
	{
		Q_OBJECT
		
		public:
			Spectrum2DGoToDialog( QWidget* parent = 0, const char* name = 0, WFlags fl = 0 );
			~Spectrum2DGoToDialog();
	    void setMinX(float minX);
	    void setMaxX(float maxX);
	    float getMinX();
	    float getMaxX();
	    void setMinY(float minY);
	    void setMaxY(float maxY);
	    float getMinY();
	    float getMaxY();
	  protected:
	    DRange<2> area_;
	    float centerX, centerY;
	  protected slots:
	    virtual void gotoButton_clicked();
	    virtual void setVisibleAreaButton_clicked();
	};

}
#endif // OPENMS_VISUAL_DIALOGS_SPECTRUM2DGOTODIALOG_H

