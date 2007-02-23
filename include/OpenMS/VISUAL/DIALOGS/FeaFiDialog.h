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


#ifndef OPENMS_VISUAL_DIALOGS_FEAFIDIALOG_H
#define OPENMS_VISUAL_DIALOGS_FEAFIDIALOG_H

#include <OpenMS/VISUAL/DIALOGS/UIC/FeaFiDialogTemplate.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinder.h>

namespace OpenMS 
{

	/**
		@brief Dialog for opening a FeaFi ini file
	
		@ingroup Dialogs
	*/
	class FeaFiDialog
		: public QDialog,
  		public Ui::FeaFiDialogTemplate
	{
		Q_OBJECT
		
		public:
			FeaFiDialog( QWidget* parent = 0 );
			~FeaFiDialog();    
	    FeatureFinder& getFeatureFinder();
	
	  public slots:
	    virtual void loadParamFile();
	    
	  protected:
	    FeatureFinder finder_;
	};

}
#endif // OPENMS_VISUAL_DIALOGS_FEAFIDIALOG_H

