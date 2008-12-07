// -*- Mode: C++; tab-width: 2; -*-
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
// $Maintainer: Marc Sturm$
// --------------------------------------------------------------------------


#ifndef OPENMS_VISUAL_DIALOGS_SPECTRUMALIGNMENTDIALOG_H
#define OPENMS_VISUAL_DIALOGS_SPECTRUMALIGNMENTDIALOG_H

#include <OpenMS/VISUAL/DIALOGS/UIC/SpectrumAlignmentDialogTemplate.h>

namespace OpenMS 
{
	/**
		@brief Dialog which allows the user to enter a tolerance and performs a spectrum alignment.
		
		Dialog which allows the user to enter a tolerance and performs a spectrum alignment. The
		alignment can only be performed if the active window is a Spectrum1DWidget with two
		canvasses (mirror mode).
		
		@ingroup Dialogs
	*/
	class SpectrumAlignmentDialog
		: public QDialog,
  		public Ui::SpectrumAlignmentDialogTemplate
	{
		Q_OBJECT
				
		public:
			
			/// Constructor
			SpectrumAlignmentDialog();
			
		protected slots:
		
		protected:
			
		private:
			
	};
	
}
#endif // OPENMS_VISUAL_DIALOGS_SPECTRUMALIGNMENTDIALOG_H
