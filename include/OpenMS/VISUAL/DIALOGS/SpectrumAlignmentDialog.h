// -*- Mode: C++; tab-width: 2; -*-
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
// $Maintainer: $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------


#ifndef OPENMS_VISUAL_DIALOGS_SPECTRUMALIGNMENTDIALOG_H
#define OPENMS_VISUAL_DIALOGS_SPECTRUMALIGNMENTDIALOG_H

#include <OpenMS/VISUAL/DIALOGS/UIC/ui_SpectrumAlignmentDialog.h>
#include <OpenMS/CONCEPT/Types.h>

namespace OpenMS 
{
	class Spectrum1DWidget;

	/**
		@brief Lets the user select two spectra and set the parameters for the spectrum alignment.
		
		@ingroup Dialogs
	*/
	class SpectrumAlignmentDialog
		: public QDialog,
  		public Ui::SpectrumAlignmentDialogTemplate
	{
		Q_OBJECT
				
		public:
			
			/// Constructor
			SpectrumAlignmentDialog(Spectrum1DWidget* parent);
			
			/// Returns the index of the selected non-flipped layer
			Int get1stLayerIndex();
			/// Returns the index of the selected flipped layer
			Int get2ndLayerIndex();
			
		protected slots:
		
		protected:
			
			/// Stores the layer indices of the layers in the left list (non-flipped layers)
			std::vector<UInt> layer_indices_1_;
			/// Stores the layer indices of the layers in the right list (flipped layers)
			std::vector<UInt> layer_indices_2_;
			
		private:
			
	};
	
}
#endif // OPENMS_VISUAL_DIALOGS_SPECTRUMALIGNMENTDIALOG_H
