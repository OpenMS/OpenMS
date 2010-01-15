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


#ifndef OPENMS_VISUAL_DIALOGS_FEATUREEDITDIALOG_H
#define OPENMS_VISUAL_DIALOGS_FEATUREEDITDIALOG_H

#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/VISUAL/DIALOGS/UIC/ui_FeatureEditDialog.h>

namespace OpenMS 
{
	/**
		@brief Dialog for editing a feature
		
		@ingroup Dialogs
	*/
	class OPENMS_DLLAPI FeatureEditDialog
		: public QDialog,
  		public Ui::FeatureEditDialogTemplate
	{
		Q_OBJECT
				
		public:
			
			/// Constructor
			FeatureEditDialog(QWidget* parent);
			/// Sets the feature
			void setFeature(const Feature& feature);
			/// Returns the feature
			const Feature& getFeature() const;
			
		protected:

			/// The feature to edit
			mutable Feature feature_;
			
		private:
			///Not implemented
			FeatureEditDialog();
			
	};
	
}
#endif // OPENMS_VISUAL_DIALOGS_FEATUREEDITDIALOG_H
