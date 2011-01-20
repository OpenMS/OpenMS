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
// $Maintainer: $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------


#ifndef OPENMS_VISUAL_DIALOGS_HISTOGRAMDIALOG_H
#define OPENMS_VISUAL_DIALOGS_HISTOGRAMDIALOG_H

#include <QtGui/QDialog>

#include <OpenMS/MATH/STATISTICS/Histogram.h>
#include <OpenMS/VISUAL/HistogramWidget.h>

namespace OpenMS
{
	/**
		@brief Dialog that show a HistogramWidget.
	
		@ingroup Dialogs
	*/
	class OPENMS_DLLAPI HistogramDialog
		: public QDialog
	{
		Q_OBJECT
		
		public:
			/// Constructor
			HistogramDialog(const Math::Histogram<>& distribution, QWidget* parent=0);
			/// Destructor
			~HistogramDialog();
			
			/// Returns the value of the left splitter
			Real getLeftSplitter();
			/// Returns the value of the right splitter
			Real getRightSplitter();
			
			/// Sets the value of the left splitter
			void setLeftSplitter(Real position);
			/// Sets the value of the right splitter
			void setRightSplitter(Real position);
			
			/// Sets the axis legend
			void setLegend(const String& legend);
			/// Sets log mode
			void setLogMode(bool log_mode);

		protected:
			HistogramWidget *mw_;
	};

} //namespace

#endif //OPENMS_VISUAL_DIALOGS_HISTOGRAMDIALOG_H

