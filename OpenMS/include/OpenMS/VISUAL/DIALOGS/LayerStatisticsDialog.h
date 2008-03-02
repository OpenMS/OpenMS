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


#ifndef OPENMS_VISUAL_DIALOGS_LAYERSTATISTICSDIALOG_H
#define OPENMS_VISUAL_DIALOGS_LAYERSTATISTICSDIALOG_H

#include <OpenMS/VISUAL/DIALOGS/UIC/LayerStatisticsDialogTemplate.h>
#include <OpenMS/VISUAL/SpectrumWidget.h>

#include <QtGui/QPushButton>

#include <utility>
#include <map>

namespace OpenMS 
{
	/**
		@brief Dialog showing statistics about the data of the current layer
		
	*/
	class LayerStatisticsDialog
		: public QDialog,
  		public Ui::LayerStatisticsDialogTemplate
	{
		Q_OBJECT
				
		public:
			
			/// constructor
			LayerStatisticsDialog(SpectrumWidget* parent);
			
			/**
				@brief struct representing the statistics about one meta information
					
			*/
			struct MetaStatsValue
			{
				MetaStatsValue(int c = 0, int mi = 0, int ma = 0, int a = 0)
				{
					count = c;
					min = mi;
					max = ma;
					avg = a;
				}
				
				unsigned long count;
				DoubleReal min, max, avg;
			};
			
			typedef MSExperiment<>::Iterator RTIterator;
			typedef MSSpectrum<>::Iterator PeakIterator;
			typedef FeatureMap<>::Iterator FeatureIterator;
			typedef std::map<String, MetaStatsValue*>::iterator MetaIterator;
			
		protected slots:
		
		protected:
						
		private:
			///Not implemented
			LayerStatisticsDialog();
			
	};
	
}
#endif // OPENMS_VISUAL_DIALOGS_LAYERSTATISTICSDIALOG_H
