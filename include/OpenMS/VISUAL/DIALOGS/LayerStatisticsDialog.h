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
// $Maintainer: Marc Sturm$
// $Authors: $
// --------------------------------------------------------------------------


#ifndef OPENMS_VISUAL_DIALOGS_LAYERSTATISTICSDIALOG_H
#define OPENMS_VISUAL_DIALOGS_LAYERSTATISTICSDIALOG_H

#include <OpenMS/VISUAL/DIALOGS/UIC/ui_LayerStatisticsDialog.h>
#include <OpenMS/VISUAL/LayerData.h>

#include <QtGui/QPushButton>

#include <utility>
#include <map>

namespace OpenMS 
{
	class SpectrumWidget;
  class SpectrumCanvas;

	/**
		@brief Dialog showing statistics about the data of the current layer
		
		@ingroup Dialogs
	*/
	class OPENMS_DLLAPI LayerStatisticsDialog
		: public QDialog,
  		public Ui::LayerStatisticsDialogTemplate
	{
		Q_OBJECT
				
		public:
			
			/// Constructor
			LayerStatisticsDialog(SpectrumWidget* parent);
			
		protected slots:
			
			/// Shows the distribution accoring to the clicked button
			void showDistribution_();
		
		protected:
			
			/**
				@brief Struct representing the statistics about one meta information
			*/
			struct MetaStatsValue_
			{
				MetaStatsValue_(unsigned long c = 0, DoubleReal mi = 0, DoubleReal ma = 0, DoubleReal a = 0)
				{
					count = c;
					min = mi;
					max = ma;
					avg = a;
				}
				
				unsigned long count;
				DoubleReal min, max, avg;
			};
			
			/// Iterates over RTs of an experiment
			typedef LayerData::ExperimentType::ConstIterator RTIterator_;
			/// Iterates over peaks of a spectrum
			typedef LayerData::ExperimentType::SpectrumType::ConstIterator PeakIterator_;
			/// Iterates over features of a feature map
			typedef LayerData::FeatureMapType::ConstIterator FeatureIterator_;
			/// Iterates over features of a feature map
			typedef LayerData::ConsensusMapType::ConstIterator ConsensusIterator_;
			/// Iterates over the meta_stats map
			typedef std::map<UInt, MetaStatsValue_>::iterator MetaIterator_;
			
			/// Computes the statistics of a peak layer
			void computePeakStats_();
			/// Computes the statistics of a feature layer
			void computeFeatureStats_();
			/// Computes the statistics of a consensus feature layer
			void computeConsensusStats_();
			/// Computes the statistics of all meta data contained in the FloatDataArray of a @p spectrum
			void computeFloatDataArrayStats_(RTIterator_ spectrum);
			/// Brings the meta values of one @p meta_interface (a peak or feature) into the statistics
			void bringInMetaStats_(const MetaInfoInterface& meta_interface);
			/// Computes the averages of all meta values stored in meta_stats and meta_array_stats
			void computeMetaAverages_();
			
			/// Map containing the statistics about all meta information of the peaks/features in the layer
			std::map<UInt,MetaStatsValue_> meta_stats_;
			/// Map containing the statistics about the FloatDataArrays of all spectra in this layer
			std::map<String, MetaStatsValue_> meta_array_stats_;
			/// The canvas of the layer
			SpectrumCanvas* canvas_;
			/// The LayerData object we compute statistics about
			LayerData layer_data_;
			/// Minimum intensity value
			DoubleReal min_intensity_;
			/// Maximum intensity value
			DoubleReal max_intensity_;
			/// Average intensity value
			DoubleReal avg_intensity_;
			/// Minimum charge value
			DoubleReal min_charge_;
			/// Maximum charge value
			DoubleReal max_charge_;
			/// Average charge value
			DoubleReal avg_charge_;
			/// Minimum quality value
			DoubleReal min_quality_;
			/// Maximum quality value
			DoubleReal max_quality_;
			/// Average quality value
			DoubleReal avg_quality_;
			/// Minimum number of elements (for consensus features only)
			DoubleReal min_elements_;
			/// Maximum number of elements (for consensus features only)
			DoubleReal max_elements_;
			/// Average number of elements (for consensus features only)
			DoubleReal avg_elements_;
			
		private:
			///Not implemented
			LayerStatisticsDialog();
			
	};
	
}
#endif // OPENMS_VISUAL_DIALOGS_LAYERSTATISTICSDIALOG_H
