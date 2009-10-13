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

#ifndef OPENMS_VISUAL_METADATABROWSER_H
#define OPENMS_VISUAL_METADATABROWSER_H

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/ConsensusMap.h>

//QT
#include <QtGui/QDialog>
#include <QtGui/QTreeWidget>
class QTreeWidgetItem;
class QPushButton;
class QStackedWidget;
class QHBoxLayout;
class QVBoxLayout;
class QGridLayout;

namespace OpenMS
{
	class BaseVisualizerGUI;
	class Acquisition;
	class AcquisitionInfo;
	class ContactPerson;
	class Digestion;
	class ExperimentalSettings;
	class Gradient;
	class HPLC;
	class PeptideIdentification;
	class Instrument;
	class IonDetector;
	class IonSource;
	class MassAnalyzer;
	class MetaInfo;
	class MetaInfoDescription;
	class MetaInfoInterface;
	class MetaInfoRegistry;
	class Modification;
	class PeptideHit;
	class Precursor;
	class DataProcessing;
	class ProteinHit;
	class ProteinIdentification;
	class Sample;
	class SampleTreatment;
	class Software;
	class SourceFile;
	class SpectrumSettings;
	class Tagging;
	class DocumentIdentifier;
	class Product;

	/**
		@brief A meta data visualization widget

		@image html MetaDataBrowser.png

		It contains a tree view showing all objects of the meta data to be viewed in hierarchical order.

		The meta info data of the tree items are shown in the right part of the viewer, when they are selected in the tree.

		If the data has beed modified exec() returns @em true . Otherwise @em false is returned.

		@improvement Add generic mechanism to add items to data vectors e.g. for Instrument - IonSource (Hiwi)

		@ingroup Visual
	*/
  class OPENMS_DLLAPI MetaDataBrowser 
  	: public QDialog
  {
    Q_OBJECT

    public:

			/// Constructor with flag for edit mode
			MetaDataBrowser(bool editable = FALSE, QWidget *parent = 0, bool modal = FALSE );

			/// Adds a peak map
			template <class PeakType>
			void add(MSExperiment<PeakType>& exp)
			{
				add(static_cast<ExperimentalSettings&>(exp));
				treeview_->expandItem( treeview_->findItems(QString::number(0),Qt::MatchExactly , 1).first() );
			}

			/// Adds a peak spectrum
			template <class PeakType>
			void add(MSSpectrum<PeakType>& spectrum)
			{
				//spectrum settings
				add(static_cast<SpectrumSettings&>(spectrum));

				//MetaInfoDescriptions
	      for (Size i=0; i<spectrum.getFloatDataArrays().size();++i)
	      {
	      	add(spectrum.getFloatDataArrays()[i]);
	      }
	      for (Size i=0; i<spectrum.getIntegerDataArrays().size();++i)
	      {
	      	add(spectrum.getIntegerDataArrays()[i]);
	      }
	      for (Size i=0; i<spectrum.getStringDataArrays().size();++i)
	      {
	      	add(spectrum.getStringDataArrays()[i]);
	      }
	      
				add(static_cast<MetaInfoInterface&>(spectrum));
				
				treeview_->expandItem( treeview_->findItems(QString::number(0),Qt::MatchExactly , 1).first() );
			}

			/// Adds a feature map
			template <class FeatureType>
			void add(FeatureMap<FeatureType>& map)
			{
				//identifier
				add(static_cast<DocumentIdentifier&>(map));
				
				//protein ids
				for (Size i=0; i<map.getProteinIdentifications().size(); ++i)
				{
					add(map.getProteinIdentifications()[i]);
				}

				//unassigned peptide ids
				for (Size i=0; i<map.getUnassignedPeptideIdentifications().size(); ++i)
				{
					add(map.getUnassignedPeptideIdentifications()[i]);
				}
				
				treeview_->expandItem( treeview_->findItems(QString::number(0),Qt::MatchExactly , 1).first() );
			}
			/// Adds a feature
			void add(Feature& feature);
			/// Adds a consensus feature
			void add(ConsensusFeature& feature);

			/// Adds a consensus map
			void add(ConsensusMap& map);

			/**
				@brief A generic function to add data.

				The meta data information of all classes that for which a visualize_ method exists can be visualized.
			*/
			template <class MetaDataType>
			void add(MetaDataType& meta_data_object)
			{
				visualize_(meta_data_object);
				treeview_->expandItem( treeview_->findItems(QString::number(0),Qt::MatchExactly , 1).first() );
			}

			///	 Check if mode is editable or not
			bool isEditable();

			/// Defines friend classess that can use the functionality of the subclasses.
			friend class ProteinIdentificationVisualizer;
			friend class PeptideIdentificationVisualizer;

    public slots:

			/// Set a list of error strings due to invalid date format.
			void setStatus(std::string status);

		protected slots:

			/// Raises the corresponding viewer from the widget stack according to the item selected in the tree.
			void showDetails_();

			/// Saves all changes and close explorer
			void saveAll_();

    protected:

    	///@name Visualizer for the different classes
    	//@{
			void visualize_(ExperimentalSettings& meta, QTreeWidgetItem* parent=0);
			void visualize_(SpectrumSettings& meta, QTreeWidgetItem* parent=0);
		 	void visualize_(MetaInfoInterface& meta, QTreeWidgetItem* parent=0);
			void visualize_(Sample& meta, QTreeWidgetItem* parent=0);
			void visualize_(HPLC& meta, QTreeWidgetItem* parent=0);
			void visualize_(Digestion& meta, QTreeWidgetItem* parent=0);
			void visualize_(Modification& meta, QTreeWidgetItem* parent=0);
			void visualize_(Tagging& meta, QTreeWidgetItem* parent=0);
			void visualize_(Gradient& meta, QTreeWidgetItem* parent=0);
			void visualize_(Software& meta, QTreeWidgetItem* parent=0);
			void visualize_(ScanWindow& meta, QTreeWidgetItem* parent=0);
			void visualize_(SourceFile& meta, QTreeWidgetItem* parent=0);
			void visualize_(ContactPerson& meta, QTreeWidgetItem* parent=0);
			void visualize_(Instrument& meta, QTreeWidgetItem* parent=0);
			void visualize_(IonSource& meta, QTreeWidgetItem* parent=0);
			void visualize_(IonDetector& meta, QTreeWidgetItem* parent=0);
			void visualize_(MassAnalyzer& meta, QTreeWidgetItem* parent=0);
			void visualize_(DataProcessing& meta, QTreeWidgetItem* parent=0);
			void visualize_(ProteinIdentification& meta, QTreeWidgetItem* parent=0);
			void visualize_(ProteinHit& meta, QTreeWidgetItem* parent=0);
			void visualize_(PeptideHit& meta, QTreeWidgetItem* parent=0);
			void visualize_(Acquisition& meta, QTreeWidgetItem* parent=0);
			void visualize_(AcquisitionInfo& meta, QTreeWidgetItem* parent=0);
			void visualize_(MetaInfoDescription& meta,  QTreeWidgetItem* parent=0);
			void visualize_(Precursor& meta, QTreeWidgetItem* parent=0);
			void visualize_(Product& meta, QTreeWidgetItem* parent=0);
			void visualize_(InstrumentSettings& meta, QTreeWidgetItem* parent=0);
			void visualize_(PeptideIdentification& meta, QTreeWidgetItem* parent=0);
			void visualize_(DocumentIdentifier& meta, QTreeWidgetItem* parent=0);
			//@}

			/// Visualizes all elements of a container
			template<typename ContainerType>
			void visualizeAll_(ContainerType& container, QTreeWidgetItem* parent)
			{
				for (typename ContainerType::iterator it=container.begin(); it!=container.end(); ++it)
				{
					visualize_(*it, parent);
				}
			}

			/// Connects the Signals of all visualier classes with Slot setStatus()
			void connectVisualizer_(BaseVisualizerGUI* ptr);

			/// Filters hits according to a score @a threshold. Takes the score orientation into account
			void filterHits_(DoubleReal threshold, bool higher_better, int tree_item_id);
			/// Shows hits.
			void showAllHits_(int tree_item_id);

			/// A list of setting errors due to invalid formats.
			std::string status_list_;

			/// Indicates the mode
			bool editable_;

			/// A widgetstack that keeps track of all widgets.
			QStackedWidget* ws_;
			/// Save button
			QPushButton* saveallbutton_;
			/// Close Button
			QPushButton* closebutton_;
			/// Cancel Button
			QPushButton* cancelbutton_;
			/// Undo Button
			QPushButton* undobutton_;

			/// The tree.
			QTreeWidget* treeview_;
	};
}
#endif
