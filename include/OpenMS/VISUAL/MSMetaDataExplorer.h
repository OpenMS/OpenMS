// -*- mode: C++; tab-width: 2; -*-
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
// $Maintainer: Marc Sturm  $
// --------------------------------------------------------------------------

#ifndef OPENMS_VISUAL_MSMETADATAEXPLORER_H
#define OPENMS_VISUAL_MSMETADATAEXPLORER_H

#include <OpenMS/DATASTRUCTURES/String.h>

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
	class BaseVisualizer;
	class Acquisition;
	class AcquisitionInfo;
	class ContactPerson;
	class Digestion;
	class ExperimentalSettings;
	class Gradient;
	class HPLC;
	class PeptideIdentification;
	class Instrument;
	class InstrumentSettings;
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
	class ProcessingMethod;
	class ProteinHit;
	class ProteinIdentification;
	class Sample;
	class SampleTreatment;
	class Software;
	class SourceFile;
	class SpectrumSettings;
	class Tagging;
	
	/**
		@brief A meta data visualization widget
		
		@image html MSMetaDataExplorer.png
		
		It contains a tree view showing all objects of the file to be viewed in hierarchical order. <BR>
		The meta info data of the tree items are shown in the right part of the viewer, when they are selected in the tree.
		
		@ingroup Visual
	*/
  class MSMetaDataExplorer 
  	: public QDialog
  {
    Q_OBJECT
  
    public: 
			/// Constructor with flag for edit mode
			MSMetaDataExplorer(bool editable = FALSE, QWidget *parent = 0, bool modal = FALSE );
		  
			
			
			/**
			@brief A template function to add classes
				
			The meta data information of many different objects can be visualized using this function. 
			The object is passed to one of the type-specific visualize_ methods, managing the visulization of the meta data.
			*/
			template <class T> void visualize(T& class_reference) 
			{                
				visualize_(class_reference);
				treeview_->expandItem(  treeview_->findItems(QString::number(0),Qt::MatchExactly , 1).first() );
			}
			
			///	 Check if mode is editable or not
			bool isEditable();	
				
			/// Defines friend classess that can use the functionality of the subclasses.
			friend class ProteinIdentificationVisualizer;
			friend class PeptideIdentificationVisualizer;
      
    public slots:
			/// Set a list of error strings due to invalid date format.
			void setStatus(std::string status);

		private slots:
			/// Raises the corresponding viewer from the widget stack according to the item selected in the tree.
			void showDetails_(QTreeWidgetItem *item, int column);
			
			/// Saves all changes and close explorer
			void saveAll_();
							 
    private:
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
			void visualize_(SourceFile& meta, QTreeWidgetItem* parent=0); 
			void visualize_(ContactPerson& meta, QTreeWidgetItem* parent=0); 
			void visualize_(Instrument& meta, QTreeWidgetItem* parent=0); 
			void visualize_(IonSource& meta, QTreeWidgetItem* parent=0);
			void visualize_(IonDetector& meta, QTreeWidgetItem* parent=0);
			void visualize_(MassAnalyzer& meta, QTreeWidgetItem* parent=0);
			void visualize_(ProcessingMethod& meta, QTreeWidgetItem* parent=0);
			void visualize_(ProteinIdentification& meta, QTreeWidgetItem* parent=0);
			void visualize_(ProteinHit& meta, QTreeWidgetItem* parent=0);
			void visualize_(PeptideHit& meta, QTreeWidgetItem* parent=0);
			void visualize_(Acquisition& meta, QTreeWidgetItem* parent=0);
			void visualize_(AcquisitionInfo& meta, QTreeWidgetItem* parent=0);
			void visualize_(MetaInfoDescription& meta,  QTreeWidgetItem* parent=0);
			void visualize_(Precursor& meta, QTreeWidgetItem* parent=0);
			void visualize_(InstrumentSettings& meta, QTreeWidgetItem* parent=0);
			void visualize_(PeptideIdentification& meta, QTreeWidgetItem* parent=0);
			//@}
			
			/// Connects the Signals of all visualier classes with Slot setStatus()
			void connectVisualizer_(BaseVisualizer*);

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
