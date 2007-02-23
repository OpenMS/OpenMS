// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2005 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: stefan_heess  $
// --------------------------------------------------------------------------

#ifndef OPENMS_VISUAL_MSMETADATAEXPLORER_H
#define OPENMS_VISUAL_MSMETADATAEXPLORER_H

#include <OpenMS/DATASTRUCTURES/String.h>

//QT
#include <QtGui/QDialog>
class QTreeWidgetItem;
class QPushButton;
class QStackedWidget;
class QTreeWidget;
class QHBoxLayout;
class QVBoxLayout;
class QGridLayout;

namespace OpenMS 
{
	class Acquisition;
	class AcquisitionInfo;
	class ContactPerson;
	class Digestion;
	class ExperimentalSettings;
	class Gradient;
	class HPLC;
	class Identification;
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
		@brief A basic viewer for metadata
		
		It contains a tree viewer showing all objects of the file to be viewed in hierarchical order. <BR>
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
				@brief A template class to add classes
				
				Simple function to add classes to the viewer. Passes the object to be displayed to the corresponding visualizer class according to its type.
			*/
			template <class T> void add(T *ptr ) 
			{                
				visualize_(*ptr);
			}
			
			///	 Check if mode is editable or not
			bool isEditable();	
				
			/// Defines friend classess that can use the functionality of the subclasses.
			friend class ProteinIdentificationVisualizer;
			friend class IdentificationVisualizer;
      
    private slots:
			/// Raises the corresponding viewer from the widget stack according to the item selected in the tree.
			void showDetails_(QTreeWidgetItem *item, int column);
			
			/// Saves all changes and close explorer
			void saveAll_();
				 
    private:
			/// Calls visualizer class for class type ExperimentalSettings.
			void visualize_(ExperimentalSettings& meta, QTreeWidgetItem* parent=0); 
			/// Calls visualizer class for class type SpectrumSettings.
			void visualize_(SpectrumSettings& meta, QTreeWidgetItem* parent=0); 
		  /// Calls visualizer class for class type MetaInfoInterface.
			void visualize_(MetaInfoInterface& meta, QTreeWidgetItem* parent=0); 
			/// Calls visualizer class for class type Sample.
			void visualize_(Sample& meta, QTreeWidgetItem* parent=0); 
			/// Calls visualizer class for class type HPLC.
			void visualize_(HPLC& meta, QTreeWidgetItem* parent=0); 
			/// Calls visualizer class for class type Digestion.
			void visualize_(Digestion& meta, QTreeWidgetItem* parent=0); 
			/// Calls visualizer class for class type Modification.
			void visualize_(Modification& meta, QTreeWidgetItem* parent=0); 
			/// Calls visualizer class for class type Tagging.
			void visualize_(Tagging& meta, QTreeWidgetItem* parent=0); 
      /// Calls visualizer class for class type Gradient.
			void visualize_(Gradient& meta, QTreeWidgetItem* parent=0); 
			/// Calls visualizer class for class type Software.
			void visualize_(Software& meta, QTreeWidgetItem* parent=0); 
			/// Calls visualizer class for class type SourceFile
			void visualize_(SourceFile& meta, QTreeWidgetItem* parent=0); 
			/// Calls visualizer class for class type ContactPerson
			void visualize_(ContactPerson& meta, QTreeWidgetItem* parent=0); 
			/// Calls visualizer class for class type Instrument
			void visualize_(Instrument& meta, QTreeWidgetItem* parent=0); 
			/// Calls visualizer class for class type IonSource
			void visualize_(IonSource& meta, QTreeWidgetItem* parent=0);
			/// Calls visualizer class for class type IonDetector
			void visualize_(IonDetector& meta, QTreeWidgetItem* parent=0);
			/// Calls visualizer class for class type MassAnalyzer
			void visualize_(MassAnalyzer& meta, QTreeWidgetItem* parent=0);
			/// Calls visualizer class for class type ProcessingMethod
			void visualize_(ProcessingMethod& meta, QTreeWidgetItem* parent=0);
			/// Calls visualizer class for class type ProteinIdentification
			void visualize_(ProteinIdentification& meta, QTreeWidgetItem* parent=0);
			/// Calls visualizer class for class type ProteinHit
			void visualize_(ProteinHit& meta, QTreeWidgetItem* parent=0);
			/// Calls visualizer class for class type PeptideHit
			void visualize_(PeptideHit& meta, QTreeWidgetItem* parent=0, std::vector<ProteinHit> *prots=0);
			/// Calls visualizer class for class type Aquisition
			void visualize_(Acquisition& meta, QTreeWidgetItem* parent=0);
			/// Calls visualizer class for class type AcquisitionInfo
			void visualize_(AcquisitionInfo& meta, QTreeWidgetItem* parent=0);
			/// Calls visualizer class for class type MetaInfoDescription
			void visualize_(MetaInfoDescription& meta,  QTreeWidgetItem* parent=0, const String& key="");
			/// Calls visualizer class for class type Precursor
			void visualize_(Precursor& meta, QTreeWidgetItem* parent=0);
			/// Calls visualizer class for class type InstrumentSettings
			void visualize_(InstrumentSettings& meta, QTreeWidgetItem* parent=0);
			/// Calls visualizer class for class type Identification
			void visualize_(Identification& meta, QTreeWidgetItem* parent=0);
			
			/** 
				@brief Update visible protein hits.
				
				Updates the protein hit information of objects of type ProteinIdentification.
				
				@param pid  Identifier of the identification object belonging to this database search.
				@param tree_item_id  Identifier of the item in the tree belonging to the ProteinIdentification object that is displayed.
			*/
			void updateProteinHits_(ProteinIdentification pid, int tree_item_id);
		 	
			 
			 
			/** 
				@brief Update visible peptide hits.
				
				Updates the peptide hit information of objects of type Identification depending on a threshold.
				
				@param id  Identifier of the identification object belonging to this database search.
				@param tree_item_id  Identifier of the item in the tree belonging to the ProteinIdentification object that is displayed.
			*/
			void updatePeptideHits_(Identification id, int tree_item_id);
			
			
			/** 
				@brief Update visible peptide hits.
				
				Updates the peptide hit information of objects of type Identification depending on referencing protein hits.
				
				@param id  Identifier of the identification object belonging to this database search.
				@param tree_item_id Identifier of the item in the tree belonging to the ProteinIdentification object that is displayed.
				@param ref_date Petide Hits referencing to protein hits with the same date of database search.
				@param ref_acc Petide Hits referencing to protein hits with the same acquisition number.
			*/
			void updateRefPeptideHits_(Identification id, int tree_item_id, String ref_date, String ref_acc);
			
			
			/** 
				@brief Update visible peptide hits.
				
				Updates the peptide hit information of objects of type Identification depending non referencing protein hits.<br>
				Only peptide hits that do not reference any protein hit will be displayed.			 			 
				
				Updates the peptide hit information of objects of type Identification depending on referencing protein hits.
				
				@param id  Identifier of the identification object belonging to this database search.
				@param tree_item_id  Identifier of the item in the tree belonging to the ProteinIdentification object that is displayed.
			*/
			void updateNonRefPeptideHits_(Identification id,  int tree_item_id);
			 
			/// Indicates the mode
			bool editable_;
			
			/// Basic layout.
			QGridLayout* basiclayout_;	
			
			/// Remove
			QVBoxLayout* vertlayout_;
			/// Remove
			QHBoxLayout* buttonlayout_;	
			
			//QSplitter 
			QStackedWidget* splitvert_;
			/// A widgetstack that keeps track of all widgets.
			QStackedWidget* ws_;
			/// Save button
			QPushButton* saveallbutton_;
			/// Close Button
			QPushButton* closebutton_;
			/// Cancel Button
			QPushButton* cancelbutton_;
			
			/// The tree.
			QTreeWidget* treeview_;
	};
}
#endif
