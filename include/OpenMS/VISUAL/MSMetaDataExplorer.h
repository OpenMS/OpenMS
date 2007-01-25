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

//OpenMS
#include <OpenMS/config.h>
#include <OpenMS/METADATA/Sample.h>
#include <OpenMS/METADATA/HPLC.h>
#include <OpenMS/METADATA/Digestion.h>
#include <OpenMS/METADATA/Modification.h>
#include <OpenMS/METADATA/Tagging.h>
#include <OpenMS/METADATA/MetaInfoInterface.h>
#include <OpenMS/METADATA/Gradient.h>
#include <OpenMS/METADATA/Software.h>
#include <OpenMS/METADATA/SourceFile.h>
#include <OpenMS/METADATA/ContactPerson.h>
#include <OpenMS/METADATA/Instrument.h>
#include <OpenMS/METADATA/IonSource.h>
#include <OpenMS/METADATA/IonDetector.h>
#include <OpenMS/METADATA/MassAnalyzer.h>
#include <OpenMS/METADATA/ProcessingMethod.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/ProteinHit.h>
#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/METADATA/ExperimentalSettings.h>
#include <OpenMS/METADATA/Acquisition.h>
#include <OpenMS/METADATA/AcquisitionInfo.h>
#include <OpenMS/METADATA/MetaInfoDescription.h>
#include <OpenMS/METADATA/Precursor.h>
#include <OpenMS/METADATA/InstrumentSettings.h>
#include <OpenMS/METADATA/Identification.h>
#include <OpenMS/METADATA/SpectrumSettings.h>

//STL
#include <iostream>

//QT
#include <qstringlist.h>
#include <qmainwindow.h>
#include <qwidget.h>
#include <qwidgetstack.h>
#include <qsplitter.h>
#include <qpushbutton.h>
#include <qlayout.h>



class QListView;
class QListViewItem;




namespace OpenMS 
{
	/**
		@brief A basic viewer for metadata
		
		It contains a tree viewer showing all objects of the file to be viewed in hierarchical order. <BR>
		The meta info data of the tree items are shown in the right part of the viewer, when they are selected in the tree.
		
		@ingroup Visual
		*/
  class MSMetaDataExplorer 
  	: public QSplitter 
  {
    Q_OBJECT
  
    public: 
		  /// Default constructor
			MSMetaDataExplorer(QWidget *parent =0, const char *name = 0);
		  
		  /**
			@brief A template class to add classes
			
			Simple function to add classes to the viewer. Passes the object to be displayed to the corresponding visualizer class according to its type.
			*/
			template <class T> void add(T *ptr ) 
			{                
				visualize_(*ptr);
			}
			
		/// Defines friend classess that can use the functionality of the subclasses.
		friend class ProteinIdentificationVisualizer;
		friend class IdentificationVisualizer;
      
    private slots:
		 /// Raises the corresponding viewer from the widget stack according to the item selected in the tree.
   	 void showDetails(QListViewItem *item);
		 /// Saves all changes 
		 void saveAll();
				 
    private:
			/// Calls visualizer class for class type ExperimentalSettings.
			void visualize_(ExperimentalSettings &es, QListViewItem* parent=0); 
			/// Calls visualizer class for class type SpectrumSettings.
			void visualize_(SpectrumSettings &es, QListViewItem* parent=0); 
		  /// Calls visualizer class for class type MetaInfoInterface.
			void visualize_(MetaInfoInterface &m, QListViewItem* parent=0); 
			/// Calls visualizer class for class type Sample.
			void visualize_(Sample &s, QListViewItem* parent=0); 
			/// Calls visualizer class for class type HPLC.
			void visualize_(HPLC &h, QListViewItem* parent=0); 
			/// Calls visualizer class for class type Digestion.
			void visualize_(Digestion &d, QListViewItem* parent=0); 
			/// Calls visualizer class for class type Modification.
			void visualize_(Modification &m, QListViewItem* parent=0); 
			/// Calls visualizer class for class type Tagging.
			void visualize_(Tagging &t, QListViewItem* parent=0); 
      /// Calls visualizer class for class type Gradient.
			void visualize_(Gradient &g, QListViewItem* parent=0); 
			/// Calls visualizer class for class type Software.
			void visualize_(Software &s, QListViewItem* parent=0); 
			/// Calls visualizer class for class type SourceFile
			void visualize_(SourceFile &s, QListViewItem* parent=0); 
			/// Calls visualizer class for class type ContactPerson
			void visualize_(ContactPerson &person, QListViewItem* parent=0); 
			/// Calls visualizer class for class type Instrument
			void visualize_(Instrument &instrument, QListViewItem* parent=0); 
			/// Calls visualizer class for class type IonSource
			void visualize_(IonSource &is, QListViewItem* parent=0);
			/// Calls visualizer class for class type IonDetector
			void visualize_(IonDetector &id, QListViewItem* parent=0);
			/// Calls visualizer class for class type MassAnalyzer
			void visualize_(MassAnalyzer &ma, QListViewItem* parent=0);
			/// Calls visualizer class for class type ProcessingMethod
			void visualize_(ProcessingMethod &ma, QListViewItem* parent=0);
			/// Calls visualizer class for class type ProteinIdentification
			void visualize_(ProteinIdentification &pid, QListViewItem* parent=0);
			/// Calls visualizer class for class type ProteinHit
			void visualize_(ProteinHit &phit, QListViewItem* parent=0);
			/// Calls visualizer class for class type PeptideHit
			void visualize_(PeptideHit &pephit, QListViewItem* parent=0, std::vector<ProteinHit> *prots=0);
			/// Calls visualizer class for class type Aquisition
			void visualize_(Acquisition &a, QListViewItem* parent=0);
			/// Calls visualizer class for class type AcquisitionInfo
			void visualize_(AcquisitionInfo &ai, QListViewItem* parent=0);
			/// Calls visualizer class for class type MetaInfoDescription
			void visualize_(MetaInfoDescription &mid,  QListViewItem* parent=0, String *key=0);
			/// Calls visualizer class for class type Precursor
			void visualize_(Precursor &pre, QListViewItem* parent=0);
			/// Calls visualizer class for class type InstrumentSettings
			void visualize_(InstrumentSettings &is, QListViewItem* parent=0);
			/// Calls visualizer class for class type Identification
			void visualize_(Identification &id, QListViewItem* parent=0);
			
			/// Creates id numbers for the widgets on the stack.
			int makeID();
			
			 /** 
			 @brief Update visible protein hits.
			 
			 Updates the protein hit information of objects of type ProteinIdentification.
			 
			 @param tree_item_id 
			 
			 Identifier of the item in the tree belonging to the ProteinIdentification object that is displayed.
			 */
			 void updateProteinHits(ProteinIdentification pid, int tree_item_id);
		 	
			 
			 
			 /** 
			 @brief Update visible peptide hits.
			 
			 Updates the peptide hit information of objects of type Identification depending on a threshold.
			 
			 @param tree_item_id 
			 
			 Identifier of the item in the tree belonging to the ProteinIdentification object that is displayed.
			 */
			 void updatePeptideHits(Identification id, int tree_item_id);
			 
			 
			 /** 
			 @brief Update visible peptide hits.
			 
			 Updates the peptide hit information of objects of type Identification depending on referencing protein hits.
			 
			 @param tree_item_id 
			 
			 Identifier of the item in the tree belonging to the ProteinIdentification object that is displayed.
			 */
			 void updateRefPeptideHits(Identification id, int tree_item_id, String ref_date, String ref_acc);
			 
			 
			 /** 
			 @brief Update visible peptide hits.
			 
			 Updates the peptide hit information of objects of type Identification depending non referencing protein hits.
			 Only peptide hits that do not reference any protein hit will be displayed.			 			 
			 
			 */
			 void updateNonRefPeptideHits(Identification id,  int tree_item_id);
			 
			 
			/// Basic layout for the two main windows.
			QHBoxLayout* vertlayout_;	
			/// Grid layout for basic layout.
			QGridLayout* glayout_;		
			/// Splitter
			QSplitter *splitvert_;
			/// A widgetstack that keeps track of all widgets.
			QWidgetStack *ws_;
			/// The selected widget.
			QWidget *wp_;
			/// Save button
			QPushButton *saveallbutton_;
			/// ID for widgets.
			int obj_id_;
			/// The tree.
			QListView *listview_;
			
			
			
			

  
};

}
#endif
