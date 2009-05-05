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
// $Maintainer: Johannes Junker $
// $Authors: Johannes Junker $
// --------------------------------------------------------------------------


#ifndef OPENMS_VISUAL_DIALOGS_TOPPASTOOLCONFIGDIALOG_H
#define OPENMS_VISUAL_DIALOGS_TOPPASTOOLCONFIGDIALOG_H

class QComboBox;
class QPushButton;
class QRadioButton;
class QString; 

#include <QtGui/QDialog>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/VISUAL/LayerData.h>

namespace OpenMS 
{
	class ParamEditor;
		
	/**
		@brief TOPP tool configuration dialog
		
		In the dialog, the user can set the parameters for the tool
	
		This information can then be used to execute the tool.
	
		@ingroup Dialogs
	*/
	class OPENMS_DLLAPI TOPPASToolConfigDialog
		: public QDialog
	{
		Q_OBJECT
		
		public:
			/**
				@brief Constructor
				
				@param parent Qt parent widget
				@param ini_file The file name of the temporary INI file created by this dialog
				@param default_dir The default directory for loading and storing INI files
				@param tool_name The name of the tool
				@param tool_type The type of the tool
				@param instance_nr The (unique) instance number
			*/
			TOPPASToolConfigDialog(QWidget* parent, String ini_file, String default_dir, String tool_name, String tool_type, UInt instance_nr);
			///Desctructor
			~TOPPASToolConfigDialog();
	
		private:		
			/// ParamEditor for reading ini-files
			ParamEditor *editor_;
			/// Param for loading the ini-file
			Param arg_param_;
			/// Param for loading configuration information in the ParamEditor
			Param vis_param_;
			/// ok-button connected with slot ok_()
			QPushButton* ok_button_;
			/// map for getting the parameter name from the full path in arg_param
			std::map<String,String> arg_map_;
			/// Location of the temporary INI file this dialog works on 
			String ini_file_;
			/// default-dir of ini-file to open
			String default_dir_;
			/// name of ini-file
			QString filename_;
			/// The name of the tool
			String tool_name_;
			/// The type of the tool
			String tool_type_;
			/// The instance number of this specific node
			UInt instance_nr_;
			
		protected slots:
			/// Slot for OK button
			void ok_();
			/// loads an ini-file into the editor_
			void loadINI_();
			/// stores an ini-file from the editor_
			void storeINI_();
			
	};

}
#endif // OPENMS_VISUAL_DIALOGS_TOOLSDIALOG_H

