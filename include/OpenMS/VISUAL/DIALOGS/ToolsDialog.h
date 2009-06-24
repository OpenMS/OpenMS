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
// $Maintainer: Marc Sturm $
// $Authors: $
// --------------------------------------------------------------------------


#ifndef OPENMS_VISUAL_DIALOGS_TOOLSDIALOG_H
#define OPENMS_VISUAL_DIALOGS_TOOLSDIALOG_H

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
		@brief TOPP tool selection dialog
		
		In the dialog, the user can 
	  - selelect a TOPP tool
	  - select the options used for the input and output file
		- and set the parameters for the tool
	
		This information can then be used to execute the tool.

		The offered tools depend on the data type set in the constructor.
	
		@ingroup Dialogs
	*/
	class OPENMS_DLLAPI ToolsDialog
		: public QDialog
	{
		Q_OBJECT
		
		public:
			/**
				@brief Constructor
				
				@param parent Qt parent widget
				@param ini_file The file name of the temporary INI file created by this dialog
				@param default_dir The default directory for loading and storing INI files
				@param type The type of data (determines that applicable tools)
			*/
			ToolsDialog(QWidget* parent, String ini_file, String default_dir, LayerData::DataType type);
			///Desctructor
			~ToolsDialog();
			
			/// to get the parameter name for output. Empty if no output was selected.
			String getOutput();
			/// to get the parameter name for input
			String getInput();
			/// to get the currently selected tool-name
			String getTool(); 
	
		private:		
			/// ParamEditor for reading ini-files
			ParamEditor *editor_;
			/// ComboBox for choosing a TOPP-tool
			QComboBox* tools_combo_;
			/// ComboBox for choosing the type of certain tools
			QComboBox* type_combo_;
			/// for choosing an input parameter
			QComboBox* input_combo_;
			/// for choosing an output parameter
			QComboBox* output_combo_;
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

			///Disables the ok button and input/output comboboxes
			void disable_();
			///Enables the ok button and input/output comboboxes
			void enable_();
			
		protected slots:
		
			/// if ok button pressed show the tool output in a new layer, a new window or standard output as messagebox 
			void ok_();
			/// Slot that handles changing of the tool
			void setTool_(int i);
			/// Slot that handles changing of the type and retrieves the defaults
			void setType_(int i);
			/// loads an ini-file into the editor_
			void loadINI_();
			/// stores an ini-file from the editor_
			void storeINI_();
			/// Updates the available types, when the tool changes
			void updateTypes_();
			
	};

}
#endif // OPENMS_VISUAL_DIALOGS_TOOLSDIALOG_H

