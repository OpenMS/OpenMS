// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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
// --------------------------------------------------------------------------


#ifndef OPENMS_VISUAL_DIALOGS_TOOLSDIALOG_H
#define OPENMS_VISUAL_DIALOGS_TOOLSDIALOG_H

class QComboBox;
class QPushButton;
class QRadioButton;
class QString; 

#include <QtGui/QDialog>
#include <OpenMS/DATASTRUCTURES/Param.h>

namespace OpenMS 
{
	class ParamEditor;
	class Param;
	class LayerData;
		
	/**
		@brief Dialog for executing a TOPP Tool
		
		In this dialog the TOPP Tools are executed by syscall
		Remember to set the PATH-variable to the OpenMS/bin directory before executing the TOPP-tools dialog!
		Before clicking ok-button you should open a file with spectrum data.
		
		@ingroup TOPPView_elements
	*/
	class ToolsDialog
		: public QDialog
	{
		Q_OBJECT
		
		public:
			/// constructor is given tmp_dir where the input-,output-files are saved
			ToolsDialog( QWidget* parent, String tmp_dir, String default_dir, const LayerData* layer);
			/// to get the parameter name for output
			String getOutput();
			/// to get the parameter name for input
			String getInput();
			/// to get the currently selected tool-name
			String getTool();
			
			/**
				@name bool functions for checking ouput action
			*/
			//@{
			bool openAsWindow();
			bool openAsLayer();
			bool noOutputAction();
			//@}
			~ToolsDialog();    
	
		private:		
			/// ParamEditor for reading ini-files
			ParamEditor *editor_;
			/// ComboBox for choosing a TOPP-tool
			QComboBox* tools_combo_;
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
			/// choosing a window as visualization of the tool-output
			QRadioButton* window_radio_;
			/// choosing a layer as visualization of the tool-output
			QRadioButton* layer_radio_;
			/// option for choosing only the output of the tool, which means it is not loaded via addSpectrum(...)
			QRadioButton* output_radio_;
			/// map for getting the parameter name from the full path in arg_param
			std::map<String,String> arg_map_;
			/// parameter chosen for input
			String input_string_;
			/// parameter chosen for output
			String output_string_;
			/// Temporary files directory
			String tmp_dir_;
			/// default-dir of ini-file to open
			String default_dir_;
			/// name of ini-file
			QString filename_;
			
		protected slots:
			/// if ok button pressed show the tool output in a new layer, a new window or standard output as messagebox 
			void ok_();
			/// get tool name from combobox
			void setTool_(int i);
			/// loads an ini-file into the editor_
			bool loadIni();
			/// stores an ini-file from the editor_
			bool storeIni();
			
	};

}
#endif // OPENMS_VISUAL_DIALOGS_TOOLSDIALOG_H_H

