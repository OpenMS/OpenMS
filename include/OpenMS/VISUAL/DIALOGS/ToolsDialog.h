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
// $Maintainer: Stefan Rink $
// --------------------------------------------------------------------------


#ifndef OPENMS_VISUAL_DIALOGS_TOOLSDIALOG_H
#define OPENMS_VISUAL_DIALOGS_TOOLSDIALOG_H

class QString;
class QStringList;
class QComboBox;
class QHBoxLayout;
class QPushButton;
class QRadioButton;
class QLabel;
class QLineEdit;

#include <QtGui/QDialog>
#include <QtCore/QDir>
#include <OpenMS/FORMAT/Param.h>

namespace OpenMS 
{
	class ParamEditor;
	class Param;
	/**
		@brief Dialog for executing a TOPP Tool
		
		TOPP Tools are executed by syscall
		Remember to set the PATH-variable to the OpenMS/bin directory before executing the TOPP-tools dialog!
		Before clicking ok-button you should open a file with spectrum data.
		
		@todo Keep only the grid layout, except for buttons layout (Stefan)
		@todo Add writing of feature pairs (Stefan)
		@todo Write docu (Stefan)
		@todo Replace QLineEdit by QComboBox (Stefan)
		@todo Show only entries below '1' (Stefan)
		
		@ingroup Dialogs
	*/
	class ToolsDialog
		: public QDialog
	{
		Q_OBJECT
		
		public:
			ToolsDialog( QWidget* parent, String tmp_dir );
			String getOutput();
			String getInput();
			String getTool();
			bool isWindow();
			~ToolsDialog();    
	
		private:
			
			/// save input argument in input_
			void setInput_();
			/// save output argument in output_
			void setOutput_();
		
			ParamEditor *editor_;
			QComboBox* tools_combo_;
			QLineEdit* input_edit_;
			QLineEdit* output_edit_;
			Param arg_param_;
			QPushButton* ok_button_;
			QRadioButton* window_radio_;
			QRadioButton* layer_radio_;
			std::map<String,String> arg_map_;
			String input_string_;
			String output_string_;
			/// Temporary files directory
			String tmp_dir_;
			
		protected slots:
			/// ok button pressed
			void ok_();
			/// cancel button pressed
			void cancel_();
			/// get tool name from combobox
			void setTool_(int i);
			
	};

}
#endif // OPENMS_VISUAL_DIALOGS_TOOLSDIALOG_H_H

