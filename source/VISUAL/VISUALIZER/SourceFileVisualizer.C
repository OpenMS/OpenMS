// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copymain (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  this library is distributed in the hope that it will be useful,
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

#include <OpenMS/VISUAL/VISUALIZER/SourceFileVisualizer.h>

#include <QtGui/QLineEdit>
#include <QtGui/QComboBox>

#include <iostream>

using namespace std;

namespace OpenMS
{
	
	SourceFileVisualizer::SourceFileVisualizer(bool editable, QWidget* parent)
		: BaseVisualizerGUI(editable, parent),
			BaseVisualizer<SourceFile>()
	{
		addLabel_("Modify source file information");
		addSeparator_();	
		addLineEdit_(name_of_file_, "Name of file" );
		addLineEdit_(path_to_file_, "Path to file" );
		addLineEdit_(file_size_, "File size (in MB)" );
		addLineEdit_(file_type_, "File type" );
		addLineEdit_(checksum_, "Checksum" );
		addComboBox_(checksum_type_, "Checksum type" );
		addLineEdit_(native_id_type_, "Native ID type of spectra");
		
		finishAdding_();
	}
	
	void SourceFileVisualizer::update_()
	{
	  name_of_file_->setText(temp_.getNameOfFile().c_str());
		path_to_file_->setText(temp_.getPathToFile().c_str() );
		file_size_->setText(String(temp_.getFileSize()).c_str());
	  file_type_->setText(temp_.getFileType().c_str());
		checksum_->setText(temp_.getChecksum().c_str());
		native_id_type_->setText(temp_.getNativeIDType().c_str());

		if(! isEditable())
		{
			fillComboBox_(checksum_type_,& temp_.NamesOfChecksumType[temp_.getChecksumType()] , 1);
		}
		else
		{
			fillComboBox_(checksum_type_, temp_.NamesOfChecksumType , SourceFile::SIZE_OF_CHECKSUMTYPE);
			checksum_type_->setCurrentIndex(temp_.getChecksumType());
		}
	}
	
	void SourceFileVisualizer::store()
	{
		ptr_->setNameOfFile(name_of_file_->text());
		ptr_->setPathToFile(path_to_file_->text());
		ptr_->setFileSize(file_size_->text().toFloat());
		ptr_->setFileType(file_type_->text());
		ptr_->setChecksum(checksum_->text(),(SourceFile::ChecksumType)checksum_type_->currentIndex());
		ptr_->setNativeIDType(native_id_type_->text());
		
		temp_=(*ptr_);
	}
	
	void SourceFileVisualizer::undo_()
	{
		update_();
	}

}
