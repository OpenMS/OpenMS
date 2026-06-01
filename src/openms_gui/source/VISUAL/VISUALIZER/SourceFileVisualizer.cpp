// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/VISUALIZER/SourceFileVisualizer.h>
#include <OpenMS/VISUAL/MISC/Qt5Port.h>

#include <QtWidgets/QLineEdit>
#include <QtWidgets/QComboBox>

#include <iostream>

using namespace std;

namespace OpenMS
{

  SourceFileVisualizer::SourceFileVisualizer(bool editable, QWidget * parent) :
    BaseVisualizerGUI(editable, parent),
    BaseVisualizer<SourceFile>()
  {
    addLabel_("Modify source file information");
    addSeparator_();
    addLineEdit_(name_of_file_, "Name of file");
    addLineEdit_(path_to_file_, "Path to file");
    addLineEdit_(file_size_, "File size (in MB)");
    addLineEdit_(file_type_, "File type");
    addLineEdit_(checksum_, "Checksum");
    addComboBox_(checksum_type_, "Checksum type");
    addLineEdit_(native_id_type_, "Native ID type of spectra");

    finishAdding_();
  }

  void SourceFileVisualizer::update_()
  {
    name_of_file_->setText(temp_.getNameOfFile().c_str());
    path_to_file_->setText(temp_.getPathToFile().c_str());
    file_size_->setText(String(temp_.getFileSize()).c_str());
    file_type_->setText(temp_.getFileType().c_str());
    checksum_->setText(temp_.getChecksum().c_str());
    native_id_type_->setText(temp_.getNativeIDType().c_str());

    if (!isEditable())
    {
      fillComboBox_(checksum_type_, &temp_.NamesOfChecksumType[static_cast<size_t>(temp_.getChecksumType())], 1);
    }
    else
    {
      fillComboBox_(checksum_type_, temp_.NamesOfChecksumType, static_cast<int>(SourceFile::ChecksumType::SIZE_OF_CHECKSUMTYPE));
      checksum_type_->setCurrentIndex(static_cast<int>(temp_.getChecksumType()));
    }
  }

  void SourceFileVisualizer::store()
  {
    ptr_->setNameOfFile(fromQString(name_of_file_->text()));
    ptr_->setPathToFile(fromQString(path_to_file_->text()));
    ptr_->setFileSize(file_size_->text().toFloat());
    ptr_->setFileType(fromQString(file_type_->text()));
    ptr_->setChecksum(fromQString(checksum_->text()), (SourceFile::ChecksumType)checksum_type_->currentIndex());
    ptr_->setNativeIDType(fromQString(native_id_type_->text()));

    temp_ = (*ptr_);
  }

  void SourceFileVisualizer::undo_()
  {
    update_();
  }

}
