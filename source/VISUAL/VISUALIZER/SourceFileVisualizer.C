// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/VISUALIZER/SourceFileVisualizer.h>

#include <QtGui/QLineEdit>
#include <QtGui/QComboBox>

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
      fillComboBox_(checksum_type_, &temp_.NamesOfChecksumType[temp_.getChecksumType()], 1);
    }
    else
    {
      fillComboBox_(checksum_type_, temp_.NamesOfChecksumType, SourceFile::SIZE_OF_CHECKSUMTYPE);
      checksum_type_->setCurrentIndex(temp_.getChecksumType());
    }
  }

  void SourceFileVisualizer::store()
  {
    ptr_->setNameOfFile(name_of_file_->text());
    ptr_->setPathToFile(path_to_file_->text());
    ptr_->setFileSize(file_size_->text().toFloat());
    ptr_->setFileType(file_type_->text());
    ptr_->setChecksum(checksum_->text(), (SourceFile::ChecksumType)checksum_type_->currentIndex());
    ptr_->setNativeIDType(native_id_type_->text());

    temp_ = (*ptr_);
  }

  void SourceFileVisualizer::undo_()
  {
    update_();
  }

}
