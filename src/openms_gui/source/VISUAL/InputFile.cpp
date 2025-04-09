// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

// OpenMS includes
#include <OpenMS/VISUAL/InputFile.h>
#include <ui_InputFile.h>

#include <OpenMS/SYSTEM/File.h>

#include <QtWidgets/QMessageBox>
#include <QtWidgets/QFileDialog>
#include <QtWidgets/QCompleter>
#include <QFileSystemModel>
#include <QDragEnterEvent>
#include <QMimeData>

namespace OpenMS
{
  InputFile::InputFile(QWidget* parent)
    : QWidget(parent),
      file_format_filter_(),
      ui_(new Ui::InputFileTemplate)
  {
    ui_->setupUi(this);
    QCompleter* completer = new QCompleter(this);
    completer->setModel(new QFileSystemModel(completer));
    ui_->line_edit->setCompleter(completer);
    connect(ui_->browse_button, SIGNAL(clicked()), this, SLOT(showFileDialog()));
  }

  InputFile::~InputFile()
  {
    delete ui_;
  }

  void InputFile::dragEnterEvent(QDragEnterEvent* e)
  {
    // file dropped from a window manager come as single URL
    if (e->mimeData()->urls().size() == 1)
    {
      e->acceptProposedAction();
    }
  }

  void InputFile::dropEvent(QDropEvent* e)
  {
    QStringList files;
    for (const QUrl& url : e->mimeData()->urls())
    {
      setFilename(url.toLocalFile());
      break;
    }
  }

  void InputFile::dragMoveEvent(QDragMoveEvent* p_event)
  {
    // TODO allow filtering?
    //if (!p_event->mimeData()->hasFormat(MY_MIMETYPE))
    //{
    //  p_event->ignore();
    //  return;
    //}
    p_event->accept();
  }

  void InputFile::setFilename(const QString& filename)
  {
    ui_->line_edit->setText(filename);
    emit updatedFile(filename);
    setCWD(File::path(filename).toQString());
  }

  QString InputFile::getFilename() const
  {
    return ui_->line_edit->text();
  }

  void InputFile::setFileFormatFilter(const QString& fff)
  {
    file_format_filter_ = fff;
  }

  const QString& InputFile::getCWD() const
  {
    return cwd_;
  }

  void InputFile::setCWD(const QString& cwd, bool force)
  {
    if (force || cwd_.isEmpty())
    {
      cwd_ = cwd;
      emit updatedCWD(cwd_);
    }
  }

  void InputFile::showFileDialog()
  {
    QFileInfo fi(getFilename()); // get path from current file as starting directory for selection
    
    QString file_name = QFileDialog::getOpenFileName(this, tr("Specify input file"), cwd_, file_format_filter_);
    if (!file_name.isEmpty())
    {
      setFilename(file_name);
    }
  }


} // namespace
