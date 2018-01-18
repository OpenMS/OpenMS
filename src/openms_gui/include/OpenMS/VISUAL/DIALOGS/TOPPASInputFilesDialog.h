// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Maintainer: Johannes Veit $
// $Authors: Johannes Junker $
// --------------------------------------------------------------------------

#ifndef OPENMS_VISUAL_DIALOGS_TOPPASINPUTFILESDIALOG_H
#define OPENMS_VISUAL_DIALOGS_TOPPASINPUTFILESDIALOG_H

// OpenMS_GUI config
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

#include <OpenMS/VISUAL/DIALOGS/UIC/ui_TOPPASInputFilesDialog.h>

namespace OpenMS
{
  /**
      @brief Dialog which allows to specify a list of input files

      @ingroup TOPPAS_elements
      @ingroup Dialogs
  */
  class OPENMS_GUI_DLLAPI TOPPASInputFilesDialog :
    public QDialog,
    public Ui::TOPPASInputFilesDialogTemplate
  {
    Q_OBJECT

public:
    /// support drag'n'drop of files from OS window manager
    void dragEnterEvent(QDragEnterEvent *e) override;
    /// support drag'n'drop of files from OS window manager
    void dropEvent(QDropEvent *e) override;

    /// Constructor
    TOPPASInputFilesDialog(const QStringList& list, const QString& cwd);

    /// Stores the list of all filenames in the list widget in @p files
    void getFilenames(QStringList& files) const;

    /// get the CWD (according to most recently added file)
    const QString& getCWD() const;

    /// support Ctrl+C to copy currently selected items to clipboard
    void keyPressEvent(QKeyEvent *e) override;

public slots:

    /// Lets the user select files via a file dialog
    void showFileDialog();
    /// Removes all currently selected files from the list
    void removeSelected();
    /// Removes all files from the list
    void removeAll();
    /// Shows a TOPPASInputFileDialog which edits the current item
    void editCurrentItem();
    /// Moves the current item up/downwards
    void moveCurrentItem();

protected:
    /// current working dir, i.e. the last position a file was added from
    QString cwd_;

  };

}
#endif // OPENMS_VISUAL_DIALOGS_TOPPASINPUTFILESDIALOG_H
