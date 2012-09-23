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
// $Maintainer: Johannes Junker $
// $Authors: Johannes Junker $
// --------------------------------------------------------------------------


#ifndef OPENMS_VISUAL_DIALOGS_TOPPASOUTPUTFILESDIALOG_H
#define OPENMS_VISUAL_DIALOGS_TOPPASOUTPUTFILESDIALOG_H

#include <OpenMS/config.h>
#include <OpenMS/VISUAL/DIALOGS/UIC/ui_TOPPASOutputFilesDialog.h>

namespace OpenMS
{
  /**
      @brief Dialog which allows to specify the directory for the output files

      @ingroup TOPPAS_elements
      @ingroup Dialogs
  */
  class OPENMS_GUI_DLLAPI TOPPASOutputFilesDialog :
    public QDialog,
    public Ui::TOPPASOutputFilesDialogTemplate
  {
    Q_OBJECT

public:

    /// Constructor
    TOPPASOutputFilesDialog(const QString & dir_name, int num_jobs);

    /// Returns the name of the directory
    QString getDirectory();

    /// Returns the maximum number of jobs in the spinbox
    int getNumJobs();

    /// Returns if the directory is valid (is a directory and writable)
    static bool dirNameValid(const QString & dir_name);

public slots:

    /// Lets the user select the directory via a file dialog
    void showFileDialog();

protected slots:

    /// Called when OK is pressed; checks if the selected file is valid
    void checkValidity_();

  };

}
#endif // OPENMS_VISUAL_DIALOGS_TOPPASOUTPUTFILESDIALOG_H
