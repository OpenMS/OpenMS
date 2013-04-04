// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
// $Maintainer: Timo Sachsenberg$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------


#ifndef OPENMS_VISUAL_DIALOGS_DBOPENDIALOG_H
#define OPENMS_VISUAL_DIALOGS_DBOPENDIALOG_H

#include <vector>
#include <QtGui/QDialog>
#include <OpenMS/CONCEPT/Types.h>

class QLineEdit;
class QTableWidget;

namespace OpenMS
{
  class DBConnection;

  /**
      @brief Dialog that allow selecting a spectrum from a DB.

      @ingroup Dialogs
  */
  class OPENMS_GUI_DLLAPI DBOpenDialog :
    public QDialog
  {
    Q_OBJECT
public:
    /**
        @brief Constructor

        The spectrum ids to load are inserted into the @p result vector.

        An external DB connection is used by handing over @p connection.
    */
    DBOpenDialog(DBConnection & connection, std::vector<UInt> & result, QWidget * parent = 0);
    /// Destructor
    ~DBOpenDialog();

private slots:
    /// Slot for accepting the selection
    void ok();
    /// Slot for refreshing the shown spectra
    void loadSpectra();

protected:
    /// DB connection
    DBConnection & connection_;
    /// reference to the result vector
    std::vector<UInt> & result_;
    /// pointer to the search string lineedit
    QLineEdit * search_string_;
    /// pointer to the table for displaying the overview
    QTableWidget * table_;
  };

} //namespace

#endif //OPENMS_VISUAL_DIALOGS_DBOPENDIALOG_H
