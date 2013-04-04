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
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/DIALOGS/DBOpenDialog.h>
#include <OpenMS/FORMAT/DB/DBConnection.h>

#include <QtGui/QPushButton>
#include <QtGui/QLayout>
#include <QtGui/QLineEdit>
#include <QtGui/QTableWidget>
#include <QtGui/QLabel>
#include <QtSql/QSqlQuery>
#include <QtGui/QHBoxLayout>
#include <QtGui/QVBoxLayout>

#include <sstream>

using namespace std;


namespace OpenMS
{
  DBOpenDialog::DBOpenDialog(DBConnection & connection, vector<UInt> & result, QWidget * parent) :
    QDialog(parent),
    connection_(connection),
    result_(result)
  {
    setWindowTitle("Select spectra from DB to open");

    //OK+Cancel button + layout
    QPushButton * cancel_button_ = new QPushButton("&Cancel", this);
    connect(cancel_button_, SIGNAL(clicked()), this, SLOT(reject()));
    QPushButton * ok_button_ = new QPushButton("&OK", this);
    connect(ok_button_, SIGNAL(clicked()), this, SLOT(ok()));

    QHBoxLayout * hbox1_ = new QHBoxLayout();
    hbox1_->addStretch(1);
    hbox1_->addWidget(ok_button_);
    hbox1_->addWidget(cancel_button_);

    //string line edit, button + layout
    QLabel * label_ = new QLabel("Description contains:", this);
    search_string_ = new QLineEdit(this);
    QPushButton * search_button_ = new QPushButton("refresh", this);
    connect(search_button_, SIGNAL(clicked()), this, SLOT(loadSpectra()));

    QHBoxLayout * hbox2_ = new QHBoxLayout();
    hbox2_->addWidget(label_);
    hbox2_->addWidget(search_string_);
    hbox2_->addWidget(search_button_);
    hbox2_->addStretch(1);
    hbox2_->setSpacing(4);
    hbox2_->setMargin(6);

    //table + layout
    table_ = new QTableWidget(this);
    table_->setColumnCount(4);
    table_->setMinimumWidth(650);
    table_->setMinimumHeight(300);
    table_->setSelectionMode(QTableWidget::NoSelection);

    QStringList header;
    header << "" << "MS Experiment id" << "description" << "type";
    table_->QTableWidget::setHorizontalHeaderLabels(header);

    QVBoxLayout * vbox1_ = new QVBoxLayout(this);
    vbox1_->addLayout(hbox2_);
    vbox1_->insertWidget(-1, table_, 1);
    vbox1_->addLayout(hbox1_);

    //resize dialog
    adjustSize();
  }

  DBOpenDialog::~DBOpenDialog()
  {

  }

  void DBOpenDialog::ok()
  {
    for (Int col = 0; col < table_->rowCount(); ++col)
    {
      if (table_->item(col, 0)->checkState() == Qt::Checked)
      {
        result_.push_back(table_->item(col, 1)->text().toInt());
      }
    }
    emit accept();
  }

  void DBOpenDialog::loadSpectra()
  {
    stringstream query;
    query << "SELECT e.id,e.Description, count(s.id) FROM META_MSExperiment e right join DATA_Spectrum s on e.id=s.fid_MSExperiment WHERE";
    if (search_string_->text() != "")
    {
      query << " e.description like '%" << search_string_->text().toStdString() << "%' and ";
    }
    query << " s.MSLevel='1' GROUP BY e.id ORDER BY e.id ASC";
    QSqlQuery result = connection_.executeQuery(query.str());
    table_->setRowCount(result.size());
    UInt row = 0;
    QTableWidgetItem * item;
    while (result.isValid())
    {
      //id, description
      for (UInt col = 0; col < 2; col++)
      {
        item = new QTableWidgetItem(result.value(col).toString());
        table_->setItem(row, col + 1, item);
      }
      //type
      if (result.value(2).toInt() == 1)
      {
        item = new QTableWidgetItem("MS");
        table_->setItem(row, 3, item);
      }
      else
      {
        item = new QTableWidgetItem("HPLC-MS");
        table_->setItem(row, 3, item);
      }
      //checkboxes
      item = new QTableWidgetItem(QTableWidgetItem::Type);
      item->setCheckState(Qt::Unchecked);
      table_->setItem(row, 0, item);

      ++row;
      result.next();
    }
  }

}
