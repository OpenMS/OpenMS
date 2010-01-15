// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: $
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
	DBOpenDialog::DBOpenDialog(DBConnection& connection, vector<UInt>& result, QWidget* parent) 
	 : QDialog(parent), 
	 	 connection_(connection), 
	 	 result_(result)
	{ 	
		setWindowTitle("Select spectra from DB to open");
		
		//OK+Cancel button + layout
		QPushButton* cancel_button_ = new QPushButton("&Cancel",this);
		connect(cancel_button_,SIGNAL(clicked()),this,SLOT(reject()));
		QPushButton* ok_button_ = new QPushButton("&OK",this);
		connect(ok_button_,SIGNAL(clicked()),this,SLOT(ok()));
	
		QHBoxLayout* hbox1_ = new QHBoxLayout();
		hbox1_->addStretch(1);
		hbox1_->addWidget(ok_button_);
		hbox1_->addWidget(cancel_button_);
		
		//string line edit, button + layout
		QLabel* label_ = new QLabel("Description contains:",this);
		search_string_ = new QLineEdit(this);
		QPushButton* search_button_ = new QPushButton("refresh",this);
		connect(search_button_,SIGNAL(clicked()),this,SLOT(loadSpectra()));
		
		QHBoxLayout* hbox2_ = new QHBoxLayout();
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
		
		QVBoxLayout* vbox1_ = new QVBoxLayout(this);
		vbox1_->addLayout(hbox2_);
		vbox1_->insertWidget(-1,table_,1);
		vbox1_->addLayout(hbox1_);
		
		//resize dialog
		adjustSize();
	}
	
	DBOpenDialog::~DBOpenDialog()
	{ 	
		
	}
	
	void DBOpenDialog::ok()
	{ 	
		for (Int col=0;col<table_->rowCount();++col)
		{
			if(table_->item(col,0)->checkState()==Qt::Checked)
			{
				result_.push_back(table_->item(col,1)->text().toInt());
			}
		}
		emit accept();
	}
	
	void DBOpenDialog::loadSpectra()
	{ 
		stringstream query;
		query << "SELECT e.id,e.Description, count(s.id) FROM META_MSExperiment e right join DATA_Spectrum s on e.id=s.fid_MSExperiment WHERE";
		if(search_string_->text()!="")
		{
			query << " e.description like '%"<<search_string_->text().toStdString() <<"%' and ";
		}
		query << " s.MSLevel='1' GROUP BY e.id ORDER BY e.id ASC";
		QSqlQuery result = connection_.executeQuery(query.str());
		table_->setRowCount(result.size());
	 	UInt row=0;
	 	QTableWidgetItem* item;
	 	while(result.isValid())
		{
			//id, description
	 		for (UInt col = 0; col < 2; col++) 
	 		{ 
	 		 	item = new QTableWidgetItem(result.value(col).toString());
     		table_->setItem(row, col+1, item);
	    }
	    //type
	    if (result.value(2).toInt()==1)
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


