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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/DIALOGS/DBSpectrumSelectorDialog.h>
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
	DBSpectrumSelectorDialog::DBSpectrumSelectorDialog(DBConnection& adapter, vector<UnsignedInt>& result,QWidget* parent) 
	 : QDialog(parent), 
	 	 adapter_(adapter), 
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
		hbox1_->setSpacing(4);
		hbox1_->setMargin(6);
		
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
		table_ = new QTableWidget(0,4,this);
		table_->setMinimumWidth(650);
		table_->setMinimumHeight(300);
		//todo sort whole rows instead of one column
		//table_->setSorting(true);
		table_->setSelectionMode(QTableWidget::NoSelection);
		table_->horizontalHeaderItem(0)->setText("");
		table_->horizontalHeaderItem(1)->setText("MS Experiment id");
		table_->horizontalHeaderItem(2)->setText("description");
		table_->horizontalHeaderItem(3)->setText("type");
		QVBoxLayout* vbox1_ = new QVBoxLayout(this);
		vbox1_->addLayout(hbox2_);
		vbox1_->insertWidget(-1,table_,1);
		vbox1_->addLayout(hbox1_);
		
		//resize dialog
		adjustSize();
	}
	
	DBSpectrumSelectorDialog::~DBSpectrumSelectorDialog()
	{ 	
		
	}
	
	void DBSpectrumSelectorDialog::ok()
	{ 	
		for (SignedInt col=0;col<table_->rowCount();++col)
		{
			if(table_->item(col,0)->checkState()==Qt::Checked)
			{
				result_.push_back(table_->item(col,1)->text().toInt());
			}
		}
		emit accept();
	}
	
	void DBSpectrumSelectorDialog::loadSpectra()
	{ 
		stringstream query;
		query << "SELECT e.id,e.Description, count(s.id) FROM META_MSExperiment e right join DATA_Spectrum s on e.id=s.fid_MSExperiment WHERE";
		if(search_string_->text()!="")
		{
			query << " e.description like '%"<<search_string_->text().toAscii().data()<<"%' and ";
		}
		query << " s.MSLevel='1' GROUP BY e.id ORDER BY e.id ASC";
		QSqlQuery result;
		adapter_.executeQuery(query.str(),result);
		table_->setRowCount(result.size());
	 	UnsignedInt row=0;
	 	while(result.isValid())
		{
			//id, description
	 		for (unsigned int col = 0; col < 3; col++) 
	 		{ 
	      table_->item(row,col+1)->setText(result.value(col).toString());
	    }
	    //type
	    if (result.value(2).toInt()==1)
	    {
	    	table_->item(row,3)->setText("MS");
	    }
	    else
	    {
	    	table_->item(row,3)->setText("HPLC-MS");
	    }
	    //checkboxes
	    table_->item(row,0)->setCheckState(Qt::Unchecked);
	    ++row;
	    result.next();
		}
	}

}


