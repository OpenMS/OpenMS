// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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
#include <qpushbutton.h>
#include <qlayout.h>
#include <qlineedit.h>
#include <qtable.h>
#include <qlabel.h>

#include <sstream>

using namespace std;


namespace OpenMS
{
	DBSpectrumSelectorDialog::DBSpectrumSelectorDialog(DBConnection& adapter, vector<UnsignedInt>& result,QWidget* parent, const char* name) 
	 : QDialog(parent,name), 
	 	 adapter_(adapter), 
	 	 result_(result)
	{ 	
		setCaption("Select spectra from DB to open");
		
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
		table_ = new QTable(0,4,this);
		table_->setMinimumWidth(650);
		table_->setMinimumHeight(300);
		//todo sort whole rows instead of one column
		//table_->setSorting(true);
		table_->setSelectionMode(QTable::NoSelection);
		table_->setRowMovingEnabled(false);
		table_->setColumnReadOnly(1,true);
		table_->setColumnReadOnly(2,true);
		table_->setColumnReadOnly(3,true);
		table_->setLeftMargin(0);
		table_->horizontalHeader()->setLabel(0, "",20);
		table_->horizontalHeader()->setLabel(1, "MS Experiment id",100);
		table_->horizontalHeader()->setLabel(2, "description",400);
		table_->horizontalHeader()->setLabel(3, "type",80);
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
		for (SignedInt col=0;col<table_->numRows();++col)
		{
			if(dynamic_cast<QCheckTableItem*>(table_->item(col,0))->isChecked())
			{
				result_.push_back(atoi(table_->text(col,1).ascii()));
			}
		}
		emit accept();
	}
	
	void DBSpectrumSelectorDialog::loadSpectra()
	{ 
		stringstream query;
		query << "SELECT e.id,e.Description, count(s.id) FROM MSExperiment e right join Spectrum s on e.id=s.MSExperiment_id WHERE";
		if(search_string_->text()!="")
		{
			query << " e.description like '%"<<search_string_->text().ascii()<<"%' and ";
		}
		query << " s.MS_Level='1' group by e.id order by e.Description";
		QSqlQuery result;
		adapter_.executeQuery(query.str(),result);
		table_->setNumRows(result.size());
	 	table_->setLeftMargin(0);
	 	UnsignedInt row=0;
	 	while(result.next())
		{
			//id, description
	 		for (unsigned int col = 0; col < 3; col++) 
	 		{ 
	      table_->setText(row,col+1,result.value(col).toString());
	    }
	    //type
	    if (result.value(2).asInt()==1)
	    {
	    	table_->setText(row,3,"MS");
	    }
	    else
	    {
	    	table_->setText(row,3,"HPLC-MS");
	    }
	    //checkboxes
	    QCheckTableItem* checker = new QCheckTableItem(table_,QString());
	    table_->setItem(row,0,checker);
	    ++row;
		}
	}

}


