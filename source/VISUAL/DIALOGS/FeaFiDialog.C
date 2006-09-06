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
// $Id: FeaFiDialog.C,v 1.4 2006/06/09 14:46:55 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

// OpenMS includes
#include <OpenMS/VISUAL/DIALOGS/FeaFiDialog.h>
#include <sstream>
#include <qfiledialog.h>
#include <qpushbutton.h>
#include <qlabel.h>

using namespace std;
using namespace OpenMS;

FeaFiDialog::FeaFiDialog( QWidget * parent, const char * name, WFlags fl):
	FeaFiDialogTemplate(parent,name,fl), finder_()
{
	resize(sizeHint());
	start_button->setEnabled(false);
}

FeaFiDialog::~FeaFiDialog()
{

}

FeatureFinder& FeaFiDialog::getFeatureFinder()
{
	return finder_;
}


void FeaFiDialog::loadParamFile()
{
	QString dir = QDir::current().path();
	dir = dir.left(dir.find("OpenMS")).append("OpenMS/source/APPLICATIONS/FEATUREFINDER/");
 	QString file = QFileDialog::getOpenFileName(dir, "Parameters (*.ini)", this,
																							"featurefinder dialog","Select file(s) to open");

	if (!file.isEmpty() && file!="")
	{
		Param param;
		param.load(file.ascii());
		bool isParamValid = finder_.setParam(param);
		if (isParamValid)
		{
			std::ostringstream label;
			label << finder_;
			param_label->setText(label.str().c_str());
			param_label->setAlignment( int( QLabel::AlignVCenter|QLabel::AlignLeft ) );
			start_button->setEnabled(true);
			start_button->setFocus();
		}
		else
		{
			param_label->setText(tr("Chosen param file was not valid.\nNo Feature-Finder loaded yet."));
		}	
	}
}
