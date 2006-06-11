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
// $Id: SpecView_SpecAnnotate.C,v 1.4 2006/05/30 15:46:40 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#include <qapplication.h>

#include "./SPECANNOTATE/VISUAL/SpectrumMDIWindowEnhanced.h"
#include <qwindowsstyle.h>


int main( int argc, char ** argv )
{
  QApplication a( argc, argv );
  OpenMS::SpectrumMDIWindowEnhanced* mw = OpenMS::SpectrumMDIWindowEnhanced::getInstance();
  a.setMainWidget(mw);
  mw->setCaption( "TOPPView with SpecAnnotate" );
  mw->show();
  mw->loadCommandLineFiles();

  a.connect( &a, SIGNAL(lastWindowClosed()), &a, SLOT(quit()) );

  int res = a.exec();

  return res;
}




