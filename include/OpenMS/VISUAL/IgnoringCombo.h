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
// $Id: IgnoringCombo.h,v 1.2 2006/03/03 17:55:35 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_VISUAL_IGNORINGCOMBO_H
#define OPENMS_VISUAL_IGNORINGCOMBO_H

#include <qcombobox.h>

namespace OpenMS
{
  
  /**
  	@brief This subclass of QComboBox ignores keyboard events
  
  
  	@ingroup Visual
  */
  class IgnoringCombo : public QComboBox
  {
    
    public:
      
      IgnoringCombo(QWidget* parent = 0, const char* name = 0);
      ~IgnoringCombo();
      
    protected:
      
      void keyReleaseEvent(QKeyEvent* e);
      void keyPressEvent(QKeyEvent* e);

  };
  
}
#endif //OPENMS_VISUAL_IGNORINGCOMBO_H

