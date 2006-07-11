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
//
// --------------------------------------------------------------------------
// $Maintainer:  $
// --------------------------------------------------------------------------


#ifndef __CUSTOMEVENTS_H__
#define __CUSTOMEVENTS_H__


namespace OpenMS
  {

  //Event that is sent by this thread to update its output-window in the gui-thread. in QT this is thread-safe
  class OutputEvent : public QCustomEvent
    {
    public:
      //initialize QCustomEvent with an ID large as the value of the "User" entry in the QEvent::Type enum.
      OutputEvent(std::string output) : QCustomEvent(65432), outp(output)
      {}
      ;
      QString output() const
        {
          return outp;
        };

    private:
      std::string outp;
    };


  //Event that is sent by this thread to notify the GUI-thread, that it is finished. in QT this is thread-safe
  class FinishEvent : public QCustomEvent
    {
    public:
      //initialize QCustomEvent with an ID large as the value of the "User" entry in the QEvent::Type enum.
      FinishEvent() : QCustomEvent(65433)
      {}
      ;
    };

}
#endif
