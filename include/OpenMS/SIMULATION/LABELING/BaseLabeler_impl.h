// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Stephan Aiche $
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_SIMULATION_LABELING_BASELABELER_IMPL_H
#define OPENMS_SIMULATION_LABELING_BASELABELER_IMPL_H

#include <OpenMS/SIMULATION/LABELING/BaseLabeler.h>

// include derived classes here
#include <OpenMS/SIMULATION/LABELING/ITRAQLabeler.h>
#include <OpenMS/SIMULATION/LABELING/LabelFreeLabeler.h>
#include <OpenMS/SIMULATION/LABELING/O18Labeler.h>

namespace OpenMS
{

    void BaseLabeler::registerChildren()
    {
        Factory< BaseLabeler >::registerProduct(LabelFreeLabeler::getProductName(), &LabelFreeLabeler::create);
        Factory< BaseLabeler >::registerProduct(O18Labeler::getProductName(), &O18Labeler::create);
        Factory< BaseLabeler >::registerProduct(ITRAQLabeler::getProductName(), &ITRAQLabeler::create);
        return;
    }

} // namespace OpenMS

#endif // OPENMS_SIMULATION_LABELING_BASELABELER_IMPL_H

