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
// $Id: Formulas.h,v 1.1 2006/06/09 14:08:43 andreas_bertsch Exp $
// $Author: andreas_bertsch $
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------
//
#ifndef OPENMS_CHEMISTRY_FORMULAS_H
#define OPENMS_CHEMISTRY_FORMULAS_H


// TODO oder array in EmpiricalFormula?

#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>

namespace OpenMS
{
	/** 
			@ingroup Chemistry
	
			@brief some often used empirical formulas
	*/
	
	namespace Formulas
	{	
		static const EmpiricalFormula& H()
		{
			static const EmpiricalFormula H("H");
			return H;
		}

		static const EmpiricalFormula& H2O()
		{
			static const EmpiricalFormula H2O("H2O");
			return H2O;
		}

		static const EmpiricalFormula& Water()
		{
			return H2O();
		}
		
		static const EmpiricalFormula& NH()
		{
			static const EmpiricalFormula NH("NH");
			return NH;
		}

		static const EmpiricalFormula& OH()
		{
			static const EmpiricalFormula OH("OH");
			return OH;
		}

		static const EmpiricalFormula& NH3()
		{
			static const EmpiricalFormula NH3("NH3");
			return NH3;
		}

	}

}

#endif
