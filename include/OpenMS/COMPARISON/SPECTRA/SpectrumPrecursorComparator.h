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
// $Id: SpectrumPrecursorComparator.h,v 1.4 2006/05/30 15:46:38 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------
//
#ifndef OPENMS_COMPARISON_SPECTRA_SPECTRUMPRECURSORCOMPARATOR_H
#define OPENMS_COMPARISON_SPECTRA_SPECTRUMPRECURSORCOMPARATOR_H
#include <OpenMS/COMPARISON/SPECTRA/CompareFunctor.h>
#include <OpenMS/KERNEL/MSSpectrum.h>

namespace OpenMS
{

  /**
  SpectrumPrecursorComparator compares just the parent mass of two spectra<br>
  \param window the mass window where spectrum pairs get similarity > 0<br>
  */
  class SpectrumPrecursorComparator : public CompareFunctor
  {
  public:
    /** @brief standard constructor <br> */
    SpectrumPrecursorComparator();

    /** @brief copy constructor <br> */
    SpectrumPrecursorComparator(const SpectrumPrecursorComparator& source);

    /** @brief destructor <br> */
    ~SpectrumPrecursorComparator();

    /** @brief assignment operator <br> */
    SpectrumPrecursorComparator& operator=(const SpectrumPrecursorComparator& source);

    static FactoryProduct* create() { return new SpectrumPrecursorComparator();}

		static const String getName()
		{
			return "SpectrumPrecursorComparator";
		}

    String info() const;

    double operator()(const ClusterSpectrum& csa,const ClusterSpectrum& csb)const;

  private:
    static const String info_;
  };

}

#endif //OPENMS_COMPARISON_SPECTRA_SPECTRUMPRECURSORCOMPARATOR_H
