// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Lars Nilse $
// $Authors: Steffen Sass, Holger Plattfaut $
// --------------------------------------------------------------------------

#ifndef OPENMS_FILTERING_DATAREDUCTION_SILACFILTER_H
#define OPENMS_FILTERING_DATAREDUCTION_SILACFILTER_H

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/DATASTRUCTURES/DataPoint.h>
#include <OpenMS/FILTERING/DATAREDUCTION/SILACFiltering.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>
#include <queue>
#include <list>

namespace OpenMS
{
   class SILACFiltering;
   class DataPoint;

   /**
   * @brief Filter to use for SILACFiltering
   *
   * A SILACFilter searches for SILAC patterns, which correspond to the defined mass shifts and charge.
   * Only peaks are taken into account, which were not blacklisted by other filters before e.g. are not part of
   * a SILAC pair yet.
   *
   * @see SILACFiltering
   */

   class OPENMS_DLLAPI SILACFilter {
     friend class SILACFiltering;

  private:

  /**
   * @brief number of peptides [i.e. number of labelled peptides +1, e.g. for SILAC triplet =3]
   */
   Int number_of_peptides_;

  /**
   * @brief charge of the ions to search for
   */
   Int charge_;

  /**
   * @brief number of peaks per peptide to search for
   */
   Int isotopes_per_peptide_;

  /**
   * @brief mass shift(s) in [Da] to search for
   */
   std::vector<DoubleReal> mass_separations_;

  /**
   * @brief peak positions of SILAC pattern
   */
   std::vector<DoubleReal> peak_positions_;

  /**
  * @brief m/z separtion between individual peptides [e.g. {0 Th, 4 Th, 5 Th}]
  */
  std::vector<DoubleReal> mz_peptide_separations_;

  /**
  * @brief m/z shifts relative to mono-isotopic peak of unlabelled peptide
  */
  std::vector<DoubleReal> expected_mz_shifts_;

  /**
   * @brief distance between isotopic peaks of a peptide in [Th]
   */
   DoubleReal isotope_distance_;

  /**
   * @brief holds the recognized features
   */
   std::vector<DataPoint> elements_;

  /**
   * @brief maximal value of which a predicted SILAC feature may deviate from the averagine model
   */
   DoubleReal model_deviation_;

  /**
   * @brief m/z at which the filter is currently applied to
   */
   DoubleReal current_mz_;

  /**
   * @brief exact m/z shift of isotopic peaks in a SILAC pattern relative to the mono-isotopic peak of the light peptide, peptides (row) x isotope (column)
   */
   std::vector<std::vector<DoubleReal> > exact_shifts_;

  /**
   * @brief intensities at mz + exact_shifts in a SILAC pattern, where mz is the m/z of the mono-isotopic peak of light peptide
   */
   std::vector<std::vector<DoubleReal> > exact_intensities_;

  /**
   * @brief expected m/z shift of isotopic peaks in a SILAC pattern relative to the mono-isotopic peak of the light peptide, peptides (row) x isotope (column)
   */
   std::vector<std::vector<DoubleReal> > expected_shifts_;

  /**
   * @brief Computes the actual m/z shift between the position mz and a region about expectedMzShift away. Returns -1 if there is no correlation between mz and signal in interval [mz + expectedMzShift - maxMzDeviation, mz + expectedMzShift + maxMzDeviation].
   * [e.g. from theoretical considerations we expect a good correlation between peaks 4.02 Th apart. But the shift between signals actually observed is 4.0189 Th.]
   * @param mz m/z position of the reference signal [e.g. mono-isotopic peak of the light peptide]
   * @param expectedMzShift poitive m/z shift at which we would expect a correlating signal [e.g. 4.02 Th]
   * @param maxMzDeviation maximum allowed deviation between expected and actual shift [In the above example the shift is 0.0011 Th.]
   */
   DoubleReal computeActualMzShift(DoubleReal mz, DoubleReal expectedMzShift, DoubleReal maxMzDeviation);

  /**
   * @brief returns true if there exists a SILAC feature at the given position, which corresponds to the filter's properties
   * @param rt RT value of the position
   * @param mz m/z value of the position
   */
   bool isSILACPattern(DoubleReal rt, DoubleReal mz);  

   public:  

  /**
   * @brief default constructor
   */
   SILACFilter();

  /**
   * @brief detailed constructor for SILAC pair filtering
   * @param mass_separations all mass shifts of the filter
   * @param charge charge of the ions to search for
   * @param model_deviation maximum deviation from the averagine model
   * @param isotopes_per_peptide number of peaks per petide to search for
   */
   SILACFilter(std::vector<DoubleReal> mass_separations, Int charge, DoubleReal model_deviation, Int isotopes_per_peptide);

  /**
   * @brief destructor
   */
   virtual ~SILACFilter();

  /**
   * @brief returns the predicted peak width at position mz
   * @param mz mz position of the peak
   */
   static DoubleReal getPeakWidth(DoubleReal mz);

  /**
   * @brief gets the m/z values of all peaks , which belong the last identified feature
   */
   std::vector<DoubleReal> getPeakPositions();

  /**
   * @brief gets the m/z shifts relative to mono-isotopic peak of unlabelled peptide
   */
   std::vector<DoubleReal> getExpectedMzShifts();

  /**
   * @brief returns all identified elements
   */
   std::vector<DataPoint> getElements();

  /**
   * @brief returns the charge of the filter
   */
   Int getCharge();

  /**
   * @brief returns the mass shifts of the filter in [Da]
   */
   std::vector<DoubleReal> getMassSeparations();

  };
}

#endif /* SILACFILTER_H_ */
