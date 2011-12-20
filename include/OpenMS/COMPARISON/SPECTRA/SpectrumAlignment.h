// -*- mode: C++; tab-width: 2; -*-
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
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------
//
#ifndef OPENMS_COMPARISON_SPECTRA_SPECTRUMALIGNMENT_H
#define OPENMS_COMPARISON_SPECTRA_SPECTRUMALIGNMENT_H

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>

#include <vector>
#include <map>
#include <utility>
#include <algorithm>

#define ALIGNMENT_DEBUG
#undef  ALIGNMENT_DEBUG

namespace OpenMS
{

  /**
	  @brief Aligns the peaks of two spectra
		
    Using a banded (width via 'tolerance' parameter) alignment.
    Scoring function is the m/z distance between peaks. Intensity does not play a role!

		@htmlinclude OpenMS_SpectrumAlignment.parameters
		
		@ingroup SpectraComparison
  */
	
  class OPENMS_DLLAPI SpectrumAlignment 
  	: public DefaultParamHandler
  {
  public:
	
		// @name Constructors and Destructors
		// @{
    /// default constructor
    SpectrumAlignment();

    /// copy constructor
    SpectrumAlignment(const SpectrumAlignment& source);

    /// destructor
    virtual ~SpectrumAlignment();

    /// assignment operator
    SpectrumAlignment& operator = (const SpectrumAlignment& source);
		// @}

		template <typename SpectrumType>
		void getSpectrumAlignment(std::vector<std::pair<Size, Size> >& alignment, const SpectrumType& s1, const SpectrumType& s2) const
		{
      if (!s1.isSorted () || !s2.isSorted())
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Input to SpectrumAlignment is not sorted!");
      }
      // clear result
      alignment.clear();

			double tolerance = (double)param_.getValue("tolerance");
			std::map<Size, std::map<Size, std::pair<Size, Size> > > traceback;
			std::map<Size, std::map<Size, double> > matrix;
			
			// init the matrix with "gap costs" tolerance
			matrix[0][0] = 0;
			for (Size i = 1; i <= s1.size(); ++i)
			{
				matrix[i][0] = i * tolerance;
				traceback[i][0]  = std::make_pair(i - 1, 0);
			}
			for (Size j = 1; j <= s2.size(); ++j)
			{
				matrix[0][j] = j * tolerance;
				traceback[0][j] = std::make_pair(0, j - 1);
			}
			
			// fill in the matrix
			Size left_ptr(1);
			Size last_i(0), last_j(0);

			//Size off_band_counter(0);
			for (Size i = 1; i <= s1.size(); ++i)
			{
        double pos1(s1[i - 1].getMZ());
        
				for (Size j = left_ptr; j <= s2.size(); ++j)
				{
					bool off_band(false);
					// find min of the three possible directions
					double pos2(s2[j - 1].getMZ());
					double diff_align = fabs(pos1 - pos2);
	
					// running off the right border of the band?
					if (pos2 > pos1 && diff_align > tolerance)
					{
						if (i < s1.size() && j < s2.size() && s1[i].getMZ() < pos2)
						{
							off_band = true;
						}
					}
	
					// can we tighten the left border of the band?
					if (pos1 > pos2 && diff_align > tolerance && j > left_ptr + 1)
					{
						++left_ptr;
					}
	
					double score_align = diff_align;
					
					if (matrix.find(i - 1) != matrix.end() && matrix[i - 1].find(j - 1) != matrix[i - 1].end())
					{
						score_align += matrix[i - 1][j - 1];
					}
					else
					{
						score_align += (i - 1 + j - 1) * tolerance;
					}
	
					double score_up = tolerance;
					if (matrix.find(i) != matrix.end() && matrix[i].find(j - 1) != matrix[i].end())
					{
						score_up += matrix[i][j - 1];
					}
					else
					{
						score_up += (i + j - 1) * tolerance;
					}
					
					double score_left = tolerance;
					if (matrix.find(i - 1) != matrix.end() && matrix[i - 1].find(j) != matrix[i - 1].end())
					{
						score_left += matrix[i - 1][j];
					}
					else
					{
						score_left += (i - 1 + j) * tolerance;
					}
			
					#ifdef ALIGNMENT_DEBUG
					cerr << i << " " << j << " " << left_ptr << " " << pos1 << " " << pos2 << " " << score_align << " " << score_left << " " << score_up << endl;
					#endif
					
					if (score_align <= score_up && score_align <= score_left && diff_align <= tolerance)
					{
						matrix[i][j] = score_align;
						traceback[i][j] = std::make_pair(i - 1, j - 1);
						last_i = i;
						last_j = j;
					}
					else
					{
						if (score_up <= score_left)
						{
							matrix[i][j] = score_up;
							traceback[i][j] = std::make_pair(i, j - 1);
						}
						else
						{
							matrix[i][j] = score_left;
							traceback[i][j] = std::make_pair(i - 1, j);
						}
					}
	
					if (off_band)
					{
						break;
					}
				}
			}
	
			//last_i = s1.size() + 1;
			//last_j = s2.size() + 1;
	
			//cerr << last_i << " " << last_j << endl;
	
#ifdef ALIGNMENT_DEBUG
#if 0
			cerr << "TheMatrix: " << endl << " \t  \t";
			for (Size j = 0; j != s2.size(); ++j)
			{
				cerr << s2[j].getPosition()[0] << " \t";
			}
			cerr << endl;
			for (Size i = 0; i <= s1.size(); ++i)
			{
				if (i != 0)
				{
					cerr << s1[i - 1].getPosition()[0] << " \t";
				}
				else
				{
					cerr << " \t";
				}
				for (Size j = 0; j <= s2.size(); ++j)
				{
					if (matrix.has(i) && matrix[i].has(j))
					{
						if (traceback[i][j].first == i - 1 && traceback[i][j].second == j - 1)
						{
							cerr << "\\";
						}
						else
						{
							if (traceback[i][j].first == i - 1 && traceback[i][j].second == j)
							{
								cerr << "|";
							}
							else
							{
								cerr << "-";
							}
						}
	
						cerr << matrix[i][j] << "  \t";
					}
					else
					{
						cerr << "-1  \t";
					}
				}
				cerr << endl;
			}
#endif
#endif
	
			// do traceback
			Size i = last_i;
			Size j = last_j;
	
			while (i >= 1 && j >= 1)
			{
				if (traceback[i][j].first == i - 1 && traceback[i][j].second == j - 1)
				{
					alignment.push_back(std::make_pair(i - 1, j - 1));
				}
				Size new_i = traceback[i][j].first;
				Size new_j = traceback[i][j].second;
	
				i = new_i;
				j = new_j;
			}
	
			std::reverse(alignment.begin(), alignment.end());
	
#ifdef ALIGNMENT_DEBUG
#if 0
			// print alignment
			cerr << "Alignment (size=" << alignment.size() << "): " << endl;
			
			Size i_s1(0), i_s2(0);
			for (vector<pair<Size, Size> >::const_reverse_iterator it = alignment.rbegin(); it != alignment.rend(); ++it, ++i_s1, ++i_s2)
			{
				while (i_s1 < it->first - 1)
				{
					cerr << i_s1 << " " << s1[i_s1].getPosition()[0] << " " << s1[i_s1].getIntensity() << endl;
					i_s1++;
				}
				while (i_s2 < it->second - 1)
				{
					cerr << " \t " <<  i_s2 << " " << s2[i_s2].getPosition()[0] << " " << s2[i_s2].getIntensity() << endl;
					i_s2++;
				}
				cerr << "(" << s1[it->first - 1].getPosition()[0] << " <-> "<< s2[it->second - 1].getPosition()[0] << ") (" 
				<< it->first << "|" << it->second << ") (" << s1[it->first - 1].getIntensity() << "|" << s2[it->second - 1].getIntensity() << ")" << endl;
			}
#endif
#endif
		
		}
  };

}
#endif //OPENMS_COMPARISON_SPECTRA_SPECTRUMALIGNMENT_H
