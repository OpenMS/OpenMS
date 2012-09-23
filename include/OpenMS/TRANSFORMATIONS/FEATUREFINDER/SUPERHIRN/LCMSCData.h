// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Florian Zeller $
// $Authors: Lukas Mueller, Markus Mueller $
// --------------------------------------------------------------------------
//
///////////////////////////////////////////////////////////////////////////
//
//  PEAK DETECTION OF FOURIER TRANSFORME MS INSTRUMENT DATA
//
//  written by Markus Mueller, markus.mueller@imsb.biol.ethz.ch
//  and Lukas Mueller, Lukas.Mueller@imsb.biol.ethz.ch
//  October 2005
//
//  Ported to OpenMS by Florian Zeller, florian.zeller@bsse.ethz.ch
//  December 2010
//
//  Group of Prof. Ruedi Aebersold, IMSB, ETH Hoenggerberg, Zurich
//
//


#ifndef _LCMSCData_h_
#define _LCMSCData_h_

namespace OpenMS
{

  typedef std::map<int, LCElutionPeak> elution_peak_list;
  typedef std::map<int, LCElutionPeak>::iterator elution_peak_list_ITERATOR;
  typedef std::map<double, elution_peak_list> MZ_LIST;
  typedef MZ_LIST::iterator MZ_LIST_ITERATOR;


  class OPENMS_DLLAPI LCMSCData
  {
public:

    LCMSCData();
    LCMSCData(const LCMSCData &);
    LCMSCData(const LCMSCData *);
    virtual ~LCMSCData();

    //////////////////////////////////////////////////
    // overload operators:
    LCMSCData & operator=(const LCMSCData &);

    ///////////////////////////////
    // start here all the get / set
    // function to access the
    // variables of the class

    // start / end of the mz list:
    MZ_LIST_ITERATOR get_DATA_start(){return DATA.begin(); }
    MZ_LIST_ITERATOR get_DATA_end(){return DATA.end(); }

    // find data of a specific m/z:
    MZ_LIST_ITERATOR get_MZ_by_iterator(double);
    // add data into the structure:
    void add_LC_elution_peak(double, LCElutionPeak *);

    // get a list of m/z observed in a scan +/- scan_tolerance:
    // return the area of the LC elution peaks
    std::vector<LCElutionPeak> get_MZ_list(int);
    std::vector<LCElutionPeak> get_MZ_list(int, int);

    // get all extracted LC peaks:
    std::vector<LCElutionPeak *> get_ALL_peak();
    // get a vector with all LC peaks ordered by their score:
    std::vector<LCElutionPeak *> get_ALL_peak_ordered();


    /*
    //////////////////////////////////////////////
    // operator which compares LC elution peaks
    // by their score:
    struct OPERATOR_LC_PEAK{
      bool operator()(LCElutionPeak* A,LCElutionPeak* B) const{
        return A->get_LC_score() < B->get_LC_score();
      }
    };
  */

    // minimal MS peak intensity allowed:
    //  static float intensity_min_threshold;



protected:

    // main data struture:
    MZ_LIST DATA;


  };

} // ns

#endif
