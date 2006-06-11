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
// $Id: ProtDigMembers.h,v 1.3 2006/03/28 08:03:27 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer:  $
// --------------------------------------------------------------------------


#ifndef __PROTDIGMEMBERS_H__
#define __PROTDIGMEMBERS_H__

//STL-includes
#include <string>
#include <vector>
#include <ext/hash_map>

//My includes
#include "ProteinDigest.h"

//forward declarations
class AminoAcid;

//namespaces
using namespace std;


namespace OpenMS
  {

  //! this class contains pointers to some member variables of \c ProteinDigest for use in \c AminoAcidModification
  class ProtDigMembers
    {
    public:

      //! pointer to \c ProteinDigest::sequence_oneletter
      string* seq_oneletter;

      //! pointer to \c ProteinDigest::sequence_threeletter
      vector<string>* seq_threeletter;

      //! pointer to \c ProteinDigest::sequence_aminoacids
      vector<AminoAcid*>* seq_aminoacids;

      //! pointer to \c ProteinDigest::sequence_overall_modifications
      vector<int>* seq_overall_modifications;

      //! pointer to \c ProteinDigest::aminoacids_occurring
      __gnu_cxx::hash_map<string, AminoAcid*>* aa_occurring;

      //! pointer to \c ProteinDigest::aminoacids_positions
      __gnu_cxx::hash_map<string, vector<int>* >* aa_positions;

      //! pointer to \c ProteinDigest::fragments
      vector<pair<int,int> >* frags;

      //! pointer to \c ProteinDigest::sequence_fragments
      vector<vector<int> >* seq_fragments;

      //! pointer to \c ProteinDigest::cleavage_positions
      list<int>* cleav_positions;

      //! default constructor
      ProtDigMembers()
      {}
      ;

      //! constructor initializable with pointer to \c ProteinDigest since this class is friend, it can fill its members
      ProtDigMembers(ProteinDigest* dig)
      {
        seq_oneletter = &(dig->sequence_oneletter);
        seq_threeletter = &(dig->sequence_threeletter);
        seq_aminoacids = &(dig->sequence_aminoacids);
        seq_overall_modifications = &(dig->sequence_overall_modifications);
        aa_occurring = &(dig->aminoacids_occurring);
        aa_positions = &(dig->aminoacids_positions);
        frags = &(dig->fragments);
        seq_fragments = &(dig->sequence_fragments);
        cleav_positions = &(dig->cleavage_positions);
      }
    };
}
#endif

