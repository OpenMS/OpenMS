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


#ifndef __PROTEINDIGEST_H__
#define __PROTEINDIGEST_H__

//STL includes
#include <string>
#include <list>
#include <vector>
#include <list>
#include <fstream>
#include <sstream>
#include <ext/hash_map>

//my includes
#include "string_hash_stl_fixes.h"
#include "../config_specannotate.h"
#include "AminoAcid.h"
#include "Modification.h"
#include "Enzyme.h"
#include "MySQLAdapter.h"


class Sample;


using namespace std;

namespace OpenMS
  {

  class ProtDigMembers;

  /** This class represents a protein, its fragments and all further relevant information that is INDEPENDENT of annotation methods
   *  Classes Sample <-> Protein Digest:
   *  - class Sample contains all information relevant for annotating a spectrum that contains partially modified proteins or protein
   *    digests using the three methods "enumerate", "improved_enumerate" and "peakwise_cormen".
   *    It therefore also contains some information about fragments, positions and so on, as long as they are closely related to one of the
   *    three annotation methods and can only be used by them. Also the methods to "theoretically modify" the contained protein PARTIALLY
   *    (meaning modifications, that are not ALLWAYS realized at the SAME POSITION, modifications, that are simply POSSIBLE modifications
   *    like: at each of the positions 1,5 and 7 there could be one of following modifications...) are contained in class Sample, because 
   *    this theoretical "partial" modification process greatly depends on respectively is the INTEGRAL PART of the three methods mentioned
   *    above.
   * - class ProteinDigest on the other hand contains all aspects of a protein and its digest, that are INDEPENDENT of the three annotation
   *   methods in class Sample. 
   *   Only OVERALL modifications are included, such as Alkylation after a protein digest, where EVERY free -SH group (Cystein) should be
   *   (here is) modified. These modifications are ALWAYS there in a protein digest and INDEPENDENT of the annotation methods.
   * Because of this design, class ProteinDigest can be used for a very broad variety of purposes without complicating things by containing
   * stuff, that is ONLY usefull for one of the annotation methods above. 
   */
  class ProteinDigest
    {
    private:

      //! database information
      std::string db_username_, db_password_, db_host_;

      //! provides connection to MySQL database
      MySQLAdapter* sql_adapter_;

      //! contains unique identifier of protein
      string protein_identifier_;

      //! contains input filename of FASTA-File
      string protein_filename;

      //! contains ID of this instance of \c ProteinDigest
      int id;

      //! contains ID of this (undigested) Protein in Database
      int protein_ID_;


      //! ------ different representations of and stored information to sequence of protein to be digested -------------------------------------

      //! contains the sequence of the protein to be digested as string: each character is an amino acid in one letter code
      string sequence_oneletter;

      //! contains the sequence of the protein to be digested as vector of strings: each element is an amino acid in three letter code
      vector<string> sequence_threeletter;

      //! contains the sequence of the protein to be digested as vector of \c AminoAcid *
      vector<AminoAcid*> sequence_aminoacids;

      //! contains for each position in sequence a vector of indices in \c fragments , indicating in what fragments each residue occurs
      vector<vector<int> > sequence_fragments;

      //! conrains for each position in the sequence an entry with an (overall)modification ID. if it is 0 position is unmodified
      vector<int> sequence_overall_modifications;


      //! --------------------------------------------------------------------------------------------------------------------------------------


      //! contains fragments of protein after digest: start position and end position
      vector<pair<int,int> > fragments;

      //! contains indices in database (table digest_fragments) of fragments stored in \c fragments (fragment index -> database index)
      vector<int> fragment_database_indices_;

      //! contains all amino acids, that occur in \c protein : they are all only once instanciated (three letter codes are hashed!)
      __gnu_cxx::hash_map<string, AminoAcid*> aminoacids_occurring;

      //! contains to each occurring a set of sequence positions indicating on what position(s) this amino acid occurs
      __gnu_cxx::hash_map<string, vector<int>* > aminoacids_positions;

      //! after applying instance of \c Enzyme this variable is filled with sequence indices indicating cleavage positions (list: sortable)
      list<int> cleavage_positions;


      //! declare \c ProtDigMembers and \c Sample as friends of this class, so that they can access its members
      friend class ProtDigMembers;
      friend class Sample;


    public:


      //! exception thrown if no fasta-filename is present in database
    class NoProteinFilename : public OpenMS::Exception::Base
        {
        public:
          NoProteinFilename(const char* file, int line, const char* function, std::string request) throw();
          ~NoProteinFilename() throw();
        };

      //! exception thrown if wrong position in filename is given
    class WrongPositionInProtein : public OpenMS::Exception::Base
        {
        public:
          WrongPositionInProtein(const char* file, int line, const char* function, std::string method) throw();
          ~WrongPositionInProtein() throw();
        };

      //! exception thrown if given Modification cannot modify residue at given Position
    class WrongModification : public OpenMS::Exception::Base
        {
        public:
          WrongModification(const char* file, int line, const char* function, int mod_ID, int pos) throw();
          ~WrongModification() throw();
        };


      //! default constructor
      ProteinDigest();

      //! constructor initializable with unique identifier of protein (must be stored in database incl filename(PDB/FASTA) absolute)
      ProteinDigest(string identifier, int ID = 0, std::string db_username="",
                    std::string db_password="", std::string db_host="");

      //! copy constructor
      ProteinDigest(const ProteinDigest& protein_digest);

      //! assignment operator
      const ProteinDigest& operator=(const ProteinDigest& protein_digest);

      /** initializes different representations of sequence of protein to be digested in synchronous way
       * fills \c protein , \c sequence_oneletter , \c sequence_threeletter and \c sequence_aminoacids with the right values.
       * sets \c sequence_modifications to be of the right size and fills positions with NULL.
       * sets \c sequence_masses and \c sequence_fragments to be of the right size.
       * this method futhermore fills \c aminoacids_occurring and \c aminacids_positions with the right values
       */
      void initialize() throw();

      //! via this method instances of \c Enzyme can be applied to contents of this class
      void digest(Enzyme* enz);

      /** via this method instances of \c Modification can be applied to contents of this class as OVERALL modifications
       *  i.e.all possible positions are see as modified.
       *
       */
      void modifyOverall(Modification* mod);

      /** the modifications are added as chained list in database Table realized_modifications
       *  the field next_realized_modification of LAST element of that list is NULL
       *  the FIRST id (primary key of database entry) of this chained list ist returned
       *  so the last entry of all the different combinations of the partial modifications can contain in its field next_realized_modification
       *  the id of first entry of this chain of overall modifications. 
       *  by doing so, the overall modifications have to be saved only once for each protein_modification_scenario,
       *  and not once for each of the modification combinations in a scenario
       *
       * important: this method does NOT check, whether actual scenario is already present in database. ensure that before calling it.
       */
      int dbStoreOverallModifications();

      //! returns pointer to fragment with index \c index
      pair<int,int> getFragment(int index)
      {
        return fragments[index];
      };

      /** returns list of vectors of ints:
       * first elements of vectors are indices in \c fragments of that fragments, that contain certain a site (given in one-letter-code)
       * second elements of vectors are positions of sites
       * third elements of vectors are ID's of this instance of \c ProteinDigest
       */
      list<vector<int> > getFragmentIndicesContaining(const vector<string>& sites);

      //! returns mass of fragment given by start and end position. start and end aa's are considered to be N-/C-terminal
      double getFragmentMonoMass(int start_pos, int end_pos);

      //! returns mass of fragment given by start and end position. start and end aa's are considered to be N-/C-terminal
      double getFragmentAverageMass(int start_pos, int end_pos);

      //! returns filename of this instance of \c ProteinDigest
      string getFilename()
      {
        return protein_filename;
      };

      //! returns name of residue specified (in three letter code)
      string getResName(int pos)
      {
        return sequence_threeletter[pos];
      };

      //! returns identifier of protein
      string getProteinIdentifier()
      {
        return protein_identifier_;
      };

      //! returns ID of this (undigested) protein in database
      int getProteinID()
      {
        return protein_ID_;
      };

      //! returns length of Protein
      int getProteinLength()
      {
        return sequence_threeletter.size();
      };

      //! returns mass of overall modified fragment, stores overall modifications and multiplicities in \c temp_o_mods
      double getFragmentOverallModifiedMass(int start_pos, int end_pos, string masstype, __gnu_cxx::hash_map<int,int>& temp_o_mods,
                                            Modification& mod);
    };

}
#endif


