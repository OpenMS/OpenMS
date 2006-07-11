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


#ifndef __SAMPLE_H__
#define __SAMPLE_H__

//STL includes
#include <vector>
#include <string>
#include <map>
#include <ext/hash_map>

//OpenMS includes
#include <OpenMS/KERNEL/DPeakArray.h>
#include <OpenMS/KERNEL/DPeak.h>
#include <OpenMS/VISUAL/Spectrum1DWidget.h>

//MyIncludes
#include "MySQLAdapter.h"
#include "Enzyme.h"
#include "Modification.h"
#include "ProteinDigest.h"
#include "Annotation.h"



namespace OpenMS
  {
  //forward declarations
  class Annotate;
  class OutputEvent;
}


using namespace OpenMS;



/** central class representing contents of protein sample and all related information
 *  IMPORTANT: this class represents ONE protein. if the sample (and therefore the spectrum) results from multiple proteins, multiple 
 *             instances of this class have to be used!
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
class Sample
  {
  private:

    /*****************************************************************************************************************************************
     * Used by more than one method 
     ****************************************************************************************************************************************/

    //VARIABLES: ------------------------------------------------------

    //! if compiled with QT support: for sending events to gui-thread
    Annotate* qannotate_;

    //! if compiled with QT support: for sending output events to gui-thread
    OutputEvent* oe_;

    //! most of the necessary data about actual sample, except for overall modifications: NEEDED FOR CALL OF "initialize_()"
    __gnu_cxx::hash_map<std::string, std::string> sample_data_;

    //! contains strings of overall modifications of actual sample: NEEDED FOR CALL OF "initialize_()"
    std::vector<std::string> overall_mod_strings_;

    //! database information
    std::string db_username_, db_password_, db_host_;

    //! provides connection to MySQL database
    MySQLAdapter* sql_adapter_;

    //! contains name of annotation_method to be used
    std::string annotation_method;

    //! contains, after execution of modify() or existsProtModScenInDB(), ID of this Sample's Protein Modification Scenario
    int prot_mod_scen_ID;

    //! contains, either after call of dbRegister() or after call of existsInDB(), this Sample's ID in table SAMPLE_TABLE
    int sample_ID;

    //! contains value for search range for search of "real" masses in calculated ones
    double range;

    //! contains either "mono" or "average", specifies what masstype should be used for annotation
    std::string masstype;

    //! contains "real" peaklist
    OpenMS::DPeakArray<1, OpenMS::DPeak<1> > peaklist_;

    //! contains iterators to instances of peaks in TOPPView, that are in peaklist
    vector<OpenMS::DPeakArray<1, OpenMS::DPeak<1> >::iterator> external_peaklist_;

    //! contains filename of "real" peaklist
    string peakfile;

    //! contains format type of peakfile
    string peakfile_format;

    //! contains output-directory
    string outputdir;

    //! contains used enzyme
    Enzyme* enzyme;

    //! flag: if true, \c calculateAnnotations has to take fragments into account, else only whole protein
    bool digested;

    //! flag: if true, overall modifications have to be applied
    bool overall_mods;

    //! flag: if true, sample can be registered in database
    bool modified;

    //! contains instance of \c ProteinDigest
    ProteinDigest* protein_digest;

    //! contains all occurring types of overall modifications
    std::vector<Modification*> overall_modifications;

    //! contains the string signifying partial modifications (must be unique for each modification scenario)
    std::string partial_modification_string;

    //! contains all overall modifications in one string (for storing in database)
    std::string overall_modification_string;

    //! contains partial modification (same datastructure as output of \c ModificationParser::parse)
    vector<pair<int, vector<Modification*> > > partial_mods;

    /** contains vectors of Annotations, each vector for one peak.
     *  indexes of this vector of vector of annotations are stored as metaValues in peaks of \c peaklist
     */
    vector<vector<Annotation> > annotation_vectors_;

    //METHODS: --------------------------------------------------------

    //! reads all necessary information out of the inifile and creates needed objects
    bool initialize_();

    //! applies member \c enzyme to member \c protein_digest
    void digest_();

    //! same as \c digest(), but if no enzyme is specified, simply nothing is done
    void tryDigest_() ;

    //! applies overall and partial modifications, specified in .ini -file, to member \c protein_digest
    void modify_();

    //! this method applies, if present, overall modifications to \c protein_digest
    void tryModifyOverall_();

    //! registers this particular instance of Sample in the Database
    void dbRegister_();

    //! returns true, if an protein_modification_scenario like in actual sample already exists. if so, prot_mod_scen_ID is set.
    bool existsProtModScenInDB_();

    //! returns true, if exactly same sample already exists in DB, meaning enzyme and prot_mod_scen. sample_ID and prot_mod_scen_ID are set
    bool existsInDB_();

    //! adds whole protein as a fragment into database (if not already present) and returns its id in table "digest_fragment"
    std::string dbAddWholeProtein_();

    //! reads "real" peaklist out of file specified in \c peaklist
    void readPeaklist_(std::string type, bool verbose=false);

    //! returns fragments (fragment_ID's) out of database-table "digest_fragment". if undigested, it ensures whole protein is in database
    vector<vector<int> > getFragments_(bool& whole_protein);


    /*****************************************************************************************************************************************
     * Method: "enumerate" 
     ****************************************************************************************************************************************/

    /** enumerates all combinations of partial modifcations
     *  each possible position is once considered to be modified, and once to be not modified
     *  the combinations are stored in database table "realized_modification". 
     *  the last entry of each combination points to first ID of overall modifications (obtained within this method by calling 
     *  \c dbStoreOverallModifications()). (they are the same for each combination)
     *  the first ID of each combination is stored in database table modification_combination
     *  (again chained: the field next_modification_combination_ID points to the entry that contains the first ID (in realized_modification)
     *                  of the next combination of actual protein_modification_scenario
     *                  after the last combination of actual scenario this field contains NULL)
     *  return value is the ID of first entry into modification_combination: to be stored in protein_modification_scenario
     ************************************************************************************************************************
     *USAGE: (for applying some overall and some partial modifications)
     *  - modifyOverall    : stores information in \c sequence_overall_modifications
     *  - modifyPartially  : enumerates (call: recursiveEnumerate), checks if position already occupied, saves everything to database
     *                                                                                       (among others: call dbStoreOverallModifications)
     ************************************************************************************************************************
     *IMPORTANT: if each positions to be modified should also be able to be unmodified, a fake modification "Unmodified" has to be used
     */
    int modifyPartiallyEnumerate_(std::vector<std::pair<int, std::vector<Modification*> > > mods, bool verbose=false);

    /** enumerates all partial modification combinations
     *  and stores them in database table "realized_modifications"
     *  the return value is a vector of pairs with the first ID and last ID in realized_midification of each combination, 
     *  first is needed to be stored in table "modification_combination"
     *  second is needed to insert pointer to already existing overall modifications
     *  datastructure "list" is used to allow effective splicing durig recursion
     */
    list<pair<int,int> > recursiveEnumerate_(vector<pair<int, vector<Modification*> > > x, vector<pair<int,int> > accu, bool verbose=false);

    //! calculates all possibly occuring masses and all their possible annotations (\c masstype: "mono" or "average")
    void calculateAnnotations_(std::string masstype = "average", bool verbose = false);

    //! that's what's special about the enumerative method
    void annotateEnumerative_();

    //! searches Peak in Database, and calls printAnnotations, if annotations found
    void annotatePeak_(OpenMS::DPeak<1>& peak);

    //! prints out file with annotations, sign. by their digest_fragment and (first) realized_modification_ID file is stored in \c outputdir
    void storeAnnotations_(vector<vector<int> > annotations, OpenMS::DPeak<1>& peak);


    /*****************************************************************************************************************************************
     * Method: "improved_enumerate" 
     ****************************************************************************************************************************************/

    //VARIABLES: ------------------------------------------------------

    //! same as \c partial_mods , just an "int" vector instead of "Modification*"
    vector<pair<int, vector<int> > > partial_mods_int;

    /** contains information about POSITIONS of partial modifications, used by improved enumerative method
     *  first "int" (map key): modification ID
     *  "int(s)" in vector   : position(s) that can be modified by this modification
     */
    __gnu_cxx::hash_map<int,vector<int> > modification_positions;

    //! same as \c partial_mods_int, only for fragment actually iterated (improved_enumerative method: isFragmentModifyable)
    vector<pair<int, vector<int> > > actual_fragment_partial_mods_int;

    //! contains modification groups, and free positions of actual fragment (improved_enumerative method: isFragmentModifyable)
    map<vector<int>, int> actual_fragment_groups;

    //! hash_map from modification-ID to a vector of group-ID's in which this mod. is (improved_enumerative method: isFragmentModifyable)
    __gnu_cxx::hash_map<int, vector<int> > actual_fragment_mod_with_groups;

    //METHODS: --------------------------------------------------------

    //! implements improved recursive enumeration method: no positions are stored, only numbers of occurrences for each modification
    int modifyPartiallyImprovedEnumerate_(std::vector<std::pair<int, std::vector<Modification*> > > mods, bool verbose=false);

    //! the heart of the improved recursive enumeration method
    list<pair<int,int> > improvedRecursiveEnumerate_(vector<int> current_possible_modification_set, int current_no_free_pos,
        vector<pair<int,int> > accu, vector<vector<int> > further_modification_sets,
        vector<int> further_free_pos);

    //! calculates annotation according to improved enumeration method
    void improvedCalculateAnnotations_(string masstype);

    //! that's what's special about the improved enumerative method
    void annotateEnumerativeImproved_();

    //! searches peaks in database, according to improved enumerative method
    void improvedAnnotatePeak_(OpenMS::DPeak<1>& peak, vector<int> fragment, double fragment_mass, double fragment_unmod_mass,
                               bool whole_protein, __gnu_cxx::hash_map<int,int> ov_mods);

    /** stores information about POSITIONS of partial modifications in \c modification_positions,
     *  to later check, wheter enough mod sites are in given frag for given modification_combination
     *  only sensible to call, if partial_mods are already parsed
     */
    void storePartialModsPosInfo_();

    //! stores information about actual fragment, needed by isFragmentModifyable (improved_enumerative method)
    void storeFragmentInfo_(int start_pos, int end_pos);

    /** returns true, if actual fragment is modifyable by given modification combination
     *  only to be used sesible after call of \c storeFragmentInfo_(int start_pos, int end_pos)
     */
    bool isActualFragmentModifyable_(map<int, int> mod_occurrences);

    /** returns true, if fragment, given by start and end, can possibly be modified by given modification in given multiplicity
     *  only to be called after \c storePartialModsPosInfo() OLD, NOT USED ANY MORE
     */
    bool isFragmentModifyable_(int start_pos, int end_pos, int mod_id, int no_of_occurrences);

    //! returns total number of modification sites contained in given fragment, OLD, NOT USED ANY MORE
    int getTotalNumberOfModSites_(int start_pos, int end_pos);


    /*****************************************************************************************************************************************
     * Method: "peakwise_cormen": modification of T.H. Cormen's "Subset Sum" - Algorithm (Cormen, "Introduction to Algorithms, p. 1045)
     ****************************************************************************************************************************************/

    //VARIABLES: ------------------------------------------------------

    //temporarily saves all combinations of modifications, that are already visited
    __gnu_cxx::hash_map<string, bool> cormen_temp_combinations;

    //indices in vector: groups; values: numbers of positions possibly modified by this group of modifications
    vector<int> cormen_groups_positions;

    /** specific modifications for actual sample:
     *  contains each modification as many times as many positions it can possibly modify at most 
     *    (eg. ID 10 can modify at pos 34 and 52 => mod with ID 10 is present two times in \c modifications
     *  the first "int" in the "pair" signifies the modification_ID, the second "int" the modification group
     */
    vector<pair<int,pair<int,double> > > cormen_modifications;

    //METHODS: --------------------------------------------------------

    //! generates a key-string out of a concrete (modification) combination, for storing in \c cormen_temp_combinations
    string generateKey_(list<int> combination_IDs);

    /** checks wheter given modification combination l does not exceed position limits for each modification group
     *  "list<int>"  : contains ID's of modifications of this combination
     *  "vector<int>": indices: groups; values: how many modifications of corresponding group are already contained in "list<int>"
     */
    bool satisfiesGroupPos_(const pair<list<int>, vector<int> >& l);

    //! adds next Modification to given list L of modification combinations (cf. cormen: "L + x_i")
    list<pair<double, pair<list<int>, vector<int> > > > addModification_(list<pair<double, pair<list<int>, vector<int> > > > L,
        double x_i, int x_i_ID, int x_i_group);

    /** this is the main function of the annotation method "peakwise_cormen"
     *  return value: "list<int>": what modifications do realize (masses sum up to) "double"
     *                "vector<int>": what groups (indices in vector), and how many mod`s per group (values in vector)
     */
    list<pair<double, pair<list<int>, vector<int> > > > exactSubsetSum_(double t, double range, bool verbose = false);

    //! fills \c cormen_modifications and cormen_groups_positions (for fragment signified by \c start_pos and \c end_pos)
    void fillCormenVariables_(int start_pos, int end_pos, Modification& mod_it);

    //! writes annotations into file, in peakwise_cormen - format
    void storeAnnotationsPeakwiseCormen_(list<pair<double, pair<list<int>, vector<int> > > > modification_combinations, OpenMS::DPeak<1>& peak,
                                         vector<int> fragment, double fragment_mass, double fragment_unmod_mass, Modification& mod,
                                         __gnu_cxx::hash_map<int,int> ov_mods);

    //! the peakwise method out of the book 'Cormen'
    void annotatePeakwiseCormen_();

    /****************************************************************************************************************************************/



  public:
    //! constructor without argument: not allowed
    Sample();

    //! constructor with argument: \c sample_data contains all necessary information about actual sample
    Sample(__gnu_cxx::hash_map<std::string, std::string> sample_data,
           std::vector<OpenMS::Spectrum1DWidget::Spectrum1D::iterator>& peaklist,
           std::vector<std::string> ov_mods=std::vector<std::string>(),
           std::string db_username="", std::string db_password="", std::string db_host="",
           Annotate* qannotate=NULL);

    //! destructor
    ~Sample();

    //! copy constructor: explicitly initialize base class, an appropriate assignment operator is generated automatically
    Sample(const Sample& sample);

    /** is the interface to the functionality of this class
     *  calls different other functions, depending on \c annotion_method
     */
    void annotate();

    //prints annotations into files. this method only can be called after call of \c annotate()
    void printAnnotations();

    // stores annotations as metadata in the spectrum from SpecAnnotate from which peaks in peaklist are.
    void storeAnnotations();

    /** exports annotated peaklist.
     *  since it is not yet possible to store pointers to objects as metavalues in OpenMS::DPeaks, the peaklist ist returned, 
     *  containing indices in a vector of vector of Annotations as metavalues, signifying what vector of annotations corresponds to
     *  what peak. 
     *  the vector of vector of Annotations is also returned
     *  IMPORTANT: like \c printAnnotations() this method can only be called, after call of \c annotate()
     */
    pair<OpenMS::DPeakArray<1, OpenMS::DPeak<1> >*, vector<vector<Annotation> > > getAnnotations();

  };

#endif
