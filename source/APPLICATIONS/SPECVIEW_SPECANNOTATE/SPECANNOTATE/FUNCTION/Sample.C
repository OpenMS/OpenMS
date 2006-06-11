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
// $Id: Sample.C,v 1.6 2006/05/30 15:46:40 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer:$
// --------------------------------------------------------------------------

#include "Sample.h"


//STL includes
#include <fstream>
#include <iomanip>
#include <sstream>

//OpenMS includes
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/String.h>

//My includes
#include "ModificationStringParser.h"
#include "string_hash_stl_fixes.h"
#include "../config_specannotate.h"

//conditional QT includes
#ifdef ANNOTATE_QT
#include "Annotate.h"
#include "CustomEvents.h"
#include <qapplication.h>
#include <qsqldatabase.h>
#endif

using namespace std;
using namespace __gnu_cxx;
using namespace OpenMS;


/*------------------------------------------------------------------------------------------------------------------------------------------
 * PUBLIC:
 *----------------------------------------------------------------------------------------------------------------------------------------*/

//! constructor without argument: not allowed
Sample::Sample()
{
  cerr << "Class Sample only can be initialized with sample_data!" << endl;
  throw OpenMS::Exception::Base(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Wrong initialization of class Sample", "No data specified!");
}



//! constructor with argument: \c sample_data contains all necessary information about actual sample
Sample::Sample(hash_map<string, string> sample_data,
               vector<OpenMS::DPeakArray<1, OpenMS::DPeak<1> >::iterator>& peaklist,
               vector<string> ov_mods, string db_username, string db_password, string db_host,
               Annotate* qannotate)
{
  db_username_ = db_username;
  db_password_ = db_password;
  db_host_     = db_host;

  qannotate_ = qannotate;
  overall_mod_strings_ = ov_mods;
  sample_data_ = sample_data;

  //fill peaklist with own objects
  if (!peaklist.empty())
    {
      for (vector<OpenMS::Spectrum1DWidget::Spectrum1D::iterator>::iterator it = peaklist.begin(); it != peaklist.end(); it++)
        {
          DPeak<1> temp_peak;

          temp_peak.getPosition()[0] = (**it).getPosition()[0];
          temp_peak.getIntensity() = (**it).getIntensity();

          peaklist_.push_back(temp_peak);

          //DEBUG
          //cout << "Peaklist (incl. borders of Spectrum): Adding Peak " << (**it).getPosition()[0] << endl;
        }

      //delete first and last element: they are not selected peaks, but borders of spectrum
      peaklist_.erase(peaklist_.begin());
      OpenMS::DPeakArray<1, OpenMS::DPeak<1> >::iterator it = peaklist_.end();
      it--;
      peaklist_.erase(it);

      //store iterators to instances of DPeak<1> of TOPPView
      external_peaklist_ = peaklist;

    }

  initialize_();
}



//! destructor
Sample::~Sample()
{
  delete enzyme;
  delete protein_digest;

  // ********************* again, this seems to be reasonable, but produces a seqfault ***
  // delete sql_adapter_;
  // *************************************************************************************

  for (vector<Modification*>::iterator it = overall_modifications.begin(); it != overall_modifications.end(); it++)
    {
      delete (*it);
    }
}



//! copy constructor
Sample::Sample(const Sample& sample)
{
  db_username_ = sample.db_username_;
  db_password_ = sample.db_password_;
  db_host_     = sample.db_host_;

  qannotate_ = sample.qannotate_;
  overall_mod_strings_ = sample.overall_mod_strings_;
  sample_data_ = sample.sample_data_;
  peaklist_ = sample.peaklist_;

  initialize_();
}



void Sample::annotate()
{
  //remove annotations of possible previous call of annotate()
  annotation_vectors_.clear();

#ifndef ANNOTATE_XML
  //if enzyme is specified in .ini-file: digest
  if (enzyme != NULL)
    {
      // Some output
#ifdef ANNOTATE_QT
      oe_ = new OutputEvent("Sample::annotate(): Digesting Sample...\n"); // Qt will delete it when done!
      QApplication::postEvent(qannotate_, oe_);
#endif
      #ifndef ANNOTATE_QT

      cout << "Sample::annotate(): Digesting Sample..." << endl;
#endif

      digest_();
    }
#endif

  //here is the point for deciding whether we have a 'peakwise' annotation method or not
  if (annotation_method.find("peakwise") != string::npos)
    {
      //read "real" peaklist in verbose mode
      readPeaklist_(peakfile_format, true);

      if (annotation_method == "peakwise_cormen")
        {
          annotatePeakwiseCormen_();
        }
    }
#ifndef ANNOTATE_XML
  else
    {

      //if sample is not already modified, do it (modify checks for itself, which \c annotation_method to use)
      if (!modified)
        {
          // Some output
#ifdef ANNOTATE_QT
          // Qt will delet it when done!
          oe_ = new OutputEvent("Sample::annotate(): Modifying Sample ("+ annotation_method +")...\n");
          QApplication::postEvent(qannotate_, oe_);
#endif
	  #ifndef ANNOTATE_QT

          cout << "Sample::annotate(): Modifying Sample ("+ annotation_method +")..." << endl;
#endif

          modify_();                      //this method is intelligently crafted, it only calculates new mod's if not already present in DB
        }

      //different behaviour for different \c annotation_methods
      if (annotation_method == "enumerate")
        {
          annotateEnumerative_();
        }
      else if (annotation_method == "improved_enumerate")
        {
          annotateEnumerativeImproved_();
        }
    }
#endif


#ifdef ANNOTATE_QT
  // Qt will delete it when done!
  oe_ = new OutputEvent("Done.\n");
  QApplication::postEvent(qannotate_, oe_);
#endif
  #ifndef ANNOTATE_QT

  cout << "Done." << endl;
#endif
}



void Sample::printAnnotations()
{
#ifdef ANNOTATE_QT
  // Qt will delete it when done!
  oe_ = new OutputEvent("Printing annotations into files:\n");
  QApplication::postEvent(qannotate_, oe_);
#endif
#ifndef ANNOTATE_QT

  cout << "Printing annotations into files:" << endl;
#endif

  if(outputdir != "")
    {
      for (DPeakArray<1, DPeak<1> >::iterator it = peaklist_.begin(); it != peaklist_.end(); it++)
        {
          //are there annotations for this peak? if yes, get them, if no, go to next peak
          int index;
          try
            {
              index = (int)(it->getMetaValue("annotations"));
            }
          catch (Exception::InvalidValue)
            {
              continue;
            }
          catch (OpenMS::Exception::ConversionError)
            {
              continue;
            }

          //generate proper filename
          ostringstream ost_file;
          ost_file << "peak_" << setw(9) << setfill('0') << setprecision(2) << fixed << it->getPosition()
            [0];
          if (annotation_method == "enumerate")
            {
              ost_file << ".enum_annot";
            }
          else if (annotation_method == "improved_enumerate")
            {
              ost_file << ".improved_enum_annot";
            }
          else if (annotation_method == "peakwise_cormen")
            {
              ost_file << ".peakw_cormen_annot";
            }
          string filename = ost_file.str();

          //open file
          ofstream ofst((outputdir + filename).c_str());
          if (!ofst)
            {
              throw OpenMS::Exception::Base(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Wrong Output Directory", "Could not create File!");
            }

          ofst << "#########################################################################################################" << endl;
          ofst << "ANNOTATIONS found for PEAK at a m/z value of " << setprecision(6) << fixed << it->getPosition()
            [0];
          ofst << ", within a search range of " << setprecision(6) << fixed << range;
        ofst << " Daltons." << endl;
        ofst << "#########################################################################################################" << endl;
        ofst << endl << endl;

        //write annotations into file
        int annotation_count = 1;
        for (vector<Annotation>::iterator annit = annotation_vectors_[index].begin(); annit != annotation_vectors_[index].end(); annit++)
            {
              annit->print(annotation_count, ofst);
              annotation_count++;
            }

          ofst.close();

#ifdef ANNOTATE_QT
          // Qt will delete it when done!
          oe_ = new OutputEvent("File " + outputdir + filename + " created.\n");
          QApplication::postEvent(qannotate_, oe_);
#endif
          #ifndef ANNOTATE_QT

          cout << "File " << + outputdir + filename + " created." << endl;
#endif

        }
    }
#ifdef ANNOTATE_QT
  // Qt will delete it when done!
  oe_ = new OutputEvent("Done.\n");
  QApplication::postEvent(qannotate_, oe_);
#endif
#ifndef ANNOTATE_QT

  cout << "Done." << endl;
#endif
}



void Sample::storeAnnotations()
{
#ifdef ANNOTATE_QT
  // Qt will delete it when done!
  oe_ = new OutputEvent("Storing annotations as meta data into spectrum :\n");
  QApplication::postEvent(qannotate_, oe_);
#endif
#ifndef ANNOTATE_QT

  cout << "Storing annotations as meta data into spectrum:" << endl;
#endif

  //register name "extended_label" for storing information to peaks in TOPPView in the MetaInfoRegistry
  unsigned int ex_lab_id = (*(external_peaklist_.begin()))->metaRegistry().registerName((string)"extended_label", (string)"annotations");

  if(outputdir == "")
    {
      //imcrement external peaklist in the same way as internal
      vector<OpenMS::DPeakArray<1, OpenMS::DPeak<1> >::iterator>::iterator ext_it = external_peaklist_.begin();
      ext_it++; //first element in external peaklist is just the border of the spectrum, no real peak;

      for (DPeakArray<1, DPeak<1> >::iterator it = peaklist_.begin(); it != peaklist_.end(); it++)
        {
          //are there annotations for this peak? if yes, get them, if no, go to next peak
          int index;
          try
            {
              index = (int)(it->getMetaValue("annotations"));
            }
          catch (Exception::InvalidValue)
            {
              ext_it++;
              continue;
            }
          catch (OpenMS::Exception::ConversionError)
            {
              ext_it++;
              continue;
            }

          //open stringstream for storing annotations as metadata with peak in spectrum in TOPPView
          ostringstream annot_ost;

          annot_ost << "#########################################################################################################" << endl;
          annot_ost << "ANNOTATIONS found for PEAK at a m/z value of " << setprecision(6) << fixed << it->getPosition()
            [0];
          annot_ost << ", within a search range of " << setprecision(6) << fixed << range;
        annot_ost << " Daltons." << endl;
        annot_ost << "#########################################################################################################" << endl;
        annot_ost << endl << endl;

        //store annotations as meta value of type string in peak referenced by corresponding iterator in external_peaklist_:
        //collect all annotations in one stringstream
        int annotation_count = 1;
        for (vector<Annotation>::iterator annit = annotation_vectors_[index].begin(); annit != annotation_vectors_[index].end(); annit++)
            {
              annit->print(annotation_count, annot_ost);
              annotation_count++;
            }

          //store the stringstream as string as meta data in spectrum
          OpenMS::DPeakArray<1, OpenMS::DPeak<1> >::iterator tmp = (*ext_it);
          tmp->setMetaValue(ex_lab_id, annot_ost.str());
          tmp->setMetaValue(3,((String)(annotation_count - 1) + " " + annotation_method + " annotations present"));

          //increment ext_it to keep track with it
          ext_it++;

#ifdef ANNOTATE_QT
          // Qt will delete it when done!
          oe_ = new OutputEvent("Annotations stored for peak at " + String(it->getPosition()[0]) + " m/z.\n");
          QApplication::postEvent(qannotate_, oe_);
#endif
          #ifndef ANNOTATE_QT

          cout << "Annotations stored for peak at " + String(it->getPosition()[0]) + " m/z.\n" << endl;
#endif

        }
    }
#ifdef ANNOTATE_QT
  // Qt will delete it when done!
  oe_ = new OutputEvent("Done.\n");
  QApplication::postEvent(qannotate_, oe_);
#endif
#ifndef ANNOTATE_QT

  cout << "Done." << endl;
#endif
}



pair<DPeakArray<1, DPeak<1> >*, vector<vector<Annotation> > > Sample::getAnnotations()
{
  /* this is a hack: peaklist_ is copied by hand into an object for which this class has no responsibility!!
     this is done, because automatic COPYING (call by value) OF DPeakArray<1, DPeak<1> > PRODUCES SEGFAULT!!!
     annotation_vector_ can be copied all right!
   */
  DPeakArray<1, DPeak<1> >* output = new DPeakArray<1, DPeak<1> >();
  for (DPeakArray<1, DPeak<1> >::iterator it = peaklist_.begin(); it != peaklist_.end(); it++)
    {
      DPeak<1> temp_peak;

      temp_peak.getPosition()[0] = (*it).getPosition()[0];
      temp_peak.getIntensity() = (*it).getIntensity();

      //register "annotations" as meta-value for peak, if already registered, nothing happens
      temp_peak.metaRegistry().registerName("annotations", "annotations found for this peak");

      try
        {
          temp_peak.setMetaValue("annotations", (int)((*it).getMetaValue("annotations")));
        }
      catch (OpenMS::Exception::ConversionError) //if no annotation is stored for (*it) simply do nothing
        {
        }

      output->push_back(temp_peak);

      //DEBUG
      cout << "Exporting peak at position: " << (*it).getPosition()[0] << endl;
    }

  return  pair<DPeakArray<1, DPeak<1> >*, vector<vector<Annotation> > >(output, annotation_vectors_);
}



/*------------------------------------------------------------------------------------------------------------------------------------------
 * PRIVATE:
 *----------------------------------------------------------------------------------------------------------------------------------------*/



// ****************************************************************************************************************************************
// * Used by more than one method 
// ***************************************************************************************************************************************


//! reads all necessary information out of the inifiles and creates needed objects of derived classes of \c Modification, as well as initial vector of \c Compound
bool Sample::initialize_()
{
  // \c enzyme not yet applied
  digested = false;

  // instance of \c Sample not yet modified
  modified = false;

  // \c sample_ID not yet known
  sample_ID = -1;

#ifndef ANNOTATE_XML
  // creating MySQL Adapter and connecting to database
  sql_adapter_ = new MySQLAdapter();
  sql_adapter_->connect(db_username_.c_str(), db_password_.c_str(), db_host_.c_str());
  sql_adapter_->selectDB(DATABASE);
#endif

  //! for printing out double values throughout the whole programm: set a fixed accuracy of 6 digits after the decimal point (default)
  cout.setf(ios_base::fixed,ios_base::floatfield)
    ;

  // set annotation_method
  annotation_method = sample_data_["annotation_method"];

  //! set filename of real peaklist
  peakfile = sample_data_["peakfile"];

  //! set format of real peaklist
  peakfile_format = sample_data_["peakfile_format"];

  //! set output directory
  outputdir = sample_data_["outputdir"];

  //! set masstype
  masstype = sample_data_["masstype"];

  //! set parameter for search of real masses in calculated masses
  range = atof((sample_data_["search_range"]).c_str());

  //! set partial_modification_string
  partial_modification_string = sample_data_["partial_modification_string"];

  //! instanciating and initializing member variables: \c Modification used as overall modifications
  if (overall_mod_strings_.size() != 0) //overall modifications present?
    {
      for (vector<string>::iterator it = overall_mod_strings_.begin(); it != overall_mod_strings_.end(); it++)
        {
          Modification* new_mod = new Modification(*it, db_username_, db_password_, db_host_ );
          overall_modifications.push_back(new_mod);

          //store ID's of modifications in one string for storing in database
          overall_modification_string += (String(new_mod->getID()) + ", ");

        }
      overall_mods = true;
    }
  else
    {
      overall_modification_string = "void";
      overall_mods = false;
    }

  //! instanciating and initializing member variables: \c Enzyme
  string enz = sample_data_["enzyme"];
  if (enz != "")
    {
      enzyme = new Enzyme(enz, db_username_, db_password_, db_host_);
    }
  else
    {
      enzyme = NULL;
    }

  //! instanciating and initializing member variables: \c ProteinDigest
  protein_digest = new ProteinDigest(sample_data_["protein"], 0, db_username_, db_password_, db_host_);

  //! up to now the boolean return-value of this procedure is not used
  return true;
}



//! applies member \c enzyme to member \c protein_digest
void Sample::digest_()
{
  if (enzyme != NULL)
    {
      protein_digest->digest(enzyme);
      digested = true;
    }
  else
    {
      throw OpenMS::Exception::Base(__FILE__, __LINE__, __PRETTY_FUNCTION__, "No Enzyme Specified", "Could not execute Sample::digest_(), no Enzyme!");
    }
}



//! same as \c digest_(), but if no enzyme is specified, simply nothing is done
void Sample::tryDigest_()
{
  if (enzyme != NULL)
    {
      protein_digest->digest(enzyme);
      digested = true;
    }
}



//! applies overall and partial modifications, specified in .ini -file, to member \c protein_digest
void Sample::modify_()
{
  //for this annotation method overall modifications have to be applied even if protein_modification_scenario already exists in DB
  if (annotation_method == "improved_enumerate")
    {
      // apply overall modifications, if present:
      tryModifyOverall_();
    }

  if (!existsProtModScenInDB_())
    {
      //parse partial modifications:
      ModificationStringParser mod_parser_(db_username_, db_password_, db_host_);
      partial_mods = mod_parser_.parse(partial_modification_string);

      //! different behaviour for different annotation methods
      if (annotation_method == "enumerate")
        {
          // apply overall modifications, if present:
          if (overall_mods)
            {
              for (vector<Modification*>::iterator it = overall_modifications.begin(); it != overall_modifications.end(); it++)
                {
                  protein_digest->modifyOverall(*it);
                }
            }

          int first_mod_comb_ID = modifyPartiallyEnumerate_(partial_mods);

          //fill database table PROT_MOD_SCEN_TABLE
          sql_adapter_->executeQuery("INSERT INTO " + (string)PROT_MOD_SCEN_TABLE + " ( `protein_ID` , `overall_modifications`"
                                     +", `annotation_method` , `partial_modifications` , `modification_combination_ID` ) "
                                     + " VALUES ( '" + String(protein_digest->getProteinID()) + "', '" + overall_modification_string
                                     + "', '" + annotation_method + "', '" + partial_modification_string + "', '"
                                     + String(first_mod_comb_ID) +"' )");
        }
      else if (annotation_method == "improved_enumerate")
        {
          int first_mod_comb_posless_ID = modifyPartiallyImprovedEnumerate_(partial_mods);

          //fill database table PROT_MOD_SCEN_TABLE
          sql_adapter_->executeQuery("INSERT INTO " + (string)PROT_MOD_SCEN_TABLE + " ( `protein_ID` , `overall_modifications`"
                                     +", `annotation_method` , `partial_modifications` , `modification_combination_positionless_ID` ) "
                                     + " VALUES ( '" + String(protein_digest->getProteinID()) + "', '" + overall_modification_string
                                     + "', '" + annotation_method + "', '" + partial_modification_string + "', '"
                                     + String(first_mod_comb_posless_ID) +"' )");
        }



      //store protein_modification_scenario_ID
      sql_adapter_->executeQuery("SELECT last_insert_id() FROM " + (string)PROT_MOD_SCEN_TABLE + " LIMIT 1");
      prot_mod_scen_ID = atoi(sql_adapter_->getUnaryResult().c_str());

      //clean up instances of partial modifications
      for (vector<pair<int, vector<Modification*> > >::iterator it = partial_mods.begin(); it != partial_mods.end(); it++)
        {
          for (vector<Modification*>::iterator itt = it->
               second.begin();
               itt != it->second.end();
               itt++)
            {
              delete *itt;
            }
        }
    }
  modified = true;
}



//! this method applies, if present, overall modifications to \c protein_digest
void Sample::tryModifyOverall_()
{
  if (overall_mods)
    {
      for (vector<Modification*>::iterator it = overall_modifications.begin(); it != overall_modifications.end(); it++)
        {
          protein_digest->modifyOverall(*it);
        }
    }
}



void Sample::dbRegister_()
{
  if (!existsInDB_())
    {
      if (modified)
        {
          if ((enzyme != NULL) && (digested))
            {
              //add to database, SAMPLE_TABLE, including enzyme
              sql_adapter_->executeQuery("INSERT INTO " + (string)SAMPLE_TABLE + " ( `enzyme_ID` , `protein_modification_scenario_ID`"
                                         + " , `annotation_method`) "
                                         + " VALUES ( '" + String(enzyme->getID()) + "', '" + String(prot_mod_scen_ID)
                                         + "', '" + annotation_method + "' )");
            }
          else
            {
              //add to database, SAMPLE_TABLE, without enzyme (no enzyme is signified by enzyme_ID = -1)
              sql_adapter_->executeQuery("INSERT INTO " + (string)SAMPLE_TABLE + " ( `enzyme_ID` , `protein_modification_scenario_ID`"
                                         + " , `annotation_method`) "
                                         + " VALUES ( '-1', '" + String(prot_mod_scen_ID)
                                         + "', '" + annotation_method +"' )");
            }

          //store sample_ID
          sql_adapter_->executeQuery("SELECT last_insert_id() FROM " + (string)SAMPLE_TABLE + " LIMIT 1");
          sample_ID = atoi(sql_adapter_->getUnaryResult().c_str());
        }
      else
        {
          throw OpenMS::Exception::Base(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Not modified",
                                        "Could not register to database: Sample::modify_() has not been executed!");
        }
    }
}



bool Sample::existsProtModScenInDB_()
{
  //! also compares \c annotation_method!!!
  sql_adapter_->executeQuery("SELECT protein_modification_scenario_ID FROM " + (string)PROT_MOD_SCEN_TABLE
                             + " WHERE protein_ID = \"" + String(protein_digest->getProteinID()) +
                             "\" AND overall_modifications = \"" + overall_modification_string +
                             "\" AND annotation_method = \"" + annotation_method +
                             "\" AND partial_modifications = \"" + partial_modification_string + "\"");

  string result_string;
  if (sql_adapter_->ifGetUnaryResult(result_string))
    {
      prot_mod_scen_ID = atoi(result_string.c_str());
      return true;
    }
  else
    {
      return false;
    }
}



//! returns true, if exactly same sample already exists in DB, meaning enzyme and prot_mod_scen. sample_ID and prot_mod_scen_ID are set
bool Sample::existsInDB_()
{
  if (existsProtModScenInDB_())
    {
      string enzyme_ID;
      if ((enzyme == NULL) || (digested == false))
        {
          enzyme_ID = "-1";
        }
      else if ((enzyme != NULL) && (digested == true))
        {
          enzyme_ID = String(enzyme->getID());
        }

      sql_adapter_->executeQuery("SELECT sample_ID FROM " + (string)SAMPLE_TABLE
                                 + " WHERE enzyme_ID = \"" + enzyme_ID +
                                 "\" AND protein_modification_scenario_ID = \"" + String(prot_mod_scen_ID) +
                                 "\" AND annotation_method = \"" + annotation_method + "\"");
      string result_string;
      if (sql_adapter_->ifGetUnaryResult(result_string))
        {
          sample_ID = atoi(result_string.c_str());
          return true;
        }
      else
        {
          return false;
        }
    }
  else
    {
      return false;
    }
}



string Sample::dbAddWholeProtein_()
{
  //get ID (in digest_fragment) of protein fragment signifying whole protein
  string frag_id;
  sql_adapter_->executeQuery("SELECT digest_fragment_ID FROM " + (string)FRAGMENT_TABLE
                             + " WHERE protein_ID = \"" + String(protein_digest->getProteinID())
                             + "\" AND d_start_pos = \"0"
                             + "\" AND enzyme_ID = \"-1"  //no enzyme: enzyme ID = -1!!
                             + "\" AND d_end_pos = \"" + String(protein_digest->getProteinLength()-1) + "\" LIMIT 0,1");
  //above, LIMIT 0,1 is set out of following reason: if first a whole protein is annotated, and no digest is
  //present in digest_fragments. then whole protein is added into it, with enzyme ID of NULL.
  //if after that, this protein is annotated with digest, whole protein is added again, but with enzyme ID
  //if after that, this protein is annotated with digest by another enzyme, whole protein is added even again
  //if after that, protein is again annotated without protein, above query would not be unique without the LIMIT

  //if such fragment does not exist yet, add it
  if (!sql_adapter_->ifGetUnaryResult(frag_id))
    {
      //add whole protein as fragment to database
      sql_adapter_->executeQuery("INSERT INTO " + (string)FRAGMENT_TABLE
                                 + " ( `protein_ID` , `enzyme_ID`,  `d_start_pos` , `d_end_pos` ) "
                                 + " VALUES ( '" + String(protein_digest->getProteinID()) + "', '-1' , '0' , '"
                                 + String(protein_digest->getProteinLength()-1) + "' )");

      //get ID of this new entry
      sql_adapter_->executeQuery("SELECT last_insert_id() FROM " + (string)FRAGMENT_TABLE + " LIMIT 1");
      frag_id = sql_adapter_->getUnaryResult();
    }
  return frag_id;
}



//! reads "real" peaklist out of file specified in \c peaklist
void Sample::readPeaklist_(string type, bool verbose)
{
  //check if peaklist was already set in constructor
  if (peaklist_.empty())
    {
      ifstream infile(peakfile.c_str());
      if (!infile)
        {
          throw OpenMS::Exception::Base(__FILE__, __LINE__, __PRETTY_FUNCTION__, "No Peakfile", ("No Valid Peakfile Specified by " + peakfile + "!").c_str());
        }

      string line;

      //is it a peakfile in michael kerber`s format?
      if (type == "kerber")
        {
          for (getline(infile,line);!infile.eof();getline(infile,line))
            {
              istringstream ist(line);
              string type;
              double mass, height, rel_height, left_width, right_width;

              ist >> type;
              ist >> mass;
              ist >> height;
              ist >> rel_height;
              ist >> left_width;
              ist >> right_width;

              //add new peak
              DPeak<1> new_peak;
              new_peak.getPosition()[0] = mass;
              new_peak.getIntensity() = height;

              // 	  new_peak.setMetaValue("RelativeHeight", rel_height);
              // 	  new_peak.setMetaValue("LeftWidth", left_width);
              // 	  new_peak.setMetaValue("RightWidth", right_width);
              // 	  new_peak.setMetaValue("PeakType", type);

              peaklist_.push_back(new_peak);

            }
        }
      //is it a peakfile in hansjoerg toll`s format?
      else if (type == "toll")
        {
          //throw the first line away
          getline(infile,line);

          //iterate rest of lines
          for (getline(infile,line);!infile.eof();getline(infile,line))
            {
              istringstream ist(line);

              double mass, height;

              //first element in each row is number of peak: throw away
              string trash;
              ist >> trash;

              ist >> mass;
              ist >> height;


              //add new peak
              DPeak<1> new_peak;
              new_peak.getPosition()[0] = mass;
              new_peak.getIntensity() = height;

              // 	  new_peak.setMetaValue("RelativeHeight", 0);
              // 	  new_peak.setMetaValue("LeftWidth", 0);
              // 	  new_peak.setMetaValue("RightWidth", 0);
              // 	  new_peak.setMetaValue("PeakType", "");

              peaklist_.push_back(new_peak);
            }
        }
      else
        {
          throw OpenMS::Exception::Base(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Wrong Peaklist Format", "Type of Peaklist Format not  known!");
        }

      infile.close();
    }

  if (verbose)
    {
      for (DPeakArray<1, DPeak<1> >::iterator it = peaklist_.begin(); it != peaklist_.end(); it++)
        {
#ifdef ANNOTATE_QT
          // Qt will delete it when done!
          oe_
          = new OutputEvent("Using peak at mass: " + String(it->getPosition()[0]) + ", height: " + String(it->getIntensity()) + "\n");
          QApplication::postEvent(qannotate_, oe_);
#endif
          #ifndef ANNOTATE_QT

          cout << "Using peak at mass: " << it->getPosition()[0] << ", height: " << it->getIntensity() << endl;
#endif

        }
    }
}


vector<vector<int> > Sample::getFragments_(bool& whole_protein)
{
  vector<vector<int> > fragments;

#ifdef ANNOTATE_XML

  vector<int> tmp_frag;
  tmp_frag.push_back(-1);                                       //digest_fragment_ID( not relevant in no database-mode)
  tmp_frag.push_back(0);                                        //start_pos
  tmp_frag.push_back((protein_digest->getProteinLength()-1));   //end_pos
  tmp_frag.push_back(-1);                                       //protein_ID ( not relevant in no database-mode )
  tmp_frag.push_back(-1);                                       //enzyme_ID

  fragments.push_back(tmp_frag);
  return fragments;

#endif


  //if no enzyme is used, whole protein has to be added into database, if already present, dbAddWholeProtein returns its id
  if (!digested)
    {
      whole_protein = true;

      string frag_id = dbAddWholeProtein_();

      vector<int> tmp_frag;
      tmp_frag.push_back(atoi(frag_id.c_str()));                    //digest_fragment_ID
      tmp_frag.push_back(0);                                        //start_pos
      tmp_frag.push_back((protein_digest->getProteinLength()-1));   //end_pos
      tmp_frag.push_back(protein_digest->getProteinID());           //protein_ID
      tmp_frag.push_back(-1);                                       //enzyme_ID

      fragments.push_back(tmp_frag);

    }
  //enzyme used: get fragments out of database
  else
    {
      whole_protein = false;

      //get protein_ID and Enzyme_ID
      int protein_ID = protein_digest->getProteinID();
      int enzyme_ID = enzyme->getID();

      //get fragments
      sql_adapter_->executeQuery("SELECT digest_fragment_ID, d_start_pos, d_end_pos FROM " + (string)FRAGMENT_TABLE
                                 + " WHERE enzyme_ID = \"" + String(enzyme_ID) +
                                 "\" AND protein_ID = \"" + String(protein_ID) + "\"");

      //acces query result directly
      QSqlDatabase* db_handle_ = QSqlDatabase::database("db_handle_");
      QSqlQuery res(db_handle_);
      res = sql_adapter_->lastResult();

      //check if result is valid
      if ((!res.isActive()) || (res.size()==0))
        {
          cerr << "Throwing exception because of query-string: " << endl << sql_adapter_->lastQuery() << endl;
          throw OpenMS::Exception::Base(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Wrong Query for Fragments", "Could not obtain fragments!");
        }

      //fill fragments into datastructure created above
      while(res.next())
        {
          vector<int> tmp_frag;
          tmp_frag.push_back(res.value(0).toInt());  //digest_fragment_ID
          tmp_frag.push_back(res.value(1).toInt());  //start_pos
          tmp_frag.push_back(res.value(2).toInt());  //end_pos
          tmp_frag.push_back(protein_ID);
          tmp_frag.push_back(enzyme_ID);

          fragments.push_back(tmp_frag);
        }
    }
  return fragments;
}




/*****************************************************************************************************************************************
 * Method: "enumerate" 
 ****************************************************************************************************************************************/


int Sample::modifyPartiallyEnumerate_(vector<pair<int, vector<Modification*> > > mods, bool verbose)
{
  vector<pair<int,int> > accu;
  list<pair<int,int> > modification_combinations = recursiveEnumerate_(mods, accu, verbose);

  //since \c sequence_modification could have been modified in \c recursiveEnumerate_() it ist not stored to database before
  int first_overall_mod_ID = protein_digest->dbStoreOverallModifications();

  //add first ID of overall modifications to last entry of each combination, only if an overall modification is present
  if (first_overall_mod_ID != -1)
    {
      for (list<pair<int,int> >::iterator it = modification_combinations.begin(); it != modification_combinations.end(); it++)
        {
          sql_adapter_->executeQuery("UPDATE " + (string)REALIZED_MOD_TABLE
                                     + " SET `next_realized_modification_ID` = " + String(first_overall_mod_ID)
                                     + " WHERE `realized_modification_ID` = " + String(it->second) + " LIMIT 1");

        }
    }

  //fill first ID of each combination into table "modification_combination"
  string last_mod_comb_ID, actual_mod_comb_ID, first_mod_comb_ID;
  for (list<pair<int,int> >::iterator it = modification_combinations.begin(); it != modification_combinations.end(); it++)
    {
      // add actual combination
      sql_adapter_->executeQuery("INSERT INTO " + (string)MOD_COMB_TABLE
                                 + " ( `first_realized_modification_ID`) "
                                 + " VALUES ( '" + String(it->first) + "' )");

      //get database id of actual combination
      sql_adapter_->executeQuery("SELECT last_insert_id() FROM " + (string)MOD_COMB_TABLE + " LIMIT 1");
      string actual_mod_comb_ID = sql_adapter_->getUnaryResult();


      // if we do not work on the first combination: add pointer to actual combination into last combination
      if (it != modification_combinations.begin())
        {
          sql_adapter_->executeQuery("UPDATE " + (string)MOD_COMB_TABLE
                                     + " SET `next_modification_combination_ID` = " + actual_mod_comb_ID
                                     + " WHERE `modification_combination_ID` = " + last_mod_comb_ID + " LIMIT 1");
        }

      //store this position
      last_mod_comb_ID = actual_mod_comb_ID;

      //if this is first combination: store in variable.
      if (it == modification_combinations.begin())
        {
          first_mod_comb_ID = actual_mod_comb_ID;
        }
    }
  return atoi(first_mod_comb_ID.c_str());
}



list<pair<int,int> > Sample::recursiveEnumerate_(vector<pair<int, vector<Modification*> > > x, vector<pair<int,int> > accu,
    bool verbose)
{
  list<pair<int,int> > result;
  if (x.size()==0)
    {
      //add accu to database (only if modified (mod_ID != 0)!)
      string actual_mod_ID, last_mod_ID = "-1", first_mod_ID = "-1";
      for (vector<pair<int,int> >::iterator it = accu.begin(); it != accu.end(); it++)
        {
          // add actual position/modification
          sql_adapter_->executeQuery("INSERT INTO " + (string)REALIZED_MOD_TABLE
                                     + " ( `m_position` , `modification_ID` ) "
                                     + " VALUES ( '" + String(it->first) + "', '" + String(it->second) + "' )");

          //get database id of actual modification
          sql_adapter_->executeQuery("SELECT last_insert_id() FROM " + (string)REALIZED_MOD_TABLE + " LIMIT 1");
          actual_mod_ID = sql_adapter_->getUnaryResult();


          // if we do not work on the first modification: add pointer to actual modification into last modification
          if (last_mod_ID != "-1")
            {
              sql_adapter_->executeQuery("UPDATE " + (string)REALIZED_MOD_TABLE
                                         + " SET `next_realized_modification_ID` = " + actual_mod_ID
                                         + " WHERE `realized_modification_ID` = " + last_mod_ID + " LIMIT 1");
            }


          //store this position
          last_mod_ID = actual_mod_ID;

          //if this is first modification: store in variable.
          if (first_mod_ID == "-1")
            {
              first_mod_ID = actual_mod_ID;
            }

          //some output
          if (verbose)
            {
              cerr << "(" << it->first << ", " << it->second << "), ";
            }
        }

      //output ctd
      if (verbose)
        {
          cerr << endl;
        }

      //put first_ID and last ID in realized_modification into result (first: for modification_combination, last: for adding overall mods)
      result.push_back(pair<int,int>(atoi(first_mod_ID.c_str()),atoi(actual_mod_ID.c_str())));
    }
  else
    {
      //enumerate different Modifications of first element of x
      for (vector<Modification*>::iterator it = x[0].second.begin(); it != x[0].second.end(); it++)
        {

          //test: is residue already modified by overall modification
          if (protein_digest->sequence_overall_modifications[x[0].first] != 0)
            {
              cout << "Site to be modified by partial modification with ID " << (*it)->getID();
              cout <<" is also to be modified by overall modification with ID ";
              cout << protein_digest->sequence_overall_modifications[x[0].first] << "." << endl;
              cout << "Using partial modification." << endl;
              protein_digest->sequence_overall_modifications[x[0].first] = 0;
            }

          //test: ist actual modification able to modify actual position?
          if (!((*it)->canModify(String(protein_digest->sequence_oneletter[x[0].first]))))
            {
              throw ProteinDigest::WrongModification(__FILE__, __LINE__, __PRETTY_FUNCTION__, (*it)->getID(), x[0].first);
            }

          //recursion with actual position and modification in accumulator-argument
          vector<pair<int,int> > temp_accu = accu;
          temp_accu.push_back(pair<int,int>(x[0].first, (*it)->getID()));

          //recursion without first element of x
          vector<pair<int, vector<Modification*> > > temp_x = x;
          temp_x.erase(temp_x.begin());                             //first element in copied list should be same as in initial list

          //construct result and recurse
          if (it == x[0].second.begin())
            {
              result = recursiveEnumerate_(temp_x,temp_accu, verbose);
            }
          else
            {
              list<pair<int,int> > recursive_result = recursiveEnumerate_(temp_x,temp_accu, verbose);
              result.splice(result.begin(), recursive_result);
            }
        }
    }
  return result;
}



void Sample::calculateAnnotations_(string masstype, bool verbose)
{
  //is sample already registered in database? if not: register
  if (sample_ID == -1)
    {
      dbRegister_();
    }

  //datastructure for fragments and their masses. both structures correspond via their indices
  vector<vector<int> > fragments;                   //first element in vector: digest_fragment_ID, second: d_start_pos, third: d_end_pos
  vector<double>       fragment_unmodified_masses;

  //datastrucuture for transiently saving a particular modification combination
  vector<pair<int,int> > mod_comb;

  //instance of \c Modification for access to different Masses of Modifications (ID 1, just dummy, will be changed later in this method)
  Modification mod(1, db_username_, db_password_, db_host_);

  //get fragments from database
  bool whole_protein;
  fragments = getFragments_(whole_protein);

  //calculate masses of unmodified fragments
  double mass;
  for (vector<vector<int> >::iterator it = fragments.begin(); it != fragments.end(); it++)
    {
      if (masstype == "mono")
        {
          mass = protein_digest->getFragmentMonoMass((*it)[1],(*it)[2]);
        }
      else
        {
          mass = protein_digest->getFragmentAverageMass((*it)[1], (*it)[2]);
        }
      fragment_unmodified_masses.push_back(mass);
    }

  //some output
  if (verbose)
    {
      cout << "Unmodified masses of fragments:" << endl;
      int abce = 0;
      for (vector<vector<int> >::iterator it = fragments.begin(); it != fragments.end(); it++)
        {
          cout << "start: " << (*it)[1] << " " << protein_digest->getResName((*it)[1]) << ", ";
          cout << "end: "   << (*it)[2] << " " << protein_digest->getResName((*it)[2]) << ", ";
          cout << "unmodified mass: " << fragment_unmodified_masses[abce] << endl;
          abce++;
        }
    }

  //flag for each fragment: true, if unmodified fragment was already added to database!!
  vector<bool> added_unmodified(fragments.size(), false);

  //create variables for modification-iteration
  string mod_comb_ID, real_mod_ID, m_position, modification_ID;

  //iterate modification_combinations
  sql_adapter_->executeQuery("SELECT modification_combination_ID FROM " + (string)PROT_MOD_SCEN_TABLE
                             + " WHERE protein_modification_scenario_ID = \"" + String(prot_mod_scen_ID) + "\"");
  mod_comb_ID = sql_adapter_->getUnaryResult();
  for( ; mod_comb_ID != "0"; )
    {
      //some output
      if (verbose)
        {
          cout << endl << "Modification Combination: " << endl;
        }

      //store always only one particular modification_combination
      mod_comb.clear();

      //iterate realized_modifications
      sql_adapter_->executeQuery("SELECT first_realized_modification_ID FROM " + (string)MOD_COMB_TABLE
                                 + " WHERE modification_combination_ID = \"" + mod_comb_ID + "\"");
      real_mod_ID = sql_adapter_->getUnaryResult();
      string first_realized_mod_ID = real_mod_ID;  //stores first_real_mod_ID of this mod_combination for adding to "annotation"
      for( ; real_mod_ID != "0"; )
        {
          //get position that is modified
          sql_adapter_->executeQuery("SELECT m_position FROM " + (string)REALIZED_MOD_TABLE
                                     + " WHERE realized_modification_ID = \"" + real_mod_ID + "\"");
          m_position = sql_adapter_->getUnaryResult();

          //get ID of Modification
          sql_adapter_->executeQuery("SELECT modification_ID FROM " + (string)REALIZED_MOD_TABLE
                                     + " WHERE realized_modification_ID = \"" + real_mod_ID + "\"");
          modification_ID = sql_adapter_->getUnaryResult();

          //add current realized modification to actual modification combination
          mod_comb.push_back(pair<int,int>(atoi(m_position.c_str()), atoi(modification_ID.c_str())));

          //some output
          if (verbose)
            {
              cout << m_position << "(" << modification_ID << "), ";
            }

          //fetch next realized modification
          sql_adapter_->executeQuery("SELECT next_realized_modification_ID FROM " + (string)REALIZED_MOD_TABLE
                                     + " WHERE realized_modification_ID = \"" + real_mod_ID + "\"");
          real_mod_ID = sql_adapter_->getUnaryResult();
        }

      //apply actual modification combination to each fragment and calculate mass of modified fragment
      int frag_count = 0;
      for (vector<vector<int> >::iterator frag_it = fragments.begin(); frag_it != fragments.end(); frag_it++)
        {
          //fragment mass starts with unmodifified mass of fragment
          double fragment_mass = fragment_unmodified_masses[frag_count];

          //iterate over actual modification_combination
          for (vector<pair<int,int> >::iterator mod_it = mod_comb.begin(); mod_it != mod_comb.end(); mod_it++)
            {
              //is actual realized modification applicable to actual fragment (i.e. is its position between start and end of fragment)?
              if ((mod_it->first >= (*frag_it)[1]) && (mod_it->first <= (*frag_it)[2]))
                {
                  //change ID of instance of \c Modification to given ID
                  mod.changeID(mod_it->second);

                  //average or mono mass?
                  if (masstype == "mono")
                    {
                      fragment_mass += mod.getMonoMass(0);
                      fragment_mass -= mod.getMonoMass(1);
                    }
                  else
                    {
                      fragment_mass += mod.getAverageMass(0);
                      fragment_mass -= mod.getAverageMass(1);
                    }
                }
            }

          //NOW WE HAVE CALCULATED THE MASS TO A PARTICULAR ANNOTATION!!!!
          //some output
          if (verbose)
            {
              cout << "Mass of modified Fragment with ID: " << (*frag_it)[0] << " is: " << fragment_mass << endl;
            }

          //check if any residue was modified
          if (fragment_mass == fragment_unmodified_masses[frag_count]) //no residue was modified
            {
              //check if unmodified fragment already present in annotations
              if (!added_unmodified[frag_count])
                {
                  //add annotation of unmodified fragment to database
                  sql_adapter_->executeQuery("INSERT INTO " + (string)ANNOTATION_TABLE
                                             + " ( `sample_ID` , `mass` , `digest_fragment_ID` , `realized_modification_ID` ) "
                                             + " VALUES ( '" + String(sample_ID) + "', '" + String(fragment_mass)
                                             + "', '" + String((*frag_it)[0]) + "', '-1' )");

                  //unmodified fragment does not have to be added again!
                  added_unmodified[frag_count] = true;
                }
            }
          else //fragment was modified
            {
              //add annotation to database
              sql_adapter_->executeQuery("INSERT INTO " + (string)ANNOTATION_TABLE
                                         + " ( `sample_ID` , `mass` , `digest_fragment_ID` , `realized_modification_ID` ) "
                                         + " VALUES ( '" + String(sample_ID) + "', '" + String(fragment_mass)
                                         + "', '" + String((*frag_it)[0]) + "', '" + first_realized_mod_ID +"' )");
            }

          frag_count++;
        }

      //fetch next modification combination
      sql_adapter_->executeQuery("SELECT next_modification_combination_ID FROM " + (string)MOD_COMB_TABLE
                                 + " WHERE modification_combination_ID = \"" + mod_comb_ID + "\"");
      mod_comb_ID = sql_adapter_->getUnaryResult();
    }
}



void Sample::annotateEnumerative_()
{
  //check: is actual sample already stored in database? if not, calculate annotations
  if (!existsInDB_())                 //after that  call, if true, sample_ID and prot_mod_scen_ID are set to right values
    {
#ifdef ANNOTATE_QT
      // Qt will delete it when done!
      oe_ = new OutputEvent("Sample::annotate(): Calculating Annotations ("+ annotation_method +")...\n");
      QApplication::postEvent(qannotate_, oe_);
#endif
      #ifndef ANNOTATE_QT

      cout << "Sample::annotate(): Calculating Annotations ("+ annotation_method +")..." << endl;
#endif

      calculateAnnotations_(masstype, false); //if above statement was false (i.e. existsInDB_() was true), annotations are already
      //generated and sample_ID and prot_mod_scen_ID are set to right values (all non verbose)
    }

  //read "real" peaklist in verbose mode
  readPeaklist_(peakfile_format, true);


  //NOW EVERYTHING IS SET FOR SEARCHING "REAL" MASSES IN THE SET OF THEORETICALLY CALCULATED MASSES, STORED IN DATABASE


  //iterate "real" peaks
#ifdef ANNOTATE_QT
  // Qt will delete it when done!
  oe_ = new OutputEvent("Sample::annotate(): Annotating Peaks...\n");
  QApplication::postEvent(qannotate_, oe_);
#endif
  #ifndef ANNOTATE_QT

  cout << "Sample::annotate(): Annotating Peaks..." << endl;
#endif

  for (DPeakArray<1, DPeak<1> >::iterator it = peaklist_.begin(); it != peaklist_.end(); it++)
    {
      annotatePeak_(*it);
    }
}



void Sample::annotatePeak_(DPeak<1>& peak)
{
  //check if conditions for use of Sample::annotatePeak_() are fulfilled
  if ((sample_ID == -1) || (!modified))
    {
      throw OpenMS::Exception::Base(__FILE__, __LINE__, __PRETTY_FUNCTION__, "No proper use of Sample::annotatePeak_()",
                                    "Sample::annotatePeak_() only can be used after calls of modify_(), and existInDB() or dbRegister_()");
    }

  //datastructure for annotations: elems in vec.: first: annotation_ID, second: digest_fragment_ID, third: (first) realized_modification_ID
  vector<vector<int> > annotations;

  //get annotations (use +/- range!!!)
  sql_adapter_->executeQuery("SELECT annotation_ID, digest_fragment_ID, realized_modification_ID FROM " + (string)ANNOTATION_TABLE
                             + " WHERE sample_ID = \"" + String(sample_ID) +
                             "\" AND mass >= \"" + String((peak.getPosition()[0]) - range) +
                             "\" AND mass <= \"" + String((peak.getPosition()[0]) + range) + "\"");

  //access query result directly
  QSqlDatabase* db_handle_ = QSqlDatabase::database("db_handle_");
  QSqlQuery res(db_handle_);
  res = sql_adapter_->lastResult();

  //check if result is valid
  if ((!res.isActive()))
    {
      cerr << "Throwing exception because of query-string: " << endl << sql_adapter_->lastQuery() << endl;
      throw OpenMS::Exception::Base(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Wrong Query for Annotations", "Could not obtain Annotations!");
    }
  //check if annotations found
  else if (res.size()==0)
    {
      return;                      //no annotations found
    }

  //fill annotations into datastructure created above
  while (res.next())
    {
      vector<int> tmp_annotation;

      tmp_annotation.push_back(res.value(0).toInt());
      tmp_annotation.push_back(res.value(1).toInt());
      tmp_annotation.push_back(res.value(2).toInt());

      annotations.push_back(tmp_annotation);
    }

  //store annotations in peaklist
  storeAnnotations_(annotations, peak);
}



void Sample::storeAnnotations_(vector<vector<int> > annotations, DPeak<1>& peak)
{
  Modification mod(1, db_username_, db_password_, db_host_); //Modification iterator (dummy ID)
  hash_map<int,int> temp_o_mods; //needed as iterator for "getFragmentOverallModifiedMass"
  vector<Annotation> annot_vec;
  for (vector<vector<int> >::iterator it = annotations.begin(); it != annotations.end(); it++)
    {
      Annotation new_annot;

      new_annot.annotation_ID = (*it)[0];
      new_annot.first_real_mod_ID = (*it)[2];
      new_annot.fragment_ID = (*it)[1];

      new_annot.annotation_method = annotation_method;
      new_annot.masstype = masstype;
      new_annot.peak_mass = peak.getPosition()[0];

      //enzyme used?
      if (enzyme == NULL)
        {
          new_annot.enzyme = "none";
        }
      else
        {
          new_annot.enzyme = enzyme->getType();
        }

      //get calculated mass information
      sql_adapter_->executeQuery("SELECT mass FROM " + (string)ANNOTATION_TABLE
                                 + " WHERE annotation_ID = \"" + String((*it)[0]) + "\"");
      string calculated_mass = sql_adapter_->getUnaryResult();

      //get digest_fragment information
      sql_adapter_->executeQuery("SELECT protein_ID FROM " + (string)FRAGMENT_TABLE
                                 + " WHERE digest_fragment_ID = \"" + String((*it)[1]) + "\"");
      string prot_ID = sql_adapter_->getUnaryResult();

      sql_adapter_->executeQuery("SELECT enzyme_ID FROM " + (string)FRAGMENT_TABLE
                                 + " WHERE digest_fragment_ID = \"" + String((*it)[1]) + "\"");
      string enz_ID = sql_adapter_->getUnaryResult();


      new_annot.calculated_annotation_mass = atof((calculated_mass).c_str());
      new_annot.protein_ID = atoi(prot_ID.c_str());


      sql_adapter_->executeQuery("SELECT d_start_pos FROM " + (string)FRAGMENT_TABLE
                                 + " WHERE digest_fragment_ID = \"" + String((*it)[1]) + "\"");
      int d_start_pos = atoi(sql_adapter_->getUnaryResult().c_str());

      sql_adapter_->executeQuery("SELECT d_end_pos FROM " + (string)FRAGMENT_TABLE
                                 + " WHERE digest_fragment_ID = \"" + String((*it)[1]) + "\"");
      int d_end_pos = atoi(sql_adapter_->getUnaryResult().c_str());

      //fragment
      new_annot.setFragment(d_start_pos, protein_digest->getResName(d_start_pos), d_end_pos, protein_digest->getResName(d_end_pos));

      if (masstype == "average")
        {
          new_annot.unmodified_fragment_mass = protein_digest->getFragmentAverageMass(d_start_pos, d_end_pos);
        }
      else if (masstype == "mono")
        {
          new_annot.unmodified_fragment_mass = protein_digest->getFragmentMonoMass(d_start_pos, d_end_pos);
        }

      new_annot.overall_modified_fragment_mass = protein_digest->getFragmentOverallModifiedMass(d_start_pos, d_end_pos, masstype,
          temp_o_mods, mod);
      new_annot.plus_mass_overall_modifications = (new_annot.overall_modified_fragment_mass
          - new_annot.unmodified_fragment_mass);
      new_annot.plus_mass_modification_combination = (new_annot.calculated_annotation_mass - new_annot.overall_modified_fragment_mass);


      //test if annotation is conform with actual instance of \c Sample
      if (prot_ID != String(protein_digest->getProteinID()))
        {
          throw OpenMS::Exception::Base(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Wrong annotation",
                                        "Wrong annotation: not the same protein ID like in actual sample");
        }
      else if ((enzyme != NULL) && (enz_ID != "-1") && (enz_ID != String(enzyme->getID())))
        {
          throw OpenMS::Exception::Base(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Wrong annotation",
                                        "Wrong annotation: not the same enzyme ID like in actual sample");
        }
      else if ((enzyme == NULL) && (enz_ID != "-1"))
        {
          throw OpenMS::Exception::Base(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Wrong annotation",
                                        "Wrong annotation: not the same enzyme ID like in actual sample (no enzyme used)");
        }

      //get protein identifier
      sql_adapter_->executeQuery("SELECT identifier FROM " + (string)PROTEIN_TABLE
                                 + " WHERE protein_ID = \"" + prot_ID + "\"");
      new_annot.protein = sql_adapter_->getUnaryResult();

      //iterate realized_modifications
      string m_position, modification_ID, real_mod_ID = String((*it)[2]);

      //if not unmodified
      if (!(real_mod_ID == "-1"))
        {
          for( ; real_mod_ID != "0"; )
            {
              //get position that is modified
              sql_adapter_->executeQuery("SELECT m_position FROM " + (string)REALIZED_MOD_TABLE
                                         + " WHERE realized_modification_ID = \"" + real_mod_ID + "\"");
              m_position = sql_adapter_->getUnaryResult();

              //only write modification into file, if it is between d_start_pos and d_end_pos
              int m_int_pos = atoi(m_position.c_str());
              if ((m_int_pos >= d_start_pos) && (m_int_pos <= d_end_pos))
                {
                  //get ID of Modification
                  sql_adapter_->executeQuery("SELECT modification_ID FROM " + (string)REALIZED_MOD_TABLE
                                             + " WHERE realized_modification_ID = \"" + real_mod_ID + "\"");
                  modification_ID = sql_adapter_->getUnaryResult();

                  //write realized modification into file
                  mod.changeID(atoi(modification_ID.c_str()));

                  //get netto plus mass of modification
                  double mod_netto_mass = 0;
                  if (masstype == "mono") //average or mono mass?
                    {
                      mod_netto_mass += mod.getMonoMass(0);
                      mod_netto_mass -= mod.getMonoMass(1);
                    }
                  else if (masstype == "average")
                    {
                      mod_netto_mass += mod.getAverageMass(0);
                      mod_netto_mass -= mod.getAverageMass(1);
                    }

                  vector<int> tmp_positions;
                  tmp_positions.push_back(atoi(m_position.c_str()));
                  new_annot.addModification(atoi(modification_ID.c_str()), mod.getType(), mod_netto_mass, 1, tmp_positions);

                }

              //fetch next realized modification
              sql_adapter_->executeQuery("SELECT next_realized_modification_ID FROM " + (string)REALIZED_MOD_TABLE
                                         + " WHERE realized_modification_ID = \"" + real_mod_ID + "\"");
              real_mod_ID = sql_adapter_->getUnaryResult();
            }
        }

      annot_vec.push_back(new_annot);
    }

  //register "annotations" as meta-value for peak, if already registered, nothing happens
  peak.metaRegistry().registerName("annotations", "annotations found for this peak");

  //already annotations stored for this peak?
  int index = -1;
  try
    {
      //get index in annotation_vectors_ under which vector of Annotations for this peak is already stored
      index = (int)(peak.getMetaValue("annotations"));

      //append existing vector of Annotations with annot_vec
      annotation_vectors_[index].insert(annotation_vectors_[index].begin(), annot_vec.begin(), annot_vec.end());
    }
  catch (OpenMS::Exception::ConversionError)
    {
      //insert annot_vec into \c annotation_vectors_
      annotation_vectors_.push_back(annot_vec);

      //insert new index from annotation_vectors_ (pointing to annot_vec) as metaValue for actual peak
      peak.setMetaValue("annotations", (int)(annotation_vectors_.size()-1));
    }
}




/*****************************************************************************************************************************************
 * Method: "improved_enumerate" 
 ****************************************************************************************************************************************/


int Sample::modifyPartiallyImprovedEnumerate_(vector<pair<int, vector<Modification*> > > mods, bool verbose)
{
  //check if positions for partial modifications are already occupied by overall modifications and if modifications all can modify their pos
  for (vector<pair<int, vector<Modification*> > >::iterator it = mods.begin(); it != mods.end(); it++)
    {
      //check occupation by overall mod
      if (protein_digest->sequence_overall_modifications[it->first] != 0)
        {
          cout << "Site " + String(it->first) + " is to be modified by partial modifications and";
          cout <<" is also to be modified by overall modification with ID " << protein_digest->sequence_overall_modifications[it->first];
          cout << "." << endl;
          cout << "Using partial modification." << endl;
          protein_digest->sequence_overall_modifications[it->first] = 0;
        }

      //check if modifications all can modify their given position
      for (vector<Modification*>::iterator mod_it = (it->second)
           .begin();
           mod_it != (it->second).end();
           mod_it++)
        {
          //test: is actual modification able to modify actual position?
          if (!((*mod_it)->canModify(String(protein_digest->sequence_oneletter[it->first]))))
            {
              throw ProteinDigest::WrongModification(__FILE__, __LINE__, __PRETTY_FUNCTION__, (*mod_it)->getID(), it->first);
            }
        }
    }

  //group modifications and determine number of positions to fill for each group
  map<vector<int>, int> temp_groups;
  for (vector<pair<int, vector<Modification*> > >::iterator it = mods.begin(); it != mods.end(); it++)
    {
      //create int-vector out of \c Modification* vector
      vector<int> tmp_curr_group;
      for (vector<Modification*>::iterator m_it = (it->second)
           .begin();
           m_it != (it->second).end();
           m_it++)
        {
          tmp_curr_group.push_back((*m_it)->getID());
        }

      //this works since a map initializes elements on first access via []-operator
      (temp_groups[tmp_curr_group])++;
    }

  //generate input variables for \c improvedRecursiveEnumerate_
  vector<vector<int> > further_modification_sets;
  vector<int> further_free_pos;
  for (map<vector<int>,int>::iterator it = temp_groups.begin(); it != temp_groups.end(); ++it)
    {
      if (verbose)
        {
          cout << "Modification group (";
          for (vector<int>::const_iterator intit = (it->first)
               .begin();
               intit != (it->first).end();
               intit++)
            {
              cout << *intit << ", ";
            }
          cout << ") can modify " << it->second << " positions." << endl;
        }

      further_modification_sets.push_back(it->first);
      further_free_pos.push_back(it->second);
    }
  vector<pair<int,int> > accu;


  //(smart, improved) enumerate all modification combinations, where position information is omitted
  list<pair<int,int> > modification_combinations_posless =
    improvedRecursiveEnumerate_(vector<int>(), 0, accu, further_modification_sets, further_free_pos);

  if (verbose)
    {
      cout << "Number of found modification combinations " << modification_combinations_posless.size() << endl;
    }

  //fill first ID of each combination into table "modification_combination_positionless"
  string last_mod_comb_ID, actual_mod_comb_ID, first_mod_comb_ID;
  for (list<pair<int,int> >::iterator it = modification_combinations_posless.begin(); it != modification_combinations_posless.end(); it++)
    {
      // add actual combination
      sql_adapter_->executeQuery("INSERT INTO " + (string)MOD_COMB_PLESS_TAB
                                 + " ( `first_realized_modification_positionless_ID`) "
                                 + " VALUES ( '" + String(it->first) + "' )");

      //get database id of actual combination
      sql_adapter_->executeQuery("SELECT last_insert_id() FROM " + (string)MOD_COMB_PLESS_TAB + " LIMIT 1");
      string actual_mod_comb_ID = sql_adapter_->getUnaryResult();

      // if we do not work on the first combination: add pointer to actual combination into last combination
      if (it != modification_combinations_posless.begin())
        {
          sql_adapter_->executeQuery("UPDATE " + (string)MOD_COMB_PLESS_TAB
                                     + " SET `next_modification_combination_positionless_ID` = " + actual_mod_comb_ID
                                     + " WHERE `modification_combination_positionless_ID` = " + last_mod_comb_ID + " LIMIT 1");
        }

      //store this position
      last_mod_comb_ID = actual_mod_comb_ID;

      //if this is first combination: store in variable.
      if (it == modification_combinations_posless.begin())
        {
          first_mod_comb_ID = actual_mod_comb_ID;
        }
    }
  return atoi(first_mod_comb_ID.c_str());
}



list<pair<int,int> > Sample::improvedRecursiveEnumerate_(vector<int> current_possible_modification_set, int current_no_free_pos,
    vector<pair<int,int> > accu, vector<vector<int> > further_modification_sets,
    vector<int> further_free_pos)
{
  list<pair<int,int> > result;

  // if all positions are full and no further sets of modifications have to be treated
  if ((current_no_free_pos == 0) && (further_modification_sets.empty()) && (further_free_pos.empty()))
    {
      //add accu to datastructure / database
      string actual_mod_ID, last_mod_ID = "-1", first_mod_ID = "-1";
      for (vector<pair<int,int> >::iterator it = accu.begin(); it != accu.end(); it++)
        {
          // add actual position/modification
          sql_adapter_->executeQuery("INSERT INTO " + (string)REAL_MOD_PLESS_TAB
                                     + " ( `modification_ID` , `no_of_occurrences` ) "
                                     + " VALUES ( '" + String(it->first) + "', '" + String(it->second) + "' )");

          //get database id of actual modification
          sql_adapter_->executeQuery("SELECT last_insert_id() FROM " + (string)REAL_MOD_PLESS_TAB + " LIMIT 1");
          actual_mod_ID = sql_adapter_->getUnaryResult();


          // if we do not work on the first modification: add pointer to actual modification into last modification
          if (last_mod_ID != "-1")
            {
              sql_adapter_->executeQuery("UPDATE " + (string)REAL_MOD_PLESS_TAB
                                         + " SET `next_realized_modification_positionless_ID` = " + actual_mod_ID
                                         + " WHERE `realized_modification_positionless_ID` = " + last_mod_ID + " LIMIT 1");
            }

          //store this position
          last_mod_ID = actual_mod_ID;

          //if this is first modification: store in variable.
          if (first_mod_ID == "-1")
            {
              first_mod_ID = actual_mod_ID;
            }
        }

      /* put first_ID and last ID in realized_modification into result
         (first: for modification_combination, last: for adding overall mods, NOT USED ANY MORE!!)
       */
      result.push_back(pair<int,int>(atoi(first_mod_ID.c_str()),atoi(actual_mod_ID.c_str())));

      return result;
    }
  // if all positions of current modification set are full, but there are still modifications to be treated
  else if (current_no_free_pos == 0)
    {
      //recurse with new set of modification as current set of possible modifications
      vector<int> new_modification_set = *(further_modification_sets.begin());
      further_modification_sets.erase(further_modification_sets.begin());

      int new_free_pos = *(further_free_pos.begin());
      further_free_pos.erase(further_free_pos.begin());

      //construct result and recurse
      list<pair<int,int> > recursive_result =
        improvedRecursiveEnumerate_(new_modification_set, new_free_pos, accu, further_modification_sets, further_free_pos);
      result.splice(result.begin(), recursive_result);
    }
  // fill all positions associated with current set of possible modifications
  else
    {
      vector<int> tmp_mods = current_possible_modification_set;

      // if current_possible_modification_set is empty: nothing is done!!! ;ELSE: iterate possible modifications
      for (vector<int>::iterator it = current_possible_modification_set.begin(); it != current_possible_modification_set.end(); it++)
        {

          //modification currently worked on, is never to be treated again (we don`t want to generate double results)
          tmp_mods.erase(find(tmp_mods.begin(), tmp_mods.end(), (*it)));
          int tmp_mod_size = tmp_mods.size();

          for (int i = current_no_free_pos; i > 0; i--)
            {
              //populate i positions with modification *it
              vector<pair<int, int> > tmp_accu = accu;
              tmp_accu.push_back(pair<int,int>(*it,i));

              //since i positions are occupied by modification *it:
              int new_no_free_pos = (current_no_free_pos - i);

              //avoid senseless calls: if no modifications would be left for rec call, dont even execute it.
              if (!((tmp_mod_size == 0) && (new_no_free_pos != 0)))
                {
                  list<pair<int,int> > recursive_result =
                    improvedRecursiveEnumerate_(tmp_mods, new_no_free_pos, tmp_accu, further_modification_sets, further_free_pos);
                  result.splice(result.begin(), recursive_result);
                }
            }
        }
    }

  return result;
}



void Sample::improvedCalculateAnnotations_(string masstype)
{
  //is sample already registered in database? if not: register
  if (sample_ID == -1)
    {
      dbRegister_();
    }

  //datastrucuture for transiently saving a particular modification combination
  vector<pair<int,int> > mod_comb_posless;

  //instance of \c Modification for access to different Masses of Modifications (ID 1, just dummy, will be changed later in this method)
  Modification mod(1, db_username_, db_password_, db_host_);

  //create variables for modification-iteration
  string mod_comb_posless_ID, real_mod_posless_ID, no_of_occurrences, modification_ID;

  //iterate modification_combinations
  sql_adapter_->executeQuery("SELECT modification_combination_positionless_ID FROM " + (string)PROT_MOD_SCEN_TABLE
                             + " WHERE protein_modification_scenario_ID = \"" + String(prot_mod_scen_ID) + "\"");
  mod_comb_posless_ID = sql_adapter_->getUnaryResult();
  for( ; ((mod_comb_posless_ID != "0") && (mod_comb_posless_ID != "0")) ; )
    {
      //store always only one particular modification_combination
      mod_comb_posless.clear();

      //iterate realized_modifications
      sql_adapter_->executeQuery("SELECT first_realized_modification_positionless_ID FROM " + (string)MOD_COMB_PLESS_TAB
                                 + " WHERE modification_combination_positionless_ID = \"" + mod_comb_posless_ID + "\"");
      real_mod_posless_ID = sql_adapter_->getUnaryResult();

      //stores first_real_mod_ID of this mod_combination for adding to "annotation"
      string first_realized_mod_posless_ID = real_mod_posless_ID;

      for( ; real_mod_posless_ID != "0"; )
        {
          //get ID of modification
          sql_adapter_->executeQuery("SELECT modification_ID FROM " + (string)REAL_MOD_PLESS_TAB
                                     + " WHERE realized_modification_positionless_ID = \"" + real_mod_posless_ID + "\"");
          modification_ID = sql_adapter_->getUnaryResult();

          //get number of occurrences of Modification
          sql_adapter_->executeQuery("SELECT no_of_occurrences FROM " + (string)REAL_MOD_PLESS_TAB
                                     + " WHERE realized_modification_positionless_ID = \"" + real_mod_posless_ID + "\"");
          no_of_occurrences = sql_adapter_->getUnaryResult();

          //add current realized modification to actual modification combination
          mod_comb_posless.push_back(pair<int,int>(atoi(modification_ID.c_str()), atoi(no_of_occurrences.c_str())));

          //fetch next realized modification
          sql_adapter_->executeQuery("SELECT next_realized_modification_positionless_ID FROM " + (string)REAL_MOD_PLESS_TAB
                                     + " WHERE realized_modification_positionless_ID = \"" + real_mod_posless_ID + "\"");
          real_mod_posless_ID = sql_adapter_->getUnaryResult();
        }


      //variable for summing up mass of this modification combination
      double mod_comb_mass = 0.0;

      //iterate over actual modification_combination
      for (vector<pair<int,int> >::iterator mod_it = mod_comb_posless.begin(); mod_it != mod_comb_posless.end(); mod_it++)
        {
          //change ID of instance of \c Modification to given ID
          mod.changeID(mod_it->first);

          //average or mono mass?
          if (masstype == "mono")
            {
              mod_comb_mass += ((mod_it->second) * (mod.getMonoMass(0)));
              mod_comb_mass -= ((mod_it->second) * (mod.getMonoMass(1)));
            }
          else
            {
              mod_comb_mass += ((mod_it->second) * (mod.getAverageMass(0)));
              mod_comb_mass -= ((mod_it->second) * (mod.getAverageMass(1)));
            }
        }


      //NOW WE HAVE CALCULATED THE MASS TO A PARTICULAR MODIFICATION COMBINATION!!!!


      //add annotation to database: IN THIS METHOD W I T H O U T digest_fragment_ID, since only mass of modification combination is stored
      sql_adapter_->executeQuery("INSERT INTO " + (string)ANNOTATION_TABLE
                                 + " ( `sample_ID` , `mass` , `realized_modification_positionless_ID` ) "
                                 + " VALUES ( '" + String(sample_ID) + "', '" + String(mod_comb_mass)
                                 + "', '" + first_realized_mod_posless_ID +"' )");

      //fetch next modification combination
      sql_adapter_->executeQuery("SELECT next_modification_combination_positionless_ID FROM " + (string)MOD_COMB_PLESS_TAB
                                 + " WHERE modification_combination_positionless_ID = \"" + mod_comb_posless_ID + "\"");
      mod_comb_posless_ID = sql_adapter_->getUnaryResult();
    }
}



void Sample::annotateEnumerativeImproved_()
{
  /* in this method we calculate an annotation only once for each modification combination. peaks are then searched for each fragment
     separately, using suitable offset. 
     so not to calculate annotations more often than needed, for storing in database always NO ENZYME is assumed 
   */
  Enzyme* tmp_enzyme = enzyme;
  enzyme = NULL;
  bool tmp_digested = digested;
  digested = false;

  //check: is actual sample (with no respect to enzyme) already stored in database? if not, calculate annotations
  if (!existsInDB_())
    {
#ifdef ANNOTATE_QT
      // Qt will delete it when done!
      oe_ = new OutputEvent("Sample::annotate(): Calculating Annotations ("+ annotation_method +")...\n");
      QApplication::postEvent(qannotate_, oe_);
#endif
      #ifndef ANNOTATE_QT

      cout << "Sample::annotate(): Calculating Annotations ("+ annotation_method +")..." << endl;
#endif

      improvedCalculateAnnotations_(masstype);
    }

  //restore true enzyme settings
  enzyme = tmp_enzyme;
  digested = tmp_digested;

  bool whole_protein; //flag, for optimizing annotation: if whole protein, isFragmentModifyable_ does not have to be executed

  //datastructure for fragments in "digest_fragments", first: ID, second: start pos, third: end pos, fourth: protein_ID, fifth: enzyme_ID
  vector<vector<int> > fragments = getFragments_(whole_protein);

  //read "real" peaklist in verbose mode
  readPeaklist_(peakfile_format, true);

  /* store information about POSITIONS of partial modifications,
     to later check, wheter enough mod sites are in given frag for given modification_combination
   */
  if (partial_mods.size()==0) // if already in database, parsing not yet done
    {
      //parse partial modifications:
      ModificationStringParser mod_parser_(db_username_, db_password_, db_host_);
      partial_mods = mod_parser_.parse(partial_modification_string);
    }
  storePartialModsPosInfo_();

  //iterate fragments
#ifdef ANNOTATE_QT
  // Qt will delete it when done!
  oe_ = new OutputEvent("Sample::annotate(): annotateEnumerativeImproved_(): Annotating Peaks...\n");
  QApplication::postEvent(qannotate_, oe_);
#endif
  #ifndef ANNOTATE_QT

  cout << "Sample::annotate(): annotateEnumerativeImproved_(): Annotating Peaks..." << endl;
#endif

  int frag_count = 0;
  Modification mod_iter(1, db_username_, db_password_, db_host_); //create instance of \c Modification for iteration in \c ProteinDigest::getFragmentOverallModifiedMass
  for (vector<vector<int> >::iterator frag_it = fragments.begin(); frag_it != fragments.end(); frag_it++)
    {
      //! store information to actual fragment
      storeFragmentInfo_((*frag_it)[1], (*frag_it)[2]);

      //! get mass of overall modified fragment
      hash_map<int,int> overall_mods;
      double fragment_overall_mod_mass =
        protein_digest->getFragmentOverallModifiedMass((*frag_it)[1], (*frag_it)[2], masstype, overall_mods, mod_iter);

      //! get mass of unmodified fragment
      double fragment_unmod_mass = 0;
      if (masstype == "average")
        {
          fragment_unmod_mass = protein_digest->getFragmentAverageMass((*frag_it)[1], (*frag_it)[2]);
        }
      else if (masstype == "mono")
        {
          fragment_unmod_mass = protein_digest->getFragmentMonoMass((*frag_it)[1], (*frag_it)[2]);
        }

      //iterate "real" peaks
      for (DPeakArray<1, DPeak<1> >::iterator it = peaklist_.begin(); it != peaklist_.end(); it++)
        {
          improvedAnnotatePeak_(*it, *frag_it, fragment_overall_mod_mass, fragment_unmod_mass, whole_protein, overall_mods);
        }
      frag_count++;
    }
}



void Sample::improvedAnnotatePeak_(DPeak<1>& peak, vector<int> fragment, double fragment_mass, double fragment_unmod_mass,
                                   bool whole_protein, hash_map<int,int> ov_mods)
{
  //check if conditions for use of Sample::improvedAnnotatePeak_() are fulfilled
  if ((sample_ID == -1) || (!modified))
    {
      throw OpenMS::Exception::Base(__FILE__, __LINE__, __PRETTY_FUNCTION__, "No proper use of Sample::improvedAnnotatePeak_()",
                                    "Sample::annotatePeak_() only can be used after calls of modify_(), and existInDB() or dbRegister_()");
    }

  //datastructure for annotations: elems in vec.: first: annotation_ID, second: mass, third: realized_modification_positionless_ID
  vector<vector<string> > annotations;

  //get annotations (use +/- range!!!) and use OFFSET ACCORDING TO FRAGMENT_MASS!!
  sql_adapter_->executeQuery("SELECT annotation_ID, mass, realized_modification_positionless_ID FROM " + (string)ANNOTATION_TABLE
                             + " WHERE sample_ID = \"" + String(sample_ID) +
                             "\" AND mass >= \"" + String((peak.getPosition()[0] - fragment_mass) - range) +
                             "\" AND mass <= \"" + String((peak.getPosition()[0] - fragment_mass) + range) + "\"");

  //access query result directly
  QSqlDatabase* db_handle_ = QSqlDatabase::database("db_handle_");
  QSqlQuery res(db_handle_);
  res = sql_adapter_->lastResult();

  //check if result is valid
  if ((!res.isActive()))
    {
      cerr << "Throwing exception because of query-string: " << endl << sql_adapter_->lastQuery() << endl;
      throw OpenMS::Exception::Base(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Wrong Query for Annotations", "Could not obtain Annotations!");
    }
  //check if annotations found
  else if (res.size()==0)
    {
      return;                      //no annotations found
    }

  //fill annotations into datastructure created above
  while(res.next())
    {
      vector<string> tmp_annotation;

      ostringstream s1, s2, s3;
      s1 << res.value(0).toString();
      s2 << res.value(1).toString();
      s3 << res.value(2).toString();
      tmp_annotation.push_back((string)(s1.str())); //annotation_ID
      tmp_annotation.push_back((string)(s2.str())); //mass
      tmp_annotation.push_back((string)(s3.str())); //realized_modification_positionless_ID

      annotations.push_back(tmp_annotation);
    }

  int annotation_count = 1;
  Modification mod(1, db_username_, db_password_, db_host_);      //Modification iterator (dummy ID)
  ostringstream ost_master; //stringstream for all annotations

  //datastructure, index of position in \c annotation_vectors_  of which is stored as metadata in currend DPeak (peak)
  vector<Annotation> annot_vec;

  for (vector<vector<string> >::iterator it = annotations.begin(); it != annotations.end(); it++)
    {
      bool annotation_valid = true;
      Annotation new_annot;

      new_annot.annotation_ID = atoi(((*it)[0]).c_str());
      new_annot.first_real_mod_pless_ID = atoi(((*it)[2]).c_str());
      new_annot.fragment_ID = fragment[0];
      new_annot.protein_ID = fragment[3];
      new_annot.annotation_method = annotation_method;

      new_annot.masstype = masstype;
      new_annot.peak_mass = peak.getPosition()[0];
      new_annot.calculated_annotation_mass = (atof(((*it)[1]).c_str()) + fragment_mass);
      new_annot.unmodified_fragment_mass = fragment_unmod_mass;
      new_annot.overall_modified_fragment_mass = fragment_mass;
      new_annot.plus_mass_overall_modifications = (fragment_mass - fragment_unmod_mass);
      new_annot.plus_mass_modification_combination = atof(((*it)[1]).c_str());

      //get protein identifier
      sql_adapter_->executeQuery("SELECT identifier FROM " + (string)PROTEIN_TABLE
                                 + " WHERE protein_ID = \"" + String(fragment[3]) + "\"");
      new_annot.protein = sql_adapter_->getUnaryResult();

      //enzyme used?
      if (enzyme == NULL)
        {
          new_annot.enzyme = "none";
        }
      else
        {
          new_annot.enzyme = enzyme->getType();
        }

      //fragment
      new_annot.setFragment(fragment[1], protein_digest->getResName(fragment[1]), fragment[2], protein_digest->getResName(fragment[2]));

      //iterate realized_modifications
      string no_of_occurrences, modification_ID, real_mod_posless_ID = (*it)[2];

      //only if modified
      if (!((real_mod_posless_ID == "-1") || (((*it)[1]) == "0"))) //is mass of modification combination == 0???
        {
          //datastructure for saving what and how many modifications occured
          map<int, int> mod_occurrences;

          for( ; real_mod_posless_ID != "0"; )
            {
              //get modification_ID
              sql_adapter_->executeQuery("SELECT modification_ID FROM " + (string)REAL_MOD_PLESS_TAB
                                         + " WHERE realized_modification_positionless_ID = \"" + real_mod_posless_ID + "\"");
              modification_ID = sql_adapter_->getUnaryResult();

              //get number of occurrences
              sql_adapter_->executeQuery("SELECT no_of_occurrences FROM " + (string)REAL_MOD_PLESS_TAB
                                         + " WHERE realized_modification_positionless_ID = \"" + real_mod_posless_ID + "\"");
              no_of_occurrences = sql_adapter_->getUnaryResult();

              //get modfication type and netto mass
              mod.changeID(atoi(modification_ID.c_str()));
              string modification_type = mod.getType();
              double mod_netto_mass = 0;
              if (masstype == "mono") //average or mono mass?
                {
                  mod_netto_mass += mod.getMonoMass(0);
                  mod_netto_mass -= mod.getMonoMass(1);
                }
              else if (masstype == "average")
                {
                  mod_netto_mass += mod.getAverageMass(0);
                  mod_netto_mass -= mod.getAverageMass(1);
                }

              new_annot.addModification(atoi(modification_ID.c_str()), modification_type, mod_netto_mass,
                                        atoi(no_of_occurrences.c_str()), vector<int>());

              if (!((modification_type == "unmodified") || (mod_netto_mass == 0)))
                {
                  //only store actual modification if actual mod is NOT the artificial modification "unmodified"
                  mod_occurrences[atoi(modification_ID.c_str())] += atoi(no_of_occurrences.c_str());
                }

              //fetch next realized modification
              sql_adapter_->executeQuery("SELECT next_realized_modification_positionless_ID FROM " + (string)REAL_MOD_PLESS_TAB
                                         + " WHERE realized_modification_positionless_ID = \"" + real_mod_posless_ID + "\"");
              real_mod_posless_ID = sql_adapter_->getUnaryResult();
            }

          /* check, if fragment is modifyable by actual modification combination
             (if it is the whole protein, this check does not have to be performed)
           */
          if (!whole_protein)
            {
              //remember: this call only makes sense, if "storeFragmentInfo" is called before!!
              if (!(isActualFragmentModifyable_(mod_occurrences)))
                {
                  annotation_valid = false;
                }
            }
        }

      if (annotation_valid)
        {
          annotation_count++;

          //add overall modifications into stringstream
          for (hash_map<int,int>::iterator it = ov_mods.begin(); it != ov_mods.end(); it++)
            {
              mod.changeID(it->first);

              //get netto plus mass of modification
              double mod_netto_mass = 0;
              if (masstype == "mono") //average or mono mass?
                {
                  mod_netto_mass += mod.getMonoMass(0);
                  mod_netto_mass -= mod.getMonoMass(1);
                }
              else if (masstype == "average")
                {
                  mod_netto_mass += mod.getAverageMass(0);
                  mod_netto_mass -= mod.getAverageMass(1);
                }

              new_annot.addModification((it->first), mod.getType(), mod_netto_mass, (it->second), vector<int>());

            }

          //add whole annotation to vector that will eventually be stored as meta info for actual peak
          annot_vec.push_back(new_annot);
        }
    }

  //only store annot_vec as meta info to this peak, if at least one found annotation was valid
  //(if == 1: ONCE annotation_count--,and once incremented: no annotation)
  if (annotation_count > 1)
    {
      //register "annotations" as meta-value for peak, if already registered, nothing happens
      peak.metaRegistry().registerName("annotations", "annotations found for this peak");

      //already annotations stored for this peak?
      int index = -1;
      try
        {
          //get index in annotation_vectors_ under which vector of Annotations for this peak is already stored
          index = (int)(peak.getMetaValue("annotations"));

          //append existing vector of Annotations with annot_vec
          annotation_vectors_[index].insert(annotation_vectors_[index].begin(), annot_vec.begin(), annot_vec.end());
        }
      catch (OpenMS::Exception::ConversionError)
        {
          //insert annot_vec into \c annotation_vectors_
          annotation_vectors_.push_back(annot_vec);

          //insert new index from annotation_vectors_ (pointing to annot_vec) as metaValue for actual peak
          peak.setMetaValue("annotations", (int)(annotation_vectors_.size()-1));
        }
    }
}



void Sample::storePartialModsPosInfo_()
{
  //just to be sure: parse partial modifiations (possibly again)
  ModificationStringParser mod_parser_(db_username_, db_password_, db_host_);
  partial_mods = mod_parser_.parse(partial_modification_string);

  //fill \c modification_positions
  for (vector<pair<int, vector<Modification*> > >::iterator x = partial_mods.begin(); x != partial_mods.end(); x++)
    {
      for (vector<Modification*>::iterator y = (x->second)
           .begin();
           y != (x->second).end();
           y++)
        {
          //store for each modification a vector of positions at which this modification may occur.
          modification_positions[(*y)->getID()].push_back(x->first);
        }
    }

  //fill \c partial_mods_int
  for (vector<pair<int, vector<Modification*> > >::iterator it = partial_mods.begin(); it != partial_mods.end(); it++)
    {
      //create int-vector out of \c Modification* vector
      vector<int> tmp_curr_group;
      for (vector<Modification*>::iterator m_it = (it->second)
           .begin();
           m_it != (it->second).end();
           m_it++)
        {
          tmp_curr_group.push_back((*m_it)->getID());
        }

      partial_mods_int.push_back(pair<int, vector<int> >(it->first, tmp_curr_group));
    }
}



void Sample::storeFragmentInfo_(int start_pos, int end_pos)
{
  //erase entries from former calls of \c storeFragmentInfo()
  actual_fragment_partial_mods_int.clear();
  actual_fragment_groups.clear();
  actual_fragment_mod_with_groups.clear();

  //create a \c partial_mods_int especially for this fragment
  for (vector<pair<int, vector<int> > >::iterator it = partial_mods_int.begin(); it != partial_mods_int.end(); it++)
    {
      if ((it->first >= start_pos) && (it->first <= end_pos))
        {
          actual_fragment_partial_mods_int.push_back(*it);
        }
    }

  //group modifications and determine number of positions to fill for each group
  for (vector<pair<int, vector<int> > >::iterator it = actual_fragment_partial_mods_int.begin();
       it != actual_fragment_partial_mods_int.end(); it++)
    {
      //this works since a map initializes elements on first access via []-operator
      (actual_fragment_groups[it->second])++;
    }

  //create hash_map from modification-ID to a vector of group-ID's in which this modification is
  int group_count = 0;  //NOTE: group_count here represents an ID for each group, corresponding to their order in \c actual_fragment_groups
  for (map<vector<int>, int>::iterator itt = actual_fragment_groups.begin(); itt != actual_fragment_groups.end(); itt++)
    {
      for (vector<int>::const_iterator intit = itt->
           first.begin();
           intit != itt->first.end();
           intit++)
        {
          actual_fragment_mod_with_groups[*intit].push_back(group_count);
        }
      group_count++;
    }
}



bool Sample::isActualFragmentModifyable_(map<int, int> mod_occurrences)
{
  //create temporary map of open positions per group in which open positions can be reduced with each modification seen
  hash_map<int, int> temp_group_pos;
  int group_count = 0;
  for (map<vector<int>, int>::iterator it = actual_fragment_groups.begin(); it != actual_fragment_groups.end(); it++)
    {
      temp_group_pos[group_count] = it->second;
      group_count++;
    }

  //iterate modifications
  for (map<int, int>::iterator it = mod_occurrences.begin(); it != mod_occurrences.end(); it++)
    {
      //iterate all groups in which modification with it it->first occurrs
      for (vector<int>::iterator intit = actual_fragment_mod_with_groups[it->first]
                                         .begin();
           intit != actual_fragment_mod_with_groups[it->first].end();
           intit++)
        {
          //is no of free positions in found group greater than multiplicity of actual modification (it->second)?
          if (temp_group_pos[*intit] >= (it->second))
            {
              temp_group_pos[*intit] -= (it->second); //assign all multiplicities of it->first to group *intit
              (it->second) = 0;
              break;                                   //no need to look into further groups, go to next modification
            }
          else
            {
              (it->second) -= temp_group_pos[*intit]; //assign as many multiplicities of it->first to group *intit as possible
              temp_group_pos[*intit] = 0;             //with the rest of the multiplicities: look into further groups
            }
        }

      //if not all multiplicities of it->first have been assigned to groups: fragment is NOT modifyable
      if ((it->second)
          > 0)
        {
          return false;
        }
    }
  return true;
}



//old, not used any more!!!
bool Sample::isFragmentModifyable_(int start_pos, int end_pos, int mod_id, int no_of_occurrences)
{
  //calculate number of positions between start_pos and end_pos, that can be modifyed with modification with id mod_id
  int no_of_pos = 0;
  for (vector<int>::iterator it = modification_positions[mod_id].begin(); it != modification_positions[mod_id].end(); it++)
    {
      if (((*it) >= start_pos) && ((*it) <= end_pos))
        {
          no_of_pos++;
        }
    }
  return (no_of_occurrences <= no_of_pos);
}



//! old, not used any more
int Sample::getTotalNumberOfModSites_(int start_pos, int end_pos)
{
  int total_no_mod_sites = 0;
  for (vector<pair<int, vector<Modification*> > >::iterator x = partial_mods.begin(); x != partial_mods.end(); x++)
    {
      if ((x->first >= start_pos) && (x->first <= end_pos))
        {
          total_no_mod_sites++;
        }
    }

  // cout << "between " << start_pos << " and " << end_pos << ": " << total_no_mod_sites << endl;

  return total_no_mod_sites;
}




// ******************************************************************************************************************************************
// * Method: "peakwise_cormen": modification of T.H. Cormen's "Subset Sum" - Algorithm (Cormen, "Introduction to Algorithms, p. 1045)
// ******************************************************************************************************************************************


//! generates a key-string out of a concrete (modification) combination, for storing in \c cormen_temp_combinations
string Sample::generateKey_(list<int> combination_IDs)
{
  combination_IDs.sort();
  ostringstream key_st;

  for (list<int>::iterator it = combination_IDs.begin(); it != combination_IDs.end(); it++)
    {
      key_st << *it;
    }

  return key_st.str();
}



/* checks wheter given modification combination l does not exceed position limits for each modification group
   "list<int>"  : contains ID's of modifications of this combination
   "vector<int>": indices: groups; values: how many modifications of corresponding group are already contained in "list<int>"
 */
bool Sample::satisfiesGroupPos_(const pair<list<int>, vector<int> >& l)
{
  for (unsigned int group_count = 0; group_count < cormen_groups_positions.size(); group_count++)
    {
      if ((l.second)[group_count] > cormen_groups_positions[group_count])
        {
          return false;
        }
    }
  return true;
}



//! adds next Modification to given list L of modification combinations (cf. cormen: "L + x_i")
list<pair<double, pair<list<int>, vector<int> > > > Sample::addModification_(list<pair<double, pair<list<int>, vector<int> > > > L,
    double x_i, int x_i_ID, int x_i_group)
{
  //for deleting elements out of L (if combination in L appended by x_i was already seen (saved in \c cormen_temp_combinations))
  vector<list<pair<double, pair<list<int>, vector<int> > > >::iterator> trash;

  //iterate L and append each combination by x_i, if not already seen (saved in \c cormen_temp_combinations)
  for (list<pair<double, pair<list<int>, vector<int> > > >::iterator it = L.begin(); it != L.end(); it++)
    {
      list<int> tmp = it->second.first;
      tmp.push_back(x_i_ID);
      string temp_key = generateKey_(tmp);
      if (!cormen_temp_combinations[temp_key])  // if this combination was never seen before, add it.
        {
          it->first += x_i;
          it->second.first.push_back(x_i_ID);
          it->second.second[x_i_group]++;

          cormen_temp_combinations[temp_key] = true;
        }
      else                                      // if it was there before, delete corresponding element out of L
        {
          trash.push_back(it);
        }
    }

  //remove objects in trash
  for(vector<list<pair<double, pair<list<int>, vector<int> > > >::iterator>::iterator it = trash.begin(); it != trash.end(); it++)
    {
      L.erase(*it);
    }

  return L;
}



/* this is the main function of the annotation method "peakwise_cormen"
   return value: "list<int>": what modifications do realize (masses sum up to) "double"
                 "vector<int>": what groups (indices in vector), and how many mod`s per group (values in vector)
 */
list<pair<double, pair<list<int>, vector<int> > > > Sample::exactSubsetSum_(double t, double range, bool verbose)
{
  //clear temporary datastructure from data of previous runs
  cormen_temp_combinations.clear();

  list<pair<double, pair<list<int>, vector<int> > > > L;

  //Initialize L with 0-element (sum 0, no modifications)
  L.push_back(pair<double, pair<list<int>, vector<int> > >(0.0, pair<list<int>, vector<int> >(list<int>(), vector<int>())));

  //to be sure all vector<int>'s signifying groups (indices: groups, values: #mods of group, already used in curr. comb.) have right size
  L.begin()->second.second.resize(cormen_groups_positions.size());

  //iterate all modifications (sample specific, therefore stored in member variable of class \c "Sample"
  for (vector<pair<int,pair<int,double> > >::iterator it = cormen_modifications.begin(); it != cormen_modifications.end(); it++)
    {
      list<pair<double, pair<list<int>, vector<int> > > > L_temp = addModification_(L, (it->second).second, it->first, (it->second).first);
      L.merge(L_temp);

      if (verbose)
        {
          cout << "L:";
        }

      /* remove every element greater than t + range and every element that does not satisfy given group position numbers from L:
         store iterators for later deletion (else the loop would not terminate)
       */
      vector<list<pair<double, pair<list<int>, vector<int> > > >::iterator> trash;
      for (list<pair<double, pair<list<int>, vector<int> > > >::iterator itt = L.begin(); itt != L.end(); itt++)
        {
          if ((itt->first > (t + range)) || (!satisfiesGroupPos_(itt->second)))
            {
              trash.push_back(itt);

              /* erase modification combination from \c cormen_temp_combinations:
                 with doing so, \c L has always one entry more than \c cormen_temp_combinations (the 0-element!)
               */
              cormen_temp_combinations.erase(generateKey_(itt->second.first));
            }

          if (verbose)
            {
              cout << itt->first << ", ";
            }
        }

      if (verbose)
        {
          cout << endl;
        }

      //actually delete elements from L.
      for(vector<list<pair<double, pair<list<int>, vector<int> > > >::iterator>::iterator it = trash.begin(); it != trash.end(); it++)
        {
          //erase from L
          L.erase(*it);
        }
    }

  if (verbose)
    {
      cout << "L.size(): " << L.size() << endl;
      cout << "cormen_temp_combination.size(): " << cormen_temp_combinations.size() << endl;
    }

  //return only elements that are greater than t-range:
  list<pair<double, pair<list<int>, vector<int> > > > result;
  for (list<pair<double, pair<list<int>, vector<int> > > >::reverse_iterator rit = L.rbegin(); rit != L.rend(); rit++)
    {
      /* since L is sorted w.r.t. the "double" component of the pairs in L (cf. cormen):
         we dont have to look furter when 1st element out of range is found
       */
      if (rit->first < (t-range))
        {
          break;
        }
      result.push_back(*rit);
    }
  return result;
}



//! fills \c cormen_modifications and cormen_groups_positions (for fragment signified by \c start_pos and \c end_pos);
void Sample::fillCormenVariables_(int start_pos, int end_pos, Modification& mod_it)
{
  //! clear variables from data of earlier calls
  cormen_modifications.clear();
  cormen_groups_positions.clear();

  //group modifications and determine number of positions to fill for each group
  map<vector<int>, int> temp_groups;
  for (vector<pair<int, vector<Modification*> > >::iterator it = partial_mods.begin(); it != partial_mods.end(); it++)
    {
      //check if actual position is within actual fragment
      if ((it->first >= start_pos) && (it->first <= end_pos))
        {
          //create int-vector out of \c Modification* vector
          vector<int> tmp_curr_group;
          for (vector<Modification*>::iterator m_it = (it->second)
               .begin();
               m_it != (it->second).end();
               m_it++)
            {
              tmp_curr_group.push_back((*m_it)->getID());
            }

          //this works since a map initializes elements on first access via []-operator
          (temp_groups[tmp_curr_group])++;
        }
    }

  //actually fill \c cormen_modifications
  int group_count = 0;
  for (map<vector<int>,int>::iterator it = temp_groups.begin(); it != temp_groups.end(); ++it)
    {
      for (vector<int>::const_iterator intit = (it->first)
           .begin();
           intit != (it->first).end();
           intit++)
        {
          mod_it.changeID(*intit);
          double mass = 0;
          if (masstype == "average")
            {
              mass += mod_it.getAverageMass(0);
              mass -= mod_it.getAverageMass(1);
            }
          else if (masstype == "mono")
            {
              mass += mod_it.getMonoMass(0);
              mass -= mod_it.getMonoMass(1);
            }
          else
            {
              throw OpenMS::Exception::Base(__FILE__, __LINE__, __PRETTY_FUNCTION__, "Unknown Masstype",
                                            ("Masstype " + masstype + " not known!").c_str());
            }

          //add each modification as often as modifyable positions in group (they possibly all could be modified by same modification)
          for (int i = 0; i < it->second; i++)
            {
              cormen_modifications.push_back(pair<int,pair<int,double> >(*intit, pair<int,double>(group_count, mass)));
            }
        }

      //fill cormen_groups_positions
      cormen_groups_positions.push_back(it->second);

      group_count++;
    }
}



void Sample::storeAnnotationsPeakwiseCormen_(list<pair<double, pair<list<int>, vector<int> > > > modification_combinations, DPeak<1>& peak,
    vector<int> fragment, double fragment_mass, double fragment_unmod_mass, Modification& mod,
    hash_map<int,int> ov_mods)
{
  int annotation_count = 1;

  //datastructure, index of position in \c annotation_vectors_  of which is stored as metadata in currend DPeak (peak)
  vector<Annotation> annot_vec;

  for (list<pair<double, pair<list<int>, vector<int> > > >::iterator it = modification_combinations.begin();
       it != modification_combinations.end(); it++)
    {
      Annotation new_annot;

      new_annot.fragment_ID = fragment[0];
      new_annot.protein_ID = fragment[3];
      new_annot.annotation_method = annotation_method;

      new_annot.masstype = masstype;
      new_annot.peak_mass = peak.getPosition()[0];
      new_annot.calculated_annotation_mass = (it->first + fragment_mass);
      new_annot.unmodified_fragment_mass = fragment_unmod_mass;
      new_annot.overall_modified_fragment_mass = fragment_mass;
      new_annot.plus_mass_overall_modifications = (fragment_mass - fragment_unmod_mass);
      new_annot.plus_mass_modification_combination = it->first;

#ifndef ANNOTATE_XML
      //get protein identifier
      sql_adapter_->executeQuery("SELECT identifier FROM " + (string)PROTEIN_TABLE
                                 + " WHERE protein_ID = \"" + String(fragment[3]) + "\"");
      new_annot.protein = sql_adapter_->getUnaryResult();
#endif

#ifdef ANNOTATE_XML

      new_annot.protein = protein_digest->getProteinIdentifier();
#endif

      //enzyme used?
      if (enzyme == NULL)
        {
          new_annot.enzyme = "none";
        }
      else
        {
          new_annot.enzyme = enzyme->getType();
        }

      //fragment
      new_annot.setFragment(fragment[1], protein_digest->getResName(fragment[1]), fragment[2], protein_digest->getResName(fragment[2]));

      //if unmodified
      if (!(it->first == 0.0)) //is mass of modification combination == 0???
        {
          string mod_type;
          int mod_id = 0;
          double mass = 0;
          int mod_count = 1;
          for (list<int>::iterator itt = it->
                                         second.first.begin();
               itt != it->second.first.end();
               itt++)
            {
              mod.changeID(*itt);
              if ((mod.getType()) == mod_type)
                {
                  mod_count++;
                }
              else if (mod_type != "")
                {
                  new_annot.addModification(mod_id, mod_type, mass, mod_count, vector<int>());

                  //calculate netto mass of modification
                  mass = 0;
                  if (masstype == "average")
                    {
                      mass += mod.getAverageMass(0);
                      mass -= mod.getAverageMass(1);
                    }
                  else if (masstype == "mono")
                    {
                      mass += mod.getMonoMass(0);
                      mass -= mod.getMonoMass(1);
                    }

                  mod_count = 1;
                  mod_type = mod.getType();
                  mod_id = *itt;
                }
              else
                {
                  //calculate netto mass of modification
                  mass = 0;
                  if (masstype == "average")
                    {
                      mass += mod.getAverageMass(0);
                      mass -= mod.getAverageMass(1);
                    }
                  else if (masstype == "mono")
                    {
                      mass += mod.getMonoMass(0);
                      mass -= mod.getMonoMass(1);
                    }

                  mod_type = mod.getType();
                  mod_id = *itt;
                }
            }

          //write last modification into file
          new_annot.addModification(mod_id, mod_type, mass, mod_count, vector<int>());

        }

      //add overall modifications into stringstream
      for (hash_map<int,int>
           ::iterator it = ov_mods.begin();
           it != ov_mods.end();
           it++)
        {
          mod.changeID(it->first);

          double mass = 0;
          if (masstype == "average")
            {
              mass += mod.getAverageMass(0);
              mass -= mod.getAverageMass(1);
            }
          else if (masstype == "mono")
            {
              mass += mod.getMonoMass(0);
              mass -= mod.getMonoMass(1);
            }

          new_annot.addModification((it->first), mod.getType(), mass, (it->second), vector<int>());
        }
      annotation_count++;

      //add annotation into annotation_vector
      annot_vec.push_back(new_annot);
    }

  //register "annotations" as meta-value for peak, if already registered, nothing happens
  peak.metaRegistry().registerName("annotations", "annotations found for this peak");

  //already annotations stored for this peak?
  int index = -1;
  try
    {
      //get index in annotation_vectors_ under which vector of Annotations for this peak is already stored
      index = (int)(peak.getMetaValue("annotations"));

      //append existing vector of Annotations with annot_vec
      annotation_vectors_[index].insert(annotation_vectors_[index].begin(), annot_vec.begin(), annot_vec.end());
    }
  catch (OpenMS::Exception::ConversionError)
    {
      //insert annot_vec into \c annotation_vectors_
      annotation_vectors_.push_back(annot_vec);

      //insert new index in annotation_vectors_ for annot_vec as metaValue for actual peak
      peak.setMetaValue("annotations", (int)(annotation_vectors_.size()-1));
    }

#ifdef ANNOTATE_QT
  // Qt will delete it when done!
  oe_ = new OutputEvent("Annotations for peak at " + String(peak.getPosition()[0]) +  " Daltons found.\n");
  QApplication::postEvent(qannotate_, oe_);
#endif
  #ifndef ANNOTATE_QT

  cout << "Annotations for peak at " + String(peak.getPosition()[0]) +  " Daltons found." << endl;
#endif

}



void Sample::annotatePeakwiseCormen_()
{
#ifdef ANNOTATE_QT
  // Qt will delete it when done!
  oe_ = new OutputEvent("Sample::annotate(): annotatePeakwiseCormen_():\n");
  QApplication::postEvent(qannotate_, oe_);
#endif
  #ifndef ANNOTATE_QT

  cout << "Sample::annotate(): annotatePeakwiseCormen_():" << endl;
#endif

  //flag for optimization: whole protein or digest? (in this method not used any more)
  bool whole_protein;

  //only overall modifying necessary, digest_() already called by annotate()
  tryModifyOverall_();

  //get digest fragments (or whole protein) out of database
  vector<vector<int> > fragments = getFragments_(whole_protein);

  //parse partial modifications:
  ModificationStringParser mod_parser_(db_username_, db_password_, db_host_);
  partial_mods = mod_parser_.parse(partial_modification_string);

  //iterate fragments
  Modification mod_iter(1, db_username_, db_password_, db_host_); //create instance of \c Modification for iteration in \c ProteinDigest::getFragmentOverallModifiedMass
  for (vector<vector<int> >::iterator frag_it = fragments.begin(); frag_it != fragments.end(); frag_it++)
    {
      //fill variables needed for this method w. r. t. actual fragment
      fillCormenVariables_((*frag_it)[1], (*frag_it)[2], mod_iter);

      //! get mass of overall modified fragment
      hash_map<int,int> temp_overall_mods;
      double fragment_overall_mod_mass =
        protein_digest->getFragmentOverallModifiedMass((*frag_it)[1], (*frag_it)[2], masstype, temp_overall_mods, mod_iter);

      //! get mass of unmodified fragment
      double fragment_unmod_mass = 0;
      if (masstype == "average")
        {
          fragment_unmod_mass = protein_digest->getFragmentAverageMass((*frag_it)[1], (*frag_it)[2]);
        }
      else if (masstype == "mono")
        {
          fragment_unmod_mass = protein_digest->getFragmentMonoMass((*frag_it)[1], (*frag_it)[2]);
        }

      //iterate "real" peaks
      for (DPeakArray<1, DPeak<1> >::iterator it = peaklist_.begin(); it != peaklist_.end(); it++)
        {
          double mass_difference = (it->getPosition()[0] - fragment_overall_mod_mass);

          if (mass_difference >= 0)
            {
              list<pair<double, pair<list<int>, vector<int> > > > modification_combinations =
                exactSubsetSum_(mass_difference, range);

              if (modification_combinations.size() != 0)
                {
                  storeAnnotationsPeakwiseCormen_(modification_combinations, *it, *frag_it, fragment_overall_mod_mass,
                                                  fragment_unmod_mass, mod_iter, temp_overall_mods);
                }
            }
          else if (mass_difference >= -range)  //found only overall modified fragment, that has a whithin range larger mass than peak
            {
              list<pair<double, pair<list<int>, vector<int> > > > modification_combinations;
              modification_combinations.push_back(pair<double, pair<list<int>, vector<int> > >(0, pair<list<int>, vector<int> >()));
              storeAnnotationsPeakwiseCormen_(modification_combinations, *it, *frag_it, fragment_overall_mod_mass,
                                              fragment_unmod_mass, mod_iter, temp_overall_mods);
            }

        }
    }
}
