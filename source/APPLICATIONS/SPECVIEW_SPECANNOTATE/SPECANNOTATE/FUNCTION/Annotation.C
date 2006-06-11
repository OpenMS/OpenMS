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
// $Id: Annotation.C,v 1.2 2006/03/28 08:03:27 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer:$
// --------------------------------------------------------------------------


#include "Annotation.h"

#include <iomanip>


using namespace std;
using namespace __gnu_cxx;
using namespace OpenMS;



Annotation::Annotation()
{}



Annotation::~Annotation()
{}



Annotation::Annotation(const Annotation& annotation)
{
  fragment_ = annotation.fragment_;
  modifications_ = annotation.modifications_;
  enzyme = annotation.enzyme;
  protein = annotation.protein;

  peak_mass = annotation.peak_mass;
  calculated_annotation_mass = annotation.calculated_annotation_mass;
  unmodified_fragment_mass = annotation.unmodified_fragment_mass;
  overall_modified_fragment_mass = annotation.overall_modified_fragment_mass;
  plus_mass_overall_modifications = annotation.plus_mass_overall_modifications;
  plus_mass_modification_combination = annotation.plus_mass_modification_combination;

  annotation_method = annotation.annotation_method;
  masstype = annotation.masstype;

  annotation_ID = annotation.annotation_ID;
  fragment_ID = annotation.fragment_ID;
  protein_ID = annotation.protein_ID;
  first_real_mod_pless_ID = annotation.first_real_mod_pless_ID;
  first_real_mod_ID = annotation.first_real_mod_ID;
}



void Annotation::setFragment(int start_pos, string start_res, int end_pos, string end_res)
{
  pair<int, string> start(start_pos, start_res);
  pair<int, string> end(end_pos, end_res);
  fragment_ = pair<pair<int, string>, pair<int, string> >(start, end);
}



void Annotation::addModification(int ID, string name, double netto_plus_mass, int no_of_occurrences, vector<int> positions)
{
  //if this modification has not been seen yet: insert it!
  if (modifications_.find(ID) == modifications_.end())
    {
      pair<string, double> name_and_mass(name, netto_plus_mass);
      pair<int, vector<int> > no_and_positions(no_of_occurrences, positions);
      modifications_[ID] = pair<pair<string, double>, pair<int, vector<int> > >(name_and_mass, no_and_positions);
    }
  //it has already been seen, add number of occurrences and positions
  else
    {
      modifications_[ID].second.first += no_of_occurrences;
      modifications_[ID].second.second.insert(modifications_[ID].second.second.begin(), positions.begin(), positions.end());
    }
}


void Annotation::print(int no, ostream& out)
{
  out.setf(ios_base::fixed,ios_base::floatfield)
    ; //set double - output to fixed decimal number of 6
  out << endl << endl;
  out << "#########################################################################################################" << endl;
  out << "ANNOTATION " << no << " (found using " << annotation_method << " method):" << endl;
  out << "-----------------------------------------------------" << endl;
  out << "(protein_ID: " << protein_ID << ", fragment_ID: " << fragment_ID;
  if (annotation_method == "peakwise_cormen")
    {
      out << ")" << endl;
    }
  else if (annotation_method == "improved_enumerate")
    {
      out << ", annotation_ID: " << annotation_ID << ", first realized_modification_positionless_ID: " << first_real_mod_pless_ID;
      out << ")" << endl;
    }
  else if (annotation_method == "enumerate")
    {
      out << ", annotation_ID: " << annotation_ID << ", first realized_modification_ID: " << first_real_mod_ID << ")" << endl;
    }
  out << endl;
  out << "MASSES (" << masstype << "):" << endl;
  out << "-----------------" << endl;
  out << setw(45) << "peak mass: "                      << setw(15) << peak_mass << endl;
  out << setw(45) << "calculated annotation mass: "     << setw(15) << calculated_annotation_mass << endl;
  out << setw(45) << "unmodified fragment mass: "       << setw(15) << unmodified_fragment_mass << endl;
  out << setw(45) << "overall modified fragment mass: " << setw(15) << overall_modified_fragment_mass << endl;
  out << setw(45) << "netto plus mass of overall modifications: " << setw(15) << plus_mass_overall_modifications << endl;
  out << setw(45) << "netto plus mass of modification combination: " << setw(15) << plus_mass_modification_combination << endl;
  out << endl;
  out << "FRAGMENT:" << endl;
  out << "---------" << endl;
  out << "residues from position " << fragment_.first.first << " (" << fragment_.first.second << ") to position ";
  out << fragment_.second.first << " (" << fragment_.second.second << "). (protein " << protein << ", digested with " << enzyme << ")";
  out << endl << endl;

  out << "MODIFICATIONS:" << endl;
  out << "--------------" <<endl;
  if (modifications_.empty())
    {
      out << "unmodified!" << endl;
    }
  else
    {
      for (hash_map<int, pair<pair<string, double>, pair<int, vector<int> > > >::iterator it = modifications_.begin();
           it != modifications_.end(); it++)
        {
          out << "ID " << setw(3) << it->first << ", netto mass " << setw(15) << it->second.first.second << ", occurring ";
          out << setw(2) << it->second.second.first << " times: " << setw(15) << it->second.first.first;
          if (annotation_method == "enumerate")
            {
              out << " at positions:" << endl << "\t";
              for (vector<int>::iterator intit = it->
                                                 second.second.second.begin();
                   intit != it->second.second.second.end();
                   intit++)
                {
                  out << *intit << ", ";
                }
              out << endl;
            }
          else
            {
              out << "." << endl;
            }
        }
    }
  out << "#########################################################################################################" << endl << endl << endl;
}
