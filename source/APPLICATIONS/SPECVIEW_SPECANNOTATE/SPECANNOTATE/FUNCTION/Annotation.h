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


#ifndef __ANNOTATION_H__
#define __ANNOTATION_H__

#include <string>
#include <map>
#include <vector>
#include <iostream>

#include <ext/hash_map>


namespace OpenMS
  {

  //! class for storing one SINGLE annotation (a peak usually is annotated by SEVERAL annotations) \todo: equality operator
  class Annotation
    {
    private:

      //main stuff
      std::pair<std::pair<int, std::string>, std::pair<int, std::string> > fragment_;

      //modifications
      __gnu_cxx::hash_map<int, std::pair<std::pair<std::string, double>, std::pair<int, std::vector<int> > > > modifications_;


    public:

      //default constructor
      Annotation();

      //default destructor
      ~Annotation();

      //copy constructor
      Annotation(const Annotation& annotation);

      //output method
      void print(int no, std::ostream& out = std::cout);

      //main stuff
      void setFragment(int start_pos, std::string start_res, int end_pos, std::string end_res);
      std::string protein;
      std::string enzyme;

      //modifications
      void addModification(int ID, std::string name, double netto_plus_mass, int no_of_occurrences, std::vector<int> positions
                           = std::vector<int>());

      //masses
      double peak_mass;
      double calculated_annotation_mass;
      double unmodified_fragment_mass;
      double overall_modified_fragment_mass;
      double plus_mass_overall_modifications;
      double plus_mass_modification_combination;

      //annotation type specific parameters
      std::string annotation_method;
      std::string masstype;

      //database-related stuff
      int annotation_ID;
      int fragment_ID;
      int protein_ID;
      int first_real_mod_pless_ID;
      int first_real_mod_ID;

    };
}
#endif // __ANNOTATION_H__
