// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2016.
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
// $Maintainer: Lukas Zimmermann $
// $Authors: Lukas Zimmermann $
// --------------------------------------------------------------------------
#ifndef OPENMS_FORMAT_CROSSLINKCLASSESFILE_H
#define OPENMS_FORMAT_CROSSLINKCLASSESFILE_H

#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/METADATA/PeptideIdentification.h>


namespace OpenMS
{
  /**
      @brief This class provides methods to handle files for cross-link class
      specifications. It can load class definitions from a class specification file
      and collect MetaValues from PeptideIdentifications and add them to the
      present classes.

  @ingroup FileIO
  */
  class OPENMS_DLLAPI CrossLinkClassesFile
  {
  public:

    /// Default Constructor
    CrossLinkClassesFile();

    /// Destruktor
    ~CrossLinkClassesFile();

    /**
      Reads the classes specification file containing the xlink class specifications and parses content
      into the classes data structure.

     * @brief loads cross-link class specifications from a file
     * @param filename path to the file containing the cross-link class specification.
     * @return true if the file @p filename could be read successfully, false otherwise.
     */
    bool load(const String & filename);

    /**
     * Groups the MetaValue @p metavalue of the peptide identification @p pep_id by cross-link class.
     *
     * @brief Groups the MetaValue @p metavalue of the peptide identification @p pep_id by class
     * @param pep_id  Input PeptideIdentification
     * @param values  Where the extracted value should be collected into
     * @param metavalue Which meta value is to be extracted from the peptide identification @p pep_id.
     */
    template<typename MetaValueType >
    void collect(const PeptideIdentification & pep_id,
                                       std::map< String, std::vector< MetaValueType > > & values,
                                       const String & metavalue)
    {
      // Inspect all classes and determine whether pep_id belongs to this cross-link class
      for (std::map< String, std::vector< std::vector< StringList > > >::const_iterator it = this->classes.begin();
           it != this->classes.end(); ++it)\
      {
         bool class_fit = false;
         std::pair< String, std::vector< std::vector< StringList > > > xlink_class = *it;

         // Go throw all the clauses of the class and see if at least one clause matches
         for (std::vector< std::vector< StringList > >::const_iterator xlink_class_it = xlink_class.second.begin();
              xlink_class_it != xlink_class.second.end(); ++xlink_class_it)
         {
             std::vector< StringList > clause = *xlink_class_it;

             bool clause_matches = true;

             // Check all criteria within the class
             for (std::vector< StringList>::const_iterator clause_it = clause.begin(); clause_it != clause.end();
                  clause_it++)
             {
                 StringList attribute = *clause_it;
                 // Decide for the MetaInfoInterface to be investigated
                 MetaInfoInterface meta_info_interface;
                 if (attribute[0] == "PEPID")
                 {
                   meta_info_interface = pep_id;
                 }
                 else
                 {
                   std::vector<PeptideHit> const &  peptide_hits = pep_id.getHits();

                   if (attribute[0] == "ALPHA")   // Assume here that Alpha always exists
                   {
                       meta_info_interface = peptide_hits[0];
                   }
                   else if (peptide_hits.size() < 2)     // Must be BETA, but there is no beta
                   {
                       clause_matches = false;
                       break;
                   }
                   else
                   {
                       meta_info_interface = peptide_hits[1];
                   }
                  }
                 // Switch on the predicate
                 String predicate = attribute[1];     // What the meta value should be tested for
                 String meta_value  = attribute[2];   // The Meta value to be tested
                 bool meta_value_exists = meta_info_interface.metaValueExists(meta_value);

                 // Each one of the following criteria throws the peptide identification out of the class
                 if (    (predicate == "HAS"    && ! meta_value_exists)
                     ||  (predicate == "HASNOT" &&   meta_value_exists)
                     ||  (predicate == "IS"     && ! meta_value_exists)
                     ||  (predicate == "IS"     &&   ((String) meta_info_interface.getMetaValue(meta_value)) != attribute[3])
                     ||  (predicate == "ISNOT"  &&   ((String) meta_info_interface.getMetaValue(meta_value)) == attribute[3]))
                 {
                    clause_matches = false;
                    break;
                 }
             }
             if(clause_matches)   // This clause matches, so the peptide_identification belongs to this classd
             {
               class_fit = true;
               break;
             }  // else:  If the clause does not match, continue to the next clause
         }
         if (class_fit)
         {
             values[xlink_class.first].push_back((MetaValueType) pep_id.getMetaValue(metavalue));
         }
    }
    }

    /**
     * @brief has determines whether class with classname @p classname has been defined.
     * @param classname Name of cross-link class to be tested.
     * @return true if the cross-link class with classname @p classname has been defined, false otherwise.
     */
    inline bool has(const String & classname) const
    {
      return this->classes.find(classname) != this->classes.end();
    }

    typedef std::map< String, std::vector< std::vector < StringList > > >::const_iterator ClassesConstIterator;
    // Allows traversal of loaded class names
    class ClassNameConstIterator : public  ClassesConstIterator
    {
      public:
        ClassNameConstIterator() : ClassesConstIterator() {}
        ClassNameConstIterator(ClassesConstIterator s) : ClassesConstIterator(s){}
        String* operator->() { return (String* const)&(ClassesConstIterator::operator->()->first); }
        String operator*() { return ClassesConstIterator::operator*().first; }
    };

    inline ClassNameConstIterator begin() const { return this->classes.begin(); }
    inline ClassNameConstIterator end()   const { return this->classes.end(); }

  private:

    // TODO Initialize with default classes, so the classes file becomes optional
    // TOOD This structure is a bit awful
    // Loaded cross-link file classes
    std::map< String, std::vector< std::vector < StringList > > >  classes;
  };
} // namespace OpenMS

#endif // OPENMS_FORMAT_CROSSLINKCLASSESFILE_H
