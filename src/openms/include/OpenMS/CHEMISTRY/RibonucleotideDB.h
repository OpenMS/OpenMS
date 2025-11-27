// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CHEMISTRY/Ribonucleotide.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <memory>
#include <unordered_map>

namespace OpenMS
{
  /** @ingroup Chemistry

      @brief Database of ribonucleotides (modified and unmodified)

      The information in this class comes primarily from the Modomics database (http://modomics.genesilico.pl/modifications/) and is read from a tab-separated text file in @p data/CHEMISTRY/Modomics.tsv.
      In addition, OpenMS-specific (as well as potentially user-supplied) modification definitions are read from the file @p data/CHEMISTRY/Custom_RNA_modifications.tsv.
  */
  class OPENMS_DLLAPI RibonucleotideDB final
  {
  public:
    using ConstRibonucleotidePtr = const Ribonucleotide *;

    /// const iterator type definition
    typedef std::vector<std::unique_ptr<Ribonucleotide>>::const_iterator ConstIterator;

    /// replacement for constructor (singleton pattern)
    static RibonucleotideDB* getInstance();

    /// destructor
    ~RibonucleotideDB() = default;

    /// copy constructor not available
    RibonucleotideDB(const RibonucleotideDB& other) = delete;

    /// assignment operator not available
    RibonucleotideDB& operator=(const RibonucleotideDB& other) = delete;

    /// Const iterator to beginning of database
    inline ConstIterator begin() const
    {
      return ribonucleotides_.begin();
    }

    /// Const iterator to end of database
    inline ConstIterator end() const
    {
      return ribonucleotides_.end();
    }

    /**
       @brief Get a ribonucleotide by its code (short name)

       @throw Exception::ElementNotFound if nothing was found
    */
    ConstRibonucleotidePtr getRibonucleotide(const std::string& code);

    /**
       @brief Get the ribonucleotide with the longest code that matches a prefix of @p seq

       @throw Exception::ElementNotFound if nothing was found
    */
    ConstRibonucleotidePtr getRibonucleotidePrefix(const std::string& seq);

    /**
       @brief Get the alternatives for an ambiguous modification code

       @throw Exception::ElementNotFound if nothing was found
    */
    std::pair<ConstRibonucleotidePtr, ConstRibonucleotidePtr> getRibonucleotideAlternatives(const std::string& code);


  protected:
    /// default constructor
    RibonucleotideDB();

    /// read (modified) nucleotides from input file
    void readFromFile_(const std::string& path);

    /// read from a newer version of Modomics that uses a JSON file
    void readFromJSON_(const std::string& path);

    /// create a (modified) nucleotide from an input row
    const std::unique_ptr<Ribonucleotide>  parseRow_(const std::string& row, Size line_count);

    /// list of known (modified) nucleotides
    std::vector<std::unique_ptr<Ribonucleotide>> ribonucleotides_;

    /// mapping of codes (short names) to indexes into @p ribonucleotides_
    std::unordered_map<std::string, Size> code_map_;

    /// mapping of ambiguity codes to the alternatives they represent
    std::map<std::string, std::pair<ConstRibonucleotidePtr, ConstRibonucleotidePtr>> ambiguity_map_;

    Size max_code_length_;
  };
}
