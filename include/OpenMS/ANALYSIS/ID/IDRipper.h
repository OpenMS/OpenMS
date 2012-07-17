// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Immanuel Luhn $
// $Authors: Immanuel Luhn $
// --------------------------------------------------------------------------
#ifndef OPENMS_ANALYSIS_ID_IDRIPPER_H
#define OPENMS_ANALYSIS_ID_IDRIPPER_H

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/ProteinIdentification.h>


namespace OpenMS
{
  /**
    @brief Ripping protein/peptide identification according their file origin.

    Helper class, which is used by @ref TOPP_ProteinQuantifier. See there for further documentation.

    @htmlinclude OpenMS_IDRipper.parameters

    @ingroup Analysis_ID
  */
  class OPENMS_DLLAPI IDRipper : public DefaultParamHandler
  {
    public:

      /// Default constructor
      IDRipper();

      /// Desctructor
      virtual ~IDRipper();

      /**
        @brief Ripping protein/peptide identification according their file origin

        Iteration over all @p PeptideIdentification. For each annotated file origin create a map entry and store the
        respective @p PeptideIdentification and @p ProteinIdentification.

        @param ripped Contains the protein identification and peptide identification for each file origin annotated in proteins and peptides
        @param proteins Peptide identification annotated with file origin
        @param peptides Protein identification
      */
      void rip(std::map<String, std::pair< std::vector<ProteinIdentification>, std::vector<PeptideIdentification> > >& ripped, std::vector<ProteinIdentification>& proteins, std::vector<PeptideIdentification>& peptides);

    private:

      //Not implemented
      /// Copy constructor
      IDRipper(const IDRipper& rhs);

      // Not implemented
      /// Assignment
      IDRipper& operator= (const IDRipper& rhs);

      /// helper function, extracts all protein hits that match the protein accession
      void getProteinHits_(std::vector<ProteinHit>& result, const std::vector<ProteinHit>& protein_hits, const std::vector<String>& protein_accessions);
      /// helper function, returns the string representation of the peptide hit accession
      void getProteinAccessions_(std::vector<String>& result, const std::vector<PeptideHit>& peptide_hits);
      /// helper function, returns the protein identification for the given peptide identification based on the same identifier
      void getProteinIdentification_(ProteinIdentification& result, PeptideIdentification pep_ident, std::vector<ProteinIdentification>& prot_idents);
  };

} // namespace OpenMS

#endif // OPENMS_ANALYSIS_ID_IDRIPPER_H
