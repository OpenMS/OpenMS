// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CHEMISTRY/Residue.h>

namespace OpenMS
{

  // MS1 (modifications) and MS2 (fragments) for the different protocols
  namespace NuXLPresets
  {

  static constexpr std::array<const char*, 10> modifications_RNA_UV
  {
    "U:", 
    "U:-H2O", 

    "C:", 
    "C:-H2O",
    "C:-NH3",

    "G:", 
    "G:-H2O", 
    "G:-NH3", 

    "A:", 
    "A:-NH3"
  };

  static constexpr std::array<const char*, 29> modifications_RNA_UV_EXTENDED
  {
    "U:", 
    "U:-H2O", 
    "U:-HPO3",
    "U:-H3PO4", 
    "U:-HPO3-C2H3NO", //RBSID

    "C:", 
    "C:-NH3", 
    "C:-H2O", 
    "C:-HPO3",
    "C:-H3PO4",
    "C:-NH3-H2O", 
    "C:-NH3-HPO3", 
    "C:-NH3-H3PO4", 

    "G:", 
    "G:-H2O", 
    "G:-NH3", 
    "G:-HPO3",
    "G:-H3PO4", 
    "G:-NH3-H2O", 
    "G:-NH3-HPO3", 
    "G:-NH3-H3PO4", 

    "A:", 
    "A:-H2O", 
    "A:-HPO3",
    "A:-H3PO4", 
    "A:-NH3",
    "A:-NH3-H2O", 
    "A:-NH3-HPO3",
    "A:-NH3-H3PO4" 
  };

  static constexpr std::array<const char*, 6> modifications_RNA_UV_4SU
  {   
    "S:",
    "S:-H2O",
    "S:-H2S",
    "S:-C5H8O4",
    "S:-C5H8O4-H2O",
    "S:-C5H8O4-H2S"
  };

  static constexpr std::array<const char*, 14> modifications_RNA_UV_4SU_EXTENDED
  {   
    "S:",
    "S:-H2O",
    "S:-H3PO4",
    "S:-HPO3",
    "S:-H2S",
    "S:+HPO3-H2S", // iTRAPP
    "S:-HPO3-H2S", //pRBSID
    "S:+HPO3", //iTRAPP

    "S:-C5H8O4",
    "S:-C5H8O4-H2O",
    "S:-C5H8O4-H3PO4",
    "S:-C5H8O4-HPO3",
    "S:-C5H8O4-HPO3-H2S", ///pRBSID
    "S:-C5H8O4-H2S"
  };

  static constexpr std::array<const char*, 6> modifications_RNA_UV_6SG
  {   
    "X:",
    "X:-H2O",
    "X:-H2S",
    "X:-C5H8O4",
    "X:-C5H8O4-H2O",
    "X:-C5H8O4-H2S"
  };

  static constexpr std::array<const char*, 11> modifications_RNA_UV_6SG_EXTENDED
  {   
    "X:",
    "X:-H2O",
    "X:-H3PO4",
    "X:-HPO3",
    "X:-H2S",
    "X:-HPO3-H2S", //pRBSID
    

    "X:-C5H8O4",
    "X:-C5H8O4-H2O",
    "X:-C5H8O4-HPO3",
    "X:-C5H8O4-HPO3-H2S", ///pRBSID
    "X:-C5H8O4-H2S"    
  };

  static constexpr std::array<const char*, 4> modifications_RNA_NM
  {
    "U:+C5H9N1",
    "G:+C5H9N1",
    "C:+C5H9N1",
    "A:+C5H9N1"
  };

  static constexpr std::array<const char*, 34> modifications_RNA_NM_EXTENDED
  {
    "U:+C5H9N1",
    "U:+C5H9N1-HPO3",
    "U:+C5H9N1-H2O",
    "U:+C5H9N1-H3PO4",
    "U:+C5H9N1-H2O-H2O",
    "U:+C5H9N1-H3PO4-H2O",

    "G:+C5H9N1",
    "G:+C5H9N1-HPO3",
    "G:+C5H9N1-H2O",
    "G:+C5H9N1-H3PO4",
    "G:+C5H9N1-H2O-H2O",
    "G:+C5H9N1-H3PO4-H2O",
    "G:+C5H9N1-NH3",
    "G:+C5H9N1-NH3-H2O",
    "G:+C5H9N1-NH3-HPO3",
    "G:+C5H9N1-NH3-H3PO4",

    "C:+C5H9N1",
    "C:+C5H9N1-HPO3",
    "C:+C5H9N1-H2O",
    "C:+C5H9N1-H3PO4",
    "C:+C5H9N1-H2O-H2O",
    "C:+C5H9N1-H3PO4-H2O",
    "C:+C5H9N1-NH3",
    "C:+C5H9N1-NH3-H2O",
    "C:+C5H9N1-NH3-HPO3",
    "C:+C5H9N1-NH3-H3PO4",

    "A:+C5H9N1",
    "A:+C5H9N1-HPO3",
    "A:+C5H9N1-H2O",
    "A:+C5H9N1-H3PO4",
    "A:+C5H9N1-NH3",
    "A:+C5H9N1-NH3-H2O",
    "A:+C5H9N1-NH3-HPO3",
    "A:+C5H9N1-NH3-H3PO4"
  };

  static constexpr std::array<const char*, 4> modifications_RNA_DEB
  { 
    "U:+C4H6O2",
    "G:+C4H6O2",
    "C:+C4H6O2",
    "A:+C4H6O2"
  };

  static constexpr std::array<const char*, 34> modifications_RNA_DEB_EXTENDED
  {
    "U:+C4H6O2",
    "U:+C4H6O2-H2O",
    "U:+C4H6O2-HPO3",
    "U:+C4H6O2-H3PO4",
    "U:+C4H6O2-H2O-H2O",
    "U:+C4H6O2-H3PO4-H2O",

    "G:+C4H6O2",
    "G:+C4H6O2-H2O",
    "G:+C4H6O2-HPO3",
    "G:+C4H6O2-H3PO4",
    "G:+C4H6O2-H2O-H2O",
    "G:+C4H6O2-H3PO4-H2O",
    "G:+C4H6O2-NH3",
    "G:+C4H6O2-NH3-H2O",
    "G:+C4H602-NH3-HPO3",
    "G:+C4H6O2-NH3-H3PO4",

    "C:+C4H6O2",
    "C:+C4H6O2-H2O",
    "C:+C4H6O2-HPO3",
    "C:+C4H6O2-H3PO4",
    "C:+C4H6O2-H2O-H2O",
    "C:+C4H6O2-H3PO4-H2O",
    "C:+C4H6O2-NH3",
    "C:+C4H6O2-NH3-H2O",
    "C:+C4H602-NH3-HPO3",
    "C:+C4H6O2-NH3-H3PO4",

    "A:+C4H6O2",
    "A:+C4H6O2-H2O",
    "A:+C4H6O2-HPO3",
    "A:+C4H6O2-H3PO4",
    "A:+C4H6O2-NH3",
    "A:+C4H6O2-NH3-H2O",
    "A:+C4H6O2-NH3-HPO3",
    "A:+C4H6O2-NH3-H3PO4"
  };  
    
  static constexpr std::array<const char*, 12> modifications_RNA_FA
  {
    "G:+C",
    "G:+C-HPO3",
    "G:+C-H3PO4",
    "G:+C-H2O",
    
    "C:+C",
    "C:+C-H2O",
    "C:+C-HPO3",
    "C:+C-H3PO4",
    
    "A:+C",
    "A:+C-HPO3",
    "A:+C-H3PO4",
    "A:+C-H2O",
  };

  static constexpr std::array<const char*, 24> modifications_RNA_FA_EXTENDED
  {
    "G:+C",
    "G:+C-HPO3",
    "G:+C-H3PO4",
    "G:+C-H2O",
    "G:+C2",
    "G:+C2-HPO3",
    "G:+C2-H2O",
    "G:+C2-H3PO4",
    
    "C:+C",
    "C:+C-H2O",
    "C:+C-HPO3",
    "C:+C-H3PO4",
    "C:+C2",
    "C:+C2-HPO3",
    "C:+C2-H2O",
    "C:+C2-H3PO4",
    
    "A:+C",
    "A:+C-HPO3",
    "A:+C-H3PO4",
    "A:+C-H2O",
    "A:+C2",
    "A:+C2-HPO3",
    "A:+C2-H2O",
    "A:+C2-H3PO4"
  };
  static constexpr std::array<const char*, 16> modifications_DNA_FA
  {
    "G:+C",
    "G:+C-HPO3",
    "G:+C-H3PO4",
    "G:+C-H2O",
    
    "C:+C",
    "C:+C-H2O",
    "C:+C-HPO3",
    "C:+C-H3PO4",
    
    "A:+C",
    "A:+C-HPO3",
    "A:+C-H2O",
    "A:+C-H3PO4",

    "d:",
    "d:-H2O",
    "d:-H3PO4",
    "d:-HPO3"
  };

  static constexpr std::array<const char*, 28> modifications_DNA_FA_EXTENDED
  {
    "G:+C",
    "G:+C-HPO3",
    "G:+C-H3PO4",
    "G:+C-H2O",
    "G:+C2",
    "G:+C2-HPO3",
    "G:+C2-H3PO4",
    "G:+C2-H2O",
    
    "C:+C",
    "C:+C-H2O",
    "C:+C-HPO3",
    "C:+C-H3PO4",
    "C:+C2",
    "C:+C2-HPO3",
    "C:+C2-H3PO4",
    "C:+C2-H2O",
    
    "A:+C",
    "A:+C-HPO3",
    "A:+C-H2O",
    "A:+C-H3PO4",
    "A:+C2",
    "A:+C2-HPO3",
    "A:+C2-H3PO4",
    "A:+C2-H2O",

    "d:",
    "d:-H2O",
    "d:-H3PO4",
    "d:-HPO3"
  };

  static constexpr std::array<const char*, 12> modifications_DNA_UV
  {
    "T:",
    "T:-H2O",

    "G:",
    "G:-H2O",
    "G:-NH3",

    "A:",
    "A:-NH3",

    "C:",
    "C:-H2O",
    "C:-NH3",

    "d:", // deoxyribosephosphate (d in lower-letter = needs to be the cross-linked nt)
    "d:-H2O"
  };

  static constexpr std::array<const char*, 32> modifications_DNA_UV_EXTENDED
  {
    "T:",
    "T:-H2O",
    "T:-H3PO4",
    "T:-HPO3",

    "C:",
    "C:-H2O",
    "C:-H3PO4",
    "C:-HPO3",
    "C:-NH3",
    "C:-NH3-H2O",
    "C:-NH3-HPO3",
    "C:-NH3-H3PO4",

    "G:",
    "G:-H2O",
    "G:-H3PO4",
    "G:-HPO3",
    "G:-NH3",
    "G:-NH3-H2O",
    "G:-NH3-HPO3",
    "G:-NH3-H3PO4",

    "A:",
    "A:-NH3",
    "A:-H2O",
    "A:-H3PO4",
    "A:-HPO3",
    "A:-NH3-H2O",
    "A:-NH3-HPO3",
    "A:-NH3-H3PO4",

    "d:", // deoxyribosephosphate (d in lower-letter = needs to be the cross-linked nt)
    "d:-H2O",
    "d:-H3PO4",
    "d:-HPO3"
  };

  static constexpr std::array<const char*, 6> modifications_DNA_DEB
  {
    "T:+C4H6O2",
    "G:+C4H6O2",
    "C:+C4H6O2",
    "A:+C4H6O2",

    "d:", // deoxyribosephosphate (d in lower-letter = needs to be the cross-linked nt)
    "d:-H2O"
  };

  static constexpr std::array<const char*, 38> modifications_DNA_DEB_EXTENDED
  {
    "T:+C4H6O2",
    "T:+C4H6O2-H2O",
    "T:+C4H6O2-HPO3",
    "T:+C4H6O2-H3PO4",
    "T:+C4H6O2-H2O-H2O",
    "T:+C4H6O2-H3PO4-H2O",

    "G:+C4H6O2",
    "G:+C4H6O2-H2O",
    "G:+C4H6O2-HPO3",
    "G:+C4H6O2-H3PO4",
    "G:+C4H6O2-H2O-H2O",
    "G:+C4H6O2-H3PO4-H2O",
    "G:+C4H6O2-NH3",
    "G:+C4H6O2-NH3-H2O",
    "G:+C4H6O2-NH3-HPO3",
    "G:+C4H6O2-NH3-H3PO4",

    "C:+C4H6O2",
    "C:+C4H6O2-H2O",
    "C:+C4H6O2-HPO3",
    "C:+C4H6O2-H3PO4",
    "C:+C4H6O2-H2O-H2O",
    "C:+C4H6O2-H3PO4-H2O",
    "C:+C4H6O2-NH3-H2O",
    "C:+C4H6O2-NH3",
    "C:+C4H6O2-NH3-HPO3",
    "C:+C4H6O2-NH3-H3PO4",

    "A:+C4H6O2",
    "A:+C4H6O2-H2O",
    "A:+C4H6O2-H3PO4",
    "A:+C4H6O2-HPO3",
    "A:+C4H6O2-NH3",
    "A:+C4H6O2-NH3-H2O",
    "A:+C4H6O2-NH3-HPO3",
    "A:+C4H6O2-NH3-H3PO4",

    "d:", // deoxyribosephosphate (d in lower-letter = needs to be the cross-linked nt)
    "d:-H2O",
    "d:-H3PO4",
    "d:-HPO3"
  };

  static constexpr std::array<const char*, 6> modifications_DNA_NM
  {
    "T:+C5H9N1",
    "G:+C5H9N1",
    "C:+C5H9N1",
    "A:+C5H9N1",
    "d:",
    "d:-H2O"
  };

  static constexpr std::array<const char*, 38> modifications_DNA_NM_EXTENDED
  {
    "T:+C5H9N1",
    "T:+C5H9N1-H2O",
    "T:+C5H9N1-HPO3",
    "T:+C5H9N1-H3PO4",
    "T:+C5H9N1-H2O-H2O",
    "T:+C5H9N1-H3PO4-H2O",

    "G:+C5H9N1",
    "G:+C5H9N1-H2O",
    "G:+C5H9N1-HPO3",
    "G:+C5H9N1-H3PO4",
    "G:+C5H9N1-H2O-H2O",
    "G:+C5H9N1-H3PO4-H2O",
    "G:+C5H9N1-NH3",
    "G:+C5H9N1-NH3-H2O",
    "G:+C5H9N1-NH3-HPO3",
    "G:+C5H9N1-NH3-H3PO4",

    "C:+C5H9N1",
    "C:+C5H9N1-H2O",
    "C:+C5H9N1-HPO3",
    "C:+C5H9N1-H3PO4",
    "C:+C5H9N1-H2O-H2O",
    "C:+C5H9N1-H3PO4-H2O",
    "C:+C5H9N1-NH3",
    "C:+C5H9N1-NH3-HPO3",
    "C:+C5H9N1-NH3-H2O",
    "C:+C5H9N1-NH3-H3PO4",

    "A:+C5H9N1",
    "A:+C5H9N1-HPO3",
    "A:+C5H9N1-H2O",
    "A:+C5H9N1-H3PO4",      
    "A:+C5H9N1-NH3",
    "A:+C5H9N1-NH3-HPO3",
    "A:+C5H9N1-NH3-H2O",
    "A:+C5H9N1-NH3-H3PO4",

    "d:",
    "d:-H2O",
    "d:-H3PO4",
    "d:-HPO3"
  };


  static constexpr std::array<const char*, 14> modifications_DNA_BrU_UV
  {
    "T:",
    "T:-H2O",

    "B:",
    "B:-H2O",

    "G:",
    "G:-H2O",
    "G:-NH3",

    "A:",
    "A:-NH3",

    "C:",
    "C:-H2O",
    "C:-NH3",

    "d:", // deoxyribosephosphate (d in lower-letter = needs to be the cross-linked nt)
    "d:-H2O"
  };


  ///////////////////////////////////////////////////////////////////////////////////
  // fragment definitions
  ///////////////////////////////////////////////////////////////////////////////////

  // shared by default and Extended
  static constexpr std::array<const char*, 39> fragments_RNA_UV
  {
    "U:C3O;C3O",
    "U:C4H4N2O2;U'",
    "U:C4H2N2O1;U'-H2O",
    "U:C9H13N2O9P1;U",
    "U:C9H11N2O8P1;U-H2O",
    "U:C9H12N2O6;U-HPO3",
    "U:C9H10N2O5;U-H3PO4",

    "C:C4H5N3O;C'",
    "C:C4H3N3;C'-H2O",
    "C:C4H2N2O;C'-NH3",
    "C:C9H14N3O8P;C",
    "C:C9H11N2O8P;C-NH3",
    "C:C9H12N3O7P;C-H2O",
    "C:C9H9N2O7P;C-NH3-H2O",
    "C:C9H13N3O5;C-HPO3",
    "C:C9H11N3O4;C-H3PO4",
    "C:C9H10N2O5;C-NH3-HPO3",
    "C:C9H8N2O4;C-NH3-H3PO4",


    "G:C5H5N5O;G'",
    "G:C5H3N5;G'-H2O",
    "G:C5H2N4O;G'-NH3",
    "G:C10H14N5O8P;G",
    "G:C10H12N5O7P;G-H2O",
    "G:C10H11N4O8P;G-NH3",
    "G:C10H9N4O7P;G-NH3-H2O",
    "G:C10H13N5O5;G-HPO3",
    "G:C10H11N5O4;G-H3PO4",
    "G:C10H10N4O5;G-NH3-HPO3",
    "G:C10H8N4O4;G-NH3-H3PO4",

    "A:C5H5N5;A'",
    "A:C5H2N4;A'-NH3",
    "A:C10H14N5O7P;A",
    "A:C10H12N5O6P;A-H2O",
    "A:C10H11N4O7P;A-NH3",
    "A:C10H9N4O6P;A-NH3-H2O",
    "A:C10H13N5O4;A-HPO3",
    "A:C10H11N5O3;A-H3PO4",
    "A:C10H10N5O4;A-NH3-HPO3",
    "A:C10H8N5O3;A-NH3-H3PO4"
  };

  // shared by default and Extended
  static constexpr std::array<const char*, 41> fragments_RNA_UV_4SU
  {
      "S:C9H10N2O5;tU-H2S", // 4SU - H2S
      "S:C4H2N2O1;tU'-H2S", // 4SU base - H2S

      "U:C3O;C3O",
      "U:C4H4N2O2;U'",
      "U:C4H2N2O1;U'-H2O",
      "U:C9H13N2O9P1;U",
      "U:C9H11N2O8P1;U-H2O",
      "U:C9H12N2O6;U-HPO3",
      "U:C9H10N2O5;U-H3PO4",

      "C:C4H5N3O;C'",
      "C:C4H3N3;C'-H2O",
      "C:C4H2N2O;C'-NH3",
      "C:C9H14N3O8P;C",
      "C:C9H11N2O8P;C-NH3",
      "C:C9H12N3O7P;C-H2O",
      "C:C9H9N2O7P;C-NH3-H2O",
      "C:C9H13N3O5;C-HPO3",
      "C:C9H11N3O4;C-H3PO4",
      "C:C9H10N2O5;C-NH3-HPO3",
      "C:C9H8N2O4;C-NH3-H3PO4",

      "G:C5H5N5O;G'",
      "G:C5H3N5;G'-H2O",
      "G:C5H2N4O;G'-NH3",
      "G:C10H14N5O8P;G",
      "G:C10H12N5O7P;G-H2O",
      "G:C10H11N4O8P;G-NH3",
      "G:C10H9N4O7P;G-NH3-H2O",
      "G:C10H13N5O5;G-HPO3",
      "G:C10H11N5O4;G-H3PO4",
      "G:C10H10N4O5;G-NH3-HPO3",
      "G:C10H8N4O4;G-NH3-H3PO4",

      "A:C5H5N5;A'",
      "A:C5H2N4;A'-NH3",
      "A:C10H14N5O7P;A",
      "A:C10H12N5O6P;A-H2O",
      "A:C10H11N4O7P;A-NH3",
      "A:C10H9N4O6P;A-NH3-H2O",
      "A:C10H13N5O4;A-HPO3",
      "A:C10H11N5O3;A-H3PO4",
      "A:C10H10N5O4;A-NH3-HPO3",
      "A:C10H8N5O3;A-NH3-H3PO4"
  };      

  static constexpr std::array<const char*, 41> fragments_RNA_UV_6SG
  {
      "X:C10H11N5O4;tG-H2S",
      "X:C5H3N5;tG'-H2S",

      "U:C3O;C3O",
      "U:C4H4N2O2;U'",
      "U:C4H2N2O1;U'-H2O",
      "U:C9H13N2O9P1;U",
      "U:C9H11N2O8P1;U-H2O",
      "U:C9H12N2O6;U-HPO3",
      "U:C9H10N2O5;U-H3PO4",

      "C:C4H5N3O;C'",
      "C:C4H3N3;C'-H2O",
      "C:C4H2N2O;C'-NH3",
      "C:C9H14N3O8P;C",
      "C:C9H11N2O8P;C-NH3",
      "C:C9H12N3O7P;C-H2O",
      "C:C9H9N2O7P;C-NH3-H2O",
      "C:C9H13N3O5;C-HPO3",
      "C:C9H11N3O4;C-H3PO4",
      "C:C9H10N2O5;C-NH3-HPO3",
      "C:C9H8N2O4;C-NH3-H3PO4",

      "G:C5H5N5O;G'",
      "G:C5H3N5;G'-H2O",
      "G:C5H2N4O;G'-NH3",
      "G:C10H14N5O8P;G",
      "G:C10H12N5O7P;G-H2O",
      "G:C10H11N4O8P;G-NH3",
      "G:C10H9N4O7P;G-NH3-H2O",
      "G:C10H13N5O5;G-HPO3",
      "G:C10H11N5O4;G-H3PO4",
      "G:C10H10N4O5;G-NH3-HPO3",
      "G:C10H8N4O4;G-NH3-H3PO4",

      "A:C5H5N5;A'",
      "A:C5H2N4;A'-NH3",
      "A:C10H14N5O7P;A",
      "A:C10H12N5O6P;A-H2O",
      "A:C10H11N4O7P;A-NH3",
      "A:C10H9N4O6P;A-NH3-H2O",
      "A:C10H13N5O4;A-HPO3",
      "A:C10H11N5O3;A-H3PO4",
      "A:C10H10N5O4;A-NH3-HPO3",
      "A:C10H8N5O3;A-NH3-H3PO4"
  };     

  // shared by default and Extended
  static constexpr std::array<const char*, 42> fragments_DNA_UV
  {
    "T:C5H6N2O2;T'",
    "T:C5H4N2O;T'-H2O",
    "T:C10H15N2O8P;T",
    "T:C10H13N2O7P;T-H2O",
    "T:C10H14N2O5;T-HPO3",
    "T:C10H12N2O4;T-H3PO4",

    "C:C9H14N3O7P;C",
    "C:C9H11N2O7P;C-NH3",
    "C:C9H12N3O6P;C-H2O",
    "C:C9H9N2O6P;C-NH3-H2O",
    "C:C9H13N3O4;C-HPO3",
    "C:C9H11N3O3;C-H3PO4",
    "C:C9H10N2O4;C-NH3-HPO3",
    "C:C9H8N2O3;C-NH3-H3PO4",
    "C:C4H5N3O;C'",
    "C:C4H3N3;C'-H2O",
    "C:C4H2N2O;C'-NH3",

    "G:C10H14N5O7P;G",
    "G:C10H12N5O6P;G-H2O",
    "G:C10H11N4O7P;G-NH3",
    "G:C10H9N4O6P;G-NH3-H2O",
    "G:C10H13N5O4;G-HPO3",
    "G:C10H10N4O4;G-NH3-HPO3",
    "G:C10H11N5O3;G-H3PO4",
    "G:C10H8N4O3;G-NH3-H3PO4",
    "G:C5H5N5O;G'",
    "G:C5H3N5;G'-H2O",
    "G:C5H2N4O;G'-NH3",

    "A:C10H14N5O6P;A",
    "A:C10H12N5O5P;A-H2O",
    "A:C10H11N4O6P;A-NH3",
    "A:C10H9N4O5P;A-NH3-H2O",
    "A:C10H13N5O3;A-HPO3",
    "A:C10H11N5O2;A-H3PO4",
    "A:C10H10N5O3;A-NH3-HPO3",
    "A:C10H8N5O2;A-NH3-H3PO4",
    "A:C5H5N5;A'",
    "A:C5H2N4;A'-NH3",

    // base was lost -> only dribose = C5H9O6P remains
    "d:C5H9O6P;C5H9O6P", 
    "d:C5H7O5P;C5H9O6P-H2O",        
    "d:C5H8O3;C5H9O6P-HPO3",
    "d:C5H6O2;C5H9O6P-H3PO4"
  };

  // shared by default and Extended
  static constexpr std::array<const char*, 48> fragments_DNA_BrU_UV
  {
    "T:C5H6N2O2;T'",
    "T:C5H4N2O;T'-H2O",
    "T:C10H15N2O8P;T",
    "T:C10H13N2O7P;T-H2O",
    "T:C10H14N2O5;T-HPO3",
    "T:C10H12N2O4;T-H3PO4",

    "B:C4H2N2O2;B'",
    "B:C4N2O;B'-H2O",
    "B:C9H11N2O8P;B",
    "B:C9H9N2O7P;B-H2O",
    "B:C9H10N2O5;B-HPO3",
    "B:C9H8N2O4;B-H3PO4",

    "C:C9H14N3O7P;C",
    "C:C9H11N2O7P;C-NH3",
    "C:C9H12N3O6P;C-H2O",
    "C:C9H9N2O6P;C-NH3-H2O",
    "C:C9H13N3O4;C-HPO3",
    "C:C9H11N3O3;C-H3PO4",
    "C:C9H10N2O4;C-NH3-HPO3",
    "C:C9H8N2O3;C-NH3-H3PO4",
    "C:C4H5N3O;C'",
    "C:C4H3N3;C'-H2O",
    "C:C4H2N2O;C'-NH3",

    "G:C10H14N5O7P;G",
    "G:C10H12N5O6P;G-H2O",
    "G:C10H11N4O7P;G-NH3",
    "G:C10H9N4O6P;G-NH3-H2O",
    "G:C10H13N5O4;G-HPO3",
    "G:C10H10N4O4;G-NH3-HPO3",
    "G:C10H11N5O3;G-H3PO4",
    "G:C10H8N4O3;G-NH3-H3PO4",
    "G:C5H5N5O;G'",
    "G:C5H3N5;G'-H2O",
    "G:C5H2N4O;G'-NH3",

    "A:C10H14N5O6P;A",
    "A:C10H12N5O5P;A-H2O",
    "A:C10H11N4O6P;A-NH3",
    "A:C10H9N4O5P;A-NH3-H2O",
    "A:C10H13N5O3;A-HPO3",
    "A:C10H11N5O2;A-H3PO4",
    "A:C10H10N5O3;A-NH3-HPO3",
    "A:C10H8N5O2;A-NH3-H3PO4",
    "A:C5H5N5;A'",
    "A:C5H2N4;A'-NH3",

    // base was lost -> only dribose = C5H9O6P remains
    "d:C5H9O6P;C5H9O6P", 
    "d:C5H7O5P;C5H9O6P-H2O",        
    "d:C5H8O3;C5H9O6P-HPO3",
    "d:C5H6O2;C5H9O6P-H3PO4"
  };

  // shared by default and Extended
  static constexpr std::array<const char*, 47> fragments_RNA_DEB
  {
    "U:C4H6O2;DEB",
    "U:C4H4O;DEB-H2O",
    "U:C7H6O3;DEB+C3O",
    "U:C8H10N2O4;DEB+U'",
    "U:C8H8N2O3;DEB+U'-H2O",        
    "U:C13H19N2O11P1;DEB+U",
    "U:C13H17N2O10P1;DEB+U-H2O",
    "U:C13H18N2O8;DEB+U-HPO3",
    "U:C13H16N2O7;DEB+U-H3PO4",

    "G:C4H6O2;DEB",
    "G:C4H4O;DEB-H2O",
    "G:C9H11N5O3;DEB+G'",
    "G:C8H9N5O3;DEB+G'-H2O",
    "G:C9H8N4O3;DEB+G'-NH3",
    "G:C14H20N5O10P1;DEB+G",
    "G:C14H18N5O9P1;DEB+G-H2O",
    "G:C14H17N4O10P1;DEB+G-NH3",
    "G:C14H19N5O7;DEB+G-HPO3",
    "G:C14H17N5O6;DEB+G-H3PO4",
    "G:C14H16N4O7;DEB+G-NH3-HPO3",
    "G:C14H14N4O6;DEB+G-NH3-H3PO4",
    "G:C14H15N4O9P1;DEB+G-NH3-H2O",

    "C:C4H6O2;DEB",
    "C:C4H4O;DEB-H2O",
    "C:C8H11N3O3;DEB+C'",		
    "C:C8H9N3O2;DEB+C'-H2O",
    "C:C8H8N2O3;DEB+C'-NH3",
    "C:C13H20N3O10P1;DEB+C",
    "C:C13H18N3O9P1;DEB+C-H2O",
    "C:C13H17N2O10P1;DEB+C-NH3",
    "C:C13H19N3O7;DEB+C-HPO3",
    "C:C13H17N3O6;DEB+C-H3PO4",
    "C:C13H16N2O7;DEB+C-NH3-HPO3",
    "C:C13H14N2O6;DEB+C-NH3-H3PO4",
    "C:C13H15N2O9P1;DEB+C-NH3-H2O",
    
    "A:C4H6O2;DEB",
    "A:C4H4O;DEB-H2O",
    "A:C9H11N5O2;DEB+A'",
    "A:C9H17N4O;DEB+A'-NH3",
    "A:C14H20N5O9P1;DEB+A",
    "A:C14H18N5O8P1;DEB+A-H2O",
    "A:C14H17N4O9P1;DEB+A-NH3",
    "A:C14H19N5O6;DEB+A-HPO3",
    "A:C14H17N5O5;DEB+A-H3PO4",
    "A:C14H16N4O6;DEB+A-NH3-HPO3",
    "A:C14H14N4O5;DEB+A-NH3-H3PO4",
    "A:C14H15N4O8P1;DEB+A-NH3-H2O"  
  };

  // shared by default and Extended
  static constexpr std::array<const char*, 50> fragments_DNA_DEB
  { 
    "T:C4H6O2;DEB",
    "T:C4H4O;DEB-H2O",
    "T:C9H12N2O4;DEB+T'",
    "T:C9H10N2O3;DEB+T'-H2O",
    "T:C14H21N2O10P1;DEB+T",   
    "T:C14H19N2O9P1;DEB+T-H2O", 
    "T:C14H20N2O7;DEB+T-HPO3",  
    "T:C14H18N2O6;DEB+T-H3PO4", 

    "C:C4H6O2;DEB",
    "C:C4H4O;DEB-H2O",
    "C:C8H11N3O3;DEB+C'",
    "C:C8H8N2O3;DEB+C'-NH3",
    "C:C8H9N3O2;DEB+C'-H2O",
    "C:C13H20N3O9P1;DEB+C",
    "C:C13H17N2O9P1;DEB+C-NH3",
    "C:C13H18N3O8P1;DEB+C-H2O",
    "C:C13H19N3O6;DEB+C-HPO3", 
    "C:C13H17N3O5;DEB+C-H3PO4",
    "C:C13H16N2O6;DEB+C-NH3-HPO3", 
    "C:C13H14N2O5;DEB+C-NH3-H3PO4",
    "C:C13H15N2O8P1;DEB+C-NH3-H2O",

    "G:C4H6O2;DEB",
    "G:C4H4O;DEB-H2O",
    "G:C9H11N5O3;DEB+G'", 
    "G:C9H8N4O3;DEB+G'-NH3", 
    "G:C9H9N5O2;DEB+G'-H2O", 
    "G:C14H20N5O9P1;DEB+G",
    "G:C14H17N4O9P1;DEB+G-NH3",
    "G:C14H18N5O8P1;DEB+G-H2O",
    "G:C14H19N5O6;DEB+G-HPO3",
    "G:C14H17N5O5;DEB+G-H3PO4",
    "G:C14H16N4O6;DEB+G-NH3-HPO3",
    "G:C14H14N4O5;DEB+G-NH3-H3PO4",
    "G:C14H15N4O8P1;DEB+G-NH3-H2O",

    "A:C4H6O2;DEB",
    "A:C4H4O;DEB-H2O",          
    "A:C9H11N5O2;DEB+A'",  
    "A:C9H8N4O2;DEB+A'-NH3",   
    "A:C14H20N5O8P1;DEB+A",   
    "A:C14H17N4O8P1;DEB+A-NH3",     
    "A:C14H18N5O7P1;DEB+A-H2O", 
    "A:C14H19N5O5;DEB+A-HPO3",  
    "A:C14H17N5O4;DEB+A-H3PO4", 
    "A:C14H16N4O5;DEB+A-NH3-HPO3",  
    "A:C14H14N4O4;DEB+A-NH3-H3PO4",
    "A:C14H15N4O7P1;DEB+A-NH3-H2O", 

    "d:C5H9O6P;C5H9O6P", 
    "d:C5H7O5P;C5H9O6P-H2O",        
    "d:C5H8O3;C5H9O6P-HPO3",
    "d:C5H6O2;C5H9O6P-H3PO4"
  };

  // shared by default and Extended
  static constexpr std::array<const char*, 43> fragments_RNA_NM
  {
    "U:C5H9N1;NM",
    "U:C8H9N1O1;NM+C3O",   
    "U:C14H22N3O9P1;NM+U",
    "U:C14H20N3O8P1;NM+U-H2O",
    "U:C14H21N3O6;NM+U-HPO3",
    "U:C14H19N3O5;NM+U-H3PO4",
    "U:C9H13N3O2;NM+U'",
    "U:C9H11N3O1;NM+U'-H2O", 

    "C:C5H9N1;NM",             
    "C:C14H23N4O8P1;NM+C",
    "C:C14H21N4O7P1;NM+C-H2O",
    "C:C14H20N3O8P1;NM+C-NH3",
    "C:C14H22N4O5;NM+C-HPO3",
    "C:C14H20N4O4;NM+C-H3PO4",   
    "C:C14H19N3O5;NM+C-NH3-HPO3",
    "C:C14H17N3O4;NM+C-NH3-H3PO4",
    "C:C9H14N4O1;NM+C'",
    "C:C9H12N4;NM+C'-H2O",
    "C:C9H11N3O1;NM+C'-NH3",
    "C:C14H18N3O7P1;NM+C-NH3-H2O",

    "G:C5H9N1;NM",             
    "G:C15H23N6O8P1;NM+G",
    "G:C15H21N6O7P1;NM+G-H2O",
    "G:C15H20N5O8P1;NM+G-NH3",
    "G:C15H22N6O5;NM+G-HPO3",
    "G:C15H20N6O4;NM+G-H3PO4",
    "G:C15H19N5O5;NM+G-NH3-HPO3",
    "G:C15H17N5O4;NM+G-NH3-H3PO4",
    "G:C10H14N6O1;NM+G'",
    "G:C10H12N6;NM+G'-H2O",
    "G:C10H14N6O1;NM+G'-NH3",
    "G:C15H18N5O7P1;NM+G-NH3-H2O",

    "A:C5H9N1;NM",        
    "A:C15H23N6O7P1;NM+A",
    "A:C15H21N6O6P1;NM+A-H2O",
    "A:C15H20N5O7P1;NM+A-NH3",
    "A:C15H22N6O4;NM+A-HPO3",
    "A:C15H20N6O3;NM+A-H3PO4", 
    "A:C15H19N5O4;NM+A-NH3-HPO3",
    "A:C15H17N5O3;NM+A-NH3-H3PO4",  		
    "A:C10H14N6;NM+A'",
    "A:C10H11N5;NM+A'-NH3",
    "A:C15H18N5O6P1;NM+A-NH3-H2O"
  };

  static constexpr std::array<const char*, 46> fragments_RNA_FA
  {
    "U:C3O;C3O",
    "U:C4H4N2O2;U'",
    "U:C4H2N2O1;U'-H2O",
    "U:C9H13N2O9P1;U",
    "U:C9H11N2O8P1;U-H2O",
    "U:C9H12N2O6;U-HPO3",
    "U:C9H10N2O5;U-H3PO4",

    "G:C;FA",
    "G:C5H5N5O;G'",
    "G:C10H14N5O8P;G",
    "G:C6H5N5O;G'+FA",
    "G:C11H14N5O8P;G+FA",
    "G:C11H12N5O7P;G+FA-H2O",
    "G:C11H11N5O4;G+FA-H3PO4",
    "G:C11H13N5O5;G+FA-HPO3",
    "G:C7H5N5O;G'+2FA",
    "G:C12H14N5O8P;G+2FA",
    "G:C12H12N5O7P;G+2FA-H2O",
    "G:C12H11N5O4;G+2FA-H3PO4",
    "G:C12H13N5O5;G+2FA-HPO3",
    
    "C:C;FA",
    "C:C4H5N3O;C'",
    "C:C9H14N3O8P;C",
    "C:C5H5N3O;C'+FA",
    "C:C10H14N3O8P;C+FA",
    "C:C10H13N3O5;C+FA-HPO3",
    "C:C10H12N3O7P;C+FA-H2O",
    "C:C10H11N3O4;C+FA-H3PO4",
    "C:C6H5N3O;C'+2FA",
    "C:C11H14N3O8P;C+2FA",
    "C:C11H13N3O5;C+2FA-HPO3",
    "C:C11H12N3O7P;C+2FA-H2O",
    "C:C11H11N3O4;C+2FA-H3PO4",
    
    "A:C;FA",
    "A:C5H5N5;A'",
    "A:C10H14N5O7P;A",
    "A:C6H5N5;A'+FA",
    "A:C11H14N5O7P;A+FA",
    "A:C11H13N5O4;A+FA-HPO3",
    "A:C11H12N5O6P;A+FA-H2O",
    "A:C11H11N5O3;A+FA-H3PO4",
    "A:C7H5N5;A'+2FA",
    "A:C12H14N5O7P;A+2FA",
    "A:C12H13N5O4;A+2FA-HPO3",
    "A:C12H12N5O6P;A+2FA-H2O",
    "A:C12H11N5O3;A+2FA-H3PO4"
  };

  static constexpr std::array<const char*, 49> fragments_DNA_FA
  {
    "T:C5H6N2O2;T'",
    "T:C5H4N2O1;T'-H2O",
    "T:C10H15N2O8P1;T",
    "T:C10H13N2O7P1;T-H2O",
    "T:C10H12N2O4;T-HPO3",
    "T:C10H10N2O3;T-H3PO4",

    "G:C;FA",
    "G:C5H5N5O;G'",
    "G:C10H14N5O7P;G",
    "G:C6H5N5O;G'+FA",
    "G:C11H14N5O7P;G+FA",
    "G:C11H13N5O4;G+FA-HPO3",
    "G:C11H12N5O6P;G+FA-H2O",
    "G:C11H11N5O3;G+FA-H3PO4",
    "G:C7H5N5O;G'+2FA",
    "G:C12H14N5O7P;G+2FA",
    "G:C12H13N5O4;G+2FA-HPO3",
    "G:C12H12N5O6P;G+2FA-H2O",
    "G:C12H11N5O3;G+2FA-H3PO4",
    
    "C:C;FA",
    "C:C4H5N3O;C'",
    "C:C9H14N3O7P;C",
    "C:C5H5N3O;C'+FA",
    "C:C10H14N3O7P;C+FA",
    "C:C10H13N3O4;C+FA-HPO3",
    "C:C10H12N3O6P;C+FA-H2O",
    "C:C10H11N3O3;C+FA-H3PO4",
    "C:C6H5N3O;C'+2FA",
    "C:C11H14N3O7P;C+2FA",
    "C:C11H13N3O4;C+2FA-HPO3",
    "C:C11H12N3O6P;C+2FA-H2O",
    "C:C11H11N3O3;C+2FA-H3PO4",
    
    "A:C;FA",
    "A:C5H5N5;A'",
    "A:C10H14N5O6P;A",
    "A:C6H5N5;A'+FA",
    "A:C11H14N5O6P;A+FA",
    "A:C11H13N5O3;A+FA-HPO3",
    "A:C11H12N5O5P;A+FA-H2O",
    "A:C11H11N5O2;A+FA-H3PO4",
    "A:C7H5N5;A'+2FA",
    "A:C12H14N5O6P;A+2FA",
    "A:C12H13N5O3;A+2FA-HPO3",
    "A:C12H12N5O5P;A+2FA-H2O",
    "A:C12H11N5O2;A+2FA-H3PO4",

    "d:C5H9O6P;C5H9O6P", 
    "d:C5H7O5P;C5H9O6P-H2O",        
    "d:C5H8O3;C5H9O6P-HPO3",
    "d:C5H6O2;C5H9O6P-H3PO4"
  };

  // shared by default and Extended
  static constexpr std::array<const char*, 46> fragments_DNA_NM
  { 
    "T:C5H9N1;NM",      
    "T:C15H24N3O8P1;NM+T",
    "T:C15H22N3O7P1;NM+T-H2O",
    "T:C15H23N3O5;NM+T-HPO3",
    "T:C15H21N3O4;NM+T-H3PO4",  
    "T:C10H15N3O2;NM+T'",
    "T:C10H13N3O1;NM+T'-H2O",

    "C:C5H9N;NM",
    "C:C9H14N4O1;NM+C'",
    "C:C9H11N3O1;NM+C'-NH3",
    "C:C9H12N4;NM+C'-H2O",
    "C:C14H23N4O7P1;NM+C",
    "C:C14H21N4O6P1;NM+C-H2O",
    "C:C14H20N3O7P1;NM+C-NH3",
    "C:C14H18N3O6P1;NM+C-NH3-H2O",
    "C:C14H20N4O3;NM+C-H3PO4",
    "C:C14H22N4O4;NM+C-HPO3",
    "C:C14H19N3O4;NM+C-NH3-HPO3",
    "C:C14H17N3O3;NM+C-NH3-H3PO4",

    "G:C5H9N1;NM",
    "G:C10H14N6O1;NM+G'",
    "G:C10H12N6;NM+G'-H2O",
    "G:C10H11N5O1;NM+G'-NH3",
    "G:C15H23N6O7P1;NM+G",
    "G:C15H21N6O6P1;NM+G-H2O",
    "G:C15H22N6O4;NM+G-HPO3",
    "G:C15H20N6O3;NM+G-H3PO4",
    "G:C15H20N5O7P1;NM+G-NH3",
    "G:C15H18N5O6P1;NM+G-NH3-H2O",
    "G:C15H19N5O4;NM+G-NH3-HPO3",
    "G:C15H17N5O3;NM+G-NH3-H3PO4",

    "A:C5H9N1;NM",
    "A:C10H14N6;NM+A'",
    "A:C10H11N5;NM+A'-NH3",
    "A:C15H23N6O6P1;NM+A",
    "A:C15H20N6O2;NM+A-H3PO4",
    "A:C15H21N6O5P1;NM+A-H2O",
    "A:C15H22N6O3;NM+A-HPO3",
    "A:C15H20N5O6P1;NM+A-NH3",
    "A:C15H18N5O5P1;NM+A-NH3-H2O",
    "A:C15H19N5O3;NM+A-NH3-HPO3",
    "A:C15H17N5O2;NM+A-NH3-H3PO4",

    "d:C5H9O6P;C5H9O6P", 
    "d:C5H7O5P;C5H9O6P-H2O",        
    "d:C5H8O3;C5H9O6P-HPO3",
    "d:C5H6O2;C5H9O6P-H3PO4"
  };

    // the nucleotides (=mono-phosphates) and the deoxyribosephosphate (for DNA)
    static constexpr std::array<const char*, 5> DNA_nucleotides {"A=C10H14N5O6P", "C=C9H14N3O7P", "G=C10H14N5O7P", "T=C10H15N2O8P", "d=C5H9O6P"};
    static constexpr std::array<const char*, 4> RNA_nucleotides {"A=C10H14N5O7P", "C=C9H14N3O8P", "G=C10H14N5O8P", "U=C9H13N2O9P"}; 
    static constexpr std::array<const char*, 5> DNA_mapping {"A->A", "C->C", "G->G", "T->T", "d->d"};
    static constexpr std::array<const char*, 4> RNA_mapping {"A->A", "C->C", "G->G", "U->U"};

    static constexpr std::array<const char*, 24> presets_names {
      "none", 
      "RNA-UV (U)", 
      "RNA-UV (UCGA)",
      "RNA-UV Extended (U)", 
      "RNA-UV Extended (UCGA)", 
      "RNA-UV (4SU)", 
      "RNA-UV Extended (4SU)",
      "RNA-UV (6SG)", 
      "RNA-UV Extended (6SG)",
      "RNA-DEB",
      "RNA-DEB Extended", 
      "RNA-NM",
      "RNA-NM Extended", 
      "DNA-UV", 
      "DNA-UV Extended", 
      "DNA-DEB", 
      "DNA-DEB Extended",
      "DNA-NM",
      "DNA-NM Extended",
      "RNA-FA",
      "RNA-FA Extended",
      "DNA-FA",
      "DNA-FA Extended",
      "DNA-UV (BrU)"
   };

  void getPresets(const String& p, 
    StringList& nucleotides, 
    StringList& mapping, 
    StringList& modifications, 
    StringList& fragment_adducts, 
    String& can_cross_link)
  {
    // construct name list from constexpr array
    const StringList names(presets_names.begin(), presets_names.end());

    // sanity check: preset name needs to be in the list of supported presets
    if (auto it = find(names.begin(), names.end(), p); it == names.end())
    {
      throw std::runtime_error("Error: unknown preset.");
    }    

    // set NTs for RNA / DNA
    if (p.hasPrefix("RNA"))
    {
      nucleotides = StringList(RNA_nucleotides.begin(), RNA_nucleotides.end());
      mapping = StringList(RNA_mapping.begin(), RNA_mapping.end());
    }
    else if (p.hasPrefix("DNA"))
    {
      nucleotides = StringList(DNA_nucleotides.begin(), DNA_nucleotides.end());
      mapping = StringList(DNA_mapping.begin(), DNA_mapping.end());
    }

    // initialize all StringLists from constexpr arrays
    // note: we do this here as this raises a logic error if e.g., size of the array doesn't match the reserved size.
    //       This can easily happen if a comma is omitted and two string literals on two lines joined
    StringList RNA_UV_modifications(modifications_RNA_UV.begin(), modifications_RNA_UV.end());
    StringList RNA_UV_EXTENDED_modifications(modifications_RNA_UV_EXTENDED.begin(), modifications_RNA_UV_EXTENDED.end());
    StringList RNA_UV_fragments(fragments_RNA_UV.begin(), fragments_RNA_UV.end());

    StringList RNA_UV_4SU_modifications(modifications_RNA_UV_4SU.begin(), modifications_RNA_UV_4SU.end());
    StringList RNA_UV_4SU_EXTENDED_modifications(modifications_RNA_UV_4SU_EXTENDED.begin(), modifications_RNA_UV_4SU_EXTENDED.end());
    StringList RNA_UV_4SU_fragments(fragments_RNA_UV_4SU.begin(), fragments_RNA_UV_4SU.end());

    StringList RNA_UV_6SG_modifications(modifications_RNA_UV_6SG.begin(), modifications_RNA_UV_6SG.end());
    StringList RNA_UV_6SG_EXTENDED_modifications(modifications_RNA_UV_6SG_EXTENDED.begin(), modifications_RNA_UV_6SG_EXTENDED.end());
    StringList RNA_UV_6SG_fragments(fragments_RNA_UV_6SG.begin(), fragments_RNA_UV_6SG.end());

    StringList DNA_UV_modifications(modifications_DNA_UV.begin(), modifications_DNA_UV.end());
    StringList DNA_UV_EXTENDED_modifications(modifications_DNA_UV_EXTENDED.begin(), modifications_DNA_UV_EXTENDED.end());
    StringList DNA_UV_fragments(fragments_DNA_UV.begin(), fragments_DNA_UV.end());

    StringList RNA_DEB_modifications(modifications_RNA_DEB.begin(), modifications_RNA_DEB.end());
    StringList RNA_DEB_EXTENDED_modifications(modifications_RNA_DEB_EXTENDED.begin(), modifications_RNA_DEB_EXTENDED.end());
    StringList RNA_DEB_fragments(fragments_RNA_DEB.begin(), fragments_RNA_DEB.end());

    StringList RNA_NM_modifications(modifications_RNA_NM.begin(), modifications_RNA_NM.end());
    StringList RNA_NM_EXTENDED_modifications(modifications_RNA_NM_EXTENDED.begin(), modifications_RNA_NM_EXTENDED.end());
    StringList RNA_NM_fragments(fragments_RNA_NM.begin(), fragments_RNA_NM.end()); 

    StringList DNA_DEB_modifications(modifications_DNA_DEB.begin(), modifications_DNA_DEB.end());
    StringList DNA_DEB_EXTENDED_modifications(modifications_DNA_DEB_EXTENDED.begin(), modifications_DNA_DEB_EXTENDED.end());
    StringList DNA_DEB_fragments(fragments_DNA_DEB.begin(), fragments_DNA_DEB.end());

    StringList DNA_NM_modifications(modifications_DNA_NM.begin(), modifications_DNA_NM.end());
    StringList DNA_NM_EXTENDED_modifications(modifications_DNA_NM_EXTENDED.begin(), modifications_DNA_NM_EXTENDED.end());
    StringList DNA_NM_fragments(fragments_DNA_NM.begin(), fragments_DNA_NM.end());
    
    StringList RNA_FA_modifications(modifications_RNA_FA.begin(), modifications_RNA_FA.end());
    StringList RNA_FA_fragments(fragments_RNA_FA.begin(), fragments_RNA_FA.end());
    StringList RNA_FA_EXTENDED_modifications(modifications_RNA_FA_EXTENDED.begin(), modifications_RNA_FA_EXTENDED.end());

    StringList DNA_FA_modifications(modifications_DNA_FA.begin(), modifications_DNA_FA.end());
    StringList DNA_FA_fragments(fragments_DNA_FA.begin(), fragments_DNA_FA.end());
    StringList DNA_FA_EXTENDED_modifications(modifications_DNA_FA_EXTENDED.begin(), modifications_DNA_FA_EXTENDED.end());
   
    StringList DNA_UV_BrU_modifications(modifications_DNA_BrU_UV.begin(), modifications_DNA_BrU_UV.end());
    StringList DNA_UV_BrU_fragments(fragments_DNA_BrU_UV.begin(), fragments_DNA_BrU_UV.end());

    const String RNA_U = "U";
    const String RNA_UCGA = "UCGA";
    const String DNA_TCGAd = "TCGAd";
    const String RNA_CGA = "CGA";
    const String DNA_CGAd = "CGAd";

    // set precursor + fragment adducts and cross-linked nucleotide
    if (p == "RNA-UV (U)" || p  == "RNA-UV (UCGA)")
    {
      modifications = RNA_UV_modifications;
      fragment_adducts = RNA_UV_fragments;
      can_cross_link = (p == "RNA-UV (U)") ? RNA_U : RNA_UCGA;
      return;
    }
    else if (p == "RNA-UV Extended (U)" || p  == "RNA-UV Extended (UCGA)")
    {
      modifications = RNA_UV_EXTENDED_modifications; 
      fragment_adducts = RNA_UV_fragments;
      can_cross_link = (p == "RNA-UV Extended (U)") ? RNA_U : RNA_UCGA ;
      return;
    }
    else if (p == "RNA-UV (4SU)")
    {
      nucleotides.push_back("S=C9H13N2O8PS"); // include 4-Thio-UMP
      mapping.push_back("S->S");
      modifications = RNA_UV_4SU_modifications;
      fragment_adducts = RNA_UV_4SU_fragments;
      can_cross_link = "S";
      return;
    }
    else if (p == "RNA-UV Extended (4SU)")
    {
      nucleotides.push_back("S=C9H13N2O8PS"); // include 4-Thio-UMP
      mapping.push_back("S->S");
      modifications = RNA_UV_4SU_EXTENDED_modifications;
      fragment_adducts = RNA_UV_4SU_fragments;
      can_cross_link = "S";
      return;
    }
    else if (p == "RNA-UV (6SG)")
    {
      nucleotides.push_back("X=C10H14N5O7PS"); // include 6-Thio-GMP
      mapping.push_back("X->X");
      modifications = RNA_UV_6SG_modifications;
      fragment_adducts = RNA_UV_6SG_fragments;
      can_cross_link = "X";
      return;
    }
    else if (p == "RNA-UV Extended (6SG)")
    {
      nucleotides.push_back("X=C10H14N5O7PS"); // include 6-Thio-GMP
      mapping.push_back("X->X");
      modifications = RNA_UV_6SG_EXTENDED_modifications;
      fragment_adducts = RNA_UV_6SG_fragments;
      can_cross_link = "X";
      return;
    }    
    else if (p == "DNA-UV")
    {
      modifications = DNA_UV_modifications;
      fragment_adducts = DNA_UV_fragments;
      can_cross_link = DNA_TCGAd;
      return;
    }
    else if (p == "DNA-UV Extended")
    {
      modifications = DNA_UV_EXTENDED_modifications;
      fragment_adducts = DNA_UV_fragments;
      can_cross_link = DNA_TCGAd;
      return;
    }
    else if (p == "DNA-UV (BrU)")
    {
      nucleotides.push_back("B=C9H11N2O8P"); // include BrU
      mapping.push_back("B->B");
      modifications = DNA_UV_BrU_modifications;
      fragment_adducts = DNA_UV_BrU_fragments;
      can_cross_link = "B";
      return;
    }    
    else if (p == "RNA-FA")
    {
      modifications = RNA_FA_modifications;     
      fragment_adducts = RNA_FA_fragments;
      can_cross_link = RNA_CGA;
      return;
    }
    else if (p == "RNA-FA Extended")
    {
      modifications = RNA_FA_EXTENDED_modifications;     
      fragment_adducts = RNA_FA_fragments;
      can_cross_link = RNA_CGA;
      return;
    }
    else if (p == "DNA-FA")
    {
      modifications = DNA_FA_modifications;     
      fragment_adducts = DNA_FA_fragments;
      can_cross_link = DNA_CGAd;
      return;
    }
    else if (p == "DNA-FA Extended")
    {
      modifications = DNA_FA_EXTENDED_modifications;     
      fragment_adducts = DNA_FA_fragments;
      can_cross_link = DNA_CGAd;
      return;
    }
    else if (p == "RNA-DEB")
    {
      // add special methionine loss
      auto r_ptr = const_cast<Residue*>(ResidueDB::getInstance()->getResidue('M'));
      r_ptr->addLossFormula(EmpiricalFormula("CH4S1"));

      modifications = RNA_DEB_modifications;     
      fragment_adducts = RNA_DEB_fragments;
      can_cross_link = RNA_UCGA;
      return;
    }
    else if (p == "RNA-DEB Extended")
    {
      // add special methionine loss
      auto r_ptr = const_cast<Residue*>(ResidueDB::getInstance()->getResidue('M'));
      r_ptr->addLossFormula(EmpiricalFormula("CH4S1"));

      modifications = RNA_DEB_EXTENDED_modifications;
      fragment_adducts = RNA_DEB_fragments;
      can_cross_link = RNA_UCGA;
      return;
    }
    else if (p == "RNA-NM")
    {
      // add special methionine loss
      auto r_ptr = const_cast<Residue*>(ResidueDB::getInstance()->getResidue('M'));
      r_ptr->addLossFormula(EmpiricalFormula("CH4S1"));

      modifications = RNA_NM_modifications;
      fragment_adducts = RNA_NM_fragments; 
      can_cross_link = RNA_UCGA;
      return;
    }
    else if (p == "RNA-NM Extended")
    {
      // add special methionine loss
      auto r_ptr = const_cast<Residue*>(ResidueDB::getInstance()->getResidue('M'));
      r_ptr->addLossFormula(EmpiricalFormula("CH4S1"));

      modifications = RNA_NM_EXTENDED_modifications;
      fragment_adducts = RNA_NM_fragments; 
      can_cross_link = RNA_UCGA;
      return;
    }    
    else if (p == "DNA-DEB")
    {
      modifications = DNA_DEB_modifications;
      fragment_adducts = DNA_DEB_fragments;
      can_cross_link = DNA_TCGAd;
      return;
    }
    else if (p == "DNA-DEB Extended")
    {
      modifications = DNA_DEB_EXTENDED_modifications;
      fragment_adducts = DNA_DEB_fragments;
      can_cross_link = DNA_TCGAd;
      return;
    }    
    else if (p == "DNA-NM")
    {
      modifications = DNA_NM_modifications;
      fragment_adducts = DNA_NM_fragments;
      can_cross_link = DNA_TCGAd;
      return;
    }
    else if (p == "DNA-NM Extended")
    {
      modifications = DNA_NM_EXTENDED_modifications;
      fragment_adducts = DNA_NM_fragments;
      can_cross_link = DNA_TCGAd;
      return;
    }
  }
  }

}
