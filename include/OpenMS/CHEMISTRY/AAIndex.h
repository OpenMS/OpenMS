// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Alexandra Scherbart $
// --------------------------------------------------------------------------

#ifndef OPENMS_CHEMISTRY_RESIDUE_INDECES_H
#define OPENMS_CHEMISTRY_RESIDUE_INDECES_H

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CHEMISTRY/AASequence.h>

namespace OpenMS
{
	/**
		@brief Represenation of selected %AAIndex properties
		
		Upper-case one-letter-code can be used to access the properties of a single amino acid.
		
		@ingroup Chemistry
	*/
	class AAIndex
	{
		public:
			
			/// Returns if the residue is aliphatic (1.0 or 0.0)
			static DoubleReal aliphatic(char aa) 
			{
				if( aa == 'A' || aa == 'G' || aa == 'F' || aa == 'I' || aa == 'M' || aa == 'L' || aa == 'P' || aa == 'V' ) 
				{
					return 1.0;
				} 
				else 
				{
					return 0.0;
				}
			}
		
			/// Returns if the residue is acidic (1.0 or 0.0)
			static DoubleReal acidic(char aa) 
			{
				if( aa == 'D' || aa == 'E' ) 
				{
					return 1.0;
				}
				else 
				{
					return 0.0;
				}
			}
			
			/// Returns if the residue is basic (1.0 or 0.0)
			static DoubleReal basic(char aa) 
			{
				if( aa == 'K' || aa == 'R' || aa == 'H' || aa == 'W' ) 
				{
					return 1.0;
				}
				else 
				{
					return 0.0;
				}
			}
	
			/// Returns if the residue is polar (1.0 or 0.0)
			static DoubleReal polar(char aa) 
			{
				if( aa == 'S' || aa == 'T' || aa == 'Y' || aa == 'H' || aa == 'C' || aa == 'N' || aa == 'Q' || aa == 'W' ) 
				{
					return 1.0;
				} 
				else
				{
					return 0.0;
				}
			}
			
			/**
				@brief TODO
				
				@exception InvalidValue is thrown if an undefined one-letter-code is used
			*/
			static DoubleReal GBsidechain(char aa)
			{
				switch(aa) 
				{
					case 'A': 
						return 0.0;
					  break;
					case 'C':
						return 0.0;
					  break;
					case 'D':
						return 784.0;
					  break;
					case 'E':
						return 790.0;
					  break;
					case 'F':
						return 0.0;
					  break;
					case 'G':
						return 0.0;
					  break;
					case 'H':
						return 927.84;
					  break;
					case 'I':
						return 0.0;
					  break;
					case 'K':
						return 926.74;
					  break;
					case 'L':
						return 0.0;
					  break;
					case 'M':
						return 830.0;
					  break;
					case 'N':
						return 864.94;
					  break;
					case 'P':
						return 0.0;
					  break;
					case 'Q':
						return 865.25;
					  break;
					case 'R':
						return 1000.0;
					  break;
					case 'S':
						return 775.0;
					  break;
					case 'T':
						return 780.0;
					  break;
					case 'V':
						return 0.0;
					  break;
					case 'W':
						return 909.53;
					  break;
					case 'Y':
						return 790.0;
					  break;
					default:			
						throw Exception::InvalidValue(__FILE__,__LINE__,__PRETTY_FUNCTION__, "Unkown amino acid one-letter-code",String(aa));
				}
			}
		
		
			/**
				@brief TODO
				
				@exception InvalidValue is thrown if an undefined one-letter-code is used
			*/
			static DoubleReal GBleft(char aa) 
			{
				switch(aa) 
				{
					case 'A':
						return 881.82;
						break;
					case 'C':
						return 881.15;
						break;
					case 'D':
						return 880.02;
						break;
					case 'E':
						return 880.10;
						break;
					case 'F':
						return 881.08;
						break;
					case 'G':
						return 881.17;
						break;
					case 'H':
						return 881.27;
						break;
					case 'I':
						return 880.99;
						break;
					case 'K':
						return 880.06;
						break;
					case 'L':
						return 881.88;
						break;
					case 'M':
						return 881.38;
						break;
					case 'N':
						return 881.18;
						break;
					case 'P':
						return 881.25;
						break;
					case 'Q':
						return 881.50;
						break;
					case 'R':
						return 882.98;
						break;
					case 'S':
						return 881.08;
						break;
					case 'T':
						return 881.14;
						break;
					case 'V':
						return 881.17;
						break;
					case 'W':
						return 881.31;
						break;
					case 'Y':
						return 881.20;
						break;
					case '>': //NH2
						return 916.84;
						break;
					default:		
						throw Exception::InvalidValue(__FILE__,__LINE__,__PRETTY_FUNCTION__, "Unkown amino acid one-letter-code",String(aa));
	
				}
			}
		
			/**
				@brief TODO
				
				@exception InvalidValue is thrown if an undefined one-letter-code is used
			*/
			static DoubleReal GBdeltaright(char aa) 
			{
				switch(aa) 
				{
					case 'A':
						return 0.0;
						break;
					case 'C':
						return -0.69;
						break;
					case 'D':
						return -0.63;
						break;
					case 'E':
						return -0.39;
						break;
					case 'F':
						return 0.03;
						break;
					case 'G':
						return 0.92;
						break;
					case 'H':
						return -0.19;
						break;
					case 'I':
						return -1.17;
						break;
					case 'K':
						return -0.71;
						break;
					case 'L':
						return -0.09;
						break;
					case 'M':
						return 0.30;
						break;
					case 'N':
						return 1.56;
						break;
					case 'P':
						return 11.75;
						break;
					case 'Q':
						return 4.10;
						break;
					case 'R':
						return 6.28;
						break;
					case 'S':
						return 0.98;
						break;
					case 'T':
						return 1.21;
						break;
					case 'V':
						return -0.90;
						break;
					case 'W':
						return 0.10;
						break;
					case 'Y':
						return -0.38;
						break;
					case '<': //COOH
						return -95.82;
						break;
					default:		
						throw Exception::InvalidValue(__FILE__,__LINE__,__PRETTY_FUNCTION__, "Unkown amino acid one-letter-code",String(aa));
				}
			}
		
			//
			//H KHAG800101
			//D The Kerr-constant increments (Khanarian-Moore, 1980)
			//R LIT:0611050b
			//A Khanarian, G. and Moore, W.J.
			//T The Kerr effect of amino acids in water
			//J Aust. J. Chem. 33, 1727-1741 (1980) (Cys Lys Tyr !)
			//C 
			//I    A/L     R/K     N/M     D/F     C/P     Q/S     E/T     G/W     H/Y     I/V
			//49.1    133.    -3.6      0.      0.     20.      0.    64.6    75.7    18.9
			//15.6      0.     6.8    54.7    43.8    44.4    31.0    70.5      0.    29.5
	
			/**
				@brief TODO
				
				@exception InvalidValue is thrown if an undefined one-letter-code is used
			*/
			static DoubleReal getKHAG800101(char aa) 
			{
				switch(aa) 
				{
					case 'A':
						return 49.1;
						break;
					case 'R':
						return 133.;
						break;
					case 'N':
						return  -3.6;
						break;
					case 'D':
						return 0.;
						break;
					case 'C':
						return 0.;
						break;
					case 'Q':
						return 20.;
						break;
					case 'E':
						return 0.;
						break;
					case 'G':
						return 64.6;
						break;
					case 'H':
						return 75.7;
						break;
					case 'I':
						return  18.9;
						break;
					case 'L':
						return 15.6;
						break;
					case 'K':
						return 0.;
						break;
					case 'M':
						return 6.8;
						break;
					case 'F':
						return 54.7;
						break;
					case 'P':
						return 43.8;
						break;
					case 'S':
						return 44.4;
						break;
					case 'T':
						return 31.0;
						break;
					case 'W':
						return 70.5;
						break;
					case 'Y':
						return 0.;
						break;
					case 'V':
						return 29.5;
						break;
					default:		
						throw Exception::InvalidValue(__FILE__,__LINE__,__PRETTY_FUNCTION__, "Unkown amino acid one-letter-code",String(aa));
				}
			}
		
			//H VASM830103
			//D Relative population of conformational state E (Vasquez et al., 1983)
			//R LIT:0908110
			//A Vasquez, M., Nemethy, G. and Scheraga, H.A.
			//T Computed conformational states of the 20 naturally occurring amino acid 
			//residues and of the prototype residue alpha-aminobutyric acid
			//J Macromolecules 16, 1043-1049 (1983) (Pro !)
			//C 
			//I    A/L     R/K     N/M     D/F     C/P     Q/S     E/T     G/W     H/Y     I/V
			//0.159   0.194   0.385   0.283   0.187   0.236   0.206   0.049   0.233   0.581
			//0.083   0.159   0.198   0.682   0.366   0.150   0.074   0.463   0.737   0.301	
	
			/**
				@brief TODO
				
				@exception InvalidValue is thrown if an undefined one-letter-code is used
			*/
			static DoubleReal getVASM830103(char aa) 
			{
				switch(aa) 
				{
					case 'A':
						return 0.159;
						break;
					case 'R':
						return 0.194;
						break;
					case 'N':
						return 0.385;
						break;
					case 'D':
						return 0.283;
						break;
					case 'C':
						return 0.187;
						break;
					case 'Q':
						return 0.236;
						break;
					case 'E':
						return 0.206;
						break;
					case 'G':
						return 0.049;
						break;
					case 'H':
						return 0.233;
						break;
					case 'I':
						return 0.581;
						break;
					case 'L':
						return 0.083;
						break;
					case 'K':
						return 0.159;
						break;
					case 'M':
						return 0.198;
						break;
					case 'F':
						return 0.682;
						break;
					case 'P':
						return 0.366;
						break;
					case 'S':
						return 0.150;
						break;
					case 'T':
						return 0.074;
						break;
					case 'W':
						return 0.463;
						break;
					case 'Y':
						return 0.737;
						break;
					case 'V':
						return 0.301;
						break;	
					default:		
						throw Exception::InvalidValue(__FILE__,__LINE__,__PRETTY_FUNCTION__, "Unkown amino acid one-letter-code",String(aa));
				}
			}
		
			
			//
			//H NADH010106
			//D Hydropathy scale based on self-information values in the two-state model (36% 
			//accessibility) (Naderi-Manesh et al., 2001)
			//R PMID:11170200
			//A Naderi-Manesh, H., Sadeghi, M., Arab, S. and Moosavi Movahedi, A.A.
			//T Prediction of protein surface accessibility with information theory
			//J Proteins. 42, 452-459 (2001)
			//C NADH010105    0.958  NADH010104    0.914  NADH010103    0.881
			//ZHOH040103    0.819  NADH010107    0.811  BAEK050101    0.809
			//NADH010102    0.808  PONP800103    0.803  VINM940103   -0.813
			//KRIW710101   -0.846  KRIW790101   -0.861
			//I    A/L     R/K     N/M     D/F     C/P     Q/S     E/T     G/W     H/Y     I/V
			//5     -57     -77      45     224     -67      -8     -47     -50      83
			//82     -38      83     117    -103     -41      79     130      27     117
	
			/**
				@brief TODO
				
				@exception InvalidValue is thrown if an undefined one-letter-code is used
			*/
			static DoubleReal getNADH010106(char aa) 
			{
				switch(aa) 
				{
					case 'A':
						return 5;
						break;
					case 'R':
						return -57;
						break;
					case 'N':
						return -77;
						break;
					case 'D':
						return 45;
						break;
					case 'C':
						return 224;
						break;
					case 'Q':
						return -67;
						break;
					case 'E':
						return -8;
						break;
					case 'G':
						return -47;
						break;
					case 'H':
						return -50;
						break;
					case 'I':
						return 83;
						break;
					case 'L':
						return 82;
						break;
					case 'K':
						return -38;
						break;
					case 'M':
						return 83;
						break;
					case 'F':
						return 117;
						break;
					case 'P':
						return -103;
						break;
					case 'S':
						return -41;
						break;
					case 'T':
						return 79;
						break;
					case 'W':
						return 130;
						break;
					case 'Y':
						return 27;
						break;
					case 'V':
						return 117;
						break;	
					default:		
						throw Exception::InvalidValue(__FILE__,__LINE__,__PRETTY_FUNCTION__, "Unkown amino acid one-letter-code",String(aa));
				}
			}
		
		
			//
			//H NADH010107
			//D Hydropathy scale based on self-information values in the two-state model (50% 
			//accessibility) (Naderi-Manesh et al., 2001)
			//R PMID:11170200
			//A Naderi-Manesh, H., Sadeghi, M., Arab, S. and Moosavi Movahedi, A.A.
			//T Prediction of protein surface accessibility with information theory
			//J Proteins. 42, 452-459 (2001)
			//C NADH010106    0.811
			//I    A/L     R/K     N/M     D/F     C/P     Q/S     E/T     G/W     H/Y     I/V
			//-2     -41     -97     248     329     -37     117     -66     -70      28
			//36     115      62     120    -132     -52     174     179      -7     114  
	
			/**
				@brief TODO
				
				@exception InvalidValue is thrown if an undefined one-letter-code is used
			*/	
			static DoubleReal getNADH010107(char aa) 
			{
				switch(aa) 
				{
					case 'A':
						return -2;
						break;
					case 'R':
						return -41;
						break;
					case 'N':
						return -97;
						break;
					case 'D':
						return 248;
						break;
					case 'C':
						return 329;
						break;
					case 'Q':
						return -37;
						break;
					case 'E':
						return 117;
						break;
					case 'G':
						return -66;
						break;
					case 'H':
						return -70;
						break;
					case 'I':
						return 28;
						break;
					case 'L':
						return 36;
						break;
					case 'K':
						return 115;
						break;
					case 'M':
						return 62;
						break;
					case 'F':
						return 120;
						break;
					case 'P':
						return -132;
						break;
					case 'S':
						return -52;
						break;
					case 'T':
						return 174;
						break;
					case 'W':
						return 179;
						break;
					case 'Y':
						return -7;
						break;
					case 'V':
						return 114;
						break;	
					default:		
						throw Exception::InvalidValue(__FILE__,__LINE__,__PRETTY_FUNCTION__, "Unkown amino acid one-letter-code",String(aa));
				}
			}
		
		
			//
			//H WILM950102
			//D Hydrophobicity coefficient in RP-HPLC, C8 with 0.1%TFA/MeCN/H2O (Wilce et al. 
			//1995)
			//R 
			//A Wilce, M.C., Aguilar, M.I. and Hearn, M.T.
			//T Physicochemical basis of amino acid hydrophobicity scales: evaluation of four 
			//new scales of amino acid hydrophobicity coefficients derived from RP-HPLC of 
			//peptides
			//J Anal Chem. 67, 1210-1219 (1995)
			//C WILM950101    0.838  MEEJ810102    0.809
			//I    A/L     R/K     N/M     D/F     C/P     Q/S     E/T     G/W     H/Y     I/V
			//2.62    1.26   -1.27   -2.84    0.73   -1.69   -0.45   -1.15   -0.74    4.38
			//6.57   -2.78   -3.12    9.14   -0.12   -1.39    1.81    5.91    1.39    2.30		
	
			/**
				@brief TODO
				
				@exception InvalidValue is thrown if an undefined one-letter-code is used
			*/
			static DoubleReal getWILM950102(char aa) 
			{
				switch(aa) 
				{
					case 'A':
						return 2.62;
						break;
					case 'R':
						return 1.26;
						break;
					case 'N':
						return -1.27;
						break;
					case 'D':
						return -2.84;
						break;
					case 'C':
						return 0.73;
						break;
					case 'Q':
						return -1.69;
						break;
					case 'E':
						return -0.45;
						break;
					case 'G':
						return -1.15;
						break;
					case 'H':
						return -0.74;
						break;
					case 'I':
						return 4.38;
						break;
					case 'L':
						return 6.57;
						break;
					case 'K':
						return -2.78;
						break;
					case 'M':
						return -3.12;
						break;
					case 'F':
						return 9.14;
						break;
					case 'P':
						return -0.12;
						break;
					case 'S':
						return -1.39;
						break;
					case 'T':
						return 1.81;
						break;
					case 'W':
						return 5.91;
						break;
					case 'Y':
						return 1.39;
						break;
					case 'V':
						return 2.30;
						break;	
					default:		
						throw Exception::InvalidValue(__FILE__,__LINE__,__PRETTY_FUNCTION__, "Unkown amino acid one-letter-code",String(aa));
				}
			}	
		
		
			//
			//H ROBB760107
			//D Information measure for extended without H-bond (Robson-Suzuki, 1976)
			//R PMID:1003471
			//A Robson, B. and Suzuki, E.
			//T Conformational properties of amino acid residues in globular proteins
			//J J. Mol. Biol. 107, 327-356 (1976)
			//C 
			//I    A/L     R/K     N/M     D/F     C/P     Q/S     E/T     G/W     H/Y     I/V
			//0.0     1.1    -2.0    -2.6     5.4     2.4     3.1    -3.4     0.8    -0.1
			//-3.7    -3.1    -2.1     0.7     7.4     1.3     0.0    -3.4     4.8     2.7
	
			/**
				@brief TODO
				
				@exception InvalidValue is thrown if an undefined one-letter-code is used
			*/
			static DoubleReal getROBB760107(char aa) 
			{
				switch(aa) 
				{
					case 'A':
						return 0.0;
						break;
					case 'R':
						return 1.1;
						break;
					case 'N':
						return -2.0;
						break;
					case 'D':
						return -2.6;
						break;
					case 'C':
						return 5.4;
						break;
					case 'Q':
						return 2.4;
						break;
					case 'E':
						return 3.1;
						break;
					case 'G':
						return -3.4;
						break;
					case 'H':
						return 0.8;
						break;
					case 'I':
						return -0.1;
						break;
					case 'L':
						return -3.7;
						break;
					case 'K':
						return -3.1;
						break;
					case 'M':
						return -2.1;
						break;
					case 'F':
						return 0.7;
						break;
					case 'P':
						return 7.4;
						break;
					case 'S':
						return 1.3;
						break;
					case 'T':
						return 0.0;
						break;
					case 'W':
						return -3.4;
						break;
					case 'Y':
						return 4.8;
						break;
					case 'V':
						return 2.7;
						break;	
					default:		
						throw Exception::InvalidValue(__FILE__,__LINE__,__PRETTY_FUNCTION__, "Unkown amino acid one-letter-code",String(aa));
				}
			}
		
		
			//
			//H OOBM850104
			//D Optimized average non-bonded energy per atom (Oobatake et al., 1985)
			//R LIT:1207075b
			//A Oobatake, M., Kubota, Y. and Ooi, T.
			//T Optimization of amino acid parameters for correspondence of sequence to 
			//tertiary structures of proteuins
			//J Bull. Inst. Chem. Res., Kyoto Univ. 63, 82-94 (1985)
			//C 
			//I    A/L     R/K     N/M     D/F     C/P     Q/S     E/T     G/W     H/Y     I/V
			//-2.49    2.55    2.27    8.86   -3.13    1.79    4.04   -0.56    4.22  -10.87
			//-7.16   -9.97   -4.96   -6.64    5.19   -1.60   -4.75  -17.84    9.25   -3.97  
	
			/**
				@brief TODO
				
				@exception InvalidValue is thrown if an undefined one-letter-code is used
			*/
			static DoubleReal getOOBM850104(char aa) 
			{
				switch(aa) 
				{
					case 'A':
						return -2.49;
						break;
					case 'R':
						return 2.55;
						break;
					case 'N':
						return 2.27;
						break;
					case 'D':
						return 8.86;
						break;
					case 'C':
						return -3.13;
						break;
					case 'Q':
						return 1.79;
						break;
					case 'E':
						return 4.04;
						break;
					case 'G':
						return -0.56;
						break;
					case 'H':
						return 4.22;
						break;
					case 'I':
						return -10.87;
						break;
					case 'L':
						return -7.16;
						break;
					case 'K':
						return -9.97;
						break;
					case 'M':
						return -4.96;
						break;
					case 'F':
						return -6.64;
						break;
					case 'P':
						return 5.19;
						break;
					case 'S':
						return -1.60;
						break;
					case 'T':
						return -4.75;
						break;
					case 'W':
						return -17.84;
						break;
					case 'Y':
						return 9.25;
						break;
					case 'V':
						return -3.97;
						break;	
					default:		
						throw Exception::InvalidValue(__FILE__,__LINE__,__PRETTY_FUNCTION__, "Unkown amino acid one-letter-code",String(aa));
				}
			}
		
		
			//
			//H FAUJ880111
			//D Positive charge (Fauchere et al., 1988)
			//R LIT:1414114 PMID:3209351
			//A Fauchere, J.L., Charton, M., Kier, L.B., Verloop, A. and Pliska, V.
			//T Amino acid side chain parameters for correlation studies in biology and 
			//pharmacology
			//J Int. J. Peptide Protein Res. 32, 269-278 (1988)
			//C ZIMJ680104    0.813
			//I    A/L     R/K     N/M     D/F     C/P     Q/S     E/T     G/W     H/Y     I/V
			//0.      1.      0.      0.      0.      0.      0.      0.      1.      0.
			//0.      1.      0.      0.      0.      0.      0.      0.      0.      0.
	
			/**
				@brief TODO
				
				@exception InvalidValue is thrown if an undefined one-letter-code is used
			*/
			static DoubleReal getFAUJ880111(char aa) 
			{
				switch(aa) 
				{
					case 'A':
						return 0.;
						break;
					case 'R':
						return 1.;
						break;
					case 'N':
						return 0.;
						break;
					case 'D':
						return 0.;
						break;
					case 'C':
						return 0.;
						break;
					case 'Q':
						return 0.;
						break;
					case 'E':
						return 0.;
						break;
					case 'G':
						return 0.;
						break;
					case 'H':
						return 1.;
						break;
					case 'I':
						return 0.;
						break;
					case 'L':
						return 0.;
						break;
					case 'K':
						return 1.;
						break;
					case 'M':
						return 0.;
						break;
					case 'F':
						return 0.;
						break;
					case 'P':
						return 0.;
						break;
					case 'S':
						return 0.;
						break;
					case 'T':
						return 0.;
						break;
					case 'W':
						return 0.;
						break;
					case 'Y':
						return 0.;
						break;
					case 'V':
						return 0.;
						break;	
					default:		
						throw Exception::InvalidValue(__FILE__,__LINE__,__PRETTY_FUNCTION__, "Unkown amino acid one-letter-code",String(aa));
				}
			}
		
		
			//
			//H FINA770101
			//D Helix-coil equilibrium constant (Finkelstein-Ptitsyn, 1977)
			//R LIT:2004052b PMID:843599
			//A Finkelstein, A.V. and Ptitsyn, O.B.
			//T Theory of protein molecule self-organization. II. A comparison of calculated 
			//thermodynamic parameters of local secondary structures with experiments
			//J Biopolymers 16, 497-524 (1977) (Pro 0.096)
			//C SUEM840101    0.883  AURR980114    0.875  AURR980113    0.849
			//PTIO830101    0.826  KANM800103    0.823  QIAN880107    0.814
			//QIAN880106    0.810  MAXF760101    0.810  AURR980109    0.802
			//I    A/L     R/K     N/M     D/F     C/P     Q/S     E/T     G/W     H/Y     I/V
			//1.08    1.05    0.85    0.85    0.95    0.95    1.15    0.55    1.00    1.05
			//1.25    1.15    1.15    1.10    0.71    0.75    0.75    1.10    1.10    0.95
	
			/**
				@brief TODO
				
				@exception InvalidValue is thrown if an undefined one-letter-code is used
			*/
			static DoubleReal getFINA770101(char aa) 
			{
				switch(aa) 
				{
					case 'A':
						return 1.08;
						break;
					case 'R':
						return 1.05;
						break;
					case 'N':
						return 0.85;
						break;
					case 'D':
						return 0.85;
						break;
					case 'C':
						return 0.95;
						break;
					case 'Q':
						return 0.95;
						break;
					case 'E':
						return 1.15;
						break;
					case 'G':
						return 0.55;
						break;
					case 'H':
						return 1.00;
						break;
					case 'I':
						return 1.05;
						break;
					case 'L':
						return 1.25;
						break;
					case 'K':
						return 1.15;
						break;
					case 'M':
						return 1.15;
						break;
					case 'F':
						return 1.10;
						break;
					case 'P':
						return 0.71;
						break;
					case 'S':
						return 0.75;
						break;
					case 'T':
						return 0.75;
						break;
					case 'W':
						return 1.10;
						break;
					case 'Y':
						return 1.10;
						break;
					case 'V':
						return 0.95;
						break;	
					default:		
						throw Exception::InvalidValue(__FILE__,__LINE__,__PRETTY_FUNCTION__, "Unkown amino acid one-letter-code",String(aa));
				}
			}
		
		
			//
			//H ARGP820102
			//D Signal sequence helical potential (Argos et al., 1982)
			//R LIT:0901079b PMID:7151796
			//A Argos, P., Rao, J.K.M. and Hargrave, P.A.
			//T Structural prediction of membrane-bound proteins
			//J Eur. J. Biochem. 128, 565-575 (1982)
			//C ARGP820103    0.961  KYTJ820101    0.803  JURD980101    0.802
			//I    A/L     R/K     N/M     D/F     C/P     Q/S     E/T     G/W     H/Y     I/V
			//1.18    0.20    0.23    0.05    1.89    0.72    0.11    0.49    0.31    1.45
			//3.23    0.06    2.67    1.96    0.76    0.97    0.84    0.77    0.39    1.08
	
			/**
				@brief TODO
				
				@exception InvalidValue is thrown if an undefined one-letter-code is used
			*/
			static DoubleReal getARGP820102(char aa) 
			{
				switch(aa) 
				{
					case 'A':
						return 1.18;
						break;
					case 'R':
						return 0.20;
						break;
					case 'N':
						return 0.23;
						break;
					case 'D':
						return 0.05;
						break;
					case 'C':
						return 1.89;
						break;
					case 'Q':
						return 0.72;
						break;
					case 'E':
						return 0.11;
						break;
					case 'G':
						return 0.49;
						break;
					case 'H':
						return 0.31;
						break;
					case 'I':
						return 1.45;
						break;
					case 'L':
						return 3.23;
						break;
					case 'K':
						return 0.06;
						break;
					case 'M':
						return 2.67;
						break;
					case 'F':
						return 1.96;
						break;
					case 'P':
						return 0.76;
						break;
					case 'S':
						return 0.97;
						break;
					case 'T':
						return 0.84;
						break;
					case 'W':
						return 0.77;
						break;
					case 'Y':
						return 0.39;
						break;
					case 'V':
						return 1.08;
						break;	
					default:		
						throw Exception::InvalidValue(__FILE__,__LINE__,__PRETTY_FUNCTION__, "Unkown amino acid one-letter-code",String(aa));
				}
			}
	
			/**
				@brief Calculates an array of properties for an amino acid sequence
				
				The array contains the following properties:
				- 0: Number of 'R' residues
				- 1: Signal sequence helical potential
				- 2: Number of 'F' residues
				- 3: Positive charge
				- 4: Helix-coil equilibrium constant
				- 5: Estimated gas-phase basicity at 500 K
				- 6: Number of 'H' residues
				- 7: Kerr-constant increments
				- 8: Number of 'M' residues
				- 9: Average amino acid weight
				- 10: Hydropathy scale (36% accessibility)
				- 11: Hydropathy scale (50% accessibility)
				- 12: Optimized average non-bonded energy per atom
				- 13: Number of 'Q' residues
				- 14: Information measure for extended without H-bond
				- 15: Relative population of conformational state E
				- 16: Hydrophobicity coefficient in RP-HPLC, C8 with 0.1%TFA/MeCN/H2 O,
				- 17: Number of 'Y' residues
				
				@exception InvalidValue is thrown if an undefined one-letter-code is used
			*/
			static std::vector<DoubleReal> getPropertyVector(const AASequence& sequence)
			{
				std::vector<DoubleReal> out(18);
				
				//for each element in sequence = residue
				for(Size pos=0; pos<sequence.size(); pos++)
				{
					char amino = sequence[pos].getOneLetterCode()[0];
					
					// numResidues of R
					out[0] += (amino=='R' ? 1.0 : 0.0);
					//The Kerr-constant increments
					out[7] += getKHAG800101(amino); 
					//Relative population of conformational state E
					out[15] += getVASM830103(amino); 
					//Hydropathy scale (36% accessibility)
					out[10] += getNADH010106(amino);
					//Hydropathy scale (50% accessibility)
					out[11] += getNADH010107(amino);
					//Hydrophobicity coefficient in RP-HPLC, C8 with 0.1%TFA/MeCN/H2 O,
					out[16] += getWILM950102(amino);
					//Information measure for extended without H-bond,
					out[14] += getROBB760107(amino);
					//Optimized average non-bonded energy per atom,
					out[12] += getOOBM850104(amino);
					//Positive charge
					out[3] += getFAUJ880111(amino);
					//Helix-coil equilibrium constant
					out[4] += getFINA770101(amino);
					//Signal sequence helical potential
					out[1] += getARGP820102(amino);
					out[8] += (amino=='M' ? 1.0 : 0.0);// numResidues of M
					out[2] += (amino=='F' ? 1.0 : 0.0);// numResidues of F
					out[6] += (amino=='H' ? 1.0 : 0.0);// numResidues of H
					out[13] += (amino=='Q' ? 1.0 : 0.0);// numResidues of Q
					out[17] += (amino=='Y' ? 1.0 : 0.0);// numResiduesof Y
				}	
				
				out[5] = calculateGB(sequence, 500.0); //Estimated gas-phase basicity at 500 K
				out[9] = sequence.getAverageWeight();
							
				return out;					
			}
			
			/**
				@brief Calculates an estimated gas-phase basicity for an amino acid sequence at a given temperature
				
				@exception InvalidValue is thrown if an undefined one-letter-code is used
			*/
			static DoubleReal calculateGB(const AASequence& seq, DoubleReal T=500.0) 
			{
						
				DoubleReal R = Constants::GAS_CONSTANT/1000.0; // ideal gas constant in kj/(K*mol)
					
				char left = '>';
				char right;
				
				DoubleReal k_app = 0.0; // apparent proton association constant
					
				// energy level E at each protonation site i is -GB(i)
				// fractional proton population of a microstate k is 
				// P_k = exp (- E_k/(RT)) / ( sum_i exp (- E_i/(RT))) 
				// the apparent proton association constant k_app:
				// k_app = sum_i GB(i)/(RT)
				// then the apparent GB is GB_app^ion = R * T * ln(k_app)
				for (Size i = 0; i <= seq.size(); i++) 
				{
					// aa left to current one
					if(i > 0) 
					{
						Residue leftchar = seq[i-1];
						left = leftchar.getOneLetterCode()[0];		
					} 
						
					// aa right to current one
					if(i == seq.size()) 
					{
						right = '<';
					} 
					else 
					{
						Residue rightchar = seq[i];
						right = rightchar.getOneLetterCode()[0];
					}
					DoubleReal contrib = exp((GBleft(left) + GBdeltaright(right))/(R*T));
					if(i > 0 && i < seq.size()) 
					{
						contrib += exp(GBsidechain(right)/(R*T));
					}
					k_app += contrib;
				}
				// calculate apparent GB
				return R * T * log(k_app)/log(2.0);
			}
		
		private:
			
			///Constructor not implemented => private
			AAIndex();
	};
	
}
#endif 
