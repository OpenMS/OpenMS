// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <cmath>

namespace OpenMS
{
  /**
    @brief Representation of selected %AAIndex properties

    The literature that describe the indices can be found with:
    @n Kawashima, S., Ogata, H., and Kanehisa, M. (1999).
    @n <em>AAindex: Amino Acid Index Database</em>,
    @n Nucleic Acids Res, 27(1), 368&ndash;369.

    The provided values are:
      - GB500      Estimated gas-phase basicity at 500 K,
      - VASM830103 Relative population of conformational state E,
      - NADH010106 Hydropathy scale (36% accessibility),
      - FAUJ880111 Positive charge,
      - WILM950102 Hydrophobicity coefficient in RP-HPLC, C8 with 0.1%TFA/MeCN/H2 O,
      - OOBM850104 Optimized average non-bonded energy per atom,
      - KHAG800101 The Kerr-constant increments,
      - NADH010107 Hydropathy scale (50% accessibility),
      - ROBB760107 Information measure for extended without H-bond,
      - FINA770101 Helix-coil equilibrium constant,
      - ARGP820102 Signal sequence helical potential.

    Upper-case one-letter-code can be used to access the properties of a single amino acid.

    @ingroup Chemistry
  */
  class OPENMS_DLLAPI AAIndex
  {

public:
    /// Constructor not implemented
    AAIndex() = delete;

    /// Returns if the residue is aliphatic (1.0 or 0.0)
    static double aliphatic(char aa)
    {
      if (aa == 'A' || aa == 'G' || aa == 'F' || aa == 'I' || aa == 'M' || aa == 'L' || aa == 'P' || aa == 'V')
      {
        return 1.0;
      }
      else
      {
        return 0.0;
      }
    }

    /// Returns if the residue is acidic (1.0 or 0.0)
    static double acidic(char aa)
    {
      if (aa == 'D' || aa == 'E')
      {
        return 1.0;
      }
      else
      {
        return 0.0;
      }
    }

    /// Returns if the residue is basic (1.0 or 0.0)
    static double basic(char aa)
    {
      if (aa == 'K' || aa == 'R' || aa == 'H' || aa == 'W')
      {
        return 1.0;
      }
      else
      {
        return 0.0;
      }
    }

    /// Returns if the residue is polar (1.0 or 0.0)
    static double polar(char aa)
    {
      if (aa == 'S' || aa == 'T' || aa == 'Y' || aa == 'H' || aa == 'C' || aa == 'N' || aa == 'Q' || aa == 'W')
      {
        return 1.0;
      }
      else
      {
        return 0.0;
      }
    }

    //I    A/L     R/K     N/M     D/F     C/P     Q/S     E/T     G/W     H/Y     I/V
    //49.1    133.    -3.6      0.      0.     20.      0.    64.6    75.7    18.9
    //15.6      0.     6.8    54.7    43.8    44.4    31.0    70.5      0.    29.5
    /**
      @brief The Kerr-constant increments (Khanarian-Moore, 1980)

      LIT:0611050b<br>
      Khanarian, G. and Moore, W.J.<br>
      The Kerr effect of amino acids in water<br>
      Aust. J. Chem. 33, 1727-1741 (1980) (Cys Lys Tyr !)

      @exception InvalidValue is thrown if an undefined one-letter-code is used
    */
    static double getKHAG800101(char aa)
    {
      switch (aa)
      {
      case 'A':
        return 49.1;

      case 'R':
        return 133.;

      case 'N':
        return -3.6;

      case 'D':
        return 0.;

      case 'C':
        return 0.;

      case 'Q':
        return 20.;

      case 'E':
        return 0.;

      case 'G':
        return 64.6;

      case 'H':
        return 75.7;

      case 'I':
        return 18.9;

      case 'L':
        return 15.6;

      case 'K':
        return 0.;

      case 'M':
        return 6.8;

      case 'F':
        return 54.7;

      case 'P':
        return 43.8;

      case 'S':
        return 44.4;

      case 'T':
        return 31.0;

      case 'W':
        return 70.5;

      case 'Y':
        return 0.;

      case 'V':
        return 29.5;

      default:
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Unknown amino acid one-letter-code", String(aa));
      }
    }

    //I    A/L     R/K     N/M     D/F     C/P     Q/S     E/T     G/W     H/Y     I/V
    //0.159   0.194   0.385   0.283   0.187   0.236   0.206   0.049   0.233   0.581
    //0.083   0.159   0.198   0.682   0.366   0.150   0.074   0.463   0.737   0.301

    /**
      @brief Relative population of conformational state E (Vasquez et al., 1983)

      LIT:0908110<br>
      Vasquez, M., Nemethy, G. and Scheraga, H.A.<br>
      Computed conformational states of the 20 naturally occurring amino acid
      residues and of the prototype residue alpha-aminobutyric acid<br>
      Macromolecules 16, 1043-1049 (1983) (Pro !)

      @exception InvalidValue is thrown if an undefined one-letter-code is used
    */
    static double getVASM830103(char aa)
    {
      switch (aa)
      {
      case 'A':
        return 0.159;

      case 'R':
        return 0.194;

      case 'N':
        return 0.385;

      case 'D':
        return 0.283;

      case 'C':
        return 0.187;

      case 'Q':
        return 0.236;

      case 'E':
        return 0.206;

      case 'G':
        return 0.049;

      case 'H':
        return 0.233;

      case 'I':
        return 0.581;

      case 'L':
        return 0.083;

      case 'K':
        return 0.159;

      case 'M':
        return 0.198;

      case 'F':
        return 0.682;

      case 'P':
        return 0.366;

      case 'S':
        return 0.150;

      case 'T':
        return 0.074;

      case 'W':
        return 0.463;

      case 'Y':
        return 0.737;

      case 'V':
        return 0.301;

      default:
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Unknown amino acid one-letter-code", String(aa));
      }
    }

    //NADH010105    0.958  NADH010104    0.914  NADH010103    0.881<br>
    //ZHOH040103    0.819  NADH010107    0.811  BAEK050101    0.809<br>
    //NADH010102    0.808  PONP800103    0.803  VINM940103   -0.813<br>
    //KRIW710101   -0.846  KRIW790101   -0.861
    //I    A/L     R/K     N/M     D/F     C/P     Q/S     E/T     G/W     H/Y     I/V
    //5     -57     -77      45     224     -67      -8     -47     -50      83
    //82     -38      83     117    -103     -41      79     130      27     117

    /**
      @brief Hydropathy scale based on self-information values in the two-state model (36% accessibility) (Naderi-Manesh et al., 2001)

      PMID:11170200<br>
      Naderi-Manesh, H., Sadeghi, M., Arab, S. and Moosavi Movahedi, A.A.<br>
      Prediction of protein surface accessibility with information theory<br>
      Proteins. 42, 452-459 (2001)

      @exception InvalidValue is thrown if an undefined one-letter-code is used
    */
    static double getNADH010106(char aa)
    {
      switch (aa)
      {
      case 'A':
        return 5;

      case 'R':
        return -57;

      case 'N':
        return -77;

      case 'D':
        return 45;

      case 'C':
        return 224;

      case 'Q':
        return -67;

      case 'E':
        return -8;

      case 'G':
        return -47;

      case 'H':
        return -50;

      case 'I':
        return 83;

      case 'L':
        return 82;

      case 'K':
        return -38;

      case 'M':
        return 83;

      case 'F':
        return 117;

      case 'P':
        return -103;

      case 'S':
        return -41;

      case 'T':
        return 79;

      case 'W':
        return 130;

      case 'Y':
        return 27;

      case 'V':
        return 117;

      default:
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Unknown amino acid one-letter-code", String(aa));
      }
    }

    //NADH010106    0.811
    //I    A/L     R/K     N/M     D/F     C/P     Q/S     E/T     G/W     H/Y     I/V
    //-2     -41     -97     248     329     -37     117     -66     -70      28
    //36     115      62     120    -132     -52     174     179      -7     114

    /**
      @brief Hydropathy scale based on self-information values in the two-state model (50% accessibility) (Naderi-Manesh et al., 2001)

      PMID:11170200<br>
      Naderi-Manesh, H., Sadeghi, M., Arab, S. and Moosavi Movahedi, A.A.<br>
      Prediction of protein surface accessibility with information theory<br>
      Proteins. 42, 452-459 (2001)

      @exception InvalidValue is thrown if an undefined one-letter-code is used
    */
    static double getNADH010107(char aa)
    {
      switch (aa)
      {
      case 'A':
        return -2;

      case 'R':
        return -41;

      case 'N':
        return -97;

      case 'D':
        return 248;

      case 'C':
        return 329;

      case 'Q':
        return -37;

      case 'E':
        return 117;

      case 'G':
        return -66;

      case 'H':
        return -70;

      case 'I':
        return 28;

      case 'L':
        return 36;

      case 'K':
        return 115;

      case 'M':
        return 62;

      case 'F':
        return 120;

      case 'P':
        return -132;

      case 'S':
        return -52;

      case 'T':
        return 174;

      case 'W':
        return 179;

      case 'Y':
        return -7;

      case 'V':
        return 114;

      default:
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Unknown amino acid one-letter-code", String(aa));
      }
    }

    //WILM950101    0.838  MEEJ810102    0.809
    //I    A/L     R/K     N/M     D/F     C/P     Q/S     E/T     G/W     H/Y     I/V
    //2.62    1.26   -1.27   -2.84    0.73   -1.69   -0.45   -1.15   -0.74    4.38
    //6.57   -2.78   -3.12    9.14   -0.12   -1.39    1.81    5.91    1.39    2.30

    /**
      @brief Hydrophobicity coefficient in RP-HPLC, C8 with 0.1%TFA/MeCN/H2O (Wilce et al. 1995)

      Wilce, M.C., Aguilar, M.I. and Hearn, M.T.<br>
      Physicochemical basis of amino acid hydrophobicity scales: evaluation of four
      new scales of amino acid hydrophobicity coefficients derived from RP-HPLC of
      peptides<br>
      Anal Chem. 67, 1210-1219 (1995)

      @exception InvalidValue is thrown if an undefined one-letter-code is used
    */
    static double getWILM950102(char aa)
    {
      switch (aa)
      {
      case 'A':
        return 2.62;

      case 'R':
        return 1.26;

      case 'N':
        return -1.27;

      case 'D':
        return -2.84;

      case 'C':
        return 0.73;

      case 'Q':
        return -1.69;

      case 'E':
        return -0.45;

      case 'G':
        return -1.15;

      case 'H':
        return -0.74;

      case 'I':
        return 4.38;

      case 'L':
        return 6.57;

      case 'K':
        return -2.78;

      case 'M':
        return -3.12;

      case 'F':
        return 9.14;

      case 'P':
        return -0.12;

      case 'S':
        return -1.39;

      case 'T':
        return 1.81;

      case 'W':
        return 5.91;

      case 'Y':
        return 1.39;

      case 'V':
        return 2.30;

      default:
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Unknown amino acid one-letter-code", String(aa));
      }
    }

    //I    A/L     R/K     N/M     D/F     C/P     Q/S     E/T     G/W     H/Y     I/V
    //0.0     1.1    -2.0    -2.6     5.4     2.4     3.1    -3.4     0.8    -0.1
    //-3.7    -3.1    -2.1     0.7     7.4     1.3     0.0    -3.4     4.8     2.7

    /**
      @brief Information measure for extended without H-bond (Robson-Suzuki, 1976)

      PMID:1003471<br>
      Robson, B. and Suzuki, E.<br>
      Conformational properties of amino acid residues in globular proteins<br>
      J. Mol. Biol. 107, 327-356 (1976)

      @exception InvalidValue is thrown if an undefined one-letter-code is used
    */
    static double getROBB760107(char aa)
    {
      switch (aa)
      {
      case 'A':
        return 0.0;

      case 'R':
        return 1.1;

      case 'N':
        return -2.0;

      case 'D':
        return -2.6;

      case 'C':
        return 5.4;

      case 'Q':
        return 2.4;

      case 'E':
        return 3.1;

      case 'G':
        return -3.4;

      case 'H':
        return 0.8;

      case 'I':
        return -0.1;

      case 'L':
        return -3.7;

      case 'K':
        return -3.1;

      case 'M':
        return -2.1;

      case 'F':
        return 0.7;

      case 'P':
        return 7.4;

      case 'S':
        return 1.3;

      case 'T':
        return 0.0;

      case 'W':
        return -3.4;

      case 'Y':
        return 4.8;

      case 'V':
        return 2.7;

      default:
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Unknown amino acid one-letter-code", String(aa));
      }
    }

    //I    A/L     R/K     N/M     D/F     C/P     Q/S     E/T     G/W     H/Y     I/V
    //-2.49    2.55    2.27    8.86   -3.13    1.79    4.04   -0.56    4.22  -10.87
    //-7.16   -9.97   -4.96   -6.64    5.19   -1.60   -4.75  -17.84    9.25   -3.97

    /**
      @brief Optimized average non-bonded energy per atom (Oobatake et al., 1985)

      LIT:1207075b<br>
      Oobatake, M., Kubota, Y. and Ooi, T.<br>
      Optimization of amino acid parameters for correspondence of sequence to
      tertiary structures of proteins<br>
      Bull. Inst. Chem. Res., Kyoto Univ. 63, 82-94 (1985)

      @exception InvalidValue is thrown if an undefined one-letter-code is used
    */
    static double getOOBM850104(char aa)
    {
      switch (aa)
      {
      case 'A':
        return -2.49;

      case 'R':
        return 2.55;

      case 'N':
        return 2.27;

      case 'D':
        return 8.86;

      case 'C':
        return -3.13;

      case 'Q':
        return 1.79;

      case 'E':
        return 4.04;

      case 'G':
        return -0.56;

      case 'H':
        return 4.22;

      case 'I':
        return -10.87;

      case 'L':
        return -7.16;

      case 'K':
        return -9.97;

      case 'M':
        return -4.96;

      case 'F':
        return -6.64;

      case 'P':
        return 5.19;

      case 'S':
        return -1.60;

      case 'T':
        return -4.75;

      case 'W':
        return -17.84;

      case 'Y':
        return 9.25;

      case 'V':
        return -3.97;

      default:
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Unknown amino acid one-letter-code", String(aa));
      }
    }

    //ZIMJ680104    0.813
    //I    A/L     R/K     N/M     D/F     C/P     Q/S     E/T     G/W     H/Y     I/V
    //0.      1.      0.      0.      0.      0.      0.      0.      1.      0.
    //0.      1.      0.      0.      0.      0.      0.      0.      0.      0.

    /**
      @brief Positive charge (Fauchere et al., 1988)

      LIT:1414114 PMID:3209351<br>
      Fauchere, J.L., Charton, M., Kier, L.B., Verloop, A. and Pliska, V.<br>
      Amino acid side chain parameters for correlation studies in biology and
      pharmacology<br>
      Int. J. Peptide Protein Res. 32, 269-278 (1988)

      @exception InvalidValue is thrown if an undefined one-letter-code is used
    */
    static double getFAUJ880111(char aa)
    {
      switch (aa)
      {
      case 'A':
        return 0.;

      case 'R':
        return 1.;

      case 'N':
        return 0.;

      case 'D':
        return 0.;

      case 'C':
        return 0.;

      case 'Q':
        return 0.;

      case 'E':
        return 0.;

      case 'G':
        return 0.;

      case 'H':
        return 1.;

      case 'I':
        return 0.;

      case 'L':
        return 0.;

      case 'K':
        return 1.;

      case 'M':
        return 0.;

      case 'F':
        return 0.;

      case 'P':
        return 0.;

      case 'S':
        return 0.;

      case 'T':
        return 0.;

      case 'W':
        return 0.;

      case 'Y':
        return 0.;

      case 'V':
        return 0.;

      default:
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Unknown amino acid one-letter-code", String(aa));
      }
    }

    //SUEM840101    0.883  AURR980114    0.875  AURR980113    0.849<br>
    //PTIO830101    0.826  KANM800103    0.823  QIAN880107    0.814<br>
    //QIAN880106    0.810  MAXF760101    0.810  AURR980109    0.802
    //I    A/L     R/K     N/M     D/F     C/P     Q/S     E/T     G/W     H/Y     I/V
    //1.08    1.05    0.85    0.85    0.95    0.95    1.15    0.55    1.00    1.05
    //1.25    1.15    1.15    1.10    0.71    0.75    0.75    1.10    1.10    0.95

    /**
      @brief Helix-coil equilibrium constant (Finkelstein-Ptitsyn, 1977)

      LIT:2004052b PMID:843599<br>
      Finkelstein, A.V. and Ptitsyn, O.B.<br>
      Theory of protein molecule self-organization. II. A comparison of calculated
      thermodynamic parameters of local secondary structures with experiments<br>
      Biopolymers 16, 497-524 (1977) (Pro 0.096)

      @exception InvalidValue is thrown if an undefined one-letter-code is used
    */
    static double getFINA770101(char aa)
    {
      switch (aa)
      {
      case 'A':
        return 1.08;

      case 'R':
        return 1.05;

      case 'N':
        return 0.85;

      case 'D':
        return 0.85;

      case 'C':
        return 0.95;

      case 'Q':
        return 0.95;

      case 'E':
        return 1.15;

      case 'G':
        return 0.55;

      case 'H':
        return 1.00;

      case 'I':
        return 1.05;

      case 'L':
        return 1.25;

      case 'K':
        return 1.15;

      case 'M':
        return 1.15;

      case 'F':
        return 1.10;

      case 'P':
        return 0.71;

      case 'S':
        return 0.75;

      case 'T':
        return 0.75;

      case 'W':
        return 1.10;

      case 'Y':
        return 1.10;

      case 'V':
        return 0.95;

      default:
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Unknown amino acid one-letter-code", String(aa));
      }
    }

    //ARGP820103    0.961  KYTJ820101    0.803  JURD980101    0.802
    //I    A/L     R/K     N/M     D/F     C/P     Q/S     E/T     G/W     H/Y     I/V
    //1.18    0.20    0.23    0.05    1.89    0.72    0.11    0.49    0.31    1.45
    //3.23    0.06    2.67    1.96    0.76    0.97    0.84    0.77    0.39    1.08

    /**
      @brief Signal sequence helical potential (Argos et al., 1982)

      LIT:0901079b PMID:7151796<br>
      Argos, P., Rao, J.K.M. and Hargrave, P.A.<br>
      Structural prediction of membrane-bound proteins<br>
      Eur. J. Biochem. 128, 565-575 (1982)

      @exception InvalidValue is thrown if an undefined one-letter-code is used
    */
    static double getARGP820102(char aa)
    {
      switch (aa)
      {
      case 'A':
        return 1.18;

      case 'R':
        return 0.20;

      case 'N':
        return 0.23;

      case 'D':
        return 0.05;

      case 'C':
        return 1.89;

      case 'Q':
        return 0.72;

      case 'E':
        return 0.11;

      case 'G':
        return 0.49;

      case 'H':
        return 0.31;

      case 'I':
        return 1.45;

      case 'L':
        return 3.23;

      case 'K':
        return 0.06;

      case 'M':
        return 2.67;

      case 'F':
        return 1.96;

      case 'P':
        return 0.76;

      case 'S':
        return 0.97;

      case 'T':
        return 0.84;

      case 'W':
        return 0.77;

      case 'Y':
        return 0.39;

      case 'V':
        return 1.08;

      default:
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Unknown amino acid one-letter-code", String(aa));
      }
    }

    /**
      @brief Calculates an estimated gas-phase basicity for an amino acid sequence at a given temperature

      Energy level E at each protonation site i is -GB(i) fractional proton population of a microstate k is <br>
      P_k = exp (- E_k/(RT)) / ( sum_i exp (- E_i/(RT))) <br>
      The apparent proton association constant K_app: K_app = sum_i GB(i)/(RT)<br>
      Then the apparent GB is GB_app^ion = R * T * ln(K_app)

      @exception InvalidValue is thrown if an undefined one-letter-code is used
    */
    static double calculateGB(const AASequence& seq, double T = 500.0)
    {

      double R = Constants::GAS_CONSTANT / 1000.0; // ideal gas constant in kj/(K*mol)

      char left = '>';
      char right;

      double k_app = 0.0; // apparent proton association constant

      // energy level E at each protonation site i is -GB(i)
      // fractional proton population of a microstate k is
      // P_k = exp (- E_k/(RT)) / ( sum_i exp (- E_i/(RT)))
      // the apparent proton association constant k_app:
      // k_app = sum_i GB(i)/(RT)
      // then the apparent GB is GB_app^ion = R * T * ln(k_app)
      for (Size i = 0; i <= seq.size(); i++)
      {
        // aa left to current one
        if (i > 0)
        {
          Residue leftchar = seq[i - 1];
          left = leftchar.getOneLetterCode()[0];
        }

        // aa right to current one
        if (i == seq.size())
        {
          right = '<';
        }
        else
        {
          Residue rightchar = seq[i];
          right = rightchar.getOneLetterCode()[0];
        }
        double contrib = exp((GBleft_(left) + GBdeltaright_(right)) / (R * T));
        if (i > 0 && i < seq.size())
        {
          contrib += exp(GBsidechain_(right) / (R * T));
        }
        k_app += contrib;
      }
      // calculate apparent GB
      return R * T * log(k_app) / log(2.0);
    }

protected:

    /**
      @brief Calculates part of the gas-phase basicity

      For a detailed description see @ref calculateGB(const AASequence&, double) .

      @exception InvalidValue is thrown if an undefined one-letter-code is used
    */
    static double GBsidechain_(char aa)
    {
      switch (aa)
      {
      case 'A':
        return 0.0;

      case 'C':
        return 0.0;

      case 'D':
        return 784.0;

      case 'E':
        return 790.0;

      case 'F':
        return 0.0;

      case 'G':
        return 0.0;

      case 'H':
        return 927.84;

      case 'I':
        return 0.0;

      case 'K':
        return 926.74;

      case 'L':
        return 0.0;

      case 'M':
        return 830.0;

      case 'N':
        return 864.94;

      case 'P':
        return 0.0;

      case 'Q':
        return 865.25;

      case 'R':
        return 1000.0;

      case 'S':
        return 775.0;

      case 'T':
        return 780.0;

      case 'V':
        return 0.0;

      case 'W':
        return 909.53;

      case 'Y':
        return 790.0;

      default:
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Unknown amino acid one-letter-code", String(aa));
      }
    }

    /**
      @brief Calculates part of the gas-phase basicity

      For a detailed description see @ref calculateGB(const AASequence&, double) .

      @exception InvalidValue is thrown if an undefined one-letter-code is used
    */
    static double GBleft_(char aa)
    {
      switch (aa)
      {
      case 'A':
        return 881.82;

      case 'C':
        return 881.15;

      case 'D':
        return 880.02;

      case 'E':
        return 880.10;

      case 'F':
        return 881.08;

      case 'G':
        return 881.17;

      case 'H':
        return 881.27;

      case 'I':
        return 880.99;

      case 'K':
        return 880.06;

      case 'L':
        return 881.88;

      case 'M':
        return 881.38;

      case 'N':
        return 881.18;

      case 'P':
        return 881.25;

      case 'Q':
        return 881.50;

      case 'R':
        return 882.98;

      case 'S':
        return 881.08;

      case 'T':
        return 881.14;

      case 'V':
        return 881.17;

      case 'W':
        return 881.31;

      case 'Y':
        return 881.20;

      case '>': //NH2
        return 916.84;

      default:
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Unknown amino acid one-letter-code", String(aa));

      }
    }

    /**
      @brief Calculates part of the gas-phase basicity

      For a detailed description see @ref calculateGB(const AASequence&, double) .

      @exception InvalidValue is thrown if an undefined one-letter-code is used
    */
    static double GBdeltaright_(char aa)
    {
      switch (aa)
      {
      case 'A':
        return 0.0;

      case 'C':
        return -0.69;

      case 'D':
        return -0.63;

      case 'E':
        return -0.39;

      case 'F':
        return 0.03;

      case 'G':
        return 0.92;

      case 'H':
        return -0.19;

      case 'I':
        return -1.17;

      case 'K':
        return -0.71;

      case 'L':
        return -0.09;

      case 'M':
        return 0.30;

      case 'N':
        return 1.56;

      case 'P':
        return 11.75;

      case 'Q':
        return 4.10;

      case 'R':
        return 6.28;

      case 'S':
        return 0.98;

      case 'T':
        return 1.21;

      case 'V':
        return -0.90;

      case 'W':
        return 0.10;

      case 'Y':
        return -0.38;

      case '<': //COOH
        return -95.82;

      default:
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Unknown amino acid one-letter-code", String(aa));
      }
    }
  };

}
