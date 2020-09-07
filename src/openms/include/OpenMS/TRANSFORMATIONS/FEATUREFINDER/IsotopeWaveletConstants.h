// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
// $Maintainer: Timo Sachsenberg$
// $Authors: Rene Hussong$
// --------------------------------------------------------------------------

#pragma once

namespace OpenMS
{
  namespace Constants
  {
#undef OPENMS_DEBUG_ISOTOPE_WAVELET

    const unsigned int DEFAULT_NUM_OF_INTERPOLATION_POINTS = 3;

    const double MASS_EPSILON = 1e-4f;

    const double MARR_WAVELET_CUTOFF = 4.f;

    const double PEPTIDE_MASS_RULE_FACTOR = 0.000507f;
    const   double PEPTIDE_MASS_RULE_BOUND =    1. / PEPTIDE_MASS_RULE_FACTOR;
    const double PEPTIDE_MASS_RULE_THEO_PPM_BOUND = 200;

    //exact
    const double IW_NEUTRON_MASS = 1.00866491578f;
    const double IW_HALF_NEUTRON_MASS = 0.5043325f;
    const double IW_QUARTER_NEUTRON_MASS = 0.252166228f;
    const double WAVELET_PERIODICITY = 6.229209734f;

    //according to Horn et al. (2000)
    /*const double IW_NEUTRON_MASS = 1.00235f;
    const double IW_HALF_NEUTRON_MASS = 0.501175f;
    const double IW_QUARTER_NEUTRON_MASS = 0.2505875f;
    const double WAVELET_PERIODICITY = 6.268454439f;*/


    const double ONEOLOG2E = 0.6931471806f;

    const double IW_PROTON_MASS = 1.00727646688f;

    //Linear Fit (standard)
    const double LAMBDA_L_0 = 0.120398590399013419f;
    const double LAMBDA_L_1 = 0.635926795694698589e-3f;

    const double CUT_LAMBDA_Q_0_A = 1.9498e+00f;
    const double CUT_LAMBDA_Q_0_B = 2.4244e-03f;
    const double CUT_LAMBDA_Q_0_C = -2.4183e-07f;
    const double CUT_LAMBDA_Q_1_A = 3.6870e+00f;
    const double CUT_LAMBDA_Q_1_B = 1.1561e-03f;
    const double CUT_LAMBDA_Q_1_C = -1.0329e-08f;
    const double CUT_LAMBDA_L_2_A = 5.7661e+00f;
    const double CUT_LAMBDA_L_2_B = 8.6301e-04f;
    const double CUT_LAMBDA_BREAK_0_1 = 2739.4f;
    const double CUT_LAMBDA_BREAK_1_2 = 1.4187e+04f;

    const int SHIFT23 = (1 << 23);
    const double SHIFT23_00 = (1.0 / (1 << 23));
    const double LOG_CONST = 0.346607f;
    const double POW_CONST = 0.33971f;
  }
}

