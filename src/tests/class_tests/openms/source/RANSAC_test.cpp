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
// $Maintainer: George Rosenberger $
// $Authors: George Rosenberger, Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/MATH/MISC/RANSAC.h>
///////////////////////////

using namespace std;
using namespace OpenMS;

///////////////////////////

START_TEST(MRMRTNormalizer, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

class RANSAC_test : public Math::RANSAC
{
  public:

    static std::pair<double, double > llsm_fit_(std::vector<std::pair<double, double> >& pairs)
    {
      return Math::RANSAC::llsm_fit_(pairs);
    }

    static double llsm_rss_(std::vector<std::pair<double, double> >& pairs, std::pair<double, double >& coefficients  )
    {
      return Math::RANSAC::llsm_rss_(pairs, coefficients);
    }

    static std::vector<std::pair<double, double> > llsm_rss_inliers_(std::vector<std::pair<double, double> >&   pairs, std::pair<double, double >& coefficients, double max_threshold)
    {
      return Math::RANSAC::llsm_rss_inliers_(pairs, coefficients, max_threshold);
    }

};

// no constructor / destructor of static class

// RANSAC() 
// ~RANSAC() 
//

START_SECTION((static std::pair<double, double > llsm_fit_(std::vector<std::pair<double, double> >& pairs)))
{
  std::vector<std::pair<double, double> > test_pairs;

  test_pairs.push_back(std::make_pair(7.66217066e+00, 3.32871078e+02));
  test_pairs.push_back(std::make_pair(1.88986378e+01, 8.41782838e+02));
  test_pairs.push_back(std::make_pair(1.43387751e+01, 6.48336013e+02));
  test_pairs.push_back(std::make_pair(1.04946477e+01, 5.30115032e+02));
  test_pairs.push_back(std::make_pair(2.40052860e+00, 1.36793947e+02));
  test_pairs.push_back(std::make_pair(2.65925164e+00, 1.38532208e+02));
  test_pairs.push_back(std::make_pair(7.00156815e+00, 3.03487855e+01));
  test_pairs.push_back(std::make_pair(1.76671412e+01, 7.67575677e+02));
  test_pairs.push_back(std::make_pair(1.02592601e+01, 5.32449429e+02));
  test_pairs.push_back(std::make_pair(1.29020672e+01,-1.74450591e+01));
  test_pairs.push_back(std::make_pair(2.66076055e-02, 1.78205080e+01));
  test_pairs.push_back(std::make_pair(1.87212750e+01, 8.59152499e+02));
  test_pairs.push_back(std::make_pair(1.81219758e+01,-5.79165989e-01));
  test_pairs.push_back(std::make_pair(5.27778174e+00, 1.88005119e+02));
  test_pairs.push_back(std::make_pair(4.56777946e+00, 1.61530045e+02));
  test_pairs.push_back(std::make_pair(2.82887267e+00, 1.64411907e+02));
  test_pairs.push_back(std::make_pair(5.77563248e+00, 2.69781852e+02));
  test_pairs.push_back(std::make_pair(1.08263921e+01, 4.65275655e+02));
  test_pairs.push_back(std::make_pair(9.61444550e+00, 3.82697907e+02));
  test_pairs.push_back(std::make_pair(5.34540857e+00, 2.56156813e+02));

  std::pair<double, double > coeff = RANSAC_test::llsm_fit_(test_pairs);
  TEST_REAL_SIMILAR( coeff.first, 46.03865245);
  TEST_REAL_SIMILAR( coeff.second, 31.20358812);

  double rss = RANSAC_test::llsm_rss_(test_pairs, coeff);
  TEST_REAL_SIMILAR( rss, 864089.67832345);

  std::vector<std::pair<double, double> > new_test_pairs;

  new_test_pairs.push_back(std::make_pair(1.20513989e+01, 5.42172984e+02));
  new_test_pairs.push_back(std::make_pair(1.68354224e+00, 1.23674095e+02));
  new_test_pairs.push_back(std::make_pair(4.64668635e+00, 2.61350113e+02));
  new_test_pairs.push_back(std::make_pair(8.13976269e+00, 3.24462812e+02));
  new_test_pairs.push_back(std::make_pair(1.04776397e+01, 4.04452477e+02));
  new_test_pairs.push_back(std::make_pair(1.56315091e+01, 6.95756737e+02));
  new_test_pairs.push_back(std::make_pair(1.27266524e+01, 6.53571377e+01));
  new_test_pairs.push_back(std::make_pair(1.33784812e+01, 3.03064682e+01));
  new_test_pairs.push_back(std::make_pair(9.73484306e+00,-1.55933991e+00));
  new_test_pairs.push_back(std::make_pair(1.29040386e+00, 4.19535249e+01));
  new_test_pairs.push_back(std::make_pair(1.36889336e+01, 5.37472495e+02));
  new_test_pairs.push_back(std::make_pair(3.37465643e+00, 1.52514434e+02));
  new_test_pairs.push_back(std::make_pair(2.86567552e+00, 5.62442618e+01));
  new_test_pairs.push_back(std::make_pair(1.63579656e+01, 8.41451166e+02));
  new_test_pairs.push_back(std::make_pair(2.01345432e+01, 8.57894838e+02));
  new_test_pairs.push_back(std::make_pair(1.62549940e+01, 7.15378774e+02));
  new_test_pairs.push_back(std::make_pair(5.79326803e+00, 2.69370208e+02));
  new_test_pairs.push_back(std::make_pair(2.04520306e+00, 8.66527618e+01));
  new_test_pairs.push_back(std::make_pair(1.16970916e+01, 6.05836392e+02));
  new_test_pairs.push_back(std::make_pair(8.68788731e+00, 9.52993526e+00));
  new_test_pairs.push_back(std::make_pair(2.79787727e+00, 1.08213952e+02));
  new_test_pairs.push_back(std::make_pair(1.95778572e+01, 1.39196902e+02));
  new_test_pairs.push_back(std::make_pair(1.69500204e-01, 3.09473207e+01));
  new_test_pairs.push_back(std::make_pair(1.17974170e+01, 2.51798532e+01));
  new_test_pairs.push_back(std::make_pair(4.67384259e+00, 2.30870376e+02));
  new_test_pairs.push_back(std::make_pair(1.41658478e+01, 5.86317425e+02));
  new_test_pairs.push_back(std::make_pair(5.00923637e+00,-1.86559595e+01));
  new_test_pairs.push_back(std::make_pair(9.87160022e+00, 4.61676941e+02));
  new_test_pairs.push_back(std::make_pair(1.14474730e+01, 4.83241860e+02));
  new_test_pairs.push_back(std::make_pair(3.79416666e+00, 1.64038065e+02));

  std::vector<std::pair<double, double> > inliers = RANSAC_test::llsm_rss_inliers_(new_test_pairs, coeff, 7e3);
  TEST_REAL_SIMILAR( inliers[0].first, 1.68354224e+00);
  TEST_REAL_SIMILAR( inliers[1].first, 4.64668635e+00);
  TEST_REAL_SIMILAR( inliers[2].first, 8.13976269e+00);
  TEST_REAL_SIMILAR( inliers[3].first, 1.04776397e+01);
  TEST_REAL_SIMILAR( inliers[4].first, 1.29040386e+00);
  TEST_REAL_SIMILAR( inliers[5].first, 1.36889336e+01);
  TEST_REAL_SIMILAR( inliers[6].first, 3.37465643e+00);
  TEST_REAL_SIMILAR( inliers[7].first, 2.86567552e+00);
  TEST_REAL_SIMILAR( inliers[8].first, 5.79326803e+00);
  TEST_REAL_SIMILAR( inliers[9].first, 2.04520306e+00);
  TEST_REAL_SIMILAR( inliers[10].first, 2.79787727e+00);
  TEST_REAL_SIMILAR( inliers[11].first, 1.69500204e-01);
  TEST_REAL_SIMILAR( inliers[12].first, 4.67384259e+00);
  TEST_REAL_SIMILAR( inliers[13].first, 1.14474730e+01);
  TEST_REAL_SIMILAR( inliers[14].first, 3.79416666e+00);
  TEST_EQUAL( inliers.size(), 15);
}
END_SECTION

START_SECTION(static double llsm_rss_(std::vector<std::pair<double, double> >& pairs, std::pair<double, double >& coefficients  ) )
{
  // tested above in llsm_fit_
  NOT_TESTABLE
}
END_SECTION

START_SECTION((static std::vector<std::pair<double, double> > llsm_rss_inliers_(std::vector<std::pair<double, double> >&   pairs, std::pair<double, double >& coefficients, double max_threshold)))
{
  // tested above in llsm_fit_
  NOT_TESTABLE
}
END_SECTION

START_SECTION((static std::vector<std::pair<double, double> > ransac(std::vector<std::pair<double, double> >& pairs, size_t n, size_t k, double t, size_t d)))
{

/*
// Python reference implementation that was used to generate the test data: http://wiki.scipy.org/Cookbook/RANSAC

import numpy
import scipy # use numpy if scipy unavailable
import scipy.linalg # use numpy if scipy unavailable

## Copyright (c) 2004-2007, Andrew D. Straw. All rights reserved.

## Redistribution and use in source and binary forms, with or without
## modification, are permitted provided that the following conditions are
## met:

##     * Redistributions of source code must retain the above copyright
##       notice, this list of conditions and the following disclaimer.

##     * Redistributions in binary form must reproduce the above
##       copyright notice, this list of conditions and the following
##       disclaimer in the documentation and/or other materials provided
##       with the distribution.

##     * Neither the name of the Andrew D. Straw nor the names of its
##       contributors may be used to endorse or promote products derived
##       from this software without specific prior written permission.

## THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
## "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
## LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
## A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
## OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
## SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
## LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
## DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
## THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
## (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
## OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

def ransac(data,model,n,k,t,d,debug=False,return_all=False):
    """fit model parameters to data using the RANSAC algorithm
    
This implementation written from pseudocode found at
http://en.wikipedia.org/w/index.php?title=RANSAC&oldid=116358182

{{{
Given:
    data - a set of observed data points
    model - a model that can be fitted to data points
    n - the minimum number of data values required to fit the model
    k - the maximum number of iterations allowed in the algorithm
    t - a threshold value for determining when a data point fits a model
    d - the number of close data values required to assert that a model fits well to data
Return:
    bestfit - model parameters which best fit the data (or nil if no good model is found)
iterations = 0
bestfit = nil
besterr = something really large
while iterations < k {
    maybeinliers = n randomly selected values from data
    maybemodel = model parameters fitted to maybeinliers
    alsoinliers = empty set
    for every point in data not in maybeinliers {
        if point fits maybemodel with an error smaller than t
             add point to alsoinliers
    }
    if the number of elements in alsoinliers is > d {
        % this implies that we may have found a good model
        % now test how good it is
        bettermodel = model parameters fitted to all points in maybeinliers and alsoinliers
        thiserr = a measure of how well model fits these points
        if thiserr < besterr {
            bestfit = bettermodel
            besterr = thiserr
        }
    }
    increment iterations
}
return bestfit
}}}
"""
    iterations = 0
    bestfit = None
    besterr = numpy.inf
    best_inlier_idxs = None
    while iterations < k:
        maybe_idxs, test_idxs = random_partition(n,data.shape[0])
        maybeinliers = data[maybe_idxs,:]
        # OPENMS_TEST: print all maybeinliers
        # print "maybeinliers"
        # print maybeinliers
        # print "end maybeinliers"
        test_points = data[test_idxs]
        maybemodel = model.fit(maybeinliers)
        test_err = model.get_error( test_points, maybemodel)
        also_idxs = test_idxs[test_err < t] # select indices of rows with accepted points
        alsoinliers = data[also_idxs,:]
        # OPENMS_TEST: print alsoinliers
        # print alsoin
        if debug:
            print 'test_err.min()',test_err.min()
            print 'test_err.max()',test_err.max()
            print 'numpy.mean(test_err)',numpy.mean(test_err)
            print 'iteration %d:len(alsoinliers) = %d'%(
                iterations,len(alsoinliers))
        if len(alsoinliers) > d:
            betterdata = numpy.concatenate( (maybeinliers, alsoinliers) )
            bettermodel = model.fit(betterdata)
            better_errs = model.get_error( betterdata, bettermodel)
            thiserr = numpy.mean( better_errs )
            if thiserr < besterr:
                bestfit = bettermodel
                besterr = thiserr
                best_inlier_idxs = numpy.concatenate( (maybe_idxs, also_idxs) )
        iterations+=1
    if bestfit is None:
        raise ValueError("did not meet fit acceptance criteria")
    if return_all:
        return bestfit, {'inliers':best_inlier_idxs}
    else:
        return bestfit

def random_partition(n,n_data):
    """return n random rows of data (and also the other len(data)-n rows)"""
    all_idxs = numpy.arange( n_data )

    # OPENMS_TEST: exclude random component
    #numpy.random.shuffle(all_idxs)
    idxs1 = all_idxs[:n]
    idxs2 = all_idxs[n:]
    return idxs1, idxs2

class LinearLeastSquaresModel:
    """linear system solved using linear least squares

    This class serves as an example that fulfills the model interface
    needed by the ransac() function.
    
    """
    def __init__(self,input_columns,output_columns,debug=False):
        self.input_columns = input_columns
        self.output_columns = output_columns
        self.debug = debug
    def fit(self, data):
        A = numpy.vstack([data[:,i] for i in self.input_columns]).T
        B = numpy.vstack([data[:,i] for i in self.output_columns]).T

        # OPENMS_TEST: make linear regression compatible
        # lstsq needs correct weights!
        W = numpy.vstack(numpy.ones(len(A)))
        A = numpy.hstack((A,W))

        x,resids,rank,s = scipy.linalg.lstsq(A,B)

        # OPENMS_TEST: print coefficients & rss
        #print x
        #print resids
        return x
    def get_error( self, data, model):
        A = numpy.vstack([data[:,i] for i in self.input_columns]).T
        B = numpy.vstack([data[:,i] for i in self.output_columns]).T

        # OPENMS_TEST: make linear regression compatible
        # lstsq needs correct weights!
        W = numpy.vstack(numpy.ones(len(A)))
        A = numpy.hstack((A,W))

        B_fit = scipy.dot(A,model)
        err_per_point = numpy.sum((B-B_fit)**2,axis=1) # sum squared error per row
        return err_per_point
        
def test():
    # generate perfect input data

    # Fix seed
    numpy.random.seed(42)

    n_samples = 50
    n_inputs = 1
    n_outputs = 1
    A_exact = 20*numpy.random.random((n_samples,n_inputs) )
    perfect_fit = 60*numpy.random.normal(size=(n_inputs,n_outputs) ) # the model
    B_exact = scipy.dot(A_exact,perfect_fit)
    assert B_exact.shape == (n_samples,n_outputs)

    # add a little gaussian noise (linear least squares alone should handle this well)
    A_noisy = A_exact + numpy.random.normal(size=A_exact.shape )
    B_noisy = B_exact + numpy.random.normal(size=B_exact.shape )

    if 1:
        # add some outliers
        n_outliers = 10
        all_idxs = numpy.arange( A_noisy.shape[0] )
        numpy.random.shuffle(all_idxs)
        outlier_idxs = all_idxs[:n_outliers]
        non_outlier_idxs = all_idxs[n_outliers:]
        A_noisy[outlier_idxs] =  20*numpy.random.random((n_outliers,n_inputs) )
        B_noisy[outlier_idxs] = 50*numpy.random.normal(size=(n_outliers,n_outputs) )

    # setup model

    all_data = numpy.hstack( (A_noisy,B_noisy) )
    input_columns = range(n_inputs) # the first columns of the array
    output_columns = [n_inputs+i for i in range(n_outputs)] # the last columns of the array
    debug = False
    model = LinearLeastSquaresModel(input_columns,output_columns,debug=debug)

    linear_fit,resids,rank,s = scipy.linalg.lstsq(all_data[:,input_columns],
                                                  all_data[:,output_columns])

    # OPENMS_TEST: print input data
    #print all_data

    # run RANSAC algorithm
    ransac_fit, ransac_data = ransac(all_data,model,
                                     20, 1, 7e3, 10, # misc. parameters
                                     debug=debug,return_all=True)


    # OPENMS_TEST: print result data
    #print all_data[ransac_data['inliers']]

    if 0:
        import pylab

        sort_idxs = numpy.argsort(A_exact[:,0])
        A_col0_sorted = A_exact[sort_idxs] # maintain as rank-2 array

        if 1:
            pylab.plot( A_noisy[:,0], B_noisy[:,0], 'k.', label='data' )
            pylab.plot( A_noisy[ransac_data['inliers'],0], B_noisy[ransac_data['inliers'],0], 'bx', label='RANSAC data' )
        else:
            pylab.plot( A_noisy[non_outlier_idxs,0], B_noisy[non_outlier_idxs,0], 'k.', label='noisy data' )
            pylab.plot( A_noisy[outlier_idxs,0], B_noisy[outlier_idxs,0], 'r.', label='outlier data' )
        pylab.plot( A_col0_sorted[:,0],
                    numpy.dot(A_col0_sorted,ransac_fit)[:,0],
                    label='RANSAC fit' )
        pylab.plot( A_col0_sorted[:,0],
                    numpy.dot(A_col0_sorted,perfect_fit)[:,0],
                    label='exact system' )
        pylab.plot( A_col0_sorted[:,0],
                    numpy.dot(A_col0_sorted,linear_fit)[:,0],
                    label='linear fit' )
        pylab.legend()
        pylab.show()

if __name__=='__main__':
    test()
    

  */

  std::vector<std::pair<double, double> > test_pairs, test_pairs_out;

  test_pairs.push_back(std::make_pair(7.66217066e+00, 3.32871078e+02));
  test_pairs.push_back(std::make_pair(1.88986378e+01, 8.41782838e+02));
  test_pairs.push_back(std::make_pair(1.43387751e+01, 6.48336013e+02));
  test_pairs.push_back(std::make_pair(1.04946477e+01, 5.30115032e+02));
  test_pairs.push_back(std::make_pair(2.40052860e+00, 1.36793947e+02));
  test_pairs.push_back(std::make_pair(2.65925164e+00, 1.38532208e+02));
  test_pairs.push_back(std::make_pair(7.00156815e+00, 3.03487855e+01));
  test_pairs.push_back(std::make_pair(1.76671412e+01, 7.67575677e+02));
  test_pairs.push_back(std::make_pair(1.02592601e+01, 5.32449429e+02));
  test_pairs.push_back(std::make_pair(1.29020672e+01,-1.74450591e+01));
  test_pairs.push_back(std::make_pair(2.66076055e-02, 1.78205080e+01));
  test_pairs.push_back(std::make_pair(1.87212750e+01, 8.59152499e+02));
  test_pairs.push_back(std::make_pair(1.81219758e+01,-5.79165989e-01));
  test_pairs.push_back(std::make_pair(5.27778174e+00, 1.88005119e+02));
  test_pairs.push_back(std::make_pair(4.56777946e+00, 1.61530045e+02));
  test_pairs.push_back(std::make_pair(2.82887267e+00, 1.64411907e+02));
  test_pairs.push_back(std::make_pair(5.77563248e+00, 2.69781852e+02));
  test_pairs.push_back(std::make_pair(1.08263921e+01, 4.65275655e+02));
  test_pairs.push_back(std::make_pair(9.61444550e+00, 3.82697907e+02));
  test_pairs.push_back(std::make_pair(5.34540857e+00, 2.56156813e+02));
  test_pairs.push_back(std::make_pair(1.20513989e+01, 5.42172984e+02));
  test_pairs.push_back(std::make_pair(1.68354224e+00, 1.23674095e+02));
  test_pairs.push_back(std::make_pair(4.64668635e+00, 2.61350113e+02));
  test_pairs.push_back(std::make_pair(8.13976269e+00, 3.24462812e+02));
  test_pairs.push_back(std::make_pair(1.04776397e+01, 4.04452477e+02));
  test_pairs.push_back(std::make_pair(1.56315091e+01, 6.95756737e+02));
  test_pairs.push_back(std::make_pair(1.27266524e+01, 6.53571377e+01));
  test_pairs.push_back(std::make_pair(1.33784812e+01, 3.03064682e+01));
  test_pairs.push_back(std::make_pair(9.73484306e+00,-1.55933991e+00));
  test_pairs.push_back(std::make_pair(1.29040386e+00, 4.19535249e+01));
  test_pairs.push_back(std::make_pair(1.36889336e+01, 5.37472495e+02));
  test_pairs.push_back(std::make_pair(3.37465643e+00, 1.52514434e+02));
  test_pairs.push_back(std::make_pair(2.86567552e+00, 5.62442618e+01));
  test_pairs.push_back(std::make_pair(1.63579656e+01, 8.41451166e+02));
  test_pairs.push_back(std::make_pair(2.01345432e+01, 8.57894838e+02));
  test_pairs.push_back(std::make_pair(1.62549940e+01, 7.15378774e+02));
  test_pairs.push_back(std::make_pair(5.79326803e+00, 2.69370208e+02));
  test_pairs.push_back(std::make_pair(2.04520306e+00, 8.66527618e+01));
  test_pairs.push_back(std::make_pair(1.16970916e+01, 6.05836392e+02));
  test_pairs.push_back(std::make_pair(8.68788731e+00, 9.52993526e+00));
  test_pairs.push_back(std::make_pair(2.79787727e+00, 1.08213952e+02));
  test_pairs.push_back(std::make_pair(1.95778572e+01, 1.39196902e+02));
  test_pairs.push_back(std::make_pair(1.69500204e-01, 3.09473207e+01));
  test_pairs.push_back(std::make_pair(1.17974170e+01, 2.51798532e+01));
  test_pairs.push_back(std::make_pair(4.67384259e+00, 2.30870376e+02));
  test_pairs.push_back(std::make_pair(1.41658478e+01, 5.86317425e+02));
  test_pairs.push_back(std::make_pair(5.00923637e+00,-1.86559595e+01));
  test_pairs.push_back(std::make_pair(9.87160022e+00, 4.61676941e+02));
  test_pairs.push_back(std::make_pair(1.14474730e+01, 4.83241860e+02));
  test_pairs.push_back(std::make_pair(3.79416666e+00, 1.64038065e+02));

  test_pairs_out = Math::RANSAC::ransac(test_pairs, 20, 1, 7e3, 10, true);

  TEST_REAL_SIMILAR( test_pairs_out[0].first, 7.66217066e+00);
  TEST_REAL_SIMILAR( test_pairs_out[1].first, 1.88986378e+01);
  TEST_REAL_SIMILAR( test_pairs_out[2].first, 1.43387751e+01);
  TEST_REAL_SIMILAR( test_pairs_out[3].first, 1.04946477e+01);
  TEST_REAL_SIMILAR( test_pairs_out[4].first, 2.40052860e+00);
  TEST_REAL_SIMILAR( test_pairs_out[5].first, 2.65925164e+00);
  TEST_REAL_SIMILAR( test_pairs_out[6].first, 7.00156815e+00);
  TEST_REAL_SIMILAR( test_pairs_out[7].first, 1.76671412e+01);
  TEST_REAL_SIMILAR( test_pairs_out[8].first, 1.02592601e+01);
  TEST_REAL_SIMILAR( test_pairs_out[9].first, 1.29020672e+01);
  TEST_REAL_SIMILAR( test_pairs_out[10].first, 2.66076055e-02);
  TEST_REAL_SIMILAR( test_pairs_out[11].first, 1.87212750e+01);
  TEST_REAL_SIMILAR( test_pairs_out[12].first, 1.81219758e+01);
  TEST_REAL_SIMILAR( test_pairs_out[13].first, 5.27778174e+00);
  TEST_REAL_SIMILAR( test_pairs_out[14].first, 4.56777946e+00);
  TEST_REAL_SIMILAR( test_pairs_out[15].first, 2.82887267e+00);
  TEST_REAL_SIMILAR( test_pairs_out[16].first, 5.77563248e+00);
  TEST_REAL_SIMILAR( test_pairs_out[17].first, 1.08263921e+01);
  TEST_REAL_SIMILAR( test_pairs_out[18].first, 9.61444550e+00);
  TEST_REAL_SIMILAR( test_pairs_out[19].first, 5.34540857e+00);
  TEST_REAL_SIMILAR( test_pairs_out[20].first, 1.68354224e+00);
  TEST_REAL_SIMILAR( test_pairs_out[21].first, 4.64668635e+00);
  TEST_REAL_SIMILAR( test_pairs_out[22].first, 8.13976269e+00);
  TEST_REAL_SIMILAR( test_pairs_out[23].first, 1.04776397e+01);
  TEST_REAL_SIMILAR( test_pairs_out[24].first, 1.29040386e+00);
  TEST_REAL_SIMILAR( test_pairs_out[25].first, 1.36889336e+01);
  TEST_REAL_SIMILAR( test_pairs_out[26].first, 3.37465643e+00);
  TEST_REAL_SIMILAR( test_pairs_out[27].first, 2.86567552e+00);
  TEST_REAL_SIMILAR( test_pairs_out[28].first, 5.79326803e+00);
  TEST_REAL_SIMILAR( test_pairs_out[29].first, 2.04520306e+00);
  TEST_REAL_SIMILAR( test_pairs_out[30].first, 2.79787727e+00);
  TEST_REAL_SIMILAR( test_pairs_out[31].first, 1.69500204e-01);
  TEST_REAL_SIMILAR( test_pairs_out[32].first, 4.67384259e+00);
  TEST_REAL_SIMILAR( test_pairs_out[33].first, 1.14474730e+01);
  TEST_REAL_SIMILAR( test_pairs_out[34].first, 3.79416666e+00);
  TEST_EQUAL( test_pairs_out.size(), 35);

}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

