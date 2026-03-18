// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/PROCESSING/BASELINE/MorphologicalFilter.h>

#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/Peak2D.h>

#include <fstream>
///////////////////////////

namespace OpenMS
{
  template < typename ValueT >
  struct SimpleTopHat
  {
    static void erosion( const std::vector<ValueT> & input, std::vector<ValueT> & output, const UInt struc_elem_length )
    {
      const Int size = Int(input.size());
      const Int struc_elem_half = struc_elem_length / 2; // yes integer division
      output.clear();
      output.resize(size);
      for ( Int index = 0; index < size; ++ index )
      {
        Int begin = std::max( 0,    index - struc_elem_half );
        Int end   = std::min( size - 1, index + struc_elem_half );
        ValueT value = std::numeric_limits<ValueT>::max();
        for ( Int i = begin; i <= end; ++i )
        {
          if ( value > input[i] ) value = input[i];
        }
        output[index] = value;
      }
      return;
    }

    static void dilation( const std::vector<ValueT> & input, std::vector<ValueT> & output, const UInt struc_elem_length )
    {
      const Int size = Int(input.size());
      const Int struc_elem_half = struc_elem_length / 2; // yes integer division
      output.clear();
      output.resize(size);
      for ( Int index = 0; index < size; ++ index )
      {
        Int begin = std::max( 0,    index - struc_elem_half );
        Int end   = std::min( size - 1, index + struc_elem_half );
        ValueT value = - std::numeric_limits<ValueT>::max();
        for ( Int i = begin; i <= end; ++i )
        {
          if ( value < input[i] ) value = input[i];
        }
        output[index] = value;
      }
      return;
    }

    static void gradient( const std::vector<ValueT> & input, std::vector<ValueT> & output, const UInt struc_elem_length )
    {
      const Int size = Int(input.size());
      output.clear();
      output.resize(size);
      std::vector<ValueT> dilation;
      std::vector<ValueT> erosion;
      SimpleTopHat::erosion(input,erosion,struc_elem_length);
      SimpleTopHat::dilation(input,dilation,struc_elem_length);
      for ( Int index = 0; index < size; ++ index )
      {
        output[index] = dilation[index] - erosion[index];
      }
      return;
    }

    static void tophat( const std::vector<ValueT> & input, std::vector<ValueT> & output, const UInt struc_elem_length )
    {
      const Int size = Int(input.size());
      std::vector<ValueT> opening;
      erosion(input,output,struc_elem_length);
      dilation(output,opening,struc_elem_length);
      for ( Int index = 0; index < size; ++ index )
      {
        output[index] = input[index] - opening[index];
      }
      return;
    }

    static void bothat( const std::vector<ValueT> & input, std::vector<ValueT> & output, const UInt struc_elem_length )
    {
      const Int size = Int(input.size());
      std::vector<ValueT> closing;
      dilation(input,output,struc_elem_length);
      erosion(output,closing,struc_elem_length);
      for ( Int index = 0; index < size; ++ index )
      {
        output[index] = input[index] - closing[index];
      }
      return;
    }
  };

}

///////////////////////////

START_TEST(MorphologicalFilter, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////


using namespace OpenMS;

Int data[] = { 1, 2, 3, -2, 0, 1, 0, 0, 0, 1, 1, 1, 4, 5, 6, 4, 3, 2, 2, 5, 5, 6, 6, 1, 0, 0, -1, 0, 0, 3, -2, -3, -1, 1, 1, 1, 1, 4, 6, 2 };
const UInt data_size = sizeof(data)/sizeof(*data);

UInt struc_elem_length = 3;
typedef SimpleTopHat<Int> STH;
std::vector<Int> input;
input.reserve(data_size);
for ( UInt i = 0; i != data_size; ++i ) input.push_back(data[i]);
std::vector<Int> erosion;
STH::erosion(input,erosion,struc_elem_length);
std::vector<Int> dilation;
STH::dilation(input,dilation,struc_elem_length);
std::vector<Int> gradient;
STH::gradient(input,gradient,struc_elem_length);
std::vector<Int> opening;
STH::dilation(erosion,opening,struc_elem_length);
std::vector<Int> closing;
STH::erosion(dilation,closing,struc_elem_length);
std::vector<Int> tophat;
STH::tophat(input,tophat,struc_elem_length);
std::vector<Int> bothat;
STH::bothat(input,bothat,struc_elem_length);

START_SECTION(([EXTRA] "struct SimpleTopHat, used as reference implementation"))
{
  std::string tmpfn;
  NEW_TMP_FILE(tmpfn);
  std::ofstream tmpf(tmpfn.c_str());

  // Using a bit of macro magic to ensure that header and numbers are always in sync (gnuplot etc.)
#ifdef SIMPLETOPHATTABLE
#error "SIMPLETOPHATTABLE already #defined, oops!"
#endif
#ifdef EnTrY
#error "EnTrY already #defined, oops!"
#endif

#define SIMPLETOPHATTABLE EnTrYfAt(input,1) EnTrY(erosion,2) EnTrY(opening,3) EnTrY(dilation,4) EnTrY(closing,5) EnTrY(gradient,6) EnTrY(tophat,7) EnTrY(bothat,8)
#define EnTrY(a,b) " "#a
#define EnTrYfAt(a,b) EnTrY(a,b)
  tmpf << "#" SIMPLETOPHATTABLE << std::endl;
#undef EnTrY
#undef EnTrYfAt
#define EnTrY(a,b) a[i] << " " <<
#define EnTrYfAt(a,b) EnTrY(a,b)
  for ( UInt i = 0; i != data_size; ++i ) tmpf << SIMPLETOPHATTABLE std::endl;
  tmpf.close();
#undef EnTrY
#undef EnTrYfAt

  TEST_FILE_EQUAL(tmpfn.c_str(),OPENMS_GET_TEST_DATA_PATH("MorphologicalFilter_test_1.txt"));

  // Documentation in MorphologicalFilter uses the following gnuplot script

  // set terminal png
  // set key outside reverse Left width 2

#define EnTrY(a,b) ", '" << tmpfn << "' using "#b" w l lw 4 lt " << b << " title '"#a"'"
#define EnTrYfAt(a,b) ", '" << tmpfn << "' using "#b" w l lw 10 lt " << b << " title '"#a"'"
  std::string tmpfn2;
  NEW_TMP_FILE(tmpfn2);
  tmpf.open(tmpfn2.c_str());
  tmpf << "set title 'morphological filtering operations  (width of structuring element: " << struc_elem_length << ")'\n";
  tmpf << "set key outside reverse Left width 2\n";
  tmpf << "set output '" << tmpfn2 << ".png'\n";
  tmpf << "set terminal png size 1000,400\n";
  tmpf << "plot [] [-5:10] 0" SIMPLETOPHATTABLE "\n";
  tmpf.close();

  // Documentation in MorphologicalFilter uses the following gnuplot script
#undef SIMPLETOPHATTABLE
#define SIMPLETOPHATTABLE EnTrYfAt(input,1) EnTrY(erosion,2) EnTrY(opening,3) /*EnTrY(dilation,4)*/ /*EnTrY(closing,5)*/ /*EnTrY(gradient,6)*/ EnTrY(tophat,7) /*EnTrY(bothat,8)*/
  NEW_TMP_FILE(tmpfn2);
  tmpf.open(tmpfn2.c_str());
  tmpf << "set title 'morphological filtering operations  (width of structuring element: " << struc_elem_length << ")'\n";
  tmpf << "set key outside reverse Left width 2\n";
  tmpf << "set output '" << tmpfn2 << ".png'\n";
  tmpf << "set terminal png size 1000,400\n";
  tmpf << "plot [] [-5:10] 0" SIMPLETOPHATTABLE "\n";
  tmpf.close();

#undef EnTrY
#undef EnTrYfAt
#undef SIMPLETOPHATTABLE


}
END_SECTION

MorphologicalFilter* tophat_ptr = nullptr;
MorphologicalFilter* tophat_nullPointer = nullptr;

START_SECTION((MorphologicalFilter()))
{
  tophat_ptr = new MorphologicalFilter;
  TEST_NOT_EQUAL(tophat_ptr, tophat_nullPointer);
}
END_SECTION

START_SECTION((virtual ~MorphologicalFilter()))
{
  delete tophat_ptr;
}
END_SECTION

typedef SimpleTopHat<Peak1D::IntensityType> STHF;
std::vector<Peak1D::IntensityType> inputf;
inputf.reserve(data_size);
for ( UInt i = 0; i != data_size; ++i ) inputf.push_back(data[i]);

START_SECTION((template < typename InputIterator, typename OutputIterator > void filterRange( InputIterator input_begin, InputIterator input_end, OutputIterator output_begin)))
{

  // This test uses increasing and decreasing sequences of numbers.  This way
  // we are likely to catch all off-by-one errors. ;-) An [EXTRA] test
  // follows, which uses more realisic data.

  using Internal::intensityIteratorWrapper;

  for ( Int data_size = 0; data_size < 50; ++ data_size )
  {
    Int offset = data_size / 2;
    std::vector<Peak1D> raw;
    raw.clear();
    Peak1D peak;
    inputf.clear();
    for ( Int i = 0; i != data_size; ++i )
    {
      peak.setIntensity(i-offset);
      peak.setPos(i);
      raw.push_back(peak);
      inputf.push_back(i-offset);
    }
    std::vector<Peak1D::IntensityType> filtered;
    std::vector<Peak1D::IntensityType> simple_filtered_1;
    MorphologicalFilter mf;

    for ( Int struc_length = 3; struc_length <= 2 * data_size + 2; struc_length += 2 )
    {
      //STATUS("data_size: " << data_size);
      //STATUS("struc_elem_length: " << struc_length);
      {
        //STATUS("erosion");
        filtered.clear();
        filtered.resize(data_size);
        simple_filtered_1.clear();
        simple_filtered_1.resize(data_size);

        Param parameters;
        parameters.setValue("method","erosion");
        parameters.setValue("struc_elem_length",(double)struc_length);
        mf.setParameters(parameters);

        mf.filterRange(  intensityIteratorWrapper(raw.begin()),
                         intensityIteratorWrapper(raw.end()),
                         filtered.begin()
                       );
        STHF::erosion( inputf,
                       simple_filtered_1,
                       struc_length
                     );
        for ( Int i = 0; i != data_size; ++i )
        {
          //STATUS(i);
          TEST_REAL_SIMILAR(filtered[i],simple_filtered_1[i]);
        }
      }

      {
        //STATUS("dilation");
        filtered.clear();
        filtered.resize(data_size);
        simple_filtered_1.clear();
        simple_filtered_1.resize(data_size);

        Param parameters;
        parameters.setValue("method","dilation");
        parameters.setValue("struc_elem_length",(double)struc_length);
        mf.setParameters(parameters);

        mf.filterRange(  intensityIteratorWrapper(raw.begin()),
                         intensityIteratorWrapper(raw.end()),
                         filtered.begin()
                       );
        STHF::dilation( inputf,
                        simple_filtered_1,
                        struc_length
                      );
        for ( Int i = 0; i != data_size; ++i )
        {
          //STATUS(i);
          TEST_REAL_SIMILAR(filtered[i],simple_filtered_1[i]);
        }
      }
    }
  }

  for ( Int data_size = 0; data_size < 50; ++ data_size )
  {
    Int offset = data_size / 2;
    std::vector<Peak1D> raw;
    raw.clear();
    Peak1D peak;
    inputf.clear();
    for ( Int i = 0; i != data_size; ++i )
    {
      peak.setIntensity(offset-i);
      peak.setPos(i);
      raw.push_back(peak);
      inputf.push_back(offset-i);
    }
    std::vector<Peak1D::IntensityType> filtered;
    std::vector<Peak1D::IntensityType> simple_filtered_1;
    MorphologicalFilter mf;

    for ( Int struc_length = 3; struc_length <= 2 * data_size + 2; struc_length += 2 )
    {
      //STATUS("data_size: " << data_size);
      //STATUS("struc_elem_length: " << struc_length);
      {
        //STATUS("erosion");
        filtered.clear();
        filtered.resize(data_size);
        simple_filtered_1.clear();
        simple_filtered_1.resize(data_size);

        Param parameters;
        parameters.setValue("method","erosion");
        parameters.setValue("struc_elem_length",(double)struc_length);
        mf.setParameters(parameters);

        mf.filterRange(  intensityIteratorWrapper(raw.begin()),
                         intensityIteratorWrapper(raw.end()),
                         filtered.begin()
                       );
        STHF::erosion( inputf,
                       simple_filtered_1,
                       struc_length
                     );
        for ( Int i = 0; i != data_size; ++i )
        {
          //STATUS(i);
          TEST_REAL_SIMILAR(filtered[i],simple_filtered_1[i]);
        }
      }

      {
        //STATUS("dilation");
        filtered.clear();
        filtered.resize(data_size);
        simple_filtered_1.clear();
        simple_filtered_1.resize(data_size);

        Param parameters;
        parameters.setValue("method","dilation");
        parameters.setValue("struc_elem_length",(double)struc_length);
        mf.setParameters(parameters);

        mf.filterRange(  intensityIteratorWrapper(raw.begin()),
                         intensityIteratorWrapper(raw.end()),
                         filtered.begin()
                       );
        STHF::dilation( inputf,
                        simple_filtered_1,
                        struc_length
                      );
        for ( Int i = 0; i != data_size; ++i )
        {
          //STATUS(i);
          TEST_REAL_SIMILAR(filtered[i],simple_filtered_1[i]);
        }
      }

    }
  }

}
END_SECTION

START_SECTION([EXTRA] (template < typename InputIterator, typename OutputIterator > void filterRange( InputIterator input_begin, InputIterator input_end, OutputIterator output_begin)))
{
  using Internal::intensityIteratorWrapper;
  std::vector<Peak1D> raw;
  Peak1D peak;
  for ( UInt i = 0; i != data_size; ++i )
  {
    peak.setIntensity(data[i]);
    peak.setPos(i);
    raw.push_back(peak);
  }
  inputf.clear();
  for ( UInt i = 0; i != data_size; ++i ) inputf.push_back(data[i]);
  std::vector<Peak1D::IntensityType> filtered;
  std::vector<Peak1D::IntensityType> simple_filtered_1;
  std::vector<Peak1D::IntensityType> simple_filtered_2;
  std::vector<Peak1D::IntensityType> simple_filtered_3;
  MorphologicalFilter mf;
  for ( UInt struc_length = 3; struc_length <= 2 * data_size + 2; struc_length += 2 )
  {
    //STATUS("struc_elem_length: " << struc_length);
    {
      //STATUS("erosion");
      filtered.clear();
      filtered.resize(data_size);
      simple_filtered_1.clear();
      simple_filtered_1.resize(data_size);

      Param parameters;
      parameters.setValue("method","erosion");
      parameters.setValue("struc_elem_length",(double)struc_length);
      mf.setParameters(parameters);

      mf.filterRange( intensityIteratorWrapper(raw.begin()),
                       intensityIteratorWrapper(raw.end()),
                       filtered.begin()
                     );
      STHF::erosion(inputf,simple_filtered_1,struc_length);
      for ( UInt i = 0; i != data_size; ++i )
      {
        //STATUS(i);
        TEST_REAL_SIMILAR(filtered[i],simple_filtered_1[i]);
      }

      //STATUS("erosion_simple");
      filtered.clear();
      filtered.resize(data_size);
      simple_filtered_1.clear();
      simple_filtered_1.resize(data_size);

      parameters.setValue("method","erosion_simple");
      parameters.setValue("struc_elem_length",(double)struc_length);
      mf.setParameters(parameters);

      mf.filterRange(  intensityIteratorWrapper(raw.begin()),
                       intensityIteratorWrapper(raw.end()),
                       filtered.begin()
                     );
      STHF::erosion(inputf,simple_filtered_1,struc_length);
      for ( UInt i = 0; i != data_size; ++i )
      {
        TEST_REAL_SIMILAR(filtered[i],simple_filtered_1[i]);
      }

      //STATUS("opening");
      filtered.clear();
      filtered.resize(data_size);
      simple_filtered_2.clear();
      simple_filtered_2.resize(data_size);

      parameters.setValue("method","opening");
      parameters.setValue("struc_elem_length",(double)struc_length);
      mf.setParameters(parameters);

      mf.filterRange(  intensityIteratorWrapper(raw.begin()),
                       intensityIteratorWrapper(raw.end()),
                       filtered.begin()
                     );
       STHF::dilation(simple_filtered_1,simple_filtered_2,struc_length);
      for ( UInt i = 0; i != data_size; ++i )
      {
        TEST_REAL_SIMILAR(filtered[i],simple_filtered_2[i]);
      }

      //STATUS("tophat");
      filtered.clear();
      filtered.resize(data_size);
      simple_filtered_3.clear();
      simple_filtered_3.resize(data_size);

      parameters.setValue("method","tophat");
      parameters.setValue("struc_elem_length",(double)struc_length);
      mf.setParameters(parameters);

      mf.filterRange(  intensityIteratorWrapper(raw.begin()),
                       intensityIteratorWrapper(raw.end()),
                       filtered.begin()
                     );
      STHF::tophat(inputf,simple_filtered_3,struc_length);
      for ( UInt i = 0; i != data_size; ++i )
      {
        TEST_REAL_SIMILAR(filtered[i],simple_filtered_3[i]);
      }
    }

    {
      //STATUS("dilation");
      filtered.clear();
      filtered.resize(data_size);
      simple_filtered_1.clear();
      simple_filtered_1.resize(data_size);

      Param parameters;
      parameters.setValue("method","dilation");
      parameters.setValue("struc_elem_length",(double)struc_length);
      mf.setParameters(parameters);

      mf.filterRange(  intensityIteratorWrapper(raw.begin()),
                       intensityIteratorWrapper(raw.end()),
                       filtered.begin()
                     );
      STHF::dilation(inputf,simple_filtered_1,struc_length);
      for ( UInt i = 0; i != data_size; ++i )
      {
        TEST_REAL_SIMILAR(filtered[i],simple_filtered_1[i]);
      }

      //STATUS("dilation_simple");
      filtered.clear();
      filtered.resize(data_size);
      simple_filtered_1.clear();
      simple_filtered_1.resize(data_size);

      parameters.setValue("method","dilation_simple");
      parameters.setValue("struc_elem_length",(double)struc_length);
      mf.setParameters(parameters);

      mf.filterRange(  intensityIteratorWrapper(raw.begin()),
                       intensityIteratorWrapper(raw.end()),
                       filtered.begin()
                     );
      STHF::dilation(inputf,simple_filtered_1,struc_length);
      for ( UInt i = 0; i != data_size; ++i )
      {
        TEST_REAL_SIMILAR(filtered[i],simple_filtered_1[i]);
      }

      //STATUS("closing");
      filtered.clear();
      filtered.resize(data_size);
      simple_filtered_2.clear();
      simple_filtered_2.resize(data_size);

      parameters.setValue("method","closing");
      parameters.setValue("struc_elem_length",(double)struc_length);
      mf.setParameters(parameters);

      mf.filterRange(  intensityIteratorWrapper(raw.begin()),
                       intensityIteratorWrapper(raw.end()),
                       filtered.begin()
                     );
      STHF::erosion(simple_filtered_1,simple_filtered_2,struc_length);
      for ( UInt i = 0; i != data_size; ++i )
      {
        TEST_REAL_SIMILAR(filtered[i],simple_filtered_2[i]);
      }

      //STATUS("bothat");
      filtered.clear();
      filtered.resize(data_size);
      simple_filtered_3.clear();
      simple_filtered_3.resize(data_size);

      parameters.setValue("method","bothat");
      parameters.setValue("struc_elem_length",(double)struc_length);
      mf.setParameters(parameters);

      mf.filterRange(  intensityIteratorWrapper(raw.begin()),
                       intensityIteratorWrapper(raw.end()),
                       filtered.begin()
                     );
      STHF::bothat(inputf,simple_filtered_3,struc_length);
      for ( UInt i = 0; i != data_size; ++i )
      {
        TEST_REAL_SIMILAR(filtered[i],simple_filtered_3[i]);
      }
    }

  }
}
END_SECTION

START_SECTION((template <typename PeakType> void filter(MSSpectrum& spectrum)))
{
  MSSpectrum raw;
  Peak1D peak;
  double spacing = 0.25;
  for ( UInt i = 0; i < data_size; ++i )
  {
    peak.setIntensity(data[i]);
    peak.setPos( double(i) * spacing );
    raw.push_back(peak);
  }
  MorphologicalFilter mf;
  for ( double struc_size = .5; struc_size <= 2; struc_size += .1 )
  {
    MSSpectrum filtered(raw);

    Param parameters;
    parameters.setValue("method","dilation");
    parameters.setValue("struc_elem_length",(double)struc_size);
    parameters.setValue("struc_elem_unit","Thomson");
    mf.setParameters(parameters);

    mf.filter(filtered);
    UInt struc_size_datapoints = UInt ( ceil ( struc_size / spacing ) );
    if ( !Math::isOdd(struc_size_datapoints) ) ++struc_size_datapoints;
    STH::dilation( input, dilation, struc_size_datapoints );
    //STATUS( "struc_size: " << struc_size << "  struc_size_datapoints: " << struc_size_datapoints );
    for ( UInt i = 0; i != data_size; ++i )
    {
      STATUS("i: " << i);
      TEST_REAL_SIMILAR(filtered[i].getIntensity(),dilation[i]);
    }
  }
}
END_SECTION

START_SECTION((template <typename PeakType > void filterExperiment(MSExperiment< PeakType > &exp)))
{
  MSSpectrum raw;
  raw.setComment("Let's see if this comment is copied by the filter.");
  Peak1D peak;
  double spacing = 0.25;
  for ( UInt i = 0; i < data_size; ++i )
  {
    peak.setIntensity(data[i]);
    peak.setPos( double(i) * spacing );
    raw.push_back(peak);
  }
  MorphologicalFilter mf;
  for ( double struc_size = .5; struc_size <= 2; struc_size += .1 )
  {
    PeakMap mse_raw;
    mse_raw.addSpectrum(raw);
    mse_raw.addSpectrum(raw);
    mse_raw.addSpectrum(raw);

    Param parameters;
    parameters.setValue("method","dilation");
    parameters.setValue("struc_elem_length",(double)struc_size);
    parameters.setValue("struc_elem_unit","Thomson");

    mf.setParameters(parameters);

    mf.filterExperiment( mse_raw );
    TEST_EQUAL(mse_raw.size(),3);
    UInt struc_size_datapoints = UInt ( ceil ( struc_size / spacing ) );
    if ( !Math::isOdd(struc_size_datapoints) ) ++struc_size_datapoints;
    STH::dilation( input, dilation, struc_size_datapoints );
    //STATUS( "struc_size: " << struc_size << "  struc_size_datapoints: " << struc_size_datapoints );
    for ( UInt scan = 0; scan < 3; ++scan )
    {
      TEST_STRING_EQUAL(mse_raw[scan].getComment(),"Let's see if this comment is copied by the filter.");
      for ( UInt i = 0; i != data_size; ++i )
      {
        //STATUS("i: " << i);
        TEST_REAL_SIMILAR(mse_raw[scan][i].getIntensity(),dilation[i]);
      }
    }
  }

}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

