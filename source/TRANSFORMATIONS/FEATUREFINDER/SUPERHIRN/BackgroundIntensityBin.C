///////////////////////////////////////////////////////////////////////////
//
//  PEAK DETECTION OF FOURIER TRANSFORME MS INSTRUMENT DATA
//
//  written by Markus Mueller, markus.mueller@imsb.biol.ethz.ch
//  ( and Lukas Mueller, Lukas.Mueller@imsb.biol.ethz.ch)
//  October 2005
//  
//  Ported to OpenMS by Florian Zeller, florian.zeller@bsse.ethz.ch
//  December 2010
//
//  Group of Prof. Ruedi Aebersold, IMSB, ETH Hoenggerberg, Zurich
// 
//

#include <string>
#include <vector>
#include <map>
#include <math.h>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/ms_peak.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/BackgroundIntensityBin.h>

using namespace std;

double BackgroundIntensityBin::TR_BINS = 2.0;
double BackgroundIntensityBin::MZ_BINS = 50.0;
double BackgroundIntensityBin::INTENS_BINS = 50;
int BackgroundIntensityBin::MIN_BIN_COUNT = 1;

////////////////////////////////////////////////
// constructor for the object BackgroundIntensityBin:
BackgroundIntensityBin::BackgroundIntensityBin(){
}

////////////////////////////////////////////////
// constructor for the object BackgroundIntensityBin:
BackgroundIntensityBin::BackgroundIntensityBin(double mz, double tr){
  mzCoord = mz;
  trCoord = tr;
  zCoord = -1;
}


////////////////////////////////////////////////
// initialization of the BackgroundIntensityBin classe:
void BackgroundIntensityBin::init(){}


////////////////////////////////////////////////
// check if a peak belongs to this intenisty bin
bool BackgroundIntensityBin::checkBelonging(ms_peak* peak){
  
  // check charge state:
  if( zCoord != -1 ){
    if( peak->get_charge_state() != zCoord ){
      return false;
    }
  }
  
  // check tr:
  double tr = peak->get_retention_time();
  if( ( tr < (trCoord - TR_BINS /2.0) ) || ( tr > (trCoord + TR_BINS /2.0) ) ){
    return false;
  }
  
  double mz = peak->get_MZ();

  if( ( mz < (mzCoord - MZ_BINS /2.0) ) || ( mz > (mzCoord + MZ_BINS /2.0) ) ){
    return false;
  }
  
  addIntensity( peak->get_intensity() );
  peak = NULL;
  return true;
}

//////////////////////////////////////////////////
// add peak to BackgroundIntensityBin
void BackgroundIntensityBin::addMSPeak( ms_peak* peak ){
  addIntensity( peak->get_intensity() );
  // cout<<mzCoord<<" "<<trCoord<<endl;
  // peak->show_info();
  peak = NULL;
}



//////////////////////////////////////////////////
// add intensity to BackgroundIntensityBin
void BackgroundIntensityBin::addIntensity( double intens ){
  IntensityMap.push_back( intens );
}


// copied from simple_math
double simple_math_WEIGHTED_AVERAGE(map< double, double >* IN){
  
  double AVERAGE = 0;
  double TOT_WEIGHT = 0;  
  
  if( IN->size() > 1 ){
    
    map< double, double >::iterator START = IN->begin();  
    while( START != IN->end() ){
      TOT_WEIGHT += (*START).second;
      AVERAGE += ( (*START).first * (*START).second );
      START++;
    }
    
    return AVERAGE / TOT_WEIGHT;  
  }
  else{
    return (*IN->begin()).first;
    IN = NULL;
  }
}


//////////////////////////////////////////////////
// process collected intensities in the map
void BackgroundIntensityBin::processIntensities( ){
  
  computeIntensityHist( );
  
  // compute histogram parameters:
  // median:
  if(!IntensityHist.empty()){
    //simple_math math;
    //mean = math.WEIGHTED_AVERAGE( &IntensityHist );
    mean = simple_math_WEIGHTED_AVERAGE( &IntensityHist );
    //mean = math.MEDIAN( &IntensityHist );  
    //mean = math.WEIGHTED_MEDIAN( &IntensityHist );
    // computeHistogramGravityPoint();
  }
  else{
    mean = 0;
  }
  
}


//////////////////////////////////////////////////
// compute the gravity of the histogram:
void BackgroundIntensityBin::computeHistogramGravityPoint(){
  
  
  double start = 0;
  double mid = 0;
  double TOTAREA = 0;
  double grav = 0;
  map<double, double>::iterator P = IntensityHist.begin( ) ;
  while( P != IntensityHist.end() ){
    mid = P->first;
    double area = (mid - start) * P->second; 
    TOTAREA += area;
    grav += ( mid * area );
    P++;
  }
  mean = grav / TOTAREA;
}


//////////////////////////////////////////////////
// copmute an intensity histogram
void BackgroundIntensityBin::computeIntensityHist( ){
  
  double constraint = BackgroundIntensityBin::INTENS_BINS;
  
  // insert into the histogram map
  vector<double>::iterator P = IntensityMap.begin();
  while( P != IntensityMap.end() ){
    
    // intensity to bin:
    double intens = (*P);
    
    // find a key:
    map<double, double>::iterator F = IntensityHist.lower_bound( intens ) ;
    if( F != IntensityHist.end() ){
      
      // check this one:
      map<double, double>::iterator check = F;      
      double mainLow = fabs( check->first - intens );
      double deltaHigh = 1000000;
      if( check != IntensityHist.begin() ){
        check--;
        deltaHigh = fabs( check->first - intens );
        if( mainLow > deltaHigh ){ 
          mainLow = deltaHigh;
          F = check;
        }
      }
      if( mainLow > constraint){
        F = IntensityHist.end();
      }
      else{
        F->second += 1.0;
      }
    }  
    
    if( F == IntensityHist.end() ){
      IntensityHist.insert( make_pair( intens, 1.0 ) );
    }
          
    P++;
  }
 
  // filter out bins of only 1 counts:
  map<double, double>::iterator F = IntensityHist.begin( ) ;
  while( F != IntensityHist.end() ){
    
    if( F->second == MIN_BIN_COUNT ){
      IntensityHist.erase( F++ );
    }
    else{
      F++;
    }
  }
  
  
  
}


///////////////////////////////////////////////////////////////
// prints the intensity data and the model to a text file
//void  BackgroundIntensityBin::writeIntensityMap(){
//  
//  // sort the intensities:
//  sort( IntensityMap.begin(), IntensityMap.end() );
//  
//  read_param* def = new read_param();  
//  // print out
//  string p_name = "ANALYSIS_" + def->search_tag("MY PROJECT NAME");
//  if( p_name[ p_name.size() - 1 ] != '/' ){
//    p_name += "/";
//  }
//  delete def;
//  def = NULL;
//  
//  char buffer[255];
//  sprintf( buffer, "%sIntensityMap_%0.0f_%0.0f.txt", p_name.c_str(), mzCoord, trCoord);
//  string file = buffer;
//  ofstream* WRITER = new ofstream();
//  WRITER->open( file.c_str(), ofstream::out); 
//  
//  if(WRITER->good()){
//    
//    ///////////////////////////////////////////////
//    string SEP = "\t";
//    int count = 0; 
//    vector<double>::iterator P = IntensityMap.begin();
//    while( P != IntensityMap.end() ){
//      // print hist data:
//      sprintf(buffer,"%d%s%0.3f\n", count, SEP.c_str(),(*P));  
//      WRITER->write(buffer,strlen(buffer));
//      count++;
//      P++;
//    }
//    
//    
//    printf("\t\t- Intensity distributions were saved to '%s'\n", file.c_str());
//  }
//  else{
//    printf("\nERROR: opening file '%s'\n", file.c_str());
//  }
//  
//  delete WRITER;
//  WRITER = NULL;
//  
//  
//}



//////////////////////////////////////////////////
// class desctructor of BackgroundIntensityBin
BackgroundIntensityBin::~BackgroundIntensityBin(){
  IntensityMap.clear();
}

//////////////////////////////////////////////////
// class copy constructor of BackgroundIntensityBin
BackgroundIntensityBin::BackgroundIntensityBin(const BackgroundIntensityBin& tmp){
  
  mzCoord = tmp.mzCoord;
  trCoord = tmp.trCoord;
  zCoord = tmp.zCoord;  
  IntensityMap = tmp.IntensityMap;
  
}

//////////////////////////////////////////////////
// class copy constructor of BackgroundIntensityBin
BackgroundIntensityBin::BackgroundIntensityBin(const BackgroundIntensityBin* tmp){
  mzCoord = tmp->mzCoord;
  trCoord = tmp->trCoord;
  zCoord = tmp->zCoord;  
  IntensityMap = tmp->IntensityMap;
  
}


//////////////////////////////////////////////////
// copy constructor:
BackgroundIntensityBin& BackgroundIntensityBin::operator=(const BackgroundIntensityBin& tmp){
  mzCoord = tmp.mzCoord;
  trCoord = tmp.trCoord;
  zCoord = tmp.zCoord;
  IntensityMap = tmp.IntensityMap;
  return *this;
}
