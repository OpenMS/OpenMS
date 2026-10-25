// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// $Maintainer: Justin Sing $
// Experimental storage and kernels only, not an OpenDIA pipeline or optimality claim.
#include "OpenDIARecords.h"
#include <OpenMS/ANALYSIS/OPENSWATH/ChromatogramExtractorAlgorithm.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessOpenMS.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/IONMOBILITY/IMTypes.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <algorithm>
#include <bit>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <malloc.h>
#include <numeric>
#include <stdexcept>
#include <sys/resource.h>
#include <unistd.h>

using namespace OpenDIA;
using Clock = std::chrono::steady_clock;
using Index = std::size_t;

void check(bool ok, const std::string& message)
{
  if (!ok) { throw std::runtime_error(message); }
}
bool bits(double a, double b) { return std::bit_cast<uint64_t>(a) == std::bit_cast<uint64_t>(b); }
bool sameMobility(const std::optional<Mobility>& a, const std::optional<Mobility>& b)
{
  return a.has_value() == b.has_value() && (!a || (bits(a->value, b->value) && a->unit == b->unit));
}
bool sameOptional(const std::optional<double>& a, const std::optional<double>& b)
{
  return a.has_value() == b.has_value() && (!a || bits(*a, *b));
}
void equal(const RawPoint& a, const RawPoint& b)
{
  const auto& x = a.coordinates;
  const auto& y = b.coordinates;
  check(bits(a.intensity,b.intensity) && bits(x.mz,y.mz) && bits(x.rt_seconds,y.rt_seconds) &&
        x.frame_id == y.frame_id && x.scan_id == y.scan_id && x.ms_level == y.ms_level &&
        bits(x.isolation_lower_mz,y.isolation_lower_mz) && bits(x.isolation_upper_mz,y.isolation_upper_mz) &&
        sameMobility(x.mobility,y.mobility) && sameOptional(x.faims_voltage,y.faims_voltage) &&
        x.native_id == y.native_id, "Bit-exact full-field roundtrip failed");
}

// Scan coordinates also preserve empty acquisitions. Numeric IDs are caller IDs,
// repeated unchanged when repeating whole fixtures, not invented physical frames.
struct Source
{
  std::vector<RawPoint> rows;
  std::vector<Coordinates> scans;
  std::vector<Index> offsets{0};
  std::vector<unsigned char> sorted;
};
MobilityUnit mobilityUnit(OpenMS::DriftTimeUnit unit)
{
  check(unit == OpenMS::DriftTimeUnit::MILLISECOND || unit == OpenMS::DriftTimeUnit::VSSC,
        "Unsupported continuous mobility unit");
  return unit == OpenMS::DriftTimeUnit::MILLISECOND ? MobilityUnit::MILLISECOND : MobilityUnit::INVERSE_REDUCED;
}
Source readSource(const OpenMS::MSExperiment& experiment)
{
  Source s;
  Index count = 0;
  for (const auto& spectrum : experiment) { count += spectrum.size(); }
  s.rows.reserve(count);
  s.scans.reserve(experiment.size());
  s.offsets.reserve(experiment.size()+1);
  s.sorted.reserve(experiment.size());
  for (Index scan = 0; scan < experiment.size(); ++scan)
  {
    const auto& spectrum = experiment[scan];
    check(spectrum.getPrecursors().size() <= 1, "Multiple isolation windows not representable by WP01 coordinates");
    Coordinates c;
    c.rt_seconds = spectrum.getRT();
    c.frame_id = scan + 1000;
    c.scan_id = scan + 2000;
    c.native_id = spectrum.getNativeID();
    c.ms_level = spectrum.getMSLevel();
    if (!spectrum.getPrecursors().empty())
    {
      const auto& p = spectrum.getPrecursors().front();
      c.isolation_lower_mz = p.getMZ()-p.getIsolationWindowLowerOffset();
      c.isolation_upper_mz = p.getMZ()+p.getIsolationWindowUpperOffset();
    }
    const auto unit = spectrum.getDriftTimeUnit();
    if (unit == OpenMS::DriftTimeUnit::FAIMS_COMPENSATION_VOLTAGE) { c.faims_voltage = spectrum.getDriftTime(); }
    else if (unit != OpenMS::DriftTimeUnit::NONE) { c.mobility = Mobility{spectrum.getDriftTime(),mobilityUnit(unit)}; }
    s.scans.push_back(c);
    s.sorted.push_back(spectrum.isSorted());
    for (Index p = 0; p < spectrum.size(); ++p)
    {
      c.mz = spectrum[p].getMZ();
      if (spectrum.containsIMData())
      {
        const auto [index, im_unit] = spectrum.getIMData();
        const auto& data = spectrum.getFloatDataArrays().at(index);
        check(data.size() == spectrum.size(), "Malformed IM array");
        c.mobility = Mobility{data[p],mobilityUnit(im_unit)};
      }
      s.rows.push_back({c,spectrum[p].getIntensity()});
    }
    s.offsets.push_back(s.rows.size());
  }
  return s;
}
Source replicate(const Source& s, Index factor)
{
  check(factor && (s.rows.empty() || factor <= std::numeric_limits<Index>::max()/s.rows.size()), "Scaling overflow");
  Source r;
  r.rows.reserve(s.rows.size()*factor);
  r.scans.reserve(s.scans.size()*factor);
  r.offsets.reserve(s.scans.size()*factor+1);
  r.sorted.reserve(s.scans.size()*factor);
  for (Index k = 0; k < factor; ++k)
  {
    const auto base = r.rows.size();
    r.rows.insert(r.rows.end(),s.rows.begin(),s.rows.end());
    r.scans.insert(r.scans.end(),s.scans.begin(),s.scans.end());
    r.sorted.insert(r.sorted.end(),s.sorted.begin(),s.sorted.end());
    for (Index i = 1; i < s.offsets.size(); ++i) { r.offsets.push_back(base+s.offsets[i]); }
  }
  return r;
}

struct HotPoint { double mz, intensity; std::optional<Mobility> mobility; };
struct ScanMetadata
{
  double rt, isolation_lower, isolation_upper;
  RecordId frame, scan;
  unsigned ms_level;
  std::optional<double> faims;
  std::string native;
  // Required even for empty scans, separate from per-peak mobility.
  std::optional<Mobility> acquisition_mobility;
  explicit ScanMetadata(const Coordinates& c) : rt(c.rt_seconds), isolation_lower(c.isolation_lower_mz),
    isolation_upper(c.isolation_upper_mz), frame(c.frame_id), scan(c.scan_id), ms_level(c.ms_level),
    faims(c.faims_voltage), native(c.native_id), acquisition_mobility(c.mobility) {}
  Coordinates coordinates() const
  {
    Coordinates c;
    c.rt_seconds=rt; c.isolation_lower_mz=isolation_lower; c.isolation_upper_mz=isolation_upper;
    c.frame_id=frame; c.scan_id=scan; c.ms_level=ms_level; c.faims_voltage=faims; c.native_id=native;
    c.mobility=acquisition_mobility;
    return c;
  }
};
// Full repeated metadata split into columns. No precision reduction or interning.
struct Columns
{
  std::vector<double> mz, intensity, rt, lower, upper;
  std::vector<RecordId> frame, scan;
  std::vector<unsigned> ms;
  std::vector<std::optional<Mobility>> mobility;
  std::vector<std::optional<double>> faims;
  std::vector<std::string> native;
  explicit Columns(const Source& s)
  {
    const auto n=s.rows.size();
    mz.reserve(n); intensity.reserve(n); rt.reserve(n); lower.reserve(n); upper.reserve(n);
    frame.reserve(n); scan.reserve(n); ms.reserve(n); mobility.reserve(n); faims.reserve(n); native.reserve(n);
    for (const auto& p : s.rows)
    {
      const auto& c=p.coordinates;
      mz.push_back(c.mz); intensity.push_back(p.intensity); rt.push_back(c.rt_seconds);
      lower.push_back(c.isolation_lower_mz); upper.push_back(c.isolation_upper_mz);
      frame.push_back(c.frame_id); scan.push_back(c.scan_id); ms.push_back(c.ms_level);
      mobility.push_back(c.mobility); faims.push_back(c.faims_voltage); native.push_back(c.native_id);
    }
  }
  RawPoint row(Index i) const
  {
    Coordinates c;
    c.mz=mz[i]; c.rt_seconds=rt[i]; c.frame_id=frame[i]; c.scan_id=scan[i]; c.ms_level=ms[i];
    c.isolation_lower_mz=lower[i]; c.isolation_upper_mz=upper[i]; c.mobility=mobility[i];
    c.faims_voltage=faims[i]; c.native_id=native[i];
    return {std::move(c),intensity[i]};
  }
};
// Compile-time variants keep dispatch out of hot loops. All own offsets and scan
// metadata, including empty spectra, because WP01 RawTile alone cannot store these.
template<int Kind> struct Layout
{
  std::vector<Index> offsets;
  std::vector<unsigned char> sorted;
  std::vector<ScanMetadata> metadata;
  std::optional<RawTile> tile;
  std::vector<RawPoint> points;
  std::optional<Columns> columns;
  std::vector<HotPoint> hot;
  std::vector<double> mz_values, intensities;
  std::vector<std::optional<Mobility>> mobilities;
  explicit Layout(const Source& s) : offsets(s.offsets), sorted(s.sorted)
  {
    metadata.reserve(s.scans.size());
    for (const auto& c : s.scans) { metadata.emplace_back(c); }
    const auto n=s.rows.size();
    if constexpr (Kind==0) { tile.emplace(RecordContext{0,71,29,503},RecordLimits{n,n},s.rows); }
    if constexpr (Kind==1) { points.reserve(n); points.insert(points.end(),s.rows.begin(),s.rows.end()); }
    if constexpr (Kind==2) { columns.emplace(s); }
    if constexpr (Kind==3)
    {
      hot.reserve(n);
      for (const auto& p : s.rows) { hot.push_back({p.coordinates.mz,p.intensity,p.coordinates.mobility}); }
    }
    if constexpr (Kind==4)
    {
      const bool has_mobility=std::any_of(s.rows.begin(),s.rows.end(),[](const auto& p){return p.coordinates.mobility.has_value();});
      mz_values.reserve(n); intensities.reserve(n);
      if(has_mobility) {mobilities.reserve(n);}
      for (const auto& p : s.rows)
      {
        mz_values.push_back(p.coordinates.mz); intensities.push_back(p.intensity);
        if(has_mobility) {mobilities.push_back(p.coordinates.mobility);}
      }
    }
  }
  Index size() const { return offsets.back(); }
  double mz(Index i) const
  {
    if constexpr (Kind==0) { return tile->point(i).coordinates.mz; }
    if constexpr (Kind==1) { return points[i].coordinates.mz; }
    if constexpr (Kind==2) { return columns->mz[i]; }
    if constexpr (Kind==3) { return hot[i].mz; }
    if constexpr (Kind==4) { return mz_values[i]; }
  }
  double intensity(Index i) const
  {
    if constexpr (Kind==0) { return tile->point(i).intensity; }
    if constexpr (Kind==1) { return points[i].intensity; }
    if constexpr (Kind==2) { return columns->intensity[i]; }
    if constexpr (Kind==3) { return hot[i].intensity; }
    if constexpr (Kind==4) { return intensities[i]; }
  }
  std::optional<Mobility> mobility(Index i) const
  {
    if constexpr (Kind==0) { return tile->point(i).coordinates.mobility; }
    if constexpr (Kind==1) { return points[i].coordinates.mobility; }
    if constexpr (Kind==2) { return columns->mobility[i]; }
    if constexpr (Kind==3) { return hot[i].mobility; }
    if constexpr (Kind==4) { return mobilities.empty() ? std::nullopt : mobilities[i]; }
  }
  HotPoint numeric(Index i) const
  {
    if constexpr (Kind==0)
    {
      const auto p=tile->point(i);
      return {p.coordinates.mz,p.intensity,p.coordinates.mobility};
    }
    else {return {mz(i),intensity(i),mobility(i)};}
  }
  RawPoint row(Index i) const
  {
    if constexpr (Kind==0) { return tile->point(i); }
    if constexpr (Kind==1) { return points[i]; }
    if constexpr (Kind==2) { return columns->row(i); }
    if constexpr (Kind>=3)
    {
      const auto scan=Index(std::upper_bound(offsets.begin(),offsets.end(),i)-offsets.begin()-1);
      auto c=metadata[scan].coordinates();
      c.mz=mz(i); c.mobility=mobility(i);
      return {std::move(c),intensity(i)};
    }
  }
  double dense() const
  {
    double sum=0;
    if constexpr (Kind==0)
    {
      for (Index i=0;i<size();++i) { const auto p=tile->point(i); sum+=p.coordinates.mz*p.intensity; }
    }
    else if constexpr (Kind==1)
    {
      for (const auto& p : points) { sum+=p.coordinates.mz*p.intensity; }
    }
    else { for (Index i=0;i<size();++i) { sum+=mz(i)*intensity(i); } }
    return sum;
  }
};

struct Target { double mz; std::optional<Mobility> mobility; };
std::vector<Target> targets(const Source& s)
{
  check(!s.rows.empty(), "No peaks in source");
  std::vector<Index> order(s.rows.size());
  std::iota(order.begin(),order.end(),0);
  std::stable_sort(order.begin(),order.end(),[&](Index a,Index b){return s.rows[a].coordinates.mz<s.rows[b].coordinates.mz;});
  std::vector<Target> result;
  result.reserve(34);
  for (Index k=0;k<32;++k)
  {
    const auto& c=s.rows[order[k*(order.size()-1)/31]].coordinates;
    result.push_back({c.mz,c.mobility});
  }
  result.push_back({s.rows[order.front()].coordinates.mz-100,std::nullopt});
  result.push_back({s.rows[order.back()].coordinates.mz+100,std::nullopt});
  std::stable_sort(result.begin(),result.end(),[](const auto& a,const auto& b){return a.mz<b.mz;});
  return result;
}
// Deliberately simple experimental lower-bound/open-edge kernel, not copied
// from the production extractor. Each result is a per-scan sum.
template<class L> double extract(const L& l, Index scan, const Target& t, double width, bool with_im)
{
  // Sortedness is verified before any timed kernel, never repaired here.
  const double lower=t.mz-width/2, upper=t.mz+width/2;
  Index lo=l.offsets[scan], hi=l.offsets[scan+1];
  while(lo<hi)
  {
    const auto mid=lo+(hi-lo)/2;
    if(l.mz(mid)<=lower) {lo=mid+1;} else {hi=mid;}
  }
  double sum=0;
  for(Index i=lo;i<l.offsets[scan+1];++i)
  {
    const auto p=l.numeric(i); // RawTile gets one actual by-value point per visited row.
    if(!(p.mz<upper)) {break;}
    if(with_im && t.mobility)
    {
      const auto& im=p.mobility;
      if(!im || im->unit!=t.mobility->unit || !(im->value>t.mobility->value-0.05 && im->value<t.mobility->value+0.05)) {continue;}
    }
    sum+=p.intensity;
  }
  return sum;
}
template<class L> double extraction(const L& l,const std::vector<Target>& ts,double width,bool im=true)
{
  double sum=0;
  for(Index scan=0;scan<l.metadata.size();++scan)
  { for(const auto& t:ts) {sum+=extract(l,scan,t,width,im);} }
  return sum;
}
// Exercise the existing public contiguous-array kernel, including adapter work.
// AoS packs each scan once. SoA passes its actual numeric vectors directly.
// IM optionals must be packed, with mixed/missing units filtered explicitly.
// Iterators restart per target, so this is not the production pipeline schedule.
template<int K> double publicKernel(const Layout<K>& l,const std::vector<Target>& ts,double width,
                                    std::vector<double>* results=nullptr)
{
  OpenMS::ChromatogramExtractorAlgorithm extractor;
  const bool need_im=std::any_of(ts.begin(),ts.end(),[](const auto& t){return t.mobility.has_value();});
  double checksum=0;
  for(Index scan=0;scan<l.metadata.size();++scan)
  {
    const auto first=l.offsets[scan], last=l.offsets[scan+1], count=last-first;
    std::vector<double> packed_mz, packed_intensity, packed_im;
    std::vector<std::optional<Mobility>> optional_im;
    const std::vector<double>* mz=nullptr;
    const std::vector<double>* intensity=nullptr;
    Index base=0;
    if constexpr(K==2) {mz=&l.columns->mz; intensity=&l.columns->intensity; base=first;}
    else if constexpr(K==4) {mz=&l.mz_values; intensity=&l.intensities; base=first;}
    else {packed_mz.reserve(count); packed_intensity.reserve(count); mz=&packed_mz; intensity=&packed_intensity;}
    if(need_im) {optional_im.reserve(count); packed_im.reserve(count);}
    for(Index i=first;i<last;++i)
    {
      if constexpr(K!=2 && K!=4)
      {
        const auto p=l.numeric(i);
        packed_mz.push_back(p.mz); packed_intensity.push_back(p.intensity);
        if(need_im) {optional_im.push_back(p.mobility);}
      }
      else if(need_im) {optional_im.push_back(l.mobility(i));}
    }
    const bool all_im=need_im && !optional_im.empty() && optional_im.front() &&
      std::all_of(optional_im.begin(),optional_im.end(),[&](const auto& im){return im && im->unit==optional_im.front()->unit;});
    if(all_im) {for(const auto& im:optional_im) {packed_im.push_back(im->value);}}
    std::vector<double> filtered_mz, filtered_intensity, filtered_im;
    for(const auto& t:ts)
    {
      auto begin=mz->cbegin()+base, end=begin+count, mz_it=begin;
      auto int_it=intensity->cbegin()+base;
      double value=0;
      if(!t.mobility) {extractor.extract_value_tophat(begin,mz_it,end,int_it,t.mz,value,width,false);}
      else if(all_im && optional_im.front()->unit==t.mobility->unit)
      {
        auto im_it=packed_im.cbegin();
        extractor.extract_value_tophat(begin,mz_it,end,int_it,im_it,t.mz,t.mobility->value,value,width,0.1,false);
      }
      else
      {
        // The public API has no optional IM or unit column. Filter only actual
        // points with a compatible unit, never fill missing IM with a number.
        filtered_mz.clear(); filtered_intensity.clear(); filtered_im.clear();
        filtered_mz.reserve(count); filtered_intensity.reserve(count); filtered_im.reserve(count);
        for(Index i=0;i<count;++i)
        {
          const auto& im=optional_im[i];
          if(im && im->unit==t.mobility->unit)
          {
            filtered_mz.push_back((*mz)[base+i]); filtered_intensity.push_back((*intensity)[base+i]); filtered_im.push_back(im->value);
          }
        }
        auto filtered_begin=filtered_mz.cbegin(), filtered_end=filtered_mz.cend(), filtered_it=filtered_begin;
        auto filtered_int=filtered_intensity.cbegin(), filtered_mobility=filtered_im.cbegin();
        extractor.extract_value_tophat(filtered_begin,filtered_it,filtered_end,filtered_int,filtered_mobility,
                                      t.mz,t.mobility->value,value,width,0.1,false);
      }
      checksum+=value;
      if(results) {results->push_back(value);}
    }
  }
  return checksum;
}
uint64_t mix(uint64_t h,uint64_t v) {return (h^v)*1099511628211ULL;}
uint64_t hashRow(const RawPoint& p)
{
  const auto& c=p.coordinates;
  uint64_t h=1469598103934665603ULL;
  for(double v:{c.mz,p.intensity,c.rt_seconds,c.isolation_lower_mz,c.isolation_upper_mz}) {h=mix(h,std::bit_cast<uint64_t>(v));}
  for(uint64_t v:{c.frame_id,c.scan_id,uint64_t(c.ms_level),uint64_t(c.mobility.has_value()),uint64_t(c.faims_voltage.has_value())}) {h=mix(h,v);}
  if(c.mobility) {h=mix(h,std::bit_cast<uint64_t>(c.mobility->value)); h=mix(h,uint64_t(c.mobility->unit));}
  if(c.faims_voltage) {h=mix(h,std::bit_cast<uint64_t>(*c.faims_voltage));}
  for(unsigned char v:c.native_id) {h=mix(h,v);}
  return h;
}
std::vector<Index> sparseIds(Index count)
{
  std::vector<Index> ids;
  ids.reserve(std::min(count,Index(4096)));
  // Deliberate sparse order, not a sequential prefix. Repeated ids are allowed.
  for(Index i=0;i<std::min(count,Index(4096));++i) {ids.push_back((i*104729ULL+17)%count);}
  return ids;
}
template<class L> double sparse(const L& l,const std::vector<Index>& ids)
{
  double sum=0;
  for(auto i:ids) {const auto p=l.numeric(i); sum+=p.mz*p.intensity;}
  return sum;
}
template<class L> uint64_t fullrow(const L& l,const std::vector<Index>& ids)
{
  uint64_t hash=0;
  for(auto i:ids) {hash=mix(hash,hashRow(l.row(i)));}
  return hash;
}
bool near(double a,double b) {return std::abs(a-b)<=1e-10*std::max({1.0,std::abs(a),std::abs(b)});}
template<class L> void validate(const L& l,const Source& s,const std::vector<Target>& ts)
{
  for (auto sorted : s.sorted) {check(sorted, "Extraction requires sorted scans, source is never silently resorted");}
  check(l.size()==s.rows.size() && l.offsets==s.offsets && l.sorted==s.sorted,"Acquisition chunks/order changed");
  for(Index scan=0;scan<s.scans.size();++scan) {equal({l.metadata[scan].coordinates(),0},{s.scans[scan],0});}
  for(Index i=0;i<s.rows.size();++i) {equal(l.row(i),s.rows[i]);}
  // A separately constructed const-reference AoS oracle, outside timing.
  Layout<1> oracle(s);
  for(double width:{0.05,10.0})
  {
    for(Index scan=0;scan<s.scans.size();++scan)
    {
      for(const auto& t:ts) {check(near(extract(l,scan,t,width,true),extract(oracle,scan,t,width,true)),"Extraction mismatch");}
    }
  }
  // Check every individual public result, not only the aggregate checksum.
  for(double width:{0.05,10.0})
  {
    std::vector<double> values;
    values.reserve(s.scans.size()*ts.size());
    publicKernel(l,ts,width,&values);
    Index index=0;
    for(Index scan=0;scan<s.scans.size();++scan)
    {
      for(const auto& t:ts) {check(near(values[index++],extract(oracle,scan,t,width,true)),"Public array kernel mismatch");}
    }
  }
  const auto ids=sparseIds(s.rows.size());
  check(bits(sparse(l,ids),sparse(oracle,ids)) && fullrow(l,ids)==fullrow(oracle,ids),"Sparse IDs/order mismatch");
}

// Actual public semantic oracle on the unreplicated fixture. The public API
// skips empty spectra and only supports IM arrays, not scalar scan IM or FAIMS.
// IM checks therefore select original spectra by actual per-peak IM unit, never
// synthesize arrays or alter peaks. mz-only checks include every nonempty scan.
void publicReference(const std::shared_ptr<OpenMS::MSExperiment>& experiment,const Source& s,const std::vector<Target>& ts)
{
  for(int mode=0;mode<3;++mode)
  {
    const bool with_im=mode!=0;
    const auto unit=mode==1 ? MobilityUnit::MILLISECOND : MobilityUnit::INVERSE_REDUCED;
    auto selected=std::make_shared<OpenMS::MSExperiment>();
    std::vector<Index> scan_ids;
    for(Index scan=0;scan<experiment->size();++scan)
    {
      const auto& spectrum=(*experiment)[scan];
      if(spectrum.empty()) {continue;}
      if(with_im && (!spectrum.containsIMData() || mobilityUnit(spectrum.getIMData().second)!=unit)) {continue;}
      selected->addSpectrum(spectrum);
      scan_ids.push_back(scan);
    }
    if(scan_ids.empty()) {continue;}
    auto access=std::make_shared<OpenMS::SpectrumAccessOpenMS>(selected);
    Layout<1> oracle(s);
    for(double width:{0.05,10.0})
    {
      std::vector<OpenMS::ChromatogramExtractorAlgorithm::ExtractionCoordinates> coords;
      std::vector<OpenSwath::ChromatogramPtr> output;
      std::vector<Index> target_ids;
      for(Index i=0;i<ts.size();++i)
      {
        if(with_im && ts[i].mobility && ts[i].mobility->unit!=unit) {continue;}
        OpenMS::ChromatogramExtractorAlgorithm::ExtractionCoordinates c;
        c.mz=ts[i].mz; c.rt_start=0; c.rt_end=-1; c.id=std::to_string(i);
        c.ion_mobility=-1; // Public API sentinel disables IM for targets without an IM coordinate.
        if(with_im && ts[i].mobility) {c.ion_mobility=ts[i].mobility->value;}
        coords.push_back(c); output.push_back(std::make_shared<OpenSwath::Chromatogram>()); target_ids.push_back(i);
      }
      OpenMS::ChromatogramExtractorAlgorithm extractor;
      extractor.extractChromatograms(access,output,coords,width,false,with_im ? 0.1 : -1,"tophat");
      for(Index t=0;t<target_ids.size();++t)
      {
        const auto& times=output[t]->getTimeArray()->data;
        const auto& values=output[t]->getIntensityArray()->data;
        check(values.size()==scan_ids.size() && times.size()==scan_ids.size(),"Public extractor scan selection changed");
        for(Index i=0;i<scan_ids.size();++i)
        {
          const auto scan=scan_ids[i];
          check(bits(times[i],s.scans[scan].rt_seconds),"Public extractor time axis differs");
          check(near(values[i],extract(oracle,scan,ts[target_ids[t]],width,with_im)),"Public extractor intensity differs");
        }
      }
    }
  }
}

Source tiny()
{
  Source s;
  Coordinates c;
  c.frame_id=9; c.scan_id=71; c.native_id=std::string(512,'x'); c.faims_voltage=0.0;
  s.scans={c,c,c}; s.offsets={0,0,1,6}; s.sorted={1,1,1};
  c.mz=2; s.rows.push_back({c,7});
  for(double mz:{1.0,2.0,2.0,3.0,4.0})
  {
    c.mz=mz;
    c.mobility=mz==3.0 ? std::nullopt : std::optional<Mobility>{{mz,MobilityUnit::INVERSE_REDUCED}};
    s.rows.push_back({c,mz});
  }
  return s;
}
template<int K> void tinyCheck(const Source& s)
{
  Layout<K> l(s);
  validate(l,s,targets(s));
  Target t{2,std::nullopt};
  check(extract(l,0,t,2,false)==0 && extract(l,1,t,2,false)==7 && extract(l,2,t,2,false)==4,"Empty/singleton/duplicate/open boundary failed");
  t.mobility=Mobility{2,MobilityUnit::INVERSE_REDUCED};
  check(extract(l,2,t,4,true)==4,"Missing IM or IM boundary failed");
  t.mobility=Mobility{2.05,MobilityUnit::INVERSE_REDUCED};
  // Exact IM open edges constructed from the same arithmetic, avoiding decimal assumptions.
  auto edge=s; edge.rows[2].coordinates.mobility->value=t.mobility->value-0.05;
  edge.rows[3].coordinates.mobility->value=t.mobility->value+0.05;
  Layout<K> e(edge);
  check(extract(e,2,t,4,true)==0,"IM exact open edges failed");
}
void selfTest()
{
  const auto s=tiny();
  tinyCheck<0>(s); tinyCheck<1>(s); tinyCheck<2>(s); tinyCheck<3>(s); tinyCheck<4>(s);
  auto absent=s;
  for(auto& row:absent.rows) {row.coordinates.mobility.reset();}
  for(auto& scan:absent.scans) {scan.mobility.reset();}
  Layout<4> no_im(absent);
  check(no_im.mobilities.empty() && no_im.mobilities.capacity()==0,"Absent IM must not allocate a factored SoA column");
  validate(no_im,absent,targets(absent));
}
uint64_t allocated()
{
  const auto m=mallinfo2();
  return uint64_t(m.uordblks)+uint64_t(m.hblkhd);
}
uint64_t rss()
{
  std::ifstream f("/proc/self/status");
  std::string line;
  while(std::getline(f,line))
  {
    if(line.starts_with("VmRSS:")) {return std::stoull(line.substr(6))*1024;}
  }
  throw std::runtime_error("/proc/self/status VmRSS unavailable");
}
uint64_t peakRss()
{
  rusage r{};
  check(getrusage(RUSAGE_SELF,&r)==0,"getrusage failed");
  return uint64_t(r.ru_maxrss)*1024;
}
struct Timing {std::vector<double> samples; double checksum=0;};
// All timing calls consume their results. No allocation accounting in hot access.
template<class F> Timing measure(Index repeats,F&& f)
{
  Timing t;
  t.samples.reserve(repeats);
  t.checksum=static_cast<double>(f()); // warm-up also consumed
  for(Index i=0;i<repeats;++i)
  {
    asm volatile("" ::: "memory"); // Prevent invariant kernel hoisting across repeats.
    const auto start=Clock::now();
    const auto value=f();
    asm volatile("" ::: "memory");
    const auto stop=Clock::now();
    t.checksum+=static_cast<double>(value);
    t.samples.push_back(std::chrono::duration<double,std::milli>(stop-start).count());
  }
  check(std::isfinite(t.checksum),"Nonfinite benchmark checksum");
  return t;
}
void printTiming(const char* name,const Timing& t,bool comma)
{
  auto sorted=t.samples;
  std::sort(sorted.begin(),sorted.end());
  const auto n=sorted.size();
  const double median=(sorted[(n-1)/2]+sorted[n/2])/2;
  std::cout << (comma ? "," : "") << '"' << name << "\":{\"checksum\":" << t.checksum
            << ",\"mean_ms\":" << std::accumulate(sorted.begin(),sorted.end(),0.0)/n
            << ",\"median_ms\":" << median << ",\"min_ms\":" << sorted.front()
            << ",\"max_ms\":" << sorted.back() << ",\"samples_ms\":[";
  for(Index i=0;i<n;++i) {std::cout << (i ? "," : "") << t.samples[i];}
  std::cout << "]}";
}
template<int K> void run(const Source& source,const Source& original,const std::string& variant,Index factor,Index repeats,const std::vector<Target>& ts)
{
  // Source, experiment, targets and repeated canonical rows remain alive equally
  // for every layout. These measurements are process/allocator totals, not payload.
  const auto baseline=rss();
  const auto before=allocated();
  const auto start=Clock::now();
  Layout<K> layout(source);
  const auto stop=Clock::now();
  const auto after=allocated();
  const auto resident=rss(), peak=peakRss();
  validate(layout,source,ts);
  const auto ids=sparseIds(source.rows.size());
  const auto dense=measure(repeats,[&]{return layout.dense();});
  const auto narrow=measure(repeats,[&]{return extraction(layout,ts,0.05);});
  const auto wide=measure(repeats,[&]{return extraction(layout,ts,10.0);});
  const auto gather=measure(repeats,[&]{return sparse(layout,ids);});
  const auto audit=measure(repeats,[&]{return fullrow(layout,ids);});
  const auto public_narrow=measure(repeats,[&]{return publicKernel(layout,ts,0.05);});
  const auto public_wide=measure(repeats,[&]{return publicKernel(layout,ts,10.0);});
  const auto final_peak=peakRss();
  std::cout << std::setprecision(17) << "LAYOUT_JSON {\"schema_version\":1,\"variant\":\"" << variant
            << "\",\"rows\":" << source.rows.size() << ",\"scans\":" << source.scans.size()
            << ",\"source_rows\":" << original.rows.size() << ",\"source_scans\":" << original.scans.size()
            << ",\"replication_factor\":" << factor << ",\"scaling\":\"repeated-fixture\",\"targets\":" << ts.size()
            << ",\"build_ms\":" << std::chrono::duration<double,std::milli>(stop-start).count()
            << ",\"allocator_before_bytes\":" << before << ",\"allocator_after_bytes\":" << after
            << ",\"allocator_delta_bytes\":" << int64_t(after)-int64_t(before)
            << ",\"baseline_rss_bytes\":" << baseline << ",\"rss_bytes\":" << resident << ",\"peak_rss_bytes\":" << peak
            << ",\"construction_rss_bytes\":" << resident << ",\"construction_peak_rss_bytes\":" << peak
            << ",\"final_peak_rss_bytes\":" << final_peak
            << ",\"memory_snapshot_stage\":\"after-construction-before-layout-validation-and-timed-kernels\""
            << ",\"build_validation\":\"" << (K==0 ? "actual-RawTile-validation-included" : "candidate-packing-only-no-RawTile-validation") << '"'
            << ",\"build_comparison\":\"not-equivalent-validation-work-do-not-rank-as-fair-construction-speed\""
            << ",\"public_kernel_adapter\":\"AoS-per-scan-pack-SoA-direct-numeric-vectors-optional-IM-pack-iterators-reset-per-target\""
            << ",\"sizeof_raw_point\":" << sizeof(RawPoint) << ",\"sizeof_coordinates\":" << sizeof(Coordinates)
            << ",\"sizeof_hot_point\":" << sizeof(HotPoint) << ",\"sizeof_scan_metadata\":" << sizeof(ScanMetadata)
            << ",\"sparse_rows\":" << ids.size() << ",\"fullrow_rows\":" << ids.size()
            << ",\"narrow_width_th\":0.05,\"wide_width_th\":10,\"im_width\":0.1"
            << ",\"memory_scope\":\"glibc-retained-allocation-and-whole-process-rss-at-construction-not-payload\""
            << ",\"public_reference\":\"unreplicated-nonempty-all-MS-levels-mz-and-perpeak-IM-unit-subsets\",\"kernel\":\"experimental-open-tophat-with-IM-where-target-present\",\"timings\":{";
  printTiming("dense",dense,false); printTiming("narrow",narrow,true); printTiming("wide",wide,true);
  printTiming("sparse",gather,true); printTiming("fullrow",audit,true);
  printTiming("public_kernel_narrow",public_narrow,true); printTiming("public_kernel_wide",public_wide,true);
  std::cout << "}}\n";
}
Index parseCount(const std::string& value)
{
  check(!value.empty() && std::all_of(value.begin(),value.end(),[](unsigned char c){return c>='0' && c<='9';}),
        "Counts must be unsigned decimal integers");
  Index consumed=0;
  const auto result=std::stoull(value,&consumed);
  check(consumed==value.size() && result<=std::numeric_limits<Index>::max(),"Count out of range");
  return static_cast<Index>(result);
}
int main(int argc,char** argv)
{
  try
  {
    selfTest(); // Always active, including -DNDEBUG.
    if(argc==2 && std::string(argv[1])=="--self-test") {std::cout << "LAYOUT_SELF_TEST_OK\n"; return 0;}
    check(argc==5,"Usage: layout_benchmark <mzML> <rawtile|aos|soa|factored-aos|factored-soa> <min_points> <repeats>");
    const std::string variant=argv[2];
    check(variant=="rawtile" || variant=="aos" || variant=="soa" || variant=="factored-aos" || variant=="factored-soa","Unknown variant");
    const Index minimum=parseCount(argv[3]), repeats=parseCount(argv[4]);
    check(repeats>0,"repeats must be positive");
    auto experiment=std::make_shared<OpenMS::MSExperiment>();
    OpenMS::MzMLFile().load(argv[1],*experiment);
    const auto original=readSource(*experiment);
    check(!original.rows.empty(),"Fixture has no points");
    for(auto sorted:original.sorted) {check(sorted,"Unsorted source is retained but sorted extraction cannot apply");}
    const auto ts=targets(original);
    publicReference(experiment,original,ts);
    const Index factor=std::max(Index(1),minimum/original.rows.size()+Index(minimum%original.rows.size()!=0));
    const auto source=replicate(original,factor);
    if(variant=="rawtile") {run<0>(source,original,variant,factor,repeats,ts);}
    if(variant=="aos") {run<1>(source,original,variant,factor,repeats,ts);}
    if(variant=="soa") {run<2>(source,original,variant,factor,repeats,ts);}
    if(variant=="factored-aos") {run<3>(source,original,variant,factor,repeats,ts);}
    if(variant=="factored-soa") {run<4>(source,original,variant,factor,repeats,ts);}
    return 0;
  }
  catch(const std::exception& e) {std::cerr << "layout_benchmark: " << e.what() << '\n'; return 1;}
}
