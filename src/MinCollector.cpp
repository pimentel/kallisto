#include "MinCollector.h"
#include <algorithm>
#include <limits>

// utility functions

std::vector<int> intersect(const std::vector<int>& x, const std::vector<int>& y) {
  std::vector<int> v;
  auto a = x.begin();
  auto b = y.begin();
  while (a != x.end() && b != y.end()) {
    if (*a < *b) {
      ++a;
    } else if (*b < *a) {
      ++b;
    } else {
      v.push_back(*a);
      ++a;
      ++b;
    }
  }
  return v;
}

void MinCollector::init_mean_fl_trunc(double mean, double sd) {
  auto tmp_trunc_fl = trunc_gaussian_fld(0, MAX_FRAG_LEN, mean, sd);
  assert( tmp_trunc_fl.size() == mean_fl_trunc.size() );

  std::copy( tmp_trunc_fl.begin(), tmp_trunc_fl.end(), mean_fl_trunc.begin() );

  mean_fl = mean_fl_trunc[ MAX_FRAG_LEN - 1 ];

  has_mean_fl = true;
  has_mean_fl_trunc = true;
}

void includeDList(Roaring& u1, Roaring& u2, const Roaring& onlist_sequences) {
  if ((u1 & onlist_sequences).cardinality() != u1.cardinality() || (u2 & onlist_sequences).cardinality() != u2.cardinality()) {
    u1.add(onlist_sequences.cardinality());
    u2.add(onlist_sequences.cardinality());
  }
}

int MinCollector::intersectKmersCFC(std::vector<std::pair<const_UnitigMap<Node>, int32_t>>& v1,
                          std::vector<std::pair<const_UnitigMap<Node>, int32_t>>& v3, 
                          std::vector<std::pair<const_UnitigMap<Node>, int32_t>>& v4, 
                          std::vector<std::pair<const_UnitigMap<Node>, int32_t>>& v5,
                          std::vector<std::pair<const_UnitigMap<Node>, int32_t>>& v6,
                          std::vector<std::pair<const_UnitigMap<Node>, int32_t>>& v7, Roaring& r) const {

  // If a single kmer matches perfectly in the host genome in any frame, throw out the entire read
  Roaring u1 = intersectECs(v1);
  if ((u1 & index.onlist_sequences).cardinality() != u1.cardinality()) return -1;
  Roaring u3 = intersectECs(v3);
  if ((u3 & index.onlist_sequences).cardinality() != u3.cardinality()) return -1;
  Roaring u4 = intersectECs(v4);
  if ((u4 & index.onlist_sequences).cardinality() != u4.cardinality()) return -1;
  Roaring u5 = intersectECs(v5);
  if ((u5 & index.onlist_sequences).cardinality() != u5.cardinality()) return -1;
  Roaring u6 = intersectECs(v6);
  if ((u6 & index.onlist_sequences).cardinality() != u6.cardinality()) return -1;
  Roaring u7 = intersectECs(v7);
  if ((u7 & index.onlist_sequences).cardinality() != u7.cardinality()) return -1;

  // Only take into account ref seqs NOT in the dlist
  u1 &= index.onlist_sequences;
  u3 &= index.onlist_sequences;
  u4 &= index.onlist_sequences;
  u5 &= index.onlist_sequences;
  u6 &= index.onlist_sequences;
  u7 &= index.onlist_sequences;

  if (u1.isEmpty() && u3.isEmpty() && u4.isEmpty() && u5.isEmpty() && u6.isEmpty() && u7.isEmpty()) {
    return -1;
  }

  std::vector<Roaring> u_vector{u1, u3, u4, u5, u6, u7};

  // non-strict intersection
  // to-do: currently the different frames are treated as if they were paired reads
  // this might not be the best way to handle them
  //bool found_non_empty = false;
  //for (const auto& u_ : u_vector) {
  //    if (!found_non_empty) {
  //      r = u_;
  //      if (!r.isEmpty()) {
  //          found_non_empty = true;
  //      }
  //    } else if (!u_.isEmpty()) {
  //        r &= u_;
  //    }
  //}

  // find best match (smallest non-zero roaring)
  // to-do/note: if two reads have the same roaring, we will use the first one (even though the second is just as valid)
  uint32_t smallest_t=-1;
  int frame_idx=0;
  int winner_frame_idx=-1;
  for (const auto& u_ : u_vector) {
    if (u_.cardinality() > 0 && u_.cardinality() < smallest_t) {
      r = u_;
      smallest_t = u_.cardinality();
      winner_frame_idx = frame_idx;
    }
    // Throw warning if new frame is as good as the current winning match
    else if (u_.cardinality() > 0 && u_.cardinality() == smallest_t) {
      // std::cerr << "[warning] Frame " << frame_idx << " has same equivalence class cardinality as winning frame, but will be dismissed." << endl;
      cardinality_clashes++;
    }
    frame_idx++;
  }

  if (r.isEmpty()) {
    return -1;
  }

  // std::cerr << "Aligned frame: " << winner_frame_idx << endl;
  return 1;
}

int MinCollector::modeKmers(std::vector<std::pair<const_UnitigMap<Node>, int32_t>>& v1,
                          std::vector<std::pair<const_UnitigMap<Node>, int32_t>>& v2, bool nonpaired, Roaring& r) const {
  Roaring u1 = intersectECs(v1);
  if (u1.isEmpty()) { u1 = modeECs(v1); }
  
  Roaring u2 = intersectECs(v2);
  if (u2.isEmpty()) { u2 = modeECs(v2); }

  if (u1.isEmpty() && u2.isEmpty()) {
    return -1;
  }

  // non-strict intersection.
  if (u1.isEmpty()) {
    if (v1.empty()) {
      r = u2;
    } else {
      return -1;
    }
  } else if (u2.isEmpty()) {
    if (v2.empty()) {
      r = u1;
    } else {
      return -1;
    }
  } else {
    if (index.dfk_onlist) { // In case we want to not intersect D-list targets
      includeDList(u1, u2, index.onlist_sequences);
    }
    r = u1 & u2;
  }

  if (r.isEmpty()) {
    return -1;
  }
  return 1;
}


int MinCollector::intersectKmers(std::vector<std::pair<const_UnitigMap<Node>, int32_t>>& v1,
                          std::vector<std::pair<const_UnitigMap<Node>, int32_t>>& v2, bool nonpaired, Roaring& r) const {
  Roaring u1 = intersectECs(v1);
  Roaring u2 = intersectECs(v2);

  if (u1.isEmpty() && u2.isEmpty()) {
    return -1;
  }

  // non-strict intersection.
  if (u1.isEmpty()) {
    if (v1.empty()) {
      r = u2;
    } else {
      return -1;
    }
  } else if (u2.isEmpty()) {
    if (v2.empty()) {
      r = u1;
    } else {
      return -1;
    }
  } else {
    if (index.dfk_onlist) { // In case we want to not intersect D-list targets
      includeDList(u1, u2, index.onlist_sequences);
    }
    r = u1 & u2;
  }

  if (r.isEmpty()) {
    return -1;
  }
  return 1;
}

int MinCollector::collect(std::vector<std::pair<const_UnitigMap<Node>, int>>& v1,
                          std::vector<std::pair<const_UnitigMap<Node>, int>>& v2, bool nonpaired) {
  Roaring u;
  int r = intersectKmers(v1, v2, nonpaired, u);
  if (r != -1) {
    return increaseCount(u);
  } else {
    return -1;
  }
}

int MinCollector::findEC(const std::vector<int32_t>& u) const {
  if (u.empty()) {
    return -1;
  }
  if (u.size() == 1) {
    return u[0];
  }

  Roaring r;
  for (int32_t i : u) {
    r.add(i);
  }
  auto search = index.ecmapinv.find(r);
  if (search != index.ecmapinv.end()) {
    return search ->second;
  } else {
    return -1;
  }
}

int MinCollector::increaseCount(const Roaring& r) {

  int ret = 0;
  if (r.isEmpty()) {
    return -1;
  } else {
    auto elem = index.ecmapinv.find(r);
    if (elem != index.ecmapinv.end()) {
      ++counts[elem->second];
      ret = elem->second;
    } else {
      size_t n_elems = index.ecmapinv.size();
      index.ecmapinv.insert({r, n_elems});
      counts.push_back(1);
      ret = n_elems;
    }
  }
  return ret;
}

int MinCollector::decreaseCount(const int ec) {
  assert(ec >= 0);
  --counts[ec];
  return ec;
}

struct ComparePairsBySecond {
  bool operator()(std::pair<KmerEntry,int> a, std::pair<KmerEntry,int> b) {
    return a.second < b.second;
  }
};

Roaring MinCollector::modeECs(std::vector<std::pair<const_UnitigMap<Node>, int32_t>>& v) const {
  Roaring mode;
  if (v.empty()) {
    return mode;
  }
/***
  sort(v.begin(), v.end(), [&](const std::pair<const_UnitigMap<Node>, int>& a, const std::pair<const_UnitigMap<Node>, int>& b)
       {
         if (a.first.isSameReferenceUnitig(b.first) &&
             a.first.getData()->ec[a.first.dist] == b.first.getData()->ec[b.first.dist]) {
           return a.second < b.second;
         } else {
           return a.first.getData()->id < b.first.getData()->id;
         }
       }); // sort by contig, and then first position
***/
  
  mode = v[0].first.getData()->ec[v[0].first.dist].getIndices();
  bool found_nonempty = !mode.isEmpty();
  bool modeMultiMapping = false; 
  Roaring lastEC = mode;
  Roaring ec;
  int modeCount = 0, curCount = 0; 
  for (int i = 1; i < v.size(); i++) {
    // Find a non-empty EC before we start taking the intersection
    if (!found_nonempty) {
      mode = v[i].first.getData()->ec[v[i].first.dist].getIndices();
      found_nonempty = !mode.isEmpty();
      if (found_nonempty && mode.cardinality() == 1 ) {
        modeMultiMapping = true; 
      }
    }
    if (!v[i].first.isSameReferenceUnitig(v[i-1].first) ||
        !(v[i].first.getData()->ec[v[i].first.dist] == v[i-1].first.getData()->ec[v[i-1].first.dist])) {
      ec = v[i].first.getData()->ec[v[i].first.dist].getIndices();
      if (ec == lastEC && !ec.isEmpty()) {
        curCount += 1; 
      } 
      
      // Don't intersect empty EC (because of thresholding)
      if (!(ec == lastEC) && !ec.isEmpty()) {
         if (index.dfk_onlist) { // In case we want to not intersect D-list targets
           includeDList(mode, ec, index.onlist_sequences);
         }
         if (curCount > modeCount && (ec.cardinality() == 1 || modeMultiMapping)) {
           if (ec.cardinality() == 1) {
             modeMultiMapping = false; 
           }
           mode = std::move(lastEC); 
           modeCount = curCount; 
           //curCount = 0; //Technically, not correct mode, but for some reason was giving better results than mode...
         }
         curCount = 0; 
         lastEC = std::move(ec);
      }
    }
  }
  // find the range of support
  int minpos = std::numeric_limits<int>::max();
  int maxpos = 0;
  for (auto& x : v) {
    minpos = std::min(minpos, x.second);
    maxpos = std::max(maxpos, x.second);
  }
  if ((maxpos-minpos + k) < min_range) {
    return {};
  }
  if (modeCount > 0) {
    return mode;
  } else {
    return {};
  }
}

Roaring MinCollector::intersectECs_long(std::vector<std::pair<const_UnitigMap<Node>, int32_t>>& v) const {
  Roaring r;
  if (v.empty()) {
    return r;
  }
  sort(v.begin(), v.end(), [&](const std::pair<const_UnitigMap<Node>, int>& a, const std::pair<const_UnitigMap<Node>, int>& b)
       {
         if (a.first.isSameReferenceUnitig(b.first) &&
             a.first.getData()->ec[a.first.dist] == b.first.getData()->ec[b.first.dist]) {
           return a.second < b.second;
         } else {
           return a.first.getData()->id < b.first.getData()->id;
         }
       }); // sort by contig, and then first position

  
  r = v[0].first.getData()->ec[v[0].first.dist].getIndices();
  bool found_nonempty = !r.isEmpty();
  Roaring lastEC = r;
  Roaring ec;
  int curCount = 0; 

  for (int i = 1; i < v.size(); i++) {

    // Find a non-empty EC before we start taking the intersection
    if (!found_nonempty) {
      r = v[i].first.getData()->ec[v[i].first.dist].getIndices();
      found_nonempty = !r.isEmpty();
    }

    if (!v[i].first.isSameReferenceUnitig(v[i-1].first) ||
        !(v[i].first.getData()->ec[v[i].first.dist] == v[i-1].first.getData()->ec[v[i-1].first.dist])) {
      curCount++; 
      ec = v[i].first.getData()->ec[v[i].first.dist].getIndices();

      // Don't intersect empty EC (because of thresholding)
      if (!(ec == lastEC) && !ec.isEmpty()) {
        if (index.dfk_onlist) { // In case we want to not intersect D-list targets
          includeDList(r, ec, index.onlist_sequences);
        }
        if (curCount > 1) {
          r &= ec;
        }
        if (r.isEmpty()) {
          return r;
        }
        lastEC = std::move(ec);
        curCount = 0; 
      }
    }
  }

  // find the range of support
  int minpos = std::numeric_limits<int>::max();
  int maxpos = 0;

  for (auto& x : v) {
    minpos = std::min(minpos, x.second);
    maxpos = std::max(maxpos, x.second);
  }

  if ((maxpos-minpos + k) < min_range) {
    return {};
  }

  return r;
}

Roaring MinCollector::intersectECs(std::vector<std::pair<const_UnitigMap<Node>, int32_t>>& v) const {
  Roaring r;
  if (v.empty()) {
    return r;
  }
  sort(v.begin(), v.end(), [&](const std::pair<const_UnitigMap<Node>, int>& a, const std::pair<const_UnitigMap<Node>, int>& b)
       {
         if (a.first.isSameReferenceUnitig(b.first) &&
             a.first.getData()->ec[a.first.dist] == b.first.getData()->ec[b.first.dist]) {
           return a.second < b.second;
         } else {
           return a.first.getData()->id < b.first.getData()->id;
         }
       }); // sort by contig, and then first position

  
  r = v[0].first.getData()->ec[v[0].first.dist].getIndices();
  bool found_nonempty = !r.isEmpty();
  Roaring lastEC = r;
  Roaring ec;

  for (int i = 1; i < v.size(); i++) {

    // Find a non-empty EC before we start taking the intersection
    if (!found_nonempty) {
      r = v[i].first.getData()->ec[v[i].first.dist].getIndices();
      found_nonempty = !r.isEmpty();
    }

    if (!v[i].first.isSameReferenceUnitig(v[i-1].first) ||
        !(v[i].first.getData()->ec[v[i].first.dist] == v[i-1].first.getData()->ec[v[i-1].first.dist])) {

      ec = v[i].first.getData()->ec[v[i].first.dist].getIndices();

      // Don't intersect empty EC (because of thresholding)
      if (!(ec == lastEC) && !ec.isEmpty()) {
        if (index.dfk_onlist) { // In case we want to not intersect D-list targets
          includeDList(r, ec, index.onlist_sequences);
        }
        r &= ec;
        if (r.isEmpty()) {
          return r;
        }
        lastEC = std::move(ec);
      }
    }
  }

  // find the range of support
  int minpos = std::numeric_limits<int>::max();
  int maxpos = 0;

  for (auto& x : v) {
    minpos = std::min(minpos, x.second);
    maxpos = std::max(maxpos, x.second);
  }

  if ((maxpos-minpos + k) < min_range) {
    return {};
  }

  return r;
}


void MinCollector::loadCounts(ProgramOptions& opt) {
  int num_ecs = counts.size();
  counts.clear();
  std::ifstream in((opt.output + "/counts.txt"));
  int i = 0;
  if (in.is_open()) {
    std::string line;
    while (getline(in, line)) {
      std::stringstream ss(line);
      int j,c;
      ss >> j;
      ss >> c;
      if (j != i) {
        std::cerr << "Error: equivalence class does not match index. Found "
                  << j << ", expected " << i << std::endl;
        exit(1);
      }
      counts.push_back(c);
      i++;
    }

    if (i != num_ecs) {
      std::cerr << "Error: number of equivalence classes does not match index. Found "
                << i << ", expected " << num_ecs << std::endl;
      exit(1);
    }
  } else {
    std::cerr << "Error: could not open file " << opt.output << "/counts.txt" << std::endl;
    exit(1);

  }

}

double MinCollector::get_mean_frag_len(bool lenient) const {
  if (has_mean_fl) {
    return mean_fl;
  }

  auto total_counts = 0;
  double total_mass = 0.0;

  for ( size_t i = 0 ; i < flens.size(); ++i ) {
    total_counts += flens[i];
    total_mass += static_cast<double>(flens[i] * i);
  }

  if (total_counts == 0) {
    if (!lenient) {
      std::cerr << "Error: could not determine mean fragment length from paired end reads, no pairs mapped to a unique transcript." << std::endl
              << "       Run kallisto quant again with a pre-specified fragment length (option -l)." << std::endl;
      exit(1);
    } else {
      return std::numeric_limits<double>::max();
    }

  }

  // cache the value
  const_cast<double&>(mean_fl) = total_mass / static_cast<double>(total_counts);
  const_cast<bool&>(has_mean_fl) = true;
  return mean_fl;
}

double MinCollector::get_sd_frag_len() const {
  double tmp = get_mean_frag_len(true);
  double m = mean_fl;

  size_t total_counts = 0;
  double total_mass = 0.0;

  for (size_t i = 0; i < flens.size(); ++i) {
    total_counts += flens[i];
    total_mass += flens[i]*(i-m)*(i-m);
  }

  double sd_fl = std::sqrt(total_mass/total_counts);
  return sd_fl;
}

void MinCollector::compute_mean_frag_lens_trunc(bool verbose)  {

  std::vector<int> counts(MAX_FRAG_LEN, 0);
  std::vector<double> mass(MAX_FRAG_LEN, 0.0);

  counts[0] = flens[0];

  for (size_t i = 1; i < MAX_FRAG_LEN; ++i) {
    // mass and counts keep track of the mass/counts up to and including index i
    mass[i] = static_cast<double>( flens[i] * i) + mass[i-1];
    counts[i] = flens[i] + counts[i-1];
    if (counts[i] > 0) {
      mean_fl_trunc[i] = mass[i] / static_cast<double>(counts[i]);
    }
  }

  has_mean_fl_trunc = true;

  if (verbose) {
    std::cerr << "[quant] estimated average fragment length: " <<
      mean_fl_trunc[MAX_FRAG_LEN - 1] << std::endl;
  }
}

int hexamerToInt(const char *s, bool revcomp) {
  int hex = 0;
  if (!revcomp) {
    for (int i = 0; i < 6; i++) {
      hex <<= 2;
      switch (*(s+i) & 0xDF) {
      case 'A': break;
      case 'C': hex += 1; break;
      case 'G': hex += 2; break;
      case 'T': hex += 3; break;
      default: return -1;
      }
    }
  } else {
    for (int i = 0; i < 6; i++) {
      switch (*(s+i) & 0xDF) {
      case 'A': hex += 3 << (2*i);break;
      case 'C': hex += 2 << (2*i); break;
      case 'G': hex += 1 << (2*i); break;
      case 'T': break;
      default: return -1;
      }
    }
  }
  return hex;
}

bool MinCollector::countBias(const char *s1, const char *s2, const std::vector<std::pair<const_UnitigMap<Node>,int>> v1, const std::vector<std::pair<const_UnitigMap<Node>,int>> v2, bool paired) {
  return countBias(s1,s2,v1,v2,paired,bias5);
}

bool MinCollector::countBias(const char *s1, const char *s2, const std::vector<std::pair<const_UnitigMap<Node>,int>> v1, const std::vector<std::pair<const_UnitigMap<Node>,int>> v2, bool paired, std::vector<int>& biasOut) const {

  const int pre = 2, post = 4;

  if (v1.empty() || (paired && v2.empty())) {
    return false;
  }

  auto getPreSeq = [&](const char *s, const_UnitigMap<Node> um, int p) -> int {
    if (s==0) {
      return -1;
    }

    size_t contig_start = 0, contig_length = um.size - um.getGraph()->getK() + 1;
    const Node* n = um.getData();
    auto mc_bounds = n->get_mc_contig(um.dist);
    contig_start += mc_bounds.first;
    contig_length = mc_bounds.second - contig_start;

    size_t pos = um.dist - contig_start;
    if (( um.strand && ((int64_t)(pos - p) >= (int64_t)pre)) ||
        (!um.strand && ((int64_t)(contig_length - 1 - pos - p) >= (int64_t)pre))) {

      int hex = -1;
      //std::cout << "  " << s << "\n";
      if (um.strand) {
        hex = hexamerToInt(um.referenceUnitigToString().c_str() + (contig_start + pos - p - pre), true);
        //std::cout << c.seq.substr(val.getPos()- p - pre,6) << "\n";
      } else {
        int pos_ = (pos + p) + k - post;
        hex = hexamerToInt(um.referenceUnitigToString().c_str() + (contig_start + pos_), false);
        //std::cout << revcomp(c.seq.substr(pos,6)) << "\n";
      }
      return hex;
    }

    return -1;
  };

  // find first contig of read
  const_UnitigMap<Node> um1 = v1[0].first;
  int p1 = v1[0].second;
  for (auto &x : v1) {
    if (x.second < p1) {
      um1 = x.first;
      p1 = x.second;
    }
  }

  bool csense1 = um1.strand; // is this in the direction of the contig?

  int hex5 = getPreSeq(s1, um1, p1);

  /*
  int hex3 = -1;
  if (paired) {
    // do the same for the second read
    KmerEntry val2 = v2[0].first;
    int p2 = v2[0].second;
    for (auto &x : v2) {
      if (x.second < p2) {
        val2 = x.first;
        p2 = x.second;
      }
    }

    Kmer km2 = Kmer((s2+p2));
    bool fw2 = (km2==km2.rep());
    bool csense2 = (fw2 == val2.isFw()); // is this in the direction of the contig?

    hex3 = getPreSeq(s2, km2, fw2, csense2, val2, p2);
  }
  */

  if (hex5 >=0) { // && (!paired || hex3 >= 0)) {
    biasOut[hex5]++;
    //bias3[hex3]++;
  } else {
    return false;
  }

  return false;
}
