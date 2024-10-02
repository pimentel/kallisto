#ifndef KALLISTO_MINCOLLECTOR_H
#define KALLISTO_MINCOLLECTOR_H

#include "common.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <unordered_map>

#include "KmerIndex.h"
#include "weights.h"
#include "Node.hpp"

const int MAX_FRAG_LEN = 1000;

struct MinCollector {

  MinCollector(KmerIndex& ind, const ProgramOptions& opt)
    :
      index(ind),
      counts(index.ecmapinv.size(), 0),
      flens(MAX_FRAG_LEN),
      unmapped_list(3000000),
      flens_lr(index.target_lens_.size()),
      flens_lr_c(index.target_lens_.size()),
      bias3(4096),
      bias5(4096),
      min_range(opt.min_range),
      k(opt.k),
      mean_fl(0.0),
      has_mean_fl(false),
      mean_fl_trunc(MAX_FRAG_LEN, 0.0),
      has_mean_fl_trunc(false),
      cardinality_clashes(0)
      // eff_len_cache(MAX_FRAG_LEN, 0.0)
       {
         if (opt.fld != 0.0) {
           mean_fl = opt.fld;
           has_mean_fl = true;
         }
       }

  int collect(std::vector<std::pair<const_UnitigMap<Node>, int>>& v1,
              std::vector<std::pair<const_UnitigMap<Node>, int>>& v2,
              bool nonpaired=false);

  int collect(std::vector<std::pair<const_UnitigMap<Node>, int>>& v1) {
    std::vector<std::pair<const_UnitigMap<Node>, int>> dummy;
    return collect(v1,dummy,true);

  }
  int increaseCount(const Roaring& u);
  int decreaseCount(const int ec);

  Roaring intersectECs(std::vector<std::pair<const_UnitigMap<Node>, int32_t>>& v) const;
  Roaring intersectECs_long(std::vector<std::pair<const_UnitigMap<Node>, int32_t>>& v) const;
  Roaring unionECs(std::vector<std::pair<const_UnitigMap<Node>, int32_t>>& v) const;
  Roaring modeECs(std::vector<std::pair<const_UnitigMap<Node>, int32_t>>& v) const;
  int intersectKmersCFC(std::vector<std::pair<const_UnitigMap<Node>, int32_t>>& v1,
                          std::vector<std::pair<const_UnitigMap<Node>, int32_t>>& v3, 
                          std::vector<std::pair<const_UnitigMap<Node>, int32_t>>& v4, 
                          std::vector<std::pair<const_UnitigMap<Node>, int32_t>>& v5,
                          std::vector<std::pair<const_UnitigMap<Node>, int32_t>>& v6,
                          std::vector<std::pair<const_UnitigMap<Node>, int32_t>>& v7, Roaring& r) const;
  int intersectKmers(std::vector<std::pair<const_UnitigMap<Node>, int32_t>>& v1,
                    std::vector<std::pair<const_UnitigMap<Node>, int32_t>>& v2, bool nonpaired, Roaring& r) const;
  int modeKmers(std::vector<std::pair<const_UnitigMap<Node>, int32_t>>& v1,
                    std::vector<std::pair<const_UnitigMap<Node>, int32_t>>& v2, bool nonpaired, Roaring& r) const;
  int findEC(const std::vector<int32_t>& u) const;


  // deprecated
  void write(std::ostream& o) {
    for (int id = 0; id < counts.size(); id++) {
      o << id << "\t" << counts[id] << "\n";
    }
  }

  void loadCounts(ProgramOptions& opt);


  bool countBias(const char *s1, const char *s2, const std::vector<std::pair<const_UnitigMap<Node>,int>> v1, const std::vector<std::pair<const_UnitigMap<Node>,int>> v2, bool paired);
  bool countBias(const char *s1, const char *s2, const std::vector<std::pair<const_UnitigMap<Node>,int>> v1, const std::vector<std::pair<const_UnitigMap<Node>,int>> v2, bool paired, std::vector<int>& biasOut) const;

  // DEPRECATED
  double get_mean_frag_len(bool lenient = false) const;
  double get_sd_frag_len() const;

  // compute the conditional mean of each target given the FLD
  void compute_mean_frag_lens_trunc(bool verbose = true);

  // this function should only be used for SE data
  void init_mean_fl_trunc(double mean, double sd);

  KmerIndex& index;
  std::vector<uint32_t> counts;
  std::vector<uint32_t> flens;
  std::vector<double> unmapped_list;
  std::vector<uint32_t> flens_lr;
  std::vector<uint32_t> flens_lr_c;
  std::vector<int32_t> bias3, bias5;
  int min_range;
  int k;

  double mean_fl;
  bool has_mean_fl;

  std::vector<double> mean_fl_trunc;
  bool has_mean_fl_trunc;

  mutable int cardinality_clashes;
};

std::vector<int> intersect(const std::vector<int>& x, const std::vector<int>& y);

int hexamerToInt(const char *s, bool revcomp);

#endif // KALLISTO_MINCOLLECTOR_H
