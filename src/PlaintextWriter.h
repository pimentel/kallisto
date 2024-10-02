#ifndef KALLISTO_PLAINTEXT_WRITER_H
#define KALLISTO_PLAINTEXT_WRITER_H

#include <assert.h>
#include <cstdlib>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "common.h"
#include "KmerIndex.h"
#include "GeneModel.h"


void plaintext_writer(
    const std::string& out_name,
    const std::vector<std::string>& targ_ids,
    const std::vector<double>& alpha,
    const std::vector<double>& eff_lens,
    const std::vector<uint32_t>& lens
    );

void plaintext_writer_gene(
    const std::string& out_name,
    const std::vector<std::string>& targ_ids,
    const std::vector<double>& alpha,
    const std::vector<double>& eff_lens,
    const Transcriptome& model
);

std::string to_json(const std::string& id, const std::string& val, bool quote,
    bool comma = true, int level = 1);

void plaintext_aux(
    const std::string& out_fname,
    const std::string& n_targs,
    const std::string& n_bootstrap,
    const std::string& n_processed,
    const std::string& n_pseudoaligned,
    const std::string& n_unique,
    const std::string& version,
    const std::string& index_v,
    const std::string& index_k,
    const std::string& start_time,
    const std::string& call,
    const std::string& cardinality_clashes="");

void writeBatchMatrix(
  const std::string &prefix,
  const KmerIndex &index,
  const std::vector<std::string> &ids,
  std::vector<std::vector<std::pair<int32_t, int32_t>>> &counts);

void writeCellIds(
  const std::string &filename,
  const std::vector<std::string> &ids);


void writeECList(
  const std::string &filename,
  const KmerIndex &index);

void writeFLD(const std::string &filename, const std::vector<std::pair<double, double>> &flds);

void writeGeneList(const std::string &filename, const Transcriptome& model, bool writeNamesOnly = false);

std::vector<double> counts_to_tpm(const std::vector<double>& est_counts,
                                  const std::vector<double>& eff_lens);

template<typename T>
void writeSparseBatchMatrix(
  const std::string &filename,
  std::vector<std::vector<std::pair<int,T>>> &data, int cols) {

  std::ofstream of;
  uint64_t num_rows = data.size();
  uint64_t num_cols = cols;
  uint64_t num_entries = 0;
  for (size_t j = 0; j < data.size(); j++) {
    const auto &v = data[j];
    for (size_t i = 0; i < v.size(); i++) {
      if (v[i].second != T(0)) {
        num_entries++;
      }
    }
  }
  of.open(filename.c_str(), std::ios::out);
  of << "%%MatrixMarket matrix coordinate real general\n";
  of << num_rows << "\t" << num_cols << "\t" << num_entries << "\n";
  if (!data.empty()) {
    for (size_t j = 0; j < data.size(); j++) {
      const auto &v = data[j];
      for (size_t i = 0; i < v.size(); i++) {
        if (v[i].second != T(0)) {
          // everything is 1-based
          of << (j+1) << "\t" << (v[i].first+1) << "\t" << v[i].second << "\n";
        }
      }
    }
  }
  of.close();
}


#endif
