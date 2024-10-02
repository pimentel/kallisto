#include "PlaintextWriter.h"
#include <iomanip>
#include <sstream>

std::vector<double> counts_to_tpm(const std::vector<double>& est_counts,
    const std::vector<double>& eff_lens) {
  assert( est_counts.size() == eff_lens.size() );
  const double MILLION {1e6};

  std::vector<double> tpm(est_counts.size(), 0.0);

  double total_mass = 0.0;

  for (size_t i = 0; i < est_counts.size(); ++i) {
    if (eff_lens[i] < 1.0) {
      std::cerr << "Why is this eff_len < 1.0? id: " << i << std::endl;
    }
    tpm[i] = (est_counts[i] / eff_lens[i]);
    total_mass += tpm[i];
  }

  for (size_t i = 0; i < est_counts.size(); ++i) {
    tpm[i] = (tpm[i] / total_mass) * MILLION;
  }

  return tpm;
}

void plaintext_writer(
    const std::string& out_name,
    const std::vector<std::string>& targ_ids,
    const std::vector<double>& alpha,
    const std::vector<double>& eff_lens,
    const std::vector<uint32_t>& lens
    ){

  std::ofstream of;
  of.open( out_name );

  if (!of.is_open()) {
    std::cerr << "Error: Couldn't open file: " << out_name << std::endl;

    exit(1);
  }

  auto tpm = counts_to_tpm(alpha, eff_lens);

  of << "target_id" << "\t"
    /* << "kallisto_id" << "\t" */
    << "length" << "\t"
    << "eff_length" << "\t"
    << "est_counts" << "\t"
    << "tpm" << std::endl;

  for (auto i = 0; i < alpha.size(); ++i) {
    of << targ_ids[i] << '\t'
      /* << i << '\t' */
      << lens[i] << '\t'
      << eff_lens[i] << '\t'
      << alpha[i] << '\t'
      << tpm[i] << std::endl;
  }

  of.close();
}

void plaintext_writer_gene(
    const std::string& out_name,
    const std::vector<std::string>& targ_ids,
    const std::vector<double>& alpha,
    const std::vector<double>& eff_lens,
    const Transcriptome& model
){

  std::ofstream of;
  of.open( out_name );

  if (!of.is_open()) {
    std::cerr << "Error: Couldn't open file: " << out_name << std::endl;

    exit(1);
  }

  std::vector<double> gc;
  std::vector<double> gc_tpm;
  gc.assign(model.genes.size(), 0.0);
  gc_tpm.assign(model.genes.size(), 0.0);
  auto tpm = counts_to_tpm(alpha, eff_lens);
  for (int i = 0; i < alpha.size(); i++) {
    if (alpha[i] > 0.0) {
      int g_id = model.transcripts[i].gene_id;
      if (g_id != -1) {
        gc[g_id] += alpha[i];
        gc_tpm[g_id] += tpm[i];
      }
    }
  }

  of << "gene_id" << "\t"
     << "gene_name" << "\t"
     << "est_counts" << "\t"
     << "tpm" << "\n";

  for (int i = 0; i < gc.size(); i++) {
    of << model.genes[i].name << '\t'
       << model.genes[i].commonName << '\t'
       << gc[i] << '\t'
       << gc_tpm[i] << "\n";
  }

  of.close();
}

std::string to_json(const std::string& id, const std::string& val, bool quote,
    bool comma, int level) {
  std::string out;

  for (auto i = 0; i < level; ++i) {
    out += "\t";
  }

  out += '"';
  out += id;
  out += "\": ";
  if (quote) {
    out += '"';
  }
  out += val;
  if (quote) {
    out += '"';
  }
  if (comma) {
    out += ',';
  }

  return out;
}


void plaintext_aux(
    const std::string& out_name,
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
    const std::string& cardinality_clashes) {
  std::ofstream of;
  of.open( out_name );

  double p_uniq =0.0;
  double p_aln = 0.0;
  std::stringstream ss;
  std::string p_aln_s, p_uniq_s;

  try {
    double nreads = std::stod(n_processed);
    double naln = std::stod(n_pseudoaligned);
    double nuniq = std::stod(n_unique);
    if (nreads > 0) {
      p_uniq = 100.0* nuniq / nreads;
      p_aln = 100.0* naln / nreads;
    }
  } catch (std::invalid_argument) {}


  ss << std::fixed << std::setprecision(1) << p_uniq;
  p_uniq_s = ss.str();
  ss.str("");
  ss << std::fixed << std::setprecision(1) << p_aln;
  p_aln_s = ss.str();

  of << "{" << std::endl <<
    to_json("n_targets", n_targs, false) << std::endl <<
    to_json("n_bootstraps", n_bootstrap, false) << std::endl <<
    to_json("n_processed", n_processed, false) << std::endl <<
    to_json("n_pseudoaligned", n_pseudoaligned, false) << std::endl <<
    to_json("n_unique", n_unique, false) << std::endl <<
    to_json("p_pseudoaligned", p_aln_s, false) << std::endl <<
    to_json("p_unique", p_uniq_s, false) << std::endl <<
    to_json("kallisto_version", version, true) << std::endl <<
    to_json("index_version", index_v, false) << std::endl <<
    to_json("k-mer length", index_k, false) << std::endl <<
    to_json("start_time", start_time, true) << std::endl;
    if (cardinality_clashes != "") {
    of << to_json("call", call, true) << std::endl;
    of <<  to_json("n_frame_clashes", cardinality_clashes, false, false) << std::endl;
    } else {
    of << to_json("call", call, true, false) << std::endl;
    }
    of << "}" << std::endl;

  of.close();
}


void writeBatchMatrix(
  const std::string &prefix,
  const KmerIndex &index,
  const std::vector<std::string> &ids,
  std::vector<std::vector<std::pair<int32_t, int32_t>>> &counts) {

  std::string ecfilename = prefix + ".ec";
  std::string countsfilename = prefix + ".tsv";
  std::string cellnamesfilename = prefix + ".cells";

  writeECList(ecfilename, index);
  writeCellIds(cellnamesfilename, ids);
  writeSparseBatchMatrix(countsfilename, counts, index.ecmapinv.size());

/*
    countsof.open(countsfilename.c_str(), std::ios::out);
    if (!counts.empty()) {
      for (int j = 0; j < counts.size(); j++) {
        const auto &v = counts[j];
        for (int i = 0; i < v.size(); i++) {
          if (v[i] != 0) {
            countsof << i << "\t" << j << "\t" << v[i] << "\n";
          }
        }
      }
    }
    countsof.close();
    */

}



void writeECList(
  const std::string &filename,
  const KmerIndex &index) {
    std::ofstream ecof;
    ecof.open(filename.c_str(), std::ios::out);

    if (!ecof.is_open()) {
      std::cerr << "Error: Couldn't open file: " << filename << std::endl;
      exit(1);
    }

    // output equivalence classes in the form "EC TXLIST";
    std::vector<Roaring> ecmap(index.ecmapinv.size());
    for (const auto& it : index.ecmapinv) {
      ecmap[it.second] = it.first;
    }
    for (int i = 0 ; i < ecmap.size(); i++) {
      ecof << i << "\t";
      // output the rest of the class
      bool first = true;
      for (auto x : ecmap[i]) {
        if (!first) {
          ecof << ",";
        } else {
          first = false;
        }
        ecof << x;
      }
      ecof << "\n";
    }
    ecof.close();
}

void writeCellIds(
  const std::string &filename,
  const std::vector<std::string> &ids) {

    std::ofstream cellsof;
    // write cell ids, one line per id
    cellsof.open(filename.c_str(), std::ios::out);

    if (!cellsof.is_open()) {
      std::cerr << "Error: Couldn't open file: " << filename << std::endl;
      exit(1);
    }

    for (int j = 0; j < ids.size(); j++) {
      cellsof << ids[j]  << "\n";
    }
    cellsof.close();
}

void writeFLD(const std::string &filename, const std::vector<std::pair<double, double>> &flds) {
  std::ofstream of;
  of.open(filename.c_str(), std::ios::out);
  if (!of.is_open()) {
    std::cerr << "Error: Couldn't open file: " << filename << std::endl;
    exit(1);
  }
  for (int j = 0; j < flds.size(); j++) {
    of << j << "\t" << flds[j].first << "\t" << flds[j].second << "\n";
  }
  of.close();
}

void writeGeneList(const std::string &filename, const Transcriptome& model, bool writeNamesOnly) {
  std::ofstream of;
  of.open(filename.c_str(), std::ios::out);
  if (!of.is_open()) {
    std::cerr << "Error: Couldn't open file: " << filename << std::endl;
    exit(1);
  }
  for (int j = 0; j < model.genes.size(); j++) {
    if (writeNamesOnly) {
      of << model.genes[j].name << "\n";
    } else {
      of << j << "\t" << model.genes[j].name << "\t" << model.genes[j].commonName << "\n";
    }
  }
  of.close();
}
