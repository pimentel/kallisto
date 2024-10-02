#include <algorithm>
#include <random>
#include <sstream>
#include <ctype.h>
#include <unordered_set>
#include <functional>
#include "common.h"
#include "KmerIndex.h"
#include "SparseVector.hpp"
#include <iostream>
#include <unordered_map>
#include <string>
#include "ColoredCDBG.hpp"

// --aa option helper functions
// first three letters of nucleotide seq -> comma-free nuc_seq
// the three stop codons will be translated to 'NNN'
constexpr const char * cfc_map(const char * nuc_seq){
    // Must use ternary operator because < C++14 allows only a single return statement (and nothing else) in a constexpr
  return (
    (nuc_seq[0] == 'T' && nuc_seq[1] == 'T' && nuc_seq[2] == 'T') ? "ACC" : 
    (nuc_seq[0] == 'T' && nuc_seq[1] == 'T' && nuc_seq[2] == 'C') ? "ACC" : 
    (nuc_seq[0] == 'T' && nuc_seq[1] == 'T' && nuc_seq[2] == 'A') ? "ACA" :
    (nuc_seq[0] == 'T' && nuc_seq[1] == 'T' && nuc_seq[2] == 'G') ? "ACA" :
    (nuc_seq[0] == 'C' && nuc_seq[1] == 'T' && nuc_seq[2] == 'T') ? "ACA" :
    (nuc_seq[0] == 'C' && nuc_seq[1] == 'T' && nuc_seq[2] == 'C') ? "ACA" :
    (nuc_seq[0] == 'C' && nuc_seq[1] == 'T' && nuc_seq[2] == 'A') ? "ACA" :
    (nuc_seq[0] == 'C' && nuc_seq[1] == 'T' && nuc_seq[2] == 'G') ? "ACA" :
    (nuc_seq[0] == 'A' && nuc_seq[1] == 'T' && nuc_seq[2] == 'T') ? "ATA" :
    (nuc_seq[0] == 'A' && nuc_seq[1] == 'T' && nuc_seq[2] == 'C') ? "ATA" :
    (nuc_seq[0] == 'A' && nuc_seq[1] == 'T' && nuc_seq[2] == 'A') ? "ATA" :
    (nuc_seq[0] == 'A' && nuc_seq[1] == 'T' && nuc_seq[2] == 'G') ? "ATC" :
    (nuc_seq[0] == 'G' && nuc_seq[1] == 'T' && nuc_seq[2] == 'T') ? "ATT" :
    (nuc_seq[0] == 'G' && nuc_seq[1] == 'T' && nuc_seq[2] == 'C') ? "ATT" :
    (nuc_seq[0] == 'G' && nuc_seq[1] == 'T' && nuc_seq[2] == 'A') ? "ATT" :
    (nuc_seq[0] == 'G' && nuc_seq[1] == 'T' && nuc_seq[2] == 'G') ? "ATT" :
    (nuc_seq[0] == 'T' && nuc_seq[1] == 'C' && nuc_seq[2] == 'T') ? "CTA" :
    (nuc_seq[0] == 'T' && nuc_seq[1] == 'C' && nuc_seq[2] == 'C') ? "CTA" :
    (nuc_seq[0] == 'T' && nuc_seq[1] == 'C' && nuc_seq[2] == 'A') ? "CTA" :
    (nuc_seq[0] == 'T' && nuc_seq[1] == 'C' && nuc_seq[2] == 'G') ? "CTA" :
    (nuc_seq[0] == 'A' && nuc_seq[1] == 'G' && nuc_seq[2] == 'T') ? "CTA" :
    (nuc_seq[0] == 'A' && nuc_seq[1] == 'G' && nuc_seq[2] == 'C') ? "CTA" :
    (nuc_seq[0] == 'C' && nuc_seq[1] == 'C' && nuc_seq[2] == 'T') ? "CTC" :
    (nuc_seq[0] == 'C' && nuc_seq[1] == 'C' && nuc_seq[2] == 'C') ? "CTC" :
    (nuc_seq[0] == 'C' && nuc_seq[1] == 'C' && nuc_seq[2] == 'A') ? "CTC" :
    (nuc_seq[0] == 'C' && nuc_seq[1] == 'C' && nuc_seq[2] == 'G') ? "CTC" :
    (nuc_seq[0] == 'A' && nuc_seq[1] == 'C' && nuc_seq[2] == 'T') ? "CTT" :
    (nuc_seq[0] == 'A' && nuc_seq[1] == 'C' && nuc_seq[2] == 'C') ? "CTT" :
    (nuc_seq[0] == 'A' && nuc_seq[1] == 'C' && nuc_seq[2] == 'A') ? "CTT" :
    (nuc_seq[0] == 'A' && nuc_seq[1] == 'C' && nuc_seq[2] == 'G') ? "CTT" :
    (nuc_seq[0] == 'G' && nuc_seq[1] == 'C' && nuc_seq[2] == 'T') ? "AGA" :
    (nuc_seq[0] == 'G' && nuc_seq[1] == 'C' && nuc_seq[2] == 'C') ? "AGA" :
    (nuc_seq[0] == 'G' && nuc_seq[1] == 'C' && nuc_seq[2] == 'A') ? "AGA" :
    (nuc_seq[0] == 'G' && nuc_seq[1] == 'C' && nuc_seq[2] == 'G') ? "AGA" :
    (nuc_seq[0] == 'T' && nuc_seq[1] == 'A' && nuc_seq[2] == 'T') ? "AGC" :
    (nuc_seq[0] == 'T' && nuc_seq[1] == 'A' && nuc_seq[2] == 'C') ? "AGC" :
    (nuc_seq[0] == 'C' && nuc_seq[1] == 'A' && nuc_seq[2] == 'T') ? "AGT" :
    (nuc_seq[0] == 'C' && nuc_seq[1] == 'A' && nuc_seq[2] == 'C') ? "AGT" :
    (nuc_seq[0] == 'C' && nuc_seq[1] == 'A' && nuc_seq[2] == 'A') ? "AGG" :
    (nuc_seq[0] == 'C' && nuc_seq[1] == 'A' && nuc_seq[2] == 'G') ? "AGG" :
    (nuc_seq[0] == 'A' && nuc_seq[1] == 'A' && nuc_seq[2] == 'T') ? "CGA" :
    (nuc_seq[0] == 'A' && nuc_seq[1] == 'A' && nuc_seq[2] == 'C') ? "CGA" :
    (nuc_seq[0] == 'A' && nuc_seq[1] == 'A' && nuc_seq[2] == 'A') ? "CGC" :
    (nuc_seq[0] == 'A' && nuc_seq[1] == 'A' && nuc_seq[2] == 'G') ? "CGC" :
    (nuc_seq[0] == 'G' && nuc_seq[1] == 'A' && nuc_seq[2] == 'T') ? "CGT" :
    (nuc_seq[0] == 'G' && nuc_seq[1] == 'A' && nuc_seq[2] == 'C') ? "CGT" :
    (nuc_seq[0] == 'G' && nuc_seq[1] == 'A' && nuc_seq[2] == 'A') ? "CGG" :
    (nuc_seq[0] == 'G' && nuc_seq[1] == 'A' && nuc_seq[2] == 'G') ? "CGG" :
    (nuc_seq[0] == 'T' && nuc_seq[1] == 'G' && nuc_seq[2] == 'T') ? "TGA" :
    (nuc_seq[0] == 'T' && nuc_seq[1] == 'G' && nuc_seq[2] == 'C') ? "TGA" :
    (nuc_seq[0] == 'T' && nuc_seq[1] == 'G' && nuc_seq[2] == 'G') ? "TGC" :
    (nuc_seq[0] == 'C' && nuc_seq[1] == 'G' && nuc_seq[2] == 'T') ? "TGT" :
    (nuc_seq[0] == 'C' && nuc_seq[1] == 'G' && nuc_seq[2] == 'C') ? "TGT" :
    (nuc_seq[0] == 'C' && nuc_seq[1] == 'G' && nuc_seq[2] == 'A') ? "TGT" :
    (nuc_seq[0] == 'C' && nuc_seq[1] == 'G' && nuc_seq[2] == 'G') ? "TGT" :
    (nuc_seq[0] == 'A' && nuc_seq[1] == 'G' && nuc_seq[2] == 'A') ? "TGT" :
    (nuc_seq[0] == 'A' && nuc_seq[1] == 'G' && nuc_seq[2] == 'G') ? "TGT" :
    (nuc_seq[0] == 'G' && nuc_seq[1] == 'G' && nuc_seq[2] == 'T') ? "TGG" :
    (nuc_seq[0] == 'G' && nuc_seq[1] == 'G' && nuc_seq[2] == 'C') ? "TGG" :
    (nuc_seq[0] == 'G' && nuc_seq[1] == 'G' && nuc_seq[2] == 'A') ? "TGG" :
    (nuc_seq[0] == 'G' && nuc_seq[1] == 'G' && nuc_seq[2] == 'G') ? "TGG" :
    "NNN"
    );
};

int countNonAA=0;
std::string AA_to_cfc (const std::string aa_string) {
  // rev translate AA sequence to comma-free code (cfc)
  std::stringstream all_stream_aa;
  int n = aa_string.size();
  for (int i = 0; i < n; i++) {
      // map amino acid to comma-free code using cfc_aa_map (in common.cpp)
      char aa = aa_string[i];
      aa = ::toupper(aa);

      auto cfc_aa_mapped = cfc_aa_map.find(aa);
      std::string cfc_aa_seq;
      // if AA not found in comma-free map, translate as "NNN"
      if (cfc_aa_mapped == cfc_aa_map.end()) {
        cfc_aa_seq = "NNN";
        ::countNonAA++;
      } else {
        cfc_aa_seq = cfc_aa_mapped->second;
      }

      // accumulate comma-free sequences into stream
      all_stream_aa << cfc_aa_seq;
  }

  // convert stream to new comma-free AA sequence in string 'str'
  std::string str = all_stream_aa.str();

  return str;
}

// int countNonNN=0;
std::string nn_to_cfc (const char * s, int l) {
  // translate nucleotide sequence s to comma-free code

  // reserve memory for new cfc string
  std::string s_cfc;
  s_cfc.reserve(1024);
  // traverse the sequence in triplets
  int incrementer = 3;
  for (int i = 0; i < l; i += incrementer) {
      if (l - i >= incrementer) {
        // add comma-free triplet to cfc string
        char bytes[3];
        memset(bytes,0,3);
        bytes[0] = ::toupper(s[i]);
        bytes[1] = ::toupper(s[i+1]);
        bytes[2] = ::toupper(s[i+2]);
        s_cfc += cfc_map(bytes);
      }
  }
  return s_cfc;
}

// other helper functions
// pre: u is sorted
bool isUnique(const std::vector<int>& u) {
  for (int j = 1; j < u.size(); j++) {
    if (u[j-1] == u[j]) {
      return false;
    }
  }
  return true;
}

std::vector<int> unique(const std::vector<int>& u) {
  std::vector<int> v;
  v.reserve(u.size());
  v.push_back(u[0]);
  for (int j = 1; j < u.size(); j++) {
    if (u[j-1] != u[j]) {
      v.push_back(u[j]);
    }
  }
  return v;
}

const char Dna(int i) {
  static const char *dna = "ACGT";
  return dna[i & 0x03];
}

int hamming(const char *a, const char *b) {
  int h = 0;
  while (*a != 0 && *b != 0) {
    if (*a != *b) {
      h++;
    }
    a++;
    b++;
  }
  return h;
}

int my_mkdir_kmer_index(const char *path, mode_t mode) {
#ifdef _WIN64
  return mkdir(path);
#else
  return mkdir(path,mode);
#endif
}

std::string generate_tmp_file(std::string seed, std::string tmp_dir) {
  struct stat stFileInfo;
  auto intStat = stat(tmp_dir.c_str(), &stFileInfo);
  if (intStat == 0) {
    // file/dir exits
    if (!S_ISDIR(stFileInfo.st_mode)) {
      std::cerr << "Error: file " << tmp_dir << " exists and is not a directory" << std::endl;
      exit(1);
    }
  } else {
    // create directory
    if (my_mkdir_kmer_index(tmp_dir.c_str(), 0777) == -1) {
      std::cerr << "Error: could not create directory " << tmp_dir << std::endl;
      exit(1);
    }
  }
  std::string base = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
  std::string tmp_file = "kallisto.";
  srand((unsigned int)std::hash<std::string>{}(seed));
  int pos;
  while(tmp_file.length() < 32) {
    pos = ((rand() % (base.size() - 1)));
    tmp_file += base.substr(pos, 1);
  }
  return tmp_dir + "/" + tmp_file;
}

std::pair<size_t,size_t> KmerIndex::getECInfo() const {
  size_t max_ec_len = 0;
  size_t cardinality_zero_encounters = 0;
  for (const auto& um : dbg) {
    const Node* n = um.getData();
    std::vector<SparseVector<uint32_t> > vs;
    n->ec.get_vals(vs);
    bool cardinality_zero = true;
    for (auto &v : vs) {
      size_t cardinality = v.cardinality();
      if (cardinality > max_ec_len) {
        max_ec_len = cardinality;
      }
      if (cardinality > 0) cardinality_zero = false;
    }
    if (cardinality_zero) cardinality_zero_encounters++;
  }
  return std::make_pair(max_ec_len, cardinality_zero_encounters);
}

void KmerIndex::BuildTranscripts(const ProgramOptions& opt, std::ofstream& out) {
  // read input
  u_set_<std::string> unique_names;

  k = opt.k;
  for (auto& fasta : opt.transfasta) {
    std::cerr << "[build] loading fasta file " << fasta
              << std::endl;
  }
  std::cerr << "[build] k-mer length: " << k << std::endl;

  // Generate random file name
  std::string tmp_file = generate_tmp_file(opt.index, opt.tmp_dir);
  std::ofstream of(tmp_file);
  num_trans = 0;

  // read fasta file using kseq (https://lh3lh3.users.sourceforge.net/kseq.shtml)
  gzFile fp = 0;
  kseq_t *seq;
  int l = 0;
  std::mt19937 gen(42);
  int countNonNucl = 0;
  int countUNuc = 0;
  int polyAcount = 0;

  for (auto& fasta : opt.transfasta) {
    fp = gzopen(fasta.c_str(), "r");
    seq = kseq_init(fp);

    if (opt.aa) {
      while (true) {
        l = kseq_read(seq);
        if (l <= 0) {
          break;
        }
        
        // Translate amino acid (AA) sequence to comma-free code (cfc)
        std::string str = AA_to_cfc (seq->seq.s);

        of << ">" << num_trans++ << "\n" << str << std::endl;
        // record length of sequence after translating to cfc (will be 3x length of AA seq)
        target_lens_.push_back(str.size());
        // record sequence name
        std::string name(seq->name.s);
        size_t p = name.find(' ');
        if (p != std::string::npos) {
          name = name.substr(0,p);
        }

        if (unique_names.find(name) != unique_names.end()) {
          if (!opt.make_unique) {
            std::cerr << "Error: repeated name in FASTA file " << fasta << "\n" << name << "\n\n" << "Run with --make-unique to replace repeated names with unique names" << std::endl;
            exit(1);
          } else {
            for (int i = 1; ; i++) { // potential bug if you have more than 2^32 repeated names
              std::string new_name = name + "_" + std::to_string(i);
              if (unique_names.find(new_name) == unique_names.end()) {
                name = new_name;
                break;
              }
            }
          }
        }

        unique_names.insert(name);
        target_names_.push_back(name);
      }
    }

    else {
      while (true) {
        l = kseq_read(seq);
        if (l <= 0) {
          break;
        }
        std::string str = seq->seq.s;
        auto n = str.size();
        for (auto i = 0; i < n; i++) {
          char c = str[i];
          c = ::toupper(c);
          if (c=='U') {
            str[i] = 'T';
            countUNuc++;
          } else if (c !='A' && c != 'C' && c != 'G' && c != 'T') {
            str[i] = Dna(gen()); // replace with pseudorandom string
            countNonNucl++;
          }
        }
        std::transform(str.begin(), str.end(),str.begin(), ::toupper);

        if (str.size() >= 10 && str.substr(str.size()-10,10) == "AAAAAAAAAA") {
          // clip off polyA tail
          //std::cerr << "[index] clipping off polyA tail" << std::endl;
          polyAcount++;
          int j;
          for (j = str.size()-1; j >= 0 && str[j] == 'A'; j--) {}
          str = str.substr(0,j+1);
        }
        of << ">" << num_trans++ << "\n" << str << std::endl;

        target_lens_.push_back(seq->seq.l);
        std::string name(seq->name.s);
        size_t p = name.find(' ');
        if (p != std::string::npos) {
          name = name.substr(0,p);
        }

        if (unique_names.find(name) != unique_names.end()) {
          if (!opt.make_unique) {
            std::cerr << "Error: repeated name in FASTA file " << fasta << "\n" << name << "\n\n" << "Run with --make-unique to replace repeated names with unique names" << std::endl;
            exit(1);
          } else {
            for (int i = 1; ; i++) { // potential bug if you have more than 2^32 repeated names
              std::string new_name = name + "_" + std::to_string(i);
              if (unique_names.find(new_name) == unique_names.end()) {
                name = new_name;
                break;
              }
            }
          }
        }
        unique_names.insert(name);
        target_names_.push_back(name);
      }
    }

    gzclose(fp);
    fp=0;
  }

  of.close();

  if (polyAcount > 0) {
    std::cerr << "[build] warning: clipped off poly-A tail (longer than 10)" << std::endl << "        from " << polyAcount << " target sequences" << std::endl;
  }

  if (countNonNucl > 0) {
    std::cerr << "[build] warning: replaced " << countNonNucl << " non-ACGUT characters in the input sequence" << std::endl << "        with pseudorandom nucleotides" << std::endl;
  }

  if (countUNuc > 0) {
    std::cerr << "[build] warning: replaced " << countUNuc << " U characters with Ts" << std::endl;
  }

  if (countNonAA > 0) {
    std::cerr << "[build] warning: found " << countNonAA << " non-standard amino acid characters in the input sequence" << std::endl << "        which were reverse translated to 'NNN'" << std::endl;
  }

  BuildDeBruijnGraph(opt, tmp_file, out);
  BuildEquivalenceClasses(opt, tmp_file);

  std::remove(tmp_file.c_str());
}

void KmerIndex::BuildDistinguishingGraph(const ProgramOptions& opt, std::ofstream& out) {
  k = opt.k;
  std::cerr << "[build] k-mer length: " << k << std::endl;
  size_t ncolors = 0;
  std::string tmp_file2 = generate_tmp_file("--" + opt.index, opt.tmp_dir);
  // Use an external input FASTA file (we'll still need to read it to determine number of targets though and we'll still write out a new FASTA file (TODO: optimize this out)
  std::cerr << "[build] Reading in FASTA file" << std::endl;
  std::ofstream of(tmp_file2); // Write external FASTA file into another (possibly temporary) file
  gzFile fp = 0;
  kseq_t *seq;
  int l = 0;
  size_t num_seqs = 0;
  int max_color = 0;
  u_set_<int> external_input_names;
  for (int i = 0; i < opt.transfasta.size(); i++) { // Currently, this should only be one file
    auto fasta = opt.transfasta[i];
    fp = opt.transfasta.size() == 1 && opt.transfasta[0] == "-" ? gzdopen(fileno(stdin), "r") : gzopen(fasta.c_str(), "r");
    seq = kseq_init(fp);
    while (true) {
      l = kseq_read(seq);
      if (l <= 0) {
        break;
      }
      std::string str = seq->seq.s;
      std::string strname = seq->name.s;
      if (strname.length() == 0) {
        continue;
      }
      int color = std::atoi(strname.c_str());
      external_input_names.insert(color);
      if (color > max_color) max_color = color;
      of << ">" << std::to_string(color) << "\n" << str << std::endl;
      num_seqs++;
    }
    gzclose(fp);
    fp=0;
  }
  of.close();
  ncolors = external_input_names.size();
  std::cerr << "[build] Read in " << num_seqs << " sequences" << std::endl;
  std::cerr << "[build] Detected " << ncolors << " used colors" << std::endl;
  ncolors = max_color+1; // +1 because color id is zero-indexed
  std::cerr << "[build] Detected " << ncolors << " total colors" << std::endl;

  // Prepare some variables:
  num_trans = ncolors;
  target_names_.clear();
  target_lens_.clear();
  for (size_t i = 0; i < num_trans; i++) {
    target_names_.push_back(std::to_string(i));
    target_lens_.push_back(k); // dummy length (k-mer size)
  }

  std::cerr << "[build] Building graph from k-mers" << std::endl;
  BuildDeBruijnGraph(opt, tmp_file2, out);
  
  std::cerr << "[build] creating equivalence classes ... " << std::endl;
  
  UnitigMap<Node> um;
  uint32_t sense = 0x80000000, missense = 0;
  
  std::vector<std::vector<TRInfo> > trinfos(dbg.size());
  std::ifstream infile_a(tmp_file2);
  int current_color = 0;
  std::string line;
  while (std::getline(infile_a, line)) {
    if (line.length() == 0) {
      continue;
    } else if (line[0] == '>') {
      if (line.size() >= 1 && strncmp(line.c_str()+1, "d_list.", 7) == 0) {
        current_color = onlist_sequences.cardinality();
      } else {
        current_color = std::atoi(line.c_str()+1);
      }
      continue;
    }
    const auto& seq = line;
    if (seq.size() < k) { continue; }
    int seqlen = seq.size() - k + 1; // number of k-mers
    size_t proc = 0;
    while (proc < seqlen) {
      um = dbg.findUnitig(seq.c_str(), proc, seq.size());

      if (um.isEmpty) {
        ++proc;
        continue;
      }

      proc += um.len;
      const Node* n = um.getData();

      TRInfo tr;
      tr.trid = current_color;
      tr.pos = (proc-um.len) | (!um.strand ? sense : missense);
      tr.start = um.dist;
      tr.stop  = um.dist + um.len;

      trinfos[n->id].push_back(tr);
    }
  }
  infile_a.close();
  PopulateMosaicECs(trinfos);
  std::remove(tmp_file2.c_str());
  
  std::cerr << "[build] target de Bruijn graph has k-mer length " << dbg.getK() << " and minimizer length "  << dbg.getG() << std::endl;
  std::cerr << "[build] target de Bruijn graph has " << dbg.size() << " contigs and contains "  << dbg.nbKmers() << " k-mers " << std::endl;

}

void KmerIndex::BuildDeBruijnGraph(const ProgramOptions& opt, const std::string& tmp_file, std::ofstream& out) {

  CDBG_Build_opt c_opt;
  c_opt.k = k;
  c_opt.nb_threads = opt.threads;
  c_opt.build = true;
  c_opt.clipTips = false;
  c_opt.deleteIsolated = false;
  c_opt.verbose = true;
  c_opt.filename_ref_in.push_back(tmp_file);

  if (opt.g > 0) { // If minimizer length supplied, override the default
    c_opt.g = opt.g;
  } else { // Define minimizer length defaults
    int g = k-8;
    if (k <= 13) {
      g = k-2;
    } else if (k <= 17) {
      g = k-4;
    } else if (k <= 19) {
      g = k-6;
    }
    c_opt.g = g;
  }
  dbg = CompactedDBG<Node>(k, c_opt.g);
  dbg.build(c_opt);

  // If off-list is supplied, add off-listed kmers flanking the common
  // sequences to the graph and append those sequences to the tmp_file
  onlist_sequences = Roaring();
  onlist_sequences.addRange(0, num_trans);
  u_set_<Kmer, KmerHash> dfks;
  DListFlankingKmers(opt, tmp_file, dfks);

  // 1. write version
  out.write((char *)&INDEX_VERSION, sizeof(INDEX_VERSION));

  // 2. serialize dBG
  auto pos1 = out.tellp();
  // Write dummy size of graph
  size_t tmp_size = 1337;
  out.write((char *)&tmp_size, sizeof(tmp_size));
  bool res = dbg.writeBinary(out, opt.threads);

  if (res == 0) {
    std::cerr << "Error: could not write de Bruijn Graph to disk." << std::endl;
    exit(1);
  }

  auto pos2 = out.tellp();
  out.seekp(pos1);

  // Write real size of graph
  tmp_size = pos2 - pos1 - sizeof(tmp_size);
  out.write((char *)&tmp_size, sizeof(tmp_size));
  out.seekp(pos2);

  std::vector<Minimizer> minz;
  dbg.clearAndGetMinimizers(minz);
  std::cerr << "[build] building MPHF" << std::endl;
  boophf_t* mphf = new boophf_t(minz.size(), std::move(minz), opt.threads, 2.0, false, 0.15);
  out.write((char *)&tmp_size, sizeof(tmp_size));
  mphf->save(out);
  pos1 = out.tellp();
  out.seekp(pos2);
  tmp_size = pos1 - pos2 - sizeof(tmp_size);
  out.write((char *)&tmp_size, sizeof(tmp_size));

  out.seekp(pos1);

  dbg.clear();
  
  uint64_t dfks_size = dfks.size();
  uint64_t dlist_overhang = opt.d_list_overhang;
  out.write((char*)&dfks_size, sizeof(dfks_size));
  out.write((char*)&dlist_overhang, sizeof(dlist_overhang));
  if (dfks.size() > 0) {
    bool found_dummy_dfk = false;
    out.write((char*)&dummy_dfk, sizeof(dummy_dfk)); // The first DFK is the dummy
    for (auto dfk : dfks) {
      if (dfk != dummy_dfk) {
        out.write((char*)&dfk, sizeof(dfk));
      } else {
        found_dummy_dfk = true;
      }
    }
    if (!found_dummy_dfk) {
      // Should never happen
      std::cerr << "Assertion failed: Dummy DFK Not Found" << std::endl;
      exit(1);
    }
  }

  std::ifstream infile;
  infile.open(opt.index, std::ios::in | std::ios::binary);
  std::istream in(0);
  in.rdbuf(infile.rdbuf());
  in.ignore(sizeof(INDEX_VERSION));
  in.ignore(sizeof(tmp_size));

  dbg = CompactedDBG<Node>(k, c_opt.g);

  dbg.readBinary(in, mphf, opt.threads);

  //infile.close();

  uint32_t running_id = 0;
  for (auto& um : dbg) {
    um.getData()->id = running_id++;
  }
}

void KmerIndex::DListFlankingKmers(const ProgramOptions& opt, const std::string& tmp_file, u_set_<Kmer, KmerHash>& kmers) {

  if (opt.d_list.empty()) return;
  
  auto kmers_special = kmers; // copy; kmers_special will be the ones we store as-is (i.e. not flanking); these are the entries where the FASTA name header is not populated

  std::cerr << "[build] extracting D-list k-mers from";
  for (std::string s : opt.d_list) std::cerr << " \"" << s << "\""; 
  std::cerr << std::endl;
  
  auto isInvalidKmer = [](const char* s, const int k) {
      int count_nonATCG = 0;
      const int max_count_nonATCG = 3;
      for (int i = 0; i < k; i++) {
        if (s[i] != 'A' && s[i] != 'T' && s[i] != 'C' && s[i] != 'G') {
          ++count_nonATCG;
          if (count_nonATCG > max_count_nonATCG) { return true; }
        }
      }
      return false;
  };

  unsigned long overhang = opt.d_list_overhang;
  // Main worker thread
  auto worker_function = [&, overhang](std::string& seq, bool is_special) {

    int lb = -1, ub = -1;
    int pos = 0;
    
    std::transform(seq.begin(), seq.end(),seq.begin(), ::toupper);
    u_set_<Kmer, KmerHash> kmers_;
    if (is_special) { // Sequence has no name; it is "special"
      int seqlen = seq.size() - k + 1; // number of k-mers
      while (pos < seqlen) {
        if (!isInvalidKmer(seq.c_str()+pos, k)) kmers_.emplace(seq.c_str()+pos);
        ++pos;
      }
      return kmers_; // Return the special k-mers
    }
    const_UnitigMap<Node> um;
    int seqlen = seq.size() - k + 1; // number of k-mers
    while (pos < seqlen) {
      um = dbg.findUnitig(seq.c_str(), pos, seq.size());

      if (um.isEmpty) {

        // Add leading kmer to set
        if (lb >= 0 && ub >= lb) {
          for (int i = 0; i < std::min<unsigned long>(lb, overhang); ++i) {
            // Add up to #overhang leading k-mers to set
            if (!isInvalidKmer(seq.c_str() + lb - i, k)) kmers_.emplace(seq.c_str() + lb - i);
          }
        }

        // Add trailing kmer to set
        if (ub > lb && ub + k < seq.length()) {
          for (int i = 0; i < std::min<unsigned long>(seq.length() - ub, overhang); ++i) {
            // Add up to #overhang trailing k-mers to set
            if (!isInvalidKmer(seq.c_str() + ub + i, k)) kmers_.emplace(seq.c_str() + ub + i);
          }
        }

        lb = -1;
        ub = -1;
        ++pos;

      } else {

        if (lb < 0 && pos > 0) {
          lb = pos - 1;
        }

        pos += um.len;
        ub = pos;
      }

    }

    // Add last leading kmer to set
    if (lb >= 0 && ub >= lb) {
      for (int i = 0; i < std::min<unsigned long>(lb, overhang); ++i) {
        // Add up to #overhang leading k-mers to set
        if (!isInvalidKmer(seq.c_str() + lb - i, k)) kmers_.emplace(seq.c_str() + lb - i);
      }
    }

    // Add last trailing kmer to set
    if (ub > lb && ub + k < seq.length()) {
      for (int i = 0; i < std::min<unsigned long>(seq.length() - ub, overhang); ++i) {
        // Add up to #overhang trailing k-mers to set
        if (!isInvalidKmer(seq.c_str() + ub + i, k)) kmers_.emplace(seq.c_str() + ub + i);
      }
    }

    return kmers_;
  };

  // FASTA reading for D-list
  std::vector<std::string> dlist_fasta_files;
  dlist_fasta_files = opt.d_list;
  gzFile fp = 0;
  kseq_t *seq;
  size_t max_threads_read = std::min(opt.threads, 8);
  size_t n = 0;
  int l = 0;
  size_t rlen = 0;
  const size_t max_rlen = 640000;

  // Translate nucleotide dlist genomes to comma-free code in all possible frames
  std::string cfc_str_f1;
  std::string cfc_str_f2;
  std::string cfc_str_f3;
  std::string cfc_str_f4;
  std::string cfc_str_f5;
  std::string cfc_str_f6;
  if (opt.aa) {
    std::string tmp_file = generate_tmp_file("aa" + opt.index, opt.tmp_dir);
    std::ofstream of(tmp_file);
    size_t i = 0;
    for (auto& fasta : dlist_fasta_files) {
      fp = gzopen(fasta.c_str(), "r");
      seq = kseq_init(fp);
      while (true) {
        l = kseq_read(seq);
        if (l <= 0) {
          break;
        }

        const char * sequence = seq->seq.s;
        bool is_special = (seq->name.l == 0);

        // Forward frame 1
        size_t seqlen1 = seq->seq.l;
        // Translate to comma-free code
        cfc_str_f1 = nn_to_cfc(sequence, seqlen1);
        // Write translated sequence to temporary file
        of << ">" << (!is_special ? std::to_string(i++) : "") << "\n" << cfc_str_f1 << std::endl;

        // Forward frame 2
        const char * seq2 = sequence+1;
        cfc_str_f2 = nn_to_cfc(seq2, seqlen1 - 1);
        of << ">" << (!is_special ? std::to_string(i++) : "") << "\n" << cfc_str_f2 << std::endl;

        // Forward frame 3
        const char * seq3 = sequence+2;
        cfc_str_f3 = nn_to_cfc(seq3, seqlen1 - 2);
        of << ">" << (!is_special ? std::to_string(i++) : "") << "\n" << cfc_str_f3 << std::endl;

        // Get reverse complement of sequence
        // const char * to string
        std::string com_seq(sequence);
        // transform comseq to its reverse complement
        com_seq = revcomp (std::move(com_seq));
        // string to const char *
        const char * com_seq_char = com_seq.c_str();

        // Rev comp frame 1
        cfc_str_f4 = nn_to_cfc(com_seq_char, seqlen1);
        of << ">" << (!is_special ? std::to_string(i++) : "") << "\n" << cfc_str_f4 << std::endl;

        // Rev comp frame 2
        const char * seq5 = com_seq_char+1;
        cfc_str_f5 = nn_to_cfc(seq5, seqlen1-1);
        of << ">" << (!is_special ? std::to_string(i++) : "") << "\n" << cfc_str_f5 << std::endl;

        // Rev comp frame 3
        const char * seq6 = com_seq_char+2;
        cfc_str_f6 = nn_to_cfc(seq6, seqlen1-2);
        of << ">" << (!is_special ? std::to_string(i++) : "") << "\n" << cfc_str_f6 << std::endl;
      }
      gzclose(fp);
      fp = 0;
    }
    of.close();
    // Overwrite dlist fastas with tmp file containing translated seqs
    dlist_fasta_files.clear();
    dlist_fasta_files.push_back(tmp_file);
  }

  for (auto& fasta : dlist_fasta_files) {
    std::vector<std::string > seqs_v(max_threads_read, "");
    std::vector<std::thread> workers;
    std::mutex mutex_file;
    std::mutex mutex_kmers;
    fp = gzopen(fasta.c_str(), "r");
    seq = kseq_init(fp);
    while (true) {
      l = kseq_read(seq);
      if (l <= 0) {
        break;
      }
      std::string sequence = seq->seq.s;
      auto slen = sequence.size();
      if (slen < max_rlen || max_threads_read <= 1 || seq->name.l == 0) { // Small enough (or "special" sequence); no need to create a new thread
        bool is_special = (seq->name.l == 0);
        auto kmers_local = worker_function(sequence, is_special);
        {
          // Preempting write access for kmer set
          std::unique_lock<std::mutex> lock(mutex_kmers);
          bool is_special = (seq->name.l == 0);
          for (auto& kmer : kmers_local) {
            if (!is_special) kmers.insert(kmer.rep());
            else kmers_special.insert(kmer.rep());
          }
        }
        continue;
      }
      seqs_v[n % seqs_v.size()] = std::move(sequence);
      workers.emplace_back(
        [&, n] {
          auto kmers_local = worker_function(seqs_v[n % seqs_v.size()], false);
          {
            // Preempting write access for kmer set
            std::unique_lock<std::mutex> lock(mutex_kmers);
            for (auto& kmer : kmers_local) {
              kmers.insert(kmer.rep());
            }
          }
        }
      );
      n++;
      if (n % seqs_v.size() == 0) {
        for (auto& t : workers) t.join();
        workers.clear();
      }
    }
    for (auto& t : workers) t.join(); // finish up all remaining threads
    gzclose(fp);
    fp = 0;
  }
  
  // Some clean ups:
  if (opt.aa) { // Get rid of the temp file for translated aa
    struct stat stFileInfo;
    std::string aa_tmp_file = generate_tmp_file("aa" + opt.index, opt.tmp_dir);
    auto intStat = stat(aa_tmp_file.c_str(), &stFileInfo);
    if (intStat == 0) {
      std::remove(aa_tmp_file.c_str());
    }
  }

  size_t N = 0;
  size_t start = num_trans;
  size_t unique_flanking_kmers = 0;
  size_t other_dlist_kmers = 0;
  std::ofstream outfile;
  outfile.open(tmp_file, std::ios_base::app);
  bool added_dummy_dfk = false;
  u_set_<Kmer, KmerHash> already_in_graph;
  for (const auto& kmer : kmers) {
    // Check whether flanking k-mer exists elsewhere in the graph
    // (may occur in the case of longer overhangs)
    //if (overhang > 1) {
    //  auto um = dbg.find(kmer);
    //  if (!um.isEmpty) continue;
    //}

    // Insert all other flanking kmers into graph
    const_UnitigMap<Node> um_ = dbg.find(kmer.rep());
    if (!um_.isEmpty) { already_in_graph.insert(kmer.rep()); continue; }
    if (!added_dummy_dfk) {
      // const_UnitigMap<Node> um_ = dbg.find(kmer);
      if (um_.isEmpty) {
        dbg.add(kmer.rep().toString());
        dummy_dfk = kmer.rep();
        added_dummy_dfk = true;
      }
    }

    // Insert all other flanking kmers into tmp_file and transcript-related member variables
    std::string tx_name = "d_list." + std::to_string(N++);

    ++num_trans;
    //target_names_.push_back(tx_name);
    //target_lens_.push_back(k);

    outfile << ">"
            << tx_name
            << std::endl
            << kmer.toString()
            << std::endl;

    ++unique_flanking_kmers;

  }
  for (const auto& kmer : kmers_special) { // Add these directly (not checking for flanking)
    if (kmers.find(kmer.rep()) != kmers.end()) continue; // Already found; just continue
    if (!added_dummy_dfk) {
      const_UnitigMap<Node> um_ = dbg.find(kmer.rep());
      if (um_.isEmpty) {
        dbg.add(kmer.rep().toString());
        dummy_dfk = kmer.rep();
        added_dummy_dfk = true;
      } else {
        dummy_dfk = kmer.rep(); // for special k-mers, it's ok to make a k-mer already in the graph the dummy
        added_dummy_dfk = true;
      }
    }
    kmers.insert(kmer.rep()); // Add to the master k-mer set
    already_in_graph.erase(kmer.rep()); // Erase from here if necessary (because these special D-listed k-mers are stringent filters)
    std::string tx_name = "d_list." + std::to_string(N++);
    ++num_trans;
    //target_names_.push_back(tx_name);
    //target_lens_.push_back(k);
    outfile << ">"
            << tx_name
            << std::endl
            << kmer.toString()
            << std::endl;
    ++other_dlist_kmers;
  }
  for (const auto& kmer : already_in_graph) {
    if (kmer != dummy_dfk) kmers.erase(kmer);
  }
  if (!added_dummy_dfk) kmers.clear();
  std::cerr << "[build] identified " << unique_flanking_kmers << " distinguishing flanking k-mers" << std::endl;
  if (other_dlist_kmers > 0) {
    std::cerr << "[build] identified " << other_dlist_kmers << " other D-list k-mers" << std::endl;
  }

  outfile.close();
}

void KmerIndex::BuildEquivalenceClasses(const ProgramOptions& opt, const std::string& tmp_file) {

  std::cerr << "[build] creating equivalence classes ... " << std::endl;

  std::vector<std::vector<TRInfo> > trinfos(dbg.size());
  UnitigMap<Node> um;
  size_t EC_THRESHOLD = 250;
  size_t EC_SOFT_THRESHOLD = 800;
  size_t EC_MAX_N_ABOVE_THRESHOLD = 6000; // Thresholding ECs to size EC_THRESHOLD will only occur if we encounter >EC_MAX_N_ABOVE_THRESHOLD number of nodes that have size >EC_SOFT_THRESHOLD
  if (opt.max_ec_size > 0) { // If max EC size is supplied, override the default thresholds
    EC_THRESHOLD = opt.max_ec_size;
    EC_SOFT_THRESHOLD = EC_THRESHOLD;
    EC_MAX_N_ABOVE_THRESHOLD = 0;
  } else if (opt.max_ec_size <= -1) { // Default: no cap
    EC_THRESHOLD = std::numeric_limits<uint32_t>::max();
    EC_SOFT_THRESHOLD = EC_THRESHOLD;
    EC_MAX_N_ABOVE_THRESHOLD = 0;
  }
  uint32_t sense = 0x80000000, missense = 0;

  std::ifstream infile(tmp_file);
  std::string line;
  size_t j = 0;
  size_t n_above_threshold = 0;
  //for (size_t i = 0; i < seqs.size(); ++i) {
  while (std::getline(infile, line)) {
    if (line[0] == '>') continue;
    const auto& seq = line;
    if (seq.size() < k) { j++; continue; }

    int seqlen = seq.size() - k + 1; // number of k-mers
    size_t proc = 0;
    while (proc < seqlen) {
      um = dbg.findUnitig(seq.c_str(), proc, seq.size());

      if (um.isEmpty) {
        ++proc;
        continue;
      }

      proc += um.len;
      const Node* n = um.getData();
      if (trinfos[n->id].size() > EC_THRESHOLD) {
        if (trinfos[n->id].size() == EC_SOFT_THRESHOLD+1) {
          n_above_threshold++;
        }
        if (n_above_threshold > EC_MAX_N_ABOVE_THRESHOLD) {
          trinfos[n->id].clear();
          std::vector<TRInfo>().swap(trinfos[n->id]); // potentially free up memory
          TRInfo tr_discard;
          tr_discard.trid = std::numeric_limits<uint32_t>::max();
          trinfos[n->id].reserve(1);
          trinfos[n->id].push_back(tr_discard);
          continue;
        }
      } else if (trinfos[n->id].size() == 1 && trinfos[n->id][0].trid == std::numeric_limits<uint32_t>::max()) {
        continue;
      }
      TRInfo tr;

      tr.trid = std::min<size_t>(j, onlist_sequences.cardinality());
      tr.pos = (proc-um.len) | (!um.strand ? sense : missense);
      tr.start = um.dist;
      tr.stop  = um.dist + um.len;

      trinfos[n->id].push_back(tr);
    }
    j++;
  }
  infile.close();

  // Threshold large ECs
  if (n_above_threshold > EC_MAX_N_ABOVE_THRESHOLD) {
    size_t n_removed = 0;
    for (auto& trinfo : trinfos) {
      if (trinfo.size() > EC_THRESHOLD) {
        trinfo.clear();
        std::vector<TRInfo>().swap(trinfo); // potentially free up memory
        ++n_removed;
      } else if (trinfo.size() == 1 && trinfo[0].trid == std::numeric_limits<uint32_t>::max()) {
        trinfo.clear();
        ++n_removed;
      }
    }
    std::cerr << "[build] discarded " << n_removed << " ECs larger than threshold." << std::endl;
  }

  PopulateMosaicECs(trinfos);

  std::cerr << "[build] target de Bruijn graph has k-mer length " << dbg.getK() << " and minimizer length "  << dbg.getG() << std::endl;
  std::cerr << "[build] target de Bruijn graph has " << dbg.size() << " contigs and contains "  << dbg.nbKmers() << " k-mers " << std::endl;
}

void KmerIndex::PopulateMosaicECs(std::vector<std::vector<TRInfo> >& trinfos) {

  for (const auto& um : dbg) {

    Node* n = um.getData();

    // Process empty ECs
    if (trinfos[n->id].size() == 0) {
      SparseVector<uint32_t> u(true);
      n->ec.insert(0, um.len, std::move(u));
      continue;
    }

    // Find the overlaps
    std::vector<int> brpoints;
    brpoints.reserve(2 * trinfos[n->id].size());
    for (const auto& x : trinfos[n->id]) {
      brpoints.push_back(x.start);
      brpoints.push_back(x.stop);
    }

    sort(brpoints.begin(), brpoints.end());
    assert(brpoints[0] == 0);
    assert(brpoints[brpoints.size()-1]==um.size-k+1);

    // Find unique break points
    if (!isUnique(brpoints)) {
      std::vector<int> u = unique(brpoints);
      swap(u,brpoints);
    }

    std::sort(trinfos[n->id].begin(), trinfos[n->id].end(),
              [](const TRInfo& lhs, const TRInfo& rhs) -> bool {
                return (lhs.trid < rhs.trid);
              });

    size_t j = 0;
    // Create a mosaic EC for the unitig, where each break point interval
    // corresponds to one set of transcripts and therefore an EC
    for (size_t i = 1; i < brpoints.size(); ++i) {

      SparseVector<uint32_t> u(true);

      for (const auto& tr : trinfos[n->id]) {
        // If a transcript encompasses the full breakpoint interval
        if (tr.start <= brpoints[i-1] && tr.stop >= brpoints[i]) {
          u.insert(tr.trid, tr.pos);
        }
      }

      assert(!u.isEmpty());
      u.runOptimize();

      // Assign mosaic EC and transcript position+sense to the corresponding part of unitig
      n->ec.insert(brpoints[i-1], brpoints[i], std::move(u));
    }
    std::vector<TRInfo>().swap(trinfos[n->id]); // potentially free up memory
  }
}

void KmerIndex::write(std::ofstream& out, int threads) {

  size_t tmp_size;

  if (!out.is_open()) {
    // TODO: better handling
    std::cerr << "Error: index output file could not be opened!";
    exit(1);
  }

  // 3. serialize nodes
  tmp_size = dbg.size();
  out.write((char *)&tmp_size, sizeof(tmp_size));
  for (const auto& um : dbg) {
    // 3.1 write head kmer to associate unitig with node
    std::string kmer = um.getUnitigHead().toString();
    out.write(kmer.c_str(), strlen(kmer.c_str()));
    // 3.2 serialize node
    std::ostringstream out_;
    um.getData()->serialize(out_);
    uint32_t s_size = out_.str().length();
    out.write((char*)&s_size, sizeof(s_size));
    out.write(out_.str().c_str(), s_size);
  }

  // 4. write number of targets
  out.write((char *)&num_trans, sizeof(num_trans));

  // 5. write out target lengths
  for (int tlen : target_lens_) {
    out.write((char *)&tlen, sizeof(tlen));
  }

  // 6. Write out target ids
  // XXX: num_trans should equal to target_names_.size(), so don't need
  // to write out again.
  assert(num_trans == target_names_.size());
  for (auto& tid : target_names_) {
    // 6.1 write out how many bytes
    tmp_size = strlen(tid.c_str());
    out.write((char *)&tmp_size, sizeof(tmp_size));

    // 6.2 write out the actual string
    out.write(tid.c_str(), tmp_size);
  }

  // 7. Write on-list
  char* buffer = new char[onlist_sequences.getSizeInBytes()];
  tmp_size = onlist_sequences.write(buffer);
  out.write((char *)&tmp_size, sizeof(tmp_size));
  out.write(buffer, tmp_size);
  delete[] buffer;
  buffer = nullptr;

}

void KmerIndex::write(const std::string& index_out, bool writeKmerTable, int threads) {
  
  if (writeKmerTable) {
    std::cerr << "KmerIndex::write() must have writeKmerTable set to false" << std::endl;
    exit(1);
  }

  std::ofstream out;
  out.open(index_out, std::ios::out | std::ios::binary);

  if (!out.is_open()) {
    // TODO: better handling
    std::cerr << "Error: index output file could not be opened!";
    exit(1);
  }

  size_t tmp_size;

  // append is true if we have already written the INDEX_VERSION and the dBG
  // when constructing the dBG
  // 1. write version
  out.write((char *)&INDEX_VERSION, sizeof(INDEX_VERSION));

  // 2. serialize dBG
  if (writeKmerTable) {
    auto pos1 = out.tellp();
    // Write dummy size of graph
    tmp_size = 1337;
    out.write((char *)&tmp_size, sizeof(tmp_size));
    bool res = dbg.writeBinary(out, threads);

    if (res == 0) {
      std::cerr << "Error: could not write de Bruijn Graph to disk." << std::endl;
      exit(1);
    }

    auto pos2 = out.tellp();
    out.seekp(pos1);

    // Write real size of graph
    tmp_size = pos2 - pos1 - sizeof(tmp_size);
    out.write((char *)&tmp_size, sizeof(tmp_size));
    out.seekp(pos2);
  } else {
    tmp_size = 0;
    out.write((char *)&tmp_size, sizeof(tmp_size));
  }
  // Don't write out a D-list:
  uint64_t dlist_size = 0;
  uint64_t dlist_overhang = 1;
  out.write((char*)&dlist_size, sizeof(dlist_size));
  out.write((char*)&dlist_overhang, sizeof(dlist_overhang));

  // 3. serialize nodes
  if (writeKmerTable) {
    tmp_size = dbg.size();
    out.write((char *)&tmp_size, sizeof(tmp_size));
    for (const auto& um : dbg) {
      // 3.1 write head kmer to associate unitig with node
      std::ostringstream out_;
      um.getData()->serialize(out_);
      uint32_t s_size = out_.str().length();
      out.write((char*)&s_size, sizeof(s_size));
      out.write(out_.str().c_str(), s_size);
    }
  } else {
    tmp_size = 0;
    out.write((char *)&tmp_size, sizeof(tmp_size));
  }

  // 4. write number of targets
  num_trans -= d_list.size();
  out.write((char *)&num_trans, sizeof(num_trans));

  // 5. write out target lengths
  for (int tlen : target_lens_) {
    out.write((char *)&tlen, sizeof(tlen));
  }

  // 6. Write out target ids
  // XXX: num_trans should equal to target_names_.size(), so don't need
  // to write out again.
  assert(num_trans == target_names_.size());
  for (auto& tid : target_names_) {
    // 6.1 write out how many bytes
    tmp_size = strlen(tid.c_str());
    out.write((char *)&tmp_size, sizeof(tmp_size));

    // 6.2 write out the actual string
    out.write(tid.c_str(), tmp_size);
  }

  // 7. Write on-list
  char* buffer = new char[onlist_sequences.getSizeInBytes()];
  tmp_size = onlist_sequences.write(buffer);
  out.write((char *)&tmp_size, sizeof(tmp_size));
  out.write(buffer, tmp_size);
  delete[] buffer;
  buffer = nullptr;

  out.flush();
  out.close();
}

void KmerIndex::load(ProgramOptions& opt, bool loadKmerTable, bool loadDlist) {

  if (opt.index.empty() && !loadKmerTable) {
    // Make an index from transcript and EC files
    loadTranscriptsFromFile(opt);
    loadECsFromFile(opt);
    return;
  }

  std::string& index_in = opt.index;
  std::ifstream in;//, in_minz;

  in.open(index_in, std::ios::in | std::ios::binary);
  //in_minz.open(index_in, std::ios::in | std::ios::binary);

  if (!in.is_open()) {
    // TODO: better handling
    std::cerr << "Error: index input file could not be opened!";
    exit(1);
  }

  // 1. read version
  size_t header_version = 0;
  in.read((char *)&header_version, sizeof(header_version));
  //in_minz.ignore(sizeof(header_version));

  if (header_version != INDEX_VERSION) {
    std::cerr << "Error: incompatible indices. Found version " << header_version << ", expected version " << INDEX_VERSION << std::endl
              << "Rerun with index to regenerate";
    exit(1);
  }

  // 2. deserialize dBG
  size_t tmp_size;
  in.read((char *)&tmp_size, sizeof(tmp_size));
  tmp_size = ((-1ULL >> 1) & tmp_size); // Mask out MSB (in case we want to use it as a toggle in some future implementation)
  if (tmp_size > 0) {

    auto pos1 = in.tellg();
    in.ignore(tmp_size);
    boophf_t* mphf = new boophf_t();
    in.read((char *)&tmp_size, sizeof(tmp_size));
    mphf->load(in);
    auto pos2 = in.tellg();
    in.seekg(pos1);

    auto success = dbg.readBinary(in, mphf, opt.threads);

    in.seekg(pos2);

    //dbg.to_static();
    k = dbg.getK();
  }
  std::cerr << "[index] k-mer length: " << std::to_string(k) << std::endl;
  
  // 2.2 deserialize distinguishing flanking k-mers for D-list
  uint64_t dlist_size;
  uint64_t dlist_overhang;
  in.read((char*)&dlist_size, sizeof(dlist_size));
  in.read((char*)&dlist_overhang, sizeof(dlist_overhang));
  d_list.reserve(dlist_size);
  for (size_t i = 0; i < dlist_size; i++) {
    Kmer dfk;     
    in.read((char*)&dfk, sizeof(dfk));
    d_list.insert(dfk);   
    if (i == 0) {    
      dummy_dfk = dfk;
      um_dummy = dbg.find(dummy_dfk);
      if (um_dummy.isEmpty) {
        std::cerr << "Error: Dummy k-mer not found in graph" << std::endl;
        exit(1);
      }                            
    }                                          
  }

  // 3. deserialize nodes
  Kmer kmer;
  size_t kmer_size = k * sizeof(char);
  char* buffer = new char[kmer_size];
  in.read((char *)&tmp_size, sizeof(tmp_size));
  const size_t max_num_nodes_buffer = 524288;
  std::vector<std::pair<char*, std::pair<Kmer, uint32_t> > > in_buf_v;
  in_buf_v.reserve(max_num_nodes_buffer);
  std::vector<std::thread> workers;
  workers.reserve(opt.threads);
  for (size_t i = 0; i < tmp_size; ++i) {
    // 3.1 read head kmer
    memset(buffer, 0, kmer_size);
    in.read(buffer, kmer_size);
    kmer = Kmer(buffer);

    // 3.2 deserialize node
    uint32_t node_size;
    in.read((char *)&node_size, sizeof(node_size));
    if (opt.threads == 1) {
      UnitigMap<Node> um;
      um = dbg.find(kmer);
      if (um.isEmpty) {
        std::cerr << "Error: Corrupted index; unitig not found: " << std::string(buffer) << std::endl;
        exit(1);
      }
      um.getData()->deserialize(in, !load_positional_info); // Just one thread; read it directly
    } else { // multi-threaded
      char* node_buf = new char[node_size];
      in.read(node_buf, node_size);
      in_buf_v.push_back(std::make_pair(node_buf, std::make_pair(kmer, node_size)));
      if (in_buf_v.size() >= in_buf_v.capacity() || i == tmp_size-1) {
        for (size_t t = 0; t < opt.threads; t++) {
          workers.emplace_back([&, t] {
            for (size_t j = t; j < in_buf_v.size(); j += opt.threads) {
              char* node_content = in_buf_v[j].first;
              Kmer km = in_buf_v[j].second.first;
              auto node_size_ = in_buf_v[j].second.second;
              UnitigMap<Node> um;
              um = dbg.find(km);
              if (um.isEmpty) {
                std::cerr << "Error: Corrupted index; unitig not found: " << std::string(buffer) << std::endl;
                exit(1);
              }
              std::stringstream oss;
              oss.write(node_content, node_size_);
              delete[] node_content; // free memory associated with node_buf
              node_content = nullptr;
              std::istringstream iss;
              iss.basic_ios<char>::rdbuf(oss.rdbuf());
              um.getData()->deserialize(iss, !load_positional_info); // No need to lock because each node is necessarily unique
            }
          });
        }
        for (auto& t : workers) t.join();
        in_buf_v.clear();
        workers.clear();
      }
    }
  }
  delete[] buffer;
  buffer = nullptr;
  std::vector<std::pair<char*, std::pair<Kmer, uint32_t> > >().swap(in_buf_v); // potentially free up memory
  std::vector<std::thread>().swap(workers);

  // 4. read number of targets
  in.read((char *)&num_trans, sizeof(num_trans));
  num_trans -= dlist_size;

  // 5. read out target lengths
  target_lens_.clear();
  target_lens_.reserve(num_trans);
  int tmp_int;
  for (size_t i = 0; i < num_trans; ++i) {
    in.read((char *)&tmp_int, sizeof(tmp_int));
    target_lens_.push_back(tmp_int);
  }

  // 6. read in target ids
  target_names_.clear();
  target_names_.reserve(num_trans);

  size_t bufsz = 1024;
  buffer = new char[bufsz];
  for (auto i = 0; i < num_trans; ++i) {

    // 6.1 read in the size
    in.read((char *)&tmp_size, sizeof(tmp_size));

    if (tmp_size +1 > bufsz) {
      delete[] buffer;
      bufsz = 2*(tmp_size+1);
      buffer = new char[bufsz];
    }

    // clear the buffer
    memset(buffer,0,bufsz);
    // 6.2 read in the character string
    in.read(buffer, tmp_size);

    target_names_.push_back(std::string(buffer));
  }
  delete[] buffer;

  // 7. Read on-list
  in.read((char *)&tmp_size, sizeof(tmp_size));
  buffer = new char[tmp_size];
  in.read(buffer, tmp_size);
  onlist_sequences = onlist_sequences.read(buffer);

  // delete the buffer
  delete[] buffer;
  buffer=nullptr;
  
  num_trans += dlist_size;

  std::cerr << "[index] number of targets: " << pretty_num(static_cast<size_t>(onlist_sequences.cardinality())) << std::endl;
  std::cerr << "[index] number of k-mers: " << pretty_num(dbg.nbKmers()) << std::endl;
  if (num_trans-onlist_sequences.cardinality() > 0) {
    std::cerr << "[index] number of D-list k-mers: " << pretty_num(static_cast<size_t>(num_trans-onlist_sequences.cardinality())) << std::endl;
  }

  in.close();

  if (!opt.ecFile.empty()) {
    loadECsFromFile(opt);
  }
  
  if (!loadDlist) { // Destroy the D-list
    if (num_trans != onlist_sequences.cardinality()) {
      std::cerr << "[index] not using the D-list k-mers" << std::endl;
      num_trans = onlist_sequences.cardinality();
      target_names_.resize(num_trans);
      target_lens_.resize(num_trans);
    }
  }
}

void KmerIndex::loadECsFromFile(const ProgramOptions& opt) {
  ecmapinv.clear();
  int32_t i = 0;
  std::ifstream in((opt.ecFile));
  if (in.is_open()) {
    std::string line;
    while (getline(in, line)) {
      std::stringstream ss(line);
      int ec;
      std::string transcripts;
      ss >> ec >> transcripts;
      if (i != ec) {
        std::cerr << "Error: equivalence class file has a misplaced equivalence class."
                  << " Found " << ec << ", expected " << i << std::endl;
        exit(1);
      }

      Roaring r;
      std::stringstream ss2(transcripts);
      while(ss2.good()) {
        std::string tmp_ecval;
        getline(ss2, tmp_ecval, ',');
        int tmp_ecval_num = std::atoi(tmp_ecval.c_str());
        if (tmp_ecval_num < 0 || tmp_ecval_num >= num_trans) {
          std::cerr << "Error: equivalence class file has invalid value: "
                    << tmp_ecval << " in " << transcripts << std::endl;
          exit(1);
        }
        r.add(tmp_ecval_num);
      }
      ecmapinv.insert({std::move(r), i}); // move
      i++;
    }
  } else {
    std::cerr << "Error: could not open file " << opt.ecFile << std::endl;
    exit(1);
  }
  std::cerr << "[index] number of equivalence classes loaded from file: "
            << pretty_num(ecmapinv.size()) << std::endl;
}

void KmerIndex::loadTranscriptsFromFile(const ProgramOptions& opt) {
  target_names_.clear();
  int i = 0;
  std::ifstream in((opt.transcriptsFile));
  if (in.is_open()) {
    std::string txp;
    while (in >> txp) {
      target_names_.push_back(txp);
      i++;
    }
  } else {
    std::cerr << "Error: could not open file " << opt.transcriptsFile << std::endl;
    exit(1);
  }
  num_trans = i;
  target_lens_.assign(num_trans, 0);
  std::cerr << "[index] number of targets loaded from file: "
            << pretty_num(num_trans) << std::endl;
}

int KmerIndex::mapPair(const char *s1, int l1, const char *s2, int l2) const {
  bool d1 = true;
  bool d2 = true;
  int p1 = -1;
  int p2 = -1;
  int c1 = -1;
  int c2 = -1;

  KmerIterator kit1(s1), kit_end;
  const_UnitigMap<Node> um1, um2;

  bool found1 = false;
  for (; kit1 != kit_end; ++kit1) {
    um1 = dbg.find(kit1->first);
    if (!um1.isEmpty) {
      found1 = true;
      if (um1.strand)
        p1 = um1.dist - kit1->second;
      else
      if (um1.strand)
        p1 = um1.dist - kit1->second;
      else
        p1 = um1.dist + k + kit1->second;
      break;
    }
  }

  if (!found1) {
    return -1;
  }

  KmerIterator kit2(s2);
  bool found2 = false;

  for (; kit2 != kit_end; ++kit2) {
    um2 = dbg.find(kit2->first);
    if (!um2.isEmpty) {
      found2 = true;
      if (um2.strand)
        p2 = um2.dist - kit2->second;
      else
        p2 = um2.dist + k + kit2->second;
      break;
    }
  }

  if (!found2) {
    return -1;
  }

  // We want the reads to map within the same EC block on the same unitig
  if (!um1.isSameReferenceUnitig(um2) ||
      !(um1.getData()->ec[um1.dist] == um2.getData()->ec[um2.dist])) {
    return -1;
  }

  // Paired reads need to map to opposite strands
  if (!(um1.strand ^ um2.strand)) {
    //std::cerr << "Reads map to same strand " << s1 << "\t" << s2 << std::endl;
    return -1;
  }

  if (um1.getData()->get_mc_contig(um1.dist).second != um2.getData()->get_mc_contig(um2.dist).second) {
    return -1; // If the mc contigs for um1 and um2 are actually not the same (despite having the same color)
  }

  if (p1>p2) {
    return p1-p2;
  } else {
    return p2-p1;
  }
}

// use:  match(s,l,v)
// pre:  v is initialized
// post: v contains all equiv classes for the k-mers in s
void KmerIndex::match(const char *s, int l, std::vector<std::pair<const_UnitigMap<Node>, int>>& v, bool partial, bool cfc) const{
  const Node* n;

  // TODO:
  // Rework KmerIndex::match() such that it uses the following type of logic
  // rather than the jumping logic below

  /*
  size_t proc = 0;
  while (proc < l - k + 1) {
    const_UnitigMap<Node> um = dbg.findUnitig(s, proc, l);
    if (um.isEmpty) {
      proc++;
      continue;
    }

    n = um.getData();
    uint32_t curr_ec = n->ec[um.dist];
    v.emplace_back(um, proc);
    // Add one entry to v for each EC that is part of the mosaic EC of the contig.
    for (size_t i = 0; i < um.len; ++i) {
      if (n->ec[um.dist + i] != curr_ec) {
        curr_ec = n->ec[um.dist + i];
        v.emplace_back(dbg.find(um.getUnitigKmer(um.dist + i)), proc + i);
      }
    }
    proc += um.len;
  }
  */

  std::string s_string;
  if (cfc) {
    // translate nucleotide sequence to comma-free code (cfc)
    s_string = nn_to_cfc(s, l);
    s = s_string.c_str();

    // if (countNonNN > 0) {
    //   std::cerr << "[warning] found " << countNonNN << " non-standard nucleotides in the input sequence" << std::endl << "           which were translated to 'NNN'" << std::endl;
    // }
    // // reset countNonNN
    // ::countNonNN = 0;
  }

  KmerIterator kit(s), kit_end;
  bool backOff = false;
  int nextPos = 0; // nextPosition to check
  Roaring rtmp;
  for (int i = 0;  kit != kit_end; ++i,++kit) {

    // need to check it
    const_UnitigMap<Node> um = dbg.find(kit->first);
    n = um.getData();

    int pos = kit->second;

    if (!um.isEmpty) {
      
      if (partial) {
        if (rtmp.isEmpty()) {
          const auto& rtmp2 = um.getData()->ec[um.dist].getIndices();
          if (!rtmp2.isEmpty()) rtmp = std::move(rtmp2);
        } else {
          const auto& rtmp2 = um.getData()->ec[um.dist].getIndices();
          if (!rtmp2.isEmpty()) rtmp &= rtmp2;
          if (rtmp.isEmpty()) {
            v.clear();
            return;
          }
        }
      }

      v.push_back({um, kit->second});

      // Find start and end of O.G. kallisto contig w.r.t. the bifrost-kallisto
      // unitig
      size_t contig_start = 0, contig_length = um.size - k + 1;
      auto p = n->get_mc_contig(um.dist);
      contig_start += p.first;
      contig_length = p.second - contig_start;

      // Looks like kallisto thinks that canonical kmer means forward strand?
      //bool forward = (um.strand == (kit->first == kit->first.rep()));
      bool forward = um.strand;
      int dist = (forward) ? (contig_length - 1 - (um.dist - contig_start)) : um.dist - contig_start;

      // see if we can skip ahead
      if (dist >= 2) {
        // where should we jump to?
        int nextPos = pos+dist; // default jump

        if (pos + dist >= l-k) {
          // if we can jump beyond the read, check the end
          nextPos = l-k;
        }

        // check next position
        KmerIterator kit2(kit);
        kit2 += nextPos-pos;
        if (kit2 != kit_end) {
          const_UnitigMap<Node> um2 = dbg.find(kit2->first);
          bool found2 = false;
          int  found2pos = pos+dist;
          if (um2.isEmpty) {
            found2=true;
            found2pos = pos;
          } else if (um.isSameReferenceUnitig(um2) &&
                     n->ec[um.dist] == um2.getData()->ec[um2.dist]) {
            // um and um2 are on the same unitig and also share the same EC
            found2=true;
            found2pos = pos+dist;
          }
          if (found2) {
            // great, a match (or nothing) see if we can move the k-mer forward
            bool is_in_dlist = um2.isEmpty && d_list.size() > 0 && (d_list.find((kit2->first).rep()) != d_list.end()); // D-list (check after first jump)
            if (found2pos >= l-k) {
              v.push_back({um, l-k}); // push back a fake position
              if (partial && is_in_dlist) { v.push_back({um_dummy, kit2->second}); return; } // D-list
              break; //
            } else {
              v.push_back({um, found2pos});
              kit = kit2; // move iterator to this new position
              if (partial && is_in_dlist) { v.push_back({um_dummy, kit2->second}); return; } // D-list
            }
          } else {
            // this is weird, let's try the middle k-mer
            bool foundMiddle = false;
            if (dist > 4) {
              int middlePos = (pos + nextPos)/2;
              int middleContig = -1;
              int found3pos = pos+dist;
              KmerIterator kit3(kit);
              kit3 += middlePos-pos;

              if (kit3 != kit_end) {
                const_UnitigMap<Node> um3 = dbg.find(kit3->first);
                if (!um3.isEmpty) {
                  if (um.isSameReferenceUnitig(um3) &&
                      n->ec[um.dist] == um3.getData()->ec[um3.dist]) {
                    foundMiddle = true;
                    found3pos = middlePos;
                  } else if (um2.isSameReferenceUnitig(um3) &&
                             um2.getData()->ec[um2.dist] == um3.getData()->ec[um3.dist]) {
                    foundMiddle = true;
                    found3pos = pos+dist;
                  }
                }

                if (foundMiddle) {
                  if (partial) {
                    if (rtmp.isEmpty()) {
                      const auto& rtmp2 = um3.getData()->ec[um3.dist].getIndices();
                      if (!rtmp2.isEmpty()) rtmp = std::move(rtmp2);
                    } else {
                      const auto& rtmp2 = um3.getData()->ec[um3.dist].getIndices();
                      if (!rtmp2.isEmpty()) rtmp &= rtmp2;
                      if (rtmp.isEmpty()) {
                        v.clear();
                        return;
                      }
                    }
                  }
                  v.push_back({um3, found3pos});
                  if (nextPos >= l-k) {
                    break;
                  } else {
                    kit = kit2;
                  }
                }
              }
            }

            if (!foundMiddle) {
              ++kit;
              backOff = true;
              goto donejumping; // sue me Dijkstra!
            }
          }
        } else {
          // the sequence is messed up at this point, let's just take the match
          //v.push_back({dbGraph.ecs[val.contig], l-k});
          break;
        }
      }
    }

donejumping:

    if (backOff) {
      // backup plan, let's play it safe and search incrementally for the rest, until nextStop
      for (int j = 0; kit != kit_end; ++kit,++j) {
        if (j==skip) {
          j=0;
        }
        if (j==0) {
          // need to check it
          const_UnitigMap<Node> um4 = dbg.find(kit->first);
          if (!um4.isEmpty) {
            // if k-mer found
            if (partial) {
              if (rtmp.isEmpty()) {
                const auto& rtmp2 = um4.getData()->ec[um4.dist].getIndices();
                if (!rtmp2.isEmpty()) rtmp = std::move(rtmp2);
              } else {
                const auto& rtmp2 = um4.getData()->ec[um4.dist].getIndices();
                if (!rtmp2.isEmpty()) rtmp &= rtmp2;
                if (rtmp.isEmpty()) {
                  v.clear();
                  return;
                }
              }
            }
            v.push_back({um4, kit->second}); // add equivalence class, and position
          }
        }

        if (kit->second >= nextPos) {
          backOff = false;
          break; // break out of backoff for loop
        }
      }
    }
  }
  
  // D-list:
  KmerIterator kit_dlist(s);
  if (d_list.size() > 0 && (v.size() > 0 || !partial)) {
    for (int i = 0;  kit_dlist != kit_end; ++i,++kit_dlist) {
      // need to check it
      bool is_in_dlist = (d_list.find((kit_dlist->first).rep()) != d_list.end());
      if (is_in_dlist) {
        v.push_back({um_dummy, kit_dlist->second});
        break;
      }
    }
  }
}

// use:  match_long(s,l,v)
// pre:  v is initialized
// post: v contains all equiv classes for the k-mers in s
double KmerIndex::match_long(const char *s, int l, std::vector<std::pair<const_UnitigMap<Node>, int>>& v, bool partial, bool cfc) const {
  const Node* n;
  double unmapped_ratio;

  std::string s_string;
  if (cfc) {
    // translate nucleotide sequence to comma-free code (cfc)
    s_string = nn_to_cfc(s, l);
    s = s_string.c_str();
  }
  int empty_count = 0; 

  Roaring rtmp;
  KmerIterator kit(s), kit_end;
  size_t proc = 0;
  while (kit != kit_end) { //should be + 2?
    if (proc < l-k){
      const_UnitigMap<Node> umU = dbg.findUnitig(s, proc, l);
      if (!umU.isEmpty) {
        int numK = 0; 
        do {
          v.push_back({umU, proc});
          numK++;
        } while (numK < int(umU.len/k));
      }
    }

    const_UnitigMap<Node> um = dbg.find(kit->first);
    //const_UnitigMap<Node> um = dbg.findUnitig(s, proc, l);     

    n = um.getData();

    int pos = kit->second;
    if (um.isEmpty) { 
      empty_count++;
    } else {      
      if (partial) {
        if (rtmp.isEmpty()) {
          const auto& rtmp2 = um.getData()->ec[um.dist].getIndices();
          if (!rtmp2.isEmpty()) rtmp = std::move(rtmp2);
        } else {
          const auto& rtmp2 = um.getData()->ec[um.dist].getIndices();
          if (!rtmp2.isEmpty()) rtmp &= rtmp2;
          if (rtmp.isEmpty()) {
            v.clear();
            unmapped_ratio = (double)empty_count/(double)(l - k);
            return unmapped_ratio;
          }
        }
      }

      v.push_back({um, kit->second});

      // Find start and end of O.G. kallisto contig w.r.t. the bifrost-kallisto
      // unitig
      size_t contig_start = 0, contig_length = um.size - k + 1;
      auto p = n->get_mc_contig(um.dist);
      contig_start += p.first;
      contig_length = p.second - contig_start;

      // Looks like kallisto thinks that canonical kmer means forward strand?
      //bool forward = (um.strand == (kit->first == kit->first.rep()));
      bool forward = um.strand;
      int dist = (forward) ? (contig_length - 1 - (um.dist - contig_start)) : um.dist - contig_start;

      // see if we can skip ahead
      if (dist >= 2) {
        // where should we jump to?
        int nextPos = pos+dist; // default jump

        if (pos + dist >= l-k) { //should be + 1?
          // if we can jump beyond the read, check the end
          nextPos = l-k; //should be +1?
        }

        // check next position
        KmerIterator kit2(kit);
        kit2 += nextPos-pos;
        if (kit2 != kit_end) { //(nextPos < l-k) { //should be +1?
          const_UnitigMap<Node> um2 = dbg.find(kit2->first); 
	        //const_UnitigMap<Node> um2 = dbg.findUnitig(s, nextPos, l); 
            
          bool found2 = false;
          int  found2pos = pos+dist;
          if (um2.isEmpty) {
            found2 = true;
            found2pos = pos;
          } else if (um.isSameReferenceUnitig(um2) &&
                       n->ec[um.dist] == um2.getData()->ec[um2.dist]) {
              // um and um2 are on the same unitig and also share the same EC
            found2=true;
            found2pos = pos+dist;
          }
            
          if (found2) {
            // great, a match (or nothing) see if we can move the k-mer forward
            if (found2pos >= l-k) { //should be +1?
              v.push_back({um, l-k}); // push back a fake position //should be +2?
              break; //
            } else {
              v.push_back({um, found2pos});
              //trying incremental search
              proc=found2pos;
              kit = kit2; // move iterator to this new position
            }
          } else {
            // this is weird, let's try the middle k-mer
            bool foundMiddle = false;
            if (dist > 4) {
              int middlePos = (pos + nextPos)/2;
              int middleContig = -1;
              int found3pos = pos + dist; //middlePos+dist; //formerly pos+dist which is same as found2pos, but I think should be middlePos+dist
              KmerIterator kit3(kit);
              kit3 += middlePos-pos;

              if (kit3 != kit_end) { //(found3pos < l-k) {
                const_UnitigMap<Node> um3 = dbg.find(kit3->first); 
                //const_UnitigMap<Node> um3 = dbg.findUnitig(s, middlePos, l); 
                if (!um3.isEmpty) {
                  if (um.isSameReferenceUnitig(um3) &&
                      n->ec[um.dist] == um3.getData()->ec[um3.dist]) {
                    foundMiddle = true;
                    found3pos = middlePos;
                  } else if (um2.isSameReferenceUnitig(um3) &&
                               um2.getData()->ec[um2.dist] == um3.getData()->ec[um3.dist]) {
                    foundMiddle = true;
                    found3pos = pos+dist;
                  }
                }

                if (foundMiddle) {
                  if (partial) {
                    if (rtmp.isEmpty()) {
                      const auto& rtmp2 = um3.getData()->ec[um3.dist].getIndices();
                      if (!rtmp2.isEmpty()) rtmp = std::move(rtmp2);
                    } else {
                      const auto& rtmp2 = um3.getData()->ec[um3.dist].getIndices();
                      if (!rtmp2.isEmpty()) {
                        rtmp &= rtmp2;
                      }

                      if (rtmp.isEmpty()) {
                        v.clear();
                        unmapped_ratio = (double)empty_count/(double)(l - k);
                        return unmapped_ratio;
                      }
                    }
                  }
                  
                  v.push_back({um3, found3pos});
                  
                  if (nextPos >= l-k) { //should be +2?
                    break;
                  } else {
                    //trying incremental search
                    proc=found2pos;
                    kit = kit2;
                  }
                }
              }
            }

            if (!foundMiddle) {
              ++proc; 
	     	      ++kit;
              // backup plan, let's play it safe and search incrementally for the rest, until nextStop
              for (int j = 0; kit != kit_end; ++kit,++j) {
                if (j == skip) {
                  j = 0;
                }
                 
                if (j == 0) {
                  // need to check it
                  const_UnitigMap<Node> um4 = dbg.find(kit->first);
                  //const_UnitigMap<Node> um4 = dbg.findUnitig(s, proc, l);;
                  if (!um4.isEmpty) {
                    // if k-mer found
                    if (partial) {
                      if (rtmp.isEmpty()) {
                        const auto& rtmp2 = um4.getData()->ec[um4.dist].getIndices();
                        if (!rtmp2.isEmpty()) {
                          rtmp = std::move(rtmp2);
                        }
                      } else {
                        const auto& rtmp2 = um4.getData()->ec[um4.dist].getIndices();
                        if (!rtmp2.isEmpty()) {
                          rtmp &= rtmp2;
                        }
                        
                        if (rtmp.isEmpty()) {
                          v.clear();
                          unmapped_ratio = (double)empty_count/(double)(l - k);
                          return unmapped_ratio;
                        }
                      }    
                    }
                  
                    v.push_back({um4, proc}); // add equivalence class, and position
                  }
                }

                if (kit->second >= nextPos) {
                  //continue;
                  //trying out the not breaking case to do a complete incremental search until kit_end
                  break; 
                }
              ++proc;   
              }
            //kit += (l-k)-nextPos;
            //proc = l; //this is checking if jump logic works how I am suspecting it is working 
            } 
          }
        }
        //} //NOTE!!! When this if is outside of the check if nextPos is at the end of the kmer, then we allow for incrementally searching which is the correct behavior? 
        else {
          // the sequence is messed up at this point, let's just take the match
          //v.push_back({dbGraph.ecs[val.contig], l-k});
          break;
        }	
      }  //adding this corresponding to NOTE!!!
    } 
    kit++;
    proc++; 
  }
  
  // D-list:
  KmerIterator kit_dlist(s);
  if (d_list.size() > 0 && (v.size() > 0 || !partial)) {
    for (int i = 0;  kit_dlist != kit_end; ++i,++kit_dlist) {
      // need to check it
      bool is_in_dlist = (d_list.find((kit_dlist->first).rep()) != d_list.end());
      if (is_in_dlist) {
        v.push_back({um_dummy, kit_dlist->second});
        break;
      }
    }
  }

  //std::cerr << "l: " << l << " empty_count: " << empty_count << std::endl;
  unmapped_ratio = (double)empty_count/(double)(l - k);
  return unmapped_ratio;
}

std::pair<int,bool> KmerIndex::findPosition(int tr, Kmer km, int p) const{
  const_UnitigMap<Node> um = dbg.find(km);
  if (!um.isEmpty) {
    return findPosition(tr, km, um, p);
  } else {
    return {-1,true};
  }
}

//use:  (pos,sense) = index.findPosition(tr,km,val,p)
//pre:  index.kmap[km] == val,
//      km is the p-th k-mer of a read
//      val.contig maps to tr
//post: km is found in position pos (1-based) on the sense/!sense strand of tr
std::pair<int,bool> KmerIndex::findPosition(int tr, Kmer km, const_UnitigMap<Node>& um, int p) const {
  bool csense = um.strand;
  int trpos = -1;
  uint32_t bitmask = 0x7FFFFFFF, rawpos;
  bool trsense = true;
  if (um.getData()->id == -1) {
    return {-1, true};
  }
  const Node* n = um.getData();
  auto mc = n->get_mc_contig(um.dist);
  auto ecs = n->ec.get_leading_vals(um.dist);
  const auto& v_ec = ecs[ecs.size() - 1];
  const Roaring& ec = v_ec.getIndices(); // transcripts
  
  rawpos = v_ec.get(tr, true).minimum();
  trpos = rawpos & bitmask;
  trsense = (rawpos == trpos);
  auto um_dist = um.dist - mc.first; // for mosaic ECs, need to start from beginning of current block, not zero
  
  if (trpos == -1) {
    return {-1,true};
  }
  
  std::pair<int,bool> ret;
  
  if (trsense) {
    if (csense) {
      size_t padding = 0;
      if (trpos == 0) {
        for (int i = ecs.size()-2; i >= 0; i--) {
          if (!ecs[i].contains(tr)) {
            padding = mc.first;
            break;
          }
          mc = n->get_mc_contig(mc.first-1);
        }
      }
      ret = {static_cast<int64_t>(trpos - p + static_cast<int64_t>(um.dist) + 1 - padding), csense}; // Case I
    } else {
      int64_t padding = um.size-mc.second+1-k;
      int64_t initial = mc.second;
      
      int right_one = 0;
      int left_one = 0;
      for (int i = ecs.size()-1; i >= 0; i--) {
        if (i == ecs.size()-1) right_one = mc.second;
        if (!ecs[i].contains(tr)) { left_one = mc.second; break; }
        else if (i == 0) left_one = 0;
        mc = n->get_mc_contig(mc.first-1); // if goes below 0, it'll be cast as the largest unsigned int and return the right-most block
      }
      padding = -(left_one+right_one-um.size+k-1);
      ret = {trpos + p + k - (um.size - k - um.dist) + initial - 1 + padding, csense}; // Case III
    }
  } else {
    if (csense) {
      int64_t start = mc.first;
      auto mc_ = n->get_mc_contig(mc.first-1);
      start = 0;
      auto ecs_ = n->ec.get_leading_vals(mc.first-1);
      int curr_mc = 0;
      int left_one = 0;
      int right_one = 0;
      int unmapped_len = 0;
      bool found_first_mapped = false;
      
      ecs_ = n->ec.get_leading_vals(-1);
      for (int i = 0; i < ecs_.size(); i++) {
        auto mc__ = n->get_mc_contig(curr_mc);
        if ((!ecs_[i].contains(tr)) && found_first_mapped) {
          if (unmapped_len == 0) left_one = mc__.first;
          right_one = mc__.second;
          unmapped_len += (mc__.second - mc__.first);
        }
        if (ecs_[i].contains(tr)) found_first_mapped = true;
        curr_mc = mc__.second;
      }
      start -= right_one-left_one;
      start += um.size-k;
      ret = {trpos + (static_cast<int64_t>(-(um.dist-start))) + k + p, !csense}; // Case IV
    } else {
      int curr_mc = 0;
      int left_one = 0;
      int right_one = 0;
      auto ecs_ = n->ec.get_leading_vals(-1);
      int unmapped_len = 0;
      bool found_first_mapped = false;
      for (int i = 0; i < ecs_.size(); i++) {
        auto mc__ = n->get_mc_contig(curr_mc);
        if ((!ecs_[i].contains(tr)) && found_first_mapped) {
          if (unmapped_len == 0) left_one = mc__.first;
          right_one = mc__.second;
          unmapped_len += (mc__.second - mc__.first);
        }
        if (ecs_[i].contains(tr)) {
          found_first_mapped = true;
        }
        curr_mc = mc__.second;//+1;
      }
      unmapped_len = right_one - left_one;
      int64_t padding = um.size-um.dist-(unmapped_len)-k+1;
      ret = {trpos+padding-p, !csense}; // case II
    }
  }
  return ret;
}

// use:  res = intersect(ec,v)
// pre:  ec is an equivalence class, v is a vector of valid targets
//       v is sorted in increasing order
// post: res contains the intersection  of ec and v
Roaring KmerIndex::intersect(const Roaring& ec, const Roaring& v) const {
  Roaring res;
  if (ec.cardinality() == 0) {
    // If the EC has no transcripts, it has been filtered due to its size
    res = v;
  } else if (v.cardinality() == 0) {
    // If transcript vector is empty, it represents an EC that has been filtered
    // due to its size
    res = ec;
  } else {
    // Do an actual intersect
    res = ec & v;
  }

  return res;
}

void KmerIndex::loadTranscriptSequences() const {
  if (target_seqs_loaded) {
    return;
  }

  bool &t = const_cast<bool&>(target_seqs_loaded);
  t = true;//target_seqs_loaded = true;
  return;
}

void KmerIndex::clear() {
  dbg.clear();

  {
    EcMapInv empty;
    std::swap(ecmapinv, empty);
  }

  target_lens_.resize(0);
  target_names_.resize(0);
  target_seqs_.resize(0);
}

void KmerIndex::writePseudoBamHeader(std::ostream &o) const {
  // write out header
  o << "@HD\tVN:1.0\n";
  for (int i = 0; i < num_trans; i++) {
    o << "@SQ\tSN:" << target_names_[i] << "\tLN:" << target_lens_[i] << "\n";
  }
  o << "@PG\tID:kallisto\tPN:kallisto\tVN:"<< KALLISTO_VERSION << "\n";
  o.flush();
}
