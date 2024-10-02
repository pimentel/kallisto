#ifndef KALLISTO_PROCESSREADS_H
#define KALLISTO_PROCESSREADS_H

#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <iostream>
#include <fstream>

#include <thread>
#include <mutex>
#include <atomic>
#include <condition_variable>

#include "common.h"
#include "MinCollector.h"
#include "KmerIndex.h"
#include "PseudoBam.h"
#include "EMAlgorithm.h"
#include "GeneModel.h"
#include "BUSData.h"
#include "BUSTools.h"

#ifndef NO_HTSLIB
#include <htslib/kstring.h>
#include <htslib/sam.h>
#endif // NO_HTSLIB

class MasterProcessor;

int64_t ProcessReads(MasterProcessor& MP, const  ProgramOptions& opt);
int64_t ProcessBatchReads(MasterProcessor& MP, const ProgramOptions& opt);
int64_t ProcessBUSReads(MasterProcessor& MP, const ProgramOptions& opt);
std::pair<const_UnitigMap<Node>, int> findFirstMappingKmer(const std::vector<std::pair<const_UnitigMap<Node>, int>> &v);

class SequenceReader {
public:

  SequenceReader(const ProgramOptions& opt) :
  readbatch_id(-1) {};
  SequenceReader() : state(false), readbatch_id(-1) {};
  virtual ~SequenceReader() {}
//  SequenceReader(SequenceReader&& o);

  virtual bool empty() = 0;
  virtual void reset();
  virtual void reserveNfiles(int n) = 0;
  virtual bool fetchSequences(char *buf, const int limit, std::vector<std::pair<const char*, int>>& seqs,
                      std::vector<std::pair<const char*, int>>& names,
                      std::vector<std::pair<const char*, int>>& quals,
                      std::vector<uint32_t>& flags,
                      std::vector<std::string>& umis, int &readbatch_id,
                      bool full=false,
                      bool comments=false) = 0;


public:
  bool state; // is the file open
  int readbatch_id = -1;
};

class FastqSequenceReader : public SequenceReader {
public:

  FastqSequenceReader(const ProgramOptions& opt) : SequenceReader(opt),
  current_file(0), paired(!opt.single_end),
  f_umi(new std::ifstream{}) {
    SequenceReader::state = false;
    files = opt.files;

    interleave_nfiles = opt.input_interleaved_nfiles;
    if (opt.bus_mode) {
      nfiles = opt.busOptions.nfiles;
    } else {
      nfiles = paired ? 2 : 1;
    }
    if (interleave_nfiles != 0) { nfiles = 1; files.clear(); files.push_back(opt.files[0]); }
    reserveNfiles(nfiles);
  }
  FastqSequenceReader() : SequenceReader(),
  paired(false),
  f_umi(new std::ifstream{}),
  current_file(0), interleave_nfiles(0) {};
  FastqSequenceReader(FastqSequenceReader &&o);
  ~FastqSequenceReader();

  bool empty();
  void reset();
  void reserveNfiles(int n);
  bool fetchSequences(char *buf, const int limit, std::vector<std::pair<const char*, int>>& seqs,
                      std::vector<std::pair<const char*, int>>& names,
                      std::vector<std::pair<const char*, int>>& quals,
                      std::vector<uint32_t>& flags,
                      std::vector<std::string>& umis, int &readbatch_id,
                      bool full=false,
                      bool comments=false);

public:
  int nfiles = 1;
  uint32_t numreads = 0;
  std::vector<gzFile> fp;
  std::vector<int> l;
  std::vector<int> nl;
  bool paired;
  std::vector<std::string> files;
  std::unique_ptr<std::ifstream> f_umi;
  int current_file;
  std::vector<kseq_t*> seq;
  int interleave_nfiles;
};

#ifndef NO_HTSLIB
class BamSequenceReader : public SequenceReader {
public:

  BamSequenceReader(const ProgramOptions& opt) :
  SequenceReader(opt) {
    SequenceReader::state = true;

    fp = bgzf_open(opt.files[0].c_str(), "rb");
    head = bam_hdr_read(fp);
    rec = bam_init1();

    err = bam_read1(fp, rec);
    eseq = bam_get_seq(rec);
    l_seq = rec->core.l_qseq;

    bc = bam_aux2Z(bam_aux_get(rec, "CR"));
    l_bc = 0;
    for (char *pbc = bc; *pbc != '\0'; ++pbc) {
      ++l_bc;
    }

    umi = bam_aux2Z(bam_aux_get(rec, "UR"));
    l_umi = 0;
    for (char *pumi = umi; *pumi != '\0'; ++pumi) {
      ++l_umi;
    }
  }
  BamSequenceReader() : SequenceReader() {};
  BamSequenceReader(BamSequenceReader &&o);
  ~BamSequenceReader();

  bool empty();
  void reset();
  void reserveNfiles(int n);
  bool fetchSequences(char *buf, const int limit, std::vector<std::pair<const char*, int>>& seqs,
                      std::vector<std::pair<const char*, int>>& names,
                      std::vector<std::pair<const char*, int>>& quals,
                      std::vector<uint32_t>& flags,
                      std::vector<std::string>& umis, int &readbatch_id,
                      bool full=false,
                      bool comments=false);

public:
  BGZF *fp;
  bam_hdr_t *head;
  bam1_t *rec;
  uint8_t *eseq;
  int32_t l_seq;
  char *bc;
  int l_bc;
  char *umi;
  int l_umi;
  int err;


private:
  static const std::string seq_enc;
};
#endif

class MasterProcessor {
public:
  MasterProcessor (KmerIndex &index, const ProgramOptions& opt, MinCollector &tc, const Transcriptome& model)
    : tc(tc), index(index), model(model), opt(opt), numreads(0), transfer_threshold(1), counter(0)
    ,nummapped(0), num_umi(0), bufsize(1ULL<<23), tlencount(0), biasCount(0), maxBiasCount((opt.bias) ? 1000000 : 0), last_pseudobatch_id (-1) {

      #ifndef NO_HTSLIB
      bamfp = nullptr;
      bamfps = nullptr;
      bamh = nullptr;
      #endif // NO_HTSLIB
      if (opt.bam) {
        #ifndef NO_HTSLIB
        SR = new BamSequenceReader(opt);
        #else
        throw std::runtime_error("HTSLIB required for bam reading but not included");
        #endif // NO_HTSLIB
      } else {
        SR = new FastqSequenceReader(opt);
      }
      
      std::vector<std::mutex> mutexes(opt.threads);
      transfer_locks.swap(mutexes);

      if (opt.batch_mode) { // Set up recording of lengths individually for each batch
        if (opt.long_read) {
          tlencounts.assign(opt.batch_ids.size(), 0);
          //std::cerr << "Is this assignment of batchFlens space causing abort?" << std::endl; std::cerr.flush();
          batchFlens_lr.assign(opt.batch_ids.size(), std::vector<uint32_t>(index.target_lens_.size(),0));
          batchFlens_lr_c.assign(opt.batch_ids.size(), std::vector<uint32_t>(index.target_lens_.size(),0));
          //std::cerr << "Passed assignment of batchFlens" << std::endl; std::cerr.flush();
          batchUnmapped_list.assign(opt.batch_ids.size(), std::vector<double>(3000000/opt.batch_ids.size(),0));
        } else {
          tlencounts.assign(opt.batch_ids.size(), 0);
          batchFlens.assign(opt.batch_ids.size(), std::vector<uint32_t>(1000,0));
        }
      } 
      if (opt.batch_ids.size() > 0) {
        std::unordered_map<std::string,int> batch_map;
        batch_id_mapping.resize(opt.batch_ids.size());
        int j = 0;
        for (int i = 0; i < opt.batch_ids.size(); i++) {
          if (batch_map.find(opt.batch_ids[i]) == batch_map.end()) {
            batch_id_mapping[i] = j;
            batch_map.insert(std::make_pair(opt.batch_ids[i], j));
            j++;
          } else {
            batch_id_mapping[i] = batch_map[opt.batch_ids[i]];
          }
        }
      }
      if (opt.fusion) {
        ofusion.open(opt.output + "/fusion.txt");
        ofusion << "TYPE\tNAME1\tSEQ1\tKPOS1\tNAME2\tSEQ2\tKPOS2\tINFO\tPOS1\tPOS2\n";
      }
      if (opt.long_read) {
        ofusion.open(opt.output + "/novel.fastq");
      }
      if (opt.pseudobam) {
        pseudobatchf_out.open(opt.output + "/pseudoaln.bin", std::ios::out | std::ios::binary);
      }
      if (opt.bus_mode || opt.batch_mode) {
        memset(&bus_bc_len[0],0,sizeof(bus_bc_len));
        memset(&bus_umi_len[0],0,sizeof(bus_umi_len));
        busf_out.open(opt.output + "/output.bus", std::ios::out | std::ios::binary);

        if (opt.batch_mode) {
          if (opt.technology.empty()) { // no_technology : No UMIs, no barcodes, so just have batch=barcode
            writeBUSHeader(busf_out, BUSFORMAT_FAKE_BARCODE_LEN, 1);
          } else if (opt.record_batch_bus_barcode) { // We want to record batch in barcode
            if (opt.busOptions.bc[0].fileno != -1) // We want to push batch at front of extracted barcode
              writeBUSHeader(busf_out, opt.busOptions.getBCLength(), opt.busOptions.getUMILength());
            else // We don't have an extracted barcode so barcode=batch
              writeBUSHeader(busf_out, BUSFORMAT_FAKE_BARCODE_LEN, opt.busOptions.getUMILength());
          } else { // We don't care about recording batch in barcode, so we follow the default technology specification
            writeBUSHeader(busf_out, opt.busOptions.getBCLength(), opt.busOptions.getUMILength());
          }
        } else { // No batches, so we follow the default technology specification
          writeBUSHeader(busf_out, opt.busOptions.getBCLength(), opt.busOptions.getUMILength());
        }
      }
    }

  ~MasterProcessor() {
    #ifndef NO_HTSLIB
    if (bamfp) {
      hts_close(bamfp);
      bamfp = nullptr;
    }
    if (bamfps) {
      for (int i = 0; i < numSortFiles; i++) {
        if (bamfps[i]) {
          hts_close(bamfps[i]);
          bamfps[i] = nullptr;
        }
      }
      bamfps = nullptr;
    }
    if (bamh) {
      bam_hdr_destroy(bamh);
      bamh = nullptr;
    }
    #endif // NO_HTSLIB
    delete SR;
  }

  std::mutex reader_lock;
  std::vector<std::mutex> parallel_bus_reader_locks;
  bool parallel_bus_read;
  std::mutex writer_lock;
  std::vector<std::mutex> transfer_locks; // one mutex per thread; for locking each thread's access to index.ecmapinv
  size_t transfer_threshold;


  SequenceReader *SR;
  std::vector<FastqSequenceReader> FSRs;
  MinCollector& tc;
  KmerIndex& index;
  const Transcriptome& model;
  const int numSortFiles = 32;

  const ProgramOptions& opt;
  int64_t numreads;
  int64_t nummapped;
  int64_t num_umi;
  int64_t counter;
  size_t bufsize;

  int bus_bc_len[33];
  int bus_umi_len[33];

  std::atomic<int> tlencount;
  std::vector<int> tlencounts;
  std::atomic<int> biasCount;
  std::vector<std::vector<double>> batchUnmapped_list; 
  std::vector<std::vector<uint32_t>> batchFlens;
  std::vector<std::vector<uint32_t>> batchFlens_lr;
  std::vector<std::vector<uint32_t>> batchFlens_lr_c;
  std::vector<std::vector<int32_t>> tmp_bc;
  const int maxBiasCount;
  u_map_<Roaring, uint32_t, RoaringHasher> newECcount;
  //  std::vector<std::pair<BUSData, std::vector<int32_t>>> newB;
  u_map_<Roaring, int32_t, RoaringHasher> bus_ecmapinv;
  std::vector<int> batch_id_mapping; // minimal perfect mapping of batch ID

  std::ofstream ofusion;
  std::ofstream onovel;
  std::ofstream pseudobatchf_out;
  std::ofstream busf_out;
  std::ifstream pseudobatchf_in;
  std::vector<PseudoAlignmentBatch> pseudobatch_stragglers;
  int last_pseudobatch_id;
  void outputFusion(const std::stringstream &o);
  void outputNovel(const std::stringstream &o);
  void processReads();
  #ifndef NO_HTSLIB
  htsFile *bamfp;
  htsFile **bamfps;
  bam_hdr_t *bamh;
  void processAln(const EMAlgorithm& em, bool useEM);
  void writePseudoBam(const std::vector<bam1_t> &bv);
  void writeSortedPseudobam(const std::vector<std::vector<bam1_t>> &bvv);
  #endif
  std::vector<uint64_t> breakpoints;
  void update(const std::vector<uint32_t>& c, const std::vector<Roaring>& newEcs, std::vector<std::pair<Roaring, std::string>>& ec_umi, std::vector<std::pair<Roaring, std::string>> &new_ec_umi, int n, std::vector<int>& flens, std::vector<double>& unmapped_list, std::vector<int>& flens_lr, std::vector<int>& flens_lr_c, std::vector<int> &bias, const PseudoAlignmentBatch& pseudobatch, std::vector<BUSData> &bv, std::vector<std::pair<BUSData, Roaring>> newB, int *bc_len, int *umi_len,   int id = -1, int local_id = -1);
};

class ReadProcessor {
public:
  ReadProcessor(const KmerIndex& index, const ProgramOptions& opt, const MinCollector& tc, MasterProcessor& mp, int id = -1, int local_id = -1);
  ReadProcessor(ReadProcessor && o);
  ~ReadProcessor();
  char *buffer;

  size_t bufsize;
  bool paired;
  const MinCollector& tc;
  std::vector<std::pair<Roaring, std::string>> ec_umi;
  std::vector<std::pair<Roaring, std::string>> new_ec_umi;
  const KmerIndex& index;
  MasterProcessor& mp;
  FastqSequenceReader batchSR;
  int64_t numreads;
  int id;
  int local_id;
  PseudoAlignmentBatch pseudobatch;

  std::vector<std::pair<const char*, int>> seqs;
  std::vector<std::pair<const char*, int>> names;
  std::vector<std::pair<const char*, int>> quals;
  std::vector<uint32_t> flags;
  std::vector<std::string> umis;
  std::vector<Roaring> newEcs;
  std::vector<double> unmapped_list; 
  std::vector<int> flens;
  std::vector<int> flens_lr;
  std::vector<int> flens_lr_c;
  std::vector<int> bias5;

  std::vector<uint32_t> counts;

  void operator()();
  void processBuffer();
  void clear();
};

class BUSProcessor {
public:
  BUSProcessor(/*const*/ KmerIndex& index, const ProgramOptions& opt, const MinCollector& tc, MasterProcessor& mp, int id = -1, int local_id = -1);
  BUSProcessor(BUSProcessor && o);
  ~BUSProcessor();
  char *buffer;

  size_t bufsize;
  bool paired;
  bool bam;
  bool num;
  const MinCollector& tc;
  const KmerIndex& index;
  MasterProcessor& mp;
  int64_t numreads;
  int id;
  int local_id;
  PseudoAlignmentBatch pseudobatch;
  FastqSequenceReader batchSR;

  int bc_len[33];
  int umi_len[33];

  std::vector<std::pair<const char*, int>> seqs;
  std::vector<std::pair<const char*, int>> names;
  std::vector<std::pair<const char*, int>> quals;
  std::vector<std::string> umis;
  std::vector<uint32_t> flags;

  std::vector<Roaring> newEcs;
  std::vector<double> unmapped_list; 
  std::vector<int> flens;
  std::vector<int> flens_lr;
  std::vector<int> flens_lr_c;
  std::vector<int> bias5;
  std::vector<uint32_t> counts;
  std::vector<BUSData> bv;
  std::vector<std::pair<BUSData, Roaring>> newB;

  void operator()();
  void processBuffer();
  void clear();
};

#ifndef NO_HTSLIB
class AlnProcessor {
public:
  AlnProcessor(const KmerIndex& index, const ProgramOptions& opt, MasterProcessor& mp, const EMAlgorithm& em, const Transcriptome& model, bool useEM, int id = -1);
  AlnProcessor(AlnProcessor && o);
  ~AlnProcessor();
  char *buffer;
  char *bambuffer;
  size_t bufsize;
  size_t bambufsize;
  bool paired;
  std::vector<std::pair<int, std::string>> ec_umi;
  const KmerIndex& index;
  const EMAlgorithm& em;
  MasterProcessor& mp;
  FastqSequenceReader batchSR;
  int64_t numreads;
  int id;
  PseudoAlignmentBatch pseudobatch;
  const Transcriptome& model;
  bool useEM;



  std::vector<std::pair<const char*, int>> seqs;
  std::vector<std::pair<const char*, int>> names;
  std::vector<std::pair<const char*, int>> quals;
  std::vector<uint32_t> flags;
  std::vector<std::string> umis;

  void operator()();
  void processBufferTrans();
  void processBufferGenome();
  void clear();
};


int fillBamRecord(bam1_t &b, uint8_t* buf, const char *seq, const char *name, const char *qual, int slen, int nlen, bool unmapped, int auxlen);
void fixCigarStringTrans(bam1_t &b, int rlen, int softclip, int overhang);
void fixCigarStringGenome(bam1_t &b, const TranscriptAlignment& tra);
void reverseComplementSeqInData(bam1_t &b);
#endif // NO_HTSLIB

#endif // KALLISTO_PROCESSREADS_H
