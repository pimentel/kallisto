#ifndef KALLISTO_INSPECTINDEX_H
#define KALLISTO_INSPECTINDEX_H

#include <iostream>
#include <unordered_set>
#include <limits>
#include "common.h"
#include "KmerIndex.h"
#include "GeneModel.h"


struct ECStruct {
  int ec;
  int chr;
  int start;
  int stop;
  std::vector<std::pair<int,int>> start_lens;
  std::vector<int> tlist;
};

std::vector<ECStruct> merge_contigs(std::vector<ECStruct> ecv) {
  if (ecv.size() <= 1) {
    return (ecv);
  }
  std::vector<ECStruct> out;
  std::sort(ecv.begin(), ecv.end(), [&](const ECStruct& a, const ECStruct& b) { return a.start < b.start;});

  // check if overlapping

  int a = 0, b = 0;
  while (b <= ecv.size()) {
    assert(a <= b);
    if (a == ecv.size()) {
      break;
    }

    // ecv[a:b] can be merged, see if ecv[a:(b+1)] can be
    if (b < ecv.size()) {
      if (a == b) {
        b++;
        continue; // ok, trivial to merge empty set
      }

      if (ecv[b-1].stop <= ecv[b].start) {
        if (b+1 == ecv.size() || ecv[b].stop <= ecv[b+1].start) {
          b++;
          continue;
        }
      }
    }
    // ok, push back ecv[a:b]
    if (a+1 == b) {
      out.push_back(ecv[a]);
    } else {
      std::unordered_set<int> tset;
      ECStruct ecs;
      ecs.ec = ecv[a].ec;
      ecs.chr = ecv[a].ec;
      ecs.start = ecv[a].start;
      ecs.stop = ecv[b-1].stop;
      int pos,len;
      for (int i = a; i < b; i++) {
        // len = 0;
        pos = ecv[i].start - ecs.start;
        auto& sp = ecs.start_lens;
        for (auto x : ecv[i].start_lens) {
          sp.push_back({pos+x.first, x.second});
          // len += x.second;
        }
        tset.insert(ecv[i].tlist.begin(), ecv[i].tlist.end());
      }
      for (auto t : tset) {
        ecs.tlist.push_back(t);
      }
      std::sort(ecs.tlist.begin(), ecs.tlist.end());
      out.push_back(ecs);
    }
    //
    a = b;
  }

  return out;
}

void printVector(const std::vector<int>& v) {
  std::cout << "[";
  int i = 0;
  for (auto x : v) {
    if (i>0) {
      std::cout << ", ";
    }
    std::cout << x;
    i++;
  }
  std::cout << "]";
}

void printHisto(const std::unordered_map<int,int>& m, const std::string& header) {
  std::cout << header << "\n";
  int mn = std::numeric_limits<int>::max();
  int mx = 0;

  for (auto kv : m) {
    mn = std::min(mn,kv.first);
    mx = std::max(mx,kv.first);
  }

  for (int i = mn; i <= mx; i++) {
    auto search = m.find(i);
    if (search == m.end()) {
      std::cout << i << "\t0\n";
    } else {
      std::cout << i << "\t" << search->second << "\n";
    }
  }
}

void InspectIndex(const KmerIndex& index, const ProgramOptions& opt) {
  std::string gfa = opt.gfa;
  std::string bed = opt.bedFile;

  static const char *dna = "ACGT";
  auto Dna = [](int i) {return dna[i & 0x03];};

  int k = index.k;
  std::cout << "[inspect] Index version number = " << index.INDEX_VERSION << std::endl;
  //std::cout << "#[inspect] k = " << index.k << std::endl;;
  //std::cout << "#[inspect] number of targets = " << index.num_trans << std::endl;

  std::cout << "[inspect] number of unitigs = " << index.dbg.size() << std::endl;
  std::cout << "[inspect] minimizer length = " << index.dbg.getG() << std::endl;

  std::pair<size_t,size_t> ec_info = index.getECInfo();

  std::cout << "[inspect] max EC size = " << ec_info.first << std::endl;
  std::cout << "[inspect] number of ECs discarded = " << ec_info.second << std::endl;


  // std::cout << "#[inspect] Number of k-mers in index = " << index.dbg.nbKmers() << std::endl;

  if (!gfa.empty()) {
    index.dbg.write(gfa);
  }

  // for every contig
  if (opt.inspect_thorough) {
    for (const auto& contig : index.dbg) {
      const Node* n = contig.getData();    
    
      auto seq = contig.referenceUnitigToString();
      KmerIterator it(seq.c_str());
      KmerIterator it_end;
      int i = 0;
      while(it != it_end) {
        const Kmer km = it->first;
        const Kmer tw = km.twin();   
        auto um = index.dbg.find(km);
        auto um_tw = index.dbg.find(tw);
        if (um.isEmpty || um_tw.isEmpty) {      
          std::cout << "contig id: " << n->id << "\n";
          std::cout << "Kmer " << km.toString() << " is empty\n";
          exit(1);
        }



        if (um.dist != i || !um.strand) {
          std::cout << "contig id: " << n->id << "\n";
          std::cout << "Kmer " << km.toString() << " is not in the correct position\n";
          std::cout << "Expected position: " << i << ", actual position: " << um.dist << "\n";
          std::cout << "seq: " << seq << "\n";
          exit(1);
        }
        if (um_tw.dist != i || um_tw.strand ) {
          std::cout << "contig id: " << n->id << "\n";
          std::cout << "Kmer " << km.toString() << " is not in the correct position\n";
          std::cout << "Expected position: " << i << ", actual position: " << um_tw.dist << "\n";
          std::cout << "seq: " << seq << "\n";
          exit(1);
        }

        ++it;
        ++i;
      }
    }
    
    

  }

  // TODO:
  // Implement bedfile output for Bifrost index
  /*
  if (!bed.empty()) {
    // export bed track with TCC information
    bool guessChromosomes = false;
    Transcriptome model;
    if (opt.genomebam) {
      if (!opt.chromFile.empty()) {
        model.loadChromosomes(opt.chromFile);
      } else {
        guessChromosomes = true;
      }
      model.parseGTF(opt.gtfFile, index, opt, guessChromosomes);
      //model.loadTranscriptome(index, in, opt);
    }

    std::ofstream out;
    out.open(bed);

    out << "track name=\"Kallisto \" gffTags=\"on\"\n";
    std::unordered_map<TranscriptAlignment, std::vector<int>> cmap;
    cmap.reserve(100);
    std::vector<std::unordered_map<int, std::vector<ECStruct>>> ec_chrom(index.ecmap.size());

    for (const auto& um : index.dbg) {
      cmap.clear();
      // structure for TRaln
      TranscriptAlignment tra;

      const Node* n = um.getData();
      int len = um.size();
      int cid = n->id;
      for (const auto& ct : c.transcripts) {
        // ct.trid, ct.pos, ct.sense
        model.translateTrPosition(ct.trid, ct.pos, len, ct.sense, tra);
        cmap[tra].push_back(ct.trid);
      }

      for (const auto& cp : cmap) {
        const auto& tra = cp.first;
        const auto& tlist = cp.second;
        if (tra.chr != -1) {
          ECStruct ecs;
          ecs.chr = tra.chr;
          ecs.ec = index.dbGraph.ecs[c.id];
          ecs.start = tra.chrpos;
          int pos = 0;
          for (uint32_t cig : tra.cigar) {
            int len = cig >> BAM_CIGAR_SHIFT;
            int type = cig & 0xF;
            if (type == BAM_CMATCH) {
              ecs.start_lens.push_back({pos, len});
              pos += len;
            } else if (type == BAM_CREF_SKIP) {
              pos += len;
            } else {
              assert(false);
            }
          }
          ecs.stop = tra.chrpos + pos;
          ecs.tlist = tlist;

          ec_chrom[ecs.ec][tra.chr].push_back(ecs);
        }
      }
    }

    int num_ecs = ec_chrom.size();
    for (int ec = 0; ec < num_ecs; ec++) {
      auto& chrmap = ec_chrom[ec];
      if (chrmap.empty()) {
        continue;
      }
      bool single_chrom = (chrmap.size() == 1);
      for (auto& x : chrmap) {
        int chr = x.first;
        auto& ecv = x.second;

        auto blocklist = merge_contigs(ecv);

        for (auto &ecs : blocklist) {
          int i;
          out << model.chr[chr].name << "\t" // chrom
              << ecs.start << "\t" << ecs.stop << "\t"; // start stop
          out << "Name=" << ec << ";Transcripts=";
          i = 0;
          for (const auto& t : index.ecmap[ec]) {
            bool ischr = (std::find(ecs.tlist.begin(), ecs.tlist.end(), t) != ecs.tlist.end());
            if (i++ > 0) {
              out << "%0A";
            }
            out << index.target_names_[t];
            if (!ischr) {
              int tchr = model.transcripts[t].chr;
              if (tchr != -1) {
                if ( tchr != chr) {
                  out << "%20(" << model.chr[tchr].name << ")";
                } else {
                  out << "%20(*)";
                }
              } else {
                out << "%20(\?\?\?)";
              }
            }
          }
          out  << ";\t" <<  (int)(1000.0 / index.ecmap[ec].size()) << "\t" // score
               << ".\t" // strand
               << ecs.start << "\t" << ecs.stop // thick part
               << "\t0\t"; // color name
          const auto& sp = ecs.start_lens;
          out << sp.size() << "\t";
          i = 0;
          for (const auto& x : sp) {
            if (i++ > 0) {
              out << ',';
            }
            out << x.second;
          }
          out << "\t";
          i = 0;
          for (const auto& x : sp) {
            if (i++ > 0) {
              out << ',';
            }
            out << x.first;
          }
          out << "\n";
        }
      }
    }
  }
  */


  auto verifySequencePositions = [&](const int tr, const std::string& sequence) {
    KmerIterator kit(sequence.c_str()), kit_end;
    int tpos = 0;
    int trlen = index.target_lens_[tr];
    // find how many As are in the end of the sequence
    int num_As = 0;
    for (int i = sequence.size()-1; i >= 0; i--) {
      if (sequence[i] == 'A') {
        num_As++;
      } else {
        break;
      }
    }
    trlen -= num_As;
    for (; kit != kit_end; ++kit) {
      Kmer km = kit->first;
      Kmer tw = km.twin();
      size_t pos = kit->second;
      if (pos + index.k - 1 >= trlen) {
        continue;
      }

      auto um = index.dbg.find(km);
      if (um.isEmpty) {
        std::cout << "Error: Kmer " << km.toString() << " from sequence " << index.target_names_[tr] 
                  << " at position " << pos << " not found in the index." << std::endl;
        exit(1);
      }
      auto um_tw = index.dbg.find(tw);
      if (um_tw.isEmpty) {
        std::cout << "Error: Kmer " << tw.toString() << " from sequence " << index.target_names_[tr] 
                  << " at position " << pos << " not found in the index." << std::endl;
        exit(1);
      }      
      auto indexPositions = index.findPositions(tr, km, um);
      bool positionMatched = false;
      for (const auto& indexPos : indexPositions) {
        // check if any of the index positions match the position in the sequence
        if (indexPos.first == pos+1 && indexPos.second) {
          positionMatched = true;
          break;
        }
      }

      if (!positionMatched) {
        std::cout << "Error: Mismatch for kmer " << km.toString() << " from sequence " << index.target_names_[tr] 
                  << " (tr: " << tr << ") at position " << pos << ". Index position not found in sequence." << std::endl;
        
        if (!indexPositions.empty()) {
          std::cout << "Index positions: ";
          for (const auto& indexPos : indexPositions) {
            std::cout << indexPos.first << " " << indexPos.second << " ";
            std::cout << " um.dist: " << um.dist << std::endl;
          }
          std::cout << std::endl;
        } else {
          std::cout << "No index positions found for kmer " << km.toString() << " from sequence " << index.target_names_[tr] 
                    << " at position " << pos << "." << std::endl;
        }
        // continue;
        exit(1);
      }
      indexPositions.clear();
      indexPositions = index.findPositions(tr, tw, um_tw);
      positionMatched = false;
      for (const auto& indexPos : indexPositions) {
        if (indexPos.first == pos+k && !indexPos.second) {
          positionMatched = true;
          break;
        }
      }
      if (!positionMatched) {
        std::cout << "Error: Mismatch for twin of kmer " << tw.toString() << " from sequence " << index.target_names_[tr] 
                  << " at position " << pos << ". Index position not found in sequence." << std::endl;
        if (!indexPositions.empty()) {
          std::cout << "Index positions: ";
          for (const auto& indexPos : indexPositions) {
            std::cout << indexPos.first << " " << indexPos.second << " ";
          }
          std::cout << std::endl;
        }
        continue;
      }
    }
  };

  
  if (opt.inspect_thorough) {
    std::string fastaFile = opt.transcriptsFile;  
    if (!fastaFile.empty()) {
      std::ifstream fasta(fastaFile);
      if (!fasta.is_open()) {
        std::cerr << "Error: Unable to open FASTA file " << fastaFile << std::endl;
        return;
      }

      std::string line, header, sequence;
      // measure the time it takes to load the transcript sequences
      
      index.loadTranscriptSequences();
      
      int tr = 0;
      while (std::getline(fasta, line)) {
        if (line[0] == '>') {
          if (!sequence.empty()) {
            if (sequence.size() >= index.k) {
              
              // Process the previous sequence
              verifySequencePositions(tr, sequence);
              if (sequence != index.target_seqs_[tr]) {
                // count polyA tail
                int num_As = 0;
                for (int i = sequence.size()-1; i >= 0; i--) {
                  if (sequence[i] == 'A') {
                    num_As++;
                  } else {
                    break;
                  }
                }
               
                // iterate over the characters and match unless the sequences has an N
                bool mismatch = false;
                for (int i = 0; i < sequence.size()-num_As; i++) {
                  if (sequence[i] != index.target_seqs_[tr][i] && sequence[i] != 'N') {
                    mismatch = true;
                    break;
                  }
                }
                if (mismatch) {
                  std::cout << "Error: Sequence mismatch for transcript " << index.target_names_[tr] << std::endl;
                  std::cout << ">expected \n" << sequence << "\n"
                            << ">found \n" << index.target_seqs_[tr] << std::endl;
                  exit(1);
                }
                
              }
            }
            sequence.clear();
            tr++;
          }
          header = line.substr(1);
        } else {
          sequence += line;
        }
      }
      // Process the last sequence
      if (!sequence.empty()) {
        verifySequencePositions(tr, sequence);
        if (sequence != index.target_seqs_[tr]) {
          std::cout << "Error: Sequence mismatch for transcript " << index.target_names_[tr] << std::endl;
          exit(1);
        }
      }
    }
  }
}



#endif // KALLISTO_INSPECTINDEX_H
