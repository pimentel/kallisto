#include "common.h"
#include "PseudoBam.h"

#ifndef NO_HTSLIB


bam_hdr_t* createPseudoBamHeaderTrans(const KmerIndex& index)  {
  bam_hdr_t *h = bam_hdr_init();
  h->n_targets = index.num_trans;
  //todo include program parameters in string
  std::string text = "@HD\tVN:1.0\n@PG\tID:kallisto\tPN:kallisto\tVN:";
  text += KALLISTO_VERSION;
  text += "\n";
  h->text = strdup(text.c_str());
  h->l_text = (uint32_t) strlen(h->text);
  h->target_len = (uint32_t *) calloc(index.num_trans, sizeof(uint32_t));
  h->target_name = (char**) calloc(index.num_trans, sizeof(char*));
  for (int i = 0; i < index.num_trans; i++) {
    h->target_len[i] = (uint32_t) index.target_lens_[i];
    h->target_name[i] = strdup(index.target_names_[i].c_str());
  }
  return h;
}

bam_hdr_t* createPseudoBamHeaderGenome(const Transcriptome& model)  {
  bam_hdr_t *h = bam_hdr_init();

  std::string text = "@HD\tVN:1.0\n@PG\tID:kallisto\tPN:kallisto\tVN:";
  text += KALLISTO_VERSION;
  text += "\n";
  int num_chr = model.chr.size();
  h->n_targets = num_chr;
  h->target_len = (uint32_t *) calloc(num_chr, sizeof(uint32_t));
  h->target_name = (char**) calloc(num_chr, sizeof(char*));
  for (int i = 0; i < num_chr; i++) {
    h->target_len[i] = (uint32_t) model.chr[i].len;
    h->target_name[i] = strdup(model.chr[i].name.c_str());
    text += "@SQ\tSN:" + model.chr[i].name + "\tLN:" + std::to_string(model.chr[i].len) + "\n";
  }
  h->text = strdup(text.c_str());
  h->l_text = (uint32_t) strlen(h->text);


  return h;
}



/** --- pseudobam functions -- **/
/*
void createBamRecord(const KmerIndex, const KmerIndex &index, const std::vector<int> &u,
    const char *s1, const char *n1, const char *q1, int slen1, int nlen1, const std::vector<std::pair<KmerEntry,int>>& v1,
    const char *s2, const char *n2, const char *q2, int slen2, int nlen2, const std::vector<std::pair<KmerEntry,int>>& v2,
    bool paired );*/



void outputPseudoBam(const KmerIndex &index, const Roaring& u,
    const char *s1, const char *n1, const char *q1, int slen1, int nlen1, const std::vector<std::pair<const_UnitigMap<Node>&, int>>& v1,
    const char *s2, const char *n2, const char *q2, int slen2, int nlen2, const std::vector<std::pair<const_UnitigMap<Node>&, int>>& v2,
    bool paired, bam_hdr_t *h, samFile *fp) {

  static char buf1[32768];
  static char buf2[32768];
  static char cig_[1000];
  char *cig = &cig_[0];
  int flag1=0,flag2=0;
  bam1_t b1,b2;

  if (nlen1 > 2 && n1[nlen1-2] == '/') {
    ((char*)n1)[nlen1-2] = 0;
    nlen1 -=2;
  }

  if (paired && nlen2 > 2 && n2[nlen2-2] == '/') {
    ((char*)n2)[nlen2-2] = 0;
    nlen2 -= 2;
  }

  /*
  b1.l_data = fillBamRecord(b1, &buf1[0], s1,n1,q1,slen1,nlen1);
  b1.data = (uint8_t*)buf1;

  if (paired) {
    b2.l_data = fillBamRecord(b2, &buf2[0], s2,n2,q2,slen2,nlen2);
    b2.data = (uint8_t*) buf2;
  }
  */
  if (u.isEmpty()) {
    // no mapping
    if (paired) {
      b1.core.tid = -1;
      b1.core.pos = -1;
      b1.core.bin = 4680; // magic bin for unmapped reads
      b1.core.qual = 0;
      b1.core.flag = BAM_FPAIRED | BAM_FREAD1 | BAM_FUNMAP | BAM_FMUNMAP;
      b1.core.mtid = -1;
      b1.core.mpos = -1;
      b1.core.isize = 0;

      b2.core.tid = -1;
      b2.core.pos = -1;
      b2.core.bin = 4680; // magic bin for unmapped reads
      b2.core.qual = 0;
      b2.core.flag = BAM_FPAIRED | BAM_FREAD2 | BAM_FUNMAP | BAM_FMUNMAP;
      b2.core.mtid = -1;
      b2.core.mpos = -1;
      b2.core.isize = 0;

      //(fp, h, &b1);
      //sam_write1(fp, h, &b2);
      //printf("%s\t77\t*\t0\t0\t*\t*\t0\t0\t%s\t%s\n", n1,s1,q1);
      //printf("%s\t141\t*\t0\t0\t*\t*\t0\t0\t%s\t%s\n", n2,s2,q2);
      //o << seq1->name.s << "" << seq1->seq.s << "\t" << seq1->qual.s << "\n";
      //o << seq2->name.s << "\t141\t*\t0\t0\t*\t*\t0\t0\t" << seq2->seq.s << "\t" << seq2->qual.s << "\n";
    } else {
      b1.core.tid = -1;
      b1.core.pos = -1;
      b1.core.bin = 4680; // magic bin for unmapped reads
      b1.core.qual = 0;
      b1.core.flag = BAM_FUNMAP ;
      b1.core.mtid = -1;
      b1.core.mpos = -1;
      //printf("%s\t4\t*\t0\t0\t*\t*\t0\t0\t%s\t%s\n", n1,s1,q1);
      //sam_write1(fp, h, &b1);
    }
  } else {
    if (paired) {
      b1.core.flag = BAM_FPAIRED | BAM_FREAD1;
      b2.core.flag = BAM_FPAIRED | BAM_FREAD2;
      int flag1 = 0x01 + 0x40;
      int flag2 = 0x01 + 0x80;

      if (v1.empty()) {
        b1.core.flag |= BAM_FUNMAP;
        b2.core.flag |= BAM_FMUNMAP;
        //flag1 += 0x04; // read unmapped
        //flag2 += 0x08; // mate unmapped
      }

      if (v2.empty()) {
        b1.core.flag |= BAM_FMUNMAP;
        b2.core.flag |= BAM_FUNMAP;
        //flag1 += 0x08; // mate unmapped
        //flag2 += 0x04; // read unmapped
      }

      if (!v1.empty() && !v2.empty()) {
        b1.core.flag |= BAM_FPROPER_PAIR;
        b2.core.flag |= BAM_FPROPER_PAIR;
        //flag1 += 0x02; // proper pair
        //flag2 += 0x02; // proper pair
      }


      int p1 = -1, p2 = -1;
      const_UnitigMap<Node> um1, um2;
      int nmap = u.cardinality();
      Kmer km1, km2;

      if (!v1.empty()) {
        um1 = v1[0].first;
        p1 = v1[0].second;
        for (auto &x : v1) {
          if (x.second < p1) {
            um1 = x.first;
            p1 = x.second;
          }
        }
        km1 = Kmer((s1+p1));
      }

      if (!v2.empty()) {
        um2 = v2[0].first;
        p2 = v2[0].second;
        for (auto &x : v2) {
          if (x.second < p2) {
            um2 = x.first;
            p2 = x.second;
          }
        }
        km2 = Kmer((s2+p2));
      }

      bool revset = false;

      // output pseudoalignments for read 1
      bool firstTr = true;

      for (auto tr : u) {
        int f1 = flag1;
        std::pair<int, bool> x1 {-1,true};
        std::pair<int, bool> x2 {-1,true};
        if (p1 != -1) {
          x1 = index.findPosition(tr, km1, p1);
          if (p2 == -1) {
            x2 = {x1.first,!x1.second};
          }
          if (!x1.second) {
            f1 += 0x10; // read reverse
            if (!revset) {
              revseq(&buf1[0], &buf2[0], s1, q1, slen1);
              revset = true;
            }
          }
        }
        if (!firstTr) {
          f1 += 0x100; // secondary alignment
        }
        if (p2 != -1) {
          x2 = index.findPosition(tr, km2 ,  p2);
          if (p1 == -1) {
            x1 = {x2.first, !x2.second};
          }
          if (!x2.second) {
            f1 += 0x20; // mate reverse
          }
        }
        firstTr = false;

        int posread = (f1 & 0x10) ? (x1.first - slen1 + 1) : x1.first;
        int posmate = (f1 & 0x20) ? (x2.first - slen2 + 1) : x2.first;
        if (v1.empty()) {
          posread = posmate;
        }
        if (v2.empty()) {
          posmate = posread;
        }

        getCIGARandSoftClip(cig, bool(f1 & 0x10), (f1 & 0x04) == 0, posread, posmate, slen1, index.target_lens_[tr]);
        int tlen = x2.first - x1.first;
        if (tlen != 0) {
          tlen += (tlen>0) ? 1 : -1;
        }

        //printf("%s\t%d\t%s\t%d\t%d\t%s\t=\t%d\t%d\t%s\t%s\tNH:i:%d\n", n1, f1 & 0xFFFF, index.target_names_[tr].c_str(), posread, (!v1.empty()) ? 255 : 0 , cig, posmate, tlen, (f1 & 0x10) ? &buf1[0] : s1, (f1 & 0x10) ? &buf2[0] : q1, nmap);
        if (v1.empty()) {
          break; // only report primary alignment
        }
      }

      revset = false;
      // output pseudoalignments for read 2
      firstTr = true;
      for (auto tr : u) {
        int f2 = flag2;
        std::pair<int, bool> x1 {-1,true};
        std::pair<int, bool> x2 {-1,true};
        if (p1 != -1) {
          x1 = index.findPosition(tr, km1, p1);
          if (p2 == -1) {
            x2 = {x1.first,!x1.second};
          }
          if (!x1.second) {
            f2 += 0x20; // mate reverse
          }
        }
        if (p2 != -1) {
          x2 = index.findPosition(tr, km2, p2);
          if (p1 == -1) {
            x1 = {x2.first, !x2.second};
          }
          if (!x2.second) {
            f2 += 0x10; // read reverse
            if (!revset) {
              revseq(&buf1[0], &buf2[0], s2, q2, slen2);
              revset = true;
            }

          }
        }
        if (!firstTr) {
          f2 += 0x100; // secondary alignment
        }

        firstTr = false;
        int posread = (f2 & 0x10) ? (x2.first - slen2 + 1) : x2.first;
        int posmate = (f2 & 0x20) ? (x1.first - slen1 + 1) : x1.first;
        if (v1.empty()) {
          posmate = posread;
        }
        if (v2.empty()) {
          posread = posmate;
        }

        getCIGARandSoftClip(cig, bool(f2 & 0x10), (f2 & 0x04) == 0, posread, posmate, slen2, index.target_lens_[tr]);
        int tlen = x1.first - x2.first;
        if (tlen != 0) {
          tlen += (tlen > 0) ? 1 : -1;
        }

        //printf("%s\t%d\t%s\t%d\t%d\t%s\t=\t%d\t%d\t%s\t%s\tNH:i:%d\n", n2, f2 & 0xFFFF, index.target_names_[tr].c_str(), posread, (!v2.empty()) ? 255 : 0, cig, posmate, tlen, (f2 & 0x10) ? &buf1[0] : s2,  (f2 & 0x10) ? &buf2[0] : q2, nmap);
        if(v2.empty()) {
          break; // only print primary alignment
        }
      }


    } else {
      // single end
      int nmap = (int) u.cardinality();
      const_UnitigMap<Node> um = v1[0].first;
      int p1 = v1[0].second;
      for (auto &x : v1) {
        if (x.second < p1) {
          um = x.first;
          p1 = x.second;
        }
      }
      Kmer km1 = Kmer((s1+p1));

      bool revset = false;
      bool firstTr = true;
      for (auto tr : u) {
        int f1 = 0;
        auto x1 = index.findPosition(tr, km1, p1);

        if (!x1.second) {
          f1 += 0x10;
          if (!revset) {
            revseq(&buf1[0], &buf2[0], s1, q1, slen1);
            revset = true;
          }
        }
        if (!firstTr) {
          f1 += 0x100; // secondary alignment
        }
        firstTr = false;
        int posread = (f1 & 0x10) ? (x1.first - slen1+1) : x1.first;
        int dummy=1;
        getCIGARandSoftClip(cig, bool(f1 & 0x10), (f1 & 0x04) == 0, posread, dummy, slen1, index.target_lens_[tr]);

        printf("%s\t%d\t%s\t%d\t255\t%s\t*\t%d\t%d\t%s\t%s\tNH:i:%d\n", n1, f1 & 0xFFFF, index.target_names_[tr].c_str(), posread, cig, 0, 0, (f1 & 0x10) ? &buf1[0] : s1, (f1 & 0x10) ? &buf2[0] : q1, nmap);
      }
    }
  }
}


void revseq(char *b1, char *b2, const char *s, const char *q, int n) {
  b1[n] = 0;
  for (int i = 0; i < n; i++) {
    switch(s[i]) {
    case 'A': b1[n-1-i] = 'T'; break;
    case 'C': b1[n-1-i] = 'G'; break;
    case 'G': b1[n-1-i] = 'C'; break;
    case 'T': b1[n-1-i] = 'A'; break;
    default:  b1[n-1-i] = 'N';
    }
  }
  b2[n] = 0;
  for (int i = 0; i < n; i++) {
    b2[n-1-i] = q[i];
  }
}



void getCIGARandSoftClip(char* cig, bool strand, bool mapped, int &posread, int &posmate, int length, int targetlength) {
  int softclip = 1 - posread;
  int overhang = (posread + length) - targetlength - 1;

  if (posread <= 0) {
    posread = 1;
  }

  if (mapped) {
    if (softclip > 0) {
      if (overhang > 0) {
        sprintf(cig, "%dS%dM%dS",softclip, (length-overhang - softclip), overhang);
      } else {
        sprintf(cig, "%dS%dM",softclip,length-softclip);
      }
    } else if (overhang > 0) {
      sprintf(cig, "%dM%dS", length-overhang, overhang);
    } else {
      sprintf(cig, "%dM",length);
    }
  } else {
    sprintf(cig, "*");
  }


  if (posmate <= 0) {
    posmate = 1;
  }
}

#endif // NO_HTSLIB

/** -- pseudoalignment info methods -- **/

void writePseudoAlignmentBatch(std::ofstream& of, const PseudoAlignmentBatch& batch) {

  of.write("BATCH=",6);
  of.write((char*)&(batch.batch_id), sizeof(int32_t));
  uint32_t bsz = batch.aln.size();
  of.write((char*)&(bsz), sizeof(uint32_t));

  // Initialize a vector with 1KB capacity
  std::vector<char> buffer;
  buffer.reserve(1024);

  for(const auto &x : batch.aln) {
    of.write((char*)&x.id, sizeof(x.id));
    uint8_t flag = 0;
    flag |= (x.paired) ? 1 : 0;
    flag |= (x.r1empty) ? 2 : 0;
    flag |= (x.r2empty) ? 4 : 0;
    of.write((char*)&flag,1);
    uint8_t k1 = (0 <= x.k1pos && x.k1pos < 255) ? x.k1pos : 255;
    uint8_t k2 = (0 <= x.k2pos && x.k2pos < 255) ? x.k2pos : 255;
    of.write((char*)&k1, 1);
    of.write((char*)&k2, 1);

    // Clear the buffer and ensure it has enough capacity
    buffer.clear();
    buffer.resize(x.ec.getSizeInBytes());

    // Write the Roaring structure to the buffer
    size_t roaring_size = x.ec.write(buffer.data());

    // Resize the buffer to the actual size used
    buffer.resize(roaring_size);

    // Write the size and content of the Roaring structure
    of.write((char*)&roaring_size, sizeof(roaring_size));
    of.write(buffer.data(), roaring_size);

    of.write((char*)&x.barcode, sizeof(uint64_t));
    of.write((char*)&x.UMI, sizeof(uint64_t));
    of.put(0); // mark the end of record
  }
}


void readPseudoAlignmentBatch(std::ifstream& in, PseudoAlignmentBatch& batch) {
  batch.aln.clear();
  char bb[7];
  char mark0;
  in.read(&bb[0], 6);
  bb[6] = 0;
  assert(strcmp(bb,"BATCH=")==0);
  in.read((char*)&(batch.batch_id), sizeof(int32_t));
  uint32_t bsz;
  in.read((char*)&bsz, sizeof(uint32_t));
  batch.aln.reserve(bsz);

  // Initialize a vector with 1KB capacity
  std::vector<char> buffer;
  buffer.reserve(1024);

  for (int i = 0; i < bsz; i++) {
    PseudoAlignmentInfo info;
    in.read((char*)&info.id, sizeof(info.id));
    uint8_t flag;
    in.read((char*)&flag, 1);
    info.paired  = (flag & 1) != 0;
    info.r1empty = (flag & 2) != 0;
    info.r2empty = (flag & 4) != 0;
    uint8_t k1,k2;
    in.read((char*)&k1,1);
    in.read((char*)&k2,1);
    info.k1pos = (k1 == 255) ? -1 : k1;
    info.k2pos = (k2 == 255) ? -1 : k2;
    size_t roaring_size;
    in.read((char*)&roaring_size, sizeof(roaring_size));

    
    buffer.resize(roaring_size);

    in.read(buffer.data(), roaring_size);
    info.ec = info.ec.read(buffer.data());

    in.read((char*)&info.barcode, sizeof(uint64_t));
    in.read((char*)&info.UMI, sizeof(uint64_t));

    mark0 = in.get();
    assert(mark0 == '\0');
    batch.aln.push_back(std::move(info));
  }
}
