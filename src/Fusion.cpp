#include <sstream>
#include <set>

#include "ProcessReads.h"
#include "Fusion.h"


typedef std::vector<std::pair<const_UnitigMap<Node>, int32_t> > MappedVector;

/** -- fusion functions -- **/

/**
FUSION output is written to OUT_DIR/fusion.txt where 'OUT_DIR' is the specified
output directory for kallisto. fusion.txt is a tab separated text file with 10 columns
TYPE  NAME1  SEQ1  KPOS1  NAME2  SEQ2  KPOS2  INFO  POS1  POS2

TYPE can be one of SPLIT or PAIR, NAME1 and NAME2 are the sequence identifiers. 
SEQ1 and SEQ2 are the raw read sequences identifie as overlapping the potential fusion.
KPOS1 and KPOS2 are the positions within the reads where the k-mers match.

PAIRs are readpairs where fragment doesn't map, but each read
maps independently, potentially the fusion breakpoint is within the fragment, but not sequenced.

SPLITs are reads where the fusion breakpoint appears within one of the reads

INFO is NA for PAIR but of the form splitat=(a,b)
where a=0,1 indicating whether the first or the second pair is split and b is the 0-based position of the
split within the read.

POS1 and POS2 are comma separated lists of locations of SEQ1 and SEQ2 respectively. Each position is of
the form (tx_name,pos,strand) where pos is 0-based from the start of the transcript tx_name and strand
is either FW or RE depending on whether the read aligns to the forward or reverse of the transcript.
**/

void printTranscripts(const KmerIndex& index, std::stringstream& o, const std::string s,
  const MappedVector& v, const Roaring& u) {

  
  Kmer km;
  KmerEntry val;
  int p;
  const_UnitigMap<Node> um;

  // find first mapping k-mer
  if (!v.empty()) {
    Roaring vtmp;
    auto res = findFirstMappingKmer(v);
    um = res.first;
    p = res.second;
    km = um.getMappedHead();
    //km = Kmer((s.c_str()+p));
  }
  

  int i = -1;
  for (const auto &tr : u) {
    i++;
  /*for (int i = 0; i < u.size(); i++) {
    int tr = u[i];*/
    if (i > 0) {
      o << ";";
    }

    std::pair<int, bool> xp = index.findPosition(tr, km, p);
    o << "(" << index.target_names_[tr] << "," << xp.first << ",";
    if (xp.second) {
      o << "FW)";
    } else {
      o << "RC)";
    }
  }
  
}

Roaring simpleIntersect(const KmerIndex& index, const MappedVector& v) {
  
  if (v.empty()) {
    return {};
  }

  // get the ec from the debruijn graph
  const auto um = v[0].first;
  
  Roaring ec = um.getData()->ec[um.dist].getIndices(); //= index.dbg.
  //int ec = index.dbGraph.ecs[v[0].first.contig];
  Roaring lastEC = ec;
  // find the roaring vector from the ecmap
  Roaring u = ec;

  for (int i = 1; i < v.size(); i++) {
    // anytime the contig changes, update the ec
    if (!v[i].first.isSameReferenceUnitig(v[i-1].first)) {
      const auto um = v[i].first;
      ec = um.getData()->ec[um.dist].getIndices();
      if (!(ec == lastEC)) {
        u &= ec;
        lastEC = ec;
        if (u.cardinality() == 0) {
          return u;
        }
      }
    }
  }
  return u;  
}




// returns true if the intersection of the union of EC classes for v1 and v2 respectively is empty
bool checkUnionIntersection(const KmerIndex& index, const std::string &s1, const std::string &s2, std::pair<int,int> &p1, std::pair<int,int> &p2) {
//const std::vector<std::pair<KmerEntry,int>>& v1, const std::vector<std::pair<KmerEntry,int>>& v2) {
  
  Roaring su1,su2;
  
  auto union_set = [&](const std::string &s, std::pair<int,int> &p) {
    p = {-1,-1};
    Roaring su;
    
    KmerIterator kit(s.c_str()), kit_end;
    Roaring lastEC;
    for (int i = 0; kit != kit_end; ++i,++kit) {
      auto um = index.dbg.find(kit->first);
      //auto search = index.kmap.find(kit->first.rep());
      //if (search != index.kmap.end()) {
      if (!um.isEmpty) {
        if (p.first == -1) {
          p.first = kit->second;
          p.second = p.first +1;
        } else {
          if (p.second + 1 == kit->second) {
            p.second++;
          }
        }

        //const KmerEntry &val = search->second;
        ///int ec = index.dbGraph.ecs[val.contig];
        Roaring ec = um.getData()->ec[um.dist].getIndices();

        if (!(ec.cardinality() == 0) && !(ec == lastEC)) {
          su |= ec;
          lastEC = ec;
        }
      }
    }    
    return su;
  };

  su1 = union_set(s1,p1);
  su2 = union_set(s2,p2);
  
  if ((su1.cardinality()==0)  || (su2.cardinality()==0)) {
    return false; // TODO, decide on this
  }
  
  Roaring su3 = su1 & su2;
  return su3.cardinality() == 0;
}


void searchFusion(const KmerIndex &index, const ProgramOptions& opt,
  const MinCollector& tc, MasterProcessor& mp, 
  const std::string &n1, const std::string &s1, 
  MappedVector &v1,
  const std::string &n2, const std::string &s2, 
  MappedVector &v2, 
  bool paired) {

  
  std::stringstream o;

  // no mapping information
  if (v1.empty() && v2.empty()) {
    return; // consider splitting in case either is empty
  }

  auto u1 = simpleIntersect(index, v1);
  auto u2 = simpleIntersect(index, v2);

  
  // discordant pairs
  if (!v1.empty() && !v2.empty()) {
    if (!(u1.cardinality()==0) && !(u2.cardinality() == 0)) {
      std::pair<int,int> p1,p2;
      if (checkUnionIntersection(index, s1, s2, p1, p2)) {
        // each pair maps to either end
        // each pair maps to either end
        o << "PAIR\t";
        o << n1 << "\t";
        o << s1 << "\t";
        o << p1.first << "," << p1.second << "\t"; 
        if (paired) {
          o << n2 << "\t";
          o << s2 << "\t";
          o << p2.first << "," << p2.second << "\t\tNA\t";
        } else {
          o << "\t\t\tNA\t";
        }
        printTranscripts(index, o, s1, v1, u1); o << "\t"; printTranscripts(index, o, s2, v2, u2);
        mp.outputFusion(o);
        return;
      }
    }
  }

  if ((v1.empty() && (u1.cardinality() != 0)) || (v2.empty() && (u2.cardinality() != 0))) {
    // read was trimmed    
    return;
  }



  // ok so ec == -1 and not both v1 and v2 are empty
  // exactly one of u1 and u2 are empty
  MappedVector vsafe, vsplit;
  // sort v1 and v2 by read position
  auto vsorter =  [&](std::pair<const_UnitigMap<Node>, int> a, std::pair<const_UnitigMap<Node>, int> b) {
    return a.second < b.second;
  };

  std::sort(v1.begin(), v1.end(), vsorter);
  std::sort(v2.begin(), v2.end(), vsorter);
  bool split1 = true;
  if (!v1.empty() && !v2.empty()) {
    if (u1.cardinality() == 0) {
      std::copy(v1.begin(), v1.end(), std::back_inserter(vsplit));
      std::copy(v2.begin(), v2.end(), std::back_inserter(vsafe));
    } else {
      std::copy(v2.begin(), v2.end(), std::back_inserter(vsplit));
      std::copy(v1.begin(), v1.end(), std::back_inserter(vsafe));
      split1 = false;
    }
  } else if (v1.empty()) {
    if (u2.cardinality() == 0) {
      std::copy(v2.begin(), v2.end(), std::back_inserter(vsplit));
    } else {
      assert(false);
      return; // can this happen?
    }
  } else if (v2.empty()) {
    if (u1.cardinality() == 0) {
      std::copy(v1.begin(), v1.end(), std::back_inserter(vsplit));
    } else {
      assert(false);
      return;
    }
  }


  // now we look for a split j s.t. vsplit[0:j] and vsplit[j:] + vsafe both have nonempty ec
  int j = vsplit.size()-1;
  while (!vsplit.empty()) {
    auto x = vsplit.back();
    vsplit.pop_back();
    vsafe.push_back(x);

    auto ut1 = simpleIntersect(index,vsplit);
    auto ut2 = simpleIntersect(index,vsafe);

    if (!(ut1.cardinality() == 0) && !(ut2.cardinality() == 0)) {
      std::pair<int,int> p1,p2;
      std::string st1,st2;
      if (u1.cardinality() == 0) {
        st1 = s1.substr(0,vsafe.back().second);
        st2 = s2;
      } else {
        st1 = s1;
        st2 = s2.substr(0,vsafe.back().second);
      }
      if (checkUnionIntersection(index, st1, st2, p1, p2)) { // need to check this more carefully
        o << "SPLIT\t";
        o << n1 << "\t" << s1 << "\t" << p1.first << "," << p1.second << "\t";
        // what to put as info?
        if (paired) {
          o << n2 << "\t" << s2 << "\t" << p2.first << "," << p2.second << "\t";
        } else {
          o << "\t\t\t";
        }
        o << "splitat=";
        if (u1.cardinality() == 0) {
          o << "(0,";
        } else {
          o << "(1,";
        }
        o << vsafe.back().second << ")\t";
        // fix this
        if (split1) {
          printTranscripts(index, o, s1, vsplit, ut1); o << "\t"; printTranscripts(index, o, s2, vsafe, ut2);
        } else {
          printTranscripts(index, o, s1, vsafe, ut2); o << "\t"; printTranscripts(index, o, s2, vsplit, ut1);
        }
        mp.outputFusion(o);
        return;
      }
    } else {
      j--;
    }
  }

  return;
  
}
