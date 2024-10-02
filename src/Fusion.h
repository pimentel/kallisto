#ifndef KALLISTO_FUSION_H
#define KALLISTO_FUSION_H

#include <sstream>
#include <set>




typedef std::vector<std::pair<const_UnitigMap<Node>, int32_t> > MappedVector;

/** -- fusion functions -- **/


void printTranscripts(const KmerIndex& index, std::stringstream& o, const std::string s,
  const MappedVector& v, const Roaring& u);
Roaring simpleIntersect(const KmerIndex& index, const MappedVector& v);
bool checkMapability(const KmerIndex& index, const std::string &s, const std::vector<std::pair<KmerEntry,int>>& v, std::vector<int> &u);
bool checkUnionIntersection(const KmerIndex& index, const std::string &s1, const std::string &s2, std::pair<int,int> &p1, std::pair<int,int> &p2);
void searchFusion(const KmerIndex &index, const ProgramOptions& opt,
  const MinCollector& tc, MasterProcessor& mp, 
  const std::string &n1, const std::string &s1, 
  MappedVector &v1,
  const std::string &n2, const std::string &s2, 
  MappedVector &v2, 
  bool paired);


#endif // KALLISTO_FUSION_H