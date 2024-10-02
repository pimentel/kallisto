#ifdef USE_HDF5
#include "H5Writer.h"

void H5Writer::init(const std::string& fname, int num_bootstrap, int num_processed,
  const std::vector<int>& fld,const std::vector<int>& preBias, const std::vector<double>& postBias,
  uint compression, size_t index_version,
  const std::string& shell_call, const std::string& start_time)
{
  primed_ = true;
  num_bootstrap_ = num_bootstrap;
  compression_ = compression;
  file_id_ = H5Fcreate(fname.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  root_ = H5Gopen(file_id_, "/", H5P_DEFAULT);
  aux_ = H5Gcreate(file_id_, "/aux", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (num_bootstrap_ > 0) {
    bs_ = H5Gcreate(file_id_, "/bootstrap", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  }
  std::vector<int> n_bs {num_bootstrap};
  vector_to_h5(n_bs, aux_, "num_bootstrap", false, compression_);

  std::vector<int> n_proc {num_processed};
  vector_to_h5(n_proc, aux_, "num_processed", false, compression_);

  vector_to_h5(fld, aux_, "fld", false, compression_);

  vector_to_h5(preBias, aux_, "bias_observed", false, compression_);

  vector_to_h5(postBias, aux_, "bias_normalized", false, compression_);

  // info about run
  std::vector<std::string> kal_version{ KALLISTO_VERSION };
  vector_to_h5(kal_version, aux_, "kallisto_version", true, compression_);

  std::vector<int> idx_version{ static_cast<int>(index_version) };
  vector_to_h5(idx_version, aux_, "index_version", false, compression_);

  std::vector<std::string> call{ shell_call };
  vector_to_h5(call, aux_, "call", true, compression_);

  std::vector<std::string> s_time{ start_time };
  vector_to_h5(s_time, aux_, "start_time", true, compression_);
}

H5Writer::~H5Writer() {
  if (!primed_) {
    return;
  }
  if (num_bootstrap_ > 0) {
    H5Gclose(bs_);
  }
  H5Gclose(aux_);
  H5Gclose(root_);
  H5Fclose(file_id_);
}

void H5Writer::write_main(const EMAlgorithm& em,
    const std::vector<std::string>& targ_ids,
    const std::vector<int>& lengths) {
  vector_to_h5(em.alpha_, root_, "est_counts", false, compression_);

  vector_to_h5(targ_ids, aux_, "ids", true, compression_);
  vector_to_h5(em.eff_lens_, aux_, "eff_lengths", false, compression_);
  vector_to_h5(lengths, aux_, "lengths", false, compression_);
}

void H5Writer::write_bootstrap(const EMAlgorithm& em, int bs_id) {
  std::string bs_id_str("bs" + std::to_string( bs_id ));
  vector_to_h5(em.alpha_, bs_, bs_id_str.c_str(), false, compression_);
}

/**********************************************************************/

H5Converter::H5Converter(const std::string& h5_fname, const std::string& out_dir) :
  out_dir_(out_dir)
{
  file_id_ = H5Fopen(h5_fname.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  root_ = H5Gopen(file_id_, "/", H5P_DEFAULT);
  aux_ = H5Gopen(file_id_, "/aux", H5P_DEFAULT);


  // <aux info>
  // read target ids
  read_dataset(aux_, "ids", targ_ids_);
  std::cerr << "[h5dump] number of targets: " << targ_ids_.size() <<
    std::endl;

  n_targs_ = targ_ids_.size();

  read_dataset(aux_, "lengths", lengths_);
  assert( n_targs_ == lengths_.size() );

  read_dataset(aux_, "eff_lengths", eff_lengths_);
  assert( n_targs_ == eff_lengths_.size() );

  // read bootstrap info
  std::vector<int> n_bs_vec;
  read_dataset(aux_, "num_bootstrap", n_bs_vec);
  n_bs_ = n_bs_vec[0];

  // read bootstrap info
  std::vector<int> n_proc_vec;
  read_dataset(aux_, "num_processed", n_proc_vec);
  n_proc_ = n_proc_vec[0];


  std::cerr << "[h5dump] number of bootstraps: " << n_bs_ << std::endl;
  // </aux info>
  if (n_bs_ > 0) {
    bs_ = H5Gopen(file_id_, "/bootstrap", H5P_DEFAULT);
  }

  std::vector<std::string> tmp;
  read_dataset(aux_, "kallisto_version", tmp);
  kal_version_ = tmp[0];
  tmp.clear();
  std::cerr << "[h5dump] kallisto version: " << kal_version_ << std::endl;

  std::vector<int> idx_version;
  read_dataset(aux_, "index_version", idx_version);
  idx_version_ = static_cast<size_t>(idx_version[0]);
  std::cerr << "[h5dump] index version: " << idx_version_ << std::endl;

  read_dataset(aux_, "start_time", tmp);
  start_time_ = tmp[0];
  tmp.clear();
  std::cerr << "[h5dump] start time: " << start_time_ << std::endl;

  read_dataset(aux_, "call", tmp);
  call_ = tmp[0];
  tmp.clear();
  std::cerr << "[h5dump] shell call: " << call_ << std::endl;

  alpha_buf_.resize( n_targs_, 0.0 );
  assert( n_targs_ == alpha_buf_.size() );

  tpm_buf_.resize( n_targs_, 0.0 );
  assert( n_targs_ == tpm_buf_.size() );
}

H5Converter::~H5Converter() {
  if (n_bs_ > 0) {
    H5Gclose(bs_);
  }
  H5Gclose(aux_);
  H5Gclose(root_);
  H5Fclose(file_id_);
}

void H5Converter::write_aux() {
  std::string out_name(out_dir_ + "/run_info.json");

  std::vector<double> alpha;
  read_dataset(root_, "est_counts", alpha);
  double s = 0.0;
  for (auto &x : alpha) {
    s += x;
  }
  n_paln_ = (int) std::round(s);

  plaintext_aux(
      out_name,
      std::string(std::to_string(n_targs_)),
      std::string(std::to_string(n_bs_)),
      std::string(std::to_string(n_proc_)),
      std::string(std::to_string(n_paln_)),
      std::string(std::to_string(-1)),
      kal_version_,
      std::string(std::to_string(idx_version_)),
      std::string("dummy k-mer length"),
      start_time_,
      call_
      );
}

void H5Converter::convert() {

  std::cerr << "[h5dump] writing abundance file: " << out_dir_ << "/abundance.tsv" << std::endl;
  rw_from_counts(root_, "est_counts", out_dir_ + "/abundance.tsv");

  if (n_bs_ > 0) {
    std::cerr << "[h5dump] writing bootstrap abundance files: " << out_dir_ << "/bs_abundance_*.tsv" << std::endl;
  }
  int i;
  for (i = 0; i < n_bs_; ++i) {
    if (i % 50 == 0 && i > 0) {
      std::cerr << std::endl;
    }
    std::cerr << ".";
    std::cerr.flush();
    std::string bs_out_fname( out_dir_ + "/bs_abundance_" + std::to_string(i) +
        ".tsv" );
    rw_from_counts(bs_, "bs" + std::to_string(i), bs_out_fname);
  }

  if (i-1 % 50 != 0 && i > 0) {
    std::cerr << std::endl;
  }
}

void H5Converter::rw_from_counts(hid_t group_id, const std::string& count_name, const std::string& out_fname) {
  std::vector<double> alpha;
  read_dataset(group_id, count_name.c_str(), alpha);

  std::vector<uint32_t> lengths_2(std::begin(lengths_), std::end(lengths_));
  plaintext_writer(out_fname, targ_ids_, alpha, eff_lengths_, lengths_2);  
}

#endif // USE_HDF5