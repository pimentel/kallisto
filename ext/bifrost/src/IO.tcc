#ifndef BIFROST_IO_CDBG_TCC
#define BIFROST_IO_CDBG_TCC

template<typename U, typename G>
bool CompactedDBG<U, G>::write(const std::string& output_fn, const size_t nb_threads, const bool GFA_output, const bool FASTA_output, const bool BFG_output, const bool write_index_file, const bool compressed_output, const bool verbose) const {

    if (invalid){

        std::cerr << "CompactedDBG::write(): Graph is invalid and cannot be written to disk" << std::endl;
        return false;
    }

    if (nb_threads <= 0){

        std::cerr << "CompactedDBG::write(): Number of threads cannot be less than 0" << std::endl;
        return false;
    }

    if (nb_threads > std::thread::hardware_concurrency()){

        std::cerr << "CompactedDBG::write(): Number of threads cannot exceed " << std::thread::hardware_concurrency() << "threads" << std::endl;
        return false;
    }

    if (!GFA_output && !FASTA_output && !BFG_output){

        std::cerr << "CompactedDBG::write(): No type of format output selected" << std::endl;
        return false;
    }

    if (static_cast<int>(GFA_output) + static_cast<int>(FASTA_output) + static_cast<int>(BFG_output) > 1){

        std::cerr << "CompactedDBG::write(): Multiple output formats selected. Please choose one." << std::endl;
        return false;
    }

    bool write_success = true;

    {
        if (verbose) std::cout << std::endl << "CompactedDBG::write(): Writing graph to disk" << std::endl;

        std::string fn = output_fn;

        // Add file extensions if missing
        if (GFA_output || FASTA_output) {

            const std::string g_ext = (GFA_output ? ".gfa" : ".fasta");
            const std::string c_ext = ".gz";
            const std::string gc_ext = g_ext + c_ext;

            const size_t pos_ext = fn.find_last_of(".");

            if (pos_ext == std::string::npos) fn.append(compressed_output ? gc_ext : g_ext);
            else if (!compressed_output && (fn.substr(pos_ext) != g_ext)) fn.append(g_ext);
            else if (compressed_output && (fn.substr(pos_ext) != c_ext)) fn.append(gc_ext);
        }
        else if (BFG_output) {

            const std::string g_ext = ".bfg";
            const size_t pos_ext = fn.find_last_of(".");

            if ((pos_ext == std::string::npos) || (fn.substr(pos_ext) != g_ext)) fn.append(g_ext);


        }

        FILE* fp = fopen(fn.c_str(), "w");

        if (fp == NULL) {

            std::cerr << "CompactedDBG::write(): Could not open file " << fn << " for writing graph" << std::endl;

            return false;
        }
        else {

            fclose(fp);

            if (std::remove(fn.c_str()) != 0) std::cerr << "CompactedDBG::write(): Could not remove temporary file " << fn << std::endl;
        }

        if (GFA_output) write_success = writeGFA(fn, nb_threads, compressed_output);
        else if (FASTA_output) write_success = writeFASTA(fn, compressed_output);
        else if (BFG_output) write_success = writeBinaryGraph(fn, nb_threads);
    }

    if (write_success && (write_index_file || BFG_output)) {

        if (verbose) std::cout << std::endl << "CompactedDBG::write(): Writing index file to disk" << std::endl;

        std::string fn = output_fn;

        // Add file extensions if missing
        {
            const std::string ext = ".bfi";

            if ((fn.length() < ext.length()) || (fn.substr(fn.length() - ext.length()) != ext)) fn.append(ext);
        }

        FILE* fp = fopen(fn.c_str(), "w");

        if (fp == NULL) {

            std::cerr << "CompactedDBG::write(): Could not open file " << fn << " for writing index file" << std::endl;

            return false;
        }
        else {

            fclose(fp);

            if (std::remove(fn.c_str()) != 0) std::cerr << "CompactedDBG::write(): Could not remove temporary file " << fn << std::endl;
        }

        const uint64_t graph_checksum = checksum();

        write_success = writeBinaryIndex(fn, graph_checksum, nb_threads);
    }

    return write_success;
}

template<typename U, typename G>
bool CompactedDBG<U, G>::read(const std::string& input_graph_fn, const size_t nb_threads, const bool verbose){

    if (verbose) std::cout << std::endl << "CompactedDBG::read(): Reading graph from disk" << std::endl;

    const int format = FileParser::getFileFormat(input_graph_fn.c_str());

    if (format == -1){

        std::cerr << "CompactedDBG::read(): Input graph file " << input_graph_fn << " does not exist, is ill-formed or is not a valid graph file format." << std::endl;

        return false;
    }
    else if ((format != 0) && (format != 2) && (format != 3)){

        std::cerr << "CompactedDBG::read(): Input graph file must be in FASTA, GFA or BFG format." << std::endl;

        return false;
    }

    std::string input_index_fn = input_graph_fn;

    // Try to open index file if available
    {
        size_t pos_ext = input_index_fn.find_last_of(".");

        std::string ext;

        if ((pos_ext != std::string::npos) && (input_index_fn.substr(pos_ext) == ".gz")) {

            input_index_fn = input_index_fn.substr(0, pos_ext);
            pos_ext = input_index_fn.find_last_of(".");
        }

        if (pos_ext != std::string::npos) ext = input_index_fn.substr(pos_ext);
        if ((pos_ext != std::string::npos) && ((ext == ".gfa")) || (ext == ".fasta") || (ext == ".bfg")) input_index_fn = input_index_fn.substr(0, pos_ext);

        input_index_fn += ".bfi"; // Add Bifrost graph index extension
    }

    if ((input_graph_fn != input_index_fn) && check_file_exists(input_index_fn)) {

        if (verbose) std::cout << "CompactedDBG::read(): Reading using Bifrost index file " << input_index_fn << "." << std::endl;

        return read(input_graph_fn, input_index_fn, nb_threads, verbose);
    }
    else {

        if (format == 0) { // FASTA input

            const int k = k_;
            const int g = g_;

            clear();

            {
                KmerStream_Build_opt kms_opt;

                kms_opt.threads = nb_threads;
                kms_opt.verbose = verbose;
                kms_opt.k = k;
                kms_opt.g = g;
                kms_opt.q = 0;

                kms_opt.files.push_back(input_graph_fn);

                KmerStream kms(kms_opt);

                MinimizerIndex hmap_min_unitigs_tmp(std::max(size_t(1), kms.MinimizerF0()) * 1.05);

                hmap_min_unitigs = std::move(hmap_min_unitigs_tmp);
            }

            setKmerGmerLength(k, g);
            makeGraphFromFASTA(input_graph_fn, nb_threads);
        }
        else if (format == 2){ // GFA format

            std::string header;

            int k = k_, g = g_;

            {
                FILE* fp = fopen(input_graph_fn.c_str(), "r");

                if (fp == NULL) {

                    std::cerr << "CompactedDBG::read(): Could not open file " << input_graph_fn << " for reading graph" << std::endl;
                    return false;
                }

                fclose(fp);
            }

            {
                size_t fn_id = 0;

                GFA_Parser gfap(input_graph_fn);

                header = gfap.open_read().first;
            }

            {
                if (header[0] != 'H'){

                    std::cerr << "CompactedDBG::read(): An error occurred while reading input GFA file." << std::endl;
                    return false;
                }

                std::stringstream hs(header.c_str() + 2); // Skip the first 2 char. of the line "H\t"
                std::string sub;

                while (hs.good()){ // Split line based on tabulation

                    getline(hs, sub, '\t');

                    const std::string tag = sub.substr(0, 5);

                    if (tag == "KL:Z:") k = atoi(sub.c_str() + 5);
                    else if (tag == "ML:Z:") g = atoi(sub.c_str() + 5);
                }

                clear();
            }

            {
                KmerStream_Build_opt kms_opt;

                kms_opt.threads = nb_threads;
                kms_opt.verbose = verbose;
                kms_opt.k = k;
                kms_opt.g = g;
                kms_opt.q = 0;

                kms_opt.files.push_back(input_graph_fn);

                KmerStream kms(kms_opt);

                MinimizerIndex hmap_min_unitigs_tmp(std::max(size_t(1), kms.MinimizerF0()) * 1.05);
                hmap_min_unitigs = std::move(hmap_min_unitigs_tmp);
            }

            setKmerGmerLength(k, g);
            makeGraphFromGFA(input_graph_fn, nb_threads);
        }
        else { // BFG format

            std::cerr << "CompactedDBG::read(): No index found for Bifrost graph file " << input_graph_fn << std::endl;

            invalid = true;
        }

        setFullCoverage(1);

        // Set coverages
        if (!invalid) for (auto& unitig : *this) unitig.setFullCoverage();

        if (verbose) std::cout << std::endl << "CompactedDBG::read(): Finished reading graph from disk" << std::endl;

        return !invalid;
    }
}

template<typename U, typename G>
bool CompactedDBG<U, G>::read(const std::string& input_graph_fn, const std::string& input_index_fn, const size_t nb_threads, const bool verbose){

    if (verbose) std::cout << std::endl << "CompactedDBG::read(): Reading graph from disk" << std::endl;

    const int format_graph = FileParser::getFileFormat(input_graph_fn.c_str());
    const int format_index = FileParser::getFileFormat(input_index_fn.c_str());

    if (format_graph == -1){

        std::cerr << "CompactedDBG::read(): Input graph file " << input_graph_fn << " does not exist, is ill-formed or is not a valid graph file format." << std::endl;

        return false;
    }
    else if ((format_graph != 0) && (format_graph != 2) && (format_graph != 3)){

        std::cerr << "CompactedDBG::read(): Input graph file must be in FASTA, GFA or BFG format." << std::endl;

        return false;
    }

    if (format_index != 4) {

        //std::cerr << format_index << std::endl;
        std::cerr << "CompactedDBG::read(): Input index file " << input_index_fn << " does not exist, is ill-formed or is not a valid index file format." << std::endl;

        return false;
    }

    std::pair<uint64_t, bool> p_readSuccess_checksum;

    if (format_graph == 0) { // FASTA input

        const int k = k_;
        const int g = g_;

        clear();
        setKmerGmerLength(k, g);

        p_readSuccess_checksum = readGraphFromIndexFASTA(input_graph_fn, input_index_fn, k_, g_);
        invalid = !p_readSuccess_checksum.second;
    }
    else if (format_graph == 2){ // GFA format

        std::string header;

        int k = k_, g = g_;

        {
            FILE* fp = fopen(input_graph_fn.c_str(), "r");

            if (fp == NULL) {

                std::cerr << "CompactedDBG::read(): Could not open file " << input_graph_fn << " for reading graph" << std::endl;
                return false;
            }

            fclose(fp);
        }

        {
            size_t fn_id = 0;

            GFA_Parser gfap(input_graph_fn);

            header = gfap.open_read().first;
        }

        {
            if (header[0] != 'H'){

                std::cerr << "CompactedDBG::read(): An error occurred while reading input GFA file." << std::endl;
                return false;
            }

            std::stringstream hs(header.c_str() + 2); // Skip the first 2 char. of the line "H\t"
            std::string sub;

            while (hs.good()){ // Split line based on tabulation

                getline(hs, sub, '\t');

                const std::string tag = sub.substr(0, 5);

                if (tag == "KL:Z:") k = atoi(sub.c_str() + 5);
                else if (tag == "ML:Z:") g = atoi(sub.c_str() + 5);
            }

            clear();
            setKmerGmerLength(k, g);
        }

        p_readSuccess_checksum = readGraphFromIndexGFA(input_graph_fn, input_index_fn, k_, g_);
        invalid = !p_readSuccess_checksum.second;
    }
    else { // BFG (binary) format

        p_readSuccess_checksum = readBinaryGraph(input_graph_fn);
        invalid = !p_readSuccess_checksum.second;
    }

    setFullCoverage(1);

    if (!invalid) {

        for (auto& unitig : *this) unitig.setFullCoverage();

        invalid = !readBinaryIndex(input_index_fn, p_readSuccess_checksum.first);
    }

    if (verbose) std::cout << std::endl << "CompactedDBG::read(): Finished reading graph from disk" << std::endl;

    return !invalid;
}

template<typename U, typename G>
template<bool is_void>
typename std::enable_if<!is_void, void>::type CompactedDBG<U, G>::writeGFA_sequence_(GFA_Parser& graph, KmerHashTable<size_t>& idmap) const {

    size_t labelA = 1;

    for (const auto& unitig : *this){

        const std::string seq(unitig.referenceUnitigToString());

        graph.write_sequence(std::to_string(labelA), seq.size(), seq, /*unitig.getData()->serialize(unitig)*/"");

        if (unitig.isAbundant) idmap.insert(Kmer(unitig.referenceUnitigToString().c_str()), labelA);

        ++labelA;
    }
}

template<typename U, typename G>
template<bool is_void>
typename std::enable_if<is_void, void>::type CompactedDBG<U, G>::writeGFA_sequence_(GFA_Parser& graph, KmerHashTable<size_t>& idmap) const {

    size_t labelA = 1;

    for (const auto& unitig : *this){

        const std::string seq(unitig.referenceUnitigToString());

        graph.write_sequence(std::to_string(labelA), seq.size(), seq, "");

        if (unitig.isAbundant) idmap.insert(Kmer(unitig.referenceUnitigToString().c_str()), labelA);

        ++labelA;
    }
}

// It is very important to write unitigs to disk in the same following order:
// 1 - All unitigs with length > k
// 2 - All unitigs with length == k which do not have abundant minimizers
// 3 - All unitigs with length == k which have abundant minimizers
// The binary graph file is written in that order
// and the checksum is stored in the index file is computed for that order
template<typename U, typename G>
bool CompactedDBG<U, G>::writeGFA(const std::string& fn, const size_t nb_threads, const bool compressed_output) const {

    const size_t v_unitigs_sz = v_unitigs.size();
    const size_t v_kmers_sz = km_unitigs.size();

    size_t labelA, labelB, id = v_unitigs_sz + v_kmers_sz + 1;

    const std::string header_tag("BV:Z:" + std::string(BFG_VERSION) + "\t" + "KL:Z:" + std::to_string(k_) + "\t" + "ML:Z:" + std::to_string(g_));

    KmerHashTable<size_t> idmap(h_kmers_ccov.size());

    GFA_Parser graph(fn);

    graph.open_write(1, header_tag, compressed_output);

    writeGFA_sequence_<std::is_void<U>::value>(graph, idmap);

    if (nb_threads == 1){

        for (labelA = 1; labelA <= v_unitigs_sz; labelA++) {

            const Unitig<U>* unitig = v_unitigs[labelA - 1];
            const Kmer head = unitig->getSeq().getKmer(0);
            const Kmer tail = unitig->getSeq().getKmer(unitig->length() - k_);

            const std::vector<const_UnitigMap<U, G>> pred = findPredecessors(head, true);
            const std::vector<const_UnitigMap<U, G>> succ = findSuccessors(tail, 4, true);

            for (const auto& unitig : pred){

                if (!unitig.isEmpty){

                    if (unitig.isAbundant) labelB = *(idmap.find(unitig.getUnitigHead().rep()));
                    else labelB = unitig.pos_unitig + 1 + ((static_cast<size_t>(!unitig.isShort) - 1) & v_unitigs_sz);

                    const std::string slabelA = std::to_string(labelA);
                    const std::string slabelB = std::to_string(labelB);

                    const size_t pos = (static_cast<size_t>(unitig.strand) - 1) & (unitig.size - k_ + 1);

                    graph.write_edge(slabelA, 0, k_-1, false, slabelB, pos, pos + k_ - 1, !unitig.strand);
                }
            }

            for (const auto& unitig : succ){

                if (!unitig.isEmpty){

                    if (unitig.isAbundant) labelB = *(idmap.find(unitig.getUnitigHead().rep()));
                    else labelB = unitig.pos_unitig + 1 + ((static_cast<size_t>(!unitig.isShort) - 1) & v_unitigs_sz);

                    const std::string slabelA = std::to_string(labelA);
                    const std::string slabelB = std::to_string(labelB);

                    const size_t pos = (static_cast<size_t>(unitig.strand) - 1) & (unitig.size - k_ + 1);

                    graph.write_edge(slabelA, unitig.size - k_ + 1, unitig.size, true, slabelB, pos, pos + k_ - 1, unitig.strand);
                }
            }
        }

        for (labelA = v_unitigs_sz + 1; labelA <= v_kmers_sz + v_unitigs_sz; labelA++) {

            const Kmer km_unitig = km_unitigs.getKmer(labelA - v_unitigs_sz - 1);

            const std::vector<const_UnitigMap<U, G>> pred = findPredecessors(km_unitig, true);
            const std::vector<const_UnitigMap<U, G>> succ = findSuccessors(km_unitig, 4, true);

            for (const auto& unitig : pred){

                if (!unitig.isEmpty){

                    if (unitig.isAbundant) labelB = *(idmap.find(unitig.getUnitigHead().rep()));
                    else labelB = unitig.pos_unitig + 1 + ((static_cast<size_t>(!unitig.isShort) - 1) & v_unitigs_sz);

                    const std::string slabelA = std::to_string(labelA);
                    const std::string slabelB = std::to_string(labelB);

                    const size_t pos = (static_cast<size_t>(unitig.strand) - 1) & (unitig.size - k_ + 1);

                    graph.write_edge(slabelA, 0, k_-1, false, slabelB, pos, pos + k_ - 1, !unitig.strand);
                }
            }

            for (const auto& unitig : succ){

                if (!unitig.isEmpty){

                    if (unitig.isAbundant) labelB = *(idmap.find(unitig.getUnitigHead().rep()));
                    else labelB = unitig.pos_unitig + 1 + ((static_cast<size_t>(!unitig.isShort) - 1) & v_unitigs_sz);

                    const std::string slabelA = std::to_string(labelA);
                    const std::string slabelB = std::to_string(labelB);

                    const size_t pos = (static_cast<size_t>(unitig.strand) - 1) & (unitig.size - k_ + 1);

                    graph.write_edge(slabelA, unitig.size - k_ + 1, unitig.size, true, slabelB, pos, pos + k_ - 1, unitig.strand);
                }
            }
        }

        for (KmerHashTable<size_t>::iterator it = idmap.begin(); it != idmap.end(); it++) {

            labelA = *it;

            const std::vector<const_UnitigMap<U, G>> pred = findPredecessors(it.getKey(), true);
            const std::vector<const_UnitigMap<U, G>> succ = findSuccessors(it.getKey(), 4, true);

            for (const auto& unitig : pred){

                if (!unitig.isEmpty){

                    if (unitig.isAbundant) labelB = *(idmap.find(unitig.getUnitigHead().rep()));
                    else labelB = unitig.pos_unitig + 1 + ((static_cast<size_t>(!unitig.isShort) - 1) & v_unitigs_sz);

                    const std::string slabelA = std::to_string(labelA);
                    const std::string slabelB = std::to_string(labelB);

                    const size_t pos = (static_cast<size_t>(unitig.strand) - 1) & (unitig.size - k_ + 1);

                    graph.write_edge(slabelA, 0, k_-1, false, slabelB, pos, pos + k_ - 1, !unitig.strand);
                }
            }

            for (const auto& unitig : succ){

                if (!unitig.isEmpty){

                    if (unitig.isAbundant) labelB = *(idmap.find(unitig.getUnitigHead().rep()));
                    else labelB = unitig.pos_unitig + 1 + ((static_cast<size_t>(!unitig.isShort) - 1) & v_unitigs_sz);

                    const std::string slabelA = std::to_string(labelA);
                    const std::string slabelB = std::to_string(labelB);

                    const size_t pos = (static_cast<size_t>(unitig.strand) - 1) & (unitig.size - k_ + 1);

                    graph.write_edge(slabelA, unitig.size - k_ + 1, unitig.size, true, slabelB, pos, pos + k_ - 1, unitig.strand);
                }
            }
        }
    }
    else {

        const size_t chunk_size = 1024;

        auto worker_v_unitigs = [v_unitigs_sz, &idmap, this](const size_t labelA_start, const size_t labelA_end,
                                                             std::vector<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>>* v_out){

            // We need to deal with the tail of long unitigs
            for (size_t labelA = labelA_start; labelA < labelA_end; ++labelA) {

                const Unitig<U>* unitig = v_unitigs[labelA - 1];

                const Kmer head = unitig->getSeq().getKmer(0);
                const Kmer tail = unitig->getSeq().getKmer(unitig->length() - k_);

                const std::vector<const_UnitigMap<U, G>> pred = this->findPredecessors(head, true);
                const std::vector<const_UnitigMap<U, G>> succ = this->findSuccessors(tail, 4, true);

                for (const auto& um : pred) {

                    if (!um.isEmpty){

                        const size_t labelB = (um.isAbundant ?  *(idmap.find(um.getUnitigHead().rep())) :
                                                                um.pos_unitig + 1 + ((static_cast<size_t>(!um.isShort) - 1) & v_unitigs_sz));

                        v_out->push_back(std::make_pair(std::make_pair(labelA, false), std::make_pair(labelB, !um.strand)));
                    }
                }

                for (const auto& um : succ) {

                    if (!um.isEmpty){

                        const size_t labelB = (um.isAbundant ?  *(idmap.find(um.getUnitigHead().rep())) :
                                                                um.pos_unitig + 1 + ((static_cast<size_t>(!um.isShort) - 1) & v_unitigs_sz));

                        v_out->push_back(std::make_pair(std::make_pair(labelA, true), std::make_pair(labelB, um.strand)));
                    }
                }
            }
        };

        auto worker_v_kmers = [v_unitigs_sz, &idmap, this](const size_t labelA_start, const size_t labelA_end,
                                                             std::vector<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>>* v_out){

            // We need to deal with the tail of long unitigs
            for (size_t labelA = labelA_start; labelA < labelA_end; ++labelA) {

                const Kmer km_unitig = km_unitigs.getKmer(labelA - v_unitigs_sz - 1);

                const std::vector<const_UnitigMap<U, G>> pred = this->findPredecessors(km_unitig, true);
                const std::vector<const_UnitigMap<U, G>> succ = this->findSuccessors(km_unitig, 4, true);

                for (const auto& um : pred) {

                    if (!um.isEmpty){

                        const size_t labelB = (um.isAbundant ?  *(idmap.find(um.getUnitigHead().rep())) :
                                                                um.pos_unitig + 1 + ((static_cast<size_t>(!um.isShort) - 1) & v_unitigs_sz));

                        v_out->push_back(std::make_pair(std::make_pair(labelA, false), std::make_pair(labelB, !um.strand)));
                    }
                }

                for (const auto& um : succ) {

                    if (!um.isEmpty){

                        const size_t labelB = (um.isAbundant ?  *(idmap.find(um.getUnitigHead().rep())) :
                                                                um.pos_unitig + 1 + ((static_cast<size_t>(!um.isShort) - 1) & v_unitigs_sz));

                        v_out->push_back(std::make_pair(std::make_pair(labelA, true), std::make_pair(labelB, um.strand)));
                    }
                }
            }
        };

        auto worker_v_abundant = [v_unitigs_sz, chunk_size, &idmap, this](  KmerHashTable<size_t>::iterator* l_it,
                                                                            std::vector<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>>* v_out){

            KmerHashTable<size_t>::iterator& it = *l_it;

            // We need to deal with the tail of long unitigs
            for (size_t i = 0; (it != idmap.end()) && (i < chunk_size); ++i, ++it) {

                const std::vector<const_UnitigMap<U, G>> pred = this->findPredecessors(it.getKey(), true);
                const std::vector<const_UnitigMap<U, G>> succ = this->findSuccessors(it.getKey(), 4, true);

                for (const auto& um : pred) {

                    if (!um.isEmpty){

                        const size_t labelB = (um.isAbundant ?  *(idmap.find(um.getUnitigHead().rep())) :
                                                                um.pos_unitig + 1 + ((static_cast<size_t>(!um.isShort) - 1) & v_unitigs_sz));

                        v_out->push_back(std::make_pair(std::make_pair(*it, false), std::make_pair(labelB, !um.strand)));
                    }
                }

                for (const auto& um : succ) {

                    if (!um.isEmpty){

                        const size_t labelB = (um.isAbundant ?  *(idmap.find(um.getUnitigHead().rep())) :
                                                                um.pos_unitig + 1 + ((static_cast<size_t>(!um.isShort) - 1) & v_unitigs_sz));

                        v_out->push_back(std::make_pair(std::make_pair(*it, true), std::make_pair(labelB, um.strand)));
                    }
                }
            }
        };

        {
            std::atomic<size_t> label(1);

            std::vector<std::vector<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>>> v_out(nb_threads);
            std::vector<std::thread> workers; // need to keep track of threads so we can join them

            std::mutex mutex_file;

            for (size_t t = 0; t < nb_threads; ++t){

                workers.emplace_back(

                    [&, t]{

                        while (true) {

                            const size_t old_labelA = label.fetch_add(chunk_size);

                            if (old_labelA <= v_unitigs_sz){

                                if (old_labelA + chunk_size <= v_unitigs_sz) worker_v_unitigs(old_labelA, old_labelA + chunk_size, &v_out[t]);
                                else worker_v_unitigs(old_labelA, v_unitigs_sz + 1, &v_out[t]);

                                {
                                    std::unique_lock<std::mutex> lock(mutex_file);

                                    for (const auto& p : v_out[t]){

                                        const std::string slabelA = std::to_string(p.first.first);
                                        const std::string slabelB = std::to_string(p.second.first);

                                        graph.write_edge(slabelA, 0, k_-1, p.first.second, slabelB, 0, k_-1, p.second.second);
                                    }
                                }

                                v_out[t].clear();
                            }
                            else return;
                        }
                    }
                );
            }

            for (auto& t : workers) t.join();
        }

        {
            const size_t v_kmers_unitigs_sz = v_kmers_sz + v_unitigs_sz;

            std::atomic<size_t> label(v_unitigs_sz + 1);

            std::vector<std::vector<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>>> v_out(nb_threads);
            std::vector<std::thread> workers; // need to keep track of threads so we can join them

            std::mutex mutex_file;

            for (size_t t = 0; t < nb_threads; ++t){

                workers.emplace_back(

                    [&, t]{

                        while (true) {

                            const size_t old_labelA = label.fetch_add(chunk_size);

                            if (old_labelA <= v_kmers_unitigs_sz){

                                if (old_labelA + chunk_size <= v_kmers_unitigs_sz) worker_v_kmers(old_labelA, old_labelA + chunk_size, &v_out[t]);
                                else worker_v_kmers(old_labelA, v_kmers_unitigs_sz + 1, &v_out[t]);

                                {
                                    std::unique_lock<std::mutex> lock(mutex_file);

                                    for (const auto& p : v_out[t]){

                                        const std::string slabelA = std::to_string(p.first.first);
                                        const std::string slabelB = std::to_string(p.second.first);

                                        graph.write_edge(slabelA, 0, k_-1, p.first.second, slabelB, 0, k_-1, p.second.second);
                                    }
                                }

                                v_out[t].clear();
                            }
                            else return;
                        }
                    }
                );
            }

            for (auto& t : workers) t.join();
        }

        {
            KmerHashTable<size_t>::iterator it = idmap.begin(), it_end = idmap.end();

            std::vector<std::vector<std::pair<std::pair<size_t, bool>, std::pair<size_t, bool>>>> v_out(nb_threads);
            std::vector<std::thread> workers; // need to keep track of threads so we can join them

            std::mutex mutex_file, mutex_it;

            for (size_t t = 0; t < nb_threads; ++t){

                workers.emplace_back(

                    [&, t]{

                        KmerHashTable<size_t>::iterator l_it;

                        bool stop;

                        while (true) {

                            {
                                std::unique_lock<std::mutex> lock(mutex_it);

                                l_it = it;

                                for (size_t i = 0; (it != it_end) && (i < chunk_size); ++i, ++it){}

                                stop = (l_it == it_end) && (it == it_end);
                            }

                            if (!stop){

                                worker_v_abundant(&l_it, &v_out[t]);

                                {
                                    std::unique_lock<std::mutex> lock(mutex_file);

                                    for (const auto& p : v_out[t]){

                                        const std::string slabelA = std::to_string(p.first.first);
                                        const std::string slabelB = std::to_string(p.second.first);

                                        graph.write_edge(slabelA, 0, k_-1, p.first.second, slabelB, 0, k_-1, p.second.second);
                                    }
                                }

                                v_out[t].clear();
                            }
                            else return;
                        }
                    }
                );
            }

            for (auto& t : workers) t.join();
        }
    }

    graph.close();

    return true;
}

// It is very important to write unitigs to disk in the same following order:
// 1 - All unitigs with length > k
// 2 - All unitigs with length == k which do not have abundant minimizers
// 3 - All unitigs with length == k which have abundant minimizers
// The binary graph file is written in that order
// and the checksum is stored in the index file is computed for that order
template<typename U, typename G>
bool CompactedDBG<U, G>::writeFASTA(const std::string& fn, const bool compressed_output) const {

    const size_t v_unitigs_sz = v_unitigs.size();
    const size_t v_kmers_sz = km_unitigs.size();

    size_t i = 0;

    bool write_success;

    if (compressed_output) {

        zstr::ofstream gout(fn, std::ios_base::out);

        for (size_t j = 0; !gout.fail() && (j < v_unitigs_sz); ++j, ++i) gout << ">" << i << "\n" << v_unitigs[j]->getSeq().toString() << "\n";
        for (size_t j = 0; !gout.fail() && (j < v_kmers_sz); ++j, ++i) gout << ">" << i << "\n" << km_unitigs.getKmer(j).toString() << "\n";

        for (typename h_kmers_ccov_t::const_iterator it = h_kmers_ccov.begin(); !gout.fail() && (it != h_kmers_ccov.end()); ++it, ++i) {

            gout << ">" << i << "\n" << it.getKey().toString() << "\n";
        }

        write_success = !gout.fail();
    }
    else {

        std::ofstream graphfile;
        std::ostream graph(0);

        graphfile.open(fn.c_str());
        graph.rdbuf(graphfile.rdbuf());

        for (size_t j = 0; !graph.fail() && (j < v_unitigs_sz); ++j, ++i) graph << ">" << i << "\n" << v_unitigs[j]->getSeq().toString() << "\n";
        for (size_t j = 0; !graph.fail() && (j < v_kmers_sz); ++j, ++i) graph << ">" << i << "\n" << km_unitigs.getKmer(j).toString() << "\n";

        for (typename h_kmers_ccov_t::const_iterator it = h_kmers_ccov.begin(); !graph.fail() && (it != h_kmers_ccov.end()); ++it, ++i) {

            graph << ">" << i << "\n" << it.getKey().toString() << "\n";
        }

        write_success = !graph.fail();

        graphfile.close();
    }

    return write_success;
}

template<typename U, typename G>
bool CompactedDBG<U, G>::readBinary(const std::string& fn, bool static_m, uint32_t threads) {

    if ((fn.length() == 0) || !check_file_exists(fn)) return false;

    std::vector<Minimizer> minz;
    if (static_m) {
        std::ifstream infile;
        std::istream in(0);

        infile.open(fn.c_str());
        in.rdbuf(infile.rdbuf());
        if (!in.fail()) {
            const std::pair<uint64_t, bool> p_readSuccess_checksum = readBinaryGraph(in);
            readBinaryMinimizers(in, p_readSuccess_checksum.first, minz, threads);
        }
        infile.close();

    }

    std::ifstream infile;
    std::istream in(0);

    infile.open(fn.c_str());
    in.rdbuf(infile.rdbuf());

    return readBinary(in, minz, threads);
}

template<typename U, typename G>
bool CompactedDBG<U, G>::readBinary(std::istream& in, std::vector<Minimizer>& minz, uint32_t threads) {

    if (!in.fail()) {

        const std::pair<uint64_t, bool> p_readSuccess_checksum = readBinaryGraph(in);

        if (p_readSuccess_checksum.second) return readBinaryIndex(in, p_readSuccess_checksum.first, minz, threads);
    }

    return false;
}

template<typename U, typename G>
bool CompactedDBG<U, G>::readBinary(std::istream& in, boophf_t* mphf, uint32_t threads) {

    if (!in.fail()) {

        const std::pair<uint64_t, bool> p_readSuccess_checksum = readBinaryGraph(in);

        if (p_readSuccess_checksum.second) return readBinaryIndex(in, p_readSuccess_checksum.first, mphf, threads);
    }

    return false;
}

template<typename U, typename G>
bool CompactedDBG<U, G>::readMinimizers(std::istream& in, std::vector<Minimizer>& minz, uint32_t threads) {

    if (!in.fail()) {
        const std::pair<uint64_t, bool> p_readSuccess_checksum = readBinaryGraph(in);
        return readBinaryMinimizers(in, p_readSuccess_checksum.first, minz, threads);
    }
    return false;
}

template<typename U, typename G>
size_t CompactedDBG<U, G>::writeMinimizers(std::ostream& out) {

    size_t n_written = 0;
    auto it = hmap_min_unitigs.begin();
    auto end = hmap_min_unitigs.end();
    while (it != end) {

        it.getKey().rep().write(out);
        ++n_written;
        ++it;
    }
    return n_written;
}

template<typename U, typename G>
void CompactedDBG<U, G>::clearAndGetMinimizers(std::vector<Minimizer>& minz) {

    hmap_min_unitigs.clearPTV();

    minz.clear();
    minz.reserve(hmap_min_unitigs.size());
    size_t n_written = 0;
    auto it = hmap_min_unitigs.begin();
    auto end = hmap_min_unitigs.end();
    while (it != end) {

        minz.push_back(std::move(it.getKey().rep()));
        ++n_written;
        ++it;
    }

    hmap_min_unitigs.clear();
}

template<typename U, typename G>
bool CompactedDBG<U, G>::readBinaryIndex(const std::string& fn, const uint64_t checksum, bool static_m, uint32_t threads) {

    if ((fn.length() == 0) || !check_file_exists(fn)) return false;

    std::ifstream infile;
    std::istream in(0);

    infile.open(fn.c_str());
    in.rdbuf(infile.rdbuf());

    std::vector<Minimizer> minz;
    return readBinaryIndex(in, checksum, minz, threads);
}

template<typename U, typename G>
bool CompactedDBG<U, G>::readBinaryIndexHead(const std::string& fn, size_t& file_format_version, size_t& v_unitigs_sz, size_t& km_unitigs_sz,
                                            size_t& h_kmers_ccov_sz, size_t& hmap_min_unitigs_sz, uint64_t& read_checksum) const {

    if ((fn.length() == 0) || !check_file_exists(fn)) return false;

    std::ifstream infile;
    std::istream in(0);

    infile.open(fn.c_str());
    in.rdbuf(infile.rdbuf());

    return readBinaryIndexHead(in, file_format_version, v_unitigs_sz, km_unitigs_sz, h_kmers_ccov_sz, hmap_min_unitigs_sz, read_checksum);
}

template<typename U, typename G>
bool CompactedDBG<U, G>::readBinaryIndexHead(std::istream& in, size_t& file_format_version, size_t& v_unitigs_sz, size_t& km_unitigs_sz,
                                            size_t& h_kmers_ccov_sz, size_t& hmap_min_unitigs_sz, uint64_t& read_checksum) const {

    if (in.fail()) return false;

    in.read(reinterpret_cast<char*>(&file_format_version), sizeof(size_t));
    in.read(reinterpret_cast<char*>(&read_checksum), sizeof(uint64_t));

    in.read(reinterpret_cast<char*>(&v_unitigs_sz), sizeof(size_t));
    in.read(reinterpret_cast<char*>(&km_unitigs_sz), sizeof(size_t));
    in.read(reinterpret_cast<char*>(&h_kmers_ccov_sz), sizeof(size_t));
    in.read(reinterpret_cast<char*>(&hmap_min_unitigs_sz), sizeof(size_t));

    return !in.fail();
}

template<typename U, typename G>
bool CompactedDBG<U, G>::readBinaryMinimizers(std::istream& in, const uint64_t checksum,
                                              std::vector<Minimizer>& minz, uint32_t threads) {


    // Build a static MPHF preemptively
    bool read_success = !in.fail();

    // 0 - Write file format version, checksum and number of minimizers
    if (read_success) {

        size_t file_format_version = 0, v_unitigs_sz = 0, km_unitigs_sz = 0, h_kmers_ccov_sz = 0, hmap_min_unitigs_sz = 0;
        uint64_t read_checksum = 0;

        read_success = readBinaryIndexHead(in, file_format_version, v_unitigs_sz, km_unitigs_sz, h_kmers_ccov_sz, hmap_min_unitigs_sz, read_checksum);

        if (!read_success || ((file_format_version >> 32) != BFG_METABIN_FORMAT_HEADER)) return false;

    }

    if (read_success) {

        size_t nb_bmp_unitigs = 0;

        in.read(reinterpret_cast<char*>(&nb_bmp_unitigs), sizeof(size_t));

        std::vector<BitContainer> v_bmp_unitigs(nb_bmp_unitigs);

        for (size_t i = 0; read_success && (i < nb_bmp_unitigs); ++i) read_success = v_bmp_unitigs[i].read(in);

        if (read_success && !v_bmp_unitigs.empty() && !v_bmp_unitigs.front().isEmpty()) {

            size_t unitig_id = 0,
                   tot_unitig_len = 0,
                   curr_unitig_len = v_unitigs[0]->getSeq().size() - g_ + 1;

            for (size_t i = 0; i < nb_bmp_unitigs; ++i) {

                const size_t id_bmp = i << 32;

                for (const auto pos_bmp : v_bmp_unitigs[i]) {

                    const size_t pos = id_bmp + pos_bmp;

                    while (pos >= (tot_unitig_len + curr_unitig_len)) {

                        ++unitig_id;
                        tot_unitig_len += curr_unitig_len;

                        curr_unitig_len = v_unitigs[unitig_id]->getSeq().size() - g_ + 1;
                    }

                    const size_t relative_pos = pos - tot_unitig_len;

                    minz.push_back(v_unitigs[unitig_id]->getSeq().getMinimizer(relative_pos).rep());
                }
                v_bmp_unitigs[i].clear();
            }
        }
    }

    if (read_success) {

        size_t nb_bmp_km = 0;

        in.read(reinterpret_cast<char*>(&nb_bmp_km), sizeof(size_t));

        std::vector<BitContainer> v_bmp_km(nb_bmp_km);

        for (size_t i = 0; (i < nb_bmp_km) && read_success; ++i) read_success = v_bmp_km[i].read(in);

        if (read_success && !v_bmp_km.empty() && !v_bmp_km.front().isEmpty()) {

            const size_t km_glen = k_ - g_ + 1;

            for (size_t i = 0; i < nb_bmp_km; ++i) {

                const size_t id_bmp = i << 32;

                for (const auto pos_bmp : v_bmp_km[i]) {

                    const size_t pos = id_bmp + pos_bmp;

                    const size_t km_id = pos / km_glen;
                    const size_t km_pos = pos % km_glen;

                    const size_t pos_id_unitig = (km_id << 32) | MASK_CONTIG_TYPE | km_pos;

                    minz.push_back(km_unitigs.getMinimizer(km_id, km_pos).rep());
                }

                v_bmp_km[i].clear();
            }
        }
    }

    if (read_success) {

        size_t nb_special_minz = 0;

        in.read(reinterpret_cast<char*>(&nb_special_minz), sizeof(size_t));

        read_success = !in.fail();

        if (read_success) {

            std::vector<BitContainer> v_bmp_abundant((nb_special_minz >> 32) + 1);
            std::vector<BitContainer> v_bmp_overcrowded((nb_special_minz >> 32) + 1);

            for (size_t i = 0; (i < v_bmp_abundant.size()) && read_success; ++i) read_success = v_bmp_abundant[i].read(in);
            for (size_t i = 0; (i < v_bmp_overcrowded.size()) && read_success; ++i) read_success = v_bmp_overcrowded[i].read(in);

            for (size_t i = 0; (i < nb_special_minz) && read_success; ++i) {

                Minimizer minz_rep;

                const bool isAbundant = v_bmp_abundant[i >> 32].contains(i & 0x00000000ffffffffULL);
                const bool isOvercrowded = v_bmp_overcrowded[i >> 32].contains(i & 0x00000000ffffffffULL);

                size_t pos_id_unitig = 0;

                read_success = minz_rep.read(in);

                if (isAbundant || isOvercrowded) {

                    pos_id_unitig = MASK_CONTIG_ID | ((static_cast<size_t>(!isOvercrowded) - 1) & MASK_CONTIG_TYPE);

                    if (isAbundant && read_success) {

                        uint32_t count_abundant = 0;

                        in.read(reinterpret_cast<char*>(&count_abundant), sizeof(uint32_t));

                        pos_id_unitig |= (static_cast<size_t>(count_abundant));
                        read_success = !in.fail();
                    }

                    if (read_success) {

                        minz.push_back(minz_rep);
                    }
                }
                else {

                    in.read(reinterpret_cast<char*>(&pos_id_unitig), sizeof(size_t));

                    read_success = !in.fail();

                    if (read_success) {

                        minz.push_back(minz_rep);
                    }
                }
            }
        }
    }

    return read_success;
}

template<typename U, typename G>
bool CompactedDBG<U, G>::readBinaryIndex(std::istream& in, const uint64_t checksum,
                                         std::vector<Minimizer>& minz, uint32_t threads) {

    bool read_success = !in.fail();

    // 0 - Write file format version, checksum and number of minimizers
    if (read_success) {

        size_t file_format_version = 0, v_unitigs_sz = 0, km_unitigs_sz = 0, h_kmers_ccov_sz = 0, hmap_min_unitigs_sz = 0;
        uint64_t read_checksum = 0;

        read_success = readBinaryIndexHead(in, file_format_version, v_unitigs_sz, km_unitigs_sz, h_kmers_ccov_sz, hmap_min_unitigs_sz, read_checksum);

        if (!read_success || ((file_format_version >> 32) != BFG_METABIN_FORMAT_HEADER)) return false;
        if (!read_success || (read_checksum != checksum)) return false;

        if (minz.empty()) {
            hmap_min_unitigs = MinimizerIndex(hmap_min_unitigs_sz);
        } else {
            hmap_min_unitigs = MinimizerIndex(2);
            hmap_min_unitigs.generate_mphf(minz);
            minz.clear();
        }
    }

    if (read_success) {

        size_t nb_bmp_unitigs = 0;

        in.read(reinterpret_cast<char*>(&nb_bmp_unitigs), sizeof(size_t));

        std::vector<BitContainer> v_bmp_unitigs(nb_bmp_unitigs);

        for (size_t i = 0; read_success && (i < nb_bmp_unitigs); ++i) read_success = v_bmp_unitigs[i].read(in);

        if (read_success && !v_bmp_unitigs.empty() && !v_bmp_unitigs.front().isEmpty()) {

            size_t unitig_id = 0,
                   tot_unitig_len = 0,
                   curr_unitig_len = v_unitigs[0]->getSeq().size() - g_ + 1;

            for (size_t i = 0; i < nb_bmp_unitigs; ++i) {

                const size_t id_bmp = i << 32;

                for (const auto pos_bmp : v_bmp_unitigs[i]) {

                    const size_t pos = id_bmp + pos_bmp;

                    while (pos >= (tot_unitig_len + curr_unitig_len)) {

                        ++unitig_id;
                        tot_unitig_len += curr_unitig_len;

                        curr_unitig_len = v_unitigs[unitig_id]->getSeq().size() - g_ + 1;
                    }

                    const size_t relative_pos = pos - tot_unitig_len;
                    const size_t pos_id_unitig = (unitig_id << 32) | relative_pos;

                    const Minimizer minz_rep = v_unitigs[unitig_id]->getSeq().getMinimizer(relative_pos).rep();

                    std::pair<MinimizerIndex::iterator, bool> p = hmap_min_unitigs.insert(minz_rep, packed_tiny_vector(), 0);

                    packed_tiny_vector& v = p.first.getVector();
                    uint8_t& flag_v = p.first.getVectorSize();

                    flag_v = v.push_back(pos_id_unitig, flag_v);
                }

                v_bmp_unitigs[i].clear();
            }
        }
    }

    if (read_success) {

        size_t nb_bmp_km = 0;

        in.read(reinterpret_cast<char*>(&nb_bmp_km), sizeof(size_t));

        std::vector<BitContainer> v_bmp_km(nb_bmp_km);

        for (size_t i = 0; (i < nb_bmp_km) && read_success; ++i) read_success = v_bmp_km[i].read(in);

        if (read_success && !v_bmp_km.empty() && !v_bmp_km.front().isEmpty()) {

            const size_t km_glen = k_ - g_ + 1;

            for (size_t i = 0; i < nb_bmp_km; ++i) {

                const size_t id_bmp = i << 32;

                for (const auto pos_bmp : v_bmp_km[i]) {

                    const size_t pos = id_bmp + pos_bmp;

                    const size_t km_id = pos / km_glen;
                    const size_t km_pos = pos % km_glen;

                    const size_t pos_id_unitig = (km_id << 32) | MASK_CONTIG_TYPE | km_pos;

                    const Minimizer minz_rep = km_unitigs.getMinimizer(km_id, km_pos).rep();

                    std::pair<MinimizerIndex::iterator, bool> p = hmap_min_unitigs.insert(minz_rep, packed_tiny_vector(), 0);

                    packed_tiny_vector& v = p.first.getVector();
                    uint8_t& flag_v = p.first.getVectorSize();

                    flag_v = v.push_back(pos_id_unitig, flag_v);
                }

                v_bmp_km[i].clear();
            }
        }
    }

    if (read_success) {

        size_t nb_special_minz = 0;

        in.read(reinterpret_cast<char*>(&nb_special_minz), sizeof(size_t));

        read_success = !in.fail();

        if (read_success) {

            std::vector<BitContainer> v_bmp_abundant((nb_special_minz >> 32) + 1);
            std::vector<BitContainer> v_bmp_overcrowded((nb_special_minz >> 32) + 1);

            for (size_t i = 0; (i < v_bmp_abundant.size()) && read_success; ++i) read_success = v_bmp_abundant[i].read(in);
            for (size_t i = 0; (i < v_bmp_overcrowded.size()) && read_success; ++i) read_success = v_bmp_overcrowded[i].read(in);

            for (size_t i = 0; (i < nb_special_minz) && read_success; ++i) {

                Minimizer minz_rep;

                const bool isAbundant = v_bmp_abundant[i >> 32].contains(i & 0x00000000ffffffffULL);
                const bool isOvercrowded = v_bmp_overcrowded[i >> 32].contains(i & 0x00000000ffffffffULL);

                size_t pos_id_unitig = 0;

                read_success = minz_rep.read(in);

                if (isAbundant || isOvercrowded) {

                    pos_id_unitig = MASK_CONTIG_ID | ((static_cast<size_t>(!isOvercrowded) - 1) & MASK_CONTIG_TYPE);

                    if (isAbundant && read_success) {

                        uint32_t count_abundant = 0;

                        in.read(reinterpret_cast<char*>(&count_abundant), sizeof(uint32_t));

                        pos_id_unitig |= (static_cast<size_t>(count_abundant));
                        read_success = !in.fail();
                    }

                    if (read_success) {

                        std::pair<MinimizerIndex::iterator, bool> p = hmap_min_unitigs.insert(minz_rep, packed_tiny_vector(), 0);

                        packed_tiny_vector& v = p.first.getVector();
                        uint8_t& flag_v = p.first.getVectorSize();

                        flag_v = v.push_back(pos_id_unitig, flag_v);
                    }
                }
                else {

                    in.read(reinterpret_cast<char*>(&pos_id_unitig), sizeof(size_t));

                    read_success = !in.fail();

                    if (read_success) {

                        std::pair<MinimizerIndex::iterator, bool> p = hmap_min_unitigs.insert(minz_rep, packed_tiny_vector(), 0);

                        packed_tiny_vector& v = p.first.getVector();
                        uint8_t& flag = p.first.getVectorSize();
                        size_t v_sz = v.size(flag);

                        if ((v_sz) == 0 || ((v(v_sz-1, flag) & MASK_CONTIG_ID) != MASK_CONTIG_ID)) flag = v.push_back(pos_id_unitig, flag);
                        else flag = v.insert(pos_id_unitig, v_sz-1, flag);
                    }
                }
            }
        }
    }

    return read_success;
}


// BEGIN TODO
// Refactor this and above function so as to minimize repetition
template<typename U, typename G>
bool CompactedDBG<U, G>::readBinaryIndex(std::istream& in, const uint64_t checksum, boophf_t* mphf, uint32_t threads) {

    bool read_success = !in.fail();

    // 0 - Write file format version, checksum and number of minimizers
    if (read_success) {

        size_t file_format_version = 0, v_unitigs_sz = 0, km_unitigs_sz = 0, h_kmers_ccov_sz = 0, hmap_min_unitigs_sz = 0;
        uint64_t read_checksum = 0;

        read_success = readBinaryIndexHead(in, file_format_version, v_unitigs_sz, km_unitigs_sz, h_kmers_ccov_sz, hmap_min_unitigs_sz, read_checksum);

        if (!read_success || ((file_format_version >> 32) != BFG_METABIN_FORMAT_HEADER)) return false;
        if (!read_success || (read_checksum != checksum)) return false;

        hmap_min_unitigs = MinimizerIndex(2);
        hmap_min_unitigs.register_mphf(mphf);
    }

    if (read_success) {

        size_t nb_bmp_unitigs = 0;

        in.read(reinterpret_cast<char*>(&nb_bmp_unitigs), sizeof(size_t));

        std::vector<BitContainer> v_bmp_unitigs(nb_bmp_unitigs);

        for (size_t i = 0; read_success && (i < nb_bmp_unitigs); ++i) read_success = v_bmp_unitigs[i].read(in);

        if (read_success && !v_bmp_unitigs.empty() && !v_bmp_unitigs.front().isEmpty()) {

            std::vector<size_t> unitig_pos_ids;
            const size_t n_elems_buffer = 33554432; // 2^25
            hmap_min_unitigs.init_threads();

            size_t unitig_id = 0,
                   tot_unitig_len = 0,
                   curr_unitig_len = v_unitigs[0]->getSeq().size() - g_ + 1;

            for (size_t i = 0; i < nb_bmp_unitigs; ++i) {

                const size_t id_bmp = i << 32;

                const auto bmp_unitigs_size = v_bmp_unitigs[i].size();
                unitig_pos_ids.reserve(n_elems_buffer > bmp_unitigs_size ? bmp_unitigs_size : n_elems_buffer);
                std::vector<std::thread> workers;
                size_t n = 0;
                for (const auto pos_bmp : v_bmp_unitigs[i]) {

                    n++;
                    const size_t pos = id_bmp + pos_bmp;

                    while (pos >= (tot_unitig_len + curr_unitig_len)) {

                        ++unitig_id;
                        tot_unitig_len += curr_unitig_len;
                        curr_unitig_len = v_unitigs[unitig_id]->getSeq().size() - g_ + 1;
                    }
                    const size_t relative_pos = pos - tot_unitig_len;
                    const size_t pos_id_unitig = (unitig_id << 32) | relative_pos;

                    unitig_pos_ids.push_back(pos_id_unitig);

                    if (unitig_pos_ids.size() >= unitig_pos_ids.capacity() || n == bmp_unitigs_size) {
                        for (size_t t = 0; t < threads; t++) {
                            workers.emplace_back([&, t]{
                                for (size_t j = t; j < unitig_pos_ids.size(); j+= threads) {
                                    const size_t relative_pos = unitig_pos_ids[j] & 0xFFFFFFFF;
                                    const size_t unitig_id_ = unitig_pos_ids[j] >> 32;
                                    const Minimizer minz_rep = v_unitigs[unitig_id_]->getSeq().getMinimizer(relative_pos).rep();
                                    hmap_min_unitigs.add_unitig_p(minz_rep, unitig_pos_ids[j]);
                                }
                           });
                        }
                        for (auto& t : workers) t.join();
                        unitig_pos_ids.clear();
                        workers.clear();
                    }
                }

            }

            hmap_min_unitigs.release_threads();
        }
    }

    if (read_success) {

        size_t nb_bmp_km = 0;

        in.read(reinterpret_cast<char*>(&nb_bmp_km), sizeof(size_t));

        std::vector<BitContainer> v_bmp_km(nb_bmp_km);

        for (size_t i = 0; (i < nb_bmp_km) && read_success; ++i) read_success = v_bmp_km[i].read(in);

        if (read_success && !v_bmp_km.empty() && !v_bmp_km.front().isEmpty()) {

            const size_t km_glen = k_ - g_ + 1;

            for (size_t i = 0; i < nb_bmp_km; ++i) {

                const size_t id_bmp = i << 32;

                for (const auto pos_bmp : v_bmp_km[i]) {

                    const size_t pos = id_bmp + pos_bmp;

                    const size_t km_id = pos / km_glen;
                    const size_t km_pos = pos % km_glen;

                    const size_t pos_id_unitig = (km_id << 32) | MASK_CONTIG_TYPE | km_pos;

                    const Minimizer minz_rep = km_unitigs.getMinimizer(km_id, km_pos).rep();

                    std::pair<MinimizerIndex::iterator, bool> p = hmap_min_unitigs.insert(minz_rep, packed_tiny_vector(), 0);

                    packed_tiny_vector& v = p.first.getVector();
                    uint8_t& flag_v = p.first.getVectorSize();

                    flag_v = v.push_back(pos_id_unitig, flag_v);
                }

                v_bmp_km[i].clear();
            }
        }
    }

    if (read_success) {

        size_t nb_special_minz = 0;

        in.read(reinterpret_cast<char*>(&nb_special_minz), sizeof(size_t));

        read_success = !in.fail();

        if (read_success) {

            std::vector<BitContainer> v_bmp_abundant((nb_special_minz >> 32) + 1);
            std::vector<BitContainer> v_bmp_overcrowded((nb_special_minz >> 32) + 1);

            for (size_t i = 0; (i < v_bmp_abundant.size()) && read_success; ++i) read_success = v_bmp_abundant[i].read(in);
            for (size_t i = 0; (i < v_bmp_overcrowded.size()) && read_success; ++i) read_success = v_bmp_overcrowded[i].read(in);

            for (size_t i = 0; (i < nb_special_minz) && read_success; ++i) {

                Minimizer minz_rep;

                const bool isAbundant = v_bmp_abundant[i >> 32].contains(i & 0x00000000ffffffffULL);
                const bool isOvercrowded = v_bmp_overcrowded[i >> 32].contains(i & 0x00000000ffffffffULL);

                size_t pos_id_unitig = 0;

                read_success = minz_rep.read(in);

                if (isAbundant || isOvercrowded) {

                    pos_id_unitig = MASK_CONTIG_ID | ((static_cast<size_t>(!isOvercrowded) - 1) & MASK_CONTIG_TYPE);

                    if (isAbundant && read_success) {

                        uint32_t count_abundant = 0;

                        in.read(reinterpret_cast<char*>(&count_abundant), sizeof(uint32_t));

                        pos_id_unitig |= (static_cast<size_t>(count_abundant));
                        read_success = !in.fail();
                    }

                    if (read_success) {

                        std::pair<MinimizerIndex::iterator, bool> p = hmap_min_unitigs.insert(minz_rep, packed_tiny_vector(), 0);

                        packed_tiny_vector& v = p.first.getVector();
                        uint8_t& flag_v = p.first.getVectorSize();

                        flag_v = v.push_back(pos_id_unitig, flag_v);
                    }
                }
                else {

                    in.read(reinterpret_cast<char*>(&pos_id_unitig), sizeof(size_t));

                    read_success = !in.fail();

                    if (read_success) {

                        std::pair<MinimizerIndex::iterator, bool> p = hmap_min_unitigs.insert(minz_rep, packed_tiny_vector(), 0);

                        packed_tiny_vector& v = p.first.getVector();
                        uint8_t& flag = p.first.getVectorSize();
                        size_t v_sz = v.size(flag);

                        if ((v_sz) == 0 || ((v(v_sz-1, flag) & MASK_CONTIG_ID) != MASK_CONTIG_ID)) flag = v.push_back(pos_id_unitig, flag);
                        else flag = v.insert(pos_id_unitig, v_sz-1, flag);
                    }
                }
            }
        }
    }

    return read_success;
}
// END TODO



template<typename U, typename G>
std::pair<uint64_t, bool> CompactedDBG<U, G>::readBinaryGraph(const std::string& fn) {

    if ((fn.length() == 0) || !check_file_exists(fn)) return {0, false};

    std::ifstream infile;
    std::istream in(0);

    infile.open(fn.c_str());
    in.rdbuf(infile.rdbuf());

    return readBinaryGraph(in);
}

template<typename U, typename G>
std::pair<uint64_t, bool> CompactedDBG<U, G>::readBinaryGraph(std::istream& in) {

    bool read_success = !in.fail();

    uint64_t graph_checksum = 0;

    clear();

    if (read_success) {

        size_t file_format_version = 0;
        int rk = 0, rg = 0;

        in.read(reinterpret_cast<char*>(&file_format_version), sizeof(size_t));
        in.read(reinterpret_cast<char*>(&rk), sizeof(int));
        in.read(reinterpret_cast<char*>(&rg), sizeof(int));

        read_success = (read_success && !in.fail());

        if (read_success) {

            const size_t k = static_cast<size_t>(rk);
            const size_t g = static_cast<size_t>(rg);

            graph_checksum = wyhash(&k, sizeof(size_t), 0, _wyp);
            graph_checksum = wyhash(&g, sizeof(size_t), graph_checksum, _wyp);
        }

        if ((file_format_version >> 32) != BFG_GRAPHBIN_FORMAT_HEADER) return {graph_checksum, false};
        if (read_success) *this = CompactedDBG<U, G>(rk, rg);
        if (invalid) return {graph_checksum, false};
    }

    if (read_success) {

        size_t v_unitigs_sz = 0;

        in.read(reinterpret_cast<char*>(&v_unitigs_sz), sizeof(size_t));

        read_success = !in.fail();

        if (read_success) {

            v_unitigs.reserve(v_unitigs_sz);

            for (size_t i = 0; (i < v_unitigs_sz) && read_success; ++i) {

                CompressedSequence cs;
                CompressedCoverage cc;

                Unitig<U>* unitig;

                read_success = cs.read(in);
                graph_checksum = cs.hash(graph_checksum);
                cc = CompressedCoverage(cs.size() - k_ + 1, false);
                unitig = new Unitig<U>(std::move(cs), std::move(cc));

                v_unitigs.push_back(unitig);
            }
        }
    }

    if (read_success) {

        read_success = km_unitigs.read(in);

        for (size_t i = 0; i < km_unitigs.size(); ++i) graph_checksum = km_unitigs.getKmer(i).hash(graph_checksum);
    }

    if (read_success) {

        const CompressedCoverage cc(1, false);

        size_t h_kmers_ccov_sz = 0;

        in.read(reinterpret_cast<char*>(&h_kmers_ccov_sz), sizeof(size_t));

        read_success = !in.fail();

        if (read_success) {

            h_kmers_ccov.reserve(h_kmers_ccov_sz);

            for (size_t i = 0; read_success && (i < h_kmers_ccov_sz); ++i) {

                Kmer km;

                read_success = km.read(in);
                graph_checksum = km.hash(graph_checksum);

                h_kmers_ccov.insert(km, cc);
            }
        }
    }

    if (read_success) {

        setFullCoverage(1);

        for (auto& unitig : *this) unitig.setFullCoverage();
    }

    return {graph_checksum, read_success};
}

template<typename U, typename G>
bool CompactedDBG<U, G>::writeBinary(const std::string& fn, const size_t nb_threads) const {

    if (fn.length() == 0) return false;

    std::ofstream outfile;
    std::ostream out(0);

    outfile.open(fn.c_str());
    out.rdbuf(outfile.rdbuf());

    return writeBinary(out, nb_threads);
}

template<typename U, typename G>
bool CompactedDBG<U, G>::writeBinary(std::ostream& out, const size_t nb_threads) const {

    if (!out.fail()) {

        const bool write_success = writeBinaryGraph(out, nb_threads);

        if (write_success) return writeBinaryIndex(out, checksum(), nb_threads);
    }

    return false;
}

template<typename U, typename G>
bool CompactedDBG<U, G>::writeBinaryGraph(const std::string& fn, const size_t nb_threads) const {

    if (fn.length() == 0) return false;

    std::ofstream outfile;
    std::ostream out(0);

    outfile.open(fn.c_str());
    out.rdbuf(outfile.rdbuf());

    return writeBinaryGraph(out, nb_threads);
}

template<typename U, typename G>
bool CompactedDBG<U, G>::writeBinaryGraph(std::ostream& out, const size_t nb_threads) const {

    bool write_success = !out.fail();

    // 0- Write file format version, k-mer and g-mer lengths
    if (write_success) {

        const size_t fileformat_version = (static_cast<size_t>(BFG_GRAPHBIN_FORMAT_HEADER) << 32) | static_cast<size_t>(BFG_GRAPHBIN_FORMAT_VERSION);

        out.write(reinterpret_cast<const char*>(&fileformat_version), sizeof(size_t));
        out.write(reinterpret_cast<const char*>(&k_), sizeof(int));
        out.write(reinterpret_cast<const char*>(&g_), sizeof(int));

        write_success = (write_success && !out.fail());
    }

    // 1 - Write unitigs longer than k
    if (write_success) {

        const size_t v_unitigs_sz = v_unitigs.size();

        out.write(reinterpret_cast<const char*>(&v_unitigs_sz), sizeof(size_t));

        write_success = !out.fail();

        for (size_t i = 0; (i < v_unitigs_sz) && write_success; ++i) write_success = v_unitigs[i]->getSeq().write(out);
    }

    // 2 - Write unitigs with length k for which minimizers are not over abundant
    if (write_success) write_success = km_unitigs.write(out);

    // 3 - Write unitigs with length k for which minimizers are over abundant
    if (write_success) {

        const size_t h_kmers_ccov_sz = h_kmers_ccov.size();

        out.write(reinterpret_cast<const char*>(&h_kmers_ccov_sz), sizeof(size_t));

        write_success = !out.fail();

        for (typename h_kmers_ccov_t::const_iterator it = h_kmers_ccov.begin(); (it != h_kmers_ccov.end()) && write_success; ++it) write_success = it.getKey().write(out);
    }

    return (write_success && !out.fail());
}

template<typename U, typename G>
bool CompactedDBG<U, G>::writeBinaryIndex(const std::string& fn, const uint64_t checksum, const size_t nb_threads) const {

    if (fn.length() == 0) return false;

    std::ofstream outfile;
    std::ostream out(0);

    outfile.open(fn.c_str());
    out.rdbuf(outfile.rdbuf());

    return writeBinaryIndex(out, checksum, nb_threads);
}

template<typename U, typename G>
bool CompactedDBG<U, G>::writeBinaryIndex(std::ostream& out, const uint64_t checksum, const size_t nb_threads) const {

    bool write_success = !out.fail();

    const size_t nb_bmp_minz = (hmap_min_unitigs.size() >> 32) + 1;

    std::atomic<size_t> nb_special_minz;

    std::vector<BitContainer> v_bmp_minz(nb_bmp_minz);

    nb_special_minz = 0;

    // 0 - Write file format version, checksum and number of minimizers
    if (write_success) {

        const size_t fileformat_version = (static_cast<size_t>(BFG_METABIN_FORMAT_HEADER) << 32) | static_cast<size_t>(BFG_METABIN_FORMAT_VERSION);

        const size_t v_unitigs_sz = v_unitigs.size();
        const size_t km_unitigs_sz = km_unitigs.size();
        const size_t h_kmers_ccov_sz = h_kmers_ccov.size();
        const size_t hmap_min_unitigs_sz = hmap_min_unitigs.size();

        out.write(reinterpret_cast<const char*>(&fileformat_version), sizeof(size_t)); // Write header for binary index file, including file format version
        out.write(reinterpret_cast<const char*>(&checksum), sizeof(uint64_t)); // Write graph checksum

        out.write(reinterpret_cast<const char*>(&v_unitigs_sz), sizeof(size_t)); // Write number of unitigs with length > k
        out.write(reinterpret_cast<const char*>(&km_unitigs_sz), sizeof(size_t)); // Write number of unitigs with length == k and non-abundant minimizer
        out.write(reinterpret_cast<const char*>(&h_kmers_ccov_sz), sizeof(size_t)); // Write number of unitigs with length == k and abundant minimizer

        out.write(reinterpret_cast<const char*>(&hmap_min_unitigs_sz), sizeof(size_t)); // Write number of minimizers

        write_success = (write_success && !out.fail());
    }

    if (write_success) {

        const size_t v_unitigs_sz = v_unitigs.size();
        const size_t nb_block_unitigs = std::max(static_cast<size_t>((v_unitigs_sz + 15) / 16) + 1, static_cast<size_t>(1));

        std::vector<size_t> v_block_len_unitigs(nb_block_unitigs, 0);

        for (size_t i = 0; i < v_unitigs_sz; ++i) v_block_len_unitigs[(i >> 4) + 1] += v_unitigs[i]->getSeq().size() - g_ + 1;
        for (size_t i = 1; i < nb_block_unitigs; ++i) v_block_len_unitigs[i] += v_block_len_unitigs[i-1];

        const size_t nb_bmp_unitigs = (v_block_len_unitigs.back() >> 32) + 1;
        const size_t nb_bmp_km_short = ((km_unitigs.size() * (k_ - g_ + 1)) >> 32) + 1;

        std::vector<BitContainer> v_bmp_unitigs(nb_bmp_unitigs), v_bmp_km_short(nb_bmp_km_short);
        std::vector<SpinLock> s_bmp_unitigs(nb_bmp_unitigs), s_bmp_km_short(nb_bmp_km_short), s_bmp_minz(nb_bmp_minz);

        auto compactMinimizers = [&](MinimizerIndex::const_iterator it, MinimizerIndex::const_iterator ite, size_t id_minz) {

            std::vector<BitContainer> lv_bmp_unitigs(nb_bmp_unitigs), lv_bmp_km_short(nb_bmp_km_short), lv_bmp_minz(nb_bmp_minz);

            size_t l_nb_special_minz = 0;

            while (it != ite) { // Annotate in bitmap the position of every minimizer in the unitigs

                const packed_tiny_vector& v = it.getVector();
                const uint8_t flag_v = it.getVectorSize();
                const int v_sz = v.size(flag_v);

                const Minimizer& minz_key = it.getKey();

                for (size_t i = 0; i < v_sz; ++i){

                    const size_t unitig_idx = v(i, flag_v);
                    const size_t unitig_id = unitig_idx >> 32;
                    const size_t unitig_pos = unitig_idx & MASK_CONTIG_POS;

                    const bool isShort = static_cast<bool>(unitig_idx & MASK_CONTIG_TYPE);

                    if (unitig_id == RESERVED_ID) {

                        lv_bmp_minz[id_minz >> 32].add(id_minz & 0x00000000ffffffffULL);
                        ++l_nb_special_minz;
                    }
                    else if (isShort) {

                        const Minimizer minz = km_unitigs.getMinimizer(unitig_id, unitig_pos).rep();

                        if (minz == minz_key) {

                            const size_t pos = (unitig_id * (k_ - g_ + 1)) + unitig_pos;

                            lv_bmp_km_short[pos >> 32].add(pos & 0x00000000ffffffffULL);
                        }
                        else {

                            lv_bmp_minz[id_minz >> 32].add(id_minz & 0x00000000ffffffffULL);

                            ++l_nb_special_minz;
                        }
                    }
                    else {

                        const Minimizer minz = v_unitigs[unitig_id]->getSeq().getMinimizer(unitig_pos).rep();

                        if (minz == minz_key) {

                            size_t pos = v_block_len_unitigs[unitig_id >> 4] + unitig_pos;

                            for (size_t j = (unitig_id & 0xfffffffffffffff0ULL); j < unitig_id; ++j) pos += v_unitigs[j]->getSeq().size() - g_ + 1;

                            lv_bmp_unitigs[pos >> 32].add(pos & 0x00000000ffffffffULL);
                        }
                        else {

                            lv_bmp_minz[id_minz >> 32].add(id_minz & 0x00000000ffffffffULL);

                            ++l_nb_special_minz;
                        }
                    }
                }

                ++it;
                ++id_minz;
            }

            {
                nb_special_minz += l_nb_special_minz;

                for (size_t i = 0; i < nb_bmp_unitigs; ++i) {

                    if (!lv_bmp_unitigs[i].isEmpty()) {

                        s_bmp_unitigs[i].acquire();

                        v_bmp_unitigs[i] |= lv_bmp_unitigs[i];

                        s_bmp_unitigs[i].release();
                    }
                }

                for (size_t i = 0; i < nb_bmp_km_short; ++i) {

                    if (!lv_bmp_km_short[i].isEmpty()) {

                        s_bmp_km_short[i].acquire();

                        v_bmp_km_short[i] |= lv_bmp_km_short[i];

                        s_bmp_km_short[i].release();
                    }
                }

                for (size_t i = 0; i < nb_bmp_minz; ++i) {

                    if (!lv_bmp_minz[i].isEmpty()) {

                        s_bmp_minz[i].acquire();

                        v_bmp_minz[i] |= lv_bmp_minz[i];

                        s_bmp_minz[i].release();
                    }
                }
            }
        };

        {
            size_t id_minz = 0;

            MinimizerIndex::const_iterator it = hmap_min_unitigs.begin();
            MinimizerIndex::const_iterator ite = hmap_min_unitigs.end();

            if (nb_threads == 1) {

                while (it != ite) {

                    MinimizerIndex::const_iterator lit = it;
                    MinimizerIndex::const_iterator lite = it;

                    size_t lid_minz = id_minz;

                    for (size_t i = 0; (i < 65536) && (lite != ite); ++i, ++id_minz) ++lite;

                    it = lite;

                    compactMinimizers(lit, lite, lid_minz);
                }
            }
            else {

                std::vector<std::thread> workers; // need to keep track of threads so we can join them

                std::mutex mutex_minz_idx;

                for (size_t t = 0; t < nb_threads; ++t){

                    workers.emplace_back(

                        [&]{

                            MinimizerIndex::const_iterator lit;
                            MinimizerIndex::const_iterator lite;

                            size_t lid_minz;

                            while (true) {

                                {
                                    std::unique_lock<std::mutex> lock(mutex_minz_idx);

                                    if (it == ite) return;

                                    lit = it;
                                    lite = it;
                                    lid_minz = id_minz;

                                    for (size_t i = 0; (i < 65536) && (lite != ite); ++i, ++id_minz) ++lite;

                                    it = lite;
                                }

                                compactMinimizers(lit, lite, lid_minz);
                            }
                        }
                    );
                }

                for (auto& t : workers) t.join();
            }
        }

        {
            for (auto& bmp : v_bmp_unitigs) bmp.runOptimize();

            out.write(reinterpret_cast<const char*>(&nb_bmp_unitigs), sizeof(size_t));

            write_success = !out.fail();

            for (size_t i = 0; (i < nb_bmp_unitigs) && write_success; ++i) write_success = v_bmp_unitigs[i].write(out);

            v_bmp_unitigs.clear();
        }

        {
            for (auto& bmp : v_bmp_km_short) bmp.runOptimize();

            out.write(reinterpret_cast<const char*>(&nb_bmp_km_short), sizeof(size_t));

            write_success = !out.fail();

            for (size_t i = 0; (i < nb_bmp_km_short) && write_success; ++i) write_success = v_bmp_km_short[i].write(out);

            v_bmp_km_short.clear();
        }

        for (auto& bmp : v_bmp_minz) bmp.runOptimize();
    }

    if (write_success) {

        const size_t nb_bmp_special_minz = (nb_special_minz >> 32) + 1;

        size_t id_special = 0;

        std::vector<BitContainer> v_bmp_abundant(nb_bmp_special_minz), v_bmp_overcrowded(nb_bmp_special_minz);

        {
            size_t id_minz = 0;

            MinimizerIndex::const_iterator it = hmap_min_unitigs.begin();
            MinimizerIndex::const_iterator ite = hmap_min_unitigs.end();

            while (it != ite) { // Annotate in bitmap the position of every minimizer in the unitigs

                if (v_bmp_minz[id_minz >> 32].contains(id_minz & 0x00000000ffffffffULL)) {

                    const packed_tiny_vector& v = it.getVector();
                    const uint8_t flag_v = it.getVectorSize();
                    const int v_sz = v.size(flag_v);

                    const Minimizer& minz_key = it.getKey();

                    for (size_t i = 0; (i < v_sz); ++i){

                        const size_t unitig_idx = v(i, flag_v);
                        const size_t unitig_id = unitig_idx >> 32;
                        const size_t unitig_pos = unitig_idx & MASK_CONTIG_POS;

                        const bool isShort = static_cast<bool>(unitig_idx & MASK_CONTIG_TYPE);

                        if (unitig_id == RESERVED_ID) {

                            const size_t id_bmp = id_special >> 32;
                            const size_t pos_bmp = id_special & 0x00000000ffffffffULL;

                            if (unitig_pos != 0) v_bmp_abundant[id_bmp].add(pos_bmp);
                            if (isShort) v_bmp_overcrowded[id_bmp].add(pos_bmp);

                               ++id_special;
                        }
                        else {

                            Minimizer minz;

                            if (isShort) minz = km_unitigs.getMinimizer(unitig_id, unitig_pos).rep();
                            else minz = v_unitigs[unitig_id]->getSeq().getMinimizer(unitig_pos).rep();

                            id_special += static_cast<size_t>(minz != minz_key);
                        }
                    }
                }

                ++it;
                ++id_minz;
            }

            for (auto& bmp : v_bmp_abundant) bmp.runOptimize();
            for (auto& bmp : v_bmp_overcrowded) bmp.runOptimize();
        }

        {
            out.write(reinterpret_cast<const char*>(&nb_special_minz), sizeof(size_t)); // Pre-reserve space

            write_success = !out.fail();

            for (size_t i = 0; (i < v_bmp_abundant.size()) && write_success; ++i) write_success = v_bmp_abundant[i].write(out);
            for (size_t i = 0; (i < v_bmp_overcrowded.size()) && write_success; ++i) write_success = v_bmp_overcrowded[i].write(out);

            v_bmp_abundant.clear();
            v_bmp_overcrowded.clear();
        }

        if (write_success) {

            size_t id_minz = 0;

            MinimizerIndex::const_iterator it = hmap_min_unitigs.begin();
            MinimizerIndex::const_iterator ite = hmap_min_unitigs.end();

            while ((it != ite) && write_success) { // Annotate in bitmap the position of every minimizer in the unitigs

                if (v_bmp_minz[id_minz >> 32].contains(id_minz & 0x00000000ffffffffULL)) {

                    const packed_tiny_vector& v = it.getVector();
                    const uint8_t flag_v = it.getVectorSize();
                    const int v_sz = v.size(flag_v);

                    const Minimizer& minz_key = it.getKey();

                    for (size_t i = 0; (i < v_sz) && write_success; ++i){

                        const size_t unitig_idx = v(i, flag_v);
                        const size_t unitig_id = unitig_idx >> 32;
                        const size_t unitig_pos = unitig_idx & MASK_CONTIG_POS;

                        if (unitig_id == RESERVED_ID) {

                            write_success = minz_key.write(out);

                            if (write_success && (unitig_pos != 0)) {

                                out.write(reinterpret_cast<const char*>(&unitig_pos), sizeof(uint32_t)); // Pre-reserve space

                                write_success = !out.fail();
                            }
                        }
                        else {

                            const bool isShort = static_cast<bool>(unitig_idx & MASK_CONTIG_TYPE);

                            Minimizer minz;

                            if (isShort) minz = km_unitigs.getMinimizer(unitig_id, unitig_pos).rep();
                            else minz = v_unitigs[unitig_id]->getSeq().getMinimizer(unitig_pos).rep();

                            if (minz != minz_key) {

                                minz_key.write(out);

                                out.write(reinterpret_cast<const char*>(&unitig_idx), sizeof(size_t)); // Pre-reserve space

                                write_success = !out.fail();
                            }
                        }
                    }
                }

                ++it;
                ++id_minz;
            }
        }
    }

    return (write_success && !out.fail());
}

template<typename U, typename G>
void CompactedDBG<U, G>::makeGraphFromGFA(const std::string& fn, const size_t nb_threads) {

    size_t graph_file_id = 0;

    bool new_file_opened = false;

    GFA_Parser graph(fn);

    graph.open_read();

    GFA_Parser::GFA_line r = graph.read(graph_file_id, new_file_opened, true);

    if (nb_threads == 1){

        while ((r.first != nullptr) || (r.second != nullptr)){

            if (r.first != nullptr) addUnitig(r.first->seq, (r.first->seq.length() == k_) ? km_unitigs.size() : v_unitigs.size());

            r = graph.read(graph_file_id, new_file_opened, true);
        }
    }
    else {

        const size_t block_sz = 1024;

        std::atomic<size_t> v_kmers_sz;
        std::atomic<size_t> v_unitigs_sz;

        bool is_first = true;
        bool stop = false;

        SpinLock lck_unitig, lck_kmer;

        std::vector<std::thread> workers; // need to keep track of threads so we can join them

        std::mutex mutex_file;

        v_kmers_sz = 0;
        v_unitigs_sz = 0;

        hmap_min_unitigs.init_threads();

        for (size_t t = 0; t < nb_threads; ++t){

            workers.emplace_back(

                [&]{

                    std::vector<std::string> seq;

                    while (true) {

                        {
                            std::unique_lock<std::mutex> lock(mutex_file);

                            if (stop) return;

                            seq.clear();

                            for (size_t i = 0; (i < block_sz) && !stop; ++i){

                                if (!is_first) r = graph.read(graph_file_id, new_file_opened, true);
                                if (r.first != nullptr) seq.push_back(r.first->seq);

                                stop = ((r.first == nullptr) && (r.second == nullptr));
                                is_first = false;
                            }
                        }

                        for (const auto& s : seq) addUnitig(s, (s.length() == k_) ? v_kmers_sz++ : v_unitigs_sz++, lck_unitig, lck_kmer);
                    }
                }
            );
        }

        for (auto& t : workers) t.join();

        hmap_min_unitigs.release_threads();

        moveToAbundant();
    }
}

template<typename U, typename G>
void CompactedDBG<U, G>::makeGraphFromFASTA(const std::string& fn, const size_t nb_threads) {

    size_t graph_file_id = 0;

    FastqFile ff(std::vector<std::string>(1, fn));

    std::string seq;

    if (nb_threads == 1){

        while (ff.read_next(seq, graph_file_id) != -1) addUnitig(seq, (seq.length() == k_) ? km_unitigs.size() : v_unitigs.size());
    }
    else {

        const size_t block_sz = 1024;

        bool stop = false;

        std::atomic<size_t> v_kmers_sz;
        std::atomic<size_t> v_unitigs_sz;

        SpinLock lck_unitig, lck_kmer;

        std::vector<std::thread> workers; // need to keep track of threads so we can join them

        std::mutex mutex_file;

        v_kmers_sz = 0;
        v_unitigs_sz = 0;

        hmap_min_unitigs.init_threads();

        for (size_t t = 0; t < nb_threads; ++t){

            workers.emplace_back(

                [&]{

                    std::vector<std::string> v_seq;

                    while (true) {

                        {
                            std::unique_lock<std::mutex> lock(mutex_file);

                            if (stop) return;

                            v_seq.clear();

                            for (size_t i = 0; (i < block_sz) && !stop; ++i){

                                stop = (ff.read_next(seq, graph_file_id) == -1);

                                if (!stop && !seq.empty()) v_seq.push_back(seq);
                            }
                        }

                        for (const auto& s : v_seq) addUnitig(s, (s.length() == k_) ? v_kmers_sz++ : v_unitigs_sz++, lck_unitig, lck_kmer);
                    }
                }
            );
        }

        for (auto& t : workers) t.join();

        hmap_min_unitigs.release_threads();

        moveToAbundant();
    }
}

template<typename U, typename G>
std::pair<uint64_t, bool> CompactedDBG<U, G>::readGraphFromIndexFASTA(const std::string& graph_fn, const std::string& index_fn, const size_t k, const size_t g) {

    FastqFile ff(std::vector<std::string>(1, graph_fn));

    std::string seq;

    size_t file_format_version = 0, v_unitigs_sz = 0, km_unitigs_sz = 0, h_kmers_ccov_sz = 0, hmap_min_unitigs_sz = 0, graph_file_id = 0;
    uint64_t read_checksum = 0, graph_checksum = 0;

    bool read_success = readBinaryIndexHead(index_fn, file_format_version, v_unitigs_sz, km_unitigs_sz, h_kmers_ccov_sz, hmap_min_unitigs_sz, read_checksum);

    {
        graph_checksum = wyhash(&k, sizeof(size_t), 0, _wyp);
        graph_checksum = wyhash(&g, sizeof(size_t), graph_checksum, _wyp);
    }

    // 1 - Read unitigs with length > k
    if (read_success) {

        size_t i = 0;

        v_unitigs.reserve(v_unitigs_sz);

        while (read_success && (i < v_unitigs_sz) && (ff.read_next(seq, graph_file_id) != -1)) {

            if (seq.length() > k_) {

                CompressedSequence cs(seq);
                CompressedCoverage cc(seq.length() - k_ + 1, false);

                Unitig<U>* unitig;

                graph_checksum = cs.hash(graph_checksum);
                unitig = new Unitig<U>(std::move(cs), std::move(cc));

                v_unitigs.push_back(unitig);
            }
            else read_success = false;

            ++i;
        }

        read_success = (read_success && (i == v_unitigs_sz));
    }

    // 2 - Read unitigs with length == k which do not have abundant or over-crowded minimizers
    if (read_success) {

        size_t i = 0;

        km_unitigs.resize(km_unitigs_sz);

        while (read_success && (i < km_unitigs_sz) && (ff.read_next(seq, graph_file_id) != -1)) {

            if (seq.length() == k_) {

                const Kmer km(seq.c_str());

                read_success = km_unitigs.set(i, km);
                graph_checksum = km.hash(graph_checksum);
            }
            else read_success = false;

            ++i;
        }

        read_success = (read_success && (i == km_unitigs_sz));
    }

    // 3 - Read unitigs with length == k which do not have abundant or over-crowded minimizers
    if (read_success) {

        const CompressedCoverage cc(1, false);

        size_t i = 0;

        h_kmers_ccov.reserve(h_kmers_ccov_sz);

        while (read_success && (i < h_kmers_ccov_sz) && (ff.read_next(seq, graph_file_id) != -1)) {

            if (seq.length() == k_) {

                const Kmer km(seq.c_str());

                graph_checksum = km.hash(graph_checksum);
                h_kmers_ccov.insert(km, cc);
            }
            else read_success = false;

            ++i;
        }

        read_success = (read_success && (i == h_kmers_ccov_sz));
    }

    // 4 - If fasta contains more sequences than what index file is saying, fail reading
    read_success = (read_success && (ff.read_next(seq, graph_file_id) == -1));

    return {graph_checksum, read_success};
}

template<typename U, typename G>
std::pair<uint64_t, bool> CompactedDBG<U, G>::readGraphFromIndexGFA(const std::string& graph_fn, const std::string& index_fn, const size_t k, const size_t g) {

    bool new_file_opened = false;

    GFA_Parser graph(graph_fn);

    GFA_Parser::GFA_line r;

    size_t file_format_version = 0, v_unitigs_sz = 0, km_unitigs_sz = 0, h_kmers_ccov_sz = 0, hmap_min_unitigs_sz = 0, graph_file_id = 0;
    uint64_t read_checksum = 0, graph_checksum = 0;

    bool read_success = readBinaryIndexHead(index_fn, file_format_version, v_unitigs_sz, km_unitigs_sz, h_kmers_ccov_sz, hmap_min_unitigs_sz, read_checksum);

    if (read_success) {

        graph_checksum = wyhash(&k, sizeof(size_t), 0, _wyp);
        graph_checksum = wyhash(&g, sizeof(size_t), graph_checksum, _wyp);

        graph.open_read();

        r = graph.read(graph_file_id, new_file_opened, true);
    }

    // 1 - Read unitigs with length > k
    if (read_success) {

        size_t i = 0;

        v_unitigs.reserve(v_unitigs_sz);

        while (read_success && (i < v_unitigs_sz) && ((r.first != nullptr) || (r.second != nullptr))) {

            if (r.first != nullptr) {

                if (r.first->seq.length() > k_) {

                    CompressedSequence cs(r.first->seq);
                    CompressedCoverage cc(r.first->seq.length() - k_ + 1, false);

                    Unitig<U>* unitig;

                    graph_checksum = cs.hash(graph_checksum);
                    unitig = new Unitig<U>(std::move(cs), std::move(cc));

                    v_unitigs.push_back(unitig);
                }
                else read_success = false;

                ++i;
            }

            r = graph.read(graph_file_id, new_file_opened, true);
        }

        read_success = (read_success && (i == v_unitigs_sz));
    }

    // 2 - Read unitigs with length == k which do not have abundant or over-crowded minimizers
    if (read_success) {

        size_t i = 0;

        km_unitigs.resize(km_unitigs_sz);

        while (read_success && (i < km_unitigs_sz) && ((r.first != nullptr) || (r.second != nullptr))) {

            if (r.first != nullptr) {

                if (r.first->seq.length() == k_) {

                    const Kmer km(r.first->seq.c_str());

                    graph_checksum = km.hash(graph_checksum);
                    read_success = km_unitigs.set(i, km);
                }
                else read_success = false;

                ++i;
            }

            r = graph.read(graph_file_id, new_file_opened, true);
        }

        read_success = (read_success && (i == km_unitigs_sz));
    }

    // 3 - Read unitigs with length == k which do not have abundant or over-crowded minimizers
    if (read_success) {

        const CompressedCoverage cc(1, false);

        size_t i = 0;

        h_kmers_ccov.reserve(h_kmers_ccov_sz);

        while (read_success && (i < h_kmers_ccov_sz) && ((r.first != nullptr) || (r.second != nullptr))) {

            if (r.first != nullptr) {

                if (r.first->seq.length() == k_) {

                    const Kmer km(r.first->seq.c_str());

                    graph_checksum = km.hash(graph_checksum);
                    h_kmers_ccov.insert(km, cc);
                }
                else read_success = false;

                ++i;
            }

            r = graph.read(graph_file_id, new_file_opened, true);
        }

        read_success = (read_success && (i == h_kmers_ccov_sz));
    }

    // 4 - If gfa contains more sequences than what index file is saying, fail reading
    while ((r.first != nullptr) || (r.second != nullptr)) {

        if (r.first != nullptr) {

            read_success = false;
            break;
        }

        r = graph.read(graph_file_id, new_file_opened, true);
    }

    return {graph_checksum, read_success};
}

#endif
