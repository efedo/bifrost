#ifndef BIFROST_COLOREDCDBG_TCC
#define BIFROST_COLOREDCDBG_TCC

template<typename U>
ColoredCDBG<U>::ColoredCDBG(int kmer_length, int minimizer_length) : CompactedDBG<DataAccessor<U>, DataStorage<U>>(kmer_length, minimizer_length){

    invalid = this->isInvalid();
}

template<typename U>
ColoredCDBG<U>::ColoredCDBG(const ColoredCDBG& o) : CompactedDBG<DataAccessor<U>, DataStorage<U>>(o), invalid(o.invalid) {}

template<typename U>
ColoredCDBG<U>::ColoredCDBG(ColoredCDBG&& o) :  CompactedDBG<DataAccessor<U>, DataStorage<U>>(std::move(o)), invalid(o.invalid) {}

template<typename U>
void ColoredCDBG<U>::clear(){

    invalid = true;

    this->getData()->clear();
    CompactedDBG<DataAccessor<U>, DataStorage<U>>::clear();
}

template<typename U>
ColoredCDBG<U>& ColoredCDBG<U>::operator=(const ColoredCDBG& o) {

    CompactedDBG<DataAccessor<U>, DataStorage<U>>::operator=(o);

    invalid = o.invalid;

    return *this;
}

template<typename U>
ColoredCDBG<U>& ColoredCDBG<U>::operator=(ColoredCDBG&& o) {

    if (this != &o) {

        CompactedDBG<DataAccessor<U>, DataStorage<U>>::operator=(std::move(o));

        invalid = o.invalid;

        o.clear();
    }

    return *this;
}

template<typename U>
ColoredCDBG<U>& ColoredCDBG<U>::operator+=(const ColoredCDBG& o) {

    if (this != &o) merge(o, 1, false);

    return *this;
}

template<typename U>
bool ColoredCDBG<U>::operator==(const ColoredCDBG& o) const {

    if (!invalid && !this->isInvalid() && !o.isInvalid() && (this->getK() == o.getK()) && (this->size() == o.size())){

        const size_t k = this->getK();

        for (const auto& unitig : *this){

            const_UnitigColorMap<U> unitig_o(o.find(unitig.getUnitigHead(), true));

            if (unitig_o.isEmpty) return false;
            else {

                unitig_o.dist = 0;
                unitig_o.len = unitig_o.size - k + 1;

                const std::string unitig_o_str = unitig_o.strand ? unitig_o.referenceUnitigToString() : reverse_complement(unitig_o.referenceUnitigToString());

                if (unitig_o_str != unitig.referenceUnitigToString()) return false;
                else {

                    const UnitigColors* uc = unitig.getData()->getUnitigColors(unitig);
                    const UnitigColors* uc_o = unitig_o.getData()->getUnitigColors(unitig_o);

                    if ((uc != nullptr) && (uc_o != nullptr)){

                        if (!uc->isEqual(unitig, *uc_o, unitig_o)) return false;
                    }
                    else if ((uc != nullptr) != (uc_o != nullptr)) return false;
                }
            }
        }

        return true;
    }

    return false;
}

template<typename U>
inline bool ColoredCDBG<U>::operator!=(const ColoredCDBG& o) const {

    return !operator==(o);
}

template<typename U>
bool ColoredCDBG<U>::merge(const ColoredCDBG& o, const size_t nb_threads, const bool verbose){

    bool ret = true;

    if (invalid){

         if (verbose) std::cerr << "ColoredCDBG::merge(): Current graph is invalid." << std::endl;
         ret = false;
    }

    if (o.invalid){

         if (verbose) std::cerr << "ColoredCDBG::merge(): Graph to merge is invalid." << std::endl;
         ret = false;
    }

    if (this->getK() != o.getK()){

         if (verbose) std::cerr << "ColoredCDBG::merge(): The graphs to merge do not have the same k-mer length." << std::endl;
         ret = false;
    }

    if (this == &o){

         if (verbose) std::cerr << "ColoredCDBG::merge(): Cannot merge graph with itself." << std::endl;
         ret = false;
    }

    if (ret){

        const size_t sz_before = this->size();

        for (auto& unitig : *this) unitig.setFullCoverage();

        ret = CompactedDBG<DataAccessor<U>, DataStorage<U>>::annotateSplitUnitigs(o, nb_threads, verbose);

        if (ret){

            const size_t sz_after = this->size();
            const std::pair<size_t, size_t> p1 = CompactedDBG<DataAccessor<U>, DataStorage<U>>::getSplitInfoAllUnitigs();

            resizeDataUC(sz_after + (p1.second - p1.first), nb_threads);

            const std::pair<size_t, size_t> p2 = CompactedDBG<DataAccessor<U>, DataStorage<U>>::splitAllUnitigs();
            const size_t joined = (p1.second != 0) ? CompactedDBG<DataAccessor<U>, DataStorage<U>>::joinUnitigs() : 0;

            if (verbose){

                std::cout << "CompactedDBG::merge(): Added " << (sz_after - sz_before) << " new unitigs." << std::endl;
                std::cout << "CompactedDBG::merge(): Split " << p1.first << " unitigs into " << p1.second << " new unitigs." << std::endl;
                std::cout << "CompactedDBG::merge(): Joined " << joined << " unitigs." << std::endl;
                std::cout << "CompactedDBG::merge(): " << this->size() << " unitigs after merging." << std::endl;
            }

            for (size_t i = 0; i < o.getNbColors(); ++i) this->getData()->color_names.push_back(o.getColorName(i));

            return CompactedDBG<DataAccessor<U>, DataStorage<U>>::mergeData(o, nb_threads, verbose);
        }
    }

    return false;
}

template<typename U>
bool ColoredCDBG<U>::merge(ColoredCDBG&& o, const size_t nb_threads, const bool verbose){

    bool ret = true;

    if (invalid){

         if (verbose) std::cerr << "ColoredCDBG::merge(): Current graph is invalid." << std::endl;
         ret = false;
    }

    if (o.invalid){

         if (verbose) std::cerr << "ColoredCDBG::merge(): Graph to merge is invalid." << std::endl;
         ret = false;
    }

    if (this->getK() != o.getK()){

         if (verbose) std::cerr << "ColoredCDBG::merge(): The graphs to merge do not have the same k-mer length." << std::endl;
         ret = false;
    }

    if (this == &o){

         if (verbose) std::cerr << "ColoredCDBG::merge(): Cannot merge graph with itself." << std::endl;
         ret = false;
    }

    if (ret){

        const size_t sz_before = this->size();

        for (auto& unitig : *this) unitig.setFullCoverage();

        ret = CompactedDBG<DataAccessor<U>, DataStorage<U>>::annotateSplitUnitigs(o, nb_threads, verbose);

        if (ret){

            const size_t sz_after = this->size();
            const std::pair<size_t, size_t> p1 = CompactedDBG<DataAccessor<U>, DataStorage<U>>::getSplitInfoAllUnitigs();

            resizeDataUC(sz_after + (p1.second - p1.first), nb_threads);

            const std::pair<size_t, size_t> p2 = CompactedDBG<DataAccessor<U>, DataStorage<U>>::splitAllUnitigs();
            const size_t joined = (p1.second != 0) ? CompactedDBG<DataAccessor<U>, DataStorage<U>>::joinUnitigs() : 0;

            if (verbose){

                std::cout << "CompactedDBG::merge(): Added " << (sz_after - sz_before) << " new unitigs." << std::endl;
                std::cout << "CompactedDBG::merge(): Split " << p1.first << " unitigs into " << p1.second << " new unitigs." << std::endl;
                std::cout << "CompactedDBG::merge(): Joined " << joined << " unitigs." << std::endl;
                std::cout << "CompactedDBG::merge(): " << this->size() << " unitigs after merging." << std::endl;
            }

            for (size_t i = 0; i < o.getNbColors(); ++i) this->getData()->color_names.push_back(o.getColorName(i));

            const bool ret = CompactedDBG<DataAccessor<U>, DataStorage<U>>::mergeData(std::move(o), nb_threads, verbose);

            o.clear();

            return ret;
        }
    }

    return false;
}

template<typename U>
bool ColoredCDBG<U>::merge(const std::vector<ColoredCDBG>& v, const size_t nb_threads, const bool verbose){

    bool ret = true;

    if (invalid){

         if (verbose) std::cerr << "ColoredCDBG::merge(): Current graph is invalid." << std::endl;
         ret = false;
    }

    for (const auto& ccdbg : v){

        if (ccdbg.invalid){

             if (verbose) std::cerr << "ColoredCDBG::merge(): One of the graph to merge is invalid." << std::endl;
             ret = false;
        }

        if (this->getK() != ccdbg.getK()){

             if (verbose) std::cerr << "ColoredCDBG::merge(): The graphs to merge do not have the same k-mer length." << std::endl;
             ret = false;
        }

        if (this == &ccdbg){

             if (verbose) std::cerr << "ColoredCDBG::merge(): Cannot merge graph with itself." << std::endl;
             ret = false;
        }
    }

    if (ret){

        const size_t sz_before = this->size();

        for (auto& unitig : *this) unitig.setFullCoverage();

        for (const auto& ccdbg : v){

            ret = CompactedDBG<DataAccessor<U>, DataStorage<U>>::annotateSplitUnitigs(ccdbg, nb_threads, verbose);

            if (!ret) break;
        }

        if (ret){

            const size_t sz_after = this->size();
            const std::pair<size_t, size_t> p1 = CompactedDBG<DataAccessor<U>, DataStorage<U>>::getSplitInfoAllUnitigs();

            resizeDataUC(sz_after + (p1.second - p1.first), nb_threads);

            const std::pair<size_t, size_t> p2 = CompactedDBG<DataAccessor<U>, DataStorage<U>>::splitAllUnitigs();
            const size_t joined = (p1.second != 0) ? CompactedDBG<DataAccessor<U>, DataStorage<U>>::joinUnitigs() : 0;

            if (verbose){

                std::cout << "CompactedDBG::merge(): Added " << (sz_after - sz_before) << " new unitigs." << std::endl;
                std::cout << "CompactedDBG::merge(): Split " << p1.first << " unitigs into " << p1.second << " new unitigs." << std::endl;
                std::cout << "CompactedDBG::merge(): Joined " << joined << " unitigs." << std::endl;
                std::cout << "CompactedDBG::merge(): " << this->size() << " unitigs after merging." << std::endl;
            }

            for (const auto& ccdbg : v){

                for (size_t i = 0; i < ccdbg.getNbColors(); ++i) this->getData()->color_names.push_back(ccdbg.getColorName(i));

                if (!CompactedDBG<DataAccessor<U>, DataStorage<U>>::mergeData(ccdbg, nb_threads, verbose)) return false;
            }

            return true;
        }
    }

    return false;
}

template<typename U>
bool ColoredCDBG<U>::merge(std::vector<ColoredCDBG>&& v, const size_t nb_threads, const bool verbose){

    bool ret = true;

    if (invalid){

         if (verbose) std::cerr << "ColoredCDBG::merge(): Current graph is invalid." << std::endl;
         ret = false;
    }

    for (const auto& ccdbg : v){

        if (ccdbg.invalid){

             if (verbose) std::cerr << "ColoredCDBG::merge(): One of the graph to merge is invalid." << std::endl;
             ret = false;
        }

        if (this->getK() != ccdbg.getK()){

             if (verbose) std::cerr << "ColoredCDBG::merge(): The graphs to merge do not have the same k-mer length." << std::endl;
             ret = false;
        }

        if (this == &ccdbg){

             if (verbose) std::cerr << "ColoredCDBG::merge(): Cannot merge graph with itself." << std::endl;
             ret = false;
        }
    }

    if (ret){

        const size_t sz_before = this->size();

        for (auto& unitig : *this) unitig.setFullCoverage();

        for (const auto& ccdbg : v){

            ret = CompactedDBG<DataAccessor<U>, DataStorage<U>>::annotateSplitUnitigs(ccdbg, nb_threads, verbose);

            if (!ret) break;
        }

        if (ret){

            const size_t sz_after = this->size();
            const std::pair<size_t, size_t> p1 = CompactedDBG<DataAccessor<U>, DataStorage<U>>::getSplitInfoAllUnitigs();

            resizeDataUC(sz_after + (p1.second - p1.first), nb_threads);

            const std::pair<size_t, size_t> p2 = CompactedDBG<DataAccessor<U>, DataStorage<U>>::splitAllUnitigs();
            const size_t joined = (p1.second != 0) ? CompactedDBG<DataAccessor<U>, DataStorage<U>>::joinUnitigs() : 0;

            if (verbose){

                std::cout << "CompactedDBG::merge(): Added " << (sz_after - sz_before) << " new unitigs." << std::endl;
                std::cout << "CompactedDBG::merge(): Split " << p1.first << " unitigs into " << p1.second << " new unitigs." << std::endl;
                std::cout << "CompactedDBG::merge(): Joined " << joined << " unitigs." << std::endl;
                std::cout << "CompactedDBG::merge(): " << this->size() << " unitigs after merging." << std::endl;
            }

            for (auto& ccdbg : v){

                for (size_t i = 0; i < ccdbg.getNbColors(); ++i) this->getData()->color_names.push_back(ccdbg.getColorName(i));

                if (!CompactedDBG<DataAccessor<U>, DataStorage<U>>::mergeData(std::move(ccdbg), nb_threads, verbose)) return false;

                ccdbg.clear();
            }

            return true;
        }
    }

    return false;
}

template<typename U>
bool ColoredCDBG<U>::buildGraph(const CCDBG_Build_opt& opt){

    if (!invalid){

        CDBG_Build_opt opt_ = opt;

        invalid = !this->build(opt_);
    }
    else std::cerr << "ColoredCDBG::buildGraph(): Graph is invalid and cannot be built." << std::endl;

    return !invalid;
}

template<typename U>
bool ColoredCDBG<U>::buildColors(const CCDBG_Build_opt& opt){

    if (!invalid){

        initUnitigColors(opt);
        buildUnitigColors(opt.nb_threads);
    }
    else std::cerr << "ColoredCDBG::buildColors(): Graph is invalid (maybe not built yet?) and colors cannot be mapped." << std::endl;

    return !invalid;
}

template<typename U>
bool ColoredCDBG<U>::write(const std::string& prefix_output_fn, const size_t nb_threads, const bool write_index_file, const bool compress_output, const bool verbose) const {

    if (!CompactedDBG<DataAccessor<U>, DataStorage<U>>::write(prefix_output_fn, nb_threads, true, false, false, write_index_file, compress_output, verbose)) return false; // Write graph

    return this->getData()->write(prefix_output_fn, verbose); // Write colors
}

template<typename U>
bool ColoredCDBG<U>::loadColors(const std::string& input_graph_fn, const std::string& input_colors_fn, const size_t nb_threads, const bool verbose) {

    if (!this->getData()->read(input_colors_fn, nb_threads, verbose)) return false; // Read colors

    if (verbose) std::cout << "ColoredCDBG::loadColors(): Joining unitigs to their color sets." << std::endl;

    GFA_Parser graph(input_graph_fn);

    graph.open_read();

    auto reading_function = [&graph](std::vector<std::pair<Kmer, uint8_t>>& unitig_tags, const size_t chunk_size) {

        size_t i = 0;
        size_t graph_file_id = 0;

        bool new_file_opened = false;

        GFA_Parser::GFA_line r = graph.read(graph_file_id, new_file_opened, true);

        while ((i < chunk_size) && ((r.first != nullptr) || (r.second != nullptr))){

            if (r.first != nullptr){ // It is a sequence

                if (r.first->tags.empty()){

                    std::cerr << "ColoredCDBG::loadColors(): One sequence line in GFA file has no DataAccessor tag. Operation aborted." << std::endl;
                    return false;
                }

                size_t i = 0;

                for (; i < r.first->tags.size(); ++i){

                    if (r.first->tags[i].substr(0, 5) == "DA:Z:") break;
                }

                if (i == r.first->tags.size()){

                    std::cerr << "ColoredCDBG::loadColors(): One sequence line in GFA file has no DataAccessor tag. Operation aborted." << std::endl;
                    return false;
                }

                unitig_tags.push_back({Kmer(r.first->seq.c_str()), atoi(r.first->tags[i].c_str() + 5)});

                ++i;
            }

            r = graph.read(graph_file_id, new_file_opened, true);
        }

        return ((r.first != nullptr) || (r.second != nullptr));
    };

    auto join_function = [this](const std::vector<std::pair<Kmer, uint8_t>>& unitig_tags) {

        for (const auto& p : unitig_tags){

            UnitigColorMap<U> ucm(this->find(p.first, true));

            if (ucm.isEmpty){

                std::cerr << "ColoredCDBG::loadColors(): Internal error, operation aborted." << std::endl;
                std::cerr << "ColoredCDBG::loadColors(): A unitig from GFA file is not found in the in-memory graph." << std::endl;
                std::cerr << "ColoredCDBG::loadColors(): Graph from GFA file possibly incorrectly compacted." << std::endl;

                return false;
            }

            DataAccessor<U>* da = ucm.getData();

            *da = DataAccessor<U>(p.second);

            if (!ucm.strand){ // Unitig has been inserted in reverse-complement, need to reverse order of color sets

                UnitigColors* uc = da->getUnitigColors(ucm);
                UnitigColors r_uc = uc->reverse(ucm);

                *uc = std::move(r_uc);
            }
        }

        return true;
    };

    {
        const size_t chunk = 10000;

        std::vector<std::thread> workers; // need to keep track of threads so we can join them

        std::mutex mutex_file;

        bool file_valid_for_read = true;

        for (size_t t = 0; t < nb_threads; ++t){

            workers.emplace_back(

                [&]{

                    std::vector<std::pair<Kmer, uint8_t>> v;

                    while (true) {

                        {
                            std::unique_lock<std::mutex> lock(mutex_file);

                            if (!file_valid_for_read) return;

                            file_valid_for_read = reading_function(v, chunk);

                        }

                        join_function(v);
                        v.clear();
                    }
                }
            );
        }

        for (auto& t : workers) t.join();
    }

    return true;
}

template<typename U>
bool ColoredCDBG<U>::read(const std::string& input_graph_fn, const std::string& input_colors_fn, const size_t nb_threads, const bool verbose) {

    bool valid_input_files = true;

    if (input_graph_fn.length() != 0){

        if (check_file_exists(input_graph_fn)){

            FILE* fp = fopen(input_graph_fn.c_str(), "r");

            if (fp == NULL) {

                std::cerr << "ColoredCDBG::read(): Could not open input graph file " << input_graph_fn << std::endl;
                valid_input_files = false;
            }
            else fclose(fp);
        }
        else {

            std::cerr << "ColoredCDBG::read(): Input graph file " << input_graph_fn << " does not exist." << std::endl;
            valid_input_files = false;
        }
    }
    else {

        std::cerr << "ColoredCDBG::read(): No input graph file provided." << std::endl;
        valid_input_files = false;
    }

    if (input_colors_fn.length() != 0){

        if (check_file_exists(input_colors_fn)){

            FILE* fp = fopen(input_colors_fn.c_str(), "rb");

            if (fp == NULL) {

                std::cerr << "ColoredCDBG::read(): Could not open input colors file " << input_colors_fn << std::endl;
                valid_input_files = false;
            }
            else fclose(fp);
        }
        else {

            std::cerr << "ColoredCDBG::read(): Input colors file " << input_colors_fn << " does not exist." << std::endl;
            valid_input_files = false;
        }
    }
    else {

        std::cerr << "ColoredCDBG::read(): No input colors file provided." << std::endl;
        valid_input_files = false;
    }

    if (valid_input_files){

        if (verbose) std::cout << "ColoredCDBG::read(): Reading graph." << std::endl;
        
        invalid = !CompactedDBG<DataAccessor<U>, DataStorage<U>>::read(input_graph_fn, nb_threads, verbose);

        if (invalid) return false; // Read graph
        if (verbose) std::cout << "ColoredCDBG::read(): Reading colors." << std::endl;
        
        invalid = !loadColors(input_graph_fn, input_colors_fn, nb_threads, verbose);

        if (invalid) return false;
    }

    return valid_input_files;
}

template<typename U>
bool ColoredCDBG<U>::read(const std::string& input_graph_fn, const std::string& input_index_fn, const std::string& input_colors_fn, const size_t nb_threads, const bool verbose) {

    bool valid_input_files = true;

    if (input_graph_fn.length() != 0){

        if (check_file_exists(input_graph_fn)){

            FILE* fp = fopen(input_graph_fn.c_str(), "r");

            if (fp == NULL) {

                std::cerr << "ColoredCDBG::read(): Could not open input graph file " << input_graph_fn << std::endl;
                valid_input_files = false;
            }
            else fclose(fp);
        }
        else {

            std::cerr << "ColoredCDBG::read(): Input graph file " << input_graph_fn << " does not exist." << std::endl;
            valid_input_files = false;
        }
    }
    else {

        std::cerr << "ColoredCDBG::read(): No input graph file provided." << std::endl;
        valid_input_files = false;
    }

    if (input_colors_fn.length() != 0){

        if (check_file_exists(input_colors_fn)){

            FILE* fp = fopen(input_colors_fn.c_str(), "rb");

            if (fp == NULL) {

                std::cerr << "ColoredCDBG::read(): Could not open input colors file " << input_colors_fn << std::endl;
                valid_input_files = false;
            }
            else fclose(fp);
        }
        else {

            std::cerr << "ColoredCDBG::read(): Input colors file " << input_colors_fn << " does not exist." << std::endl;
            valid_input_files = false;
        }
    }
    else {

        std::cerr << "ColoredCDBG::read(): No input colors file provided." << std::endl;
        valid_input_files = false;
    }

    if (input_index_fn.length() != 0){

        if (check_file_exists(input_index_fn)){

            FILE* fp = fopen(input_index_fn.c_str(), "rb");

            if (fp == NULL) {

                std::cerr << "ColoredCDBG::read(): Could not open input index file " << input_index_fn << std::endl;
                valid_input_files = false;
            }
            else fclose(fp);
        }
        else {

            std::cerr << "ColoredCDBG::read(): Input index file " << input_index_fn << " does not exist." << std::endl;
            valid_input_files = false;
        }
    }
    else {

        std::cerr << "ColoredCDBG::read(): No input index file provided." << std::endl;
        valid_input_files = false;
    }

    if (valid_input_files){

        if (verbose) std::cout << "ColoredCDBG::read(): Reading graph." << std::endl;
        
        invalid = !CompactedDBG<DataAccessor<U>, DataStorage<U>>::read(input_graph_fn, input_index_fn, nb_threads, verbose);

        if (invalid) return false; // Read graph
        if (verbose) std::cout << "ColoredCDBG::read(): Reading colors." << std::endl;
        
        invalid = !loadColors(input_graph_fn, input_colors_fn, nb_threads, verbose);

        if (invalid) return false;
    }

    return valid_input_files;
}

template<typename U>
bool ColoredCDBG<U>::readGraph(const std::string& input_graph_fn, const size_t nb_threads, const bool verbose) {

    bool valid_input_files = true;

    if (input_graph_fn.length() != 0){

        if (check_file_exists(input_graph_fn)){

            FILE* fp = fopen(input_graph_fn.c_str(), "r");

            if (fp == NULL) {

                std::cerr << "ColoredCDBG::readGraph(): Could not open input graph file " << input_graph_fn << std::endl;
                valid_input_files = false;
            }
            else fclose(fp);
        }
        else {

            std::cerr << "ColoredCDBG::readGraph(): Input graph file " << input_graph_fn << " does not exist." << std::endl;
            valid_input_files = false;
        }
    }
    else {

        std::cerr << "ColoredCDBG::readGraph(): No input graph file provided." << std::endl;
        valid_input_files = false;
    }

    if (valid_input_files){

        if (verbose) std::cout << "ColoredCDBG::readGraph(): Reading graph." << std::endl;
        
        invalid = !CompactedDBG<DataAccessor<U>, DataStorage<U>>::read(input_graph_fn, nb_threads, verbose);

        if (invalid) return false; // Read graph
    }

    return valid_input_files;
}

template<typename U>
bool ColoredCDBG<U>::readGraph(const std::string& input_graph_fn, const std::string& input_index_fn, const size_t nb_threads, const bool verbose) {

    bool valid_input_files = true;

    if (input_graph_fn.length() != 0){

        if (check_file_exists(input_graph_fn)){

            FILE* fp = fopen(input_graph_fn.c_str(), "r");

            if (fp == NULL) {

                std::cerr << "ColoredCDBG::readGraph(): Could not open input graph file " << input_graph_fn << std::endl;
                valid_input_files = false;
            }
            else fclose(fp);
        }
        else {

            std::cerr << "ColoredCDBG::readGraph(): Input graph file " << input_graph_fn << " does not exist." << std::endl;
            valid_input_files = false;
        }
    }
    else {

        std::cerr << "ColoredCDBG::readGraph(): No input graph file provided." << std::endl;
        valid_input_files = false;
    }

    if (input_index_fn.length() != 0){

        if (check_file_exists(input_index_fn)){

            FILE* fp = fopen(input_index_fn.c_str(), "rb");

            if (fp == NULL) {

                std::cerr << "ColoredCDBG::readGraph(): Could not open input index file " << input_index_fn << std::endl;
                valid_input_files = false;
            }
            else fclose(fp);
        }
        else {

            std::cerr << "ColoredCDBG::readGraph(): Input index file " << input_index_fn << " does not exist." << std::endl;
            valid_input_files = false;
        }
    }
    else {

        std::cerr << "ColoredCDBG::readGraph(): No input index file provided." << std::endl;
        valid_input_files = false;
    }

    if (valid_input_files){

        if (verbose) std::cout << "ColoredCDBG::readGraph(): Reading graph." << std::endl;
        
        invalid = !CompactedDBG<DataAccessor<U>, DataStorage<U>>::read(input_graph_fn, input_index_fn, nb_threads, verbose);

        if (invalid) return false; // Read graph
    }

    return valid_input_files;
}

template<typename U>
void ColoredCDBG<U>::initUnitigColors(const CCDBG_Build_opt& opt, const size_t max_nb_hash){

    std::vector<std::string> v_files(opt.filename_seq_in);

    v_files.insert(v_files.end(), opt.filename_ref_in.begin(), opt.filename_ref_in.end());

    DataStorage<U>* ds = this->getData();
    DataStorage<U> new_ds(max_nb_hash, this->size(), v_files);

    *ds = std::move(new_ds);

    v_files.clear();

    const size_t chunk = 1000;

    std::vector<std::thread> workers; // need to keep track of threads so we can join them

    typename ColoredCDBG<U>::iterator g_a = this->begin();
    typename ColoredCDBG<U>::iterator g_b = this->end();

    std::mutex mutex_it;

    for (size_t t = 0; t < opt.nb_threads; ++t){

        workers.emplace_back(

            [&]{

                typename ColoredCDBG<U>::iterator l_a, l_b;

                while (true) {

                    {
                        std::unique_lock<std::mutex> lock(mutex_it);

                        if (g_a == g_b) return;

                        l_a = g_a;
                        l_b = g_a;

                        for (size_t cpt = 0; (cpt < chunk) && (l_b != g_b); ++cpt, ++l_b){}

                        g_a = l_b;
                    }

                    for (auto& it_unitig = l_a; it_unitig != l_b; ++it_unitig) {

                        *(it_unitig->getData()) = ds->insert(*it_unitig).first;
                    }
                }
            }
        );
    }

    for (auto& t : workers) t.join();

    //cout << "Number of unitigs not hashed is " << ds->overflow.size() << " on " << ds->nb_cs << " unitigs." << std::endl;
}

template<typename U>
void ColoredCDBG<U>::resizeDataUC(const size_t sz, const size_t nb_threads, const size_t max_nb_hash){

    DataStorage<U>* ds = this->getData();

    DataStorage<U> new_ds(max_nb_hash, sz, ds->color_names);

    const size_t chunk = 100;

    std::vector<std::thread> workers; // need to keep track of threads so we can join them

    typename ColoredCDBG<U>::iterator g_a = this->begin();
    typename ColoredCDBG<U>::iterator g_b = this->end();

    std::mutex mutex_it;

    for (size_t t = 0; t < nb_threads; ++t){

        workers.emplace_back(

            [&]{

                typename ColoredCDBG<U>::iterator l_a, l_b;

                while (true) {

                    {
                        std::unique_lock<std::mutex> lock(mutex_it);

                        if (g_a == g_b) return;

                        l_a = g_a;
                        l_b = g_a;

                        for (size_t cpt = 0; (cpt < chunk) && (l_b != g_b); ++cpt, ++l_b){}

                        g_a = l_b;
                    }

                    for (auto& it_unitig = l_a; it_unitig != l_b; ++it_unitig) {

                        UnitigColors* uc = ds->getUnitigColors(*it_unitig);
                        U* data = ds->getData(*it_unitig);

                        if ((uc != nullptr) || (data != nullptr)){

                            const std::pair<DataAccessor<U>, std::pair<UnitigColors*, U*>> p  = new_ds.insert(*it_unitig);

                            *(it_unitig->getData()) = p.first;

                            if (uc != nullptr) *(p.second.first) = std::move(*uc);
                            if (data != nullptr) *(p.second.second) = std::move(*data);
                        }
                    }
                }
            }
        );
    }

    for (auto& t : workers) t.join();

    *ds = std::move(new_ds);

    //cout << "Number of unitigs not hashed is " << ds->overflow.size() << " on " << ds->nb_cs << " unitigs." << std::endl;
}

template<>
inline void ColoredCDBG<void>::resizeDataUC(const size_t sz, const size_t nb_threads, const size_t max_nb_hash){

    DataStorage<void>* ds = this->getData();

    DataStorage<void> new_ds(max_nb_hash, sz, ds->color_names);

    const size_t chunk = 100;

    std::vector<std::thread> workers; // need to keep track of threads so we can join them

    typename ColoredCDBG<void>::iterator g_a = this->begin();
    typename ColoredCDBG<void>::iterator g_b = this->end();

    std::mutex mutex_it;

    for (size_t t = 0; t < nb_threads; ++t){

        workers.emplace_back(

            [&]{

                typename ColoredCDBG<void>::iterator l_a, l_b;

                while (true) {

                    {
                        std::unique_lock<std::mutex> lock(mutex_it);

                        if (g_a == g_b) return;

                        l_a = g_a;
                        l_b = g_a;

                        for (size_t cpt = 0; (cpt < chunk) && (l_b != g_b); ++cpt, ++l_b){}

                        g_a = l_b;
                    }

                    for (auto& it_unitig = l_a; it_unitig != l_b; ++it_unitig) {

                        UnitigColors* uc = ds->getUnitigColors(*it_unitig);

                        if (uc != nullptr){

                            const std::pair<DataAccessor<void>, std::pair<UnitigColors*, void*>> p  = new_ds.insert(*it_unitig);

                            *(it_unitig->getData()) = p.first;
                            *(p.second.first) = std::move(*uc);
                        }
                    }
                }
            }
        );
    }

    for (auto& t : workers) t.join();

    *ds = std::move(new_ds);

    //cout << "Number of unitigs not hashed is " << ds->overflow.size() << " on " << ds->nb_cs << " unitigs." << std::endl;
}

template<typename U>
void ColoredCDBG<U>::buildUnitigColors(const size_t nb_threads){

    DataStorage<U>* ds = this->getData();

    const int k_ = this->getK();

    const size_t nb_locks = nb_threads * 1024;
    const size_t chunk_size = 64;
    const size_t max_len_seq = rndup(static_cast<size_t>(1024 + k_ - 1));
    const size_t thread_seq_buf_sz = BUFFER_SIZE;
    const size_t thread_col_buf_sz = (thread_seq_buf_sz / (k_ + 1)) + 1;

    size_t prev_file_id = 0;

    size_t pos_read = 0;
    size_t len_read = 0;

    bool next_file = true;

    std::string s;

    FileParser fp(ds->color_names);

    std::atomic_flag* cs_locks = new std::atomic_flag[nb_locks];

    for (size_t i = 0; i < nb_locks; ++i) cs_locks[i].clear();

    // Main worker thread
    auto worker_function = [&](char* seq_buf, const size_t seq_buf_sz, const size_t* col_buf) {

        char* str = seq_buf;
        const char* str_end = &seq_buf[seq_buf_sz];

        size_t c_id = 0;

        while (str < str_end) { // for each input

            const int len = strlen(str);

            for (char* s = str; s != &str[len]; ++s) *s &= 0xDF;

            for (size_t i = 0; i < len - k_ + 1; i += max_len_seq - k_ + 1){

                const int curr_len = std::min(len - i, max_len_seq);
                const char saved_char = str[i + curr_len];
                const char* str_tmp = &str[i];

                str[i + curr_len] = '\0';

                for (KmerIterator it_km(str_tmp), it_km_end; it_km != it_km_end; ++it_km) {

                    UnitigColorMap<U> um = this->find(it_km->first);

                    if (!um.isEmpty) {

                        if (um.strand || (um.dist != 0)){

                            um.len = 1 + um.lcp(str_tmp, it_km->second + k_, um.strand ? um.dist + k_ : um.dist - 1, !um.strand);

                            //if ((um.size != k_) && !um.strand) um.dist -= um.len - 1;
                            um.dist -= (um.len - 1) & (static_cast<size_t>((um.size == k_) || um.strand) - 1);

                            it_km += um.len - 1;
                        }

                        const uint64_t id_lock = ds->getHash(um) % nb_locks;
                        UnitigColors* uc = ds->getUnitigColors(um);

                        while (cs_locks[id_lock].test_and_set(std::memory_order_acquire)); // Set the corresponding lock

                        uc->add(um, col_buf[c_id]);

                        cs_locks[id_lock].clear(std::memory_order_release);
                    }
                }

                str[i + curr_len] = saved_char;
            }

            str += len + 1;
            ++c_id;
        }
    };

    auto reading_function = [&](char* seq_buf, size_t& seq_buf_sz, size_t* col_buf) {

        size_t file_id = prev_file_id;
        size_t i = 0;

        const size_t sz_buf = thread_seq_buf_sz - k_;

        const char* s_str = s.c_str();

        seq_buf_sz = 0;

        while (seq_buf_sz < sz_buf) {

            const bool new_reading = (pos_read >= len_read);

            if (!new_reading || fp.read(s, file_id)) {

                //pos_read = (new_reading ? 0 : pos_read);
                pos_read &= static_cast<size_t>(new_reading) - 1;

                len_read = s.length();
                s_str = s.c_str();

                if (len_read >= k_){

                    if ((thread_seq_buf_sz - seq_buf_sz - 1) < (len_read - pos_read)){

                        strncpy(&seq_buf[seq_buf_sz], &s_str[pos_read], thread_seq_buf_sz - seq_buf_sz - 1);

                        seq_buf[thread_seq_buf_sz - 1] = '\0';
                        col_buf[i++] = file_id;

                        pos_read += sz_buf - seq_buf_sz;
                        seq_buf_sz = thread_seq_buf_sz;

                        break;
                    }
                    else {

                        strcpy(&seq_buf[seq_buf_sz], &s_str[pos_read]);

                        col_buf[i++] = file_id;

                        seq_buf_sz += (len_read - pos_read) + 1;
                        pos_read = len_read;
                    }
                }
                else pos_read = len_read;
            }
            else {

                next_file = false;

                return true;
            }
        }

        const bool ret = (file_id != prev_file_id);

        next_file = true;
        prev_file_id = file_id;

        return ret;
    };

    {
        bool stop = false;

        std::vector<std::thread> workers; // need to keep track of threads so we can join them

        std::mutex mutex_file;

        size_t prev_uc_sz = getCurrentRSS();

        while (next_file){

            stop = false;

            for (size_t t = 0; t < nb_threads; ++t){

                workers.emplace_back(

                    [&, t]{

                        char* buffer_seq = new char[thread_seq_buf_sz];
                        size_t* buffer_col = new size_t[thread_col_buf_sz];

                        size_t buffer_seq_sz = 0;

                        while (true) {

                            {
                                std::unique_lock<std::mutex> lock(mutex_file);

                                if (stop) {

                                    delete[] buffer_seq;
                                    delete[] buffer_col;

                                    return;
                                }

                                stop = reading_function(buffer_seq, buffer_seq_sz, buffer_col);
                            }

                            worker_function(buffer_seq, buffer_seq_sz, buffer_col);
                        }

                        delete[] buffer_seq;
                        delete[] buffer_col;
                    }
                );
            }

            for (auto& t : workers) t.join();

            workers.clear();

            const size_t curr_uc_sz = getCurrentRSS();

            if ((curr_uc_sz - prev_uc_sz) >= 1073741824ULL){

                const size_t chunk = 1000;

                typename ColoredCDBG<U>::iterator g_a = this->begin();
                typename ColoredCDBG<U>::iterator g_b = this->end();

                std::mutex mutex_it;

                for (size_t t = 0; t < nb_threads; ++t){

                    workers.emplace_back(

                        [&]{

                            typename ColoredCDBG<U>::iterator l_a, l_b;

                            while (true) {

                                {
                                    std::unique_lock<std::mutex> lock(mutex_it);

                                    if (g_a == g_b) return;

                                    l_a = g_a;
                                    l_b = g_a;

                                    for (size_t cpt = 0; (cpt < chunk) && (l_b != g_b); ++cpt, ++l_b){}

                                    g_a = l_b;
                                }

                                while (l_a != l_b){

                                    l_a->getData()->getUnitigColors(*l_a)->optimizeFullColors(*l_a);
                                    ++l_a;
                                }
                            }
                        }
                    );
                }

                for (auto& t : workers) t.join();

                workers.clear();

                prev_uc_sz = getCurrentRSS();
            }
        }
    }

    fp.close();

    //checkColors(ds->color_names);

    /*typedef std::unordered_map<uint64_t, pair<int64_t, size_t>> uc_unordered_map;

    uc_unordered_map u_map;

    mutex mutex_u_map;

    std::vector<std::thread> workers;

    auto add_hash_function = [&](typename ColoredCDBG<U>::iterator it_a, typename ColoredCDBG<U>::iterator it_b) {

        while (it_a != it_b) {

            const const_UnitigColorMap<U> unitig(*it_a);
            const UnitigColors* uc = unitig.getData()->getUnitigColors(unitig);
            const UnitigColors uc_full = uc->makeFullColors(unitig);
            const UnitigColors* uc_full_array = uc_full.getFullColorsPtr();

            if (uc_full_array[0].size() != 0){

                const std::pair<int64_t, size_t> pv(0 - static_cast<int64_t>(uc_full_array[0].getSizeInBytes()) - static_cast<int64_t>(sizeof(size_t)), 0);

                const int64_t to_add = static_cast<int64_t>(uc->getSizeInBytes());
                const int64_t to_rm = (static_cast<int64_t>(uc_full_array[1].getSizeInBytes() + 2 * sizeof(UnitigColors)));

                {
                    unique_lock<mutex> lock(mutex_u_map);

                    pair<uc_unordered_map::iterator, bool> p = u_map.insert(make_pair(uc_full_array[0].hash(), pv));

                    p.first->second.first += to_add - to_rm;
                }
            }

            ++it_a;
        }
    };

    auto add_shared_function = [&](typename ColoredCDBG<U>::iterator it_a, typename ColoredCDBG<U>::iterator it_b) {

        while (it_a != it_b) {

            const UnitigColorMap<U> unitig(*it_a);

            UnitigColors* uc = unitig.getData()->getUnitigColors(unitig);
            UnitigColors uc_full = uc->makeFullColors(unitig);
            UnitigColors* uc_full_array = uc_full.getFullColorsPtr();

            if (uc_full_array[0].size() != 0){

                uc_unordered_map::const_iterator it = u_map.find(uc_full_array[0].hash());

                if (it->second.first > 0){ // If there is some sharing

                    const size_t id_shared = it->second.second;
                    const uint64_t id_lock = id_shared % nb_locks;

                    bool move_full = false;

                    while (cs_locks[id_lock].test_and_set(std::memory_order_acquire)); // Set the corresponding lock

                    if (ds->shared_color_sets[id_shared].second == 0){

                        ds->shared_color_sets[id_shared].first = std::move(uc_full_array[0]);
                        ds->shared_color_sets[id_shared].second = 0;

                        uc_full_array[0] = ds->shared_color_sets[id_shared];
                        move_full = true;
                    }
                    else if (uc_full_array[0] == ds->shared_color_sets[id_shared]) {

                        uc_full_array[0] = ds->shared_color_sets[id_shared];
                        move_full = true;
                    }

                    cs_locks[id_lock].clear(std::memory_order_release);

                    if (move_full) *uc = std::move(uc_full);
                }
            }

            ++it_a;
        }
    };

    {
        const size_t chunk = 1000;

        typename ColoredCDBG<U>::iterator g_a = this->begin();
        typename ColoredCDBG<U>::iterator g_b = this->end();

        std::mutex mutex_it;

        for (size_t t = 0; t < nb_threads; ++t){

            workers.emplace_back(

                [&, t]{

                    typename ColoredCDBG<U>::iterator l_a, l_b;

                    while (true) {

                        {
                            std::unique_lock<std::mutex> lock(mutex_it);

                            if (g_a == g_b) return;

                            l_a = g_a;
                            l_b = g_a;

                            for (size_t cpt = 0; (cpt < chunk) && (l_b != g_b); ++cpt, ++l_b){}

                            g_a = l_b;
                        }

                        add_hash_function(l_a, l_b);
                    }
                }
            );
        }

        for (auto& t : workers) t.join();

        workers.clear();
    }

    for (auto& p : u_map){

        if (p.second.first > 0){

            p.second.second = ds->sz_shared_cs;
            ++(ds->sz_shared_cs);
        }
    }

    ds->shared_color_sets = new UnitigColors::SharedUnitigColors[ds->sz_shared_cs];

    {
        const size_t chunk = 1000;

        typename ColoredCDBG<U>::iterator g_a = this->begin();
        typename ColoredCDBG<U>::iterator g_b = this->end();

        std::mutex mutex_it;

        for (size_t t = 0; t < nb_threads; ++t){

            workers.emplace_back(

                [&, t]{

                    typename ColoredCDBG<U>::iterator l_a, l_b;

                    while (true) {

                        {
                            std::unique_lock<std::mutex> lock(mutex_it);

                            if (g_a == g_b) return;

                            l_a = g_a;
                            l_b = g_a;

                            for (size_t cpt = 0; (cpt < chunk) && (l_b != g_b); ++cpt, ++l_b){}

                            g_a = l_b;
                        }

                        add_shared_function(l_a, l_b);
                    }
                }
            );
        }

        for (auto& t : workers) t.join();

        workers.clear();
    }*/

    delete[] cs_locks;
}

template<typename U>
std::string ColoredCDBG<U>::getColorName(const size_t color_id) const {

    if (invalid){

        std::cerr << "ColoredCDBG::getColorName(): Graph is invalid or colors are not yet mapped to unitigs." << std::endl;
        return std::string();
    }

    const DataStorage<U>* ds = this->getData();

    if (color_id >= ds->color_names.size()){

        std::cerr << "ColoredCDBG::getColorName(): Color ID " << color_id << " is invalid, graph only has " <<
        ds->color_names.size() << " colors." << std::endl;

        return std::string();
    }

    return ds->color_names[color_id];
}

template<typename U>
std::vector<std::string> ColoredCDBG<U>::getColorNames() const {

    if (invalid){

        std::cerr << "ColoredCDBG::getColorNames(): Graph is invalid or colors are not yet mapped to unitigs." << std::endl;
        return std::vector<std::string>();
    }

    return this->getData()->color_names;
}

template<typename U>
bool ColoredCDBG<U>::searchMinRatioKmer(const std::vector<std::string>& query_filenames, const std::string& out_filename_prefix,
                                        const double min_ratio_kmers,
                                        const bool inexact_search, const bool files_as_queries,
                                        const size_t nb_threads, const bool verbose) const {

    const std::string out_tmp = out_filename_prefix + ".tsv";

    {
        FILE* fp_tmp = fopen(out_tmp.c_str(), "w");

        if (fp_tmp == NULL) {

            std::cerr << "ColoredCDBG::searchMinRatioKmer(): Could not open file " << out_tmp << " for writing." << std::endl;
            return false;
        }
        else {

            fclose(fp_tmp);

            if (std::remove(out_tmp.c_str()) != 0) std::cerr << "ColoredCDBG::searchMinRatioKmer(): Could not remove temporary file " << out_tmp << std::endl;
        }
    }

    std::ofstream outfile;
    std::ostream out(0);

    outfile.open(out_tmp.c_str());
    out.rdbuf(outfile.rdbuf());

    const bool ret = this->searchMinRatioKmer(query_filenames, out, min_ratio_kmers, inexact_search, files_as_queries, nb_threads, verbose);

    outfile.close();

    return ret;
}

template<typename U>
bool ColoredCDBG<U>::searchMinRatioKmer(const std::vector<std::string>& query_filenames, const std::string& out_filename_prefix,
                                        const double min_ratio_kmers, const size_t nb_min_colors,
                                        const bool inexact_search, const bool files_as_queries,
                                        const size_t nb_threads, const bool verbose) const {

    const std::string out_tmp = out_filename_prefix + ".tsv";

    {
        FILE* fp_tmp = fopen(out_tmp.c_str(), "w");

        if (fp_tmp == NULL) {

            std::cerr << "ColoredCDBG::searchMinRatioKmer(): Could not open file " << out_tmp << " for writing." << std::endl;
            return false;
        }
        else {

            fclose(fp_tmp);

            if (std::remove(out_tmp.c_str()) != 0) std::cerr << "ColoredCDBG::searchMinRatioKmer(): Could not remove temporary file " << out_tmp << std::endl;
        }
    }

    std::ofstream outfile;
    std::ostream out(0);

    outfile.open(out_tmp.c_str());
    out.rdbuf(outfile.rdbuf());

    const bool ret = this->searchMinRatioKmer(query_filenames, out, min_ratio_kmers, nb_min_colors, inexact_search, files_as_queries, nb_threads, verbose);

    outfile.close();

    return ret;
}

template<typename U>
bool ColoredCDBG<U>::search(const std::vector<std::string>& query_filenames, const std::string& out_filename_prefix,
                                    const bool found_km_ratio_out,
                                    const bool inexact_search, const bool files_as_queries,
                                    const size_t nb_threads, const bool verbose) const {

    const std::string out_tmp = out_filename_prefix + ".tsv";

    {
        FILE* fp_tmp = fopen(out_tmp.c_str(), "w");

        if (fp_tmp == NULL) {

            std::cerr << "ColoredCDBG::search(): Could not open file " << out_tmp << " for writing." << std::endl;
            return false;
        }
        else {

            fclose(fp_tmp);

            if (std::remove(out_tmp.c_str()) != 0) std::cerr << "ColoredCDBG::search(): Could not remove temporary file " << out_tmp << std::endl;
        }
    }

    std::ofstream outfile;
    std::ostream out(0);

    outfile.open(out_tmp.c_str());
    out.rdbuf(outfile.rdbuf());

    const bool ret = this->search(query_filenames, out, found_km_ratio_out, inexact_search, files_as_queries, nb_threads, verbose);

    outfile.close();

    return ret;
}

template<typename U>
bool ColoredCDBG<U>::searchMinRatioKmer(const std::vector<std::string>& query_filenames, std::ostream& out,
                                        const double min_ratio_kmers,
                                        const bool inexact_search, const bool files_as_queries,
                                        const size_t nb_threads, const bool verbose) const {

    return searchMinRatioKmer_( query_filenames, out, min_ratio_kmers, 0,
                                inexact_search, files_as_queries, nb_threads, verbose);
}

template<typename U>
bool ColoredCDBG<U>::searchMinRatioKmer(const std::vector<std::string>& query_filenames, std::ostream& out,
                                        const double min_ratio_kmers, const size_t min_nb_colors,
                                        const bool inexact_search, const bool files_as_queries,
                                        const size_t nb_threads, const bool verbose) const {

    if (min_nb_colors == 0) {

        std::cerr << "ColoredCDBG::searchMinRatioKmer(): Minimum number of required colors is 0." << std::endl;
        return false;
    }

    return searchMinRatioKmer_( query_filenames, out, min_ratio_kmers, min_nb_colors,
                                inexact_search, files_as_queries, nb_threads, verbose);
}

template<typename U>
bool ColoredCDBG<U>::searchMinRatioKmer_(   const std::vector<std::string>& query_filenames, std::ostream& out,
                                            const double min_ratio_kmers, const size_t min_nb_colors,
                                            const bool inexact_search, const bool files_as_queries,
                                            const size_t nb_threads, const bool verbose) const {

    if (invalid){

        std::cerr << "ColoredCDBG::searchMinRatioKmer(): Graph is invalid and cannot be searched" << std::endl;
        return false;
    }

    if (nb_threads > std::thread::hardware_concurrency()){

        std::cerr << "ColoredCDBG::searchMinRatioKmer(): Number of threads cannot be greater than or equal to " << std::thread::hardware_concurrency() << "." << std::endl;
        return false;
    }

    if (nb_threads <= 0){

        std::cerr << "ColoredCDBG::searchMinRatioKmer(): Number of threads cannot be less than or equal to 0." << std::endl;
        return false;
    }

    if (min_ratio_kmers <= 0.0){

        std::cerr << "ColoredCDBG::searchMinRatioKmer(): Ratio of k-mers is less than or equal to 0.0." << std::endl;
        return false;
    }

    if (min_ratio_kmers > 1.0){

        std::cerr << "ColoredCDBG::searchMinRatioKmer(): Ratio of k-mers is greater than 1.0." << std::endl;
        return false;
    }

    if (min_nb_colors > getNbColors()) {

        std::cerr << "ColoredCDBG::searchMinRatioKmer(): Minimum number of required colors is larger than total number of colors in graph." << std::endl;
        return false;
    }

    if (out.fail()) {

        std::cerr << "CompactedDBG::searchMinRatioKmer(): Output stream is in a failed state and cannot be written to." << std::endl;
        return false;
    }

    if (verbose) std::cout << "ColoredCDBG::searchMinRatioKmer(): Querying graph." << std::endl;

    const size_t k = this->getK();
    const size_t thread_seq_buf_sz = BUFFER_SIZE;

    const char query_pres[2] = {'\t', '1'};
    const char query_abs[2] = {'\t', '0'};

    const char eol = '\n';

    const bool nb_color_found_filtering = (min_nb_colors != 0);

    const size_t l_query_res = 2;
    const size_t nb_colors = getNbColors();
    const size_t sz_binary_color_query_out = nb_colors * l_query_res + 1;

    std::string s;

    bool write_success = true;
    bool query_success = true;

    size_t file_id = 0;
    size_t prev_file_id = 0xffffffffffffffffULL;

    auto processCounts = [&](const std::vector<std::pair<size_t, const_UnitigColorMap<U>>>& v_um, Roaring* color_occ_r, uint32_t* color_occ_u){

        struct hash_pair {

            size_t operator() (const std::pair<size_t, std::pair<Kmer, size_t>>& p) const {

                return wyhash(&p, sizeof(std::pair<size_t, std::pair<Kmer, size_t>>), 0, _wyp);
            }
        };

        std::unordered_set<std::pair<size_t, std::pair<Kmer, size_t>>, hash_pair> s_um;

        typename std::unordered_set<std::pair<size_t, std::pair<Kmer, size_t>>, hash_pair>::const_iterator it;

        for (const auto& p : v_um){

            s_um.insert({p.first, {p.second.strand ? p.second.getUnitigHead() : p.second.getUnitigTail().twin(), p.second.dist}});
        }

        for (const auto& p : v_um){

            const_UnitigColorMap<U> um = p.second;

            size_t pos_query = p.first;

            if (um.strand) {

                const Kmer head = um.getUnitigHead();

                size_t pos_unitig = um.dist;

                it = s_um.find({pos_query, {head, um.dist}});

                if (it != s_um.end()) {

                    s_um.erase(it);

                    while ((pos_unitig + k) < um.size){

                        ++pos_query;
                        ++pos_unitig;

                        it = s_um.find({pos_query, {head, pos_unitig}});

                        if (it == s_um.end()) break;
                        else {

                            ++(um.len);

                            s_um.erase(it);
                        }
                    }

                    const UnitigColors* uc = um.getData()->getUnitigColors(um);

                    UnitigColors::const_iterator it_uc = uc->begin(um);
                    UnitigColors::const_iterator it_uc_end = uc->end();

                    if (inexact_search){

                        for (; it_uc != it_uc_end; ++it_uc) color_occ_r[it_uc.getColorID()].add(it_uc.getKmerPosition() - um.dist + p.first);
                    }
                    else {

                        for (; it_uc != it_uc_end; ++it_uc) color_occ_u[it_uc.getColorID()] += 1;
                    }
                }
            }
            else {

                const Kmer head = um.getUnitigTail().twin();

                it = s_um.find({pos_query, {head, um.dist}});

                if (it != s_um.end()) {

                    s_um.erase(it);

                    while (um.dist > 0){

                        ++pos_query;

                        it = s_um.find({pos_query, {head, um.dist - 1}});

                        if (it == s_um.end()) break;
                        else {

                            --(um.dist);
                            ++(um.len);

                            s_um.erase(it);
                        }
                    }

                    const UnitigColors* uc = um.getData()->getUnitigColors(um);

                    UnitigColors::const_iterator it_uc = uc->begin(um);
                    UnitigColors::const_iterator it_uc_end = uc->end();

                    if (inexact_search){

                        const size_t max_pos_um = um.dist + um.len - 1;

                        for (; it_uc != it_uc_end; ++it_uc) color_occ_r[it_uc.getColorID()].add(max_pos_um - it_uc.getKmerPosition() + p.first);
                    }
                    else {

                        for (; it_uc != it_uc_end; ++it_uc) color_occ_u[it_uc.getColorID()] += 1;
                    }
                }
            }
        }

        if (inexact_search){

            for (size_t i = 0; i < nb_colors; ++i) color_occ_u[i] = color_occ_r[i].cardinality();
        }
    };

    auto searchQuery = [&](const std::string& query, Roaring* color_occ_r, uint32_t* color_occ_u, const size_t nb_km_min){

        const std::vector<std::pair<size_t, const_UnitigColorMap<U>>> v_um_e = this->searchSequence(query, true, false, false, false, false);

        processCounts(v_um_e, color_occ_r, color_occ_u); // Extract k-mer occurrences for each color

        if (inexact_search){

            if (!files_as_queries) {

                size_t nb_color_pres = 0;

                for (size_t j = 0; j < nb_colors; ++j) nb_color_pres += static_cast<size_t>(color_occ_u[j] >= nb_km_min);

                if (nb_color_found_filtering && (nb_color_pres >= min_nb_colors)) return;
                else if (!nb_color_found_filtering && (nb_color_pres == nb_colors)) return;
            }

            const std::vector<std::pair<size_t, const_UnitigColorMap<U>>> v_um_d = this->searchSequence(query, false, false, true, false, false);

            processCounts(v_um_d, color_occ_r, color_occ_u); // Extract k-mer occurrences for each color

            if (!files_as_queries) {
                
                size_t nb_color_pres = 0;

                for (size_t j = 0; j < nb_colors; ++j) nb_color_pres += static_cast<size_t>(color_occ_u[j] >= nb_km_min);

                if (nb_color_found_filtering && (nb_color_pres >= min_nb_colors)) return;
                else if (!nb_color_found_filtering && (nb_color_pres == nb_colors)) return;
            }

            const std::vector<std::pair<size_t, const_UnitigColorMap<U>>> v_um_m = this->searchSequence(query, false, false, false, true, false);

            processCounts(v_um_m, color_occ_r, color_occ_u); // Extract k-mer occurrences for each color

            if (!files_as_queries) {
                
                size_t nb_color_pres = 0;

                for (size_t j = 0; j < nb_colors; ++j) nb_color_pres += static_cast<size_t>(color_occ_u[j] >= nb_km_min);

                if (nb_color_found_filtering && (nb_color_pres >= min_nb_colors)) return;
                else if (!nb_color_found_filtering && (nb_color_pres == nb_colors)) return;
            }

            const std::vector<std::pair<size_t, const_UnitigColorMap<U>>> v_um_i = this->searchSequence(query, false, true, false, false, false);

            processCounts(v_um_i, color_occ_r, color_occ_u); // Extract k-mer occurrences for each color
        }
    };

    auto writeBinaryOutput = [&](   const char* query_name, const size_t len_query_name, const uint32_t* color_occ,
                                    char* buffer_res, size_t& pos_buffer_out, const size_t nb_km_min, std::mutex& mtx){

        bool is_found = false;

        if (nb_color_found_filtering) {

            size_t nb_colors_found = 0;

            for (size_t i = 0; i < nb_colors; ++i) nb_colors_found += static_cast<size_t>(color_occ[i] >= nb_km_min);

            is_found = (nb_colors_found >= min_nb_colors);

            if ((pos_buffer_out + len_query_name + l_query_res + 1) > thread_seq_buf_sz){

                std::unique_lock<std::mutex> lock(mtx);

                if (pos_buffer_out > 0) {

                    out.write(buffer_res, pos_buffer_out); // Write result buffer

                    pos_buffer_out = 0; // Reset position to 0;
                    write_success = (write_success && !out.fail());
                }
            }

            // Copy new result to buffer
            std::memcpy(buffer_res + pos_buffer_out, query_name, len_query_name * sizeof(char));
            std::memcpy(buffer_res + pos_buffer_out + len_query_name, is_found ? query_pres : query_abs, l_query_res * sizeof(char));

            pos_buffer_out += len_query_name + l_query_res;

            buffer_res[pos_buffer_out++] = eol;
        }
        else {

            if ((pos_buffer_out + len_query_name + sz_binary_color_query_out) > thread_seq_buf_sz){

                std::unique_lock<std::mutex> lock(mtx);

                if (pos_buffer_out > 0) {

                    out.write(buffer_res, pos_buffer_out); // Write result buffer

                    pos_buffer_out = 0; // Reset position to 0;
                    write_success = (write_success && !out.fail());
                }

                out.write(query_name, len_query_name * sizeof(char)); // Write title
                write_success = (write_success && !out.fail());

                for (size_t i = 0; i < nb_colors; ++i) {

                    if (color_occ[i] >= nb_km_min) {

                        out.write(query_pres, l_query_res * sizeof(char));
                        is_found = true;
                    }
                    else out.write(query_abs, l_query_res * sizeof(char));

                    write_success = (write_success && !out.fail());
                }

                out.write(&eol, sizeof(char));

                write_success = (write_success && !out.fail());
            }
            else {

                // Copy new result to buffer
                std::memcpy(buffer_res + pos_buffer_out, query_name, len_query_name * sizeof(char));

                pos_buffer_out += len_query_name;

                for (size_t i = 0; i < nb_colors; ++i, pos_buffer_out += l_query_res){

                    if (color_occ[i] >= nb_km_min) {

                        std::memcpy(buffer_res + pos_buffer_out, query_pres, l_query_res * sizeof(char));
                        is_found = true;
                    }
                    else std::memcpy(buffer_res + pos_buffer_out, query_abs, l_query_res * sizeof(char));
                }

                buffer_res[pos_buffer_out++] = eol;
            }
        }

        return is_found;
    };

    // Write header to TSV file
    if (write_success) {

        out << "query_name";

        if (nb_color_found_filtering) out << "\tpresence_query\n";
        else {

            const std::vector<std::string> color_names = getColorNames();

            for (const auto& name : color_names) out << '\t' << name;

            out << eol;
        }

        write_success = (write_success && !out.fail());
    }

    if (write_success) {

        FileParser fp(query_filenames);

        if (nb_threads == 1){

            const char* query_name = nullptr;

            char* buffer_res = new char[thread_seq_buf_sz];

            uint32_t* color_occ_u = new uint32_t[nb_colors]();
            uint32_t* color_occ_u_fid = files_as_queries ? new uint32_t[nb_colors]() : nullptr;

            Roaring* color_occ_r = inexact_search ? new Roaring[nb_colors] : nullptr;

            std::mutex mtx_file_out; // Dummy

            size_t pos_buffer_out = 0;
            size_t nb_queries_found = 0;
            size_t nb_queries_processed = 0;

            size_t nb_km_query = 0;

            while (write_success && fp.read(s, file_id)){

                // Write previous results if there are any
                {
                    const size_t nb_km_min = std::max(static_cast<size_t>(1), static_cast<size_t>(round(static_cast<double>(nb_km_query) * min_ratio_kmers)));

                    if (files_as_queries) {

                        for (size_t j = 0; j < nb_colors; ++j) color_occ_u_fid[j] += color_occ_u[j];

                        if (file_id != prev_file_id) {

                            if (prev_file_id != 0xffffffffffffffffULL) { // Push results to buffer, write buffer if overflow

                                nb_queries_found += static_cast<size_t>(writeBinaryOutput(query_name, strlen(query_name), color_occ_u_fid, buffer_res, pos_buffer_out, nb_km_min, mtx_file_out));
                                ++nb_queries_processed;
                            }

                            query_name = query_filenames[file_id].c_str(); // Query name is the filename
                            nb_km_query = 0;

                            std::memset(color_occ_u_fid, 0, nb_colors * sizeof(uint32_t));
                        }

                        nb_km_query += s.length() - k + 1;
                    }
                    else {

                        // Push results to buffer, write buffer if overflow
                        if (prev_file_id != 0xffffffffffffffffULL) {

                            nb_queries_found += static_cast<size_t>(writeBinaryOutput(query_name, strlen(query_name), color_occ_u, buffer_res, pos_buffer_out, nb_km_min, mtx_file_out));
                            ++nb_queries_processed;
                        }

                        query_name = fp.getNameString(); // Query name is the record name
                        nb_km_query = s.length() - k + 1;
                    }
                }

                // Clean-up
                {
                    std::memset(color_occ_u, 0, nb_colors * sizeof(uint32_t));

                    if (inexact_search){

                        for (size_t j = 0; j < nb_colors; ++j) color_occ_r[j] = Roaring(); // Reset k-mer occurences for each color
                    }
                }

                for (auto& c : s) c &= 0xDF;

                if (files_as_queries) searchQuery(s, color_occ_r, color_occ_u, s.length() - k + 1);
                else searchQuery(s, color_occ_r, color_occ_u, std::max(static_cast<size_t>(1), static_cast<size_t>(round(static_cast<double>(s.length() - k + 1) * min_ratio_kmers))));

                prev_file_id = file_id;
            }

            // Flush rest of buffer result to final output
            if (write_success && (prev_file_id != 0xffffffffffffffffULL)) {

                const size_t nb_km_min = std::max(static_cast<size_t>(1), static_cast<size_t>(round(static_cast<double>(nb_km_query) * min_ratio_kmers)));

                if (files_as_queries) {

                    for (size_t j = 0; j < nb_colors; ++j) color_occ_u_fid[j] += color_occ_u[j];

                    nb_queries_found += static_cast<size_t>(writeBinaryOutput(query_name, strlen(query_name), color_occ_u_fid, buffer_res, pos_buffer_out, nb_km_min, mtx_file_out));
                }
                else nb_queries_found += static_cast<size_t>(writeBinaryOutput(query_name, strlen(query_name), color_occ_u, buffer_res, pos_buffer_out, nb_km_min, mtx_file_out));

                ++nb_queries_processed;

                if (pos_buffer_out > 0) out.write(buffer_res, pos_buffer_out);
            }

            delete[] buffer_res;
            delete[] color_occ_u;

            if (color_occ_r != nullptr) delete[] color_occ_r;
            if (color_occ_u_fid != nullptr) delete[] color_occ_u_fid;

            if (write_success && verbose) {

                std::cout << "ColoredCDBG::searchMinRatioKmer(): Processed " << nb_queries_processed << " queries." << std::endl;
                std::cout << "ColoredCDBG::searchMinRatioKmer(): Found " << nb_queries_found << " queries in at least " << ((min_nb_colors == 0) ? 1 : min_nb_colors) << " color(s)." << std::endl;
            }
        }
        else {

            struct ResultFileQuery {

                uint32_t* color_occ;
                size_t nb_km_queries;
                size_t nb_queries;
                bool is_read;

                ResultFileQuery() : color_occ(nullptr), nb_km_queries(0), nb_queries(0), is_read(false) {}
            };

            bool stop = false;

            std::vector<std::thread> workers; // need to keep track of threads so we can join them

            std::atomic<size_t> nb_queries_found;
            std::atomic<size_t> nb_queries_processed;

            std::mutex mtx_files_in, mtx_file_out, mtx_file_id;

            std::unordered_map<size_t, ResultFileQuery> um_file_id;

            nb_queries_found = 0;
            nb_queries_processed = 0;

            for (size_t t = 0; t < nb_threads; ++t){

                workers.emplace_back(

                    [&]{

                        size_t pos_buffer_out = 0;

                        char* buffer_res = new char[thread_seq_buf_sz];
                        uint32_t* color_occ_u = new uint32_t[nb_colors]();
                        Roaring* color_occ_r = inexact_search ? new Roaring[nb_colors] : nullptr;

                        std::vector<std::string> buffer_seq;
                        std::vector<std::string> buffer_name;
                        std::vector<size_t> buffer_file_id;

                        while (true) {

                            {
                                size_t buffer_sz = 0;

                                std::unique_lock<std::mutex> lock(mtx_files_in);

                                if (stop) break; // Exit loop

                                while (buffer_sz < thread_seq_buf_sz){

                                    stop = !fp.read(s, file_id);

                                    if (!stop) {

                                        buffer_sz += s.length();

                                        buffer_seq.push_back(std::move(s));

                                        if (files_as_queries) buffer_file_id.push_back(file_id);
                                        else buffer_name.push_back(std::string(fp.getNameString()));
                                    }
                                    else break;
                                }

                                if (files_as_queries) {

                                    std::unique_lock<std::mutex> lock(mtx_file_id);

                                    std::pair<typename std::unordered_map<size_t, ResultFileQuery>::iterator, bool> p_it_um_file_id;

                                    size_t prev_file_id_local = 0xffffffffffffffffULL;

                                    for (const auto file_id_local : buffer_file_id) {

                                        if (file_id_local != prev_file_id_local) {

                                            p_it_um_file_id = um_file_id.insert(std::pair<size_t, ResultFileQuery>(file_id_local, ResultFileQuery()));
                                            prev_file_id_local = file_id_local;

                                            if (p_it_um_file_id.first->second.color_occ == nullptr) p_it_um_file_id.first->second.color_occ = new uint32_t[nb_colors]();
                                        }

                                        p_it_um_file_id.first->second.nb_queries += 1;

                                        if ((file_id_local != prev_file_id) && (prev_file_id != 0xffffffffffffffffULL)) {

                                            typename std::unordered_map<size_t, ResultFileQuery>::iterator it_um_file_id = um_file_id.find(prev_file_id);

                                            if (it_um_file_id == um_file_id.end()) {

                                                query_success = false;
                                                break;
                                            }
                                            else it_um_file_id->second.is_read = true;
                                        }

                                        prev_file_id = file_id_local;
                                    }

                                    if (query_success && stop && (prev_file_id != 0xffffffffffffffffULL)) { // This thread is the last one to read from input

                                        typename std::unordered_map<size_t, ResultFileQuery>::iterator it_um_file_id = um_file_id.find(prev_file_id);

                                        if (it_um_file_id == um_file_id.end()) {

                                            query_success = false;
                                            break;
                                        }
                                        else it_um_file_id->second.is_read = true;
                                    }
                                }
                            }

                            size_t l_nb_queries_found = 0;
                            size_t l_nb_queries_processed = 0;

                            std::vector<std::pair<size_t, ResultFileQuery>> v_rfq_out;

                            for (size_t i = 0; i < buffer_seq.size(); ++i){

                                const size_t nb_km_query = buffer_seq[i].length() - k + 1;

                                for (auto& c : buffer_seq[i]) c &= 0xDF;

                                if (files_as_queries){

                                    searchQuery(buffer_seq[i], color_occ_r, color_occ_u, nb_km_query);

                                    {
                                        std::unique_lock<std::mutex> lock(mtx_file_id);

                                        typename std::unordered_map<size_t, ResultFileQuery>::iterator it_um_file_id = um_file_id.find(buffer_file_id[i]);

                                        if (it_um_file_id == um_file_id.end()) {

                                            query_success = false;
                                            break;
                                        }
                                        else {

                                            ResultFileQuery& rfq = it_um_file_id->second;

                                            rfq.nb_queries -= 1;
                                            rfq.nb_km_queries += nb_km_query;

                                            for (size_t j = 0; j < nb_colors; ++j) rfq.color_occ[j] += color_occ_u[j];

                                            if (rfq.is_read && (rfq.nb_queries == 0)) { // All records for this file have been 1/ read from input 2/ queried

                                                v_rfq_out.push_back(std::pair<size_t, ResultFileQuery>(buffer_file_id[i], rfq)); // Makes a copy
                                                rfq.color_occ = nullptr;
                                                um_file_id.erase(it_um_file_id);
                                            }
                                        }
                                    }
                                }
                                else {

                                    const size_t nb_km_min = std::max(static_cast<size_t>(1), static_cast<size_t>(std::round(static_cast<double>(nb_km_query) * min_ratio_kmers)));

                                    searchQuery(buffer_seq[i], color_occ_r, color_occ_u, nb_km_min);

                                    const bool is_found = writeBinaryOutput(buffer_name[i].c_str(), buffer_name[i].length(), color_occ_u, buffer_res, pos_buffer_out, nb_km_min, mtx_file_out);

                                    l_nb_queries_found += static_cast<size_t>(is_found);
                                    l_nb_queries_processed += 1;
                                }

                                // Clean up
                                {
                                    std::memset(color_occ_u, 0, nb_colors * sizeof(uint32_t));

                                    if (inexact_search){

                                        for (size_t j = 0; j < nb_colors; ++j) color_occ_r[j] = Roaring(); // Reset k-mer occurences for each color
                                    }
                                }
                            }

                            {
                                for (const auto& p : v_rfq_out) {

                                    const ResultFileQuery& rfq = p.second;
                                    const std::string& q_name = query_filenames[p.first];
                                    const size_t nb_km_min = std::max(static_cast<size_t>(1), static_cast<size_t>(std::round(static_cast<double>(rfq.nb_km_queries) * min_ratio_kmers)));
                                    const bool is_found = writeBinaryOutput(q_name.c_str(), q_name.length(), rfq.color_occ, buffer_res, pos_buffer_out, nb_km_min, mtx_file_out);

                                    delete[] rfq.color_occ;

                                    l_nb_queries_found += static_cast<size_t>(is_found);
                                    ++l_nb_queries_processed;
                                }

                                v_rfq_out.clear();
                            }

                            nb_queries_found += l_nb_queries_found;
                            nb_queries_processed += l_nb_queries_processed;

                            // Clear buffers for next round
                            buffer_seq.clear();
                            buffer_name.clear();
                            buffer_file_id.clear();
                        }

                        // Flush rest of the thread buffer to output
                        if (write_success && (pos_buffer_out > 0)) {

                            std::unique_lock<std::mutex> lock(mtx_file_out);

                            out.write(buffer_res, pos_buffer_out);
                        }

                        delete[] buffer_res;
                        delete[] color_occ_u;

                        if (color_occ_r != nullptr) delete[] color_occ_r;
                    }
                );
            }

            for (auto& t : workers) t.join();

            if (files_as_queries && !um_file_id.empty()) {

                // This (rarely) happens when a thread overflows its input buffer just after reading the last record.
                // At that point, the thread doesn't know yet that the last record was read ("stop" is still false)
                // so the last file hasn't been annotated has fully read ("is_read" still false).
                // Subsequent threads attempting to read from input will set is_read=true but will not get
                // any additional sequences to search, hence the result for that file never gets pushed to output

                size_t pos_buffer_out = 0;

                char* buffer_res = new char[thread_seq_buf_sz];

                for (auto& p : um_file_id) {

                    ResultFileQuery& rfq = p.second;

                    if (!rfq.is_read || (rfq.nb_queries != 0)) {

                        query_success = false;
                        break;
                    }
                    else {

                        const std::string& q_name = query_filenames[p.first];
                        const size_t nb_km_min = std::max(static_cast<size_t>(1), static_cast<size_t>(std::round(static_cast<double>(rfq.nb_km_queries) * min_ratio_kmers)));
                        const bool is_found = writeBinaryOutput(q_name.c_str(), q_name.length(), rfq.color_occ, buffer_res, pos_buffer_out, nb_km_min, mtx_file_out);

                        delete[] rfq.color_occ;

                        nb_queries_found += static_cast<size_t>(is_found);
                        ++nb_queries_processed;
                    }
                }

                // Flush rest of the thread buffer to output
                if (write_success && (pos_buffer_out > 0)) {

                    std::unique_lock<std::mutex> lock(mtx_file_out);

                    out.write(buffer_res, pos_buffer_out);
                }

                delete[] buffer_res;
            }

            if (write_success && query_success && verbose) {

                std::cout << "ColoredCDBG::searchMinRatioKmer(): Processed " << nb_queries_processed << " queries." << std::endl;
                std::cout << "ColoredCDBG::searchMinRatioKmer(): Found " << nb_queries_found << " queries in at least " << ((min_nb_colors == 0) ? 1 : min_nb_colors) << " color(s)." << std::endl;
            }
        }

        fp.close();
    }

    if (!query_success) std::cerr << "ColoredCDBG::searchMinRatioKmer(): Unexpected error encountered. Please file an issue. Operation aborted." << std::endl;
    if (!write_success) std::cerr << "ColoredCDBG::searchMinRatioKmer(): Output stream is in a failed state and cannot be written to. Operation aborted." << std::endl;

    return write_success && query_success;
}

template<typename U>
bool ColoredCDBG<U>::search(const std::vector<std::string>& query_filenames, std::ostream& out,
                                    const bool found_km_ratio_out, const bool inexact_search,
                                    const bool files_as_queries, const size_t nb_threads, const bool verbose) const {

    if (invalid){

        std::cerr << "ColoredCDBG::search(): Graph is invalid and cannot be searched" << std::endl;
        return false;
    }

    if (nb_threads > std::thread::hardware_concurrency()){

        std::cerr << "ColoredCDBG::search(): Number of threads cannot be greater than or equal to " << std::thread::hardware_concurrency() << "." << std::endl;
        return false;
    }

    if (nb_threads <= 0){

        std::cerr << "ColoredCDBG::search(): Number of threads cannot be less than or equal to 0." << std::endl;
        return false;
    }

    if (out.fail()) {

        std::cerr << "CompactedDBG::search(): Output stream is in a failed state and cannot be written to." << std::endl;
        return false;
    }

    if (verbose) std::cout << "ColoredCDBG::search(): Querying graph." << std::endl;

    const size_t k = this->getK();
    const size_t thread_seq_buf_sz = BUFFER_SIZE;
    const size_t nb_colors = getNbColors();

    const char eol = '\n';

    std::string s;

    bool write_success = true;
    bool query_success = true;

    size_t file_id = 0;
    size_t prev_file_id = 0xffffffffffffffffULL;

    auto processCounts = [&](const std::vector<std::pair<size_t, const_UnitigColorMap<U>>>& v_um, Roaring* color_occ_r, uint32_t* color_occ_u){

        struct hash_pair {

            size_t operator() (const std::pair<size_t, std::pair<Kmer, size_t>>& p) const {

                return wyhash(&p, sizeof(std::pair<size_t, std::pair<Kmer, size_t>>), 0, _wyp);
            }
        };

        std::unordered_set<std::pair<size_t, std::pair<Kmer, size_t>>, hash_pair> s_um;

        typename std::unordered_set<std::pair<size_t, std::pair<Kmer, size_t>>, hash_pair>::const_iterator it;

        for (const auto& p : v_um){

            s_um.insert({p.first, {p.second.strand ? p.second.getUnitigHead() : p.second.getUnitigTail().twin(), p.second.dist}});
        }

        for (const auto& p : v_um){

            const_UnitigColorMap<U> um = p.second;

            size_t pos_query = p.first;

            if (um.strand) {

                const Kmer head = um.getUnitigHead();

                size_t pos_unitig = um.dist;

                it = s_um.find({pos_query, {head, um.dist}});

                if (it != s_um.end()) {

                    s_um.erase(it);

                    while ((pos_unitig + k) < um.size){

                        ++pos_query;
                        ++pos_unitig;

                        it = s_um.find({pos_query, {head, pos_unitig}});

                        if (it == s_um.end()) break;
                        else {

                            ++(um.len);

                            s_um.erase(it);
                        }
                    }

                    const UnitigColors* uc = um.getData()->getUnitigColors(um);

                    UnitigColors::const_iterator it_uc = uc->begin(um);
                    UnitigColors::const_iterator it_uc_end = uc->end();

                    if (inexact_search){

                        for (; it_uc != it_uc_end; ++it_uc) color_occ_r[it_uc.getColorID()].add(it_uc.getKmerPosition() - um.dist + p.first);
                    }
                    else {

                        for (; it_uc != it_uc_end; ++it_uc) color_occ_u[it_uc.getColorID()] += 1;
                    }
                }
            }
            else {

                const Kmer head = um.getUnitigTail().twin();

                it = s_um.find({pos_query, {head, um.dist}});

                if (it != s_um.end()) {

                    s_um.erase(it);

                    while (um.dist > 0){

                        ++pos_query;

                        it = s_um.find({pos_query, {head, um.dist - 1}});

                        if (it == s_um.end()) break;
                        else {

                            --(um.dist);
                            ++(um.len);

                            s_um.erase(it);
                        }
                    }

                    const UnitigColors* uc = um.getData()->getUnitigColors(um);

                    UnitigColors::const_iterator it_uc = uc->begin(um);
                    UnitigColors::const_iterator it_uc_end = uc->end();

                    if (inexact_search){

                        const size_t max_pos_um = um.dist + um.len - 1;

                        for (; it_uc != it_uc_end; ++it_uc) color_occ_r[it_uc.getColorID()].add(max_pos_um - it_uc.getKmerPosition() + p.first);
                    }
                    else {

                        for (; it_uc != it_uc_end; ++it_uc) color_occ_u[it_uc.getColorID()] += 1;
                    }
                }
            }
        }

        if (inexact_search){

            for (size_t i = 0; i < nb_colors; ++i) color_occ_u[i] = color_occ_r[i].cardinality();
        }
    };

    auto searchQuery = [&](const std::string& query, Roaring* color_occ_r, uint32_t* color_occ_u, const size_t nb_km_min){

        const std::vector<std::pair<size_t, const_UnitigColorMap<U>>> v_um_e = this->searchSequence(query, true, false, false, false, false);

        processCounts(v_um_e, color_occ_r, color_occ_u); // Extract k-mer occurrences for each color

        if (inexact_search){

            const std::vector<std::pair<size_t, const_UnitigColorMap<U>>> v_um_d = this->searchSequence(query, false, false, true, false, false);

            processCounts(v_um_d, color_occ_r, color_occ_u); // Extract k-mer occurrences for each color

            const std::vector<std::pair<size_t, const_UnitigColorMap<U>>> v_um_m = this->searchSequence(query, false, false, false, true, false);

            processCounts(v_um_m, color_occ_r, color_occ_u); // Extract k-mer occurrences for each color

            const std::vector<std::pair<size_t, const_UnitigColorMap<U>>> v_um_i = this->searchSequence(query, false, true, false, false, false);

            processCounts(v_um_i, color_occ_r, color_occ_u); // Extract k-mer occurrences for each color
        }
    };

    auto writeQuantifiedOutput = [&](   const char* query_name, const size_t len_query_name, const size_t nb_km_query,
                                        const uint32_t* color_occ, char* buffer_res, size_t& pos_buffer_out, std::mutex& mtx){

        std::string color_query_out = "";

        for (size_t i = 0; i < nb_colors; ++i) {

            color_query_out += '\t';

            if (!found_km_ratio_out) color_query_out += std::to_string(color_occ[i]);
            else color_query_out += std::to_string(static_cast<double>(color_occ[i]) / static_cast<double>(nb_km_query));
        }

        const size_t l_color_query_out = color_query_out.length();

        if ((pos_buffer_out + len_query_name + l_color_query_out + 1) > thread_seq_buf_sz){

            std::unique_lock<std::mutex> lock(mtx);

            if (pos_buffer_out > 0) {

                out.write(buffer_res, pos_buffer_out); // Write result buffer
                pos_buffer_out = 0; // Reset position to 0;
            }

            out.write(query_name, len_query_name * sizeof(char)); // Write query name
            out.write(color_query_out.c_str(), l_color_query_out * sizeof(char)); // Write query name
            out.write(&eol, sizeof(char));

            write_success = (write_success && !out.fail());
        }
        else {

            // Copy new result to buffer
            std::memcpy(buffer_res + pos_buffer_out, query_name, len_query_name * sizeof(char));
            std::memcpy(buffer_res + pos_buffer_out + len_query_name, color_query_out.c_str(), l_color_query_out * sizeof(char));

            pos_buffer_out += len_query_name + l_color_query_out;

            buffer_res[pos_buffer_out++] = eol;
        }
    };

    // Write header to TSV file
    if (write_success) {

        const std::vector<std::string> color_names = getColorNames();

        out << "query_name";

        for (const auto& name : color_names) out << '\t' << name;

        out << '\n';

        write_success = (write_success && !out.fail());
    }

    if (write_success) {

        FileParser fp(query_filenames);

        if (nb_threads == 1){

            const char* query_name = nullptr;

            char* buffer_res = new char[thread_seq_buf_sz];

            uint32_t* color_occ_u = new uint32_t[nb_colors]();
            uint32_t* color_occ_u_fid = files_as_queries ? new uint32_t[nb_colors]() : nullptr;

            Roaring* color_occ_r = inexact_search ? new Roaring[nb_colors] : nullptr;

            std::mutex mtx_file_out; // Dummy

            size_t pos_buffer_out = 0;
            size_t nb_queries_processed = 0;

            size_t nb_km_query = 0;

            while (write_success && fp.read(s, file_id)){

                if (files_as_queries) {

                    for (size_t j = 0; j < nb_colors; ++j) color_occ_u_fid[j] += color_occ_u[j];

                    if (file_id != prev_file_id) {

                        if (prev_file_id != 0xffffffffffffffffULL) { // Push results to buffer, write buffer if overflow

                            writeQuantifiedOutput(query_name, strlen(query_name), nb_km_query, color_occ_u_fid, buffer_res, pos_buffer_out, mtx_file_out);

                            ++nb_queries_processed;
                        }

                        query_name = query_filenames[file_id].c_str(); // Query name is the filename
                        nb_km_query = 0;

                        std::memset(color_occ_u_fid, 0, nb_colors * sizeof(uint32_t));
                    }

                    nb_km_query += s.length() - k + 1;
                }
                else {

                    // Push results to buffer, write buffer if overflow
                    if (prev_file_id != 0xffffffffffffffffULL) {

                        writeQuantifiedOutput(query_name, strlen(query_name), nb_km_query, color_occ_u, buffer_res, pos_buffer_out, mtx_file_out);

                        ++nb_queries_processed;
                    }

                    query_name = fp.getNameString(); // Query name is the record name
                    nb_km_query = s.length() - k + 1;
                }

                // Clean-up
                {
                    std::memset(color_occ_u, 0, nb_colors * sizeof(uint32_t));

                    if (inexact_search){

                        for (size_t j = 0; j < nb_colors; ++j) color_occ_r[j] = Roaring(); // Reset k-mer occurences for each color
                    }
                }

                for (auto& c : s) c &= 0xDF;

                searchQuery(s, color_occ_r, color_occ_u, s.length() - k + 1);

                prev_file_id = file_id;
            }

            // Flush rest of buffer result to final output
            if (write_success && (prev_file_id != 0xffffffffffffffffULL)) {

                if (files_as_queries) {

                    for (size_t j = 0; j < nb_colors; ++j) color_occ_u_fid[j] += color_occ_u[j];

                    writeQuantifiedOutput(query_name, strlen(query_name), nb_km_query, color_occ_u_fid, buffer_res, pos_buffer_out, mtx_file_out);
                }
                else writeQuantifiedOutput(query_name, strlen(query_name), nb_km_query, color_occ_u, buffer_res, pos_buffer_out, mtx_file_out);

                ++nb_queries_processed;

                if (pos_buffer_out > 0) out.write(buffer_res, pos_buffer_out);
            }

            delete[] buffer_res;
            delete[] color_occ_u;

            if (color_occ_r != nullptr) delete[] color_occ_r;
            if (color_occ_u_fid != nullptr) delete[] color_occ_u_fid;

            if (write_success && verbose) std::cout << "ColoredCDBG::search(): Processed " << nb_queries_processed << " queries. " << std::endl;
        }
        else {

            struct ResultFileQuery {

                uint32_t* color_occ;
                size_t nb_km_queries;
                size_t nb_queries;
                bool is_read;

                ResultFileQuery() : color_occ(nullptr), nb_km_queries(0), nb_queries(0), is_read(false) {}
            };

            bool stop = false;

            std::vector<std::thread> workers; // need to keep track of threads so we can join them

            std::atomic<size_t> nb_queries_processed;

            std::mutex mtx_files_in, mtx_file_out, mtx_file_id;

            std::unordered_map<size_t, ResultFileQuery> um_file_id;

            nb_queries_processed = 0;

            for (size_t t = 0; t < nb_threads; ++t){

                workers.emplace_back(

                    [&]{

                        size_t pos_buffer_out = 0;

                        char* buffer_res = new char[thread_seq_buf_sz];
                        uint32_t* color_occ_u = new uint32_t[nb_colors]();
                        Roaring* color_occ_r = inexact_search ? new Roaring[nb_colors] : nullptr;

                        std::vector<std::string> buffer_seq;
                        std::vector<std::string> buffer_name;
                        std::vector<size_t> buffer_file_id;

                        while (true) {

                            {
                                size_t buffer_sz = 0;

                                std::unique_lock<std::mutex> lock(mtx_files_in);

                                if (stop) break; // Exit loop

                                while (buffer_sz < thread_seq_buf_sz){

                                    stop = !fp.read(s, file_id);

                                    if (!stop) {

                                        buffer_sz += s.length();

                                        buffer_seq.push_back(std::move(s));

                                        if (files_as_queries) buffer_file_id.push_back(file_id);
                                        else buffer_name.push_back(std::string(fp.getNameString()));
                                    }
                                    else break;
                                }

                                if (files_as_queries) {

                                    std::unique_lock<std::mutex> lock(mtx_file_id);

                                    std::pair<typename std::unordered_map<size_t, ResultFileQuery>::iterator, bool> p_it_um_file_id;

                                    size_t prev_file_id_local = 0xffffffffffffffffULL;

                                    for (const auto file_id_local : buffer_file_id) {

                                        if (file_id_local != prev_file_id_local) {

                                            p_it_um_file_id = um_file_id.insert(std::pair<size_t, ResultFileQuery>(file_id_local, ResultFileQuery()));
                                            prev_file_id_local = file_id_local;

                                            if (p_it_um_file_id.first->second.color_occ == nullptr) p_it_um_file_id.first->second.color_occ = new uint32_t[nb_colors]();
                                        }

                                        p_it_um_file_id.first->second.nb_queries += 1;

                                        if ((file_id_local != prev_file_id) && (prev_file_id != 0xffffffffffffffffULL)) {

                                            typename std::unordered_map<size_t, ResultFileQuery>::iterator it_um_file_id = um_file_id.find(prev_file_id);

                                            if (it_um_file_id == um_file_id.end()) {

                                                query_success = false;
                                                break;
                                            }
                                            else it_um_file_id->second.is_read = true;
                                        }

                                        prev_file_id = file_id_local;
                                    }

                                    if (query_success && stop && (prev_file_id != 0xffffffffffffffffULL)) { // This thread is the last one to read from input

                                        typename std::unordered_map<size_t, ResultFileQuery>::iterator it_um_file_id = um_file_id.find(prev_file_id);

                                        if (it_um_file_id == um_file_id.end()) {

                                            query_success = false;
                                            break;
                                        }
                                        else it_um_file_id->second.is_read = true;
                                    }
                                }
                            }

                            size_t l_nb_queries_found = 0;
                            size_t l_nb_queries_processed = 0;

                            std::vector<std::pair<size_t, ResultFileQuery>> v_rfq_out;

                            for (size_t i = 0; i < buffer_seq.size(); ++i) {

                                const size_t nb_km_query = buffer_seq[i].length() - k + 1;

                                for (auto& c : buffer_seq[i]) c &= 0xDF;

                                searchQuery(buffer_seq[i], color_occ_r, color_occ_u, nb_km_query);

                                if (files_as_queries) {

                                    std::unique_lock<std::mutex> lock(mtx_file_id);

                                    typename std::unordered_map<size_t, ResultFileQuery>::iterator it_um_file_id = um_file_id.find(buffer_file_id[i]);

                                    if (it_um_file_id == um_file_id.end()) {

                                        query_success = false;
                                        break;
                                    }
                                    else {

                                        ResultFileQuery& rfq = it_um_file_id->second;

                                        rfq.nb_queries -= 1;
                                        rfq.nb_km_queries += nb_km_query;

                                        for (size_t j = 0; j < nb_colors; ++j) rfq.color_occ[j] += color_occ_u[j];

                                        if (rfq.is_read && (rfq.nb_queries == 0)) { // All records for this file have been 1/ read from input 2/ queried

                                            v_rfq_out.push_back(std::pair<size_t, ResultFileQuery>(buffer_file_id[i], rfq)); // Makes a copy
                                            rfq.color_occ = nullptr;
                                            um_file_id.erase(it_um_file_id);
                                        }
                                    }
                                }
                                else {

                                    writeQuantifiedOutput(buffer_name[i].c_str(), buffer_name[i].length(), nb_km_query, color_occ_u, buffer_res, pos_buffer_out, mtx_file_out);

                                    ++l_nb_queries_processed;
                                }

                                // Clean up
                                {
                                    std::memset(color_occ_u, 0, nb_colors * sizeof(uint32_t));

                                    if (inexact_search){

                                        for (size_t j = 0; j < nb_colors; ++j) color_occ_r[j] = Roaring(); // Reset k-mer occurences for each color
                                    }
                                }
                            }

                            {
                                for (const auto& p : v_rfq_out) {

                                    const ResultFileQuery& rfq = p.second;
                                    const std::string& q_name = query_filenames[p.first];

                                    writeQuantifiedOutput(q_name.c_str(), q_name.length(), rfq.nb_km_queries, rfq.color_occ, buffer_res, pos_buffer_out, mtx_file_out);

                                    delete[] rfq.color_occ;

                                    ++l_nb_queries_processed;
                                }

                                v_rfq_out.clear();
                            }

                            nb_queries_processed += l_nb_queries_processed;

                            // Clear buffers for next round
                            buffer_seq.clear();
                            buffer_name.clear();
                            buffer_file_id.clear();
                        }

                        // Flush rest of the thread buffer to output
                        if (write_success && (pos_buffer_out > 0)) {

                            std::unique_lock<std::mutex> lock(mtx_file_out);

                            out.write(buffer_res, pos_buffer_out);
                        }

                        delete[] buffer_res;
                        delete[] color_occ_u;

                        if (color_occ_r != nullptr) delete[] color_occ_r;
                    }
                );
            }

            for (auto& t : workers) t.join();

            if (files_as_queries && !um_file_id.empty()) {

                // This (rarely) happens when a thread overflows its input buffer just after reading the last record.
                // At that point, the thread doesn't know yet that the last record was read ("stop" is still false)
                // so the last file hasn't been annotated has fully read ("is_read" still false).
                // Subsequent threads attempting to read from input will set is_read=true but will not get
                // any additional sequences to search, hence the result for that file never gets pushed to output

                size_t pos_buffer_out = 0;

                char* buffer_res = new char[thread_seq_buf_sz];

                for (auto& p : um_file_id) {

                    ResultFileQuery& rfq = p.second;

                    if (!rfq.is_read || (rfq.nb_queries != 0)) {

                        query_success = false;
                        break;
                    }
                    else {

                        const std::string& q_name = query_filenames[p.first];

                        writeQuantifiedOutput(q_name.c_str(), q_name.length(), rfq.nb_km_queries, rfq.color_occ, buffer_res, pos_buffer_out, mtx_file_out);

                        delete[] rfq.color_occ;

                        ++nb_queries_processed;
                    }
                }

                // Flush rest of the thread buffer to output
                if (write_success && (pos_buffer_out > 0)) {

                    std::unique_lock<std::mutex> lock(mtx_file_out);

                    out.write(buffer_res, pos_buffer_out);
                }

                delete[] buffer_res;
            }

            if (write_success && query_success && verbose) std::cout << "ColoredCDBG::search(): Processed " << nb_queries_processed << " queries. " << std::endl;
        }

        fp.close();
    }

    if (!query_success) std::cerr << "ColoredCDBG::search(): Unexpected error encountered. Please file an issue. Operation aborted." << std::endl;
    if (!write_success) std::cerr << "ColoredCDBG::search(): Output stream is in a failed state and cannot be written to. Operation aborted." << std::endl;

    return write_success && query_success;
}

template<typename U>
void ColoredCDBG<U>::checkColors(const std::vector<std::string>& filename_seq_in) const {

    std::cout << "ColoredCDBG::checkColors(): Start" << std::endl;

    size_t file_id = 0;

    std::string s;

    KmerHashTable<tiny_vector<size_t, 1>> km_h;

    FastqFile FQ(filename_seq_in);

    while (FQ.read_next(s, file_id) >= 0){

        for (KmerIterator it_km(s.c_str()), it_km_end; it_km != it_km_end; ++it_km) {

            std::pair<KmerHashTable<tiny_vector<size_t, 1>>::iterator, bool> it = km_h.insert(it_km->first.rep(), tiny_vector<size_t, 1>());

            tiny_vector<size_t, 1>& tv = *(it.first);

            const size_t id = file_id / 64;

            while (tv.size() < (id + 1)) tv.push_back(0);

            tv[id] |= (1ULL << (file_id % 64));
        }
    }

    FQ.close();

    std::cout << "ColoredCDBG::checkColors(): All k-mers in the hash table with their colors" << std::endl;

    for (typename KmerHashTable<tiny_vector<size_t, 1>>::const_iterator it_km = km_h.begin(), it_km_end = km_h.end(); it_km != it_km_end; ++it_km){

        const Kmer km = it_km.getKey();
        const const_UnitigColorMap<U> ucm = this->find(km);

        if (ucm.isEmpty){

            std::cerr << "ColoredCDBG::checkColors(): K-mer " << km.toString() << " is not found in the graph" << std::endl;
            exit(1);
        }

        const UnitigColors* cs = ucm.getData()->getUnitigColors(ucm);

        if (cs == nullptr){

            std::cerr << "ColoredCDBG::checkColors(): K-mer " << km.toString() << " has no color set associated" << std::endl;
            exit(1);
        }

        const tiny_vector<size_t, 1>& tv = *it_km;
        const size_t tv_nb_max_elem = tv.size() * 64;

        for (size_t i = 0; i < std::min(filename_seq_in.size(), tv_nb_max_elem); ++i){

            const bool color_pres_graph = cs->contains(ucm, i);
            const bool color_pres_hasht = ((tv[i/64] >> (i%64)) & 0x1) == 0x1;

            if (color_pres_graph != color_pres_hasht){

                std::cerr << "ColoredCDBG::checkColors(): Current color is " << i << ": " << filename_seq_in[i] << std::endl;
                std::cerr << "ColoredCDBG::checkColors(): K-mer " << km.toString() << " for color " << i << ": " << filename_seq_in[i] << std::endl;
                std::cerr << "ColoredCDBG::checkColors(): Size unitig: " << ucm.size << std::endl;
                std::cerr << "ColoredCDBG::checkColors(): Mapping position: " << ucm.dist << std::endl;
                std::cerr << "ColoredCDBG::checkColors(): Mapping strand: " << ucm.strand << std::endl;
                std::cerr << "ColoredCDBG::checkColors(): Present in graph: " << color_pres_graph << std::endl;
                std::cerr << "ColoredCDBG::checkColors(): Present in hash table: " << color_pres_hasht << std::endl;

                exit(1);
            }
        }
    }

    std::cout << "ColoredCDBG::checkColors(): Checked all colors of all k-mers: everything is fine" << std::endl;
    std::cout << "ColoredCDBG::checkColors(): Number of k-mers in the graph: " << km_h.size() << std::endl;
}

#endif
