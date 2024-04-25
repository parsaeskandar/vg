//
// Created by Parsa Eskandar on 3/4/24.
//
#include "r-index/internal/r_index.hpp"
#include "r-index/internal/utils.hpp"
#include <getopt.h>
#include <random>
#include <string>
#include <vector>
#include <iostream>
#include <sys/stat.h>
#include "subcommand.hpp"
//#include <spdlog/spdlog.h>
#include <vg/io/vpkg.hpp>
#include <gbwtgraph/gbwtgraph.h>
#include "../gbwtgraph_helper.hpp"
#include "../gbwt_helper.hpp"
#include "../index_registry.hpp"
#include <gbwtgraph/index.h>
#include <../BPlusTree.hpp>


#include <string>

using namespace std;
using namespace ri;
using namespace vg;
using namespace vg::subcommand;
using namespace gbwtgraph;

void help_panindexer(char **argv) {
    cerr << "usage: " << argv[0] << " panindexer [options]" << endl
         << endl
         << "options:" << endl
         << "    -g, --graph FILE              input graph file" << endl
         << "    -k, --kmer-length N            length of the kmers in the index [default 29]" << endl
         << "    -t, --threads N          number of threads to use [1]" << endl;
}

// This function is a version of canonical_kmers function from gbwtgraph library but in this version we separate the kmer
// and its reverse complement. So this function returns the sorted vector of kmers
template<class KeyType>
std::vector<Kmer<KeyType>>
forward_strand_kmers(std::string::const_iterator begin, std::string::const_iterator end, size_t k) {
    std::vector<Kmer<KeyType>> result;
    if (k == 0 || k > KeyType::KMER_MAX_LENGTH) {
        std::cerr << "forward_strand_kmers(): k must be between 1 and " << KeyType::KMER_MAX_LENGTH << std::endl;
        return result;
    }

    size_t valid_chars = 0, offset = 0;
    KeyType forward_key;
    std::string::const_iterator iter = begin;
    while (iter != end) {
        forward_key.forward(k, *iter, valid_chars);
        if (valid_chars >= k) {
            result.push_back({forward_key, forward_key.hash(), offset_type(offset - (k - 1)), false});
        }
        ++iter;
        offset++;
    }

    std::sort(result.begin(), result.end());
    return result;
}

template<class KeyType>
std::vector<Kmer<KeyType>>
forward_strand_kmers(const std::string &seq, size_t k) {
    return forward_strand_kmers<KeyType>(seq.begin(), seq.end(), k);
}


// This function returns the unique kmers in the graph and stores them in the index
template<class KeyType>
void
unique_kmers_parallel(const GBWTGraph &graph, vg::hash_map<gbwtgraph::Key64::value_type, gbwtgraph::Position> &index,
                      size_t k) {
    typedef KeyType key_type;
    typedef Kmer<key_type> kmer_type;
    constexpr size_t KMER_CACHE_SIZE = 1024;  // Adjust cache size as needed

    int threads = omp_get_max_threads();
    std::vector<std::vector<std::pair<gbwtgraph::Key64::value_type, gbwtgraph::Position>>> cache(threads);
    // make a duplicates cache to handle the duplicates faster
    hash_set<gbwtgraph::Key64::value_type> duplicates;


    // Lambda to flush the thread-local cache to the shared index
    auto flush_cache = [&](int thread_number) {
        auto current_cache = cache[thread_number];

        // Sort the cache by key
        std::sort(current_cache.begin(), current_cache.end(), [](const auto &a, const auto &b) {
            return a.first < b.first;
        });
        // Remove duplicates
        auto last = std::unique(current_cache.begin(), current_cache.end(), [](const auto &a, const auto &b) {
            return a.first == b.first;
        });
        current_cache.erase(last, current_cache.end());
#pragma omp critical
        {
            for (auto entry: current_cache) {
                if (duplicates.find(entry.first) != duplicates.end()) {
                    // Key is a known duplicate, skip it
                    continue;
                }
                auto it = index.find(entry.first);
                if (it == index.end()) {
                    index[entry.first] = entry.second;
                } else if (it->second != entry.second) {
                    // remove the kmer from the index
                    index.erase(it); // TODO: check the time complexity of this
                    duplicates.insert(entry.first);
                }
            }
        }

        cache[thread_number].clear();  // Clear the cache after flushing
    };

    // Main lambda function to process kmers
    auto hash_kmers = [&](const std::vector<handle_t> &traversal, const std::string &seq) {
        int thread_id = omp_get_thread_num();
        std::vector<kmer_type> kmers = forward_strand_kmers<key_type>(seq, k);
        auto iter = traversal.begin();
        size_t node_start = 0;

        for (auto kmer: kmers) {
            if (kmer.empty()) continue;

            size_t node_length = graph.get_length(*iter);
            while (node_start + node_length <= kmer.offset) {
                node_start += node_length;
                ++iter;
                node_length = graph.get_length(*iter);
            }

            vg::pos_t pos{graph.get_id(*iter), graph.get_is_reverse(*iter), kmer.offset - node_start};

            if (!gbwtgraph::Position::valid_offset(pos)) {
#pragma omp critical (cerr)
                {
                    std::cerr << "index_haplotypes(): Node offset " << vg::offset(pos) << " is too large" << std::endl;
                }
                std::exit(EXIT_FAILURE);
            }
            // Use thread-local cache to reduce contention on the shared index
            cache[thread_id].emplace_back(kmer.key.get_key(), gbwtgraph::Position::encode(pos));



            // Flush the cache if it reaches the size limit
            if (cache[thread_id].size() >= KMER_CACHE_SIZE) {
                flush_cache(thread_id);
            }
        }

    };

    // Parallel execution
    gbwtgraph::for_each_nonredundant_window(graph, k, hash_kmers, true);
    for (int thread_id = 0; thread_id < threads; thread_id++) { flush_cache(thread_id); }

}


// This function input is the OCC vector of the end of the sequences (#) and it returns the sorted end_of_seq vector
// which is the sorted vector of pairs (i, SA[i]) for the end of each sequence which is (ISA[j], j)
vector<pair<uint64_t, uint64_t>> sort_end_of_seq(vector<uint64_t> &OCC) {
    vector<pair<uint64_t, uint64_t>> end_of_seq;
    auto number_of_sequences = OCC.size();
    for (int64_t i = 0; i <= number_of_sequences - 1; ++i) {
        end_of_seq.emplace_back(number_of_sequences - i - 1, OCC[i]);
    }

    // Sort the end_of_seq vector by the second element of each pair
    sort(end_of_seq.begin(), end_of_seq.end(),
         [](const pair<uint64_t, uint64_t> &a, const pair<uint64_t, uint64_t> &b) {
             return a.second < b.second;
         });

    return end_of_seq;
}

// Create a Run struct to store the run data structure for the BPlusTree
struct Run {
    size_t start_position;
    size_t graph_position;
    char run_char; // the character of the run TODO: can handle this with 3 bits instead of 8 bits

    // Operators for the struct
    bool operator<(const Run &other) const {
        return start_position < other.start_position;
    }

    bool operator==(const Run &other) const {
        return start_position == other.start_position;
    }

    Run &operator=(const Run &other) {
        if (this == &other) {
            return *this;
        }
        start_position = other.start_position;
        graph_position = other.graph_position;
        run_char = other.run_char;
        return *this;
    }

    // Set the Run struct to a zero one when it is assigned to 0
    Run& operator=(int zero) {
        if (zero == 0) {  // Ensures it only responds to 0
            start_position = 0;
            graph_position = 0;
            run_char = '\0';

        }
        return *this;
    }

    // print the Run struct using cout
    friend std::ostream &operator<<(std::ostream &os, const Run &run) {
        if (run.run_char != '\0') {
            os << "start_position: " << run.start_position
               << ", graph_position: " << run.graph_position
               << ", run_char: " << run.run_char
               << '\n';
        }

        return os;
    }

};


int main_panindexer(int argc, char **argv) {

    if (argc <= 2) {
        help_panindexer(argv);
        return 1;
    }

    int threads = 1;
    size_t k = 29;
    string graph_file = argv[2];

    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
                {
                        {"graph",       required_argument, 0, 'g'},
                        {"kmer-length", required_argument, 0, 'k'},
                        {"threads",     required_argument, 0, 't'},
                        {0, 0,                             0, 0}
                };
        int option_index = 0;
        int c = getopt_long(argc, argv, "g:k:t:", long_options, &option_index);

        if (c == -1) { break; }

        switch (c) {
            case 'g':
                graph_file = optarg;
                break;
            case 'k':
                k = parse<size_t>(optarg);
                break;
            case 't':
                threads = parse<int>(optarg);
                omp_set_num_threads(threads);
                break;
            case 'h':
            case '?':
            default:
                help_panindexer(argv);
                exit(1);
        }

    }

    // read the ri index file and store it in the ri index data structure
//    string index_file = "/Users/seeskand/Documents/pangenome-index/test_data/test.txt.ri";
//    std::ifstream in(index_file);
//    bool fast;
//
//    //fast or small index?
//    in.read((char*)&fast,sizeof(fast));
//
//    string pattern = "#";
//
////    locate<r_index<> >(in, pattern);
//    r_index<> idx;
//    idx.load(in);
//
//    // get the ISA values for the end of each sequence by searching for pattern # and locating the results
//    auto OCC = idx.locate_all(pattern);
//    cout << OCC << endl;
//
//    // print the bwt of the index
//    cout << idx.get_bwt() << endl;
//
//    // using the sort_end_of_seq function
//    auto end_of_seq = sort_end_of_seq(OCC);
//
//    // print the sorted end_of_seq vector
//    for (auto& i : end_of_seq) {
//        cout << i.first << " " << i.second << endl;
//    }
//
//    // Find the first sequence end in the bwt (using end_of_seq vector)
//    auto first_seq_end = end_of_seq[0].second;
//    auto temp =idx.F_at(5);
//    // print temp
//
//
//    cout << "First sequence end in text: " << end_of_seq[0].second << " First sequence end in First column " << end_of_seq[0].first << " " << temp<< endl;


    // read the gbz file and store it in the gbz data structure

    unique_ptr<gbwtgraph::GBZ> gbz;
    auto input = vg::io::VPKG::try_load_first<gbwtgraph::GBZ, gbwtgraph::GBWTGraph, HandleGraph>(graph_file);
    gbz = std::move(get<0>(input));

    bool progress = true;
    auto threshold = 0;
    auto space_efficient_counting = false;


    typedef gbwtgraph::Key64::value_type kmer_type;

    hash_map<kmer_type, gbwtgraph::Position> index;
//    unique_kmers<gbwtgraph::Key64>(gbz->graph, index, 29);
    unique_kmers_parallel<gbwtgraph::Key64>(gbz->graph, index, k);

    BPlusTree<Run> bptree(8);
    Run R = {1, 2, 'A'};
    Run a = {5, 3, 'A'};
    bptree.insert(R, 4);
    bptree.insert(a, 3);
    // insert 10 nodes to the bptree
//    for (size_t i = 0; i < 10; i++) {
//        Run temp = {i, 5 * i + 1, 'B'};
//        bptree.insert(temp, 2 * i + 1);
//    }


    bptree.bpt_print();






//
//    // print one element of the index
//    cout << "Key " << index.begin()->first << " " << index.begin()->second.id() << "Offset " << index.begin()->second.offset() << endl;
//


    return 0;
}


static vg::subcommand::Subcommand vg_panindexer("panindexer", "index pangenomes", PIPELINE, main_panindexer);
