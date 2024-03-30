//
// Created by Parsa Eskandar on 3/4/24.
//
#include "r-index/internal/r_index.hpp"
#include "r-index/internal/utils.hpp"
#include <random>
#include <string>
#include <vector>
#include <iostream>
#include <sys/stat.h>
#include "subcommand.hpp"
#include <spdlog/spdlog.h>
#include <vg/io/vpkg.hpp>
#include <gbwtgraph/gbwtgraph.h>
#include "../gbwtgraph_helper.hpp"
#include "../gbwt_helper.hpp"
#include "../index_registry.hpp"
#include <gbwtgraph/index.h>


#include <string>

using namespace std;
using namespace ri;
using namespace vg;
using namespace vg::subcommand;
using namespace gbwtgraph;

// This function is a version of canonical_kmers function from gbwtgraph library but in this version we separate the kmer
// and its reverse complement. So this function returns the sorted vector of kmers
template<class KeyType>
std::vector<Kmer<KeyType>>
forward_strand_kmers(std::string::const_iterator begin, std::string::const_iterator end, size_t k) {
    std::vector<Kmer<KeyType>> result;
    if (k == 0 || k > KeyType::KMER_MAX_LENGTH) {
        std::cerr << "canonical_kmers(): k must be between 1 and " << KeyType::KMER_MAX_LENGTH << std::endl;
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
forward_strand_kmers(const std::string& seq, size_t k)
{
    return forward_strand_kmers<KeyType>(seq.begin(), seq.end(), k);
}


// return an index of unique kmers and their position in the graph
// if the kmer is not unique, the special value pos_t (0, 0, 0) TODO: check this

std::unordered_map<gbwtgraph::Key64::value_type, int> keyOccurrences;
template<class KeyType>
void
unique_kmers(const GBWTGraph& graph, vg::hash_map<gbwtgraph::Key64::value_type, vg::pos_t>& index, size_t k){
    typedef KeyType key_type;
    typedef Kmer<key_type> kmer_type;

    // add the kmers to the hash_map and if the kmer is not unique, set the value to std::numeric_limits<size_t> max()
    auto hash_kmers = [&](const std::vector<handle_t>& traversal, const std::string& seq){
        // get the kmers from the forward strand of the sequence
        std::vector<kmer_type> kmers = forward_strand_kmers<key_type>(seq, k);
        auto iter = traversal.begin();
        size_t node_start = 0;
        for (auto kmer: kmers){
            if(kmer.empty()) { continue; }

            auto it = index.find(kmer.key.get_key());
            // if the kmer is not in the index, calculate the position of the kmer in the graph and add it to the index
            if(it == index.end()) {
                size_t node_length = graph.get_length(*iter);
                // find the node that contains the kmer
                while (node_start + node_length <= kmer.offset) {
                    node_start += node_length;
                    ++iter;
                    node_length = graph.get_length(*iter);
                }
                // get the position of the kmer in the graph
                vg::pos_t pos{graph.get_id(*iter), graph.get_is_reverse(*iter), kmer.offset - node_start};
                index[kmer.key.get_key()] = pos;
            } else {
                // if the kmer is already in the index, set the value to the special value pos_t (0, 0, 0)
                index[kmer.key.get_key()] = vg::pos_t(0, 0, 0);
            }

        }
    };




    // TODO: make this parallel if possible
    gbwtgraph::for_each_haplotype_window(graph, k, hash_kmers, false);
    // print all the keyoccurrences
    for (auto& i : keyOccurrences) {
        cout << i.first << " " << i.second << endl;
    }

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


int main_panindexer(int argc, char **argv) {

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
    string gbz_file = "/Users/seeskand/Documents/pangenome-index/test_data/x.giraffe.gbz";
    unique_ptr<gbwtgraph::GBZ> gbz;
    auto input = vg::io::VPKG::try_load_first<gbwtgraph::GBZ, gbwtgraph::GBWTGraph, HandleGraph>(gbz_file);
    gbz = std::move(get<0>(input));

    bool progress = true;
    auto threshold = 0;
    auto space_efficient_counting = false;


    typedef gbwtgraph::Key64::value_type kmer_type;




    hash_map<kmer_type, vg::pos_t > index;
    unique_kmers<gbwtgraph::Key64>(gbz->graph, index, 29);

    // print the index values that are not the special value pos_t (0, 0, 0)
    for (auto& i : index) {
        if (i.second != vg::pos_t(0, 0, 0)) {
            cout << i.first << " " << vg::id(i.second) << " " << vg::offset(i.second) << endl;
        }
    }




    return 0;
}


static vg::subcommand::Subcommand vg_panindexer("panindexer", "index pangenomes", PIPELINE, main_panindexer);
