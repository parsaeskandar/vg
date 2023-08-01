/*
 * Define the "vg remove_duplicate" subcommand, which remove the duplicate PCRs
 *
 */
#include <omp.h>
#include <unistd.h>
#include <getopt.h>

#include <string>
#include <vector>
#include <set>
#include <cstdint>
#include <ctime>

#include <subcommand.hpp>

#include "../algorithms/alignment_path_offsets.hpp"
#include "../multipath_alignment.hpp"
#include "../alignment.hpp"
#include "../vg.hpp"
#include <vg/io/stream.hpp>
#include <vg/io/vpkg.hpp>
#include <bdsg/overlays/overlay_helper.hpp>
#include "../stream_sorter.hpp"
#include <BooPHF.h>
#include <bitset>
#include "../hash_map.hpp"
#include "../hts_alignment_emitter.hpp"

using namespace std;
using namespace vg;
using namespace vg::subcommand;

#include <iostream>
#include <vector>
#include <algorithm>


class Custom_string_hasher {
public:


    uint64_t operator()(const string &s, uint64_t seed = 0) const {
        size_t hash = seed;
        hash_combine(hash, s);
        return hash;


    }
};

vector <pair<long long, long long>> make_coalesced_sorted_intervals(const Alignment &aln) {
    /// return a vector of pairs of id-ts of the node ids that are themselves sorted and properly coalesced
    vector<long long> sorted_node_ids;
    vector <pair<long long, long long>> intervals;
    for (const Mapping mapping: aln.path().mapping()) {
        sorted_node_ids.push_back(mapping.position().node_id());
    }

    if (sorted_node_ids.empty()) {
        return intervals;
    }
    // We have a sorted vector of node_ids
    sort(sorted_node_ids.begin(), sorted_node_ids.end());


    long long start = sorted_node_ids[0];
    long long end = sorted_node_ids[0];
    for (long long i = 1; i < sorted_node_ids.size(); i++) {
        if (sorted_node_ids[i] == end + 1) {
            end = sorted_node_ids[i];
        } else {
            intervals.push_back(make_pair(start, end));
            start = end = sorted_node_ids[i];
        }
    }

    intervals.push_back(make_pair(start, end));
    return intervals;
}

void remove_rep_elements(std::vector<long long> &v) {
    // removing repetitive elements from a vector
    std::unordered_set<long long> s;
    auto end = std::remove_if(v.begin(), v.end(),
                              [&s](long long const &i) {
                                  return !s.insert(i).second;
                              });

    v.erase(end, v.end());
}

vector <pair<long long, long long>> make_coalesced_sorted_intervals_batch(vector <Alignment> batch) {
    /// return a vector of pairs of id-ts of the node ids that are themselves sorted and properly coalesced
    vector<long long> sorted_node_ids;
    vector <pair<long long, long long>> intervals;
    for (const Alignment alignment: batch) {
        for (const Mapping mapping: alignment.path().mapping()) {
            sorted_node_ids.push_back(mapping.position().node_id());
        }
    }
    if (sorted_node_ids.empty()) {
        return intervals;
    }

    remove_rep_elements(sorted_node_ids);

    // We have a sorted vector of node_ids
    sort(sorted_node_ids.begin(), sorted_node_ids.end());


    long long start = sorted_node_ids[0];
    long long end = sorted_node_ids[0];
    for (long long i = 1; i < sorted_node_ids.size(); i++) {
        if (sorted_node_ids[i] == end + 1) {
            end = sorted_node_ids[i];
        } else {
            intervals.push_back(make_pair(start, end));
            start = end = sorted_node_ids[i];
        }
    }

    intervals.push_back(make_pair(start, end));
    return intervals;
}


vector <pos_t> alignment_position(const Alignment &aln) {
    /// return a vector of pos_t variables which the index of the vector points to the read offset and the value points to the node pos_t
    vector <pos_t> positions;
    for (const Mapping &mapping: aln.path().mapping()) {
        long long node_offset = mapping.position().offset(); // This is the position of the first Edit
        long long node_id = mapping.position().node_id(); // The id of the node we are working with
        bool is_reverse = mapping.position().is_reverse();

        for (const Edit &edit: mapping.edit()) {

            // Match or Mismatch
            if (edit.from_length() == edit.to_length()) {

                for (size_t i = 0; i < edit.to_length(); ++i) {
                    positions.push_back(make_pos_t(node_id, is_reverse, node_offset + i));
                }
                node_offset += edit.from_length();

            } else if (edit.to_length() == 0 && edit.from_length() > edit.to_length()) { // Deletion
                node_offset += edit.from_length();

            } else if (edit.from_length() < edit.to_length()) { // insertion
                node_offset += edit.from_length();
                for (size_t i = 0; i < edit.to_length(); ++i) {
                    positions.push_back(make_pos_t(node_id, is_reverse, node_offset));
                }
            }


        }
    }
    return positions;

}

bool check_duplicate_batch(vector <pos_t> first_alignment_pos, vector <pos_t> second_alignment_pos) {

    if (first_alignment_pos.empty() || second_alignment_pos.empty()) {
        return false;
    }
//
//    if (id(first_alignment_pos[0]) == id(second_alignment_pos[0]) &&
//        offset(first_alignment_pos[0]) == offset(second_alignment_pos[0]) &&
//        (is_rev(first_alignment_pos[0]) == is_rev(second_alignment_pos[0]))) {


    for (size_t i = 0; i < first_alignment_pos.size(); ++i) {
        if (id(first_alignment_pos[i]) == id(second_alignment_pos[i]) &&
            offset(first_alignment_pos[i]) == offset(second_alignment_pos[i]) &&
            (is_rev(first_alignment_pos[i]) == is_rev(second_alignment_pos[i]))) {
//            size_t t = first_alignment_pos.size() - 1;
//            // We check if the alignments end on the same position
//            if (id(first_alignment_pos[t]) == id(second_alignment_pos[t]) &&
//                offset(first_alignment_pos[t]) == offset(second_alignment_pos[t]) &&
//                (is_rev(first_alignment_pos[t]) == is_rev(second_alignment_pos[t])))
            return true;
        }

    }

//    }


    return false;

}


// Assuming v1 and v2 are sorted vectors
bool have_common(const vector<long long> &v1, const vector<long long> &v2) {

    int i = 0, j = 0;

    while (i < v1.size() && j < v2.size()) {
        if ((v1[i] == v2[j])) {
            return true;
        } else if (v1[i] < v2[j]) {
            i++;
        } else {
            j++;
        }
    }

    return false;
}


bool check_duplicate_on_fly(const Alignment aln1, const Alignment aln2) {
    /// This is a function that try to open the mapping while checking if they are duplicates


    if (aln1.path().mapping_size() == 0 || aln2.path).mapping_size() == 0){
        return false;
    }

    vector<long long> aln1_mapping_nodes;
    for (const Mapping &mapping: aln1.path().mapping()) {
        aln1_mapping_nodes.push_back(mapping.position().node_id());
    }
    vector<long long> aln2_mapping_nodes;
    for (const Mapping &mapping: aln2.path().mapping()) {
        aln2_mapping_nodes.push_back(mapping.position().node_id());
    }


    // find the indexes of common nodes
    vector<int> idx1;
    vector<int> idx2;

    unordered_map<long long, int> m1;
    for (int i = 0; i < aln1_mapping_nodes.size(); i++) {
        m1[aln1_mapping_nodes[i]] = i;
    }

    for (int j = 0; j < aln2_mapping_nodes.size(); j++) {
        if (m1.count(aln2_mapping_nodes[j]) > 0) {
            idx1.push_back(m1[aln2_mapping_nodes[j]]);
            idx2.push_back(j);
        }
    }

    for (int i = 0; i < idx1.size(); i++) {
        int first_index = idx1[i];
        int second_index = idx2[i];

        // TODO: Check is_reverse equality here
        long long node_offset1 = aln1.path().mapping(first_index).position().offset();
        bool is_rev1 = aln1.path().mapping(first_index).position().is_reverse();
        long long node_offset2 = aln2.path().mapping(second_index).position().offset();
        bool is_rev2 = aln2.path().mapping(second_index).position().is_reverse();


        if (is_rev1 == is_rev2){


            vector<long long> positions1;
            vector<long long> positions2;

            for (const Edit &edit: aln1.path().mapping(first_index).edit()) {

                // Match or Mismatch
                if (edit.from_length() == edit.to_length()) {

                    for (size_t index = 0; index < edit.to_length(); ++index) {
                        positions1.push_back(node_offset1 + index);
                    }
                    node_offset1 += edit.from_length();

                } else if (edit.to_length() == 0 && edit.from_length() > edit.to_length()) { // Deletion
                    node_offset1 += edit.from_length();

                    // TODO: check if this is a true statement (For this way I don't add one per seq)
                } else if (edit.from_length() < edit.to_length()) { // insertion
                    node_offset1 += edit.from_length();
                    positions1.push_back(node_offset1);
//                for (size_t ind = 0; ind < edit.to_length(); ++ind) {
//                    positions1.push_back(node_offset1);
//                    reverse1.push_back(is_rev1);
//                }
                }
            }

            for (const Edit &edit: aln2.path().mapping(second_index).edit()) {

                // Match or Mismatch
                if (edit.from_length() == edit.to_length()) {

                    for (size_t index2 = 0; index2 < edit.to_length(); ++index2) {
                        positions2.push_back(node_offset2 + index2);
                    }
                    node_offset2 += edit.from_length();

                } else if (edit.to_length() == 0 && edit.from_length() > edit.to_length()) { // Deletion
                    node_offset2 += edit.from_length();

                    // TODO: check if this is a true statement (For this way I don't add one per seq)
                } else if (edit.from_length() < edit.to_length()) { // insertion
                    node_offset2 += edit.from_length();
                    positions2.push_back(node_offset2);
//                for (size_t ind2 = 0; ind2 < edit.to_length(); ++ind2) {
//                    positions2.push_back(node_offset2);
//                    reverse2.push_back(is_rev2);
//                }
                }
            }

            if (have_common(positions1, positions2)) {
                return true;
            }

        }


    }

    return false;

}


bool check_duplicate(const Alignment aln1, const Alignment aln2) {
    /// This function check if the two input alignments are duplicate or not, They are duplicate even if have one equal position on one base
    vector <pos_t> first_alignment_pos = alignment_position(aln1);
    vector <pos_t> second_alignment_pos = alignment_position(aln2);


    if (first_alignment_pos.empty() || second_alignment_pos.empty()) {
        return false;
    }
//
//    if (id(first_alignment_pos[0]) == id(second_alignment_pos[0]) &&
//        offset(first_alignment_pos[0]) == offset(second_alignment_pos[0]) &&
//        (is_rev(first_alignment_pos[0]) == is_rev(second_alignment_pos[0]))) {


    for (size_t i = 0; i < first_alignment_pos.size(); ++i) {
        if (id(first_alignment_pos[i]) == id(second_alignment_pos[i]) &&
            offset(first_alignment_pos[i]) == offset(second_alignment_pos[i]) &&
            (is_rev(first_alignment_pos[i]) == is_rev(second_alignment_pos[i]))) {
//            size_t t = first_alignment_pos.size() - 1;
//            // We check if the alignments end on the same position
//            if (id(first_alignment_pos[t]) == id(second_alignment_pos[t]) &&
//                offset(first_alignment_pos[t]) == offset(second_alignment_pos[t]) &&
//                (is_rev(first_alignment_pos[t]) == is_rev(second_alignment_pos[t])))
            return true;
        }

    }

//    }


    return false;


}

bool check_pair_duplicate(const Alignment aln1, const Alignment aln2) {
    /// This function check if the two input alignments are duplicate or not, They are duplicate even if have same begining and end
    vector <pos_t> first_alignment_pos = alignment_position(aln1);
    vector <pos_t> second_alignment_pos = alignment_position(aln2);
    if (id(first_alignment_pos.front()) == id(second_alignment_pos.front()) &&
        offset(first_alignment_pos.front()) == offset(second_alignment_pos.front()) &&
        (is_rev(first_alignment_pos.front()) == is_rev(second_alignment_pos.front()))) {
//        size_t t = first_alignment_pos.size() - 1;
        // We check if the alignments end on the same position
        if (id(first_alignment_pos.back()) == id(second_alignment_pos.back()) &&
            offset(first_alignment_pos.back()) == offset(second_alignment_pos.back()) &&
            (is_rev(first_alignment_pos.back()) == is_rev(second_alignment_pos.back())))
            return true;
    }

    return false;


}

string name_id(const Alignment &aln) {
    /// Make a unique name for each alignment base on the name and if they have next/prev fragment
    string id;
    if (aln.has_fragment_next()) {
        id = aln.name() + "/1";
    } else if (aln.has_fragment_prev()) {
        id = aln.name() + "/2";
    } else {
        id = aln.name();
    }
    return id;

}

vector<int> find_indexes(const std::vector<std::vector<long long>> &batch_nodes, long long node_id) {
    std::vector<int> indexes;

    for (int i = 0; i < batch_nodes.size(); i++) {
        const std::vector<long long> &vec = batch_nodes[i];

        if (std::find(vec.begin(), vec.end(), node_id) != vec.end()) {
            indexes.push_back(i);
        }
    }

    return indexes;
}


void help_rmvdup(char **argv) {
    cerr << "usage: " << argv[0] << " rmvdup [options] inputfile.gam > output.gam " << endl
         << "Remove duplicate PCRs from the input file. A gam index file (.gam.gai) must exists." << endl
         << "  -s --single_end                  set this flag if the input gam file is single-end read" << endl
         << "  -o, --output_type               prints the pairs of (sequence, name) as strings in the output" << endl
         << "    -t, --threads N            number of threads to use" << endl
         << "   -g --get_duplicates            prints the duplicate candidates in the output(can use with the -o)"
         << endl;

}


typedef boomphf::mphf <string, Custom_string_hasher> boophf_t;

int main_rmvdup(int argc, char *argv[]) {
    string filename;
    bool output_t = false;
    bool print_duplicates = false;
    int threads = 1;
    omp_set_num_threads(threads);
    bool single_ended = false;

    int c;
    optind = 2;  // force optind past command positional argument
    while (true) {
        static struct option long_options[] = {
                {"help",           no_argument,       0, 'h'},
                {"output_type",    no_argument,       0, 'o'},
                {"get_duplicates", no_argument,       0, 'g'},
                {"threads",        required_argument, 0, 't'},
                {"single_end",     no_argument,       0, 's'},
                {0,                0,                 0, 0}
        };
        int option_index = 0;
        c = getopt_long(argc, argv, "hogt:s",
                        long_options, &option_index);

        if (c == -1) break;

        switch (c) {

            case 'o':
                output_t = true;
                break;

            case 't':
                threads = parse<int>(optarg);
                omp_set_num_threads(threads);
                break;

            case 'g':
                print_duplicates = true;
                break;

            case 's':
                single_ended = true;
                break;

            case 'h':
            case '?':
                help_rmvdup(argv);
                exit(1);
                break;

            default:
                help_rmvdup(argv);
                abort();
        }
    }

    if (optind >= argc) {
        cerr << "[vg view] error: no filename given" << endl;
        exit(1);
    }


    string sorted_gam_name = get_input_file_name(optind, argc, argv);
    unique_ptr <GAMIndex> gam_index;
    if (!sorted_gam_name.empty()) {
        // Load the GAM index
        gam_index = unique_ptr<GAMIndex>(new GAMIndex());
        get_input_file(sorted_gam_name + ".gai", [&](istream &in) {
            // We get it form the appropriate .gai, which must exist
            gam_index->load(in);
        });
    }

    vector <vector<string>> fill_hash_memory(get_thread_count());
    vector <string> keys;

    function<void(Alignment & )> fill_hash = [&](const Alignment &aln) {
        int thread_number = omp_get_thread_num();
        fill_hash_memory[thread_number].push_back(name_id(aln));
    };
    get_input_file(sorted_gam_name, [&](istream &in) {
        vg::io::for_each_parallel(in, fill_hash);
    });

    // Making hashing multi-threaded

    size_t totalSize = 0;

    // Calculate the total size of all vectors
    for (const auto &vec: fill_hash_memory) {
        totalSize += fill_hash_memory.size();
    }

    // Reserve memory to avoid frequent reallocation
    keys.reserve(totalSize);

    // Merge the vectors
    for (const auto &vec: fill_hash_memory) {
        keys.insert(keys.end(), vec.begin(), vec.end());
    }


    vector<bool> checked(keys.size(), false);

    boophf_t *bphf = new boomphf::mphf<string, Custom_string_hasher>(keys.size(), keys, threads, 2.0, false, false);

//    std::unique_ptr<vg::io::ProtobufEmitter<Alignment>> emitter;
    vector <tuple<path_handle_t, size_t, size_t>> paths;
    PathPositionHandleGraph *path_position_graph = nullptr;
//    paths = get_sequence_dictionary("", {}, *path_position_graph);
    unique_ptr <AlignmentEmitter> alignment_emitter;
    alignment_emitter = get_alignment_emitter("-", "GAM",
                                              paths, get_thread_count());
//    emitter = std::unique_ptr<vg::io::ProtobufEmitter<Alignment>>(new vg::io::ProtobufEmitter<Alignment>(cout));
    function<void(Alignment & )> pcr_removal = [&](Alignment &aln) {
        if (gam_index.get() != nullptr) {
            // I mark all the reads that have to get deleted as duplicates. This means one read from each duplicate set remains unmarked.
            // This way we can remove all reads with duplicate flag and not worry about deleting them all
            if (!checked[bphf->lookup(name_id(aln))]) {
                // make the alignment nodes list that can be use as input if .find function of the gam_index
                vector <pair<long long, long long>> intervals = make_coalesced_sorted_intervals(aln);
                // Find all alignments that share at least one node with the current working alignment
                get_input_file(sorted_gam_name, [&](istream &input_gam) {
                    vg::io::ProtobufIterator<Alignment> gam_cursor(input_gam);
                    // find all sharing nodes alignments and call the function to handle the result
                    gam_index->find(gam_cursor, intervals, [&](const Alignment &share_aln) {

                        if (!checked[bphf->lookup(name_id(share_aln))]) {
                            if (name_id(aln) != name_id(share_aln)) {
                                if (check_duplicate(aln, share_aln)) {
                                    if (print_duplicates) {
//#pragma omp critical (cerr)           int a = 0;
//                                        if (output_t)
//                                            cout << share_aln.sequence() << endl;
//                                        else
//                                            emitter->write(std::move(const_cast<Alignment &>(share_aln)));
                                    }
                                    checked[bphf->lookup(name_id(share_aln))] = true;
                                }

                            }

                        }
                    });

                });
                // Writing the non-duplicate PCRs to the file
                if (!checked[bphf->lookup(name_id(aln))]) {
#pragma omp critical (cerr)
                    checked[bphf->lookup(name_id(aln))] = true;
                    if (!print_duplicates) {
//                        if (output_t)
//                            cout << aln.sequence() << "\t" << aln.name() << endl;
//                        else
//                            emitter->write(std::move(aln));
                    }
                }


            }

        }
    };


    vector <vector<Alignment>> batch(get_thread_count());
    int batch_size = 30000;
    int tmp = 0;
    function<void(Alignment & , Alignment & , bool)> batch_pair_removal = [&](Alignment &aln, Alignment &aln_pair,
                                                                              bool is_last) {
        /// This function handle the alignments in batches. It calls all the .find functions needed in one .find function

        int thread_number = omp_get_thread_num();
        if ((batch[thread_number].size() == batch_size - 2) || is_last) {

            clock_t first_time = clock();

            // is_last is for handling the last batch, which size could be less than the batch_size
            if (!is_last) {
                batch[thread_number].push_back(aln);
                batch[thread_number].push_back(aln_pair);
            } else {
                for (size_t ind = 0; ind < batch.size(); ++ind) {
                    // find the thread that still has more alignments in the batch
                    if (batch[ind].size() > 0) {
                        thread_number = ind;
                        break;
                    }
                }

            }


            vector <Alignment> output;
            vector <Alignment> filtered_batch;


            bool is_different = false;
            if (batch[thread_number].size() >= 2) {
                filtered_batch.push_back(batch[thread_number][0]);
                filtered_batch.push_back(batch[thread_number][1]);


                for (int i = 2; i < batch[thread_number].size(); i += 2) {
//                    filtered_batch.push_back(batch[thread_number][i]);
//                    filtered_batch.push_back(batch[thread_number][i + 1]);

                    Alignment first_pair, second_pair;
                    first_pair = std::move(batch[thread_number][i]);
                    second_pair = std::move(batch[thread_number][i + 1]);

//                    Alignment res_pair1, res_pair2;
//                    res_pair1 = filtered_batch[filtered_batch.size() - 2];
//                    res_pair2 = filtered_batch[filtered_batch.size() - 1];
//
//
                    if (first_pair.path().mapping_size() == 0 || second_pair.path().mapping_size() == 0) {
                        output.push_back(first_pair);
                        output.push_back(second_pair);
                    } else {
                        filtered_batch.push_back(first_pair);
                        filtered_batch.push_back(second_pair);
                    }
//                    else {
//
//                        if (!checked[bphf->lookup(name_id(first_pair))]) {
//
//                            if (!check_duplicate(first_pair, res_pair1) || !check_duplicate(second_pair, res_pair2)) {
//
//                                if (is_different) {
//
//                                    Alignment aln_dup1, aln_dup2;
//                                    aln_dup1 = std::move(batch[thread_number][i - 2]);
//                                    aln_dup2 = std::move(batch[thread_number][i - 1]);
//
//                                    if (!check_duplicate(first_pair, aln_dup1) ||
//                                        !check_duplicate(second_pair, aln_dup2)) {
//                                        filtered_batch.push_back(first_pair);
//                                        filtered_batch.push_back(second_pair);
//                                        is_different = false;
//
//                                    } else {
//                                        is_different = true;
//                                        checked[bphf->lookup(name_id(first_pair))] = true;
//                                        checked[bphf->lookup(name_id(second_pair))] = true;
//                                    }
//                                } else {
//                                    filtered_batch.push_back(first_pair);
//                                    filtered_batch.push_back(second_pair);
//                                    is_different = false;
//
//                                }
//
//                            } else {
//                                is_different = true;
//                                checked[bphf->lookup(name_id(first_pair))] = true;
//                                checked[bphf->lookup(name_id(second_pair))] = true;
//                            }
//
//                        }
//
//                    }


                }

            }


            batch[thread_number].clear();


//            batch[thread_number].push_back(aln);
//            batch[thread_number].push_back(aln_pair);





            vector <pair<long long, long long>> intervals_batch = make_coalesced_sorted_intervals_batch(
                    filtered_batch);



            // This is a map from node ids to the vector of indexes of alignments
            unordered_map<long long, vector<int>> find_results;
            vector <Alignment> find_alignments;
            get_input_file(sorted_gam_name, [&](istream &input_gam) {

                vg::io::ProtobufIterator<Alignment> gam_cursor(input_gam);


                // This finds all other alns that share node/nodes with the one of the alignments in the batch
                gam_index->find(gam_cursor, intervals_batch, [&](const Alignment &share_aln) {
                    if (!checked[bphf->lookup(name_id(share_aln))])
                        find_alignments.push_back(share_aln);

                });

            });

            vector<long long> node_ids;
            vector <vector<long long>> batch_nodes;
            for (const Alignment alignment: find_alignments) {
                vector<long long> temp;
                for (const Mapping mapping: alignment.path().mapping()) {
                    node_ids.push_back(mapping.position().node_id());
                    temp.push_back(mapping.position().node_id());


                }
                batch_nodes.push_back(temp);
            }


//            cerr << node_ids.size() << "____" << find_alignments.size() << " " << thread_number << endl;
            // this function remove duplicate node ids and keep one of them
            remove_rep_elements(node_ids);


            for (long long ids: node_ids) {
                find_results[ids] = find_indexes(batch_nodes, ids);
            }

            clock_t second_time = clock();

            if (!intervals_batch.empty()) {
                cerr << node_ids.size() << "____" << find_alignments.size() << "____________" << find_results.size()
                     << " " << thread_number << "First " << intervals_batch.front().first << "last "
                     << intervals_batch.back().second << "size of interval " << intervals_batch.size()
                     << "size of filtered batch "
                     << filtered_batch.size() << endl;


            }


            int checking = 0;
            for (int i = 0; i < filtered_batch.size(); i += 2) {
                // handle two pair alignments using the in-mem finder we have
                Alignment pair1, pair2;
                pair1 = filtered_batch[i];
                pair2 = filtered_batch[i + 1];


                if (!checked[bphf->lookup(name_id(pair1))]) {
                    checked[bphf->lookup(name_id(pair1))] = true;
                    checked[bphf->lookup(name_id(pair2))] = true;


                    unordered_set<int> cumulative_node_ids_pair1;
                    for (const Mapping mapping: pair1.path().mapping()) {
                        for (int find_index_pair1: find_results[mapping.position().node_id()]) {
                            cumulative_node_ids_pair1.insert(find_index_pair1);
                        }
                    }


                    // This is an in-mem version of .find function
//                        for (const Mapping mapping: pair1.path().mapping()) {
//                            vector<int> find_index_pair1_list = find_results[mapping.position().node_id()];
//                            for (int find_index_pair1: find_index_pair1_list) {
                    unordered_set <string> pair1_duplicates_set;
//                    vector <pos_t> pair1_alignment_pos = alignment_position(pair1);
                    for (int find_index_pair1: cumulative_node_ids_pair1) {
                        Alignment share_aln = find_alignments[find_index_pair1];
                        if (!checked[bphf->lookup(name_id(share_aln))]) {
//                            vector <pos_t> share_aln_alignment_pos = alignment_position(share_aln);
                            if (check_duplicate_on_fly(pair1, share_aln)) {
                                pair1_duplicates_set.insert(
                                        share_aln.has_fragment_prev() ? share_aln.fragment_prev().name()
                                                                      : share_aln.fragment_next().name());
                            }
                        }
                    }
//                        }




                    if (pair1_duplicates_set.size() > 1) {
//                        cerr << ++checking << " mapping size pair1 " << pair1.path().mapping_size() << " mapping size pair2 " << pair2.path().mapping_size() << " pair 1 duplicate set size " << pair1_duplicates_set.size() << " pair1_name " << pair1.sequence() << "pair2 seq " << pair2.sequence() << " --- " << thread_number << endl;
                        unordered_set<int> cumulative_node_ids;
                        for (const Mapping mapping: pair2.path().mapping()) {
                            for (int find_index_pair2: find_results[mapping.position().node_id()]) {
                                cumulative_node_ids.insert(find_index_pair2);
                            }
                        }



//                        for (const Mapping mapping: pair2.path().mapping()) {
//                            vector<int> find_index_pair2_list = find_results[mapping.position().node_id()];
//                            for (int find_index_pair2: find_index_pair2_list) {

                        // Handling check duplicate manually to make it faster!!!
//                        vector <pos_t> first_alignment_pos = alignment_position(pair2);

//                        cout << pair1_alignment_pos.size() << " " << first_alignment_pos.size() << " tt " << thread_number << endl;

                        for (int find_index_pair2: cumulative_node_ids) {

                            Alignment share_aln = find_alignments[find_index_pair2];

                            if (!checked[bphf->lookup(name_id(share_aln))]) {


//                                vector <pos_t> second_alignment_pos = alignment_position(share_aln);

                                if (check_duplicate_on_fly(pair2, share_aln)) {


                                    auto it = pair1_duplicates_set.find(share_aln.name());

                                    // The pair is in the set
                                    if (it != pair1_duplicates_set.end()) {
                                        // This means we find pairs that are duplicates with the main pairs so we check them for being duplicate
                                        pair1_duplicates_set.erase(it);


                                        // we add /1 to the prev one because we are naming the other pair
                                        string pair_name_id = share_aln.has_fragment_prev() ?
                                                              share_aln.fragment_prev().name() + "/1" :
                                                              share_aln.fragment_next().name() + "/2";
                                        checked[bphf->lookup(name_id(share_aln))] = true;
                                        checked[bphf->lookup(pair_name_id)] = true;


                                    }


                                }
                            }


                        }

                    }

                    pair1_duplicates_set.clear();
//                    }


//                        cout << name_id(aln) << endl;
//                    checked[bphf->lookup(name_id(pair1))] = true;
//                    checked[bphf->lookup(name_id(pair2))] = true;
                    if (!print_duplicates) {
                        if (output_t) {
//#pragma omp critical (cerr)
                            cout << "Pair1 " << aln.sequence() << "\t" << pair1.name() << "Pair2 "
                                 << aln_pair.sequence() << "\t" << pair2.name() << endl;
                        } else {
//                            cerr << "we reach here" << name_id(pair1) << " " << thread_number << endl;

                            // TODO: can change the emit from pairs to batches and check the speed
//                            vector <Alignment> output = {pair1, pair2};
                            output.push_back(pair1);
                            output.push_back(pair2);
//                            alignment_emitter->emit_singles(std::move(output));
//                            emitter->write(std::move(pair1));
//                            emitter->write(std::move(pair2));
                        }

                    }


                }


            }

            clock_t third_time = clock();
            double first_part = (second_time - first_time) / (double) 1000000;
            double second_part = (third_time - second_time) / (double) 1000000;

            cerr << "first part time " << first_part << " second part time " << second_part << " thread "
                 << thread_number << endl;



//            cerr << "have to clear batch " << endl;
            cerr << ++tmp << " " << thread_number << endl;
            alignment_emitter->emit_singles(std::move(output));
            // clear the batch of this thread to use it again with another batch
//            batch[thread_number].clear();
        } else {
            batch[thread_number].push_back(aln);
            batch[thread_number].push_back(aln_pair);
        }


    };


    function<void(Alignment & , Alignment & )> batch_pair_two_var = [&](Alignment &o1, Alignment &o2) {
        batch_pair_removal(o1, o2, false);
    };
    if (gam_index.get() != nullptr) {

        if (single_ended) {

            get_input_file(sorted_gam_name, [&](istream &in) {
                vg::io::for_each_parallel(in, pcr_removal);

            });
        } else {
            get_input_file(sorted_gam_name, [&](istream &in) {
                vg::io::for_each_parallel_shuffled_double(in, batch_pair_two_var);

            });


            Alignment a, b;
            batch_pair_removal(a, b, true);
        }


    }


    return 0;

}


// Register subcommand
static Subcommand vg_removeduplicate("rmvdup", "Remove duplicate PCRs from the input file", WIDGET, main_rmvdup);