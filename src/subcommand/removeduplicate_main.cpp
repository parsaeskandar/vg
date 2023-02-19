/*
 * Define the "vg remove_duplicate" subcommand, which remove the duplicate PCRs
 *
 *
 */
#include <omp.h>
#include <unistd.h>
#include <getopt.h>

#include <string>
#include <vector>
#include <set>
#include <cstdint>

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

using namespace std;
using namespace vg;
using namespace vg::subcommand;

#include <iostream>
#include <vector>
#include <algorithm>


//const uint64_t SEED = 0x9747b28c;
//const uint64_t M = 0xc6a4a7935bd1e995;
//const int R = 47;
//
//uint64_t murmurHash3(const void* key, size_t len, uint64_t seed) {
//    const uint64_t *data = (const uint64_t *)key;
//    const size_t nblocks = len / 8;
//
//    uint64_t h = seed;
//
//    for (size_t i = 0; i < nblocks; i++) {
//        uint64_t k = data[i];
//        k *= M;
//        k ^= k >> R;
//        k *= M;
//        h ^= k;
//        h *= M;
//    }
//
//    const uint8_t *tail = (const uint8_t *)(data + nblocks);
//    uint64_t k1 = 0;
//
//    switch (len & 7) {
//        case 7: k1 ^= uint64_t(tail[6]) << 48;
//        case 6: k1 ^= uint64_t(tail[5]) << 40;
//        case 5: k1 ^= uint64_t(tail[4]) << 32;
//        case 4: k1 ^= uint64_t(tail[3]) << 24;
//        case 3: k1 ^= uint64_t(tail[2]) << 16;
//        case 2: k1 ^= uint64_t(tail[1]) << 8;
//        case 1: k1 ^= uint64_t(tail[0]);
//            k1 *= M;
//            k1 ^= k1 >> R;
//            k1 *= M;
//            h ^= k1;
//    }
//
//    h ^= len;
//    h ^= h >> 33;
//    h *= 0xff51afd7ed558ccd;
//    h ^= h >> 33;std::hash<string>
//    h *= 0xc4ceb9fe1a85ec53;
//    h ^= h >> 33;
//
//    return h;
//}
//
class Custom_string_hasher {
public:


    uint64_t operator() (const string& s, uint64_t seed=0) const {
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

bool check_duplicate(const Alignment aln1, const Alignment aln2) {
    /// This function check if the two input alignments are duplicate or not, They are duplicate even if have one equal position on one base
    vector <pos_t> first_alignment_pos = alignment_position(aln1);
    vector <pos_t> second_alignment_pos = alignment_position(aln2);
    for (size_t i = 0; i < first_alignment_pos.size(); ++i) {
        if (id(first_alignment_pos[i]) == id(second_alignment_pos[i]) &&
            offset(first_alignment_pos[i]) == offset(second_alignment_pos[i]) &&
            (is_rev(first_alignment_pos[i]) == is_rev(second_alignment_pos[i])))
            return true;

    }
    return false;


}


void help_rmvdup(char **argv) {
// TODO: add whatever option is needed to this list. Change long_option and getopt_long if want to add an option.
// TODO: see what input and output file formats is possible
    cerr << "usage: " << argv[0] << " rmvdup [options] inputfile.gam > output.gam " << endl
         << "Remove duplicate PCRs from the input file. A gam index file (.gam.gai) must exists." << endl
         << "  -p, --progress               Show progress." << endl
         << "    -t, --threads N            number of threads to use" << endl;

}

typedef boomphf::mphf <string, Custom_string_hasher> boophf_t;

int main_rmvdup(int argc, char *argv[]) {
    string filename;
    bool show_progress = false;
    int threads = 8;

    int c;
    optind = 2;  // force optind past command positional argument
    while (true) {
        static struct option long_options[] = {
                {"help",     no_argument,       0, 'h'},
                {"progress", no_argument,       0, 'p'},
                {"threads",  required_argument, 0, 't'},
                {0,          0,                 0, 0}
        };

        int option_index = 0;
        c = getopt_long(argc, argv, "hpt",
                        long_options, &option_index);

        if (c == -1) break;

        switch (c) {

            case 'p':
                show_progress = true;
                break;

            case 't':
                threads = parse<int>(optarg);
                omp_set_num_threads(threads);
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

    vector <string> keys;
//    int i = 0;
    function<void(Alignment & )> fill_hash = [&](const Alignment &aln) {

        keys.push_back(aln.name());
//        cout << aln.name() << ++i << endl;
    };

    get_input_file(sorted_gam_name, [&](istream &in) {
        vg::io::for_each(in, fill_hash);
    });

//    bitset<keys.size()> checked;
    cout << keys.size() << endl;
    boophf_t *bphf = new boomphf::mphf<string, Custom_string_hasher>(keys.size(), keys, threads);


    function<void(Alignment & )> test = [&](const Alignment &aln) {
        if (gam_index.get() != nullptr) {
            cout << bphf->lookup(aln.name()) << " " << aln.name() << endl;

            // This is a schema of what I am going to do
            // I mark all the reads that have to get deleted as duplicates. This means one read from each duplicate set remains unmarked.
            // This way we can remove all reads with duplicate flag and not worry about deleting them all
            // TODO: check if the above algorithm is logical
//            if (!checked.test(bphf->lookup(aln.name()))){
//                // make the alignment nodes list that can be use as input if .find function of the gam_index
//                vector <pair<long long, long long>> intervals = make_coalesced_sorted_intervals(aln);
//                // Find all alignments that share at least one node with the current working alignment
//                get_input_file(sorted_gam_name, [&](istream &input_gam) {
//                    vg::io::ProtobufIterator<Alignment> gam_cursor(input_gam);
//                    vector <Alignment> alns;
//                    // find all sharing nodes alignments and call the function to handle the result
//                    gam_index->find(gam_cursor, intervals, [&](const Alignment &share_aln) {
//
//                        if (!checked.test(bphf->lookup(share_aln.name()))){
//                        if (check_duplicate(aln, share_aln)) {
//                            checked.set(bphf->lookup(share_aln.name()));
//                            // these alignments are duplicate if we are here
//
//                        }
//                    }
//                    });
//                });
//
//
//            }


        }
    };


    if (gam_index.get() != nullptr) {
        get_input_file(sorted_gam_name, [&](istream &in) {
            vg::io::for_each(in, test);
//            vg::io::for_each_parallel(in, test);

        });
    }


    return 0;

}


// Register subcommand
static Subcommand vg_removeduplicate("rmvdup", "Remove duplicate PCRs from the input file", WIDGET, main_rmvdup);