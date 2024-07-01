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
#include <iostream>
#include <../bplus_tree.hpp>
//#include <../bplus_tree_light.hpp>
#include <omp.h>
#include <queue>
#include <mutex>
#include <condition_variable>

#include <string>


bool debug = false;

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
         << "    -i, --index FILE              input index file" << endl
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
        std::cerr << "The maximum kmer size is " << KeyType::KMER_MAX_LENGTH << std::endl;
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
//        cout << kmers.size() << endl;
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
//            cout << " The pos is: " << gbwtgraph::id(pos) << " " << gbwtgraph::offset(pos) << " " << gbwtgraph::is_rev(pos) << endl;

            if (debug){
                if (gbwtgraph::Key64::encode("GACAAATCTGGGTTCAAATCCTCACTTTG") == kmer.key){
                    cout << "The key is: " << kmer.key << " " << kmer.offset << " id " << graph.get_id(*iter) << " rev " << graph.get_is_reverse(*iter) << " kmer offset " << kmer.offset << " node_start " << node_start << " node_length " << node_length  << endl;
                }

            }

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
//    gbwtgraph::for_each_nonredundant_window(graph, k, hash_kmers, true);
    gbwtgraph::for_each_haplotype_window(graph, k, hash_kmers, true);
    for (int thread_id = 0; thread_id < threads; thread_id++) { flush_cache(thread_id); }

}


// This function input is the OCC vector of the end of the sequences (#) and it returns the sorted end_of_seq vector
// which is the sorted vector of pairs (i, SA[i]) for the end of each sequence which is (ISA[j], j)
vector<pair<uint64_t, uint64_t>> sort_end_of_seq(vector< pair<uint64_t, uint64_t> > &OCC) {


    // Sort the end_of_seq vector by the second element of each pair
    sort(OCC.begin(), OCC.end(),
         [](const pair<uint64_t, uint64_t> &a, const pair<uint64_t, uint64_t> &b) {
             return a.second < b.second;
         });

    return OCC;
}

// Create a Run struct to store the run data structure for the BPlusTree
struct Run {
    size_t start_position;
    gbwtgraph::Position graph_position;
//    char run_char; // the character of the run TODO: can handle this with 3 bits instead of 8 bits

    // Operators for the struct
    bool operator<(const Run &other) const {
        return start_position < other.start_position;
    }

//    bool operator==(const Run &other) const {
//        return start_position == other.start_position;
//    }

    bool operator==(const Run &other) const {
        return (start_position == other.start_position && graph_position == other.graph_position);
    }



    Run &operator=(const Run &other) {
        if (this == &other) {
            return *this;
        }
        start_position = other.start_position;
        graph_position = other.graph_position;
//        run_char = other.run_char;
        return *this;
    }

    // Set the Run struct to a zero one when it is assigned to 0
    Run &operator=(int zero) {
        if (zero == 0) {  // Ensures it only responds to 0
            start_position = 0;
            graph_position = gbwtgraph::Position::no_value();
        }
        return *this;
    }

    // print the Run struct
    friend std::ostream &operator<<(std::ostream &os, const Run &run) {
        os << "start_position: " << run.start_position
           << ", graph_position: " << run.graph_position.value;

        return os;
    }

};


// this function iterate over all of the kmers in the r-index and add the runs to the BPlusTree recursively
// TODO: make this function parallel
void kmers_to_bplustree(r_index<> &idx, BplusTree<Run> &bptree,
                        vg::hash_map<gbwtgraph::Key64::value_type, gbwtgraph::Position> &index, size_t k,
                        range_t interval, const string current_kmer) {
    if (current_kmer.length() == k && interval.first <= interval.second) {
        if (debug) {
            cout << "The current kmer is: " << current_kmer << endl;
            cout << "The interval is: " << interval.first << " " << interval.second << endl;
        }
        // creating the kmer with the key type
        gbwtgraph::Key64 kmer_key = gbwtgraph::Key64::encode(current_kmer);
        // check if the kmer_key is in the index and if it is add the run to the BPlusTree
        auto it = index.find(kmer_key.get_key());
        if (it != index.end()) {
            Run run = {interval.first, it->second};
            if (debug)
                cout << "The adding run is: " << run << " with len " << interval.second - interval.first + 1 << endl;
            bptree.insert(run, interval.second - interval.first + 1); // CHECK
        }
        return;
    }

    for (char base: {'A', 'C', 'G', 'T'}) {

        if (interval.first <= interval.second) {
            kmers_to_bplustree(idx, bptree, index, k, idx.LF(interval, base), base + current_kmer);
        }
    }
}


// Thread-safe queue implementation
template<typename T>
class ThreadSafeQueue {
public:
    void push(const T &item) {
        std::unique_lock<std::mutex> lock(mutex_);
        queue_.push(item);
        lock.unlock();
        cond_var_.notify_one();
    }

    bool try_pop(T &item) {
        std::unique_lock<std::mutex> lock(mutex_);
        if (queue_.empty()) {
            return false;
        }
        item = queue_.front();
        queue_.pop();
        return true;
    }

private:
    std::queue<T> queue_;
    std::mutex mutex_;
    std::condition_variable cond_var_;
};

void kmers_to_bplustree_worker(r_index<> &idx, ThreadSafeQueue<std::pair<Run, size_t>> &queue,
                               vg::hash_map<gbwtgraph::Key64::value_type, gbwtgraph::Position> &index, size_t k,
                               range_t interval, const string &current_kmer) {
    if (current_kmer.length() == k && interval.first <= interval.second) {
        if (debug) {
            cout << "The current kmer is: " << current_kmer << endl;
            cout << "The interval is: " << interval.first << " " << interval.second << endl;
        }
        // creating the kmer with the key type
        gbwtgraph::Key64 kmer_key = gbwtgraph::Key64::encode(current_kmer);
        // check if the kmer_key is in the index and if it is add the run to the queue
        auto it = index.find(kmer_key.get_key());
        if (it != index.end()) {
            Run run = {interval.first, it->second};
            if (debug)
                cout << "The adding run is: " << run << " with len " << interval.second - interval.first + 1 << endl;
            queue.push({run, interval.second - interval.first + 1});
        }
        return;
    }

    for (char base: {'A', 'C', 'G', 'T'}) {
        if (interval.first <= interval.second) {
            kmers_to_bplustree_worker(idx, queue, index, k, idx.LF(interval, base), base + current_kmer);
        }
    }
}

void parallel_kmers_to_bplustree(r_index<> &idx, BplusTree<Run> &bptree,
                                 vg::hash_map<gbwtgraph::Key64::value_type, gbwtgraph::Position> &index, size_t k,
                                 range_t interval) {
    // Thread-safe queue to collect results
    ThreadSafeQueue<std::pair<Run, size_t>> queue;

    // number of threads
    int threads = omp_get_max_threads();
    // splitting the starting range (0, idx.bwt_size() - 1) into parts and call the kmers_to_bplustree_worker function
    // for each part
    size_t part_size = (idx.bwt_size() - 1) / threads;
#pragma omp parallel for
    for (int i = 0; i < threads; i++) {
        size_t start = i * part_size;
        size_t end = (i + 1) * part_size - 1;
        if (i == threads - 1) {
            end = idx.bwt_size() - 1;
        }
        kmers_to_bplustree_worker(idx, queue, index, k, {start, end}, "");
    }

    // Single-threaded insertion into BPlusTree
    std::pair<Run, size_t> result;
    while (queue.try_pop(result)) {
        bptree.insert(result.first, result.second);
    }
}


Run generate_random_run() {
    static std::mt19937_64 rng(std::random_device{}());
    static std::uniform_int_distribution<size_t> pos_dist(1, 100000); // Adjust range as needed
    static std::uniform_int_distribution<uint64_t> graph_pos_dist(1, 15); // Adjust range as needed

    Run run;
    run.start_position = pos_dist(rng);
    vg::pos_t pos{graph_pos_dist(rng), false, 10};
    run.graph_position = gbwtgraph::Position::encode(pos);
    return run;
}

// Function to perform a randomized test on BPlusTree
void randomized_test_bplustree(int num_tests) {
    BplusTree<Run> bptree(15); // Adjust the degree of BPlusTree as needed
    BplusTree<Run> bptree2(16);
    BplusTree<Run> bptree3(num_tests + 15);


    // Insert random Run objects into the BPlusTree
    for (int i = 0; i < num_tests; ++i) {
        Run run = generate_random_run();
        size_t run_length = 10; // Fixed run length; adjust if necessary
        std::cout << "Inserting Run: " << run << " with length " << run_length << std::endl;
        bptree.insert(run, run_length);
        bptree2.insert(run, run_length);
        bptree3.insert(run, run_length);
    }

    // Check the integrity of the BPlusTree
//    bptree.bpt_check_items();

    // Additional checks to ensure the tree is structured correctly
    // For example: ensure items are in sorted order, no duplicates, etc.
    auto it = bptree.begin();
    auto prev_it = it;
    if (it != bptree.end()) ++it;
    while (it != bptree.end()) {
//        cout << "Checking " << *prev_it << " and " << *it << endl;
        assert((*prev_it).start_position <= (*it).start_position);
        assert((*prev_it).graph_position.value != (*it).graph_position.value);
        ++prev_it;
        ++it;
    }

    cout << "elements are sorted correctly" << endl;

    auto it1 = bptree.begin();
    auto it2 = bptree2.begin();
    auto it3 = bptree3.begin();
    while (it1 != bptree.end() && it2 != bptree2.end() && it3 != bptree3.end()) {
        cout << "Checking " << *it1 << " and " << *it2 << " and " << *it3 << endl;
        assert(*it1 == *it2);
        assert(*it1 == *it3);
        ++it1;
        ++it2;
        ++it3;
    }
    assert(it1 == bptree.end()); // Ensure both iterators reached the end of bptree
    assert(it2 == bptree.end());
    assert(it3 == bptree.end());

    std::cout << "Randomized test completed successfully!" << std::endl;
}




// This function iterate the bptree and try to extend the kmers backward (from the start position) on the graph
vector< pair<Run, size_t> > extend_kmers_on_graph(GBWTGraph &graph, r_index<> &idx, BplusTree<Run> &bptree) {

    vector< pair<Run, size_t> > extension_candidates;

    // first want to iterate over the BPlusTree and get the graph position for each interval
    for (auto it = bptree.begin(); it != bptree.end(); ++it) {
        Run current_item = *it;
        if (current_item.graph_position.value != 0) { // Check if the current item is not a gap
            auto next_it = it;
            ++next_it; // Move to the next element
            if (next_it != bptree.end()) { // Check if the next element is not the end
                Run next_item = *next_it; // Get the next item

                auto current_starting_pos = current_item.start_position;
                auto next_starting_pos = next_item.start_position;
                auto current_graph_pos = current_item.graph_position;

                // decode the Position to a pos_t
                pos_t current_pos = current_graph_pos.decode();

                // get the traversal of the graph position
                handle_t current_handle = graph.get_handle(vg::id(current_pos), false);


                // the bwt interval of the current kmer
                range_t current_interval = {current_starting_pos, next_starting_pos - 1};



                // trying to see if this kmer is extendable
                if (vg::is_rev(current_pos)) {
                    if (debug) cout << "The current position is reversed" << endl;
                } else {
                    // trying to see if the previous node base is unique in the graph
                    if (vg::offset(current_pos) > 0) {
                        // meaning that there is room for extension
                        // get the previous base
                        auto prev_base = graph.get_base(graph.get_handle(vg::id(current_pos), vg::is_rev(current_pos)),
                                                        vg::offset(current_pos) - 1);
                        auto prev_graph_pos = vg::pos_t{vg::id(current_pos), vg::is_rev(current_pos),
                                                        vg::offset(current_pos) - 1};

                        auto new_range = idx.LF(current_interval, prev_base);

                        if (new_range.first <= new_range.second){
                            extension_candidates.push_back({{new_range.first, gbwtgraph::Position::encode(prev_graph_pos)}, new_range.second - new_range.first + 1});

                        }


                    } else {
                        int prev_bases_num = 0;
                        handle_t prev_node;
                        char prev_base;
                        pos_t prev_graph_pos;

                        graph.follow_edges(current_handle, true, [&](const handle_t& prev) {
                            prev_bases_num++;
                            if (prev_bases_num != 1) return false;
                            prev_base = graph.get_base(prev, graph.get_length(prev) - 1);
                            prev_node = prev;
                            return true;
                        });

                        if (prev_bases_num == 1){
                            // add to the candidates
                            prev_graph_pos = vg::pos_t{graph.get_id(prev_node), graph.get_is_reverse(prev_node),
                                                       graph.get_length(prev_node) - 1};
                            auto new_range = idx.LF(current_interval, prev_base);
                            if (new_range.first <= new_range.second){
                                extension_candidates.push_back({{new_range.first, gbwtgraph::Position::encode(prev_graph_pos)}, new_range.second - new_range.first + 1});
                            }
                        }
                    }

                }
            }
        }
    }

    return extension_candidates;
}

vector< pair<Run, size_t> > extend_kmers_bfs(GBWTGraph &graph, r_index<> &idx, BplusTree<Run> &bptree) {
    bool succ = true;
    vector< pair<Run, size_t> > extension_candidates;

    // Queue for BFS
    std::queue< pair<Run, size_t> > bfs_queue;

    // Add initial runs to the queue
    for (auto it = bptree.begin(); it != bptree.end(); ++it) {
        Run current_item = *it;
        if (current_item.graph_position.value != 0) {
            auto next_it = it;
            ++next_it;
            if (next_it != bptree.end()) {
                Run next_item = *next_it;
                bfs_queue.push(make_pair(current_item, next_item.start_position));

//                auto current_starting_pos = current_item.start_position;
//                auto next_starting_pos = next_item.start_position;
//                auto current_graph_pos = current_item.graph_position;

            }

        }

    }

    // BFS to extend kmers
    while (!bfs_queue.empty()) {
        pair <Run, size_t> current_pair = bfs_queue.front();
        Run current_item = current_pair.first;
        size_t interval_end = current_pair.second;
        bfs_queue.pop();

        if (current_item.graph_position.value != 0) { // Check if the current item is not a gap
            auto current_starting_pos = current_item.start_position;
            auto current_graph_pos = current_item.graph_position;

            // Decode the Position to a pos_t
            pos_t current_pos = current_graph_pos.decode();

            // Get the traversal of the graph position
            handle_t current_handle = graph.get_handle(vg::id(current_pos), false);

            // The BWT interval of the current kmer
            range_t current_interval = {current_starting_pos, interval_end - 1}; // Start with single interval

            if (vg::is_rev(current_pos)) {
                if (debug) cout << "The current position is reversed" << endl;
            } else {
                // trying to see if the previous node base is unique in the graph
                if (vg::offset(current_pos) > 0) {
                    // meaning that there is room for extension
                    // get the previous base
                    auto prev_base = graph.get_base(graph.get_handle(vg::id(current_pos), vg::is_rev(current_pos)),
                                                    vg::offset(current_pos) - 1);
                    auto prev_graph_pos = vg::pos_t{vg::id(current_pos), vg::is_rev(current_pos),
                                                    vg::offset(current_pos) - 1};

                    auto new_range = idx.LF(current_interval, prev_base);

                    if (new_range.first <= new_range.second){
                        Run temp_run = {new_range.first, gbwtgraph::Position::encode(prev_graph_pos)};
//                        extension_candidates.push_back({temp_run, new_range.second - new_range.first + 1});
//                        bptree.insert(temp_run, new_range.second - new_range.first + 1);
//                        bfs_queue.push(make_pair(temp_run, new_range.second + 1));
                        if (bptree.insert_success(temp_run, new_range.second - new_range.first + 1)){
                            bfs_queue.push(make_pair(temp_run, new_range.second + 1));
                        }
//                        bfs_queue.push(temp_run, new_range.second + 1);



                    }


                } else {
                    int prev_bases_num = 0;
                    handle_t prev_node;
                    char prev_base;
                    pos_t prev_graph_pos;

                    graph.follow_edges(current_handle, true, [&](const handle_t& prev) {
                        prev_bases_num++;
                        if (prev_bases_num != 1) return false;
                        prev_base = graph.get_base(prev, graph.get_length(prev) - 1);
                        prev_node = prev;
                        return true;
                    });

                    if (prev_bases_num == 1){
                        // add to the candidates
                        prev_graph_pos = vg::pos_t{graph.get_id(prev_node), graph.get_is_reverse(prev_node),
                                                   graph.get_length(prev_node) - 1};
                        auto new_range = idx.LF(current_interval, prev_base);
                        if (new_range.first <= new_range.second){
                            Run new_run = {new_range.first, gbwtgraph::Position::encode(prev_graph_pos)};
//                            extension_candidates.push_back({new_run, new_range.second - new_range.first + 1});

                            if (bptree.insert_success(new_run, new_range.second - new_range.first + 1)){
                                bfs_queue.push(make_pair(new_run, new_range.second + 1));
                            }
//                            bptree.insert(new_run, new_range.second - new_range.first + 1);
//                            bfs_queue.push(make_pair(new_run, new_range.second + 1));


                        }
                    }
                }

            }
        }
    }

    return extension_candidates;
}




vector< pair<Run, size_t> > extend_kmers_bfs_parallel(GBWTGraph &graph, r_index<> &idx, BplusTree<Run> &bptree, int batch_size) {
    vector< pair<Run, size_t> > extension_candidates;

    int num_threads = omp_get_max_threads();
    vector<std::queue< pair<Run, size_t> >> bfs_queues(num_threads);
    vector<vector<pair<Run, size_t>>> batches(num_threads);

    // Mutex for synchronizing access to bptree
    std::mutex bptree_mutex;

    // Add initial runs to the queues
    int item_count = 0;
    for (auto it = bptree.begin(); it != bptree.end(); ++it) {
        Run current_item = *it;
        if (current_item.graph_position.value != 0) {
            auto next_it = it;
            ++next_it;
            if (next_it != bptree.end()) {
                Run next_item = *next_it;
                bfs_queues[item_count % num_threads].push(make_pair(current_item, next_item.start_position));
                ++item_count;
            }
        }
    }

    // Parallel BFS to extend kmers
#pragma omp parallel
    {
        int thread_num = omp_get_thread_num();
        while (!bfs_queues[thread_num].empty()) {
            pair<Run, size_t> current_pair = bfs_queues[thread_num].front();
            bfs_queues[thread_num].pop();

            Run current_item = current_pair.first;
            size_t interval_end = current_pair.second;

            if (current_item.graph_position.value != 0) {
                auto current_starting_pos = current_item.start_position;
                auto current_graph_pos = current_item.graph_position;

                // Decode the Position to a pos_t
                pos_t current_pos = current_graph_pos.decode();

                // Get the traversal of the graph position
                handle_t current_handle = graph.get_handle(vg::id(current_pos), false);

                // The BWT interval of the current kmer
                range_t current_interval = {current_starting_pos, interval_end - 1};

                if (!vg::is_rev(current_pos)) {
                    if (vg::offset(current_pos) > 0) {
                        auto prev_base = graph.get_base(graph.get_handle(vg::id(current_pos), vg::is_rev(current_pos)),
                                                        vg::offset(current_pos) - 1);
                        auto prev_graph_pos = vg::pos_t{vg::id(current_pos), vg::is_rev(current_pos),
                                                        vg::offset(current_pos) - 1};

                        auto new_range = idx.LF(current_interval, prev_base);

                        if (new_range.first <= new_range.second) {
                            Run temp_run = {new_range.first, gbwtgraph::Position::encode(prev_graph_pos)};
                            batches[thread_num].push_back(make_pair(temp_run, new_range.second - new_range.first + 1));

                            if (batches[thread_num].size() >= batch_size) {
                                std::lock_guard<std::mutex> lock(bptree_mutex);
                                for (const auto& item : batches[thread_num]) {
//                                    bptree.insert(item.first, item.second);
//                                    bfs_queues[thread_num].push(make_pair(item.first, item.second + item.first.start_position));
                                    if (bptree.insert_success(item.first, item.second)) {
                                        bfs_queues[thread_num].push(make_pair(item.first, item.second + item.first.start_position));
                                    }
                                }
                                batches[thread_num].clear();
                            }
                        }
                    } else {
                        int prev_bases_num = 0;
                        handle_t prev_node;
                        char prev_base;
                        pos_t prev_graph_pos;

                        graph.follow_edges(current_handle, true, [&](const handle_t& prev) {
                            prev_bases_num++;
                            if (prev_bases_num != 1) return false;
                            prev_base = graph.get_base(prev, graph.get_length(prev) - 1);
                            prev_node = prev;
                            return true;
                        });

                        if (prev_bases_num == 1) {
                            prev_graph_pos = vg::pos_t{graph.get_id(prev_node), graph.get_is_reverse(prev_node),
                                                       graph.get_length(prev_node) - 1};
                            auto new_range = idx.LF(current_interval, prev_base);
                            if (new_range.first <= new_range.second) {
                                Run new_run = {new_range.first, gbwtgraph::Position::encode(prev_graph_pos)};
                                batches[thread_num].push_back({new_run, new_range.second - new_range.first + 1});

                                if (batches[thread_num].size() >= batch_size) {
                                    std::lock_guard<std::mutex> lock(bptree_mutex);
                                    for (const auto& item : batches[thread_num]) {
//                                        bptree.insert(item.first, item.second);
//                                        bfs_queues[thread_num].push(make_pair(item.first, item.second + item.first.start_position));
                                        if (bptree.insert_success(item.first, item.second)) {
                                            bfs_queues[thread_num].push(make_pair(item.first, item.second + item.first.start_position));
                                        }
                                    }
                                    batches[thread_num].clear();
                                }
                            }
                        }
                    }
                }
            }
        }

        // Process remaining items in the local batch
        if (!batches[thread_num].empty()) {
            std::lock_guard<std::mutex> lock(bptree_mutex);
            for (const auto& item : batches[thread_num]) {
//                bptree.insert(item.first, item.second);
//                bfs_queues[thread_num].push(make_pair(item.first, item.second + item.first.start_position));
                if (bptree.insert_success(item.first, item.second)) {
                    bfs_queues[thread_num].push(make_pair(item.first, item.second + item.first.start_position));
                }
            }
            batches[thread_num].clear();
        }
    }

    return extension_candidates;
}












int main_panindexer(int argc, char **argv) {

    if (argc <= 3) {
        help_panindexer(argv);
        return 1;
    }

    int threads = 1;
    size_t k = 29; // TODO: change the default value
    string graph_file = argv[2];
    string index_file;

    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
                {
                        {"graph",       required_argument, 0, 'g'},
                        {"kmer-length", required_argument, 0, 'k'},
                        {"index",       required_argument, 0, 'i'},
                        {"threads",     required_argument, 0, 't'},
                        {0, 0,                             0, 0}
                };
        int option_index = 0;
        int c = getopt_long(argc, argv, "g:k:i:t:", long_options, &option_index);

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
            case 'i':
                index_file = optarg;
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
//    r_index<> idx;
//    idx.load(in);


//
//    string pattern = "#";
//
////    locate<r_index<> >(in, pattern);

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

//    int num_tests = 1000000; // Number of random tests
//    randomized_test_bplustree(num_tests);
//    return 0;



    unique_ptr<gbwtgraph::GBZ> gbz;
    cout << "Loading the graph file" << endl;
    auto input = vg::io::VPKG::try_load_first<gbwtgraph::GBZ, gbwtgraph::GBWTGraph, HandleGraph>(graph_file);
    gbz = std::move(get<0>(input));

//    auto temp = gbz->index.extract(0);
//    for (auto& i : temp) {
//        auto x = GBWTGraph::node_to_handle(i);
//        cout << gbz->graph.get_sequence(x) << endl;
//
//    }
//
//    cout << "second           ASFSŠ" << endl;
//
//    auto temp1 = gbz->index.extract(2);
//    for (auto& i : temp1) {
//        auto x = GBWTGraph::node_to_handle(i);
//        cout << gbz->graph.get_sequence(x) << endl;
//
//    }


//    cout << "Reading the rindex file" << endl;
//    std::ifstream in(index_file);
//    bool fast;
//    //fast or small index?
//    in.read((char *) &fast, sizeof(fast));
//    r_index<> idx;
//    idx.load(in);





    bool progress = true;
    auto threshold = 0;
    auto space_efficient_counting = false;

    typedef gbwtgraph::Key64::value_type kmer_type;

    hash_map<kmer_type, gbwtgraph::Position> index;
    cout << "Computing the unique kmers in the graph" << endl;

    unique_kmers_parallel<gbwtgraph::Key64>(gbz->graph, index, k);


    BplusTree<Run> bptree(15); // TODO: determine the BPlusTree degree

    cout << "Reading the rindex file" << endl;
    std::ifstream in(index_file);
    bool fast;
    //fast or small index?
    in.read((char *) &fast, sizeof(fast));
    r_index<> idx;
    idx.load(in);





    cout << "Adding the kmers to the BPlusTree" << endl;
//    kmers_to_bplustree(idx, bptree, index, k, {0, idx.bwt_size() - 1}, "");
    parallel_kmers_to_bplustree(idx, bptree, index, k, {0, idx.bwt_size() - 1});


    // computing some statistics
    auto unique_kmers_size = index.size();
    auto bptree_items = bptree.get_bpt_size();
    auto bwt_size = idx.bwt_size();
    size_t tag_arrays_covered = 0;

    cout << "The number of unique kmers in the index is: " << unique_kmers_size << endl;
    cout << "The number of items in the BPlusTree is: " << bptree_items << endl;
    cout << "The size of the BWT is: " << bwt_size << endl;


    cout << "calculating the fraction of the tag arrays covered" << endl;
    // calculating the fraction of the tag arrays covered
    for (auto it = bptree.begin(); it != bptree.end(); ++it) {
        Run current_item = *it;
        if (current_item.graph_position.value != 0) { // Check if the current item is not a gap
            auto next_it = it;
            ++next_it; // Move to the next element
            if (next_it != bptree.end()) { // Check if the next element is not the end
                Run next_item = *next_it; // Get the next item
                tag_arrays_covered += (next_item.start_position - current_item.start_position);
            }
        }
    }



    cout << "The fraction of the tag arrays covered is: " << tag_arrays_covered << " / " << bwt_size << " = "
         << (double) tag_arrays_covered / bwt_size << endl;


    // calling the extension function
    for (int i = 0; i < 1; i++){
        cout << "Extending the kmers on the graph" << endl;
        auto extension_candidates = extend_kmers_bfs_parallel(gbz->graph, idx, bptree, 1024);
//        auto extension_candidates = extend_kmers_bfs(gbz->graph, idx, bptree);
        cout << "The extension candidates are: " << endl;
        cout << extension_candidates.size() << endl;

        // adding the extension candidates to the BPlusTree
//        for (auto &candidate: extension_candidates) {
//            bptree.insert(candidate.first, candidate.second);
//        }
    }




    cout << "The number of items in the BPlusTree is: " << bptree.get_bpt_size() << endl;
    tag_arrays_covered = 0;

    cout << "calculating the fraction of the tag arrays covered after one extension" << endl;
    // calculating the fraction of the tag arrays covered
    for (auto it = bptree.begin(); it != bptree.end(); ++it) {
        Run current_item = *it;
        if (current_item.graph_position.value != 0) { // Check if the current item is not a gap
            auto next_it = it;
            ++next_it; // Move to the next element
            if (next_it != bptree.end()) { // Check if the next element is not the end
                Run next_item = *next_it; // Get the next item
                tag_arrays_covered += (next_item.start_position - current_item.start_position);
//                cout << "The current item is: " << current_item << " NEXT " << next_item << endl;
//                if (next_item.start_position - current_item.start_position > 500){
//                    cout << "The current item is: " << current_item << " The next item is: " << next_item << endl;
//                }
            }
        }
    }

    cout << "The fraction of the tag arrays covered is: " << tag_arrays_covered << " / " << bwt_size << " = "
         << (double) tag_arrays_covered / bwt_size << endl;



    string pattern = "$";


    // get the ISA values for the end of each sequence by searching for pattern # and locating the results
    auto OCC = idx.ISA(pattern);

    // print all the ISA values
    for (auto& i : OCC) {
        cout << i.first << " " << i.second << endl;
    }



    // print the bwt of the index
//    cout << idx.get_bwt() << endl;

    // using the sort_end_of_seq function
    // the first item in the nth element of this is the end of the nth sequence in the bwt
    auto end_of_seq = sort_end_of_seq(OCC);

    for (auto& i : end_of_seq) {
        cout << i.first << " " << i.second << endl;
    }

    // traversing the sequenece in the RLBWT and the graph
    int seq_num = 0;
    auto seq_graph_nodes = gbz->index.extract(seq_num);
    auto bwt_index = end_of_seq[seq_num].first;



//    // traversing the oriented path of a sequence backward on the graph
//    for (int i = seq_graph_nodes.size() - 1; i >= 0; --i){
//
//        //traverse the graph backwards
//        auto node = GBWTGraph::node_to_handle(seq_graph_nodes[i]);
//
//        // search the
//
////        cout << gbz->graph.get_length(node) << endl;
//
//
//        for (int j = gbz->graph.get_length(node) - 1; j >= 0; --j){
//            bwt_index = idx.LF(bwt_index);
//            cout << idx.F_at(bwt_index) << " " << gbz->graph.get_base(node, j) << endl;
//            cout << bwt_index << endl;
//
//        }
//
//
////        cout << gbz->graph.get_sequence(node) << endl;
////        cout << gbz->graph.get_length(node) << endl;
//
//    }

    auto current_nodes_index = seq_graph_nodes.size() - 1;
    auto current_node = GBWTGraph::node_to_handle(seq_graph_nodes[current_nodes_index]);
    auto in_node_index = gbz->graph.get_length(current_node) - 1;

    // traversing the RLBWT of a sequence
    while (true){
        // moving backwards
        bwt_index = idx.LF(bwt_index);
        auto first = idx.F_at(bwt_index);
        if (first == '$'){
            break;
        }

        // make sure there is still room to travese on the current node
        if (in_node_index == -1){
            if (current_nodes_index == 0){
                break;
            }
            current_nodes_index--;
            current_node = GBWTGraph::node_to_handle(seq_graph_nodes[current_nodes_index]);
            in_node_index = gbz->graph.get_length(current_node) - 1;
        }


        assert(first == gbz->graph.get_base(current_node, in_node_index));
        //traverse the nodes on the graph to get the same base


//        cout << "F " << first << " graph base " << gbz->graph.get_base(current_node, in_node_index) << endl;
//        cout << bwt_index << endl;
        in_node_index--;


        // have to check if the current bwt position is in the bptree or not
        // calling search function if it return a gap run then this position is not currently in the tree and have to add it
        // if it returns a non-gap run then this position is already in the tree and we can continue traversing the bwt


        auto bptree_search = bptree.search(bwt_index);
//        cout << "tree search " << bptree_search << endl;

        // the current position is not in the tree
        if (bptree_search.graph_position.value == 0){

            // adding the current position to the tree


            pos_t current_pos = vg::pos_t{gbz->graph.get_id(current_node), gbz->graph.get_is_reverse(current_node),
                                          in_node_index};

            Run current_run = {bwt_index, gbwtgraph::Position::encode(current_pos)};
            bptree.insert(current_run, 1);




        } else {
            // not adding the current position to the tree however checking if the tree position and the current graph
            // positions are the same
            pos_t current_pos = vg::pos_t{gbz->graph.get_id(current_node), gbz->graph.get_is_reverse(current_node),
                                          in_node_index + 1};


            if (gbwtgraph::Position::encode(current_pos).value != bptree_search.graph_position.value) {
                cout << "reverse? " << gbz->graph.get_is_reverse(current_node) << endl;
                cout << "node len " << gbz->graph.get_length(current_node) << endl;
                cout << "The graph position in the tree " << vg::id(bptree_search.graph_position.decode()) << " "
                     << vg::offset(bptree_search.graph_position.decode()) << " the actual brute force graph position "
                     << vg::id(current_pos) << " " << vg::offset(current_pos) << endl;
            }

            assert(gbwtgraph::Position::encode(current_pos).value == bptree_search.graph_position.value);
        }








    }


    cout << "The number of items in the BPlusTree is: " << bptree.get_bpt_size() << endl;
    tag_arrays_covered = 0;

    cout << "calculating the fraction of the tag arrays covered after one extension" << endl;
    // calculating the fraction of the tag arrays covered
    for (auto it = bptree.begin(); it != bptree.end(); ++it) {
        Run current_item = *it;
        if (current_item.graph_position.value != 0) { // Check if the current item is not a gap
            auto next_it = it;
            ++next_it; // Move to the next element
            if (next_it != bptree.end()) { // Check if the next element is not the end
                Run next_item = *next_it; // Get the next item
                tag_arrays_covered += (next_item.start_position - current_item.start_position);
//                cout << "The current item is: " << current_item << " NEXT " << next_item << endl;
//                if (next_item.start_position - current_item.start_position > 500){
//                    cout << "The current item is: " << current_item << " The next item is: " << next_item << endl;
//                }
            }
        }
    }

    cout << "The fraction of the tag arrays covered is: " << tag_arrays_covered << " / " << bwt_size << " = "
         << (double) tag_arrays_covered / bwt_size << endl;











//    BplusTree<Run> bptree(10);
//    bptree.insert({5,1}, 3);
//
//    bptree.insert({11,5}, 2);
//
//    bptree.insert({10,5}, 7);
//    bptree.print_whole();
//
//    cout << "-------------------------------1" << endl;
//
//    BplusTree<Run> bptree1(10);
//    bptree1.insert({5,1}, 3);
//
//    bptree1.insert({11,5}, 2);
//
//    bptree1.insert({4,1}, 6);
//    bptree1.print_whole();
//
//    cout << "-------------------------------2" << endl;
//
//    BplusTree<Run> bptree3(10);
//    bptree3.insert({5,1}, 3);
//
//    bptree3.insert({11,5}, 2);
//
//    bptree3.insert({5,1}, 4);
//    bptree3.print_whole();
//
//
//    cout << "-------------------------------3" << endl;
//
//    BplusTree<Run> bptree2(10);
//    bptree2.insert({5,1}, 3);
//
//    bptree2.insert({11,5}, 2);
//
//    bptree2.insert({4,1}, 7);
//    bptree2.print_whole();
//
//    cout << "-------------------------------4" << endl;
//
//    BplusTree<Run> bptree4(10);
//    bptree4.insert({5,1}, 3);
//
//    bptree4.insert({11,5}, 2);
//
//    bptree4.insert({6,1}, 1);
//    bptree4.print_whole();






//
//        for (size_t i = 3; i < 10; i+=2){
//        bptree.insert({i, 3}, 1);
//        cout << "The B+ tree is: " << endl;
//        bptree.print_whole();
//    }
//    for (size_t i = 4; i < 9; i+=4){
//        bptree.insert({i, 3}, 1);
//        cout << "The B+ tree is: " << endl;
//        bptree.print_whole();
//    }
//
//    for (size_t i = 30; i > 9; i-=2){
//        bptree.insert({i, 3}, 1);
//        cout << "The B+ tree is: " << endl;
//        bptree.print_whole();
//    }
//    for (size_t i = 13; i < 30; i+=2){
//        bptree.insert({i, 3}, 1);
//        cout << "The B+ tree is: " << endl;
//        bptree.print_whole();
//    }


//    for (size_t i = 10; i < 30; i+=2){
//        bptree.insert({i, 3}, 1);
//        cout << "The B+ tree is: " << endl;
//        bptree.print_whole();
//    }
//    for (size_t i = 29; i > 12; i-=2){
//        bptree.insert({i, 3}, 1);
//        cout << "The B+ tree is: " << endl;
//        bptree.print_whole();
//    }
//

////    bptree.insert({2,4}, 2);
//    bptree.insert({7,8}, 1);
//    bptree.insert({5,5}, 1);
//    bptree.insert({9,5}, 1);
//    bptree.insert({2,7}, 1);
//    bptree.print_whole();


//    cout << node.is_leaf() << endl;
//    node.insert({5,2}, 2);
//    node.print();
//    node.insert({9,2}, 5);
//    node.print();
//    node.insert({7,2}, 2);
//
//
//    node.print();


    return 0;
}


static vg::subcommand::Subcommand vg_panindexer("panindexer", "index pangenomes", PIPELINE, main_panindexer);
