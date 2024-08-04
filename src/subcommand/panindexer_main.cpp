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
#include <tag_arrays.hpp>

#include <string>

#define TIME 1
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
         << "    -t, --threads N          number of threads to use [1]" << endl
         << "    -b, --build-index             build the index" << endl;
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
//        std::sort(current_cache.begin(), current_cache.end(), [](const auto &a, const auto &b) {
//            return a.first < b.first;
//        });
//        // Remove duplicates
//        auto last = std::unique(current_cache.begin(), current_cache.end(), [](const auto &a, const auto &b) {
//            return a.first == b.first;
//        });
//        current_cache.erase(last, current_cache.end());
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
                    index.erase(it);
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
//        cerr << kmers.size() << endl;
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
//            cerr << " The pos is: " << gbwtgraph::id(pos) << " " << gbwtgraph::offset(pos) << " " << gbwtgraph::is_rev(pos) << endl;

            if (debug) {
                if (gbwtgraph::Key64::encode("GACAAATCTGGGTTCAAATCCTCACTTTG") == kmer.key) {
                    cerr << "The key is: " << kmer.key << " " << kmer.offset << " id " << graph.get_id(*iter) << " rev "
                         << graph.get_is_reverse(*iter) << " kmer offset " << kmer.offset << " node_start "
                         << node_start << " node_length " << node_length << endl;
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
vector<pair<uint64_t, uint64_t>> sort_end_of_seq(vector<pair<uint64_t, uint64_t> > &OCC) {


    // Sort the end_of_seq vector by the second element of each pair
    sort(OCC.begin(), OCC.end(),
         [](const pair<uint64_t, uint64_t> &a, const pair<uint64_t, uint64_t> &b) {
             return a.second < b.second;
         });

    return OCC;
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
            cerr << "The current kmer is: " << current_kmer << endl;
            cerr << "The interval is: " << interval.first << " " << interval.second << endl;
        }
        // creating the kmer with the key type
        gbwtgraph::Key64 kmer_key = gbwtgraph::Key64::encode(current_kmer);
        // check if the kmer_key is in the index and if it is add the run to the queue
        auto it = index.find(kmer_key.get_key());
        if (it != index.end()) {
            Run run = {interval.first, it->second};
            if (debug)
                cerr << "The adding run is: " << run << " with len " << interval.second - interval.first + 1 << endl;
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
        std::cerr << "Inserting Run: " << run << " with length " << run_length << std::endl;
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
//        cerr << "Checking " << *prev_it << " and " << *it << endl;
        assert((*prev_it).start_position <= (*it).start_position);
        assert((*prev_it).graph_position.value != (*it).graph_position.value);
        ++prev_it;
        ++it;
    }

    cerr << "elements are sorted correctly" << endl;

    auto it1 = bptree.begin();
    auto it2 = bptree2.begin();
    auto it3 = bptree3.begin();
    while (it1 != bptree.end() && it2 != bptree2.end() && it3 != bptree3.end()) {
        cerr << "Checking " << *it1 << " and " << *it2 << " and " << *it3 << endl;
        assert(*it1 == *it2);
        assert(*it1 == *it3);
        ++it1;
        ++it2;
        ++it3;
    }
    assert(it1 == bptree.end()); // Ensure both iterators reached the end of bptree
    assert(it2 == bptree.end());
    assert(it3 == bptree.end());

    std::cerr << "Randomized test completed successfully!" << std::endl;
}


// This function iterate the bptree and try to extend the kmers backward (from the start position) on the graph
vector<pair<Run, size_t> > extend_kmers_on_graph(GBWTGraph &graph, r_index<> &idx, BplusTree<Run> &bptree) {

    vector<pair<Run, size_t> > extension_candidates;

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
                    if (debug) cerr << "The current position is reversed" << endl;
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

                        if (new_range.first <= new_range.second) {
                            extension_candidates.push_back(
                                    {{new_range.first, gbwtgraph::Position::encode(prev_graph_pos)},
                                     new_range.second - new_range.first + 1});

                        }


                    } else {
                        int prev_bases_num = 0;
                        handle_t prev_node;
                        char prev_base;
                        pos_t prev_graph_pos;

                        graph.follow_edges(current_handle, true, [&](const handle_t &prev) {
                            prev_bases_num++;
                            if (prev_bases_num != 1) return false;
                            prev_base = graph.get_base(prev, graph.get_length(prev) - 1);
                            prev_node = prev;
                            return true;
                        });

                        if (prev_bases_num == 1) {
                            // add to the candidates
                            prev_graph_pos = vg::pos_t{graph.get_id(prev_node), graph.get_is_reverse(prev_node),
                                                       graph.get_length(prev_node) - 1};
                            auto new_range = idx.LF(current_interval, prev_base);
                            if (new_range.first <= new_range.second) {
                                extension_candidates.push_back(
                                        {{new_range.first, gbwtgraph::Position::encode(prev_graph_pos)},
                                         new_range.second - new_range.first + 1});
                            }
                        }
                    }

                }
            }
        }
    }

    return extension_candidates;
}

vector<pair<Run, size_t> > extend_kmers_bfs(GBWTGraph &graph, r_index<> &idx, BplusTree<Run> &bptree) {
    bool succ = true;
    vector<pair<Run, size_t> > extension_candidates;

    // Queue for BFS
    std::queue<pair<Run, size_t> > bfs_queue;

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
        pair<Run, size_t> current_pair = bfs_queue.front();
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
                if (debug) cerr << "The current position is reversed" << endl;
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

                    if (new_range.first <= new_range.second) {
                        Run temp_run = {new_range.first, gbwtgraph::Position::encode(prev_graph_pos)};
//                        extension_candidates.push_back({temp_run, new_range.second - new_range.first + 1});
//                        bptree.insert(temp_run, new_range.second - new_range.first + 1);
//                        bfs_queue.push(make_pair(temp_run, new_range.second + 1));
                        if (bptree.insert_success(temp_run, new_range.second - new_range.first + 1)) {
                            bfs_queue.push(make_pair(temp_run, new_range.second + 1));
                        }
//                        bfs_queue.push(temp_run, new_range.second + 1);



                    }


                } else {
                    int prev_bases_num = 0;
                    handle_t prev_node;
                    char prev_base;
                    pos_t prev_graph_pos;

                    graph.follow_edges(current_handle, true, [&](const handle_t &prev) {
                        prev_bases_num++;
                        if (prev_bases_num != 1) return false;
                        prev_base = graph.get_base(prev, graph.get_length(prev) - 1);
                        prev_node = prev;
                        return true;
                    });

                    if (prev_bases_num == 1) {
                        // add to the candidates
                        prev_graph_pos = vg::pos_t{graph.get_id(prev_node), graph.get_is_reverse(prev_node),
                                                   graph.get_length(prev_node) - 1};
                        auto new_range = idx.LF(current_interval, prev_base);
                        if (new_range.first <= new_range.second) {
                            Run new_run = {new_range.first, gbwtgraph::Position::encode(prev_graph_pos)};
//                            extension_candidates.push_back({new_run, new_range.second - new_range.first + 1});

                            if (bptree.insert_success(new_run, new_range.second - new_range.first + 1)) {
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


vector<pair<Run, size_t> >
extend_kmers_bfs_parallel(GBWTGraph &graph, r_index<> &idx, BplusTree<Run> &bptree, int batch_size) {
    vector<pair<Run, size_t> > extension_candidates;

    int num_threads = omp_get_max_threads();
    vector<std::queue<pair<Run, size_t> >> bfs_queues(num_threads);
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
                                for (const auto &item: batches[thread_num]) {
//                                    bptree.insert(item.first, item.second);
//                                    bfs_queues[thread_num].push(make_pair(item.first, item.second + item.first.start_position));
                                    if (bptree.insert_success(item.first, item.second)) {
                                        bfs_queues[thread_num].push(
                                                make_pair(item.first, item.second + item.first.start_position));
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

                        graph.follow_edges(current_handle, true, [&](const handle_t &prev) {
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
                                    for (const auto &item: batches[thread_num]) {
//                                        bptree.insert(item.first, item.second);
//                                        bfs_queues[thread_num].push(make_pair(item.first, item.second + item.first.start_position));
                                        if (bptree.insert_success(item.first, item.second)) {
                                            bfs_queues[thread_num].push(
                                                    make_pair(item.first, item.second + item.first.start_position));
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
            for (const auto &item: batches[thread_num]) {
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


void traverse_sequences(gbwtgraph::GBZ &gbz, BplusTree<Run> &bptree, r_index<> &idx,
                        vector<pair<uint64_t, uint64_t> > &end_of_seq) {
    auto number_of_sequences = end_of_seq.size();
    int traverse = 0;

    vector<int> tmp;

    vector<Run> tmp1;

    for (int seq_num = 0; seq_num < number_of_sequences; ++seq_num) {
        cerr << "running for sequence number " << seq_num << endl;
        cerr << traverse << endl;

        auto seq_graph_nodes = gbz.index.extract(seq_num * 2);
        auto bwt_index = end_of_seq[seq_num].first;


        auto current_nodes_index = seq_graph_nodes.size() - 1;
        auto current_node = GBWTGraph::node_to_handle(seq_graph_nodes[current_nodes_index]);
        auto in_node_index = gbz.graph.get_length(current_node) - 1;

        // traversing the RLBWT of a sequence
        while (true) {
            traverse++;

            // moving backwards
            bwt_index = idx.LF(bwt_index);
            auto first = idx.F_at(bwt_index);
            if (first == '$') {
                cerr << "The end of the sequence at bwt index " << bwt_index << endl;
                break;
            }


            // make sure there is still room to travese on the current node
            if (in_node_index == -1) {
                if (current_nodes_index == 0) {
                    break;
                }
                current_nodes_index--;
                current_node = GBWTGraph::node_to_handle(seq_graph_nodes[current_nodes_index]);
                in_node_index = gbz.graph.get_length(current_node) - 1;
            }


            //traverse the nodes on the graph to get the same base
            assert(first == gbz.graph.get_base(current_node, in_node_index));


            in_node_index--;


            // have to check if the current bwt position is in the bptree or not
            // calling search function if it return a gap run then this position is not currently in the tree and have to add it
            // if it returns a non-gap run then this position is already in the tree and we can continue traversing the bwt


            auto bptree_search = bptree.search(bwt_index);


            // the current position is not in the tree
            if (bptree_search.graph_position.value == 0) {


                // adding the current position to the tree
                pos_t current_pos = vg::pos_t{gbz.graph.get_id(current_node), gbz.graph.get_is_reverse(current_node),
                                              in_node_index + 1};

                Run current_run = {bwt_index, gbwtgraph::Position::encode(current_pos)};
//                bptree.insert(current_run, 1);
                tmp1.push_back(current_run);


            } else {
                // not adding the current position to the tree however checking if the tree position and the current graph
                // positions are the same
                pos_t current_pos = vg::pos_t{gbz.graph.get_id(current_node), gbz.graph.get_is_reverse(current_node),
                                              in_node_index + 1};


//                if (gbwtgraph::Position::encode(current_pos).value != bptree_search.graph_position.value) {
//                    cerr << "tree search " << bptree_search << endl;
//                    cerr << "brute " << gbwtgraph::Position::encode(current_pos).value << endl;
//                    cerr << "reverse? " << gbz->graph.get_is_reverse(current_node) << endl;
//                    cerr << "node len " << gbz->graph.get_length(current_node) << endl;
//                    cerr << "The graph position in the tree " << vg::id(bptree_search.graph_position.decode()) << " "
//                         << vg::offset(bptree_search.graph_position.decode()) << " the actual brute force graph position "
//                         << vg::id(current_pos) << " " << vg::offset(current_pos) << endl;
//                }

                assert(gbwtgraph::Position::encode(current_pos).value == bptree_search.graph_position.value);
            }


        }


    }

    cerr << "The traverse completed" << endl;

    // sorting the tmp1
    sort(tmp1.begin(), tmp1.end());

    // insert all the items in the tmp1 to the bptree
    for (auto &i: tmp1) {
        bptree.insert(i, 1);
    }
}

void traverse_sequences_parallel(gbwtgraph::GBZ &gbz, BplusTree<Run> &bptree, r_index<> &idx,
                                 vector<pair<uint64_t, uint64_t>> &end_of_seq) {
    auto number_of_sequences = end_of_seq.size();
    int traverse = 0;

    vector<int> tmp;

    vector<Run> tmp1;
    omp_lock_t lock;
    omp_init_lock(&lock);

    cerr << "Filling the gaps on the bptree" << endl;

#pragma omp parallel for
    for (int seq_num = 0; seq_num < number_of_sequences; ++seq_num) {
//        cerr << "running for sequence number " << seq_num << endl;
//        cerr << traverse << endl;

        auto seq_graph_nodes = gbz.index.extract(seq_num * 2);
        auto bwt_index = end_of_seq[seq_num].first;

        auto current_nodes_index = seq_graph_nodes.size() - 1;
        auto current_node = GBWTGraph::node_to_handle(seq_graph_nodes[current_nodes_index]);
        auto in_node_index = gbz.graph.get_length(current_node) - 1;

        vector<Run> local_tmp1;

        // traversing the RLBWT of a sequence
        while (true) {
#pragma omp atomic
            traverse++;

            // moving backwards
            bwt_index = idx.LF(bwt_index);
            auto first = idx.F_at(bwt_index);
            if (first == '$') {
//                cerr << "The end of the sequence at bwt index " << bwt_index << endl;
                break;
            }

            // make sure there is still room to traverse on the current node
            if (in_node_index == -1) {
                if (current_nodes_index == 0) {
                    break;
                }
                current_nodes_index--;
                current_node = GBWTGraph::node_to_handle(seq_graph_nodes[current_nodes_index]);
                in_node_index = gbz.graph.get_length(current_node) - 1;
            }

            // traverse the nodes on the graph to get the same base
            assert(first == gbz.graph.get_base(current_node, in_node_index));

            in_node_index--;

            // have to check if the current bwt position is in the bptree or not
            // calling search function if it return a gap run then this position is not currently in the tree and have to add it
            // if it returns a non-gap run then this position is already in the tree and we can continue traversing the bwt

            auto bptree_search = bptree.search(bwt_index);

            // the current position is not in the tree
            if (bptree_search.graph_position.value == 0) {
                // adding the current position to the tree
                pos_t current_pos = vg::pos_t{gbz.graph.get_id(current_node), gbz.graph.get_is_reverse(current_node),
                                              in_node_index + 1};

                Run current_run = {bwt_index, gbwtgraph::Position::encode(current_pos)};
                local_tmp1.push_back(current_run);

            } else {
                // not adding the current position to the tree however checking if the tree position and the current graph
                // positions are the same
                pos_t current_pos = vg::pos_t{gbz.graph.get_id(current_node), gbz.graph.get_is_reverse(current_node),
                                              in_node_index + 1};

                assert(gbwtgraph::Position::encode(current_pos).value == bptree_search.graph_position.value);
            }
        }

        // Lock to update shared data
        omp_set_lock(&lock);
        tmp1.insert(tmp1.end(), local_tmp1.begin(), local_tmp1.end());
        omp_unset_lock(&lock);
    }

    cerr << "The traverse completed" << endl;

    // Sort before inserting into the bptree
    sort(tmp1.begin(), tmp1.end());

    // insert all the items in the tmp1 to the bptree
    for (auto &i: tmp1) {
        bptree.insert(i, 1);
    }

    omp_destroy_lock(&lock);
}


template<typename T>
int count_bits_using_bitset(T value) {
    std::bitset<std::numeric_limits<T>::digits> bits(value);
    return bits.count();
}

void bin(unsigned int n) {
    /* step 1 */
    if (n > 1)
        bin(n / 2);

    /* step 2 */
    cerr << n % 2;
}


int main_panindexer(int argc, char **argv) {

    if (argc <= 3) {
        help_panindexer(argv);
        return 1;
    }

    int threads = 1;
    size_t k = 31; // TODO: change the default value
    string graph_file = argv[2];
    string index_file;
    bool build_index = false;

    optind = 2; // force optind past command positional argument
    while (true) {
        static struct option long_options[] =
                {
                        {"graph",       required_argument, 0, 'g'},
                        {"kmer-length", required_argument, 0, 'k'},
                        {"index",       required_argument, 0, 'i'},
                        {"threads",     required_argument, 0, 't'},
                        {"build-index", no_argument, 0, 'b'},
                        {0, 0,                             0, 0}
                };
        int option_index = 0;
        int c = getopt_long(argc, argv, "g:k:i:t:b", long_options, &option_index);

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
            case 'b':
                build_index = true;
                break;
            case 'h':
            case '?':
            default:
                help_panindexer(argv);
                exit(1);
        }

    }


//    int num_tests = 1000000; // Number of random tests
//    randomized_test_bplustree(num_tests);
//    return 0;

    cerr << "Reading the rindex file" << endl;
    std::ifstream in(index_file);
    bool fast;
    //fast or small index?
    in.read((char *) &fast, sizeof(fast));
    r_index<> idx;
    idx.load(in);


    if (build_index) {


        unique_ptr<gbwtgraph::GBZ> gbz;
        cerr << "Loading the graph file" << endl;
        auto input = vg::io::VPKG::try_load_first<gbwtgraph::GBZ, gbwtgraph::GBWTGraph, HandleGraph>(graph_file);
        gbz = std::move(get<0>(input));


        bool progress = true;
        auto threshold = 0;
        auto space_efficient_counting = false;

        typedef gbwtgraph::Key64::value_type kmer_type;

        hash_map<kmer_type, gbwtgraph::Position> index;
        cerr << "Computing the unique kmers in the graph" << endl;

#if TIME
        auto time1 = chrono::high_resolution_clock::now();
#endif

        unique_kmers_parallel<gbwtgraph::Key64>(gbz->graph, index, k);

#if TIME
        auto time2 = chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration1 = time2 - time1;
        std::cerr << "Indexing unique kmers took " << duration1.count() << " seconds" << std::endl;
#endif


        BplusTree<Run> bptree(15); // TODO: determine the BPlusTree degree




        cerr << "Adding the kmers to the BPlusTree" << endl;
//    kmers_to_bplustree(idx, bptree, index, k, {0, idx.bwt_size() - 1}, "");
        parallel_kmers_to_bplustree(idx, bptree, index, k, {0, idx.bwt_size() - 1});


#if TIME
        auto time4 = chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration3 = time4 - time2;
        std::cerr << "Adding kmers to bptree took " << duration3.count() << " seconds" << std::endl;
#endif


        // computing some statistics
        auto unique_kmers_size = index.size();
        auto bptree_items = bptree.get_bpt_size();
        auto bwt_size = idx.bwt_size();
        size_t tag_arrays_covered = 0;

        cerr << "The number of unique kmers in the index is: " << unique_kmers_size << endl;
        cerr << "The number of items in the BPlusTree is: " << bptree_items << endl;
        cerr << "The size of the BWT is: " << bwt_size << endl;


        cerr << "calculating the fraction of the tag arrays covered" << endl;
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


        cerr << "The fraction of the tag arrays covered by unique kmers is: " << tag_arrays_covered << " / " << bwt_size
             << " = "
             << (double) tag_arrays_covered / bwt_size << endl;


#if TIME
        auto time5 = chrono::high_resolution_clock::now();
#endif
        cerr << "Extending the kmers on the graph" << endl;
        auto extension_candidates = extend_kmers_bfs_parallel(gbz->graph, idx, bptree, 1024);

#if TIME
        auto time6 = chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration4 = time6 - time5;
        std::cerr << "Extending kmers took " << duration4.count() << " seconds" << std::endl;
#endif


//    cerr << "The number of items in the BPlusTree is: " << bptree.get_bpt_size() << endl;
        tag_arrays_covered = 0;

        cerr << "calculating the fraction of the tag arrays covered after one extension" << endl;
        // calculating the fraction of the tag arrays covered
        for (auto it = bptree.begin(); it != bptree.end(); ++it) {
            Run current_item = *it;
            if (current_item.graph_position.value != 0) { // Check if the current item is not a gap
                auto next_it = it;
                ++next_it; // Move to the next element
                if (next_it != bptree.end()) { // Check if the next element is not the end
                    Run next_item = *next_it; // Get the next item
                    tag_arrays_covered += (next_item.start_position - current_item.start_position);
//                cerr << "The current item is: " << current_item << " NEXT " << next_item << endl;
//                if (next_item.start_position - current_item.start_position > 500){
//                    cerr << "The current item is: " << current_item << " The next item is: " << next_item << endl;
//                }
                }
            }
        }

        cerr << "The fraction of the tag arrays covered after extending the kmers is: " << tag_arrays_covered << " / "
             << bwt_size << " = "
             << (double) tag_arrays_covered / bwt_size << endl;


        string pattern = "$";


        // get the ISA values for the end of each sequence by searching for pattern # and locating the results
        auto OCC = idx.ISA(pattern);

        // using the sort_end_of_seq function
        // the first item in the nth element of this is the end of the nth sequence in the bwt
        auto end_of_seq = sort_end_of_seq(OCC);

#if TIME
        auto time7 = chrono::high_resolution_clock::now();
#endif

        traverse_sequences_parallel(*gbz, bptree, idx, end_of_seq);
#if TIME
        auto time8 = chrono::high_resolution_clock::now();
        std::chrono::duration<double> duration5 = time8 - time7;
        std::cerr << "Traversing all paths and fill all the gaps took " << duration5.count() << " seconds" << std::endl;
#endif


        cerr << "The final number of items in the BPlusTree is: " << bptree.get_bpt_size() << endl;
        tag_arrays_covered = 0;

        cerr << "calculating the fraction of the tag arrays covered " << endl;
        // calculating the fraction of the tag arrays covered
        for (auto it = bptree.begin(); it != bptree.end(); ++it) {
            Run current_item = *it;
            if (current_item.graph_position.value != 0) { // Check if the current item is not a gap
                auto next_it = it;
                ++next_it; // Move to the next element
                if (next_it != bptree.end()) { // Check if the next element is not the end
                    Run next_item = *next_it; // Get the next item
                    tag_arrays_covered += (next_item.start_position - current_item.start_position);
                    // print the decoded items of current_item and next item
                    pos_t t1 = current_item.graph_position.decode();
//                cerr << "The current item is: " << current_item << endl;
//                cerr << "node id: " << vg::id(t1) << " offset " << vg::offset(t1) << " rev? " << vg::is_rev(t1) << endl;

//                if (next_item.start_position - current_item.start_position > 500){
//                    cerr << "The current item is: " << current_item << " The next item is: " << next_item << endl;
//                }
                }
            }
        }

        cerr << "The fraction of the tag arrays covered after filling the gaps is: " << tag_arrays_covered << " / "
             << bwt_size << " = "
             << (double) tag_arrays_covered / bwt_size << endl;

        vg::TagArray tag_array;
        tag_array.load_bptree(bptree, idx.bwt_size());

        tag_array.serialize(std::cout);
    } else {
        vg::TagArray tag_array;
        std::ifstream in_ds("/Users/seeskand/Documents/pangenome-index/test_data/1mb.ser");
        tag_array.load(in_ds);

//    tag_array.query(4, 9);
        std::default_random_engine generator(static_cast<unsigned>(std::time(0)));
        std::uniform_int_distribution<int> distribution(0, 3);
        auto generate_random_kmer = [&](int k, const std::string &alphabet) {
            std::string kmer;

            for (int i = 0; i < k; ++i) {
                kmer += alphabet[distribution(generator)];
            }

            return kmer;
        };

        int num_kmers = 1000000;
        std::string alphabet = "ACGT";

        std::vector<std::string> kmers;
        kmers.reserve(num_kmers);
        for (int i = 0; i < num_kmers; ++i) {
            kmers.push_back(generate_random_kmer(10, alphabet));
        }

        for (auto &kmer: kmers) {


            //find the interval in the bwt
            auto range = idx.count(kmer);
//        cout << range.first << " " << range.second << endl;



            if (range.first <= range.second) {
//            cerr << "The kmer is: " << kmer << endl;
//            cout << "bwt range is " << range.first << " " << range.second << endl;
                tag_array.query(range.first, range.second);
            }

        }
    }
//    exit(0);



    return 0;
}


static vg::subcommand::Subcommand vg_panindexer("panindexer", "index pangenomes", PIPELINE, main_panindexer);
