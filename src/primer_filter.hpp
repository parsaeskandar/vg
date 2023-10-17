// primer_filter.hpp
//
// Contains class Primer_finder for storing and filtering primers predicted 
// using Primer3. Also contains Primer struct and Primer_pair struct that stores
// information on primers and primer pairs.

#ifndef VG_PRIMER_FILTER_HPP_INCLUDED
#define VG_PRIMER_FILTER_HPP_INCLUDED

// Not sure what to include.. Will include everything from the unittest for now
#include <stdio.h>
#include <iostream>
#include <regex>
#include <fstream>
#include <vector>
#include <sstream>
#include <set>
#include <vg/vg.pb.h>
#include "utility.hpp"
#include "snarl_distance_index.hpp"
#include "integrated_snarl_finder.hpp"
#include "genotypekit.hpp"
#include "traversal_finder.hpp"
#include <vg/io/protobuf_emitter.hpp>
#include <vg/io/vpkg.hpp>

using namespace std;

namespace vg {

/**
 * Primer struct contains primer attributes, including sequence, left/right primer,
 * position on the reference genome, length, index offset in corresponding node on
 * sequence graph, and vector of corresponding nodes on the sequence graph. Everything
 * is in the positive/forward orientation.
 */
struct Primer {
    string sequence;
    bool left = true;
    size_t position = numeric_limits<size_t>::max();
    size_t length = numeric_limits<size_t>::max();
    size_t offset = numeric_limits<size_t>::max();
    vector<size_t> mapped_nodes_ids;
};

/**
 * Primer_pair struct contains primer pair attributesm including left primer, right primer,
 * linear product size, minimum and maximum product size on the sequence graph, and boolean on
 * whether the primers locate in low variation region of the sequence graph.
 */
struct PrimerPair {
    Primer left_primer;
    Primer right_primer;
    size_t linear_product_size = numeric_limits<size_t>::max();
    size_t min_product_size = numeric_limits<size_t>::max();
    size_t max_product_size = numeric_limits<size_t>::max();
    bool no_variation = false;
};

class PrimerFinder {

private:
    vector<PrimerPair> primer_pairs;
    vector<PrimerPair> selected_primer_pairs;
    const PathPositionHandleGraph* graph;
    const SnarlDistanceIndex* distance_index;
    path_handle_t reference_path_handle; 

public:
    PrimerFinder() = default;
    
    /**
     * Construct Primer finder given PathPositionHandleGraph, reference graph name
     * and pointer to SnarlDistanceIndex
     */
    PrimerFinder(const unique_ptr<handlegraph::PathPositionHandleGraph>& graph_param,
                const string& reference_path_name, const SnarlDistanceIndex* distance_index_param);

    /**
     * Destructor
     */
    ~PrimerFinder();

    /**
     * Add a Primer_pair object given primers' starting node id, offset relative
     * to the starting node, and length, all in the POSTIVE orientation. The new
     * primer_pair object is automatically added to primer_pairs vector - and
     * selected_primer_pairs if conditions are met. Mainly used for unit testing.
     */
    void add_primer_pair(const size_t& left_primer_starting_node_id,
                    const size_t& left_primer_offset, const size_t& left_primer_length,
                    const size_t& right_primer_starting_node_id,
                    const size_t& right_primer_offset, const size_t& right_primer_length);

    /**
     * Read the path to the primer3 output. Primers information is parsed,
     * processed, and  stored in primer_pairs vector - and selected_primer_pairs
     * if conditions are met.
     */
    void load_primers(const string& path_to_primers);

    /**
     * return vector of Primer pairs
     */
    const vector<PrimerPair>& get_primer_pairs() const;

    /**
     * return vector selected primer pairs
     */
    const vector<PrimerPair>& get_selected_primer_pairs() const;

private:
    /**
     * Private functions used by public or private functions.
     */

    /**
     * Update minimum and maximum prodcut to a primer pair object.
     * Used in: add_primer_pair
     *          load_primers
     */
    void update_min_max_product_size(PrimerPair& primer_pair);

    /**
     * Update a Primer object given starting node id, offset relative to the starting node,
     * and the length of primer.
     * Used in: add_primer_pair
     */
    void make_primer(Primer& primer, const size_t& starting_node_id,
            const size_t& offset, const size_t& length, const bool& is_left);

    /**
     * Find and store corresponding node ids to Primer object.
     * Used in: make_primer
     *          load_primers
     */
    void map_to_nodes(Primer& primer);

    /**
     * Find the length of the longest match between two sequences. Also find and
     * store offset in Primer object.
     * Used in: map_to_nodes
     */
    size_t longest_match_len(Primer& primer, const string& left_seq, const string& right_seq,
       const bool& first_node);

    /**
     * Strip empty spaces on the right side of a string.
     * Used in: load_primers
     */
    const string rstrip(const string& s) const;

    /**
     * Check if primers in a primer_pair object have variations on the pangenome.
     * Used in: add_primer_node
     *          load_primers
     */
    const bool no_variation(const PrimerPair& primer_pair) const;

    /**
     * Split a string into vectors.
     * Used in: load_priemrs
     */
    const vector<string> split(string str, const string& delim) const;

};

}

#endif /* primder_filter_hpp */