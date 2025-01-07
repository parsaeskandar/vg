/**
 * \file haplotype_indexer.cpp: implementations of haplotype indexing with the GBWT
 */

#include "haplotype_indexer.hpp"

#include <iostream>
#include <mutex>
#include <vector>
#include <string>

#include <vg/io/stream.hpp>
#include <gbwtgraph/algorithms.h>
#include <gbwtgraph/path_cover.h>

#include "gbwt_helper.hpp"

#include "alignment.hpp"
#include "hash_map.hpp"
#include "path.hpp"

using namespace std;

namespace vg {

HaplotypeIndexer::HaplotypeIndexer() {
    // Build the GBWTs silently.
    gbwt::Verbosity::set(gbwt::Verbosity::SILENT);
}

std::vector<std::string> HaplotypeIndexer::parse_vcf(const std::string& filename, const PathHandleGraph& graph, const std::string& job_name) const {

    // Parse all non-alt paths.
    std::vector<path_handle_t> path_handles;
    graph.for_each_path_handle([&](path_handle_t path_handle) {
        std::string path_name = graph.get_path_name(path_handle);
        if (!Paths::is_alt(path_name)) {
            path_handles.push_back(path_handle);
        }
    });

    return this->parse_vcf(filename, graph, path_handles, job_name);
}

std::vector<std::string> HaplotypeIndexer::parse_vcf(const std::string& filename, const PathHandleGraph& graph, const std::vector<path_handle_t>& paths, const std::string& job_name) const {

    // Open the VCF file.
    vcflib::VariantCallFile variant_file;
    variant_file.parseSamples = false; // vcflib parsing is very slow if there are many samples.
    std::string temp_filename = filename;
    variant_file.open(temp_filename);
    if (!variant_file.is_open()) {
        std::cerr << "error: [HaplotypeIndexer::parse_vcf] could not open " << filename << std::endl;
        std::exit(EXIT_FAILURE);
    }

    // How many samples are there?
    size_t num_samples = variant_file.sampleNames.size();
    if (num_samples == 0) {
        std::cerr << "error: [HaplotypeIndexer::parse_vcf] variant file '" << filename << "' does not contain phasings" << std::endl;
        std::exit(EXIT_FAILURE);
    }

    // Determine the samples we want to index.
    std::pair<size_t, size_t> sample_range = this->sample_range;
    sample_range.second = std::min(sample_range.second, num_samples);
    std::vector<std::string> sample_names(variant_file.sampleNames.begin() + sample_range.first, variant_file.sampleNames.begin() + sample_range.second);
    if (this->show_progress) {
        #pragma omp critical
        {
            std::cerr << job_name << ": Parsing VCF file " << filename << " with options";
            if (!this->phase_homozygous) {
                std::cerr << " --actual-phasing";
            }
            if (this->force_phasing) {
                std::cerr << " --force-phasing";
            }
            if (this->discard_overlaps) {
                std::cerr << " --discard-overlaps";
            }
            if (!this->rename_variants) {
                std::cerr << "  --vcf-variants";
            }
            std::cerr << std::endl;
            std::cerr << job_name << ": Samples " << sample_range.first << " to " << (sample_range.second - 1) << ", batch size " << samples_in_batch << std::endl;
        }
    }

    // Parse the contigs we are interested in.
    std::vector<std::string> result;
    size_t total_variants_processed = 0;
    std::mt19937 rng(0xDEADBEEF);
    std::uniform_int_distribution<std::mt19937::result_type> random_bit(0, 1);
    size_t found_missing_variants = 0;
    for (size_t path_id = 0; path_id < paths.size(); path_id++) {
        std::string path_name = graph.get_path_name(paths[path_id]);
        std::string vcf_contig_name = (this->path_to_vcf.count(path_name) > 0 ? this->path_to_vcf.at(path_name) : path_name);

        // Set the VCF region or process the entire contig.
        if (this->regions.count(vcf_contig_name)) {
            std::pair<size_t, size_t> region = this->regions.at(vcf_contig_name);
            variant_file.setRegion(vcf_contig_name, region.first, region.second);
        } else {
            variant_file.setRegion(vcf_contig_name);
        }

        // Check that the VCF file contains this contig.
        vcflib::Variant var(variant_file);
        if (!(variant_file.is_open() && variant_file.getNextVariant(var) && var.sequenceName == vcf_contig_name)) {
            std::cerr << "warning: [HaplotypeIndexer::parse_vcf] contig " << vcf_contig_name << " not present in file " << filename << std::endl;
            continue;
        }
        if (this->show_progress) {
            #pragma omp critical
            {
                std::cerr << job_name << ": Path " << path_name << " matches VCF contig " << vcf_contig_name;
                if (this->regions.count(vcf_contig_name)) {
                    std::pair<size_t, size_t> region = this->regions.at(vcf_contig_name);
                    std::cerr << ", region " << region.first << " to " << region.second;
                }
                std::cerr << std::endl;
            }
        }

        // Structures to parse the VCF file into.
        std::string parse_file = (this->batch_file_prefix.empty() ? gbwt::TempFile::getName("parse") : this->batch_file_prefix + '_' + vcf_contig_name);
        gbwt::VariantPaths variants(graph.get_step_count(paths[path_id]));
        variants.setSampleNames(sample_names);
        variants.setContigName(path_name);
        std::vector<gbwt::PhasingInformation> phasings;

        // Add the reference to VariantPaths.
        for (handle_t handle : graph.scan_path(paths[path_id])) {
            variants.appendToReference(gbwt::Node::encode(graph.get_id(handle), graph.get_is_reverse(handle)));
        }
        variants.indexReference();

        // Create a PhasingInformation file for each batch.
        for (size_t batch_start = sample_range.first; batch_start < sample_range.second; batch_start += samples_in_batch) {
            size_t batch_size = std::min(samples_in_batch, sample_range.second - batch_start);
            if (!this->batch_file_prefix.empty()) {
                // Use a permanent file.
                phasings.emplace_back(parse_file, batch_start, batch_size);
            } else {
                // Use a temporary file that persists until the program exits.
                phasings.emplace_back(batch_start, batch_size);
                phasings.back().makePersistent();
            }
            variants.addFile(phasings.back().name(), phasings.back().offset(), phasings.back().size());
        }

        // Parse the variants and the phasings.
        size_t variants_processed = 0;
        std::vector<bool> was_diploid(sample_range.second, true); // Was the sample diploid at the previous site?
        do {
            // Skip variants with non-DNA sequence, as they are not included in the graph.
            bool isDNA = allATGC(var.ref);
            for (std::vector<std::string>::iterator a = var.alt.begin(); a != var.alt.end(); ++a) {
                if (!allATGC(*a)) isDNA = false;
            }
            if (!isDNA) {
                continue;
            }
            
            if (this->rename_variants) {
                // We need to move the variant over to the contig name
                // used in the graph, in order to get the right id for
                // it in the graph.
                var.sequenceName = path_name;
            }
            
            // Determine the reference nodes for the current variant and create a variant site.
            // If the variant is not an insertion, there should be a path for the ref allele.
            std::string var_name = make_variant_id(var);
            std::string ref_path_name = "_alt_" + var_name + "_0";
            gbwt::vector_type ref_path = extract_as_gbwt_path(graph, ref_path_name);
            
            size_t ref_pos = variants.invalid_position();
            if (!ref_path.empty()) {
                ref_pos = variants.firstOccurrence(ref_path.front());
                if (ref_pos == variants.invalid_position()) {
                    #pragma omp critical
                    {
                        std::cerr << "warning: [HaplotypeIndexer::parse_vcf] invalid ref path for " << var_name << " at "
                            << var.sequenceName << ":" << var.position << std::endl;
                    }
                    continue;
                }
            } else {
                // Try using the alternate alleles instead.
                bool found = false;
                for (size_t alt_index = 1; alt_index < var.alleles.size(); alt_index++) {
                    std::string alt_path_name = "_alt_" + var_name + "_" + std::to_string(alt_index);
                    size_t candidate_pos = 0;
                    bool candidate_found = false;
                    gbwt::vector_type pred_nodes = path_predecessors(graph, alt_path_name);
                    if (!pred_nodes.empty()) {
                        for (auto node : pred_nodes) {
                            size_t pred_pos = variants.firstOccurrence(node);
                            if (pred_pos != variants.invalid_position()) {
                                candidate_pos = std::max(candidate_pos, pred_pos + 1);
                                candidate_found = true;
                                found = true;
                            }
                        }
                        // For each alternate allele, find the rightmost reference node among
                        // its predecessors. If multiple alleles have candidates for the
                        // reference position, choose the leftmost one.
                        if (candidate_found) {
                            ref_pos = std::min(ref_pos, candidate_pos);
                        }
                    }
                }
                if (!found) {
                    // This variant from the VCF is just not in the graph, so skip it.
                    found_missing_variants++;
                    if (this->warn_on_missing_variants && found_missing_variants <= this->max_missing_variant_warnings) {
                        #pragma omp critical
                        {
                            // The user might not know it. Warn them in case they mixed up their VCFs.
                            std::cerr << "warning: [HaplotypeIndexer::parse_vcf] alt and ref paths for " << var_name
                                    << " at " << var.sequenceName << ":" << var.position
                                    << " missing/empty! Was the variant skipped during construction?" << std::endl;
                            if (found_missing_variants == this->max_missing_variant_warnings) {
                                std::cerr << "warning: [HaplotypeIndexer::parse_vcf] suppressing further missing variant warnings" << std::endl;
                            }
                        }
                    }
                    continue;
                }
            }
            variants.addSite(ref_pos, ref_pos + ref_path.size());

            // Add alternate alleles to the site.
            for (size_t alt_index = 1; alt_index < var.alleles.size(); alt_index++) {
                std::string alt_path_name = "_alt_" + var_name + "_" + std::to_string(alt_index);
                variants.addAllele(extract_as_gbwt_path(graph, alt_path_name));
            }

            // Store the phasings in PhasingInformation structures.
            std::vector<std::string> genotypes = parseGenotypes(var.originalLine, num_samples);
            for (size_t batch = 0; batch < phasings.size(); batch++) {
                std::vector<gbwt::Phasing> current_phasings;
                for (size_t sample = phasings[batch].offset(); sample < phasings[batch].limit(); sample++) {
                    current_phasings.emplace_back(genotypes[sample], was_diploid[sample], this->phase_homozygous);
                    was_diploid[sample] = current_phasings.back().diploid;
                    if(this->force_phasing) {
                        current_phasings.back().forcePhased([&]() {
                            return random_bit(rng);
                        });
                    }
                }
                phasings[batch].append(current_phasings);
            }
            variants_processed++;
        }
        while (variant_file.is_open() && variant_file.getNextVariant(var) && var.sequenceName == vcf_contig_name); // End of variants.
        if (this->show_progress) {
            size_t phasing_bytes = 0;
            for (size_t batch = 0; batch < phasings.size(); batch++) {
                phasing_bytes += phasings[batch].bytes();
            }
            #pragma omp critical
            {
                std::cerr << job_name << ": Processed " << variants_processed << " variants on path " << path_name << ", " << gbwt::inMegabytes(phasing_bytes) << " MiB phasing information" << std::endl;
                std::cerr << job_name << ": Saving the VCF parse for path " << path_name << " to " << parse_file << std::endl;
            }
        }

        // Save the VCF parse.
        if (!sdsl::store_to_file(variants, parse_file)) {
            std::cerr << "error: [HaplotypeIndexer::parse_vcf] cannot write parse file " << parse_file << std::endl;
            std::exit(EXIT_FAILURE);
        }
        result.push_back(parse_file);

        // End of haplotype generation for the current contig.
        total_variants_processed += variants_processed;
    } // End of contigs.
        
    if (this->warn_on_missing_variants && found_missing_variants > 0) {
        #pragma omp critical
        {
            std::cerr << "warning: [HaplotypeIndexer::parse_vcf] Found " << found_missing_variants << "/" << total_variants_processed
                << " variants in phasing VCF but not in graph! Do your graph and VCF match?" << std::endl;
        }
    }

    return result;
}

std::unique_ptr<gbwt::DynamicGBWT> HaplotypeIndexer::build_gbwt(const std::vector<std::string>& vcf_parse_files,
                                                                const std::string& job_name,
                                                                const PathHandleGraph* graph,
                                                                const std::unordered_set<path_handle_t>* paths,
                                                                bool skip_unvisited_paths) const {

    // GBWT index.
    std::unique_ptr<gbwt::DynamicGBWT> index(new gbwt::DynamicGBWT());
    index->addMetadata();
    if (vcf_parse_files.empty() && !graph) {
        return index;
    }
    
    {
        // New stuff we're adding to the GBWT metadata
        std::set<gbwt::range_type> haplotypes;
        std::vector<std::string> sample_names, contig_names;
        
        // Haplotype construction for each contig from VCF.
        for (const std::string& filename : vcf_parse_files) {
            gbwt::VariantPaths variants;
            if (!sdsl::load_from_file(variants, filename)) {
                std::cerr << "error: [HaplotypeIndexer::build_gbwt] cannot load VCF parse from " << filename << std::endl;
                std::exit(EXIT_FAILURE);
            }
            if (!variants.hasContigName() || !variants.hasSampleNames()) {
                std::cerr << "error: [HaplotypeIndexer::build_gbwt] VCF parse file " << filename << " does not contain sample/contig names" << std::endl;
                std::exit(EXIT_FAILURE);
            }
            if (sample_names.empty()) {
                sample_names = variants.getSampleNames();
            } else if (sample_names != variants.getSampleNames()) {
                std::cerr << "error: [HaplotypeIndexer::build_gbwt] invalid sample names in VCF parse file " << filename << std::endl;
                std::exit(EXIT_FAILURE);
            }
            contig_names.emplace_back(variants.getContigName());
            if (this->show_progress) {
                #pragma omp critical
                {
                    std::cerr << job_name << ": Generating haplotypes for path " << variants.getContigName() << " from file " << filename << std::endl;
                }
            }
            gbwt::GBWTBuilder builder(variants.nodeWidth(true), this->gbwt_buffer_size * gbwt::MILLION, this->id_interval);
            builder.swapIndex(*index);
            gbwt::generateHaplotypes(variants, std::set<std::string>(),
            [&](gbwt::size_type sample_id) -> bool {
                return (this->excluded_samples.find(sample_names[sample_id]) == this->excluded_samples.end());
            }, [&](const gbwt::Haplotype& haplotype) {
                builder.insert(haplotype.path, true); // Insert in both orientations.
                builder.index.metadata.addPath(haplotype.sample, contig_names.size() - 1, haplotype.phase, haplotype.count);
                haplotypes.insert(gbwt::range_type(haplotype.sample, haplotype.phase));
            }, [&](gbwt::size_type, gbwt::size_type) -> bool {
                // For each overlap, discard it if our global flag is set.
                return this->discard_overlaps;
            });
            builder.finish();
            builder.swapIndex(*index);
        }
        
        // Now count the haplotypes we added
        index->metadata.setHaplotypes(index->metadata.haplotypes() + haplotypes.size());
        // And remember the samples and contigs we created
        index->metadata.setSamples(sample_names);
        index->metadata.setContigs(contig_names);
    }
    
    if (graph) {
        // Also include graph named paths from the graph
        
        // Actual work.
        if (show_progress) {
            #pragma omp critical
            {
                std::cerr << "Indexing embedded paths" << std::endl;
            }
        }
        
        // GBWT construction, into existing cumulative index.
        gbwt::GBWTBuilder builder(gbwt_node_width(*graph), this->gbwt_buffer_size * gbwt::MILLION, this->id_interval);
        builder.swapIndex(*index);
        
        std::unordered_set<std::string> visited_contig_names;
        if (skip_unvisited_paths) {
            // Make a set of the contigs visited so we can filter down to them.
            for (size_t i = 0; i < builder.index.metadata.contig_names.size(); i++) {
                visited_contig_names.insert(builder.index.metadata.contig(i));
            }
        }
        
        // Set up a filter to drop alt paths and other paths we don't want.
        std::function<bool(const path_handle_t&)> path_filter = [&](const path_handle_t& path_handle) {
            std::string path_name = graph->get_path_name(path_handle);
            if (Paths::is_alt(path_name)) {
                return false;
            }
            if (skip_unvisited_paths && !visited_contig_names.count(graph->get_locus_name(path_handle))) {
                // This path is on an unvisited contig
                return false;
            }
            if (paths && !paths->count(path_handle)) {
                return false;
            }
            return true;
        };
        // And add every path that passes the filter (including haplotype paths) from the source graph.
        gbwtgraph::store_paths(builder, *graph, {PathSense::GENERIC, PathSense::REFERENCE, PathSense::HAPLOTYPE}, &path_filter);

        // Finish the construction for this set of paths and put the index back.
        builder.finish();
        builder.swapIndex(*index);
    }

    // Finish the construction overall.
    if (this->show_progress) {
        std::cerr << job_name << ": ";
        gbwt::operator<<(std::cerr, index->metadata);
        std::cerr << std::endl;
    }
    return index;
}

std::unique_ptr<gbwt::DynamicGBWT> HaplotypeIndexer::build_gbwt(const PathHandleGraph& graph) const {
    // Fall back to the general vcf-and-graph implementation
    return build_gbwt({}, "GBWT", &graph);
}

std::unique_ptr<gbwt::GBWT> HaplotypeIndexer::build_gbwt(const HandleGraph& graph,
    const std::vector<std::string>& aln_filenames, const std::string& aln_format, size_t parallel_jobs) const {

    // Handle multithreading and parallel jobs.
    parallel_jobs = std::max<size_t>(1, parallel_jobs);
    int old_threads = omp_get_max_threads();
    omp_set_num_threads(parallel_jobs);

    // Partition the graph into construction jobs.
    if (this->show_progress) {
        #pragma omp critical
        {
            std::cerr << "Partitioning the graph into GBWT construction jobs" << std::endl;
        }
    }
    size_t target_size = graph.get_node_count() / parallel_jobs;
    gbwtgraph::ConstructionJobs jobs = gbwtgraph::gbwt_construction_jobs(graph, target_size);
    if (this->show_progress) {
        #pragma omp critical
        {
            std::cerr << "Created " << jobs.size() << " parallel construction jobs" << std::endl;
        }
    }

    // GBWT construction.
    std::vector<std::mutex> builder_mutexes(jobs.size());
    std::vector<std::unique_ptr<gbwt::GBWTBuilder>> builders(jobs.size());
    // This is a bit inefficient, as read names are often longer than the SSO threshold for GCC (but not for Clang).
    // TODO: Maybe use concatenated 0-terminated names?
    std::vector<std::vector<std::string>> read_names(jobs.size());
    for (size_t i = 0; i < jobs.size(); i++) {
        builders[i].reset(new gbwt::GBWTBuilder(gbwt_node_width(graph), this->gbwt_buffer_size * gbwt::MILLION, this->id_interval));
    }

    // Actual work.
    if (this->show_progress) {
        #pragma omp critical
        {
            std::cerr << "Converting " << aln_format << " to GBWT paths" << std::endl;
        }
    }
    std::function<void(Alignment&)> lambda = [&](Alignment& aln) {
        gbwt::vector_type buffer;
        for (auto& m : aln.path().mapping()) {
            buffer.push_back(mapping_to_gbwt(m));
        }
        size_t job_id = 0;
        if (buffer.size() > 0) {
            job_id = jobs.job(gbwt::Node::id(buffer.front()));
            if (job_id >= jobs.size()) {
                job_id = 0;
            }
        }
        {
            // Insert the path into the appropriate builder and record the read name.
            std::lock_guard<std::mutex> lock(builder_mutexes[job_id]);
            builders[job_id]->insert(buffer, true);
            read_names[job_id].push_back(aln.name());
        }
    };
    for (auto& file_name : aln_filenames) {
        if (aln_format == "GAM") {
            get_input_file(file_name, [&](istream& in) {
                vg::io::for_each_parallel(in, lambda);
            });
        } else {
            assert(aln_format == "GAF");
            vg::io::gaf_unpaired_for_each_parallel(graph, file_name, lambda);
        }
    }

    // Finish the construction and convert to compressed GBWT.
    std::vector<gbwt::GBWT> partial_indexes(jobs.size());
    #pragma omp parallel for schedule(dynamic, 1)
    for (size_t i = 0; i < jobs.size(); i++) {
        builders[i]->finish();
        partial_indexes[i] = gbwt::GBWT(builders[i]->index);
        builders[i].reset();
    }

    // Merge the partial indexes.
    if (this->show_progress) {
        #pragma omp critical
        {
            std::cerr << "Merging the partial GBWTs" << std::endl;
        }
    }
    std::unique_ptr<gbwt::GBWT> result(new gbwt::GBWT(partial_indexes));
    partial_indexes.clear();

    // Create the metadata.
    if (this->show_progress) {
        #pragma omp critical
        {
            std::cerr << "Creating metadata" << std::endl;
        }
    }
    result->addMetadata();
    result->metadata.setContigs({ "unknown" });
    {
        // We can use 32-bit values, as GBWT metadata uses them as well.
        string_hash_map<std::string, std::pair<std::uint32_t, std::uint32_t>> read_info; // name -> (sample id, fragment count)
        for (auto& names : read_names) {
            for (const std::string& name : names) {
                std::uint32_t sample_id = 0, fragment_count = 0;
                auto iter = read_info.find(name);
                if (iter == read_info.end()) {
                    sample_id = read_info.size();
                    read_info[name] = std::make_pair(sample_id, fragment_count);
                } else {
                    sample_id = iter->second.first;
                    fragment_count = iter->second.second;
                    iter->second.second++;
                }
                result->metadata.addPath(sample_id, 0, 0, fragment_count);
            }
            names = std::vector<std::string>();
        }
        std::vector<std::string> sample_names(read_info.size());
        for (auto& p : read_info) {
            sample_names[p.second.first] = p.first;
        }
        read_info = string_hash_map<std::string, std::pair<std::uint32_t, std::uint32_t>>();
        result->metadata.setSamples(sample_names);
        result->metadata.setHaplotypes(sample_names.size());
    }
    if (this->show_progress) {
        #pragma omp critical
        {
            std::cerr << "GBWT: ";
            gbwt::operator<<(std::cerr, result->metadata);
            std::cerr << std::endl;
        }
    }

    // Restore the number of threads.
    omp_set_num_threads(old_threads);
    return result;
}

}
