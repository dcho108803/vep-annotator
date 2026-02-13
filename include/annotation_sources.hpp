/**
 * Annotation Source Factory Functions
 *
 * Provides factory functions to create various annotation sources.
 */

#ifndef ANNOTATION_SOURCES_HPP
#define ANNOTATION_SOURCES_HPP

#include "annotation_source.hpp"
#include <string>
#include <memory>
#include <set>

namespace vep {

// Forward declarations
class ReferenceGenome;

// ============================================================================
// Pathogenicity Sources
// ============================================================================

/**
 * Create dbNSFP annotation source
 * @param path Path to tabix-indexed dbNSFP file (.txt.gz + .tbi)
 * @param fields Field specification ("essential", "all", or comma-separated list)
 */
std::shared_ptr<AnnotationSource> create_dbnsfp_source(
    const std::string& path,
    const std::string& fields = "essential"
);

// ============================================================================
// Splice Prediction Sources
// ============================================================================

/**
 * Create SpliceAI annotation source
 * @param path Path to SpliceAI VCF file (.vcf.gz + .tbi)
 */
std::shared_ptr<AnnotationSource> create_spliceai_source(const std::string& path);

/**
 * Create MaxEntScan annotation source (algorithmic, no data file needed)
 */
std::shared_ptr<AnnotationSource> create_maxentscan_source();

/**
 * Create MaxEntScan annotation source with reference genome for actual scoring
 */
std::shared_ptr<AnnotationSource> create_maxentscan_source(const ReferenceGenome* ref);

/**
 * Create dbscSNV annotation source
 * @param path Path to dbscSNV VCF file (.vcf.gz + .tbi)
 */
std::shared_ptr<AnnotationSource> create_dbscsnv_source(const std::string& path);

// ============================================================================
// Conservation Sources
// ============================================================================

/**
 * Create PhyloP conservation score source
 * @param path Path to PhyloP bigWig file (.bw)
 */
std::shared_ptr<AnnotationSource> create_phylop_source(const std::string& path);

/**
 * Create PhastCons conservation score source
 * @param path Path to PhastCons bigWig file (.bw)
 */
std::shared_ptr<AnnotationSource> create_phastcons_source(const std::string& path);

/**
 * Create GERP++ conservation score source
 * @param path Path to GERP++ bigWig file (.bw)
 */
std::shared_ptr<AnnotationSource> create_gerp_source(const std::string& path);

// ============================================================================
// Regulatory Sources
// ============================================================================

/**
 * Create Ensembl Regulatory Build annotation source
 * @param path Path to regulatory GFF3 file
 * @param feature_types Feature types to include (empty = all)
 */
std::shared_ptr<AnnotationSource> create_regulatory_source(
    const std::string& path,
    const std::set<std::string>& feature_types = {}
);

// ============================================================================
// Protein Domain Sources
// ============================================================================

/**
 * Create Pfam domain annotation source
 * @param path Path to Pfam TSV file
 */
std::shared_ptr<AnnotationSource> create_pfam_source(const std::string& path);

/**
 * Create InterPro domain annotation source
 * @param path Path to InterPro TSV file
 */
std::shared_ptr<AnnotationSource> create_interpro_source(const std::string& path);

// ============================================================================
// LoF Sources
// ============================================================================

/**
 * Create LOFTEE-style LoF annotation source
 */
std::shared_ptr<AnnotationSource> create_loftee_source();

/**
 * Create NMD prediction source
 */
std::shared_ptr<AnnotationSource> create_nmd_source();

/**
 * Create LoFtool gene tolerance source
 * @param path Path to LoFtool scores file
 */
std::shared_ptr<AnnotationSource> create_loftool_source(const std::string& path);

} // namespace vep

#endif // ANNOTATION_SOURCES_HPP
