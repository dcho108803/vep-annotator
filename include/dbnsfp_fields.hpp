/**
 * dbNSFP Field Definitions
 *
 * Defines all available fields from dbNSFP database and their column mappings.
 * dbNSFP includes 35+ pathogenicity prediction scores and conservation metrics.
 */

#ifndef DBNSFP_FIELDS_HPP
#define DBNSFP_FIELDS_HPP

#include <string>
#include <vector>
#include <map>
#include <set>
#include <sstream>

namespace vep {

/**
 * dbNSFP field definition
 */
struct DbNSFPField {
    std::string name;           // Output field name
    std::string column;         // Column name in dbNSFP file
    std::string description;    // Field description
    bool is_score;              // True if numeric score
    bool higher_is_damaging;    // For scores: true if higher = more damaging
    double damaging_threshold;  // Threshold for "damaging" prediction (-1 = no threshold)
};

/**
 * Standard dbNSFP fields organized by category
 */

// Pathogenicity Prediction Scores
inline const std::vector<DbNSFPField> DBNSFP_PATHOGENICITY_FIELDS = {
    // SIFT
    {"SIFT_score", "SIFT_score", "SIFT score (0-1, lower = more damaging)", true, false, 0.05},
    {"SIFT_pred", "SIFT_pred", "SIFT prediction (D=Damaging, T=Tolerated)", false, false, -1},
    {"SIFT4G_score", "SIFT4G_score", "SIFT4G score", true, false, 0.05},
    {"SIFT4G_pred", "SIFT4G_pred", "SIFT4G prediction", false, false, -1},

    // PolyPhen-2
    {"Polyphen2_HDIV_score", "Polyphen2_HDIV_score", "PolyPhen2 HDIV score (0-1, higher = more damaging)", true, true, 0.453},
    {"Polyphen2_HDIV_pred", "Polyphen2_HDIV_pred", "PolyPhen2 HDIV prediction (D/P/B)", false, false, -1},
    {"Polyphen2_HVAR_score", "Polyphen2_HVAR_score", "PolyPhen2 HVAR score (0-1, higher = more damaging)", true, true, 0.446},
    {"Polyphen2_HVAR_pred", "Polyphen2_HVAR_pred", "PolyPhen2 HVAR prediction (D/P/B)", false, false, -1},

    // CADD
    {"CADD_raw", "CADD_raw", "CADD raw score", true, true, -1},
    {"CADD_phred", "CADD_phred", "CADD phred-scaled score (>20 = top 1%)", true, true, 20.0},

    // REVEL
    {"REVEL_score", "REVEL_score", "REVEL ensemble score (0-1, higher = more damaging)", true, true, 0.5},

    // AlphaMissense
    {"AlphaMissense_score", "AlphaMissense_score", "AlphaMissense score (0-1)", true, true, 0.564},
    {"AlphaMissense_pred", "AlphaMissense_class", "AlphaMissense prediction (likely_pathogenic/ambiguous/likely_benign)", false, false, -1},

    // MetaSVM/MetaLR
    {"MetaSVM_score", "MetaSVM_score", "MetaSVM score", true, true, 0.0},
    {"MetaSVM_pred", "MetaSVM_pred", "MetaSVM prediction (D/T)", false, false, -1},
    {"MetaLR_score", "MetaLR_score", "MetaLR score", true, true, 0.5},
    {"MetaLR_pred", "MetaLR_pred", "MetaLR prediction (D/T)", false, false, -1},
    {"MetaRNN_score", "MetaRNN_score", "MetaRNN score", true, true, 0.5},
    {"MetaRNN_pred", "MetaRNN_pred", "MetaRNN prediction", false, false, -1},

    // VEST4
    {"VEST4_score", "VEST4_score", "VEST4 score (0-1, higher = more damaging)", true, true, 0.5},

    // PROVEAN
    {"PROVEAN_score", "PROVEAN_score", "PROVEAN score (lower = more damaging)", true, false, -2.5},
    {"PROVEAN_pred", "PROVEAN_pred", "PROVEAN prediction (D/N)", false, false, -1},

    // FATHMM
    {"FATHMM_score", "FATHMM_score", "FATHMM score (lower = more damaging)", true, false, -1.5},
    {"FATHMM_pred", "FATHMM_pred", "FATHMM prediction (D/T)", false, false, -1},

    // MutationTaster
    {"MutationTaster_score", "MutationTaster_score", "MutationTaster probability", true, true, 0.5},
    {"MutationTaster_pred", "MutationTaster_pred", "MutationTaster prediction (A/D/N/P)", false, false, -1},

    // MutationAssessor
    {"MutationAssessor_score", "MutationAssessor_score", "MutationAssessor score", true, true, 1.935},
    {"MutationAssessor_pred", "MutationAssessor_pred", "MutationAssessor prediction (H/M/L/N)", false, false, -1},

    // LRT
    {"LRT_score", "LRT_score", "LRT score", true, false, -1},
    {"LRT_pred", "LRT_pred", "LRT prediction (D/N/U)", false, false, -1},

    // DANN
    {"DANN_score", "DANN_score", "DANN score (0-1, higher = more damaging)", true, true, 0.95},

    // Eigen
    {"Eigen_raw", "Eigen-raw_coding", "Eigen raw score", true, true, -1},
    {"Eigen_phred", "Eigen-phred_coding", "Eigen phred score", true, true, -1},
    {"Eigen_PC_raw", "Eigen-PC-raw_coding", "Eigen-PC raw score", true, true, -1},
    {"Eigen_PC_phred", "Eigen-PC-phred_coding", "Eigen-PC phred score", true, true, -1},

    // M-CAP
    {"M_CAP_score", "M-CAP_score", "M-CAP score", true, true, 0.025},
    {"M_CAP_pred", "M-CAP_pred", "M-CAP prediction (D/T)", false, false, -1},

    // MPC
    {"MPC_score", "MPC_score", "MPC (Missense badness, PolyPhen-2, Constraint) score", true, true, 2.0},

    // PrimateAI
    {"PrimateAI_score", "PrimateAI_score", "PrimateAI score (0-1)", true, true, 0.803},
    {"PrimateAI_pred", "PrimateAI_pred", "PrimateAI prediction (D/T)", false, false, -1},

    // BayesDel
    {"BayesDel_addAF_score", "BayesDel_addAF_score", "BayesDel score with AF", true, true, 0.0692},
    {"BayesDel_addAF_pred", "BayesDel_addAF_pred", "BayesDel prediction with AF", false, false, -1},
    {"BayesDel_noAF_score", "BayesDel_noAF_score", "BayesDel score without AF", true, true, -0.0570},
    {"BayesDel_noAF_pred", "BayesDel_noAF_pred", "BayesDel prediction without AF", false, false, -1},

    // ClinPred
    {"ClinPred_score", "ClinPred_score", "ClinPred score", true, true, 0.5},
    {"ClinPred_pred", "ClinPred_pred", "ClinPred prediction (D/T)", false, false, -1},

    // LIST-S2
    {"LIST_S2_score", "LIST-S2_score", "LIST-S2 score", true, true, 0.85},
    {"LIST_S2_pred", "LIST-S2_pred", "LIST-S2 prediction (D/T)", false, false, -1},

    // DEOGEN2
    {"DEOGEN2_score", "DEOGEN2_score", "DEOGEN2 score", true, true, 0.5},
    {"DEOGEN2_pred", "DEOGEN2_pred", "DEOGEN2 prediction (D/T)", false, false, -1},

    // MVP
    {"MVP_score", "MVP_score", "MVP (Missense Variant Pathogenicity) score", true, true, 0.7},
    {"MVP_pred", "MVP_rankscore", "MVP rank score", true, true, -1},

    // gMVP
    {"gMVP_score", "gMVP_score", "gMVP score", true, true, 0.5},
};

// Conservation Scores (from dbNSFP)
inline const std::vector<DbNSFPField> DBNSFP_CONSERVATION_FIELDS = {
    {"phyloP100way_vertebrate", "phyloP100way_vertebrate", "PhyloP 100-way vertebrate score", true, true, 1.6},
    {"phyloP470way_mammalian", "phyloP470way_mammalian", "PhyloP 470-way mammalian score", true, true, 1.6},
    {"phyloP17way_primate", "phyloP17way_primate", "PhyloP 17-way primate score", true, true, 0.5},
    {"phastCons100way_vertebrate", "phastCons100way_vertebrate", "PhastCons 100-way vertebrate score (0-1)", true, true, 0.5},
    {"phastCons470way_mammalian", "phastCons470way_mammalian", "PhastCons 470-way mammalian score (0-1)", true, true, 0.5},
    {"phastCons17way_primate", "phastCons17way_primate", "PhastCons 17-way primate score (0-1)", true, true, 0.5},
    {"GERP_NR", "GERP++_NR", "GERP++ neutral rate", true, true, -1},
    {"GERP_RS", "GERP++_RS", "GERP++ rejected substitutions score", true, true, 2.0},
    {"SiPhy_29way_logOdds", "SiPhy_29way_logOdds", "SiPhy 29-way log odds score", true, true, 10.0},
    {"bStatistic", "bStatistic", "Background selection statistic", true, true, -1},
};

// Splicing Predictions (from dbNSFP)
inline const std::vector<DbNSFPField> DBNSFP_SPLICE_FIELDS = {
    {"GERP_RS_rankscore", "GERP++_RS_rankscore", "GERP++ RS rank score", true, true, -1},
    {"Interpro_domain", "Interpro_domain", "InterPro domain", false, false, -1},
    {"GTEx_V8_gene", "GTEx_V8_gene", "GTEx V8 gene expression", false, false, -1},
    {"GTEx_V8_tissue", "GTEx_V8_tissue", "GTEx V8 tissue expression", false, false, -1},
};

// Population Frequency Fields
inline const std::vector<DbNSFPField> DBNSFP_FREQUENCY_FIELDS = {
    {"gnomAD_exomes_AF", "gnomAD_exomes_AF", "gnomAD exomes allele frequency", true, false, -1},
    {"gnomAD_exomes_NFE_AF", "gnomAD_exomes_NFE_AF", "gnomAD exomes NFE allele frequency", true, false, -1},
    {"gnomAD_genomes_AF", "gnomAD_genomes_AF", "gnomAD genomes allele frequency", true, false, -1},
    {"gnomAD_genomes_NFE_AF", "gnomAD_genomes_NFE_AF", "gnomAD genomes NFE allele frequency", true, false, -1},
    {"ExAC_AF", "ExAC_AF", "ExAC allele frequency", true, false, -1},
    {"ExAC_NFE_AF", "ExAC_NFE_AF", "ExAC NFE allele frequency", true, false, -1},
    {"1000Gp3_AF", "1000Gp3_AF", "1000 Genomes phase 3 allele frequency", true, false, -1},
    {"1000Gp3_EUR_AF", "1000Gp3_EUR_AF", "1000 Genomes phase 3 EUR allele frequency", true, false, -1},
    {"ESP6500_EA_AF", "ESP6500_EA_AF", "ESP6500 European American allele frequency", true, false, -1},
    {"ESP6500_AA_AF", "ESP6500_AA_AF", "ESP6500 African American allele frequency", true, false, -1},
};

// Clinical Significance Fields
inline const std::vector<DbNSFPField> DBNSFP_CLINICAL_FIELDS = {
    {"clinvar_id", "clinvar_id", "ClinVar variation ID", false, false, -1},
    {"clinvar_clnsig", "clinvar_clnsig", "ClinVar clinical significance", false, false, -1},
    {"clinvar_trait", "clinvar_trait", "ClinVar trait/condition", false, false, -1},
    {"clinvar_review", "clinvar_review", "ClinVar review status", false, false, -1},
    {"clinvar_hgvs", "clinvar_hgvs", "ClinVar HGVS notation", false, false, -1},
    {"Uniprot_acc", "Uniprot_acc", "UniProt accession", false, false, -1},
    {"Uniprot_entry", "Uniprot_entry", "UniProt entry name", false, false, -1},
    {"HGVSc_snpEff", "HGVSc_snpEff", "HGVSc from snpEff", false, false, -1},
    {"HGVSp_snpEff", "HGVSp_snpEff", "HGVSp from snpEff", false, false, -1},
};

/**
 * Get all standard dbNSFP fields
 */
inline std::vector<DbNSFPField> get_all_dbnsfp_fields() {
    std::vector<DbNSFPField> all_fields;
    all_fields.insert(all_fields.end(), DBNSFP_PATHOGENICITY_FIELDS.begin(), DBNSFP_PATHOGENICITY_FIELDS.end());
    all_fields.insert(all_fields.end(), DBNSFP_CONSERVATION_FIELDS.begin(), DBNSFP_CONSERVATION_FIELDS.end());
    all_fields.insert(all_fields.end(), DBNSFP_SPLICE_FIELDS.begin(), DBNSFP_SPLICE_FIELDS.end());
    all_fields.insert(all_fields.end(), DBNSFP_FREQUENCY_FIELDS.begin(), DBNSFP_FREQUENCY_FIELDS.end());
    all_fields.insert(all_fields.end(), DBNSFP_CLINICAL_FIELDS.begin(), DBNSFP_CLINICAL_FIELDS.end());
    return all_fields;
}

/**
 * Get field names for a category
 */
inline std::set<std::string> get_dbnsfp_field_names(const std::vector<DbNSFPField>& fields) {
    std::set<std::string> names;
    for (const auto& f : fields) {
        names.insert(f.name);
    }
    return names;
}

/**
 * Common field presets
 */
inline std::vector<DbNSFPField> get_dbnsfp_preset(const std::string& preset) {
    if (preset == "essential") {
        // Most commonly used scores
        return {
            DBNSFP_PATHOGENICITY_FIELDS[0],  // SIFT_score
            DBNSFP_PATHOGENICITY_FIELDS[4],  // Polyphen2_HDIV_score
            DBNSFP_PATHOGENICITY_FIELDS[9],  // CADD_phred
            DBNSFP_PATHOGENICITY_FIELDS[10], // REVEL_score
            DBNSFP_PATHOGENICITY_FIELDS[11], // AlphaMissense_score
        };
    } else if (preset == "pathogenicity") {
        return DBNSFP_PATHOGENICITY_FIELDS;
    } else if (preset == "conservation") {
        return DBNSFP_CONSERVATION_FIELDS;
    } else if (preset == "frequency") {
        return DBNSFP_FREQUENCY_FIELDS;
    } else if (preset == "clinical") {
        return DBNSFP_CLINICAL_FIELDS;
    }
    return get_all_dbnsfp_fields();
}

/**
 * Parse field specification string (comma-separated)
 * Returns requested fields, or all fields if "all" is specified
 */
inline std::vector<DbNSFPField> parse_dbnsfp_fields(const std::string& field_spec) {
    std::vector<DbNSFPField> result;

    if (field_spec.empty() || field_spec == "all") {
        return get_all_dbnsfp_fields();
    }

    // Check for preset names
    if (field_spec == "essential" || field_spec == "pathogenicity" ||
        field_spec == "conservation" || field_spec == "frequency" ||
        field_spec == "clinical" || field_spec == "splicing") {
        return get_dbnsfp_preset(field_spec);
    }

    // Build lookup map
    std::map<std::string, DbNSFPField> lookup;
    for (const auto& f : get_all_dbnsfp_fields()) {
        lookup[f.name] = f;
    }

    // Parse comma-separated list
    std::istringstream iss(field_spec);
    std::string field;
    while (std::getline(iss, field, ',')) {
        // Trim whitespace
        size_t start = field.find_first_not_of(" \t");
        size_t end = field.find_last_not_of(" \t");
        if (start != std::string::npos) {
            field = field.substr(start, end - start + 1);
        }

        auto it = lookup.find(field);
        if (it != lookup.end()) {
            result.push_back(it->second);
        }
    }

    return result;
}

} // namespace vep

#endif // DBNSFP_FIELDS_HPP
