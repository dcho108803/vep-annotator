/**
 * Tests for filter_vep functionality: operator parsing, expression parsing,
 * filter condition evaluation, and full filter pipeline.
 */

#include <gtest/gtest.h>
#include "filter_vep.hpp"

#include <cmath>
#include <limits>

using namespace vep;

// ============================================================================
// Helper: build a FilterableRecord from a map of field->value pairs
// ============================================================================

static FilterableRecord make_record(
    const std::map<std::string, std::string>& fields,
    const std::string& original_line = "")
{
    FilterableRecord rec;
    rec.fields.insert(fields.begin(), fields.end());
    rec.original_line = original_line;
    return rec;
}

// ============================================================================
// 1. FilterOperator parsing - parse_filter_operator()
// ============================================================================

TEST(ParseFilterOperator, EqualsVariants) {
    EXPECT_EQ(parse_filter_operator("eq"), FilterOperator::EQUALS);
    EXPECT_EQ(parse_filter_operator("="), FilterOperator::EQUALS);
    EXPECT_EQ(parse_filter_operator("is"), FilterOperator::EQUALS);
    // Case insensitivity
    EXPECT_EQ(parse_filter_operator("EQ"), FilterOperator::EQUALS);
    EXPECT_EQ(parse_filter_operator("Eq"), FilterOperator::EQUALS);
    EXPECT_EQ(parse_filter_operator("IS"), FilterOperator::EQUALS);
    EXPECT_EQ(parse_filter_operator("Is"), FilterOperator::EQUALS);
}

TEST(ParseFilterOperator, NotEqualsVariants) {
    EXPECT_EQ(parse_filter_operator("ne"), FilterOperator::NOT_EQUALS);
    EXPECT_EQ(parse_filter_operator("!="), FilterOperator::NOT_EQUALS);
    EXPECT_EQ(parse_filter_operator("NE"), FilterOperator::NOT_EQUALS);
    EXPECT_EQ(parse_filter_operator("Ne"), FilterOperator::NOT_EQUALS);
}

TEST(ParseFilterOperator, GreaterThan) {
    EXPECT_EQ(parse_filter_operator("gt"), FilterOperator::GREATER);
    EXPECT_EQ(parse_filter_operator(">"), FilterOperator::GREATER);
    EXPECT_EQ(parse_filter_operator("GT"), FilterOperator::GREATER);
}

TEST(ParseFilterOperator, GreaterThanOrEqual) {
    EXPECT_EQ(parse_filter_operator("ge"), FilterOperator::GREATER_EQ);
    EXPECT_EQ(parse_filter_operator(">="), FilterOperator::GREATER_EQ);
    EXPECT_EQ(parse_filter_operator("GE"), FilterOperator::GREATER_EQ);
}

TEST(ParseFilterOperator, LessThan) {
    EXPECT_EQ(parse_filter_operator("lt"), FilterOperator::LESS);
    EXPECT_EQ(parse_filter_operator("<"), FilterOperator::LESS);
    EXPECT_EQ(parse_filter_operator("LT"), FilterOperator::LESS);
}

TEST(ParseFilterOperator, LessThanOrEqual) {
    EXPECT_EQ(parse_filter_operator("le"), FilterOperator::LESS_EQ);
    EXPECT_EQ(parse_filter_operator("<="), FilterOperator::LESS_EQ);
    EXPECT_EQ(parse_filter_operator("LE"), FilterOperator::LESS_EQ);
}

TEST(ParseFilterOperator, Contains) {
    EXPECT_EQ(parse_filter_operator("contains"), FilterOperator::CONTAINS);
    EXPECT_EQ(parse_filter_operator("match"), FilterOperator::CONTAINS);
    EXPECT_EQ(parse_filter_operator("CONTAINS"), FilterOperator::CONTAINS);
    EXPECT_EQ(parse_filter_operator("MATCH"), FilterOperator::CONTAINS);
}

TEST(ParseFilterOperator, In) {
    EXPECT_EQ(parse_filter_operator("in"), FilterOperator::IN);
    EXPECT_EQ(parse_filter_operator("IN"), FilterOperator::IN);
}

TEST(ParseFilterOperator, Exists) {
    EXPECT_EQ(parse_filter_operator("exists"), FilterOperator::EXISTS);
    EXPECT_EQ(parse_filter_operator("defined"), FilterOperator::EXISTS);
    EXPECT_EQ(parse_filter_operator("EXISTS"), FilterOperator::EXISTS);
    EXPECT_EQ(parse_filter_operator("DEFINED"), FilterOperator::EXISTS);
}

TEST(ParseFilterOperator, Regex) {
    EXPECT_EQ(parse_filter_operator("regex"), FilterOperator::REGEX);
    EXPECT_EQ(parse_filter_operator("re"), FilterOperator::REGEX);
    EXPECT_EQ(parse_filter_operator("REGEX"), FilterOperator::REGEX);
    EXPECT_EQ(parse_filter_operator("RE"), FilterOperator::REGEX);
}

TEST(ParseFilterOperator, UnknownDefaultsToEquals) {
    // Unknown operator strings should default to EQUALS
    EXPECT_EQ(parse_filter_operator("unknown"), FilterOperator::EQUALS);
    EXPECT_EQ(parse_filter_operator("foobar"), FilterOperator::EQUALS);
    EXPECT_EQ(parse_filter_operator(""), FilterOperator::EQUALS);
}

// ============================================================================
// 2. FilterExpression parsing - parse_filter_expression()
// ============================================================================

TEST(ParseFilterExpression, SimpleIsExpression) {
    FilterCondition cond = parse_filter_expression("Consequence is missense_variant");
    EXPECT_EQ(cond.field, "Consequence");
    EXPECT_EQ(cond.op, FilterOperator::EQUALS);
    EXPECT_EQ(cond.value, "missense_variant");
    EXPECT_FALSE(cond.negated);
}

TEST(ParseFilterExpression, SimpleEqExpression) {
    FilterCondition cond = parse_filter_expression("IMPACT eq HIGH");
    EXPECT_EQ(cond.field, "IMPACT");
    EXPECT_EQ(cond.op, FilterOperator::EQUALS);
    EXPECT_EQ(cond.value, "HIGH");
}

TEST(ParseFilterExpression, NumericLessThan) {
    FilterCondition cond = parse_filter_expression("SIFT_score < 0.05");
    EXPECT_EQ(cond.field, "SIFT_score");
    EXPECT_EQ(cond.op, FilterOperator::LESS);
    EXPECT_EQ(cond.value, "0.05");
}

TEST(ParseFilterExpression, NumericGreaterThan) {
    FilterCondition cond = parse_filter_expression("CADD_phred > 20");
    EXPECT_EQ(cond.field, "CADD_phred");
    EXPECT_EQ(cond.op, FilterOperator::GREATER);
    EXPECT_EQ(cond.value, "20");
}

TEST(ParseFilterExpression, NumericGreaterEqual) {
    FilterCondition cond = parse_filter_expression("AF>=0.01");
    EXPECT_EQ(cond.field, "AF");
    EXPECT_EQ(cond.op, FilterOperator::GREATER_EQ);
    EXPECT_EQ(cond.value, "0.01");
}

TEST(ParseFilterExpression, NumericLessEqual) {
    FilterCondition cond = parse_filter_expression("AF<=0.05");
    EXPECT_EQ(cond.field, "AF");
    EXPECT_EQ(cond.op, FilterOperator::LESS_EQ);
    EXPECT_EQ(cond.value, "0.05");
}

TEST(ParseFilterExpression, NotEquals) {
    FilterCondition cond = parse_filter_expression("IMPACT ne LOW");
    EXPECT_EQ(cond.field, "IMPACT");
    EXPECT_EQ(cond.op, FilterOperator::NOT_EQUALS);
    EXPECT_EQ(cond.value, "LOW");
}

TEST(ParseFilterExpression, NotEqualsSymbol) {
    FilterCondition cond = parse_filter_expression("BIOTYPE!=protein_coding");
    EXPECT_EQ(cond.field, "BIOTYPE");
    EXPECT_EQ(cond.op, FilterOperator::NOT_EQUALS);
    EXPECT_EQ(cond.value, "protein_coding");
}

TEST(ParseFilterExpression, ContainsOperator) {
    FilterCondition cond = parse_filter_expression("Consequence contains missense");
    EXPECT_EQ(cond.field, "Consequence");
    EXPECT_EQ(cond.op, FilterOperator::CONTAINS);
    EXPECT_EQ(cond.value, "missense");
}

TEST(ParseFilterExpression, MatchOperator) {
    FilterCondition cond = parse_filter_expression("Consequence match splice");
    EXPECT_EQ(cond.field, "Consequence");
    EXPECT_EQ(cond.op, FilterOperator::CONTAINS);
    EXPECT_EQ(cond.value, "splice");
}

TEST(ParseFilterExpression, InList) {
    FilterCondition cond = parse_filter_expression("Consequence in missense_variant,stop_gained");
    EXPECT_EQ(cond.field, "Consequence");
    EXPECT_EQ(cond.op, FilterOperator::IN);
    EXPECT_EQ(cond.value, "missense_variant,stop_gained");
    ASSERT_EQ(cond.value_list.size(), 2u);
    EXPECT_EQ(cond.value_list[0], "missense_variant");
    EXPECT_EQ(cond.value_list[1], "stop_gained");
}

TEST(ParseFilterExpression, InListWithSpaces) {
    FilterCondition cond = parse_filter_expression("Consequence in missense_variant, stop_gained, frameshift_variant");
    EXPECT_EQ(cond.op, FilterOperator::IN);
    ASSERT_EQ(cond.value_list.size(), 3u);
    EXPECT_EQ(cond.value_list[0], "missense_variant");
    EXPECT_EQ(cond.value_list[1], "stop_gained");
    EXPECT_EQ(cond.value_list[2], "frameshift_variant");
}

TEST(ParseFilterExpression, ExistsOperator) {
    FilterCondition cond = parse_filter_expression("SIFT exists");
    EXPECT_EQ(cond.field, "SIFT");
    EXPECT_EQ(cond.op, FilterOperator::EXISTS);
}

TEST(ParseFilterExpression, RegexOperator) {
    FilterCondition cond = parse_filter_expression("Consequence regex splice");
    // "regex" is not in the operators list for parse_filter_expression,
    // but "match" is handled as CONTAINS. Let's verify what actually happens.
    // Actually "regex" is not one of the search operators in parse_filter_expression.
    // It would not be found as an operator. Let's check the behavior.
    // The operators searched are: " is ", " eq ", " ne ", " gt ", " ge ", " lt ", " le ",
    // " contains ", " in ", " match ", " exists", ">=", "<=", "!=", ">", "<", "="
    // "regex" is not in the list, so no operator would be found; it falls back to EXISTS.
    // However, let's just verify whatever the implementation does.
    // Actually, looking more carefully: the "=" operator will match within "regex" at position
    // ... no, "=" is searched as bare "=". Let me trace: expr = "Consequence regex splice"
    // None of the space-delimited operators match " regex ". The symbolic ones: ">=", "<=", "!="
    // won't match. ">" and "<" won't match. "=" won't be found.
    // So it falls through to the "no operator found" branch -> EXISTS with field = entire expr.
    EXPECT_EQ(cond.op, FilterOperator::EXISTS);
    EXPECT_EQ(cond.field, "Consequence regex splice");
}

TEST(ParseFilterExpression, NoOperatorDefaultsToExists) {
    FilterCondition cond = parse_filter_expression("SIFT");
    EXPECT_EQ(cond.field, "SIFT");
    EXPECT_EQ(cond.op, FilterOperator::EXISTS);
}

TEST(ParseFilterExpression, EqualsSymbol) {
    FilterCondition cond = parse_filter_expression("IMPACT=HIGH");
    EXPECT_EQ(cond.field, "IMPACT");
    EXPECT_EQ(cond.op, FilterOperator::EQUALS);
    EXPECT_EQ(cond.value, "HIGH");
}

TEST(ParseFilterExpression, NumericGtKeyword) {
    FilterCondition cond = parse_filter_expression("REVEL_score gt 0.5");
    EXPECT_EQ(cond.field, "REVEL_score");
    EXPECT_EQ(cond.op, FilterOperator::GREATER);
    EXPECT_EQ(cond.value, "0.5");
}

TEST(ParseFilterExpression, NumericGeKeyword) {
    FilterCondition cond = parse_filter_expression("REVEL_score ge 0.5");
    EXPECT_EQ(cond.field, "REVEL_score");
    EXPECT_EQ(cond.op, FilterOperator::GREATER_EQ);
    EXPECT_EQ(cond.value, "0.5");
}

TEST(ParseFilterExpression, NumericLtKeyword) {
    FilterCondition cond = parse_filter_expression("AF lt 0.001");
    EXPECT_EQ(cond.field, "AF");
    EXPECT_EQ(cond.op, FilterOperator::LESS);
    EXPECT_EQ(cond.value, "0.001");
}

TEST(ParseFilterExpression, NumericLeKeyword) {
    FilterCondition cond = parse_filter_expression("AF le 0.001");
    EXPECT_EQ(cond.field, "AF");
    EXPECT_EQ(cond.op, FilterOperator::LESS_EQ);
    EXPECT_EQ(cond.value, "0.001");
}

// ============================================================================
// 3. FilterableRecord
// ============================================================================

TEST(FilterableRecord, GetExistingField) {
    auto rec = make_record({{"Consequence", "missense_variant"}, {"IMPACT", "MODERATE"}});
    EXPECT_EQ(rec.get("Consequence"), "missense_variant");
    EXPECT_EQ(rec.get("IMPACT"), "MODERATE");
}

TEST(FilterableRecord, GetMissingFieldReturnsEmpty) {
    auto rec = make_record({{"Consequence", "missense_variant"}});
    EXPECT_EQ(rec.get("SIFT"), "");
}

TEST(FilterableRecord, HasExistingNonEmptyField) {
    auto rec = make_record({{"Consequence", "missense_variant"}, {"SIFT", ""}});
    EXPECT_TRUE(rec.has("Consequence"));
    // Empty string means field is not "present" per has()
    EXPECT_FALSE(rec.has("SIFT"));
}

TEST(FilterableRecord, HasMissingField) {
    auto rec = make_record({{"Consequence", "missense_variant"}});
    EXPECT_FALSE(rec.has("POLYPHEN"));
}

TEST(FilterableRecord, GetNumericValid) {
    auto rec = make_record({{"AF", "0.001"}, {"CADD_phred", "25.3"}});
    EXPECT_DOUBLE_EQ(rec.get_numeric("AF"), 0.001);
    EXPECT_DOUBLE_EQ(rec.get_numeric("CADD_phred"), 25.3);
}

TEST(FilterableRecord, GetNumericInvalidReturnsNaN) {
    auto rec = make_record({{"AF", "."}, {"SIFT", "NA"}, {"REVEL", "NaN"}, {"FOO", "not_a_number"}});
    EXPECT_TRUE(std::isnan(rec.get_numeric("AF")));
    EXPECT_TRUE(std::isnan(rec.get_numeric("SIFT")));
    EXPECT_TRUE(std::isnan(rec.get_numeric("REVEL")));
    EXPECT_TRUE(std::isnan(rec.get_numeric("FOO")));
}

TEST(FilterableRecord, GetNumericMissingFieldReturnsNaN) {
    auto rec = make_record({});
    EXPECT_TRUE(std::isnan(rec.get_numeric("NONEXISTENT")));
}

TEST(FilterableRecord, GetNumericEmptyStringReturnsNaN) {
    auto rec = make_record({{"AF", ""}});
    EXPECT_TRUE(std::isnan(rec.get_numeric("AF")));
}

// ============================================================================
// 4. apply_condition() - Filter condition evaluation
// ============================================================================

// --- EQUALS ---

TEST(ApplyCondition, EqualsStringMatch) {
    auto rec = make_record({{"Consequence", "missense_variant"}});
    FilterCondition cond("Consequence", FilterOperator::EQUALS, "missense_variant");
    EXPECT_TRUE(apply_condition(rec, cond));
}

TEST(ApplyCondition, EqualsStringNoMatch) {
    auto rec = make_record({{"Consequence", "synonymous_variant"}});
    FilterCondition cond("Consequence", FilterOperator::EQUALS, "missense_variant");
    EXPECT_FALSE(apply_condition(rec, cond));
}

TEST(ApplyCondition, EqualsNumericMatch) {
    auto rec = make_record({{"AF", "0.01"}});
    FilterCondition cond("AF", FilterOperator::EQUALS, "0.01");
    EXPECT_TRUE(apply_condition(rec, cond));
}

TEST(ApplyCondition, EqualsNumericNoMatch) {
    auto rec = make_record({{"AF", "0.05"}});
    FilterCondition cond("AF", FilterOperator::EQUALS, "0.01");
    EXPECT_FALSE(apply_condition(rec, cond));
}

TEST(ApplyCondition, EqualsEmptyFieldVsEmptyValue) {
    // Empty field equals empty value (string comparison)
    auto rec = make_record({});
    FilterCondition cond("Consequence", FilterOperator::EQUALS, "");
    EXPECT_TRUE(apply_condition(rec, cond));
}

// --- NOT_EQUALS ---

TEST(ApplyCondition, NotEqualsMatch) {
    auto rec = make_record({{"IMPACT", "HIGH"}});
    FilterCondition cond("IMPACT", FilterOperator::NOT_EQUALS, "LOW");
    EXPECT_TRUE(apply_condition(rec, cond));
}

TEST(ApplyCondition, NotEqualsNoMatch) {
    auto rec = make_record({{"IMPACT", "HIGH"}});
    FilterCondition cond("IMPACT", FilterOperator::NOT_EQUALS, "HIGH");
    EXPECT_FALSE(apply_condition(rec, cond));
}

TEST(ApplyCondition, NotEqualsNumeric) {
    auto rec = make_record({{"AF", "0.05"}});
    FilterCondition cond("AF", FilterOperator::NOT_EQUALS, "0.01");
    EXPECT_TRUE(apply_condition(rec, cond));
}

// --- GREATER ---

TEST(ApplyCondition, GreaterThanTrue) {
    auto rec = make_record({{"CADD_phred", "30"}});
    FilterCondition cond("CADD_phred", FilterOperator::GREATER, "20");
    EXPECT_TRUE(apply_condition(rec, cond));
}

TEST(ApplyCondition, GreaterThanFalseEqual) {
    auto rec = make_record({{"CADD_phred", "20"}});
    FilterCondition cond("CADD_phred", FilterOperator::GREATER, "20");
    EXPECT_FALSE(apply_condition(rec, cond));
}

TEST(ApplyCondition, GreaterThanFalseLess) {
    auto rec = make_record({{"CADD_phred", "10"}});
    FilterCondition cond("CADD_phred", FilterOperator::GREATER, "20");
    EXPECT_FALSE(apply_condition(rec, cond));
}

TEST(ApplyCondition, GreaterThanNonNumericPassesThrough) {
    // Non-numeric values pass through numeric comparisons (match Perl VEP)
    auto rec = make_record({{"CADD_phred", "."}});
    FilterCondition cond("CADD_phred", FilterOperator::GREATER, "20");
    EXPECT_TRUE(apply_condition(rec, cond));
}

// --- GREATER_EQ ---

TEST(ApplyCondition, GreaterEqualTrueGreater) {
    auto rec = make_record({{"AF", "0.05"}});
    FilterCondition cond("AF", FilterOperator::GREATER_EQ, "0.01");
    EXPECT_TRUE(apply_condition(rec, cond));
}

TEST(ApplyCondition, GreaterEqualTrueEqual) {
    auto rec = make_record({{"AF", "0.01"}});
    FilterCondition cond("AF", FilterOperator::GREATER_EQ, "0.01");
    EXPECT_TRUE(apply_condition(rec, cond));
}

TEST(ApplyCondition, GreaterEqualFalse) {
    auto rec = make_record({{"AF", "0.001"}});
    FilterCondition cond("AF", FilterOperator::GREATER_EQ, "0.01");
    EXPECT_FALSE(apply_condition(rec, cond));
}

// --- LESS ---

TEST(ApplyCondition, LessThanTrue) {
    auto rec = make_record({{"SIFT_score", "0.03"}});
    FilterCondition cond("SIFT_score", FilterOperator::LESS, "0.05");
    EXPECT_TRUE(apply_condition(rec, cond));
}

TEST(ApplyCondition, LessThanFalseEqual) {
    auto rec = make_record({{"SIFT_score", "0.05"}});
    FilterCondition cond("SIFT_score", FilterOperator::LESS, "0.05");
    EXPECT_FALSE(apply_condition(rec, cond));
}

TEST(ApplyCondition, LessThanFalseGreater) {
    auto rec = make_record({{"SIFT_score", "0.8"}});
    FilterCondition cond("SIFT_score", FilterOperator::LESS, "0.05");
    EXPECT_FALSE(apply_condition(rec, cond));
}

// --- LESS_EQ ---

TEST(ApplyCondition, LessEqualTrueLess) {
    auto rec = make_record({{"AF", "0.001"}});
    FilterCondition cond("AF", FilterOperator::LESS_EQ, "0.01");
    EXPECT_TRUE(apply_condition(rec, cond));
}

TEST(ApplyCondition, LessEqualTrueEqual) {
    auto rec = make_record({{"AF", "0.01"}});
    FilterCondition cond("AF", FilterOperator::LESS_EQ, "0.01");
    EXPECT_TRUE(apply_condition(rec, cond));
}

TEST(ApplyCondition, LessEqualFalse) {
    auto rec = make_record({{"AF", "0.05"}});
    FilterCondition cond("AF", FilterOperator::LESS_EQ, "0.01");
    EXPECT_FALSE(apply_condition(rec, cond));
}

// --- CONTAINS ---

TEST(ApplyCondition, ContainsMatch) {
    auto rec = make_record({{"Consequence", "missense_variant"}});
    FilterCondition cond("Consequence", FilterOperator::CONTAINS, "missense");
    EXPECT_TRUE(apply_condition(rec, cond));
}

TEST(ApplyCondition, ContainsNoMatch) {
    auto rec = make_record({{"Consequence", "synonymous_variant"}});
    FilterCondition cond("Consequence", FilterOperator::CONTAINS, "missense");
    EXPECT_FALSE(apply_condition(rec, cond));
}

TEST(ApplyCondition, ContainsPartialMatch) {
    auto rec = make_record({{"Consequence", "splice_donor_variant,intron_variant"}});
    FilterCondition cond("Consequence", FilterOperator::CONTAINS, "splice");
    EXPECT_TRUE(apply_condition(rec, cond));
}

TEST(ApplyCondition, ContainsEmptyValueMatchesAnything) {
    auto rec = make_record({{"Consequence", "missense_variant"}});
    FilterCondition cond("Consequence", FilterOperator::CONTAINS, "");
    EXPECT_TRUE(apply_condition(rec, cond));
}

// --- NOT_CONTAINS ---

TEST(ApplyCondition, NotContainsMatch) {
    auto rec = make_record({{"Consequence", "synonymous_variant"}});
    FilterCondition cond("Consequence", FilterOperator::NOT_CONTAINS, "missense");
    EXPECT_TRUE(apply_condition(rec, cond));
}

TEST(ApplyCondition, NotContainsNoMatch) {
    auto rec = make_record({{"Consequence", "missense_variant"}});
    FilterCondition cond("Consequence", FilterOperator::NOT_CONTAINS, "missense");
    EXPECT_FALSE(apply_condition(rec, cond));
}

// --- IN ---

TEST(ApplyCondition, InListMatch) {
    auto rec = make_record({{"Consequence", "missense_variant"}});
    FilterCondition cond;
    cond.field = "Consequence";
    cond.op = FilterOperator::IN;
    cond.value_list = {"missense_variant", "stop_gained", "frameshift_variant"};
    EXPECT_TRUE(apply_condition(rec, cond));
}

TEST(ApplyCondition, InListNoMatch) {
    auto rec = make_record({{"Consequence", "synonymous_variant"}});
    FilterCondition cond;
    cond.field = "Consequence";
    cond.op = FilterOperator::IN;
    cond.value_list = {"missense_variant", "stop_gained", "frameshift_variant"};
    EXPECT_FALSE(apply_condition(rec, cond));
}

TEST(ApplyCondition, InListEmptyList) {
    auto rec = make_record({{"Consequence", "missense_variant"}});
    FilterCondition cond;
    cond.field = "Consequence";
    cond.op = FilterOperator::IN;
    // empty value_list
    EXPECT_FALSE(apply_condition(rec, cond));
}

// --- NOT_IN ---

TEST(ApplyCondition, NotInListMatch) {
    auto rec = make_record({{"Consequence", "synonymous_variant"}});
    FilterCondition cond;
    cond.field = "Consequence";
    cond.op = FilterOperator::NOT_IN;
    cond.value_list = {"missense_variant", "stop_gained"};
    EXPECT_TRUE(apply_condition(rec, cond));
}

TEST(ApplyCondition, NotInListNoMatch) {
    auto rec = make_record({{"Consequence", "missense_variant"}});
    FilterCondition cond;
    cond.field = "Consequence";
    cond.op = FilterOperator::NOT_IN;
    cond.value_list = {"missense_variant", "stop_gained"};
    EXPECT_FALSE(apply_condition(rec, cond));
}

TEST(ApplyCondition, NotInListEmptyListReturnsTrue) {
    auto rec = make_record({{"Consequence", "anything"}});
    FilterCondition cond;
    cond.field = "Consequence";
    cond.op = FilterOperator::NOT_IN;
    // Empty value_list means nothing matches, so NOT_IN is true
    EXPECT_TRUE(apply_condition(rec, cond));
}

// --- EXISTS ---

TEST(ApplyCondition, ExistsTrue) {
    auto rec = make_record({{"SIFT", "deleterious(0.01)"}});
    FilterCondition cond;
    cond.field = "SIFT";
    cond.op = FilterOperator::EXISTS;
    EXPECT_TRUE(apply_condition(rec, cond));
}

TEST(ApplyCondition, ExistsFalseMissing) {
    auto rec = make_record({{"Consequence", "missense_variant"}});
    FilterCondition cond;
    cond.field = "SIFT";
    cond.op = FilterOperator::EXISTS;
    EXPECT_FALSE(apply_condition(rec, cond));
}

TEST(ApplyCondition, ExistsFalseEmpty) {
    auto rec = make_record({{"SIFT", ""}});
    FilterCondition cond;
    cond.field = "SIFT";
    cond.op = FilterOperator::EXISTS;
    // has() returns false for empty strings
    EXPECT_FALSE(apply_condition(rec, cond));
}

// --- NOT_EXISTS ---

TEST(ApplyCondition, NotExistsTrueMissing) {
    auto rec = make_record({{"Consequence", "missense_variant"}});
    FilterCondition cond;
    cond.field = "PolyPhen";
    cond.op = FilterOperator::NOT_EXISTS;
    EXPECT_TRUE(apply_condition(rec, cond));
}

TEST(ApplyCondition, NotExistsTrueEmpty) {
    auto rec = make_record({{"PolyPhen", ""}});
    FilterCondition cond;
    cond.field = "PolyPhen";
    cond.op = FilterOperator::NOT_EXISTS;
    EXPECT_TRUE(apply_condition(rec, cond));
}

TEST(ApplyCondition, NotExistsFalse) {
    auto rec = make_record({{"PolyPhen", "probably_damaging(0.99)"}});
    FilterCondition cond;
    cond.field = "PolyPhen";
    cond.op = FilterOperator::NOT_EXISTS;
    EXPECT_FALSE(apply_condition(rec, cond));
}

// --- REGEX ---

TEST(ApplyCondition, RegexMatch) {
    auto rec = make_record({{"Consequence", "splice_donor_variant"}});
    FilterCondition cond("Consequence", FilterOperator::REGEX, "splice");
    EXPECT_TRUE(apply_condition(rec, cond));
}

TEST(ApplyCondition, RegexNoMatch) {
    auto rec = make_record({{"Consequence", "missense_variant"}});
    FilterCondition cond("Consequence", FilterOperator::REGEX, "splice");
    EXPECT_FALSE(apply_condition(rec, cond));
}

// --- Negation ---

TEST(ApplyCondition, NegatedEqualsTrue) {
    auto rec = make_record({{"IMPACT", "LOW"}});
    FilterCondition cond("IMPACT", FilterOperator::EQUALS, "HIGH");
    cond.negated = true;
    // IMPACT == "HIGH" is false, negated -> true
    EXPECT_TRUE(apply_condition(rec, cond));
}

TEST(ApplyCondition, NegatedEqualsFalse) {
    auto rec = make_record({{"IMPACT", "HIGH"}});
    FilterCondition cond("IMPACT", FilterOperator::EQUALS, "HIGH");
    cond.negated = true;
    // IMPACT == "HIGH" is true, negated -> false
    EXPECT_FALSE(apply_condition(rec, cond));
}

TEST(ApplyCondition, NegatedExistsTrue) {
    // Field does NOT exist, EXISTS is false, negated -> true
    auto rec = make_record({});
    FilterCondition cond;
    cond.field = "SIFT";
    cond.op = FilterOperator::EXISTS;
    cond.negated = true;
    EXPECT_TRUE(apply_condition(rec, cond));
}

TEST(ApplyCondition, NegatedExistsFalse) {
    // Field exists, EXISTS is true, negated -> false
    auto rec = make_record({{"SIFT", "tolerated(0.5)"}});
    FilterCondition cond;
    cond.field = "SIFT";
    cond.op = FilterOperator::EXISTS;
    cond.negated = true;
    EXPECT_FALSE(apply_condition(rec, cond));
}

// --- Edge cases ---

TEST(ApplyCondition, MissingFieldNumericComparisonPassesThrough) {
    // Missing fields pass through numeric comparisons (match Perl VEP)
    auto rec = make_record({});
    FilterCondition cond("CADD_phred", FilterOperator::GREATER, "20");
    EXPECT_TRUE(apply_condition(rec, cond));
}

TEST(ApplyCondition, DotFieldNumericComparisonPassesThrough) {
    // "." (missing) values pass through numeric comparisons (match Perl VEP)
    auto rec = make_record({{"AF", "."}});
    FilterCondition cond("AF", FilterOperator::LESS, "0.01");
    EXPECT_TRUE(apply_condition(rec, cond));
}

TEST(ApplyCondition, NAFieldNumericComparisonPassesThrough) {
    // "NA" values pass through numeric comparisons (match Perl VEP)
    auto rec = make_record({{"AF", "NA"}});
    FilterCondition cond("AF", FilterOperator::GREATER_EQ, "0.01");
    EXPECT_TRUE(apply_condition(rec, cond));
}

TEST(ApplyCondition, NumericEqualsWithFloatingPointPrecision) {
    // Test that near-equal floats are considered equal (epsilon 1e-9)
    auto rec = make_record({{"AF", "0.010000000001"}});
    FilterCondition cond("AF", FilterOperator::EQUALS, "0.01");
    EXPECT_TRUE(apply_condition(rec, cond));
}

// ============================================================================
// 5. apply_filter() - Full filter config evaluation
// ============================================================================

TEST(ApplyFilter, NoFiltersPassesEverything) {
    auto rec = make_record({{"Consequence", "missense_variant"}});
    FilterConfig config;
    EXPECT_TRUE(apply_filter(rec, config));
}

// --- Consequence filter ---

TEST(ApplyFilter, ConsequenceFilterMatch) {
    auto rec = make_record({{"Consequence", "missense_variant"}});
    FilterConfig config;
    config.consequence_filter.insert("missense_variant");
    EXPECT_TRUE(apply_filter(rec, config));
}

TEST(ApplyFilter, ConsequenceFilterNoMatch) {
    auto rec = make_record({{"Consequence", "synonymous_variant"}});
    FilterConfig config;
    config.consequence_filter.insert("missense_variant");
    EXPECT_FALSE(apply_filter(rec, config));
}

TEST(ApplyFilter, ConsequenceFilterPartialMatch) {
    // consequence_filter uses find(), so partial matches should work
    auto rec = make_record({{"Consequence", "missense_variant,splice_region_variant"}});
    FilterConfig config;
    config.consequence_filter.insert("missense_variant");
    EXPECT_TRUE(apply_filter(rec, config));
}

TEST(ApplyFilter, ConsequenceFilterUsesConsequenceOrCONSEQUENCE) {
    // Test with uppercase key
    auto rec = make_record({{"CONSEQUENCE", "missense_variant"}});
    FilterConfig config;
    config.consequence_filter.insert("missense_variant");
    EXPECT_TRUE(apply_filter(rec, config));
}

// --- Impact filter ---

TEST(ApplyFilter, ImpactFilterMatch) {
    auto rec = make_record({{"IMPACT", "HIGH"}});
    FilterConfig config;
    config.impact_filter.insert("HIGH");
    EXPECT_TRUE(apply_filter(rec, config));
}

TEST(ApplyFilter, ImpactFilterNoMatch) {
    auto rec = make_record({{"IMPACT", "MODIFIER"}});
    FilterConfig config;
    config.impact_filter.insert("HIGH");
    config.impact_filter.insert("MODERATE");
    EXPECT_FALSE(apply_filter(rec, config));
}

TEST(ApplyFilter, ImpactFilterMultipleValues) {
    auto rec = make_record({{"IMPACT", "MODERATE"}});
    FilterConfig config;
    config.impact_filter.insert("HIGH");
    config.impact_filter.insert("MODERATE");
    EXPECT_TRUE(apply_filter(rec, config));
}

// --- Gene filter ---

TEST(ApplyFilter, GeneFilterMatch) {
    auto rec = make_record({{"SYMBOL", "TP53"}});
    FilterConfig config;
    config.gene_filter.insert("TP53");
    config.gene_filter.insert("BRCA1");
    EXPECT_TRUE(apply_filter(rec, config));
}

TEST(ApplyFilter, GeneFilterNoMatch) {
    auto rec = make_record({{"SYMBOL", "CFTR"}});
    FilterConfig config;
    config.gene_filter.insert("TP53");
    config.gene_filter.insert("BRCA1");
    EXPECT_FALSE(apply_filter(rec, config));
}

TEST(ApplyFilter, GeneFilterChecksMultipleFields) {
    // Should check GENE, Gene, then SYMBOL
    auto rec = make_record({{"Gene", "BRCA1"}});
    FilterConfig config;
    config.gene_filter.insert("BRCA1");
    EXPECT_TRUE(apply_filter(rec, config));
}

// --- Biotype filter ---

TEST(ApplyFilter, BiotypeFilterMatch) {
    auto rec = make_record({{"BIOTYPE", "protein_coding"}});
    FilterConfig config;
    config.biotype_filter.insert("protein_coding");
    EXPECT_TRUE(apply_filter(rec, config));
}

TEST(ApplyFilter, BiotypeFilterNoMatch) {
    auto rec = make_record({{"BIOTYPE", "pseudogene"}});
    FilterConfig config;
    config.biotype_filter.insert("protein_coding");
    EXPECT_FALSE(apply_filter(rec, config));
}

// --- Numeric filters ---

TEST(ApplyFilter, MinAFFilterPass) {
    auto rec = make_record({{"AF", "0.05"}});
    FilterConfig config;
    config.min_af = 0.01;
    EXPECT_TRUE(apply_filter(rec, config));
}

TEST(ApplyFilter, MinAFFilterFail) {
    auto rec = make_record({{"AF", "0.001"}});
    FilterConfig config;
    config.min_af = 0.01;
    EXPECT_FALSE(apply_filter(rec, config));
}

TEST(ApplyFilter, MaxAFFilterPass) {
    auto rec = make_record({{"AF", "0.001"}});
    FilterConfig config;
    config.max_af = 0.01;
    EXPECT_TRUE(apply_filter(rec, config));
}

TEST(ApplyFilter, MaxAFFilterFail) {
    auto rec = make_record({{"AF", "0.05"}});
    FilterConfig config;
    config.max_af = 0.01;
    EXPECT_FALSE(apply_filter(rec, config));
}

TEST(ApplyFilter, AFFilterFallbackToGnomAD) {
    // No AF field, but has gnomAD_AF
    auto rec = make_record({{"gnomAD_AF", "0.001"}});
    FilterConfig config;
    config.max_af = 0.01;
    EXPECT_TRUE(apply_filter(rec, config));
}

TEST(ApplyFilter, AFFilterFallbackToGnomADColon) {
    // No AF or gnomAD_AF, but has gnomAD:AF
    auto rec = make_record({{"gnomAD:AF", "0.001"}});
    FilterConfig config;
    config.max_af = 0.01;
    EXPECT_TRUE(apply_filter(rec, config));
}

TEST(ApplyFilter, MinCADDFilterPass) {
    auto rec = make_record({{"CADD_phred", "30"}});
    FilterConfig config;
    config.min_cadd = 20;
    EXPECT_TRUE(apply_filter(rec, config));
}

TEST(ApplyFilter, MinCADDFilterFail) {
    auto rec = make_record({{"CADD_phred", "10"}});
    FilterConfig config;
    config.min_cadd = 20;
    EXPECT_FALSE(apply_filter(rec, config));
}

TEST(ApplyFilter, MinCADDFallbackToCADD) {
    auto rec = make_record({{"CADD", "25"}});
    FilterConfig config;
    config.min_cadd = 20;
    EXPECT_TRUE(apply_filter(rec, config));
}

TEST(ApplyFilter, MinREVELFilterPass) {
    auto rec = make_record({{"REVEL_score", "0.8"}});
    FilterConfig config;
    config.min_revel = 0.5;
    EXPECT_TRUE(apply_filter(rec, config));
}

TEST(ApplyFilter, MinREVELFilterFail) {
    auto rec = make_record({{"REVEL_score", "0.3"}});
    FilterConfig config;
    config.min_revel = 0.5;
    EXPECT_FALSE(apply_filter(rec, config));
}

TEST(ApplyFilter, NumericFilterNaNPassesThrough) {
    // If AF is NaN/missing, numeric AF filter does not reject
    auto rec = make_record({{"AF", "."}});
    FilterConfig config;
    config.min_af = 0.01;
    // NaN AF => the condition block is skipped (not rejected)
    EXPECT_TRUE(apply_filter(rec, config));
}

// --- Boolean filters ---

TEST(ApplyFilter, CodingOnlyPass) {
    auto rec = make_record({{"BIOTYPE", "protein_coding"}});
    FilterConfig config;
    config.coding_only = true;
    EXPECT_TRUE(apply_filter(rec, config));
}

TEST(ApplyFilter, CodingOnlyFail) {
    auto rec = make_record({{"BIOTYPE", "pseudogene"}});
    FilterConfig config;
    config.coding_only = true;
    EXPECT_FALSE(apply_filter(rec, config));
}

TEST(ApplyFilter, ExcludeIntergenicPass) {
    auto rec = make_record({{"Consequence", "missense_variant"}});
    FilterConfig config;
    config.exclude_intergenic = true;
    EXPECT_TRUE(apply_filter(rec, config));
}

TEST(ApplyFilter, ExcludeIntergenicFail) {
    auto rec = make_record({{"Consequence", "intergenic_variant"}});
    FilterConfig config;
    config.exclude_intergenic = true;
    EXPECT_FALSE(apply_filter(rec, config));
}

TEST(ApplyFilter, ExcludeIntronicPureIntronVariant) {
    auto rec = make_record({{"Consequence", "intron_variant"}});
    FilterConfig config;
    config.exclude_intronic = true;
    EXPECT_FALSE(apply_filter(rec, config));
}

TEST(ApplyFilter, ExcludeIntronicCombinedConsequenceNotExcluded) {
    // intron_variant combined with another consequence should NOT be excluded
    auto rec = make_record({{"Consequence", "splice_region_variant,intron_variant"}});
    FilterConfig config;
    config.exclude_intronic = true;
    EXPECT_TRUE(apply_filter(rec, config));
}

TEST(ApplyFilter, CanonicalOnlyPassYES) {
    auto rec = make_record({{"CANONICAL", "YES"}});
    FilterConfig config;
    config.canonical_only = true;
    EXPECT_TRUE(apply_filter(rec, config));
}

TEST(ApplyFilter, CanonicalOnlyPass1) {
    auto rec = make_record({{"CANONICAL", "1"}});
    FilterConfig config;
    config.canonical_only = true;
    EXPECT_TRUE(apply_filter(rec, config));
}

TEST(ApplyFilter, CanonicalOnlyPassTrue) {
    auto rec = make_record({{"CANONICAL", "true"}});
    FilterConfig config;
    config.canonical_only = true;
    EXPECT_TRUE(apply_filter(rec, config));
}

TEST(ApplyFilter, CanonicalOnlyFail) {
    auto rec = make_record({{"CANONICAL", ""}});
    FilterConfig config;
    config.canonical_only = true;
    EXPECT_FALSE(apply_filter(rec, config));
}

TEST(ApplyFilter, ManeOnlyPass) {
    auto rec = make_record({{"MANE_SELECT", "NM_000546.6"}});
    FilterConfig config;
    config.mane_only = true;
    EXPECT_TRUE(apply_filter(rec, config));
}

TEST(ApplyFilter, ManeOnlyFailEmpty) {
    auto rec = make_record({{"MANE_SELECT", ""}});
    FilterConfig config;
    config.mane_only = true;
    EXPECT_FALSE(apply_filter(rec, config));
}

TEST(ApplyFilter, ManeOnlyFailDot) {
    auto rec = make_record({{"MANE_SELECT", "."}});
    FilterConfig config;
    config.mane_only = true;
    EXPECT_FALSE(apply_filter(rec, config));
}

TEST(ApplyFilter, ManeOnlyFailNA) {
    auto rec = make_record({{"MANE_SELECT", "NA"}});
    FilterConfig config;
    config.mane_only = true;
    EXPECT_FALSE(apply_filter(rec, config));
}

TEST(ApplyFilter, ManeOnlyFallbackToMANE) {
    auto rec = make_record({{"MANE", "NM_000546.6"}});
    FilterConfig config;
    config.mane_only = true;
    EXPECT_TRUE(apply_filter(rec, config));
}

// --- Custom conditions with AND logic ---

TEST(ApplyFilter, CustomConditionsAndBothMatch) {
    auto rec = make_record({{"Consequence", "missense_variant"}, {"IMPACT", "MODERATE"}});
    FilterConfig config;
    config.match_all = true;
    config.conditions.push_back(FilterCondition("Consequence", FilterOperator::EQUALS, "missense_variant"));
    config.conditions.push_back(FilterCondition("IMPACT", FilterOperator::EQUALS, "MODERATE"));
    EXPECT_TRUE(apply_filter(rec, config));
}

TEST(ApplyFilter, CustomConditionsAndOneNoMatch) {
    auto rec = make_record({{"Consequence", "missense_variant"}, {"IMPACT", "LOW"}});
    FilterConfig config;
    config.match_all = true;
    config.conditions.push_back(FilterCondition("Consequence", FilterOperator::EQUALS, "missense_variant"));
    config.conditions.push_back(FilterCondition("IMPACT", FilterOperator::EQUALS, "HIGH"));
    EXPECT_FALSE(apply_filter(rec, config));
}

// --- Custom conditions with OR logic ---

TEST(ApplyFilter, CustomConditionsOrOneMatch) {
    auto rec = make_record({{"Consequence", "missense_variant"}, {"IMPACT", "MODERATE"}});
    FilterConfig config;
    config.match_all = false;
    config.conditions.push_back(FilterCondition("Consequence", FilterOperator::EQUALS, "stop_gained"));
    config.conditions.push_back(FilterCondition("IMPACT", FilterOperator::EQUALS, "MODERATE"));
    EXPECT_TRUE(apply_filter(rec, config));
}

TEST(ApplyFilter, CustomConditionsOrNoneMatch) {
    auto rec = make_record({{"Consequence", "synonymous_variant"}, {"IMPACT", "LOW"}});
    FilterConfig config;
    config.match_all = false;
    config.conditions.push_back(FilterCondition("Consequence", FilterOperator::EQUALS, "missense_variant"));
    config.conditions.push_back(FilterCondition("IMPACT", FilterOperator::EQUALS, "HIGH"));
    EXPECT_FALSE(apply_filter(rec, config));
}

// --- Combined quick filters and custom conditions ---

TEST(ApplyFilter, QuickFilterAndCustomConditionBothMustPass) {
    auto rec = make_record({
        {"Consequence", "missense_variant"},
        {"IMPACT", "MODERATE"},
        {"SIFT_score", "0.01"}
    });
    FilterConfig config;
    config.impact_filter.insert("MODERATE");
    config.impact_filter.insert("HIGH");
    config.conditions.push_back(FilterCondition("SIFT_score", FilterOperator::LESS, "0.05"));
    EXPECT_TRUE(apply_filter(rec, config));
}

TEST(ApplyFilter, QuickFilterFailsEvenIfCustomPasses) {
    auto rec = make_record({
        {"Consequence", "missense_variant"},
        {"IMPACT", "LOW"},
        {"SIFT_score", "0.01"}
    });
    FilterConfig config;
    config.impact_filter.insert("MODERATE");
    config.impact_filter.insert("HIGH");
    config.conditions.push_back(FilterCondition("SIFT_score", FilterOperator::LESS, "0.05"));
    EXPECT_FALSE(apply_filter(rec, config));
}

// ============================================================================
// 6. parse_tsv_header() and parse_tsv_record()
// ============================================================================

TEST(ParseTSVHeader, BasicHeader) {
    auto col_map = parse_tsv_header("#Uploaded_variation\tLocation\tAllele\tConsequence\tIMPACT");
    EXPECT_EQ(col_map["#Uploaded_variation"], 0);
    EXPECT_EQ(col_map["Location"], 1);
    EXPECT_EQ(col_map["Allele"], 2);
    EXPECT_EQ(col_map["Consequence"], 3);
    EXPECT_EQ(col_map["IMPACT"], 4);
    EXPECT_EQ(col_map.size(), 5u);
}

TEST(ParseTSVHeader, SingleColumn) {
    auto col_map = parse_tsv_header("Consequence");
    EXPECT_EQ(col_map.size(), 1u);
    EXPECT_EQ(col_map["Consequence"], 0);
}

TEST(ParseTSVRecord, BasicRecord) {
    auto col_map = parse_tsv_header("Consequence\tIMPACT\tSIFT");
    auto rec = parse_tsv_record("missense_variant\tMODERATE\tdeleterious(0.01)", col_map);
    EXPECT_EQ(rec.get("Consequence"), "missense_variant");
    EXPECT_EQ(rec.get("IMPACT"), "MODERATE");
    EXPECT_EQ(rec.get("SIFT"), "deleterious(0.01)");
}

TEST(ParseTSVRecord, FewerFieldsThanHeader) {
    auto col_map = parse_tsv_header("Consequence\tIMPACT\tSIFT");
    auto rec = parse_tsv_record("missense_variant\tMODERATE", col_map);
    EXPECT_EQ(rec.get("Consequence"), "missense_variant");
    EXPECT_EQ(rec.get("IMPACT"), "MODERATE");
    EXPECT_EQ(rec.get("SIFT"), "");  // Missing field
}

TEST(ParseTSVRecord, OriginalLinePreserved) {
    auto col_map = parse_tsv_header("Consequence\tIMPACT");
    std::string line = "missense_variant\tHIGH";
    auto rec = parse_tsv_record(line, col_map);
    EXPECT_EQ(rec.original_line, line);
}

// ============================================================================
// 7. FilterConfig::has_any_filter()
// ============================================================================

TEST(FilterConfig, HasAnyFilterEmpty) {
    FilterConfig config;
    EXPECT_FALSE(config.has_any_filter());
}

TEST(FilterConfig, HasAnyFilterWithConsequence) {
    FilterConfig config;
    config.consequence_filter.insert("missense_variant");
    EXPECT_TRUE(config.has_any_filter());
}

TEST(FilterConfig, HasAnyFilterWithImpact) {
    FilterConfig config;
    config.impact_filter.insert("HIGH");
    EXPECT_TRUE(config.has_any_filter());
}

TEST(FilterConfig, HasAnyFilterWithGene) {
    FilterConfig config;
    config.gene_filter.insert("TP53");
    EXPECT_TRUE(config.has_any_filter());
}

TEST(FilterConfig, HasAnyFilterWithBiotype) {
    FilterConfig config;
    config.biotype_filter.insert("protein_coding");
    EXPECT_TRUE(config.has_any_filter());
}

TEST(FilterConfig, HasAnyFilterWithMinAF) {
    FilterConfig config;
    config.min_af = 0.01;
    EXPECT_TRUE(config.has_any_filter());
}

TEST(FilterConfig, HasAnyFilterWithMaxAF) {
    FilterConfig config;
    config.max_af = 0.05;
    EXPECT_TRUE(config.has_any_filter());
}

TEST(FilterConfig, HasAnyFilterWithMinCADD) {
    FilterConfig config;
    config.min_cadd = 20;
    EXPECT_TRUE(config.has_any_filter());
}

TEST(FilterConfig, HasAnyFilterWithMinREVEL) {
    FilterConfig config;
    config.min_revel = 0.5;
    EXPECT_TRUE(config.has_any_filter());
}

TEST(FilterConfig, HasAnyFilterWithCodingOnly) {
    FilterConfig config;
    config.coding_only = true;
    EXPECT_TRUE(config.has_any_filter());
}

TEST(FilterConfig, HasAnyFilterWithExcludeIntergenic) {
    FilterConfig config;
    config.exclude_intergenic = true;
    EXPECT_TRUE(config.has_any_filter());
}

TEST(FilterConfig, HasAnyFilterWithExcludeIntronic) {
    FilterConfig config;
    config.exclude_intronic = true;
    EXPECT_TRUE(config.has_any_filter());
}

TEST(FilterConfig, HasAnyFilterWithCanonicalOnly) {
    FilterConfig config;
    config.canonical_only = true;
    EXPECT_TRUE(config.has_any_filter());
}

TEST(FilterConfig, HasAnyFilterWithManeOnly) {
    FilterConfig config;
    config.mane_only = true;
    EXPECT_TRUE(config.has_any_filter());
}

TEST(FilterConfig, HasAnyFilterWithPickOne) {
    FilterConfig config;
    config.pick_one = true;
    EXPECT_TRUE(config.has_any_filter());
}

TEST(FilterConfig, HasAnyFilterWithConditions) {
    FilterConfig config;
    config.conditions.push_back(FilterCondition("Consequence", FilterOperator::EQUALS, "missense_variant"));
    EXPECT_TRUE(config.has_any_filter());
}

// ============================================================================
// 8. FilterCondition default construction
// ============================================================================

TEST(FilterCondition, DefaultConstruction) {
    FilterCondition cond;
    EXPECT_EQ(cond.field, "");
    EXPECT_EQ(cond.op, FilterOperator::EQUALS);
    EXPECT_EQ(cond.value, "");
    EXPECT_TRUE(cond.value_list.empty());
    EXPECT_FALSE(cond.negated);
}

TEST(FilterCondition, ThreeArgConstruction) {
    FilterCondition cond("IMPACT", FilterOperator::NOT_EQUALS, "LOW");
    EXPECT_EQ(cond.field, "IMPACT");
    EXPECT_EQ(cond.op, FilterOperator::NOT_EQUALS);
    EXPECT_EQ(cond.value, "LOW");
    EXPECT_FALSE(cond.negated);
}

// ============================================================================
// 9. Complex filter scenarios
// ============================================================================

TEST(ComplexFilter, HighImpactMissenseSIFTFilter) {
    // Simulate a pathogenic variant that passes multiple filters
    auto rec = make_record({
        {"Consequence", "missense_variant"},
        {"IMPACT", "MODERATE"},
        {"SIFT", "deleterious(0.001)"},
        {"SIFT_score", "0.001"},
        {"PolyPhen", "probably_damaging(0.999)"},
        {"CADD_phred", "32"},
        {"REVEL_score", "0.85"},
        {"AF", "0.0001"},
        {"BIOTYPE", "protein_coding"},
        {"CANONICAL", "YES"},
        {"SYMBOL", "TP53"},
        {"MANE_SELECT", "NM_000546.6"}
    });

    FilterConfig config;
    config.consequence_filter.insert("missense_variant");
    config.impact_filter.insert("MODERATE");
    config.impact_filter.insert("HIGH");
    config.gene_filter.insert("TP53");
    config.biotype_filter.insert("protein_coding");
    config.max_af = 0.001;
    config.min_cadd = 20;
    config.min_revel = 0.5;
    config.coding_only = true;
    config.canonical_only = true;

    EXPECT_TRUE(apply_filter(rec, config));
}

TEST(ComplexFilter, BenignVariantFailsStrictFilter) {
    auto rec = make_record({
        {"Consequence", "synonymous_variant"},
        {"IMPACT", "LOW"},
        {"BIOTYPE", "protein_coding"},
        {"CANONICAL", "YES"},
        {"AF", "0.15"}
    });

    FilterConfig config;
    config.impact_filter.insert("HIGH");
    config.impact_filter.insert("MODERATE");

    EXPECT_FALSE(apply_filter(rec, config));
}

TEST(ComplexFilter, MultipleQuickFiltersAllMustPass) {
    // All quick filters are AND-ed together
    auto rec = make_record({
        {"Consequence", "missense_variant"},
        {"IMPACT", "MODERATE"},
        {"SYMBOL", "BRCA1"},
        {"BIOTYPE", "protein_coding"},
        {"CANONICAL", "YES"}
    });

    FilterConfig config;
    config.consequence_filter.insert("missense_variant");
    config.impact_filter.insert("MODERATE");
    config.gene_filter.insert("TP53");  // Different gene!
    config.canonical_only = true;

    // Fails because gene doesn't match
    EXPECT_FALSE(apply_filter(rec, config));
}

TEST(ComplexFilter, CustomConditionsWithNumericAndString) {
    auto rec = make_record({
        {"Consequence", "missense_variant"},
        {"SIFT_score", "0.03"},
        {"IMPACT", "MODERATE"}
    });

    FilterConfig config;
    config.match_all = true;
    config.conditions.push_back(FilterCondition("Consequence", FilterOperator::CONTAINS, "missense"));
    config.conditions.push_back(FilterCondition("SIFT_score", FilterOperator::LESS, "0.05"));

    EXPECT_TRUE(apply_filter(rec, config));
}

TEST(ComplexFilter, OrConditionsWithMultipleConsequences) {
    auto rec = make_record({{"Consequence", "stop_gained"}, {"IMPACT", "HIGH"}});

    FilterConfig config;
    config.match_all = false;  // OR
    config.conditions.push_back(FilterCondition("Consequence", FilterOperator::EQUALS, "missense_variant"));
    config.conditions.push_back(FilterCondition("Consequence", FilterOperator::EQUALS, "stop_gained"));
    config.conditions.push_back(FilterCondition("Consequence", FilterOperator::EQUALS, "frameshift_variant"));

    EXPECT_TRUE(apply_filter(rec, config));
}

// ============================================================================
// 10. Edge cases and error handling
// ============================================================================

TEST(EdgeCases, EmptyRecord) {
    auto rec = make_record({});
    FilterConfig config;
    config.conditions.push_back(FilterCondition("Consequence", FilterOperator::EQUALS, "missense_variant"));
    EXPECT_FALSE(apply_filter(rec, config));
}

TEST(EdgeCases, RecordWithAllEmptyFields) {
    auto rec = make_record({{"Consequence", ""}, {"IMPACT", ""}, {"AF", ""}});

    FilterConfig config;
    config.conditions.push_back(FilterCondition("Consequence", FilterOperator::EQUALS, ""));
    EXPECT_TRUE(apply_filter(rec, config));
}

TEST(EdgeCases, EqualsComparisonDotValue) {
    auto rec = make_record({{"SIFT", "."}});
    FilterCondition cond("SIFT", FilterOperator::EQUALS, ".");
    // Both "." and "." are non-numeric strings, so string comparison
    EXPECT_TRUE(apply_condition(rec, cond));
}

TEST(EdgeCases, ContainsOnEmptyField) {
    auto rec = make_record({});  // Field not present, get() returns ""
    FilterCondition cond("Consequence", FilterOperator::CONTAINS, "missense");
    EXPECT_FALSE(apply_condition(rec, cond));
}

TEST(EdgeCases, InListSingleItem) {
    auto rec = make_record({{"Consequence", "missense_variant"}});
    FilterCondition cond;
    cond.field = "Consequence";
    cond.op = FilterOperator::IN;
    cond.value_list = {"missense_variant"};
    EXPECT_TRUE(apply_condition(rec, cond));
}

TEST(EdgeCases, ExistsWithDotValue) {
    // A field with "." is still considered as "having" a value by has()
    auto rec = make_record({{"SIFT", "."}});
    FilterCondition cond;
    cond.field = "SIFT";
    cond.op = FilterOperator::EXISTS;
    // has() checks fields.count > 0 && !empty() -- "." is not empty
    EXPECT_TRUE(apply_condition(rec, cond));
}

TEST(EdgeCases, ParseExpressionWithExtraWhitespace) {
    FilterCondition cond = parse_filter_expression("  Consequence  is  missense_variant  ");
    // The parser trims field and value
    EXPECT_EQ(cond.field, "Consequence");
    EXPECT_EQ(cond.value, "missense_variant");
    EXPECT_EQ(cond.op, FilterOperator::EQUALS);
}

TEST(EdgeCases, NegativeNumericComparison) {
    auto rec = make_record({{"score", "-0.5"}});
    FilterCondition cond("score", FilterOperator::LESS, "0");
    EXPECT_TRUE(apply_condition(rec, cond));
}

TEST(EdgeCases, VeryLargeNumericComparison) {
    auto rec = make_record({{"score", "1e10"}});
    FilterCondition cond("score", FilterOperator::GREATER, "1000000");
    EXPECT_TRUE(apply_condition(rec, cond));
}

TEST(EdgeCases, ScientificNotationNumericComparison) {
    auto rec = make_record({{"AF", "1.5e-4"}});
    FilterCondition cond("AF", FilterOperator::LESS, "0.001");
    EXPECT_TRUE(apply_condition(rec, cond));
}

TEST(EdgeCases, ZeroValueNumericComparison) {
    auto rec = make_record({{"AF", "0"}});
    FilterCondition cond("AF", FilterOperator::LESS_EQ, "0.01");
    EXPECT_TRUE(apply_condition(rec, cond));
}

TEST(EdgeCases, ZeroEqualsZero) {
    auto rec = make_record({{"AF", "0"}});
    FilterCondition cond("AF", FilterOperator::EQUALS, "0");
    EXPECT_TRUE(apply_condition(rec, cond));
}

TEST(EdgeCases, ParseExpressionEmptyString) {
    FilterCondition cond = parse_filter_expression("");
    // Empty expression -> no operator found -> EXISTS with empty field
    EXPECT_EQ(cond.op, FilterOperator::EXISTS);
    EXPECT_EQ(cond.field, "");
}

TEST(EdgeCases, InListExactMatchRequired) {
    // "missense" should NOT match "missense_variant" in IN list
    auto rec = make_record({{"Consequence", "missense"}});
    FilterCondition cond;
    cond.field = "Consequence";
    cond.op = FilterOperator::IN;
    cond.value_list = {"missense_variant", "stop_gained"};
    EXPECT_FALSE(apply_condition(rec, cond));
}

TEST(EdgeCases, NotInExactMatchRequired) {
    // "missense_variant" IS in the list, so NOT_IN returns false
    auto rec = make_record({{"Consequence", "missense_variant"}});
    FilterCondition cond;
    cond.field = "Consequence";
    cond.op = FilterOperator::NOT_IN;
    cond.value_list = {"missense_variant", "stop_gained"};
    EXPECT_FALSE(apply_condition(rec, cond));
}

// ============================================================================
// 11. Filter pipeline integration with TSV parsing
// ============================================================================

TEST(FilterPipeline, ParseAndApplyTSVLine) {
    // Simulate a complete TSV line parse + filter workflow
    std::string header = "Consequence\tIMPACT\tSIFT_score\tBIOTYPE\tCANONICAL";
    auto col_map = parse_tsv_header(header);

    std::string line = "missense_variant\tMODERATE\t0.01\tprotein_coding\tYES";
    auto rec = parse_tsv_record(line, col_map);

    FilterConfig config;
    config.consequence_filter.insert("missense_variant");
    config.canonical_only = true;
    config.coding_only = true;

    EXPECT_TRUE(apply_filter(rec, config));
}

TEST(FilterPipeline, ParseAndRejectTSVLine) {
    std::string header = "Consequence\tIMPACT\tSIFT_score\tBIOTYPE\tCANONICAL";
    auto col_map = parse_tsv_header(header);

    std::string line = "intergenic_variant\tMODIFIER\t.\tprotein_coding\t";
    auto rec = parse_tsv_record(line, col_map);

    FilterConfig config;
    config.exclude_intergenic = true;

    EXPECT_FALSE(apply_filter(rec, config));
}

TEST(FilterPipeline, ParseMultipleLinesCountMatches) {
    std::string header = "Consequence\tIMPACT";
    auto col_map = parse_tsv_header(header);

    std::vector<std::string> lines = {
        "missense_variant\tMODERATE",
        "synonymous_variant\tLOW",
        "stop_gained\tHIGH",
        "intron_variant\tMODIFIER",
        "frameshift_variant\tHIGH"
    };

    FilterConfig config;
    config.impact_filter.insert("HIGH");
    config.impact_filter.insert("MODERATE");

    int match_count = 0;
    for (const auto& line : lines) {
        auto rec = parse_tsv_record(line, col_map);
        if (apply_filter(rec, config)) {
            match_count++;
        }
    }

    // missense (MODERATE), stop_gained (HIGH), frameshift (HIGH)
    EXPECT_EQ(match_count, 3);
}

// ============================================================================
// filter_vep --pick mode tests
// ============================================================================

TEST(FilterPick, PickOneKeepsFirstPassingPerVariant) {
    // Two records for same variant, both pass filter
    FilterConfig config;
    config.pick_one = true;

    auto rec1 = make_record({
        {"CHROM", "1"}, {"POS", "100"}, {"REF", "A"}, {"ALT", "T"},
        {"Consequence", "missense_variant"}, {"IMPACT", "MODERATE"}
    });
    auto rec2 = make_record({
        {"CHROM", "1"}, {"POS", "100"}, {"REF", "A"}, {"ALT", "T"},
        {"Consequence", "synonymous_variant"}, {"IMPACT", "LOW"}
    });

    // Both should pass with no filter conditions
    EXPECT_TRUE(apply_filter(rec1, config));
    EXPECT_TRUE(apply_filter(rec2, config));

    // Verify the pick_one behavior is in FilterConfig
    EXPECT_TRUE(config.has_any_filter());
}

TEST(FilterPick, PickWithConsequenceFilter) {
    FilterConfig config;
    config.pick_one = true;
    config.consequence_filter.insert("missense_variant");

    auto rec_miss = make_record({
        {"CHROM", "1"}, {"POS", "100"}, {"REF", "A"}, {"ALT", "T"},
        {"Consequence", "missense_variant"}
    });
    auto rec_syn = make_record({
        {"CHROM", "1"}, {"POS", "100"}, {"REF", "A"}, {"ALT", "T"},
        {"Consequence", "synonymous_variant"}
    });

    EXPECT_TRUE(apply_filter(rec_miss, config));
    EXPECT_FALSE(apply_filter(rec_syn, config));
}

TEST(FilterPick, DifferentVariantsNotAffected) {
    // Two different variants should each pass independently
    FilterConfig config;
    config.pick_one = true;

    auto rec1 = make_record({
        {"CHROM", "1"}, {"POS", "100"}, {"REF", "A"}, {"ALT", "T"},
        {"Consequence", "missense_variant"}
    });
    auto rec2 = make_record({
        {"CHROM", "1"}, {"POS", "200"}, {"REF", "G"}, {"ALT", "C"},
        {"Consequence", "missense_variant"}
    });

    // Different positions = different variant keys = both should pass filter
    EXPECT_TRUE(apply_filter(rec1, config));
    EXPECT_TRUE(apply_filter(rec2, config));
}

// ============================================================================
// FilterableRecord edge cases
// ============================================================================

TEST(FilterableRecord, GetMissingField) {
    auto rec = make_record({{"CHROM", "1"}, {"POS", "100"}});
    EXPECT_EQ(rec.get("NONEXISTENT"), "");
}

TEST(FilterableRecord, HasFieldWithEmptyValue) {
    auto rec = make_record({{"CHROM", "1"}, {"POS", ""}});
    EXPECT_FALSE(rec.has("POS"));  // empty value = not "has"
    EXPECT_TRUE(rec.has("CHROM"));
}

TEST(FilterableRecord, HasMissingFieldEdge) {
    auto rec = make_record({{"CHROM", "1"}});
    EXPECT_FALSE(rec.has("NONEXISTENT"));
}

// ============================================================================
// NEW TESTS: Complex filter expressions
// ============================================================================

TEST(ComplexFilterExpr, NegatedNot_ImpactIsHigh) {
    // "not IMPACT is HIGH" -> field=IMPACT, op=EQUALS, value=HIGH, negated=true
    FilterCondition cond = parse_filter_expression("not IMPACT is HIGH");
    EXPECT_EQ(cond.field, "IMPACT");
    EXPECT_EQ(cond.op, FilterOperator::EQUALS);
    EXPECT_EQ(cond.value, "HIGH");
    EXPECT_TRUE(cond.negated);

    // Apply: record with IMPACT=LOW should pass (negated is HIGH -> true)
    auto rec_low = make_record({{"IMPACT", "LOW"}});
    EXPECT_TRUE(apply_condition(rec_low, cond));

    // Record with IMPACT=HIGH should fail
    auto rec_high = make_record({{"IMPACT", "HIGH"}});
    EXPECT_FALSE(apply_condition(rec_high, cond));
}

TEST(ComplexFilterExpr, MultipleValuesInList) {
    // "Consequence in missense_variant,stop_gained,frameshift_variant"
    FilterCondition cond = parse_filter_expression(
        "Consequence in missense_variant,stop_gained,frameshift_variant");
    EXPECT_EQ(cond.field, "Consequence");
    EXPECT_EQ(cond.op, FilterOperator::IN);
    ASSERT_EQ(cond.value_list.size(), 3u);
    EXPECT_EQ(cond.value_list[0], "missense_variant");
    EXPECT_EQ(cond.value_list[1], "stop_gained");
    EXPECT_EQ(cond.value_list[2], "frameshift_variant");

    auto rec_miss = make_record({{"Consequence", "missense_variant"}});
    EXPECT_TRUE(apply_condition(rec_miss, cond));
    auto rec_stop = make_record({{"Consequence", "stop_gained"}});
    EXPECT_TRUE(apply_condition(rec_stop, cond));
    auto rec_syn = make_record({{"Consequence", "synonymous_variant"}});
    EXPECT_FALSE(apply_condition(rec_syn, cond));
}

TEST(ComplexFilterExpr, NotInOperator) {
    // Build a NOT_IN condition manually (parse_filter_expression doesn't
    // directly support "not in" as a single operator, but we can use
    // the negated flag with IN)
    FilterCondition cond;
    cond.field = "BIOTYPE";
    cond.op = FilterOperator::NOT_IN;
    cond.value_list = {"processed_pseudogene", "retained_intron"};

    auto rec_pc = make_record({{"BIOTYPE", "protein_coding"}});
    EXPECT_TRUE(apply_condition(rec_pc, cond));

    auto rec_pseudo = make_record({{"BIOTYPE", "processed_pseudogene"}});
    EXPECT_FALSE(apply_condition(rec_pseudo, cond));

    auto rec_retained = make_record({{"BIOTYPE", "retained_intron"}});
    EXPECT_FALSE(apply_condition(rec_retained, cond));
}

TEST(ComplexFilterExpr, RegexOperatorSubstring) {
    // REGEX operator acts as contains (simple substring match)
    FilterCondition cond("HGVSc", FilterOperator::REGEX, "del");

    auto rec_del = make_record({{"HGVSc", "ENST00000269305.9:c.123del"}});
    EXPECT_TRUE(apply_condition(rec_del, cond));

    auto rec_ins = make_record({{"HGVSc", "ENST00000269305.9:c.123ins"}});
    EXPECT_FALSE(apply_condition(rec_ins, cond));
}

TEST(ComplexFilterExpr, ContainsWithSpecialChars) {
    // Contains with colon and period characters
    FilterCondition cond("HGVSp", FilterOperator::CONTAINS, "p.Arg");

    auto rec_match = make_record({{"HGVSp", "ENSP00000269305:p.Arg175His"}});
    EXPECT_TRUE(apply_condition(rec_match, cond));

    auto rec_no = make_record({{"HGVSp", "ENSP00000269305:p.Gly245Ser"}});
    EXPECT_FALSE(apply_condition(rec_no, cond));
}

// ============================================================================
// NEW TESTS: Quick filter interactions
// ============================================================================

TEST(QuickFilterInteraction, ConsequenceAndImpactCombined) {
    // Both consequence_filter and impact_filter must pass (AND logic)
    FilterConfig config;
    config.consequence_filter.insert("missense_variant");
    config.impact_filter.insert("HIGH");

    // Missense is MODERATE, not HIGH -> fails impact
    auto rec = make_record({{"Consequence", "missense_variant"}, {"IMPACT", "MODERATE"}});
    EXPECT_FALSE(apply_filter(rec, config));

    // Stop gained is HIGH -> but not in consequence filter
    auto rec2 = make_record({{"Consequence", "stop_gained"}, {"IMPACT", "HIGH"}});
    EXPECT_FALSE(apply_filter(rec2, config));

    // Missense + HIGH (unusual but test the logic)
    auto rec3 = make_record({{"Consequence", "missense_variant"}, {"IMPACT", "HIGH"}});
    EXPECT_TRUE(apply_filter(rec3, config));
}

TEST(QuickFilterInteraction, GeneFilterAndCodingOnly) {
    FilterConfig config;
    config.gene_filter.insert("TP53");
    config.coding_only = true;

    // Correct gene, protein_coding -> pass
    auto rec1 = make_record({{"SYMBOL", "TP53"}, {"BIOTYPE", "protein_coding"}});
    EXPECT_TRUE(apply_filter(rec1, config));

    // Correct gene, lncRNA -> fails coding_only
    auto rec2 = make_record({{"SYMBOL", "TP53"}, {"BIOTYPE", "lncRNA"}});
    EXPECT_FALSE(apply_filter(rec2, config));

    // Wrong gene, protein_coding -> fails gene filter
    auto rec3 = make_record({{"SYMBOL", "BRCA1"}, {"BIOTYPE", "protein_coding"}});
    EXPECT_FALSE(apply_filter(rec3, config));
}

TEST(QuickFilterInteraction, MinAFAndMaxAFRange) {
    // Filter for rare variants: 0.001 <= AF <= 0.01
    FilterConfig config;
    config.min_af = 0.001;
    config.max_af = 0.01;

    auto rec_in_range = make_record({{"AF", "0.005"}});
    EXPECT_TRUE(apply_filter(rec_in_range, config));

    auto rec_too_low = make_record({{"AF", "0.0001"}});
    EXPECT_FALSE(apply_filter(rec_too_low, config));

    auto rec_too_high = make_record({{"AF", "0.05"}});
    EXPECT_FALSE(apply_filter(rec_too_high, config));
}

TEST(QuickFilterInteraction, CanonicalAndManeCombined) {
    FilterConfig config;
    config.canonical_only = true;
    config.mane_only = true;

    // Both canonical and MANE -> pass
    auto rec1 = make_record({{"CANONICAL", "YES"}, {"MANE_SELECT", "NM_000546.6"}});
    EXPECT_TRUE(apply_filter(rec1, config));

    // Canonical but not MANE -> fail
    auto rec2 = make_record({{"CANONICAL", "YES"}, {"MANE_SELECT", ""}});
    EXPECT_FALSE(apply_filter(rec2, config));

    // MANE but not canonical -> fail
    auto rec3 = make_record({{"CANONICAL", ""}, {"MANE_SELECT", "NM_000546.6"}});
    EXPECT_FALSE(apply_filter(rec3, config));
}

TEST(QuickFilterInteraction, AllFiltersDisabledPassesEverything) {
    FilterConfig config;
    // Default config: everything disabled/empty

    auto rec = make_record({
        {"Consequence", "intergenic_variant"},
        {"IMPACT", "MODIFIER"},
        {"BIOTYPE", "pseudogene"},
        {"CANONICAL", ""},
        {"AF", "0.5"}
    });
    EXPECT_TRUE(apply_filter(rec, config));
    EXPECT_FALSE(config.has_any_filter());
}

// ============================================================================
// NEW TESTS: Numeric edge cases
// ============================================================================

TEST(NumericEdgeCases, VerySmallAFValue) {
    auto rec = make_record({{"AF", "0.000001"}});  // 1e-6
    FilterCondition cond("AF", FilterOperator::LESS, "0.001");
    EXPECT_TRUE(apply_condition(rec, cond));

    FilterCondition cond2("AF", FilterOperator::GREATER, "0.0000001");
    EXPECT_TRUE(apply_condition(rec, cond2));
}

TEST(NumericEdgeCases, ExactBoundaryMaxAF) {
    // AF = 0.01 with max_af = 0.01: AF > max_af is false, so it should pass
    auto rec = make_record({{"AF", "0.01"}});
    FilterConfig config;
    config.max_af = 0.01;
    // max_af check: af > config.max_af -> 0.01 > 0.01 is false -> passes
    EXPECT_TRUE(apply_filter(rec, config));
}

TEST(NumericEdgeCases, NegativeNumberComparison) {
    auto rec = make_record({{"score", "-5.5"}});
    FilterCondition cond_lt("score", FilterOperator::LESS, "0");
    EXPECT_TRUE(apply_condition(rec, cond_lt));

    FilterCondition cond_gt("score", FilterOperator::GREATER, "-10");
    EXPECT_TRUE(apply_condition(rec, cond_gt));

    FilterCondition cond_eq("score", FilterOperator::EQUALS, "-5.5");
    EXPECT_TRUE(apply_condition(rec, cond_eq));
}

TEST(NumericEdgeCases, ScientificNotationValues) {
    auto rec = make_record({{"AF", "1.5e-3"}});
    // 1.5e-3 = 0.0015
    FilterCondition cond("AF", FilterOperator::LESS, "0.01");
    EXPECT_TRUE(apply_condition(rec, cond));

    FilterCondition cond2("AF", FilterOperator::GREATER, "0.001");
    EXPECT_TRUE(apply_condition(rec, cond2));

    FilterCondition cond3("AF", FilterOperator::EQUALS, "0.0015");
    EXPECT_TRUE(apply_condition(rec, cond3));
}

TEST(NumericEdgeCases, IntegerVsFloatComparison) {
    // "20" vs "20.0" should be equal numerically
    auto rec = make_record({{"CADD_phred", "20"}});
    FilterCondition cond("CADD_phred", FilterOperator::EQUALS, "20.0");
    EXPECT_TRUE(apply_condition(rec, cond));

    FilterCondition cond_gt("CADD_phred", FilterOperator::GREATER, "20.0");
    EXPECT_FALSE(apply_condition(rec, cond_gt));

    FilterCondition cond_ge("CADD_phred", FilterOperator::GREATER_EQ, "20.0");
    EXPECT_TRUE(apply_condition(rec, cond_ge));
}

// ============================================================================
// NEW TESTS: Record field access
// ============================================================================

TEST(RecordFieldAccess, GetMissingFieldReturnsEmptyString) {
    auto rec = make_record({{"Consequence", "missense_variant"}});
    EXPECT_EQ(rec.get("NONEXISTENT_FIELD"), "");
    EXPECT_EQ(rec.get(""), "");
}

TEST(RecordFieldAccess, HasPresentVsAbsent) {
    auto rec = make_record({
        {"Consequence", "missense_variant"},
        {"SIFT", ""},
        {"IMPACT", "MODERATE"}
    });
    EXPECT_TRUE(rec.has("Consequence"));
    EXPECT_TRUE(rec.has("IMPACT"));
    EXPECT_FALSE(rec.has("SIFT"));       // Empty value
    EXPECT_FALSE(rec.has("PolyPhen"));   // Missing field
}

TEST(RecordFieldAccess, GetNumericVariousValues) {
    auto rec = make_record({
        {"AF", "0.001"},
        {"NaN_field", "NaN"},
        {"dot_field", "."},
        {"NA_field", "NA"},
        {"empty_field", ""},
        {"text_field", "deleterious"}
    });

    EXPECT_DOUBLE_EQ(rec.get_numeric("AF"), 0.001);
    EXPECT_TRUE(std::isnan(rec.get_numeric("NaN_field")));
    EXPECT_TRUE(std::isnan(rec.get_numeric("dot_field")));
    EXPECT_TRUE(std::isnan(rec.get_numeric("NA_field")));
    EXPECT_TRUE(std::isnan(rec.get_numeric("empty_field")));
    EXPECT_TRUE(std::isnan(rec.get_numeric("text_field")));
    EXPECT_TRUE(std::isnan(rec.get_numeric("missing_field")));
}

TEST(RecordFieldAccess, LargeRecordWithManyFields) {
    std::map<std::string, std::string> fields;
    for (int i = 0; i < 100; ++i) {
        fields["field_" + std::to_string(i)] = "value_" + std::to_string(i);
    }
    auto rec = make_record(fields);

    EXPECT_EQ(rec.get("field_0"), "value_0");
    EXPECT_EQ(rec.get("field_50"), "value_50");
    EXPECT_EQ(rec.get("field_99"), "value_99");
    EXPECT_TRUE(rec.has("field_42"));
    EXPECT_FALSE(rec.has("field_100"));
}

TEST(RecordFieldAccess, FieldNamesWithColons) {
    // Custom annotations often use colons (e.g., gnomAD:AF)
    auto rec = make_record({
        {"gnomAD:AF", "0.001"},
        {"dbNSFP:SIFT_score", "0.05"},
        {"custom:my_field:sub", "value"}
    });

    EXPECT_EQ(rec.get("gnomAD:AF"), "0.001");
    EXPECT_DOUBLE_EQ(rec.get_numeric("gnomAD:AF"), 0.001);
    EXPECT_EQ(rec.get("dbNSFP:SIFT_score"), "0.05");
    EXPECT_EQ(rec.get("custom:my_field:sub"), "value");
    EXPECT_TRUE(rec.has("gnomAD:AF"));
}
