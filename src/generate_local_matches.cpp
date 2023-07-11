#include <random>

#include <seqan3/argument_parser/all.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/io/sequence_file/output.hpp>
#include <seqan3/alphabet/views/complement.hpp>

struct my_traits : seqan3::sequence_file_input_default_traits_dna
{
    using sequence_alphabet = seqan3::dna4;
};

struct cmd_arguments
{
    std::filesystem::path ref_path{};
    std::filesystem::path out_path{};
    double max_error_rate{0.01};
    uint32_t min_match_length{50u};
    uint32_t max_match_length{200u};
    uint32_t total_num_matches{1ULL<<20};
    uint64_t ref_len;
    bool reverse{false};
    uint32_t seed{42};

    bool verbose_ids{false};
};

void run_program(cmd_arguments const & arguments)
{
    std::random_device rd;
    // std::mt19937_64 rng(rd());
    std::mt19937_64 rng(arguments.seed);
    std::uniform_int_distribution<uint8_t> dna4_rank_dis(0, 3);
    std::uniform_int_distribution<> match_len_dis(arguments.min_match_length, arguments.max_match_length);

    uint32_t match_counter{};

    // Immediately invoked initialising lambda expession (IIILE).
    std::filesystem::path const out_file = [&]
                                               {
                                                   std::filesystem::path out_file = arguments.out_path;
                                                   out_file /= arguments.ref_path.stem();
                                                   out_file += ".fastq";
                                                   return out_file;
                                               }();

    seqan3::sequence_file_input<my_traits, seqan3::fields<seqan3::field::seq, seqan3::field::id>> fin{arguments.ref_path};
    seqan3::sequence_file_output fout{out_file};

    uint32_t num_matches = arguments.total_num_matches;
    if (arguments.reverse)
        num_matches = arguments.total_num_matches / 2;


    auto sample_matches = [&](auto const & seq, auto const & reference_name, uint32_t const & num_matches, bool const reverse = false)
    {
        uint64_t const seq_len = std::ranges::size(seq);
        uint32_t per_seq_matches = num_matches;
        if (seq.size() < arguments.ref_len - 1)  // if more than one chromosome
            per_seq_matches = std::round(num_matches * (double) seq_len / (double) arguments.ref_len);

        std::uniform_int_distribution<uint64_t> match_start_dis(0, seq_len - arguments.max_match_length);
        for (uint32_t current_match_number = 0; current_match_number < num_matches; ++current_match_number, ++match_counter)
        {
            uint32_t match_length = match_len_dis(rng);
            std::vector<seqan3::phred42> const quality(match_length, seqan3::assign_rank_to(40u, seqan3::phred42{}));
            std::uniform_int_distribution<uint32_t> match_error_position_dis(0, match_length - 1);
            uint64_t const match_start_pos = match_start_dis(rng);
            std::vector<seqan3::dna4> match = seq |
                                              seqan3::views::slice(match_start_pos, match_start_pos + match_length) |
                                              seqan3::views::to<std::vector>;

            uint32_t max_errors = (int) (match_length * arguments.max_error_rate);  // convert error rate to error count and round down
            for (uint8_t error_count = 0; error_count < max_errors; ++error_count)
            {
                uint32_t const error_pos = match_error_position_dis(rng);
                seqan3::dna4 const current_base = match[error_pos];
                seqan3::dna4 new_base = current_base;
                while (new_base == current_base)
                    seqan3::assign_rank_to(dna4_rank_dis(rng), new_base);
                match[error_pos] = new_base;
            }

            std::string query_id = std::to_string(match_counter);
            std::string meta_info{};

            if (arguments.verbose_ids)
            {
                meta_info += ' ';
                if (reverse)
                    meta_info += "reverse,";
                meta_info += "start_position=" + std::to_string(match_start_pos);
                meta_info += ",length=" + std::to_string(match_length);
                meta_info += ",errors=" + std::to_string(max_errors);
                meta_info += ",reference_id='" + reference_name + "'";
                meta_info += ",reference_file='" + std::string{arguments.ref_path} + "'";
            }
            fout.emplace_back(match, query_id + meta_info, quality);
        }
    };

    for (auto const & [seq, reference_name] : fin)
        sample_matches(seq, reference_name, num_matches);

    if (arguments.reverse)
    {
        seqan3::sequence_file_input<my_traits, seqan3::fields<seqan3::field::seq, seqan3::field::id>> fin{arguments.ref_path};
        for (auto const & [seq, reference_name] : fin)
        {
            std::vector<seqan3::dna4> compl_seq;
            for (auto const & c : seq | std::views::reverse | seqan3::views::complement)
                compl_seq.emplace_back(c);
            sample_matches(compl_seq, reference_name, num_matches, true);
        }
    }
}

void initialise_argument_parser(seqan3::argument_parser & parser, cmd_arguments & arguments)
{
    parser.info.author = "Enrico Seiler";
    parser.info.author = "enrico.seiler@fu-berlin.de";
    parser.info.short_description = "Generate local matches from reference.";
    parser.info.version = "0.0.1";
    parser.info.examples = {"./generate_local_matches --output ./matches_e2 ./big_dataset/64/bins/bin_{00..63}.fasta"};
    parser.add_positional_option(arguments.ref_path,
                                 "Provide path to reference sequence.");
    parser.add_option(arguments.out_path,
                      '\0',
                      "output",
                      "Provide the base dir where the matchs should be written to.",
                      seqan3::option_spec::required);
    parser.add_option(arguments.max_error_rate,
                      '\0',
                      "max-error-rate",
                      "The maximum number of errors.",
                      seqan3::option_spec::required,
		              seqan3::arithmetic_range_validator{0, 1});

    parser.add_option(arguments.min_match_length,
                      '\0',
                      "min-match-length",
                      "The minimum match length.");
    parser.add_option(arguments.max_match_length,
                      '\0',
                      "max-match-length",
                      "The maximum match length.");
    parser.add_option(arguments.total_num_matches,
                      '\0',
                      "num-matches",
                      "The number of matches.");
    parser.add_option(arguments.ref_len,
                      '\0',
                      "ref-len",
                      "Length of the reference.",
		              seqan3::option_spec::required);
    parser.add_flag(arguments.reverse,
                    'r',
                    "reverse",
                    "Also simulate matches from reverse strand.");
    parser.add_option(arguments.seed,
                    '\0',
                    "seed",
                    "Seed for random generator.");
    parser.add_flag(arguments.verbose_ids,
                    '\0',
                    "verbose-ids",
                    "Puts position information into the ID (where the match was sampled from)");
}

int main(int argc, char ** argv)
{
    seqan3::argument_parser myparser{"build_ibf", argc, argv, seqan3::update_notifications::off};
    cmd_arguments arguments{};
    initialise_argument_parser(myparser, arguments);
    try
    {
         myparser.parse();
    }
    catch (seqan3::argument_parser_error const & ext)
    {
        std::cout << "[Error] " << ext.what() << "\n";
        return -1;
    }

    run_program(arguments);

    return 0;
}
