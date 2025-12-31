#include <algorithm>
#include <random>
#include <ranges>

#include <seqan3/argument_parser/all.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/io/sequence_file/output.hpp>
#include <seqan3/alphabet/views/complement.hpp>
#include <seqan3/core/debug_stream.hpp>

struct my_traits : seqan3::sequence_file_input_default_traits_dna
{
    using sequence_alphabet = seqan3::dna4;
};

struct cmd_arguments
{
    std::filesystem::path ref_path{};
    std::filesystem::path matches_out_path{};
    std::filesystem::path genome_out_path{};
    std::filesystem::path query_path{};
    double max_error_rate{0.01};
    uint32_t min_match_length{50u};
    uint32_t max_match_length{200u};
    uint32_t total_num_matches{1ULL<<20};
    uint64_t ref_len;
    bool uniform{true};
    bool normal{false};
    double std_dev_fraction{0.1};

    bool reverse{false};
    uint32_t seed{42};

    bool verbose_ids{false};
};

template <typename func_t>
void runtime_to_compile_time(func_t const & func, bool b1)
{
    if (b1) 
        func.template operator()<true>();
    else
        func.template operator()<false>();
}

template <typename dist_T>
dist_T get_match_start_dis(uint64_t const seq_size, cmd_arguments const & arguments)
{
	if (std::is_same<dist_T, std::uniform_int_distribution<uint64_t>>::value)
		return dist_T{0, seq_size - arguments.max_match_length};
	else
	{
		const uint64_t seq_mean = std::round((seq_size - arguments.max_match_length) / 2.0);
        	// seqan3::debug_stream << "seq_mean\t" << seq_mean << '\n'; 
		return dist_T{(double) seq_mean, seq_mean * arguments.std_dev_fraction};
	}
}

template <bool is_uniform>
void run_program(cmd_arguments const & arguments)
{
    std::random_device rd;
    // std::mt19937_64 rng(rd());
    std::mt19937_64 rng(arguments.seed);
    std::uniform_int_distribution<uint8_t> dna4_rank_dis(0, 3);
    std::uniform_int_distribution<> match_len_dis(arguments.min_match_length, arguments.max_match_length);

    uint32_t num_matches = arguments.total_num_matches;
    if (arguments.reverse)
        num_matches = arguments.total_num_matches / 2;

    seqan3::sequence_file_input<my_traits, seqan3::fields<seqan3::field::seq, seqan3::field::id>> fref{arguments.ref_path};
    using record_type = typename decltype(fref)::record_type;
    std::vector<record_type> matches;
    matches.reserve(arguments.total_num_matches);

    uint32_t match_counter{};
    auto sample_matches = [&](auto const & seq, auto const & reference_name, uint32_t const & num_matches, bool const reverse = false)
    {
        uint32_t per_seq_matches = num_matches;
        if (seq.size() < arguments.ref_len - 1)  // if more than one chromosome
	{
	    //seqan3::debug_stream << "seq.size()\t" << seq.size() << '\n';
	    //seqan3::debug_stream << "arguments.ref_len\t" << arguments.ref_len  << '\n';
	    //seqan3::debug_stream << "num_matches\t" << num_matches << '\n';
	    per_seq_matches = std::round(num_matches * (double) seq.size() / (double) arguments.ref_len);
	}
        //seqan3::debug_stream << "Simulating " << per_seq_matches << " matches from sequence " << reference_name << '\n'; 
	using start_dis_T = std::conditional<is_uniform, 
	      				     std::uniform_int_distribution<uint64_t>, 
	      				     std::normal_distribution<double>>::type;
	start_dis_T match_start_dis = get_match_start_dis<start_dis_T>(seq.size(), arguments);

	for (uint32_t current_match_number = 0; current_match_number < num_matches; ++current_match_number, ++match_counter)
        {
            uint32_t match_length = match_len_dis(rng);
            std::uniform_int_distribution<uint32_t> match_error_position_dis(0, match_length - 1);
            // uint64_t const match_start_pos = match_start_dis(rng);
            uint64_t const match_start_pos = std::clamp<uint64_t>(std::round(match_start_dis(rng)), 0,  seq.size() - arguments.max_match_length);
            std::vector<seqan3::dna4> match;
            for (auto & n : seq | seqan3::views::slice(match_start_pos, match_start_pos + match_length))
                match.push_back(n);

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
            matches.emplace_back(match, query_id + meta_info);
        }

	//seqan3::debug_stream << "emplaced_back matches\n";
    };

    for (auto const & [seq, reference_name] : fref)
        sample_matches(seq, reference_name, num_matches);

    if (arguments.reverse)
    {
        seqan3::sequence_file_input<my_traits, seqan3::fields<seqan3::field::seq, seqan3::field::id>> fref{arguments.ref_path};
        for (auto const & [seq, reference_name] : fref)
        {
            std::vector<seqan3::dna4> compl_seq;
            for (auto const & c : seq | std::views::reverse | seqan3::views::complement)
                compl_seq.emplace_back(c);
            sample_matches(compl_seq, reference_name, num_matches, true);
        }
    }

    //seqan3::debug_stream << "reversed\n";

    seqan3::sequence_file_output fout_matches{arguments.matches_out_path};
    for (auto & match : matches)
        fout_matches.push_back(match);

    //seqan3::debug_stream << "pushed_back\n";

    if (!arguments.query_path.empty())
    {
        seqan3::sequence_file_input<my_traits, seqan3::fields<seqan3::field::seq, seqan3::field::id>> fquery{arguments.query_path};
        uint64_t total_query_len{0};

        std::vector<seqan3::dna4_vector> query_sequences;
        std::vector<std::string> query_ids;    
        for (auto & [seq, query_name] : fquery)
        {
            total_query_len += seq.size();
            query_sequences.emplace_back(std::move(seq));
            query_ids.emplace_back(std::move(query_name));
        }
        std::vector<bool> mut_mask(total_query_len, 0);

        std::uniform_int_distribution<> match_insertion_loc_dis(0, total_query_len - arguments.max_match_length);
        for (uint32_t i = 0; i < arguments.total_num_matches; i++)
        {
            uint8_t failure_counter{0};
            uint64_t elapsed_length{0};
            size_t query_ind{0};

            auto loc = match_insertion_loc_dis(rng);
            auto [match, match_id] = matches[i];

            while (std::any_of(mut_mask.begin() + loc, 
                               mut_mask.begin() + loc + match.size(), 
                               [](bool x) { return x; }))
            {
                loc = match_insertion_loc_dis(rng);
                failure_counter++;
                if (failure_counter > 10)
                    throw std::runtime_error("Can not insert another local match without overlapping existing ones.");
            }
            
            while (loc - elapsed_length + match.size() >= query_sequences[query_ind].size())
            {
                elapsed_length += query_sequences[query_ind].size();
                query_ind++;
            }
            
            auto & seq = query_sequences[query_ind];
            for (auto & nuc : match)
                std::cout << nuc.to_char();

            std::cout << '\n' << match_id << '\n';
            for (size_t l{0}; l < match.size(); l++)
            {
                seq[loc - elapsed_length + l] = match[l];
                mut_mask[loc + l] = 1;
                /*std::cout << std::to_string(loc) << '\t' << std::to_string(elapsed_length) << '\t' 
                          << std::to_string(l) << '\t' <<  std::to_string(loc - elapsed_length + l) << '\t' 
                          << match[l].to_char() << '\t' << seq[loc - elapsed_length + l].to_char() << '\t' << query_sequences[query_ind][loc - elapsed_length + query_ind].to_char() << '\n';
                */
            }
        }

        seqan3::sequence_file_output fout_genome{arguments.genome_out_path};
        for (size_t i{0}; i < query_sequences.size(); i++)
            fout_genome.emplace_back(query_sequences[i], query_ids[i]);


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
    parser.add_option(arguments.matches_out_path,
                      '\0',
                      "matches-out",
                      "Provide the path to the local alignment FASTA output.",
                      seqan3::option_spec::required);
    parser.add_option(arguments.genome_out_path,
                      '\0',
                      "genome-out",
                      "Provide the path to the local alignment FASTA output.");
    parser.add_option(arguments.query_path,
                      '\0',
                      "query",
                      "Provide the query sequence where the local matches should be inserted into.");
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
    parser.add_flag(arguments.normal,
                    '\0',
                    "normal",
                    "Local alignment positions are normally distributed across the sequence. \n"
		    "Default: uniform distribution of start positions.");
    parser.add_option(arguments.std_dev_fraction,
		    '\0',
		    "std-dev-frac",
		    "Set the standard deviation of the normal distribution as a fraction of the mean."
		    "Default: std dev 10% of the mean.");
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

    if (arguments.normal)
        arguments.uniform = false;

    runtime_to_compile_time([&]<bool is_uniform>()
    {
        run_program<is_uniform>(arguments);
    }, arguments.uniform);

    return 0;
}
