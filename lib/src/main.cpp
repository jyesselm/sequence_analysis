#include <iostream>
#include <fstream>

#include <args.hxx>

#include <plog/Log.h>
#include <plog/Appenders/ColorConsoleAppender.h>

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// logging
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

namespace plog {
class CustomFormatter {
public:
    static util::nstring header() { return util::nstring(); }

    static util::nstring format(const Record & record) {
        tm t;
        util::localtime_s(&t, &record.getTime().time);

        util::nostringstream ss;
        ss << std::setfill(PLOG_NSTR('0')) << std::setw(2) << t.tm_hour << PLOG_NSTR(":")
           << std::setfill(PLOG_NSTR('0')) << std::setw(2) << t.tm_min << PLOG_NSTR(":") << std::setfill(PLOG_NSTR('0'))
           << std::setw(2) << t.tm_sec << PLOG_NSTR(" ");
        ss << std::setfill(PLOG_NSTR(' ')) << std::setw(5) << std::left << severityToString(record.getSeverity())
           << PLOG_NSTR(" ");
        ss << PLOG_NSTR("[") << record.getFunc() << PLOG_NSTR("@") << record.getLine() << PLOG_NSTR("] ");
        ss << record.getMessage() << PLOG_NSTR("\n");

        return ss.str();
    }
};
}

void
init_logging(
        plog::Severity const & log_level) {
    static plog::ColorConsoleAppender<plog::CustomFormatter> consoleAppender;
    plog::init(log_level, &consoleAppender);
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// command line parsing
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct Parameters {
    std::string nomod_fastq;
    std::string dms_fasta;
    std::string source_seqs;
};

void parse_command_line(
        int argc,
        char **argv,
        Parameters & parameters) {
    args::ArgumentParser parser("This is a test program.", "This goes after the options.");
    args::HelpFlag help(parser, "help", "Display this help menu", {'h', "help"});
    args::ValueFlag<std::string> nomod(parser, "", "path to nomod fastq", {'n', "nomod"});
    args::ValueFlag<std::string> dms(parser, "", "path to dms fastq", {'d', "dms"});
    args::ValueFlag<std::string> source_seqs(parser, "", "path to source seqs", {'s', "source"});

    try {
        parser.ParseCLI(argc, argv);
    }
    catch (args::Help const &) {
        std::cout << parser;
        exit(0);
    }
    catch (args::ParseError const & e) {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        exit(1);
    }
    catch (args::ValidationError const & e) {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        exit(1);
    }

    if(! (nomod && dms && source_seqs)) {
        std::cout << parser;
        LOG_ERROR << "both -n, -d and -s are required" << std::endl;
        exit(1);
    }

    parameters.nomod_fastq = args::get(nomod);
    parameters.dms_fasta   = args::get(dms);
    parameters.source_seqs = args::get(source_seqs);

}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// classes
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

class Sequence {
public:
    explicit
    Sequence(
            std::string & sequence):
            _sequence(std::move(sequence)) {
        _generate_res_codes();
    }

public:
    size_t
    length() const { return _sequence.length(); }

    inline
    int
    get_res_code(
            int pos) const {
        return _res_codes[pos];
    }

    inline
    std::string const &
    str() const {
        return _sequence;
    }

private:
    void
    _generate_res_codes() {
        for(auto const & e : _sequence) {
            if     (e == 'A') { _res_codes.push_back(0); }
            else if(e == 'C') { _res_codes.push_back(1); }
            else if(e == 'G') { _res_codes.push_back(2); }
            else if(e == 'T') { _res_codes.push_back(3); }
            else if(e == 'U') { _res_codes.push_back(3); }
            else              { _res_codes.push_back(4); }
        }
    }

private:
    std::string _sequence;
    std::vector<int> _res_codes;
};

typedef std::vector<Sequence> Sequences;

struct Alignment {
    int pos;
    int mismatches;
};

typedef std::vector<Alignment> Alignments;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// main functions
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::string & ltrim(std::string & s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(),
                                    std::ptr_fun<int, int>(std::isgraph)));
    return s;
}

std::string & rtrim(std::string & s) {
    s.erase(std::find_if(s.rbegin(), s.rend(),
                         std::ptr_fun<int, int>(std::isgraph)).base(), s.end());
    return s;
}

std::string & trim(std::string & s) {
    return ltrim(rtrim(s));
}

void
parse_reads_from_fastq(
        std::string const & filename,
        Sequences & sequences /*return*/) {
    auto previous_line = std::string();
    auto line = std::string();
    auto in = std::ifstream();
    in.open(filename);
    getline(in, previous_line);
    while(in.good()) {
        getline(in, line);
        if (previous_line[0] == '@' && line[0] != '@') {
            sequences.push_back(Sequence(trim(line)));
        }
        previous_line = line;
    }
    in.close();
    LOG_INFO << sequences.size() << " read(s) loaded from fastq file: " << filename;
}

void
get_source_sequences(
        std::string const & filename,
        Sequences & sequences /*return*/) {
    auto line = std::string();
    auto in = std::ifstream();
    in.open(filename);
    while(in.good()) {
        getline(in, line);
        if(line.length() == 0) { continue; }
        sequences.push_back(Sequence(trim(line)));
    }
    LOG_INFO << sequences.size() << " source sequence(s) from file: " << filename;
}

bool
get_best_pairwise_alignment(
        Sequence const & seq_1,
        Sequence const & seq_2,
        int threshold,
        Alignment & alignment /*return*/) {
    alignment.mismatches = threshold;
    auto best_score = threshold;
    auto best_pos = 0;
    auto score = 0;
    for(int i = 0; i < seq_1.length(); i++) {
        score = 0;
        if(seq_2.length() >= seq_1.length() - i) { break; }
        for(int j = 0; j < seq_2.length(); j++) {
            if(seq_1.get_res_code(i+j) != seq_2.get_res_code(j)) {
                score += 1;
                if(score > threshold) { break; }
            }
        }
        if(score < best_score) {
            best_pos = i;
            best_score = score;
        }
    }
    if(best_score < threshold) {
        alignment.pos = best_pos;
        alignment.mismatches = best_score;
        return true;
    }
    return false;
}

void
get_seq_reads(
        Sequences const & source_seqs,
        Sequences const & reads,
        int mismatch_threshold,
        std::vector<Sequences> & seq_reads /*return*/,
        std::vector<Alignments> & seq_alignments /*return*/) {

    auto best_i = 0;
    auto i = 0;
    auto aligned = false;
    auto alignment = Alignment();
    auto best_alignment = Alignment();
    auto j = -1;
    for(auto const & read: reads) {
        j++;
        best_alignment.mismatches = 100;
        best_i = 0;
        i = -1;
        for(auto const & seq : source_seqs) {
            i++;
            aligned = get_best_pairwise_alignment(read, seq, 20, alignment);
            if(!aligned) { continue; }
            if(alignment.mismatches < best_alignment.mismatches) {
                best_alignment = alignment;
                best_i = i;
            }
        }
        if(best_alignment.mismatches > mismatch_threshold) { continue; }
        seq_reads[best_i].push_back(read);
        seq_alignments[best_i].push_back(alignment);
    }
}


int main(
        int argc,
        char **argv) {
    auto parameters = Parameters();
    init_logging(plog::Severity::info);
    parse_command_line(argc, argv, parameters);

    auto source_seqs = Sequences();
    auto dms_reads = Sequences();
    auto nomod_reads = Sequences();

    get_source_sequences(parameters.source_seqs, source_seqs);
    parse_reads_from_fastq(parameters.dms_fasta, dms_reads);
    parse_reads_from_fastq(parameters.nomod_fastq, nomod_reads);
    auto dms_seq_reads = std::vector<Sequences>(source_seqs.size());
    auto dms_alignments = std::vector<Alignments>(source_seqs.size());

    get_seq_reads(source_seqs, dms_reads, 3, dms_seq_reads, dms_alignments);
    int i = 0;
    for(auto const & seq : source_seqs) {
        LOG_INFO << seq.str() << " " << dms_seq_reads[i].size();
        i++;
    }


    return 0;
}

