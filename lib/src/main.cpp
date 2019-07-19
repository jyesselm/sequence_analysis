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
    length() { return _sequence.length(); }

    inline
    int
    get_res_code(
            int pos) {
        return _res_codes[pos];
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

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
// main functions
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
parse_reads_from_fastq(
        std::string const & filename,
        Sequences & sequences /*return*/) {
    auto line = std::string();
    auto in = std::ifstream();
    in.open(filename);
    while(in.good()) {
        getline(in, line);
        if(line.length() == 0) { continue; }
        if(line[0] == '@') {
            getline(in, line);
            if(line.length() == 0) { continue; }
            if(line[0] != '@') {
                sequences.push_back(Sequence(line));
            }
        }
    }
    LOG_INFO << sequences.size() << " reads loaded from fastq file: " << filename;
}

void
get_source_sequences(
        std::string const & filename,
        Sequences & sequences /*return*/) {

}


int main(
        int argc,
        char **argv) {
    auto parameters = Parameters();
    init_logging(plog::Severity::info);
    parse_command_line(argc, argv, parameters);

    auto dms_reads = Sequences();
    auto nomod_reads = Sequences();

    parse_reads_from_fastq(parameters.dms_fasta, dms_reads);
    parse_reads_from_fastq(parameters.nomod_fastq, nomod_reads);

    return 0;
}

