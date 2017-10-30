#include <cmath>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <functional>
#include <getopt.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdint.h>
#include <string>
#include <sys/stat.h>
#include <utility>
#include <vector>

#include <thread>
#include <atomic>

#include <jemalloc/jemalloc.h>

#include "CompactedDBG.hpp"
#include "Common.hpp"

using namespace std;

struct ProgramOptions {

    bool ref;
    bool verbose;
    bool clipTips;
    bool deleteIsolated;

    size_t threads;
    size_t k, g;
    size_t read_chunksize;
    size_t contig_size;
    size_t bf, bf2;
    size_t nkmers, nkmers2;

    string prefixFilenameGFA;
    string filenameGFA;
    string inFilenameBBF;
    string outFilenameBBF;

    vector<string> files;

    ProgramOptions() :  threads(1), k(31), g(23), nkmers(0), nkmers2(0), bf(14), bf2(14), read_chunksize(10000), \
                        contig_size(1000000), ref(false), verbose(false), clipTips(false), deleteIsolated(false), \
                        inFilenameBBF(""), outFilenameBBF("") {}
};

// use:  PrintVersion();
// post: The version of the program has been printed to cout
void PrintVersion() {
    cout <<  BFG_VERSION << endl;
}

// use:  PrintCite();
// post: Information of how to cite this software has been printed to cerr
void PrintCite() {
    cout << "The paper describing this software has not been published." << endl;
}

void PrintUsage() {

    cout << endl << "BFGraph " << BFG_VERSION << endl << endl;
    cout << "Highly Parallel and Memory Efficient Compacted de Bruijn Graph Construction" << endl << endl;
    cout << "Usage: BFGraph [Parameters] FAST(A|Q)_file_1 ..." << endl << endl;
    cout << "Parameters with required argument:" << endl << endl <<
    "  -n, --num-kmers          [MANDATORY] Estimated number (upper bound) of different k-mers in the FASTA/FASTQ files" << endl <<
    "  -N, --num-kmer2          [MANDATORY] Estimated number (upper bound) of different k-mers occurring twice or more in the FASTA/FASTQ files" << endl <<
    "  -o, --output             [MANDATORY] Prefix for output GFA file" << endl <<
    "  -t, --threads            Number of threads (default is 1)" << endl <<
    "  -k, --kmer-length        Length of k-mers (default is 31)" << endl <<
    "  -g, --min-length         Length of minimizers (default is 23)" << endl <<
    "  -b, --bloom-bits         Number of Bloom filter bits per k-mer occurring at least once in the FASTA/FASTQ files (default is 14)" << endl <<
    "  -B, --bloom-bits2        Number of Bloom filter bits per k-mer occurring at least twice in the FASTA/FASTQ files (default is 14)" << endl <<
    "  -l, --load               Filename for input Blocked Bloom Filter, skips filtering step (default is no input)" << endl <<
    "  -f, --output2            Filename for output Blocked Bloom Filter (default is no output)" << endl <<
    "  -s, --chunk-size         Read chunksize to split between threads (default is 10000)" << endl <<
    endl << "Parameters with no argument:" << endl << endl <<
    "      --ref                Reference mode, no filtering" << endl <<
    "  -c, --clip-tips          Clip tips shorter than k k-mers in length" << endl <<
    "  -r, --rm-isolated        Delete isolated contigs shorter than k k-mers in length" << endl <<
    "  -v, --verbose            Print information messages during construction" << endl <<
    endl;
}

void parse_ProgramOptions(int argc, char **argv, ProgramOptions& opt) {

    const char* opt_string = "n:N:o:t:k:g:b:B:l:f:s:crv";

    static struct option long_options[] = {

        {"num-kmers",       required_argument,  0, 'n'},
        {"num-kmers2",      required_argument,  0, 'N'},
        {"output",          required_argument,  0, 'o'},
        {"threads",         required_argument,  0, 't'},
        {"kmer-length",     required_argument,  0, 'k'},
        {"min-length",      required_argument,  0, 'g'},
        {"bloom-bits",      required_argument,  0, 'b'},
        {"bloom-bits2",     required_argument,  0, 'B'},
        {"load",            required_argument,  0, 'l'},
        {"output2",         required_argument,  0, 'f'},
        {"chunk-size",      required_argument,  0, 's'},
        {"ref",             no_argument,        0,  0 },
        {"clip-tips",       no_argument,        0, 'c'},
        {"rm-isolated",     no_argument,        0, 'r'},
        {"verbose",         no_argument,        0, 'v'},
        {0,                 0,                  0,  0 }
    };

    int option_index = 0, c;

    while ((c = getopt_long(argc, argv, opt_string, long_options, &option_index)) != -1) {

        switch (c) {

            case 0:
                if (strcmp(long_options[option_index].name, "ref") == 0) opt.ref = true;
                break;
            case 'v':
                opt.verbose = true;
                break;
            case 'c':
                opt.clipTips = true;
                break;
            case 'r':
                opt.deleteIsolated = true;
                break;
            case 't':
                opt.threads = atoi(optarg);
                break;
            case 'k':
                opt.k = atoi(optarg);
                break;
            case 'g':
                opt.g = atoi(optarg);
                break;
            case 's':
                opt.read_chunksize = atoi(optarg);
                break;
            case 'n':
                opt.nkmers = atoi(optarg);
                break;
            case 'N':
                opt.nkmers2 = atoi(optarg);
                break;
            case 'b':
                opt.bf = atoi(optarg);
                break;
            case 'B':
                opt.bf2 = atoi(optarg);
                break;
            case 'o':
                opt.prefixFilenameGFA = optarg;
                break;
            case 'f':
                opt.outFilenameBBF = optarg;
                break;
            case 'l':
                opt.inFilenameBBF = optarg;
                break;
            default: break;
        }
    }

    // all other arguments are fast[a/q] files to be read
    while (optind < argc) opt.files.push_back(argv[optind++]);
}

bool check_ProgramOptions(ProgramOptions& opt) {

    bool ret = true;

    size_t max_threads = std::thread::hardware_concurrency();

    if (opt.threads <= 0){

        cerr << "Error: Number of threads cannot be less than or equal to 0" << endl;
        ret = false;
    }
    else if (opt.threads > max_threads){

        cerr << "Error: Number of threads cannot be greater than or equal to " << max_threads << endl;
        ret = false;
    }

    if (opt.read_chunksize <= 0) {

        cerr << "Error: Chunk size of reads to share among threads cannot be less than or equal to 0" << endl;
        ret = false;
    }

    if (opt.k <= 0){

        cerr << "Error: Length k of k-mers cannot be less than or equal to 0" << endl;
        ret = false;
    }

    if (opt.k >= MAX_KMER_SIZE){

        cerr << "Error: Length k of k-mers cannot exceed or be equal to " << MAX_KMER_SIZE << endl;
        ret = false;
    }

    if (opt.g <= 0){

        cerr << "Error: Length g of minimizers cannot be less than or equal to 0" << endl;
        ret = false;
    }

    if (opt.g >= MAX_KMER_SIZE){

        cerr << "Error: Length g of minimizers cannot exceed or be equal to " << MAX_KMER_SIZE << endl;
        ret = false;
    }

    if (opt.nkmers <= 0){

        cerr << "Error: Number of Bloom filter bits per unique k-mer cannot be less than or equal to 0" << endl;
        ret = false;
    }

    if (!opt.ref && (opt.nkmers2 <= 0)){

        cerr << "Error: Number of Bloom filter bits per non unique k-mer cannot be less than or equal to 0" << endl;
        ret = false;
    }

    if (!opt.ref && (opt.nkmers2 > opt.nkmers)){

        cerr << "Error: The estimated number of non unique k-mers ";
        cerr << "cannot be greater than the estimated number of unique k-mers" << endl;
        ret = false;
    }

    if (opt.bf <= 0){

        cerr << "Error: Number of Bloom filter bits per unique k-mer cannot be less than or equal to 0" << endl;
        ret = false;
    }

    if (!opt.ref && (opt.bf2 <= 0)){

        cerr << "Error: Number of Bloom filter bits per non unique k-mer cannot be less than or equal to 0" << endl;
        ret = false;
    }

    if (opt.ref) {

        opt.bf2 = 0;
        opt.nkmers2 = 0;
    }

    if (opt.outFilenameBBF.length() != 0){

        FILE* fp = fopen(opt.outFilenameBBF.c_str(), "wb");

        if (fp == NULL) {

            cerr << "Error: Could not open file for writing output Blocked Bloom filter: " << opt.outFilenameBBF << endl;
            ret = false;
        }
        else fclose(fp);
    }

    if (opt.inFilenameBBF.length() != 0){

        FILE* fp = fopen(opt.inFilenameBBF.c_str(), "rb");

        if (fp == NULL) {

            cerr << "Error: Could not read file input Blocked Bloom filter: " << opt.inFilenameBBF << endl;
            ret = false;
        }
        else fclose(fp);
    }

    opt.filenameGFA = opt.prefixFilenameGFA + ".gfa";

    FILE* fp = fopen(opt.filenameGFA.c_str(), "w");

    if (fp == NULL) {

        cerr << "Error: Could not open file for writing output graph in GFA format: " << opt.filenameGFA << endl;
        ret = false;
    }
    else fclose(fp);

    if (opt.files.size() == 0) {

        cerr << "Error: Missing FASTA/FASTQ input files" << endl;
        ret = false;
    }
    else {

        struct stat stFileInfo;
        vector<string>::const_iterator it;
        int intStat;

        for(it = opt.files.begin(); it != opt.files.end(); ++it) {

            intStat = stat(it->c_str(), &stFileInfo);

            if (intStat != 0) {
                cerr << "Error: File not found, " << *it << endl;
                ret = false;
            }
        }
    }

    return ret;
}

/*class myInt {

    public:

        myInt(int int_init = 1) : Int(int_init) {}

        static myInt* joinData(const UnitigMap& um_tail, const UnitigMap& um_head, const CompactedDBG<myInt>& cdbg){

            myInt* join = new myInt;
            join->Int = cdbg.getData(um_tail)->Int + cdbg.getData(um_head)->Int;

            return join;
        }

        static vector<myInt*> splitData(const UnitigMap& a, const vector<pair<int,int>>& split_a, const CompactedDBG<myInt>& cdbg){

            vector<myInt*> v_return;

            for (size_t i = 0; i < split_a.size(); i++){

                v_return.push_back(new myInt);
                v_return[i]->Int = cdbg.getData(a)->Int / split_a.size();
            }

            return v_return;
        }

        int Int;
};*/

int main(int argc, char **argv){

    if (argc < 2) PrintUsage();
    else {

        ProgramOptions opt;

        parse_ProgramOptions(argc, argv, opt);

        if (check_ProgramOptions(opt)){ //Program options are valid

            /*int i = 1;
            for (unsigned char c = 0, c_tmp; i <= 256; i++, c++){

                c_tmp = (c << 6) | ((c << 2) & 0x30) | ((c >> 2) & 0xc) | (c >> 6);

                cout << "0x" << hex << (int) c_tmp << ",";

                if (i%16 == 0) cout << endl;
                else cout << " ";
            }

            exit(1);*/

            CompactedDBG<> cdbg;

            cdbg.setK(opt.k, opt.g);

            cdbg.build(opt.files, opt.nkmers, opt.nkmers2, opt.ref, opt.threads, nullptr,
                       opt.bf, opt.bf2, opt.inFilenameBBF, opt.outFilenameBBF,
                       opt.read_chunksize, 1000000, opt.verbose);

            cdbg.simplify(opt.deleteIsolated, opt.clipTips, opt.verbose);

            cdbg.write(opt.prefixFilenameGFA, opt.verbose);


            /*CompactedDBG<myInt> cdbg(&myInt::joinData, &myInt::splitData);
            cdbg.setK(opt.k, opt.g);

            string seq1 = "CAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
            string seq2 = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAC";

            cdbg.add(seq1, true);

            UnitigMap um1 = cdbg.find(Kmer(seq1.c_str()));

            const myInt* data1 = cdbg.getData(um1);

            cout << "Before: data1 = " << data1->Int << endl;


            cdbg.add(seq2, true);

            UnitigMap um2 = cdbg.find(Kmer(seq2.c_str()));

            const myInt* data2 = cdbg.getData(um2);

            cout << "data2 = " << data2->Int << endl;

            for (CompactedDBG<myInt>::iterator it = cdbg.begin(), it_end; it != it_end; it++){

                cout << cdbg.toString(*it) << endl;
            }*/
        }
    }
}

/*
// use:  PrintUsage();
// post: How to run BFGraph has been printed to cout
void PrintUsage() {
    cout << "BFGraph " << BFG_VERSION << endl << endl;
    cout << "Highly Parallel and Memory Efficient Compacted de Bruijn Graph Construction" << endl << endl;
    cout << "Usage: BFGraph <cmd> [options] ..." << endl << endl;
    cout << "Where <cmd> can be one of:" << endl;
    cout <<
        "    filter       Filters errors from reads" << endl <<
        "    contigs      Builds a compacted de Bruijn graph" << endl <<
        "    cite         Prints information for citing the paper" << endl <<
        "    version      Displays version number" << endl << endl;
    ;
}


int main(int argc, char **argv) {

    if (argc < 2) PrintUsage();
    else if (strcmp(argv[1], "cite") == 0) PrintCite();
    else if (strcmp(argv[1], "version") == 0) PrintVersion();
    else if (strcmp(argv[1], "filter") == 0) FilterReads(argc-1,argv+1);
    else if (strcmp(argv[1], "contigs") == 0) BuildContigs(argc-1,argv+1);
    else {

        cout << "Did not understand command " << argv[1] << endl;
        PrintUsage();
    }
}*/
