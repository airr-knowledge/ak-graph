#include <iostream>
#include <fstream>
#include <pqxx/pqxx>
#include <compairr/api.h>
#include <vector>
#include <string>
#include "utils/timer.h"

using std::cout;
using std::endl;
using std::ofstream;

struct mem_file {
    FILE *fp;
    char *mem;
};

mem_file query_and_stream(void) {
    // Connect to the database
    pqxx::connection cx{"postgresql://postgres:example@ak-db/airrkb_v1"};
    pqxx::work tx{cx};

    std::string buffer = "repertoire_id\tsequence_id\tduplicate_count\tjunction_aa\n";
    std::string repertoire_id = "rep1"; // same for all rows
    for (auto [int_index, junction_aa] :
        tx.stream<std::string, std::string>(
            "SELECT int_index, junction_aa "
            "FROM unique_junctions"
        )
    ) {
        buffer += repertoire_id + "\t";         // repertoire_id
        buffer += int_index + "\t";             // sequence_id
        buffer += "1\t";                        // duplicate_count
        buffer += junction_aa;                  // junction_aa
        buffer += "\n";
    }

    tx.commit();

    // construct a file in memory
    char *mem = (char *)malloc(buffer.size());
    memcpy(mem, buffer.data(), buffer.size());

    FILE *fp = fmemopen(
        mem,
        buffer.size(),
        "r"
    );

    if (!fp) {
        perror("fmemopen");
        exit(1);
    }

    return {fp, mem};
}

void free_mem_file(mem_file mf) {
    fclose(mf.fp);
    free(mf.mem);
}

void run_compairr(mem_file mf) {
    std::vector<std::string> argv_strs;
    argv_strs.push_back("filler for prog_name, doesn't matter but needs to be here");
    argv_strs.push_back("-m");                  // compute overlap matrix
    argv_strs.push_back("repertoire.tsv");      // filename 1
    argv_strs.push_back("repertoire.tsv");      // filename 2
    argv_strs.push_back("-g");                  // ignore genes
    argv_strs.push_back("-u");                  // ignore sequences with unknowns (X)
    argv_strs.push_back("-d");                  // differences
    argv_strs.push_back("1");
    argv_strs.push_back("-i");                  // indels
    argv_strs.push_back("-l");                  // log filename
    argv_strs.push_back("compairr.log");
    argv_strs.push_back("-o");                  // matrix output filename
    argv_strs.push_back("output_matrix.tsv");
    argv_strs.push_back("--no-matrix");         // no matrix output
    argv_strs.push_back("-p");                  // pairs output filename
    argv_strs.push_back("pairs_with_indels.tsv");
    argv_strs.push_back("-t");                  // threads
    argv_strs.push_back("8");

    std::vector<char*> argv_vec;
    for (auto &s : argv_strs) {
        argv_vec.push_back(&s[0]);
    }

    int argc = argv_vec.size();
    char **argv = argv_vec.data();
    
    logfile = stderr;

    arch_srandom(1);

    args_init(argc, argv);

    open_files();

    db_init();
    db * d1 = db_create();

    Timer db1_read_time;
    db1_read_time.start();
    db_read(d1, mf.fp, true, "1");
    std::cout << "\tdb_read() for db1:   ";
    db1_read_time.view(std::cout); 
    std::cout << std::endl;
    
    db * d2 = db_create();
    d2 = d1;

    Timer co_timer;
    co_timer.start();
    compairr::computeOverlap(d1, d2, true);
    std::cout << "\tcomputeOverlap():  ";
    co_timer.view(std::cout);
    std::cout << std::endl;

    close_files();
}

int main(void) {
    Timer qs_time, rc_time;
    
    qs_time.start();
    mem_file mf = query_and_stream();
    std::cout << "query_and_stream():  ";
    qs_time.view(std::cout); 
    std::cout << std::endl;

    rc_time.start();
    run_compairr(mf);
    std::cout << "run_compairr():      ";
    rc_time.view(std::cout);
    std::cout << std::endl;

    return 0;
}
// time docker run -v $(pwd):/data compairr -m -g -t 100 -d 1 --distance /data/dedup_repertoire.tsv -p /data/pairs_2.tsv