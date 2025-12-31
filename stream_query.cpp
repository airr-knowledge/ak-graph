#include <iostream>
#include <fstream>
#include <pqxx/pqxx>
#include <compairr/api.h>
#include <vector>
#include <string>

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

    // Query akc_id and junction_aa
    // pqxx::result r{tx.exec(
    //     "SELECT akc_id, junction_aa, v_call, j_call "
    //     "FROM \"unique_junctions\" "
    //     // "INNER JOIN \"Chain\" c ON tcr.trb_chain = c.akc_id "
    //     "LIMIT 3000"
    // )};
    std::string buffer = "repertoire_id\tsequence_id\tduplicate_count\tjunction_aa\n";
    std::string repertoire_id = "rep1"; // same for all rows
    for (auto [int_index, junction_aa] :
        tx.stream<std::string, std::string>(
            "SELECT int_index, junction_aa "
            "FROM unique_junctions"
        )
    ) {
        buffer += repertoire_id + "\t";         // repertoire_id
        buffer += int_index + "\t"; // sequence_id
        buffer += "1\t";                         // duplicate_count
        buffer += junction_aa;    // junction_aa
        buffer += "\n";
    }

    tx.commit();

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
    argv_strs.push_back("output_pairs.tsv");

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
    // std::cout << "d1_create ---------------------" << std::endl;
    // db_debug_print(d1, std::cout);
    // std::cout << 
    db_read_fp(d1, mf.fp, true, "1");
    // std::cout << "d1_fill ---------------------" << std::endl;
    // db_debug_print(d1, std::cout);
    
    db * d2 = db_create();
    // db_read(d2, "repertoire.tsv", true, "2");
    d2 = d1;
    // std::cout << "d2_fill ---------------------" << std::endl;
    // db_debug_print(d2, std::cout);

    compairr::computeOverlap(d1, d2, true);

    close_files();
}

int main(void) {
    mem_file mf = query_and_stream();
    run_compairr(mf);
    // free_mem_file(mf);

    return 0;
}
// time docker run -v $(pwd):/data compairr -m -g -t 100 -d 1 --distance /data/dedup_repertoire.tsv -p /data/pairs_2.tsv