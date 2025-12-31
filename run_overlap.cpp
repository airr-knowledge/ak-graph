#include <iostream>
#include <fstream>
#include <pqxx/pqxx>
#include <compairr/api.h>
#include <vector>
#include <string>

using std::cout;
using std::endl;
using std::ofstream;

pqxx::result query_and_write_tsv(void) {
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
    pqxx::result r{tx.exec(
        "SELECT tcr.akc_id, c.junction_aa "
        "FROM \"TCellReceptor\" tcr "
        "INNER JOIN \"Chain\" c ON tcr.trb_chain = c.akc_id "
        "LIMIT 5"
    )};
    cout << "Query Done." << endl;
    if (r.empty()) {
        cout << "No rows found!" << endl;
        return r;
    }

    // Open TSV file
    ofstream tsv_file("repertoire.tsv");
    if (!tsv_file.is_open()) {
        cout << "Failed to open TSV file for writing." << endl;
        return r;
    }

    // Write TSV header
    tsv_file << "repertoire_id\tsequence_id\tduplicate_count\tjunction_aa\n";

    std::string repertoire_id = "rep1"; // same for all rows
    for (const auto& row : r) {
        tsv_file 
                 << repertoire_id << "\t"         // repertoire_id
                 << row["akc_id"].c_str() << "\t" // sequence_id
                 << "1\t"                         // duplicate_count
                 << row["junction_aa"].c_str()    // junction_aa
                 << "\n";
    }

    cout << "TSV file written successfully with " << r.size() << " rows." << endl;

    tx.commit();
    
    return r;
}

void run_compairr(void) {
    std::vector<std::string> argv_strs;
    argv_strs.push_back("filler for prog_name, doesn't matter but needs to be here");
    argv_strs.push_back("-m");                  // compute overlap matrix
    argv_strs.push_back("repertoire.tsv");      // filename 1
    argv_strs.push_back("repertoire.tsv");      // filename 2
    argv_strs.push_back("-g");                  // ignore genes
    argv_strs.push_back("-d");                  // differences
    argv_strs.push_back("1");
    argv_strs.push_back("-i");                  // indels
    argv_strs.push_back("-l");                  // log filename
    argv_strs.push_back("compairr.log");
    argv_strs.push_back("-o");                  // matrix output filename
    argv_strs.push_back("output_matrix.tsv");
    // argv_strs.push_back("--no-matrix");         // no matrix output
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
    db_read(d1, "repertoire.tsv", true, "1");
    
    db * d2 = db_create();
    // db_read(d2, "repertoire.tsv", true, "2");
    d2 = d1;

    compairr::computeOverlap(d1, d2, true);

    close_files();
}

int main(void) {
    query_and_write_tsv();
    run_compairr();

    return 0;
}
// time docker run -v $(pwd):/data compairr -m -g -t 100 -d 1 --distance /data/dedup_repertoire.tsv -p /data/pairs_2.tsv