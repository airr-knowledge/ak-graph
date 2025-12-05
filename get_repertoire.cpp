#include <iostream>
#include <fstream>
#include <pqxx/pqxx>

using std::cout;
using std::endl;
using std::ofstream;

pqxx::result query_and_write_tsv() {
    // Connect to the database
    pqxx::connection cx{"postgresql://postgres:example@ak-db/airrkb_v1"};
    pqxx::work tx{cx};

    // Query akc_id and junction_aa
    pqxx::result r{tx.exec(
        "SELECT tcr.akc_id, c.junction_aa "
        "FROM \"TCellReceptor\" tcr "
        "INNER JOIN \"Chain\" c ON tcr.trb_chain = c.akc_id "
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
        tsv_file << repertoire_id << "\t"          // repertoire_id
                 << row["akc_id"].c_str() << "\t" // sequence_id
                 << "1\t"                          // duplicate_count
                 << row["junction_aa"].c_str()    // junction_aa
                 << "\n";
    }

    cout << "TSV file written successfully with " << r.size() << " rows." << endl;

    tx.commit();
    return r;
}

int main() {
    query_and_write_tsv();
    return 0;
}
// time docker run -v $(pwd):/data compairr -m -g -t 100 -d 1 --distance /data/dedup_repertoire.tsv -p /data/pairs_2.tsv