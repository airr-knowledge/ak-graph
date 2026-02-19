#include <iostream>
#include <fstream>
#include <pqxx/pqxx>
#include <compairr/api.h>
#include <compairr/db.h>
#include <vector>
#include <string>

using std::cout;
using std::endl;
using std::ofstream;

struct mem_file {
    FILE *fp;
    char *mem;
};

void query_stream_and_read(struct db * d) {
    // optional to add funcionality for in future
    bool require_sequence_id = true;
    const char * default_repertoire_id = "1";
    
    // Setup data
    struct stat fs;
    
    int fd = -1;
    bool is_regular = false;
    uint64_t filesize = 0;
    uint64_t fileread = 0;
    
    if (! is_regular)
    fprintf(logfile, "Waiting for data from SQL query...\n");

    size_t line_alloc = 4096;
    char * line = (char *) xmalloc(line_alloc);
    uint64_t lineno = 0;
    ssize_t linelen = 0;

    // d->longest = 0;
    db_set_longest(d, 0);
    // d->shortest = UINT_MAX;
    db_set_shortest(d, UINT_MAX);
    // d->ignored_unknown = 0;
    db_set_ignored_unknown(d, 0);
    // d->ignored_empty = 0;
    db_set_ignored_empty(d, 0);
    
    int state = 0;
    progress_init("Streaming sequences:", filesize);
    
    // Connect to the database
    pqxx::connection cx{"postgresql://postgres:example@ak-db/airrkb_v1"};
    pqxx::work tx{cx};
    
    // leave off \n at end 
    std::string input = "repertoire_id\tsequence_id\tduplicate_count\tjunction_aa";
    std::string buffer = input;
    
    // size_t line_alloc = input.size() + 1;   // +1 for '\0'
    // char *line = (char *) xmalloc(line_alloc);
    std::memcpy(line, input.c_str(), line_alloc);    
    linelen = input.length();

    parse_airr_tsv_header(line, d, require_sequence_id);

    if (is_regular)
        progress_update(fileread);

    lineno++;

    std::string repertoire_id = "rep1"; // same for all rows

    for (auto [int_index, junction_aa] :
        tx.stream<std::string, std::string>(
            "SELECT int_index, junction_aa "
            "FROM unique_junctions"
        )
    ) {
        input  = repertoire_id + "\t";   // repertoire_id
        input += int_index + "\t";       // sequence_id
        input += "1\t";                  // duplicate_count
        input += junction_aa;            // junction_aa
        // input += "\n";                   // leave off

        buffer += input;
        std::memcpy(line, input.c_str(), line_alloc);

        // if the \n is left out then we don't need to use this
        // std::cout << line << std::endl;
        // linelen = input.length();
        // if ((linelen > 0) && (line[linelen-1] == '\n'))
        // {
        //     line[linelen-1] = 0;
        //     linelen--;
        // }

        // if ((linelen > 0) && (line[linelen-1] == '\r'))
        // {
        //     line[linelen-1] = 0;
        //     linelen--;
        // }

        lineno++;

        parse_airr_tsv_line(line,
                            lineno,
                            d,
                            require_sequence_id,
                            default_repertoire_id);
        
        // if (is_regular)
        //     progress_update(fileread);

        fileread += linelen;
    }

    tx.commit();

    // make mem_file struct
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

    // progress_done();

    if (line)
        xfree(line);
    line = nullptr;

    fclose(fp);

    // d->repertoire_count = d->repertoire_id_vector.size();
    db_set_repertoire_count(d, db_get_repertoire_id_vector(d).size());

    // if (d->ignored_unknown > 0)
    if (db_get_ignored_unknown(d) > 0)
        // fprintf(logfile, "%" PRIu64 " sequences with unknown symbols ignored.\n", d->ignored_unknown);
        fprintf(logfile, "%" PRIu64 " sequences with unknown symbols ignored.\n", db_get_ignored_unknown(d));

    // if (d->ignored_empty > 0)
    if (db_get_ignored_empty(d) > 0)
        // fprintf(logfile, "%" PRIu64 " empty sequences ignored.\n", d->ignored_empty);
        fprintf(logfile, "%" PRIu64 " empty sequences ignored.\n", db_get_ignored_empty(d));

    // if (d->sequences > 0)
    if (db_getsequencecount(d) > 0)
        {
        fprintf(logfile,
                "Repertoires:       %" PRIu64 "\n"
                "Sequences:         %" PRIu64 "\n"
                "Residues:          %" PRIu64 "\n"
                "Shortest:          %u\n"
                "Longest:           %u\n"
                "Average length:    %.1lf\n"
                "Total dupl. count: %" PRIu64 "\n",
                // d->repertoire_count,
                db_get_repertoire_count(d),
                // d->sequences,
                db_getsequencecount(d),
                // d->residues_count,
                db_getresiduescount(d),
                // d->shortest,
                db_get_shortest(d),
                // d->longest,
                db_get_longest(d),
                // 1.0 * d->residues_count / d->sequences,
                1.0 * db_getresiduescount(d) / db_getsequencecount(d),
                // d->total_duplicate_count);
                db_get_total_duplicate_count(d));
        }
    else
        {
        fprintf(logfile,
                "Repertoires:       %" PRIu64 "\n"
                "Sequences:         %" PRIu64 "\n"
                "Residues:          %" PRIu64 "\n"
                "Shortest:          -\n"
                "Longest:           -\n"
                "Average length:    -\n"
                "Total dupl. count: %" PRIu64 "\n",
                // d->repertoire_count,
                db_get_repertoire_count(d),
                // d->sequences,
                db_getsequencecount(d),
                // d->residues_count,
                db_getresiduescount(d),
                // d->total_duplicate_count);
                db_get_total_duplicate_count(d));
        }

    /* add sequence pointers to index table */

    // progress_init("Indexing:         ", d->sequences);
    progress_init("Indexing:         ", db_getsequencecount(d));
    // char * r = d->residues_p;
    char * r = db_get_residues_p(d);
    for(uint64_t i = 0; i < db_getsequencecount(d); i++)
        {
        // seqinfo_s * p = d->seqindex + i;
        // p->seq = r;
        // r += p->seqlen;    
        db_set_seqinfo_s(d, r, i);
        progress_update(i+1);
        }
    progress_done();
}

void free_mem_file(mem_file mf) {
    fclose(mf.fp);
    free(mf.mem);
}

void run_compairr(void) {
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
    argv_strs.push_back("-t");                  // threads
    argv_strs.push_back("64");

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

    query_stream_and_read(d1);

    db * d2 = db_create();
    d2 = d1;

    compairr::computeOverlap(d1, d2, true);

    close_files();
}

int main(void) {
    run_compairr();

    return 0;
}
// time docker run -v $(pwd):/data compairr -m -g -t 100 -d 1 --distance /data/dedup_repertoire.tsv -p /data/pairs_2.tsv