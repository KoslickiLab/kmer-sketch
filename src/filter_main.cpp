#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <stdexcept>
// (optional) #include <cerrno>

#include "SketchIO.hpp"

static void usage() {
    std::cerr << "Usage:\n"
              << "  filter --query Q.sketch "
              << "--refs R1.sketch [R2.sketch ...] "
              << "[--refs @list.txt] [--refs-filelist list.txt] "
              << "--metric {jaccard|containment|cosine} --threshold X [--output OUT.tsv]\n";
}

static std::string get_arg(std::vector<std::string>& args, const std::string& key, const std::string& def="") {
    for (size_t i=0;i<args.size();++i) if (args[i]==key && i+1<args.size()) return args[i+1];
    return def;
}

// --- NEW: helpers for filelists ---------------------------------------------

static std::string trim_ws(const std::string& s) {
    const char* ws = " \t\r\n";
    const auto a = s.find_first_not_of(ws);
    if (a == std::string::npos) return "";
    const auto b = s.find_last_not_of(ws);
    return s.substr(a, b - a + 1);
}

static void add_refs_from_filelist(std::vector<std::string>& out, const std::string& listpath) {
    std::ifstream fin(listpath);
    if (!fin) throw std::runtime_error("cannot open refs list: " + listpath);

    // base directory of list file for resolving relative entries
    std::string basedir;
    size_t pos = listpath.find_last_of("/\\");
    if (pos != std::string::npos) basedir = listpath.substr(0, pos + 1);

    std::string line;
    while (std::getline(fin, line)) {
        line = trim_ws(line);
        if (line.empty() || line[0] == '#') continue;

        // unquote if "..." or '...'
        if (line.size() >= 2 &&
            ((line.front()=='"' && line.back()=='"') ||
             (line.front()=='\'' && line.back()=='\''))) {
            line = line.substr(1, line.size()-2);
        }

        // if relative path, resolve relative to the list file directory
        bool is_abs_unix = (!line.empty() && (line[0] == '/'));
        bool is_abs_win  = (line.size() > 1 && line[1] == ':'); // C:\...
        if (!basedir.empty() && !is_abs_unix && !is_abs_win) {
            out.push_back(basedir + line);
        } else {
            out.push_back(line);
        }
    }
}

static void add_ref_token(std::vector<std::string>& out, const std::string& tok) {
    if (!tok.empty() && tok[0] == '@') {
        // response file expansion
        add_refs_from_filelist(out, tok.substr(1));
    } else {
        out.push_back(tok);
    }
}

// ---------------------------------------------------------------------------

int main(int argc, char** argv) {
    if (argc < 2) { usage(); return 1; }
    std::vector<std::string> args(argv+1, argv+argc);

    std::string qpath     = get_arg(args, "--query");
    std::string metric    = get_arg(args, "--metric", "jaccard");
    std::string outpath   = get_arg(args, "--output", "");
    std::string reflist   = get_arg(args, "--refs-filelist", ""); // NEW
    double threshold      = std::stod(get_arg(args, "--threshold", "0.0"));

    // Collect refs:
    std::vector<std::string> refs;

    // 1) Everything after `--refs` (until next flag), supporting @file expansion
    bool after_refs = false;
    for (size_t i=0;i<args.size();++i) {
        if (args[i] == "--refs") { after_refs = true; continue; }
        if (after_refs) {
            if (!args[i].empty() && args[i].rfind("--",0)!=0) {
                add_ref_token(refs, args[i]); // handles @list.txt too
            } else {
                after_refs = false;
            }
        }
    }

    // 2) Optional dedicated file list: --refs-filelist list.txt
    if (!reflist.empty()) {
        try { add_refs_from_filelist(refs, reflist); }
        catch (const std::exception& e) {
            std::cerr << "Error reading --refs-filelist '" << reflist << "': " << e.what() << "\n";
            return 1;
        }
    }

    // de-dup
    std::sort(refs.begin(), refs.end());
    refs.erase(std::unique(refs.begin(), refs.end()), refs.end());

    if (qpath.empty() || refs.empty()) { usage(); return 1; }

    // Load query
    VariantSketch q = load_sketch(qpath);

    // Load refs and compute scores
    struct Row { std::string path; double score; size_t inter; size_t size1; size_t size2; size_t uni; };
    std::vector<Row> rows;

    for (const auto& rp : refs) {
        if (rp == qpath) continue;
        try {
            VariantSketch r = load_sketch(rp);
            std::string why;
            if (!compatible(q, r, why)) {
                std::cerr << "Skipping incompatible reference '" << rp << "': " << why << "\n";
                continue;
            }
            double s = 0.0; size_t inter=0, uni=0, size1=0, size2=0;
            if (q.maxgeom && r.maxgeom) {
                if (metric=="cosine") s = q.maxgeom->cosine(*r.maxgeom);
                else if (metric=="containment") s = r.maxgeom->containment_in(*q.maxgeom);
                else s = q.maxgeom->jaccard(*r.maxgeom);

                for (auto& kv: q.maxgeom->buckets()) size1 += kv.second.size();
                for (auto& kv: r.maxgeom->buckets()) size2 += kv.second.size();
            } else if (q.alphamaxgeom && r.alphamaxgeom) {
                if (metric=="cosine") s = q.alphamaxgeom->cosine(*r.alphamaxgeom);
                else if (metric=="containment") s = r.alphamaxgeom->containment_in(*q.alphamaxgeom);
                else s = q.alphamaxgeom->jaccard(*r.alphamaxgeom);
                for (auto& kv: q.alphamaxgeom->buckets()) size1 += kv.second.size();
                for (auto& kv: r.alphamaxgeom->buckets()) size2 += kv.second.size();
            } else if (q.fracmh && r.fracmh) {
                if (metric=="cosine") s = FracMinHash::cosine(*q.fracmh, *r.fracmh);
                else if (metric=="containment") s = FracMinHash::containment_in(*r.fracmh, *q.fracmh);
                else s = FracMinHash::jaccard(*q.fracmh, *r.fracmh);
                size1 = q.fracmh->hashes().size(); size2 = r.fracmh->hashes().size();
                for (auto x: q.fracmh->hashes()) if (r.fracmh->hashes().count(x)) ++inter;
                uni = size1 + size2 - inter;
            } else if (q.minhash && r.minhash) {
                s = (metric=="cosine") ? MinHash::cosine(*q.minhash,*r.minhash) : MinHash::jaccard(*q.minhash,*r.minhash);
                size1 = q.minhash->num_perm(); size2 = r.minhash->num_perm();
                for (size_t i=0;i<q.minhash->num_perm();++i) if (q.minhash->mins()[i]==r.minhash->mins()[i]) ++inter;
                uni = q.minhash->num_perm();
            } else if (q.bottomk && r.bottomk) {
                if (metric=="cosine") s = BottomK::cosine(*q.bottomk, *r.bottomk);
                else if (metric=="containment") s = BottomK::containment_in(*r.bottomk, *q.bottomk);
                else s = BottomK::jaccard(*q.bottomk, *r.bottomk);
                size1 = q.bottomk->hashes().size(); size2 = r.bottomk->hashes().size();
                for (auto x: q.bottomk->hashes()) if (r.bottomk->hashes().count(x)) ++inter;
                uni = size1 + size2 - inter;
            } else {
                std::cerr << "Skipping '" << rp << "': internal type mismatch\n";
                continue;
            }
            if (s >= threshold) rows.push_back({rp, s, inter, size1, size2, uni});
        } catch (const std::exception& e) {
            std::cerr << "Error processing '" << rp << "': " << e.what() << "\n";
        }
    }

    std::sort(rows.begin(), rows.end(), [](const Row& a, const Row& b){ return a.score > b.score; });

    std::ostream* pout = &std::cout;
    std::ofstream fout;
    if (!outpath.empty()) {
        fout.open(outpath);
        if (!fout) { std::cerr << "Cannot open --output '" << outpath << "'\n"; return 2; }
        pout = &fout;
    }
    *pout << "reference\t" << metric << "_score\tintersection\tsize_query\tsize_ref\tunion\n";
    for (auto& r : rows) {
        *pout << r.path << "\t" << r.score << "\t" << r.inter << "\t" << r.size1 << "\t" << r.size2 << "\t" << r.uni << "\n";
    }
    return 0;
}
