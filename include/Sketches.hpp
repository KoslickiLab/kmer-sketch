#ifndef SKETCHES_HPP
#define SKETCHES_HPP

#include <cstdint>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <queue>
#include <string>
#include <stdexcept>
#include <cmath>
#include <algorithm>
#include <limits>
#include <iostream>
#include <sstream>

#include "Hash.hpp"

static constexpr uint64_t HASH_MAX = std::numeric_limits<uint64_t>::max();

struct SketchInfo {
    std::string algo;
    size_t kmer_size = 0;
    uint64_t seed = 42;
    std::string hash_function = "MurmurHash3_x64_128_low64";
};

class MaxGeomSample {
public:
    struct Entry {
        uint64_t h;
        uint64_t hprime;
        uint64_t freq;
    };

    MaxGeomSample(size_t b, size_t w=64, uint64_t seed=42)
        : b_(b), w_(w), seed_(seed) {
        if (w_ < 1 || w_ > 64) throw std::runtime_error("w must be in [1,64]");
        if (b_ == 0) throw std::runtime_error("b must be positive");
    }

    void add_hash(uint64_t h) {
        size_t i = zpl_plus_one(h);
        uint64_t hprime = tail_after_leftmost_one(h, i);
        auto& bucket = buckets_[i];
        auto& heap = heaps_[i];
        auto it = bucket.find(h);
        if (it != bucket.end()) {
            it->second.freq += 1;
            return;
        }
        if (bucket.size() < b_) {
            bucket.emplace(h, Entry{h, hprime, 1});
            heap.emplace(hprime, h);
        } else {
            if (!heap.empty() && hprime > heap.top().first) {
                bucket.emplace(h, Entry{h, hprime, 1});
                heap.emplace(hprime, h);
                evict_smallest(i);
            }
        }
    }

    double jaccard(const MaxGeomSample& other) const {
        if (w_ != other.w_ || b_ != other.b_) throw std::runtime_error("Incompatible MaxGeomSample for Jaccard");
        size_t union_size = 0, inter_size = 0;
        for (const auto& kv : buckets_) {
            size_t i = kv.first;
            auto it2 = other.buckets_.find(i);
            if (it2 == other.buckets_.end()) continue;
            std::unordered_set<uint64_t> self_hprimes;
            self_hprimes.reserve(kv.second.size()*2);
            for (const auto& e : kv.second) self_hprimes.insert(e.second.hprime);
            std::unordered_set<uint64_t> other_hprimes;
            other_hprimes.reserve(it2->second.size()*2);
            for (const auto& e : it2->second) other_hprimes.insert(e.second.hprime);
            std::vector<uint64_t> uni;
            uni.reserve(self_hprimes.size()+other_hprimes.size());
            for (auto x: self_hprimes) uni.push_back(x);
            for (auto x: other_hprimes) uni.push_back(x);
            std::sort(uni.begin(), uni.end(), std::greater<uint64_t>());
            uni.erase(std::unique(uni.begin(), uni.end()), uni.end());
            if (uni.size() > b_) uni.resize(b_);
            union_size += uni.size();
            std::unordered_set<uint64_t> union_set(uni.begin(), uni.end());
            size_t c = 0;
            for (auto x : self_hprimes) if (other_hprimes.count(x) && union_set.count(x)) ++c;
            inter_size += c;
        }
        if (union_size == 0) return 1.0;
        return double(inter_size) / double(union_size);
    }


    double containment_in(const MaxGeomSample& other) const {
        if (w_ != other.w_ || b_ != other.b_)
            throw std::runtime_error("Incompatible MaxGeomSample for Containment (w or b differ)");
        size_t inter_sum = 0, denom_sum = 0;

        for (const auto& kv : buckets_) {
            size_t i = kv.first;
            auto it2 = other.buckets_.find(i);
            if (it2 == other.buckets_.end()) continue;

            // Build sets of h' (dedup within bucket)
            std::unordered_set<uint64_t> self_hprimes; self_hprimes.reserve(kv.second.size()*2);
            for (const auto& e : kv.second) self_hprimes.insert(e.second.hprime);

            std::unordered_set<uint64_t> other_hprimes; other_hprimes.reserve(it2->second.size()*2);
            for (const auto& e : it2->second) other_hprimes.insert(e.second.hprime);

            // Capacity-limited union (descending, then unique, capped at k_)
            std::vector<uint64_t> uni; uni.reserve(self_hprimes.size() + other_hprimes.size());
            for (auto x: self_hprimes) uni.push_back(x);
            for (auto x: other_hprimes) uni.push_back(x);
            std::sort(uni.begin(), uni.end(), std::greater<uint64_t>());
            uni.erase(std::unique(uni.begin(), uni.end()), uni.end());
            if (uni.size() > b_) uni.resize(b_);

            std::unordered_set<uint64_t> allow(uni.begin(), uni.end());

            // Denominator = how many of self's kept elements survive capacity cap
            size_t denom_bucket = 0;
            for (auto x : self_hprimes) if (allow.count(x)) ++denom_bucket;

            // Intersection = overlap that also survives capacity cap
            size_t inter_bucket = 0;
            for (auto x : self_hprimes) if (allow.count(x) && other_hprimes.count(x)) ++inter_bucket;

            denom_sum += denom_bucket;
            inter_sum += inter_bucket;
        }
        if (denom_sum == 0) return 1.0; // both empty under caps ⇒ define containment as 1
        return double(inter_sum) / double(denom_sum);
    }


    double cosine(const MaxGeomSample& other) const {
        if (w_ != other.w_ || b_ != other.b_) throw std::runtime_error("Incompatible MaxGeomSample for Cosine");
        double dot = 0.0, n1 = 0.0, n2 = 0.0;
        for (const auto& kv : buckets_) {
            size_t i = kv.first;
            auto it2 = other.buckets_.find(i);
            if (it2 == other.buckets_.end()) continue;
            std::unordered_map<uint64_t, uint64_t> f1, f2;
            for (const auto& e : kv.second) f1[e.second.hprime] = e.second.freq;
            for (const auto& e : it2->second) f2[e.second.hprime] = e.second.freq;
            std::vector<uint64_t> uni;
            uni.reserve(f1.size()+f2.size());
            for (auto& p : f1) uni.push_back(p.first);
            for (auto& p : f2) uni.push_back(p.first);
            std::sort(uni.begin(), uni.end(), std::greater<uint64_t>());
            uni.erase(std::unique(uni.begin(), uni.end()), uni.end());
            if (uni.size() > b_) uni.resize(b_);
            for (auto hpr : uni) {
                double a = double(f1.count(hpr) ? f1[hpr] : 0);
                double b = double(f2.count(hpr) ? f2[hpr] : 0);
                dot += a*b;
                n1 += a*a;
                n2 += b*b;
            }
        }
        if (n1 == 0.0 || n2 == 0.0) return 0.0;
        return dot / (std::sqrt(n1) * std::sqrt(n2));
    }

    size_t w() const { return w_; }
    size_t b() const { return b_; }
    uint64_t seed() const { return seed_; }
    const std::unordered_map<size_t, std::unordered_map<uint64_t, Entry>>& buckets() const { return buckets_; }

    void write(std::ostream& out, size_t kmer_size) const {
        out << "# sketch_version=1\n";
        out << "# algo=MaxGeom\n";
        out << "# kmer_size=" << kmer_size << "\n";
        out << "# hash_seed=" << seed_ << "\n";
        out << "# hash_function=MurmurHash3_x64_128_low64\n";
        out << "# params.b=" << b_ << "\n";
        out << "# params.w=" << w_ << "\n";
        out << "# fields: bucket_index,h,hprime,freq\n";
        for (const auto& kv : buckets_) {
            size_t i = kv.first;
            for (const auto& e : kv.second) {
                out << i << "," << e.second.h << "," << e.second.hprime << "," << e.second.freq << "\n";
            }
        }
    }

    static MaxGeomSample read(std::istream& in, size_t& kmer_size, uint64_t& seed_out) {
        std::string line;
        size_t b=0, w=64;
        while (std::getline(in, line)) {
            if (line.empty()) continue;
            if (line[0] != '#') break;
            auto pos = line.find('=');
            if (pos==std::string::npos) continue;
            std::string key = line.substr(2, pos-2);
            std::string val = line.substr(pos+1);
            if (key=="kmer_size") kmer_size = std::stoull(val);
            else if (key=="hash_seed") seed_out = std::stoull(val);
            else if (key=="params.b") b = std::stoull(val);
            else if (key=="params.w") w = std::stoull(val);
        }
        if (b==0) throw std::runtime_error("MaxGeomSample.read: missing params.b");
        MaxGeomSample s(b,w,seed_out);
        if (!line.empty() && line[0] != '#') {
            size_t i; uint64_t h, hp, f;
            std::replace(line.begin(), line.end(), ',', ' ');
            std::stringstream ls(line);
            if (!(ls >> i >> h >> hp >> f)) throw std::runtime_error("Malformed data line in MaxGeom sample");
            s.buckets_[i][h] = Entry{h,hp,f};
            s.heaps_[i].emplace(hp, h);
        }
        while (std::getline(in, line)) {
            if (line.empty() || line[0]=='#') continue;
            size_t i; uint64_t h, hp, f;
            std::replace(line.begin(), line.end(), ',', ' ');
            std::stringstream ls(line);
            if (!(ls >> i >> h >> hp >> f)) continue;
            s.buckets_[i][h] = Entry{h,hp,f};
            s.heaps_[i].emplace(hp, h);
        }
        return s;
    }

private:
    size_t zpl_plus_one(uint64_t h) const {
        if (w_ == 64) {
            if (h == 0) return w_;
            int lz = __builtin_clzll(h);
            return size_t(lz + 1);
        } else {
            uint64_t topw = (h >> (64 - w_)) & ((w_==64? ~0ULL : ((1ULL<<w_) - 1)));
            if (topw == 0) return w_;
            int lz_full = __builtin_clzll(topw);
            int lz = lz_full - (64 - (int)w_);
            return size_t(lz + 1);
        }
    }

    uint64_t tail_after_leftmost_one(uint64_t h, size_t i) const {
        size_t lower_bits = 64 - w_;
        uint64_t topw = (h >> (64 - w_)) & (w_==64? ~0ULL : ((1ULL<<w_) - 1));
        size_t rem_len = (w_ > i ? w_ - i : 0);
        uint64_t rem_top = (rem_len==64? topw : (topw & ((rem_len==64? ~0ULL : ((1ULL<<rem_len) - 1)))));
        uint64_t low = (lower_bits==64? 0 : (h & ((lower_bits==64? ~0ULL : ((1ULL<<lower_bits) - 1)))));
        return (rem_top << lower_bits) | low;
    }

    void evict_smallest(size_t i) {
        auto& bucket = buckets_[i];
        auto& heap = heaps_[i];
        while (bucket.size() > b_) {
            if (heap.empty()) break;
            auto [hp, h] = heap.top();
            heap.pop();
            auto it = bucket.find(h);
            if (it != bucket.end() && it->second.hprime == hp) {
                bucket.erase(it);
                break;
            }
        }
    }

    size_t b_;
    size_t w_;
    uint64_t seed_;
    std::unordered_map<size_t, std::unordered_map<uint64_t, Entry>> buckets_;
    struct MinCmp { bool operator()(const std::pair<uint64_t,uint64_t>& a,
                                    const std::pair<uint64_t,uint64_t>& b) const {
                        return a.first > b.first;
                    } };
    std::unordered_map<size_t, std::priority_queue<std::pair<uint64_t,uint64_t>,
                                                   std::vector<std::pair<uint64_t,uint64_t> >,
                                                   MinCmp>> heaps_;
};

class AlphaMaxGeomSample {
public:
    struct Entry { uint64_t h; uint64_t hprime; uint64_t freq; };

    AlphaMaxGeomSample(double alpha, size_t w=64, uint64_t seed=42)
        : alpha_(alpha), w_(w), seed_(seed) {
        if (!(alpha > 0.0 && alpha < 1.0)) throw std::runtime_error("alpha must be in (0,1)");
        if (w_ < 1 || w_ > 64) throw std::runtime_error("w must be in [1,64]");
        b_sizes_.resize(w_+1, 1);
        double beta = alpha_/(1.0-alpha_);
        for (size_t i=0; i<=w_; ++i) {
            double val = std::pow(2.0, beta * (double)i);
            // take ceiling of val
            size_t b = (size_t)std::ceil(val);
            if (b < 1) b = 1;
            b_sizes_[i] = b;
        }
    }

    void add_hash(uint64_t h) {
        size_t i = zpl_plus_one(h);
        uint64_t hprime = tail_after_leftmost_one(h, i);
        auto& bucket = buckets_[i];
        auto& heap = heaps_[i];
        auto it = bucket.find(h);
        if (it != bucket.end()) {
            it->second.freq += 1;
            return;
        }
        if (bucket.size() < b_sizes_[i]) {
            bucket.emplace(h, Entry{h,hprime,1});
            heap.emplace(hprime, h);
        } else {
            if (!heap.empty() && hprime > heap.top().first) {
                bucket.emplace(h, Entry{h,hprime,1});
                heap.emplace(hprime, h);
                evict_smallest(i);
            }
        }
    }

    double jaccard(const AlphaMaxGeomSample& other) const {
        if (w_ != other.w_ || alpha_ != other.alpha_) throw std::runtime_error("Incompatible AlphaMaxGeomSample for Jaccard");
        size_t union_size = 0, inter_size = 0;
        for (const auto& kv : buckets_) {
            size_t i = kv.first;
            auto it2 = other.buckets_.find(i);
            if (it2 == other.buckets_.end()) continue;
            std::unordered_set<uint64_t> self_hprimes;
            self_hprimes.reserve(kv.second.size()*2);
            for (const auto& e : kv.second) self_hprimes.insert(e.second.hprime);
            std::unordered_set<uint64_t> other_hprimes;
            other_hprimes.reserve(it2->second.size()*2);
            for (const auto& e : it2->second) other_hprimes.insert(e.second.hprime);
            std::vector<uint64_t> uni;
            uni.reserve(self_hprimes.size()+other_hprimes.size());
            for (auto x: self_hprimes) uni.push_back(x);
            for (auto x: other_hprimes) uni.push_back(x);
            std::sort(uni.begin(), uni.end(), std::greater<uint64_t>());
            uni.erase(std::unique(uni.begin(), uni.end()), uni.end());
            size_t cap = b_sizes_[i];
            if (uni.size() > cap) uni.resize(cap);
            union_size += uni.size();
            std::unordered_set<uint64_t> union_set(uni.begin(), uni.end());
            size_t c = 0;
            for (auto x : self_hprimes) if (other_hprimes.count(x) && union_set.count(x)) ++c;
            inter_size += c;
        }
        if (union_size == 0) return 1.0;
        return double(inter_size) / double(union_size);
    }

    double containment_in(const AlphaMaxGeomSample& other) const {
        if (w_ != other.w_ || alpha_ != other.alpha_)
            throw std::runtime_error("Incompatible AlphaMaxGeomSample for Containment (w or alpha differ)");
        size_t inter_sum = 0, denom_sum = 0;

        for (const auto& kv : buckets_) {
            size_t i = kv.first;
            auto it2 = other.buckets_.find(i);
            if (it2 == other.buckets_.end()) continue;

            std::unordered_set<uint64_t> self_hprimes; self_hprimes.reserve(kv.second.size()*2);
            for (const auto& e : kv.second) self_hprimes.insert(e.second.hprime);

            std::unordered_set<uint64_t> other_hprimes; other_hprimes.reserve(it2->second.size()*2);
            for (const auto& e : it2->second) other_hprimes.insert(e.second.hprime);

            // Capacity-limited union (cap depends on bucket i)
            std::vector<uint64_t> uni; uni.reserve(self_hprimes.size() + other_hprimes.size());
            for (auto x: self_hprimes) uni.push_back(x);
            for (auto x: other_hprimes) uni.push_back(x);
            std::sort(uni.begin(), uni.end(), std::greater<uint64_t>());
            uni.erase(std::unique(uni.begin(), uni.end()), uni.end());
            size_t cap = b_sizes_[i];
            if (uni.size() > cap) uni.resize(cap);

            std::unordered_set<uint64_t> allow(uni.begin(), uni.end());

            size_t denom_bucket = 0;
            for (auto x : self_hprimes) if (allow.count(x)) ++denom_bucket;

            size_t inter_bucket = 0;
            for (auto x : self_hprimes) if (allow.count(x) && other_hprimes.count(x)) ++inter_bucket;

            denom_sum += denom_bucket;
            inter_sum += inter_bucket;
        }
        if (denom_sum == 0) return 1.0;
        return double(inter_sum) / double(denom_sum);
    }

    double cosine(const AlphaMaxGeomSample& other) const {
        if (w_ != other.w_ || alpha_ != other.alpha_) throw std::runtime_error("Incompatible AlphaMaxGeomSample for Cosine");
        double dot = 0.0, n1 = 0.0, n2 = 0.0;
        for (const auto& kv : buckets_) {
            size_t i = kv.first;
            auto it2 = other.buckets_.find(i);
            if (it2 == other.buckets_.end()) continue;
            std::unordered_map<uint64_t,uint64_t> f1, f2;
            for (const auto& e : kv.second) f1[e.second.hprime] = e.second.freq;
            for (const auto& e : it2->second) f2[e.second.hprime] = e.second.freq;
            std::vector<uint64_t> uni;
            uni.reserve(f1.size()+f2.size());
            for (auto& p : f1) uni.push_back(p.first);
            for (auto& p : f2) uni.push_back(p.first);
            std::sort(uni.begin(), uni.end(), std::greater<uint64_t>());
            uni.erase(std::unique(uni.begin(), uni.end()), uni.end());
            size_t cap = b_sizes_[i];
            if (uni.size() > cap) uni.resize(cap);
            for (auto hpr : uni) {
                double a = double(f1.count(hpr) ? f1[hpr] : 0);
                double b = double(f2.count(hpr) ? f2[hpr] : 0);
                dot += a*b;
                n1 += a*a;
                n2 += b*b;
            }
        }
        if (n1 == 0.0 || n2 == 0.0) return 0.0;
        return dot / (std::sqrt(n1) * std::sqrt(n2));
    }

    size_t w() const { return w_; }
    double alpha() const { return alpha_; }
    uint64_t seed() const { return seed_; }
    const std::unordered_map<size_t, std::unordered_map<uint64_t, Entry>>& buckets() const { return buckets_; }

    void write(std::ostream& out, size_t kmer_size) const {
        out << "# sketch_version=1\n";
        out << "# algo=AlphaMaxGeom\n";
        out << "# kmer_size=" << kmer_size << "\n";
        out << "# hash_seed=" << seed_ << "\n";
        out << "# hash_function=MurmurHash3_x64_128_low64\n";
        out << "# params.alpha=" << alpha_ << "\n";
        out << "# params.w=" << w_ << "\n";
        out << "# fields: bucket_index,h,hprime,freq\n";
        for (const auto& kv : buckets_) {
            size_t i = kv.first;
            for (const auto& e : kv.second) {
                out << i << "," << e.second.h << "," << e.second.hprime << "," << e.second.freq << "\n";
            }
        }
    }

    static AlphaMaxGeomSample read(std::istream& in, size_t& kmer_size, uint64_t& seed_out, double& alpha_out) {
        std::string line;
        size_t w=64;
        while (std::getline(in, line)) {
            if (line.empty()) continue;
            if (line[0] != '#') break;
            auto pos = line.find('=');
            if (pos==std::string::npos) continue;
            std::string key = line.substr(2, pos-2);
            std::string val = line.substr(pos+1);
            if (key=="kmer_size") kmer_size = std::stoull(val);
            else if (key=="hash_seed") seed_out = std::stoull(val);
            else if (key=="params.alpha") alpha_out = std::stod(val);
            else if (key=="params.w") w = std::stoull(val);
        }
        AlphaMaxGeomSample s(alpha_out, w, seed_out);
        if (!line.empty() && line[0] != '#') {
            size_t i; uint64_t h, hp, f;
            std::replace(line.begin(), line.end(), ',', ' ');
            std::stringstream ls(line);
            if (!(ls >> i >> h >> hp >> f)) throw std::runtime_error("Malformed data line in AlphaMaxGeom sample");
            s.buckets_[i][h] = Entry{h,hp,f};
            s.heaps_[i].emplace(hp, h);
        }
        while (std::getline(in, line)) {
            if (line.empty() || line[0]=='#') continue;
            size_t i; uint64_t h, hp, f;
            std::replace(line.begin(), line.end(), ',', ' ');
            std::stringstream ls(line);
            if (!(ls >> i >> h >> hp >> f)) continue;
            s.buckets_[i][h] = Entry{h,hp,f};
            s.heaps_[i].emplace(hp, h);
        }
        return s;
    }

private:
    size_t zpl_plus_one(uint64_t h) const {
        if (w_ == 64) {
            if (h == 0) return w_;
            int lz = __builtin_clzll(h);
            return size_t(lz + 1);
        } else {
            uint64_t topw = (h >> (64 - w_)) & ((w_==64? ~0ULL : ((1ULL<<w_) - 1)));
            if (topw == 0) return w_;
            int lz_full = __builtin_clzll(topw);
            int lz = lz_full - (64 - (int)w_);
            return size_t(lz + 1);
        }
    }

    uint64_t tail_after_leftmost_one(uint64_t h, size_t i) const {
        size_t lower_bits = 64 - w_;
        uint64_t topw = (h >> (64 - w_)) & (w_==64? ~0ULL : ((1ULL<<w_) - 1));
        size_t rem_len = (w_ > i ? w_ - i : 0);
        uint64_t rem_top = (rem_len==64? topw : (topw & ((rem_len==64? ~0ULL : ((1ULL<<rem_len) - 1)))));
        uint64_t low = (lower_bits==64? 0 : (h & ((lower_bits==64? ~0ULL : ((1ULL<<lower_bits) - 1)))));
        return (rem_top << lower_bits) | low;
    }

    void evict_smallest(size_t i) {
        auto& bucket = buckets_[i];
        auto& heap = heaps_[i];
        while (bucket.size() > b_sizes_[i]) {
            if (heap.empty()) break;
            auto [hp, h] = heap.top();
            heap.pop();
            auto it = bucket.find(h);
            if (it != bucket.end() && it->second.hprime == hp) {
                bucket.erase(it);
                break;
            }
        }
    }

    double alpha_;
    size_t w_;
    uint64_t seed_;
    std::vector<size_t> b_sizes_;
    std::unordered_map<size_t, std::unordered_map<uint64_t, Entry>> buckets_;
    struct MinCmp { bool operator()(const std::pair<uint64_t,uint64_t>& a,
                                    const std::pair<uint64_t,uint64_t>& b) const {
                        return a.first > b.first;
                    } };
    std::unordered_map<size_t, std::priority_queue<std::pair<uint64_t,uint64_t>,
                                                   std::vector<std::pair<uint64_t,uint64_t> >,
                                                   MinCmp>> heaps_;
};

class FracMinHash {
public:
    explicit FracMinHash(double scale, uint64_t seed=42)
        : scale_(scale), seed_(seed) {
        if (!(scale_ > 0.0 && scale_ <= 1.0)) throw std::runtime_error("scale must be in (0,1]");
        threshold_ = (uint64_t)(scale_ * double(HASH_MAX));
    }

    void add_hash(uint64_t h) {
        if (h <= threshold_) hashes_.insert(h);
    }

    size_t size() const { return hashes_.size(); }

    static double jaccard(const FracMinHash& a, const FracMinHash& b) {
        size_t inter = 0;
        if (a.hashes_.size() < b.hashes_.size()) {
            for (auto& x : a.hashes_) if (b.hashes_.count(x)) ++inter;
        } else {
            for (auto& x : b.hashes_) if (a.hashes_.count(x)) ++inter;
        }
        size_t uni = a.hashes_.size() + b.hashes_.size() - inter;
        if (uni == 0) return 1.0;
        return double(inter) / double(uni);
    }
    // containment_in(self, other): fraction of self's sketch elements found in other.
    // Estimates C(self ⊆ other) = |self ∩ other| / |self|.
    // Returns 1.0 when self is empty (vacuously contained).
    static double containment_in(const FracMinHash& self, const FracMinHash& other) {
        if (self.hashes_.empty()) return 1.0;
        size_t inter = 0;
        for (auto x : self.hashes_) if (other.hashes_.count(x)) ++inter;
        return double(inter) / double(self.hashes_.size());
    }

    static double cosine(const FracMinHash& a, const FracMinHash& b) {
        size_t inter = 0;
        if (a.hashes_.size() < b.hashes_.size()) {
            for (auto& x : a.hashes_) if (b.hashes_.count(x)) ++inter;
        } else {
            for (auto& x : b.hashes_) if (a.hashes_.count(x)) ++inter;
        }
        if (a.hashes_.empty() || b.hashes_.empty()) return 1.0;
        return double(inter) / std::sqrt(double(a.hashes_.size()) * double(b.hashes_.size()));
    }

    void write(std::ostream& out, size_t kmer_size) const {
        out << "# sketch_version=1\n";
        out << "# algo=FracMinHash\n";
        out << "# kmer_size=" << kmer_size << "\n";
        out << "# hash_seed=" << seed_ << "\n";
        out << "# hash_function=MurmurHash3_x64_128_low64\n";
        out << "# params.scale=" << scale_ << "\n";
        out << "# params.threshold=" << threshold_ << "\n";
        out << "# fields: hash\n";
        // sort the hashes for output
        std::vector<uint64_t> v(hashes_.begin(), hashes_.end());
        std::sort(v.begin(), v.end());
        for (auto h : v) out << h << "\n";
    }

    static FracMinHash read(std::istream& in, size_t& kmer_size, uint64_t& seed_out, double& scale_out) {
        std::string line;
        uint64_t threshold=0;
        while (std::getline(in, line)) {
            if (line.empty()) continue;
            if (line[0] != '#') break;
            auto pos = line.find('=');
            if (pos==std::string::npos) continue;
            std::string key = line.substr(2, pos-2);
            std::string val = line.substr(pos+1);
            if (key=="kmer_size") kmer_size = std::stoull(val);
            else if (key=="hash_seed") seed_out = std::stoull(val);
            else if (key=="params.scale") scale_out = std::stod(val);
            else if (key=="params.threshold") threshold = std::stoull(val);
        }
        FracMinHash s(scale_out, seed_out);
        s.threshold_ = threshold;
        if (!line.empty() && line[0] != '#') {
            uint64_t h = 0;
            std::stringstream ls(line);
            if (ls >> h) s.hashes_.insert(h);
        }
        while (std::getline(in, line)) {
            if (line.empty() || line[0]=='#') continue;
            uint64_t h = 0; std::stringstream ls(line); if (ls >> h) s.hashes_.insert(h);
        }
        return s;
    }

    const std::unordered_set<uint64_t>& hashes() const { return hashes_; }
    double scale() const { return scale_; }
    uint64_t threshold() const { return threshold_; }
    uint64_t seed() const { return seed_; }

private:
    double scale_;
    uint64_t seed_;
    uint64_t threshold_;
    std::unordered_set<uint64_t> hashes_;
};

class BottomK {
public:
    explicit BottomK(size_t k, uint64_t seed=42) : k_(k), seed_(seed) {
        if (k_ == 0) throw std::runtime_error("k must be positive");
    }
    void add_hash(uint64_t h) {
        if (set_.count(h)) return;
        if (set_.size() < k_) {
            set_.insert(h);
            heap_.push(h);
        } else if (!heap_.empty() && h < heap_.top()) {
            uint64_t old = heap_.top(); heap_.pop();
            set_.erase(old);
            set_.insert(h);
            heap_.push(h);
        }
    }
    size_t size() const { return set_.size(); }
    size_t k() const { return k_; }
    const std::unordered_set<uint64_t>& hashes() const { return set_; }

    static double jaccard(const BottomK& a, const BottomK& b) {
        size_t inter = 0;
        if (a.size() < b.size()) { for (auto x: a.set_) if (b.set_.count(x)) ++inter; }
        else { for (auto x: b.set_) if (a.set_.count(x)) ++inter; }
        size_t uni = a.size() + b.size() - inter;
        if (uni == 0) return 1.0;
        return double(inter)/double(uni);
    }
    // containment_in(self, other): estimates C(self ⊆ other) = |self ∩ other| / |self|.
    // Returns 1.0 when self is empty (vacuously contained).
    //
    // For bottom-k sketches the naive estimator |self ∩ other| / |self| is biased
    // when the two genomes differ in size: if other is larger, its sketch's k-th
    // minimum (τ) is smaller than self's, so many of self's sketch elements lie
    // above τ and are invisible to other's sketch — causing the naive estimator to
    // undercount containment.
    //
    // We correct for this by restricting the denominator to only those elements of
    // self that fall within other's sketch window (hash ≤ τ = max(other)):
    //
    //   C(self ⊆ other) ≈ |{x ∈ self : x ∈ other AND x ≤ τ}|
    //                      ─────────────────────────────────────
    //                         |{x ∈ self : x ≤ τ}|
    //
    // This is unbiased when both genomes have ≥ k k-mers.  When other has fewer
    // than k k-mers its sketch is exhaustive (τ = max of other's hashes), so the
    // estimator reduces to the naive form, which is exact in that case.
    // Returns 0.0 when none of self's sketch elements fall within other's window.
    static double containment_in(const BottomK& self, const BottomK& other) {
        if (self.set_.empty()) return 1.0;
        if (other.set_.empty()) return 0.0;
        uint64_t tau = *std::max_element(other.set_.begin(), other.set_.end());
        size_t denom = 0, inter = 0;
        for (auto x : self.set_) {
            if (x <= tau) {
                ++denom;
                if (other.set_.count(x)) ++inter;
            }
        }
        if (denom == 0) return 0.0;
        return double(inter) / double(denom);
    }

    static double cosine(const BottomK& a, const BottomK& b) {
        size_t inter = 0;
        if (a.size() < b.size()) { for (auto x: a.set_) if (b.set_.count(x)) ++inter; }
        else { for (auto x: b.set_) if (a.set_.count(x)) ++inter; }
        if (a.size()==0 || b.size()==0) return 1.0;
        return double(inter) / std::sqrt(double(a.size())*double(b.size()));
    }

    void write(std::ostream& out, size_t kmer_size) const {
        out << "# sketch_version=1\n";
        out << "# algo=BottomK\n";
        out << "# kmer_size=" << kmer_size << "\n";
        out << "# hash_seed=" << seed_ << "\n";
        out << "# hash_function=MurmurHash3_x64_128_low64\n";
        out << "# params.k=" << k_ << "\n";
        out << "# fields: hash\n";
        std::vector<uint64_t> v(set_.begin(), set_.end());
        std::sort(v.begin(), v.end());
        for (auto h : v) out << h << "\n";
    }

    static BottomK read(std::istream& in, size_t& kmer_size, uint64_t& seed_out, size_t& k_out) {
        std::string line;
        while (std::getline(in, line)) {
            if (line.empty()) continue;
            if (line[0] != '#') break;
            auto pos = line.find('=');
            if (pos==std::string::npos) continue;
            std::string key = line.substr(2, pos-2);
            std::string val = line.substr(pos+1);
            if (key=="kmer_size") kmer_size = std::stoull(val);
            else if (key=="hash_seed") seed_out = std::stoull(val);
            else if (key=="params.k") k_out = std::stoull(val);
        }
        BottomK s(k_out, seed_out);
        if (!line.empty() && line[0] != '#') {
            uint64_t h = 0;
            std::stringstream ls(line);
            if (ls >> h) { s.set_.insert(h); s.heap_.push(h); }
        }
        while (std::getline(in, line)) {
            if (line.empty() || line[0]=='#') continue;
            uint64_t h = 0; std::stringstream ls(line); if (ls >> h) { s.set_.insert(h); s.heap_.push(h); }
        }
        return s;
    }

private:
    size_t k_;
    uint64_t seed_;
    std::unordered_set<uint64_t> set_;
    std::priority_queue<uint64_t> heap_;
};

class MinHash {
public:
    MinHash(size_t num_perm, uint64_t seed=42) : num_perm_(num_perm), seed_(seed) {
        if (num_perm_ == 0) throw std::runtime_error("num_perm must be positive");
        mins_.assign(num_perm_, HASH_MAX);
        hashutil::SplitMix64 gen(seed_);
        a_.resize(num_perm_); b_.resize(num_perm_);
        for (size_t i=0;i<num_perm_;++i) {
            uint64_t ai = gen.next() | 1ULL;
            uint64_t bi = gen.next();
            a_[i] = ai;
            b_[i] = bi;
        }
    }

    void add_hash(uint64_t h) {
        for (size_t i=0;i<num_perm_;++i) {
            uint64_t x = hashutil::mix64(a_[i] * h + b_[i]);
            if (x < mins_[i]) mins_[i] = x;
        }
    }

    size_t num_perm() const { return num_perm_; }
    const std::vector<uint64_t>& mins() const { return mins_; }
    uint64_t seed() const { return seed_; }
    size_t size() const { return num_perm_; }

    static double jaccard(const MinHash& a, const MinHash& b) {
        if (a.num_perm_ != b.num_perm_ || a.seed_ != b.seed_) throw std::runtime_error("Incompatible MinHash for comparison");
        size_t matches = 0;
        for (size_t i=0;i<a.num_perm_;++i) if (a.mins_[i] == b.mins_[i]) ++matches;
        return double(matches) / double(a.num_perm_);
    }
    static double cosine(const MinHash& a, const MinHash& b) {
        return jaccard(a,b);
    }

    void write(std::ostream& out, size_t kmer_size) const {
        out << "# sketch_version=1\n";
        out << "# algo=MinHash\n";
        out << "# kmer_size=" << kmer_size << "\n";
        out << "# hash_seed=" << seed_ << "\n";
        out << "# hash_function=MurmurHash3_x64_128_low64\n";
        out << "# params.num_perm=" << num_perm_ << "\n";
        out << "# fields: index,hash\n";
        for (size_t i=0;i<mins_.size();++i) out << i << "," << mins_[i] << "\n";
    }

    static MinHash read(std::istream& in, size_t& kmer_size, uint64_t& seed_out, size_t& num_perm_out) {
        std::string line;
        while (std::getline(in, line)) {
            if (line.empty()) continue;
            if (line[0] != '#') break;
            auto pos = line.find('=');
            if (pos==std::string::npos) continue;
            std::string key = line.substr(2, pos-2);
            std::string val = line.substr(pos+1);
            if (key=="kmer_size") kmer_size = std::stoull(val);
            else if (key=="hash_seed") seed_out = std::stoull(val);
            else if (key=="params.num_perm") num_perm_out = std::stoull(val);
        }
        MinHash s(num_perm_out, seed_out);
        s.mins_.assign(num_perm_out, HASH_MAX);
        if (!line.empty() && line[0] != '#') {
            size_t idx; uint64_t h;
            std::replace(line.begin(), line.end(), ',', ' ');
            std::stringstream ls(line);
            if (ls >> idx >> h) {
                if (idx < s.mins_.size()) s.mins_[idx] = h;
            }
        }
        while (std::getline(in, line)) {
            if (line.empty() || line[0]=='#') continue;
            std::replace(line.begin(), line.end(), ',', ' ');
            std::stringstream ls(line);
            size_t idx; uint64_t h;
            if (ls >> idx >> h) {
                if (idx < s.mins_.size()) s.mins_[idx] = h;
            }
        }
        return s;
    }

private:
    size_t num_perm_;
    uint64_t seed_;
    std::vector<uint64_t> a_, b_;
    std::vector<uint64_t> mins_;
};

#endif // SKETCHES_HPP
