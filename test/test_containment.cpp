#include <cassert>
#include <cmath>
#include <iostream>
#include <sstream>
#include "Sketches.hpp"

static void near(double actual, double expected) {
    assert(std::isfinite(actual) && std::abs(actual - expected) < 1e-12);
}

template <typename F> static void throws(F f) {
    bool caught = false;
    try { f(); } catch (const std::runtime_error&) { caught = true; }
    assert(caught);
}

template <typename Factory> static void exact_tests(Factory factory) {
    auto a = factory(), b = factory(), empty = factory();
    const uint64_t h1 = (1ULL << 63) + 11, h2 = (1ULL << 62) + 22;
    a.add_hash(h1); a.add_hash(h2); a.add_hash(h1); // repeated input is still a set
    b.add_hash(h1);
    near(a.containment_in(b), 0.5); // includes an A-only bucket
    near(b.containment_in(a), 1.0); // directional, not Jaccard
    near(a.containment_in(a), 1.0);
    near(a.containment_in(empty), 0.0);
    b = factory(); b.add_hash(1ULL << 61);
    near(a.containment_in(b), 0.0); // entirely disjoint buckets
}

template <typename Factory> static void tie_tests(Factory factory) {
    auto a = factory(), b = factory();
    // Distinct full hashes can have the same hprime in the terminal bucket.
    const uint64_t lo = 5, hi = (1ULL << 63) + lo;
    a.add_hash(lo); b.add_hash(hi);
    const double disjoint = a.containment_in(b);
    assert(std::isnan(disjoint) || disjoint == 0.0); // never a false match
    // Equal priorities are broken consistently, independent of insertion order.
    a.add_hash(hi); b.add_hash(lo); a.add_hash(hi); b.add_hash(hi);
    near(a.containment_in(b), 1.0);
    for (const auto& [i, bucket] : a.buckets())
        for (const auto& [h, e] : bucket) {
            assert(b.buckets().at(i).count(h));
            assert(e.freq == b.buckets().at(i).at(h).freq);
        }
}

template <typename S, typename Reader> static void roundtrip(S a, S b, Reader read) {
    for (uint64_t i = 0; i < 100; ++i) a.add_hash(hashutil::mix64(i));
    for (uint64_t i = 50; i < 400; ++i) b.add_hash(hashutil::mix64(i));
    std::stringstream sa, sb;
    a.write(sa, 31); b.write(sb, 31);
    const auto aa = read(sa), bb = read(sb);
    const double before = a.containment_in(b), after = aa.containment_in(bb);
    assert((std::isnan(before) && std::isnan(after)) || before == after);
}

int main() {
    exact_tests([] { return MaxGeomSample(10); });
    exact_tests([] { return AlphaMaxGeomSample(0.45); });
    exact_tests([] { return BottomK(10); });
    exact_tests([] { return FracMinHash(1.0); });
    tie_tests([] { return MaxGeomSample(1, 1); });
    tie_tests([] { return AlphaMaxGeomSample(0.45, 1); });

    MaxGeomSample ma(1), mb(1);
    ma.add_hash((1ULL << 63) + 1); mb.add_hash((1ULL << 63) + 2);
    assert(std::isnan(ma.containment_in(mb)));
    AlphaMaxGeomSample aa(0.45), ab(0.45);
    for (uint64_t i : {1, 2}) aa.add_hash((1ULL << 63) + i);
    for (uint64_t i : {3, 4}) ab.add_hash((1ULL << 63) + i);
    assert(std::isnan(aa.containment_in(ab)));
    BottomK ba(2), bb(2);
    ba.add_hash(100); ba.add_hash(101); bb.add_hash(1); bb.add_hash(2);
    assert(std::isnan(ba.containment_in(bb)));
    near(MaxGeomSample(1).containment_in(ma), 1.0);
    near(AlphaMaxGeomSample(0.45).containment_in(aa), 1.0);
    near(BottomK(2).containment_in(ba), 1.0);

    // A common threshold is essential when FracMinHash scales differ.
    FracMinHash fa(0.5), fb(0.25), empty(0.5), full(1.0);
    for (uint64_t h : {1ULL, 2ULL, 1ULL << 63}) fa.add_hash(h);
    fb.add_hash(1); fb.add_hash(3);
    near(fa.containment_in(fb), 0.5);
    assert(std::isnan(empty.containment_in(fb)));
    full.add_hash(HASH_MAX); full.add_hash(0);
    assert(full.size() == 2 && full.threshold() == HASH_MAX);
    // Large alpha previously caused out-of-range floating-to-integer casts.
    AlphaMaxGeomSample large_alpha(0.9);
    large_alpha.add_hash(0); large_alpha.add_hash(1);
    near(large_alpha.containment_in(large_alpha), 1.0);

    throws([&] { ma.containment_in(MaxGeomSample(1, 64, 7)); });
    throws([&] { ma.containment_in(MaxGeomSample(2)); });
    throws([&] { ma.containment_in(MaxGeomSample(1, 63)); });
    throws([&] { aa.containment_in(AlphaMaxGeomSample(0.45, 64, 7)); });
    throws([&] { aa.containment_in(AlphaMaxGeomSample(0.4)); });
    throws([&] { ba.containment_in(BottomK(2, 7)); });
    throws([&] { fa.containment_in(FracMinHash(0.5, 7)); });

    roundtrip(MaxGeomSample(10), MaxGeomSample(10), [](std::istream& in) {
        size_t k = 0; uint64_t seed = 0; return MaxGeomSample::read(in, k, seed);
    });
    roundtrip(AlphaMaxGeomSample(0.45), AlphaMaxGeomSample(0.45), [](std::istream& in) {
        size_t k = 0; uint64_t seed = 0; double a = 0;
        return AlphaMaxGeomSample::read(in, k, seed, a);
    });
    roundtrip(BottomK(100), BottomK(100), [](std::istream& in) {
        size_t k = 0, cap = 0; uint64_t seed = 0; return BottomK::read(in, k, seed, cap);
    });
    roundtrip(FracMinHash(0.1), FracMinHash(0.1), [](std::istream& in) {
        size_t k = 0; uint64_t seed = 0; double scale = 0;
        return FracMinHash::read(in, k, seed, scale);
    });
    std::cout << "Containment unit tests passed\n";
}
