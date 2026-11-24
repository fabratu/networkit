#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <memory>
#include <thread>

#if defined(__APPLE__)
#include <mach/mach.h>
#elif defined(__linux__)
#include <unistd.h>
#endif

#include <gtest/gtest.h>

#include <networkit/auxiliary/Timer.hpp>
#include <networkit/centrality/HyperKatzCentrality.hpp>
#include <networkit/centrality/HyperKatzCentralityDeltaInPlace.hpp>
#include <networkit/centrality/HyperKatzCentralityDeltaInPlaceTight.hpp>
#include <networkit/centrality/HyperKatzCentralityDeltaSum.hpp>
#include <networkit/centrality/HyperKatzCentralityGraph.hpp>
#include <networkit/centrality/HyperKatzCentralityNaiveSum.hpp>
#include <networkit/centrality/HyperKatzCentralityNaiveSumTight.hpp>
#include <networkit/centrality/HyperKatzCentralityNaiveSumPerLevelEps.hpp>
#include <networkit/io/HMETISHypergraphReader.hpp>

namespace NetworKit {
namespace {

// constexpr auto inputPath =
//     "input/hyperlap-sbm-n500000-d2.5-5-50-c1.5-x1.0-q2.0-2-20-wmajority.hmetis";
constexpr auto inputPath =
    "input/hyperlap-sbm-n500000-d2.5-20-50-c1.5-x1.0-q2.0-2-20-wmajority.hmetis";
// constexpr auto inputPath =
//     "input/hyperlap-sbm-n1000000-d2.5-5-50-c1.5-x1.0-q2.0-2-20-wmajority.hmetis";
// constexpr auto inputPath =
//     "input/hyperlap-sbm-n1000000-d2.5-20-50-c1.5-x1.0-q2.0-2-20-wmajority.hmetis";

// constexpr auto inputPath = "input/abcdh-n10000-d2.5-5-50-c1.5-x0.5-q-2.0-2-20-wstrict.hmetis";
// constexpr auto inputPath = "input/edge-walmart-edge-label.hmetis";
// constexpr auto inputPath = "input/gene-disease.hmetis";
constexpr count topK = 100;

template <typename Centrality>
uint64_t topKHash(const Centrality &centrality) {
    uint64_t hash = 1469598103934665603ULL;
    for (count i = 0; i < topK; ++i) {
        hash ^= centrality.top(i);
        hash *= 1099511628211ULL;
    }
    return hash;
}

uint64_t residentSetBytes() noexcept {
#if defined(__APPLE__)
    mach_task_basic_info_data_t info{};
    mach_msg_type_number_t count = MACH_TASK_BASIC_INFO_COUNT;
    const kern_return_t result = task_info(mach_task_self(), MACH_TASK_BASIC_INFO,
                                           reinterpret_cast<task_info_t>(&info), &count);
    return result == KERN_SUCCESS ? static_cast<uint64_t>(info.resident_size) : 0;
#elif defined(__linux__)
    std::ifstream statm{"/proc/self/statm"};
    uint64_t totalPages = 0;
    uint64_t residentPages = 0;
    statm >> totalPages >> residentPages;
    const long pageSize = sysconf(_SC_PAGESIZE);
    return statm && pageSize > 0 ? residentPages * static_cast<uint64_t>(pageSize) : 0;
#else
    return 0;
#endif
}

class PeakResidentMemorySampler final {
public:
    PeakResidentMemorySampler()
        : baseline{residentSetBytes()}, peak{baseline}, worker{[this] { sample(); }} {}

    ~PeakResidentMemorySampler() { stop(); }

    void stop() {
        if (!worker.joinable())
            return;
        stopped.store(true, std::memory_order_release);
        worker.join();
    }

    uint64_t additionalBytes() const noexcept { return peak > baseline ? peak - baseline : 0; }

    uint64_t peakBytes() const noexcept { return peak; }

    bool isSupported() const noexcept { return baseline != 0; }

private:
    void sample() {
        do {
            peak = std::max(peak, residentSetBytes());
            std::this_thread::sleep_for(std::chrono::milliseconds{1});
        } while (!stopped.load(std::memory_order_acquire));
        peak = std::max(peak, residentSetBytes());
    }

    const uint64_t baseline;
    uint64_t peak;
    std::atomic<bool> stopped{false};
    std::thread worker;
};

} // namespace

TEST(HyperKatzCentralityBenchmark, benchmarkCodex) {
#if !defined(__APPLE__) && !defined(__linux__)
    GTEST_SKIP() << "Resident-memory sampling is only available on macOS and Linux";
#endif

    const Hypergraph hGraph = HMETISHypergraphReader{}.read(inputPath);
    // ASSERT_EQ(hGraph.numberOfNodes(), 88836);
    // ASSERT_EQ(hGraph.numberOfEdges(), 65898);

    PeakResidentMemorySampler memory;
    ASSERT_TRUE(memory.isSupported());

    Aux::Timer timer;
    timer.start();
    auto centrality = std::make_unique<HyperKatzCentrality>(hGraph, topK, true);
    centrality->run();
    timer.stop();
    memory.stop();

    constexpr double bytesPerMiB = 1024.0 * 1024.0;
    std::cout << "HyperKatzCentrality benchmark: " << hGraph.numberOfNodes() << " nodes, "
              << hGraph.numberOfEdges() << " hyperedges, top-k=" << topK << '\n'
              << "  runtime: " << timer.elapsedMilliseconds() << " ms\n"
              << "  peak resident memory: " << static_cast<double>(memory.peakBytes()) / bytesPerMiB
              << " MiB\n"
              << "  additional resident memory after loading: "
              << static_cast<double>(memory.additionalBytes()) / bytesPerMiB << " MiB\n"
              << "  iterations: " << centrality->iterationReached
              << " top 3 edges: " << centrality->top(0) << "," << centrality->top(1) << ","
              << centrality->top(2) << " top 3 score: " << centrality->score(centrality->top(0))
              << "," << centrality->score(centrality->top(1)) << ","
              << centrality->score(centrality->top(2)) << std::endl;

    EXPECT_TRUE(centrality->hasFinished());
}

TEST(HyperKatzCentralityBenchmark, benchmarkNaiveSum) {
#if !defined(__APPLE__) && !defined(__linux__)
    GTEST_SKIP() << "Resident-memory sampling is only available on macOS and Linux";
#endif

    const Hypergraph hGraph = HMETISHypergraphReader{}.read(inputPath);
    // ASSERT_EQ(hGraph.numberOfNodes(), 88836);
    // ASSERT_EQ(hGraph.numberOfEdges(), 65898);

    PeakResidentMemorySampler memory;
    ASSERT_TRUE(memory.isSupported());

    Aux::Timer timer;
    timer.start();
    auto centrality = std::make_unique<HyperKatzCentralityNaiveSum>(hGraph, topK, true);
    centrality->run();
    timer.stop();
    memory.stop();

    constexpr double bytesPerMiB = 1024.0 * 1024.0;
    std::cout << "HyperKatzCentrality benchmark: " << hGraph.numberOfNodes() << " nodes, "
              << hGraph.numberOfEdges() << " hyperedges, top-k=" << topK << '\n'
              << "  runtime: " << timer.elapsedMilliseconds() << " ms\n"
              << "  peak resident memory: " << static_cast<double>(memory.peakBytes()) / bytesPerMiB
              << " MiB\n"
              << "  additional resident memory after loading: "
              << static_cast<double>(memory.additionalBytes()) / bytesPerMiB << " MiB\n"
              << "  iterations: " << centrality->iterationReached
              << " top-100 hash: " << topKHash(*centrality)
              << " top 3 edges: " << centrality->top(0) << "," << centrality->top(1) << ","
              << centrality->top(2) << " top 3 score: " << centrality->score(centrality->top(0))
              << "," << centrality->score(centrality->top(1)) << ","
              << centrality->score(centrality->top(2)) << std::endl;

    EXPECT_TRUE(centrality->hasFinished());
}

TEST(HyperKatzCentralityBenchmark, benchmarkNaiveSumTight) {
#if !defined(__APPLE__) && !defined(__linux__)
    GTEST_SKIP() << "Resident-memory sampling is only available on macOS and Linux";
#endif

    const Hypergraph hGraph = HMETISHypergraphReader{}.read(inputPath);
    PeakResidentMemorySampler memory;
    ASSERT_TRUE(memory.isSupported());

    Aux::Timer timer;
    timer.start();
    const auto constructionStart = std::chrono::steady_clock::now();
    auto centrality = std::make_unique<HyperKatzCentralityNaiveSumTight>(hGraph, topK, true);
    const auto constructionEnd = std::chrono::steady_clock::now();
    centrality->run();
    const auto runEnd = std::chrono::steady_clock::now();
    timer.stop();
    memory.stop();

    constexpr double bytesPerMiB = 1024.0 * 1024.0;
    std::cout << "HyperKatzCentralityNaiveSumTight benchmark: " << hGraph.numberOfNodes()
              << " nodes, " << hGraph.numberOfEdges() << " hyperedges, top-k=" << topK << '\n'
              << "  runtime: " << timer.elapsedMilliseconds() << " ms\n"
              << "  construction: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(constructionEnd
                                                                       - constructionStart)
                     .count()
              << " ms\n"
              << "  iterations: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(runEnd - constructionEnd)
                     .count()
              << " ms\n"
              << "  peak resident memory: " << static_cast<double>(memory.peakBytes()) / bytesPerMiB
              << " MiB\n"
              << "  additional resident memory after loading: "
              << static_cast<double>(memory.additionalBytes()) / bytesPerMiB << " MiB\n"
              << "  iterations: " << centrality->iterationReached
              << " top-100 hash: " << topKHash(*centrality)
              << " top 3 edges: " << centrality->top(0) << "," << centrality->top(1) << ","
              << centrality->top(2) << " top 3 score: " << centrality->score(centrality->top(0))
              << "," << centrality->score(centrality->top(1)) << ","
              << centrality->score(centrality->top(2)) << std::endl;

    EXPECT_TRUE(centrality->hasFinished());
}

TEST(HyperKatzCentralityBenchmark, compareNaiveSumTightTop100) {
    const Hypergraph hGraph = HMETISHypergraphReader{}.read(inputPath);
    std::vector<edgeid> naiveTop(topK);
    count naiveIterations = 0;
    {
        HyperKatzCentralityNaiveSum naive(hGraph, topK, true);
        naive.run();
        naiveIterations = naive.iterationReached;
        for (count rank = 0; rank < topK; ++rank)
            naiveTop[rank] = naive.top(rank);
    }

    HyperKatzCentralityNaiveSumTight tight(hGraph, topK, true);
    tight.run();
    count shared = 0;
    count sameRank = 0;
    for (count rank = 0; rank < topK; ++rank) {
        const edgeid eid = tight.top(rank);
        shared += std::find(naiveTop.begin(), naiveTop.end(), eid) != naiveTop.end();
        sameRank += eid == naiveTop[rank];
    }

    std::cout << "NaiveSum versus NaiveSumTight top-100: shared=" << shared
              << " same-rank=" << sameRank << " Jaccard="
              << static_cast<double>(shared) / (2 * topK - shared)
              << " naive iterations=" << naiveIterations
              << " tight iterations=" << tight.iterationReached
              << "\n  NaiveSum top 10:";
    for (count rank = 0; rank < 10; ++rank)
        std::cout << ' ' << naiveTop[rank];
    std::cout << "\n  NaiveSumTight top 10:";
    for (count rank = 0; rank < 10; ++rank)
        std::cout << ' ' << tight.top(rank);
    std::cout << std::endl;
}

TEST(HyperKatzCentralityBenchmark, compareDeltaInPlaceTightWithNaiveSumTight) {
    const Hypergraph hGraph = HMETISHypergraphReader{}.read(inputPath);
    const count dimension = hGraph.upperEdgeIdBound();
    std::vector<double> naiveScores(dimension, 0.0);
    std::vector<double> naiveBounds(dimension, 0.0);
    std::vector<edgeid> naiveTop(topK);
    count naiveIterations = 0;
    std::chrono::milliseconds naiveTime{0};
    {
        const auto start = std::chrono::steady_clock::now();
        HyperKatzCentralityNaiveSumTight naive(hGraph, topK, true);
        naive.run();
        naiveTime = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::steady_clock::now() - start);
        naiveIterations = naive.iterationReached;
        for (count rank = 0; rank < topK; ++rank)
            naiveTop[rank] = naive.top(rank);
        hGraph.forEdges([&](edgeid eid) {
            naiveScores[eid] = naive.score(eid);
            naiveBounds[eid] = naive.bound(eid);
        });
    }

    const auto start = std::chrono::steady_clock::now();
    HyperKatzCentralityDeltaInPlaceTight inPlace(hGraph, topK, true);
    inPlace.run();
    const auto inPlaceTime = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::steady_clock::now() - start);

    count sameRank = 0;
    for (count rank = 0; rank < topK; ++rank)
        sameRank += naiveTop[rank] == inPlace.top(rank);
    double maxScoreDifference = 0.0;
    double maxBoundDifference = 0.0;
    hGraph.forEdges([&](edgeid eid) {
        maxScoreDifference =
            std::max(maxScoreDifference, std::abs(naiveScores[eid] - inPlace.score(eid)));
        maxBoundDifference =
            std::max(maxBoundDifference, std::abs(naiveBounds[eid] - inPlace.bound(eid)));
    });

    std::cout << "NaiveSumTight versus DeltaInPlaceTight: naive=" << naiveTime.count()
              << " ms in-place=" << inPlaceTime.count() << " ms iterations=" << naiveIterations
              << "/" << inPlace.iterationReached << " same-rank top-100=" << sameRank
              << " maximum score difference=" << maxScoreDifference
              << " maximum bound difference=" << maxBoundDifference << std::endl;
    EXPECT_EQ(inPlace.iterationReached, naiveIterations);
    EXPECT_EQ(sameRank, topK);
    EXPECT_LE(maxScoreDifference, 1e-9);
    EXPECT_LE(maxBoundDifference, 1e-9);
}

TEST(HyperKatzCentralityBenchmark, benchmarkGraph) {
#if !defined(__APPLE__) && !defined(__linux__)
    GTEST_SKIP() << "Resident-memory sampling is only available on macOS and Linux";
#endif

    const Hypergraph hGraph = HMETISHypergraphReader{}.read(inputPath);
    PeakResidentMemorySampler memory;
    ASSERT_TRUE(memory.isSupported());

    Aux::Timer timer;
    timer.start();
    const auto constructionStart = std::chrono::steady_clock::now();
    auto centrality = std::make_unique<HyperKatzCentralityGraph>(hGraph, topK, true);
    const auto constructionEnd = std::chrono::steady_clock::now();
    centrality->run();
    const auto runEnd = std::chrono::steady_clock::now();
    timer.stop();
    memory.stop();

    constexpr double bytesPerMiB = 1024.0 * 1024.0;
    std::cout
        << "HyperKatzCentralityGraph benchmark: " << hGraph.numberOfNodes() << " nodes, "
        << hGraph.numberOfEdges() << " hyperedges, top-k=" << topK << '\n'
        << "  runtime: " << timer.elapsedMilliseconds() << " ms\n"
        << "  construction: "
        << std::chrono::duration_cast<std::chrono::milliseconds>(constructionEnd
                                                                 - constructionStart)
               .count()
        << " ms\n"
        << "  iterations: "
        << std::chrono::duration_cast<std::chrono::milliseconds>(runEnd - constructionEnd).count()
        << " ms\n"
        << "  peak resident memory: " << static_cast<double>(memory.peakBytes()) / bytesPerMiB
        << " MiB\n"
        << "  additional resident memory after loading: "
        << static_cast<double>(memory.additionalBytes()) / bytesPerMiB << " MiB\n"
        << "  iterations: " << centrality->iterationReached
        << " top-100 hash: " << topKHash(*centrality) << " top 3 edges: " << centrality->top(0)
        << "," << centrality->top(1) << "," << centrality->top(2)
        << " top 3 score: " << centrality->score(centrality->top(0)) << ","
        << centrality->score(centrality->top(1)) << "," << centrality->score(centrality->top(2))
        << std::endl;

    EXPECT_TRUE(centrality->hasFinished());
}

TEST(HyperKatzCentralityBenchmark, benchmarkNaiveSumPerLevelEps) {
#if !defined(__APPLE__) && !defined(__linux__)
    GTEST_SKIP() << "Resident-memory sampling is only available on macOS and Linux";
#endif

    const Hypergraph hGraph = HMETISHypergraphReader{}.read(inputPath);
    // ASSERT_EQ(hGraph.numberOfNodes(), 88836);
    // ASSERT_EQ(hGraph.numberOfEdges(), 65898);

    PeakResidentMemorySampler memory;
    ASSERT_TRUE(memory.isSupported());

    Aux::Log::setLogLevel("INFO");

    Aux::Timer timer;
    timer.start();
    auto centrality = std::make_unique<HyperKatzCentralityNaiveSumPerLevelEps>(hGraph, topK, true);
    centrality->run();
    timer.stop();
    memory.stop();

    constexpr double bytesPerMiB = 1024.0 * 1024.0;
    std::cout << "HyperKatzCentrality benchmark: " << hGraph.numberOfNodes() << " nodes, "
              << hGraph.numberOfEdges() << " hyperedges, top-k=" << topK << '\n'
              << "  runtime: " << timer.elapsedMilliseconds() << " ms\n"
              << "  peak resident memory: " << static_cast<double>(memory.peakBytes()) / bytesPerMiB
              << " MiB\n"
              << "  additional resident memory after loading: "
              << static_cast<double>(memory.additionalBytes()) / bytesPerMiB << " MiB\n"
              << "  iterations: " << centrality->iterationReached
              << " top 3 edges: " << centrality->top(0) << "," << centrality->top(1) << ","
              << centrality->top(2) << " top 3 score: " << centrality->score(centrality->top(0))
              << "," << centrality->score(centrality->top(1)) << ","
              << centrality->score(centrality->top(2)) << "\n"
              << "Skipped SpMV: " << centrality->skippedSpMV << std::endl;

    EXPECT_TRUE(centrality->hasFinished());
}

TEST(HyperKatzCentralityBenchmark, benchmarkDeltaSum) {
#if !defined(__APPLE__) && !defined(__linux__)
    GTEST_SKIP() << "Resident-memory sampling is only available on macOS and Linux";
#endif

    const Hypergraph hGraph = HMETISHypergraphReader{}.read(inputPath);
    // ASSERT_EQ(hGraph.numberOfNodes(), 88836);
    // ASSERT_EQ(hGraph.numberOfEdges(), 65898);

    PeakResidentMemorySampler memory;
    ASSERT_TRUE(memory.isSupported());

    Aux::Log::setLogLevel("INFO");

    Aux::Timer timer;
    timer.start();
    auto centrality = std::make_unique<HyperKatzCentralityDeltaSum>(hGraph, topK, true);
    centrality->run();
    timer.stop();
    memory.stop();

    constexpr double bytesPerMiB = 1024.0 * 1024.0;
    std::cout << "HyperKatzCentrality benchmark: " << hGraph.numberOfNodes() << " nodes, "
              << hGraph.numberOfEdges() << " hyperedges, top-k=" << topK << '\n'
              << "  runtime: " << timer.elapsedMilliseconds() << " ms\n"
              << "  peak resident memory: " << static_cast<double>(memory.peakBytes()) / bytesPerMiB
              << " MiB\n"
              << "  additional resident memory after loading: "
              << static_cast<double>(memory.additionalBytes()) / bytesPerMiB << " MiB\n"
              << "  iterations: " << centrality->iterationReached
              << " top 3 edges: " << centrality->top(0) << "," << centrality->top(1) << ","
              << centrality->top(2) << " top 3 score: " << centrality->score(centrality->top(0))
              << "," << centrality->score(centrality->top(1)) << ","
              << centrality->score(centrality->top(2)) << std::endl;

    EXPECT_TRUE(centrality->hasFinished());
}

TEST(HyperKatzCentralityBenchmark, benchmarkDeltaInPlace) {
#if !defined(__APPLE__) && !defined(__linux__)
    GTEST_SKIP() << "Resident-memory sampling is only available on macOS and Linux";
#endif

    const Hypergraph hGraph = HMETISHypergraphReader{}.read(inputPath);
    // ASSERT_EQ(hGraph.numberOfNodes(), 88836);
    // ASSERT_EQ(hGraph.numberOfEdges(), 65898);

    PeakResidentMemorySampler memory;
    ASSERT_TRUE(memory.isSupported());

    Aux::Log::setLogLevel("INFO");

    Aux::Timer timer;
    timer.start();
    auto centrality = std::make_unique<HyperKatzCentralityDeltaInPlace>(hGraph, topK, true);
    centrality->run();
    timer.stop();
    memory.stop();

    constexpr double bytesPerMiB = 1024.0 * 1024.0;
    std::cout << "HyperKatzCentralityDeltaInPlace benchmark: " << hGraph.numberOfNodes()
              << " nodes, "
              << hGraph.numberOfEdges() << " hyperedges, top-k=" << topK << '\n'
              << "  runtime: " << timer.elapsedMilliseconds() << " ms\n"
              << "  peak resident memory: " << static_cast<double>(memory.peakBytes()) / bytesPerMiB
              << " MiB\n"
              << "  additional resident memory after loading: "
              << static_cast<double>(memory.additionalBytes()) / bytesPerMiB << " MiB\n"
              << "  iterations: " << centrality->iterationReached
              << " top-100 hash: " << topKHash(*centrality)
              << " top 3 edges: " << centrality->top(0) << "," << centrality->top(1) << ","
              << centrality->top(2) << " top 3 score: " << centrality->score(centrality->top(0))
              << "," << centrality->score(centrality->top(1)) << ","
              << centrality->score(centrality->top(2)) << std::endl;

    EXPECT_TRUE(centrality->hasFinished());
}

TEST(HyperKatzCentralityBenchmark, benchmarkDeltaInPlaceTight) {
#if !defined(__APPLE__) && !defined(__linux__)
    GTEST_SKIP() << "Resident-memory sampling is only available on macOS and Linux";
#endif

    const Hypergraph hGraph = HMETISHypergraphReader{}.read(inputPath);
    PeakResidentMemorySampler memory;
    ASSERT_TRUE(memory.isSupported());

    Aux::Timer timer;
    timer.start();
    const auto constructionStart = std::chrono::steady_clock::now();
    auto centrality = std::make_unique<HyperKatzCentralityDeltaInPlaceTight>(hGraph, topK, true);
    const auto constructionEnd = std::chrono::steady_clock::now();
    centrality->run();
    const auto runEnd = std::chrono::steady_clock::now();
    timer.stop();
    memory.stop();

    constexpr double bytesPerMiB = 1024.0 * 1024.0;
    std::cout << "HyperKatzCentralityDeltaInPlaceTight benchmark: " << hGraph.numberOfNodes()
              << " nodes, " << hGraph.numberOfEdges() << " hyperedges, top-k=" << topK << '\n'
              << "  runtime: " << timer.elapsedMilliseconds() << " ms\n"
              << "  construction: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(constructionEnd
                                                                       - constructionStart)
                     .count()
              << " ms\n"
              << "  iterations: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(runEnd - constructionEnd)
                     .count()
              << " ms\n"
              << "  peak resident memory: " << static_cast<double>(memory.peakBytes()) / bytesPerMiB
              << " MiB\n"
              << "  additional resident memory after loading: "
              << static_cast<double>(memory.additionalBytes()) / bytesPerMiB << " MiB\n"
              << "  iterations: " << centrality->iterationReached
              << " top-100 hash: " << topKHash(*centrality)
              << " top 3 edges: " << centrality->top(0) << "," << centrality->top(1) << ","
              << centrality->top(2) << " top 3 score: " << centrality->score(centrality->top(0))
              << "," << centrality->score(centrality->top(1)) << ","
              << centrality->score(centrality->top(2)) << std::endl;

    EXPECT_TRUE(centrality->hasFinished());
}

} // namespace NetworKit
