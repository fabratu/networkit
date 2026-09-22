#include <algorithm>
#include <atomic>
#include <chrono>
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
#include <networkit/io/HMETISHypergraphReader.hpp>

namespace NetworKit {
namespace {

constexpr auto inputPath = "input/abcdh-n10000-d2.5-5-50-c1.5-x0.5-q-2.0-2-20-wstrict.hmetis";
constexpr count topK = 100;

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

TEST(HyperKatzCentralityBenchmark, benchmarkRuntimeAndMemory) {
#if !defined(__APPLE__) && !defined(__linux__)
    GTEST_SKIP() << "Resident-memory sampling is only available on macOS and Linux";
#endif

    const Hypergraph hGraph = HMETISHypergraphReader{}.read(inputPath);
    ASSERT_EQ(hGraph.numberOfNodes(), 10000);
    ASSERT_EQ(hGraph.numberOfEdges(), 7867);

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
              << "  iterations: " << centrality->iterationReached << std::endl;

    EXPECT_TRUE(centrality->hasFinished());
}

} // namespace NetworKit
