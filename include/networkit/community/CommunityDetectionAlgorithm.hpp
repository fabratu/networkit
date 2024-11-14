/*
 * CommunityDetectionAlgorithm.hpp
 *
 *  Created on: 30.10.2012
 *      Author: Christian Staudt
 */

#ifndef NETWORKIT_COMMUNITY_COMMUNITY_DETECTION_ALGORITHM_HPP_
#define NETWORKIT_COMMUNITY_COMMUNITY_DETECTION_ALGORITHM_HPP_

#include <type_traits>

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/structures/Partition.hpp>

namespace NetworKit {

/**
 * @ingroup community
 * Abstract base class for community detection/graph clustering algorithms.
 */
template <typename GraphType>
class CommunityDetectionAlgorithm : public Algorithm {
public:
    /**
     * A community detection algorithm operates on a graph, so the constructor expects a graph.
     *
     * @param[in] G input graph
     */
    CommunityDetectionAlgorithm(const GraphType &G) : Algorithm(), G(&G), result(0) {
        // currently our community detection methods are not defined on directed graphs
        if constexpr (std::is_same<GraphType, Graph>::value) {
            if (G.isDirected()) throw std::runtime_error("This community detection method is undefined on directed graphs");
        }
    }

    /**
     * A community detection algorithm operates on a graph, so the constructor expects a graph.
     *
     * @param[in] G input graph
     * @param[in] baseClustering optional; the algorithm will start from the given clustering.
     */
    CommunityDetectionAlgorithm(const GraphType &G, Partition baseClustering) : Algorithm(), G(&G), result(std::move(baseClustering)) {}

    /** Default destructor */
    ~CommunityDetectionAlgorithm() override = default;

    /**
     * Apply algorithm to graph
     */
    void run() override = 0;

    /**
     * Returns the result of the run method or throws an error, if the algorithm hasn't run yet.
     * @return partition of the node set
     */
    const Partition &getPartition() const {
        if (!hasRun) {
            throw std::runtime_error("Call run()-function first.");
        }
        return result;
    };

protected:
    const GraphType *G;
    Partition result;
};

} /* namespace NetworKit */
#endif // NETWORKIT_COMMUNITY_COMMUNITY_DETECTION_ALGORITHM_HPP_
