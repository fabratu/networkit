
/*
 * File:   Correspondences.hpp
 * Author: Roland Glantz (roland.glantz@kit.edu)
 *
 */

// networkit-format

#ifndef NETWORKIT_CORRESPONDENCES_CORRESPONDENCES_HPP_
#define NETWORKIT_CORRESPONDENCES_CORRESPONDENCES_HPP_

#include <networkit/structures/Partition.hpp>

namespace NetworKit {

class Correspondences {

public:
    count numberOfElements;

    count cardPartition1, cardPartition2;

    std::map<index, count> cardinalityOfCluster1;

    std::map<index, count> cardinalityOfCluster2;

    std::vector<std::vector<count>> distributions;

    void normalize(const Partition &partitionA, const Partition &partitionB,
                   Partition &normalPartitionA, Partition &normalPartitionB);

    void getDistributions(Partition &partition1, Partition &partition2);

    count peak(count cardCluster2, count overlap);

    count getBoundPeak(std::vector<count> &distriB_s, std::vector<count> &distriB_t);

    count potFracs(std::vector<int> &belongs, count &sFracPotNew, count &tFracPotNew,
                   std::vector<count> &distriB_s, std::vector<count> &distriB_t);

    count greedyDescent(index s, index t, count &bestFrac, count &currentFrac,
                        std::vector<index> &insertedAt, count &position, std::vector<int> &belongs,
                        std::vector<int> &bestBelongs, std::vector<count> &distriB_s,
                        std::vector<count> &distriB_t);

    count greedyBB(index s, index t, count bestFrac, count currentFrac, std::vector<int> &belongs,
                   std::vector<int> &bestBelongs, std::vector<count> &distriB_s,
                   std::vector<count> &distriB_t);

    void getBestBelongsPrime(std::vector<int> &bestBelongs, std::vector<int> &bestBelongsPrime);

    void evaluateCorrespondence(std::vector<int> &bestBelongs, std::vector<int> &bestBelongsPrime,
                                count &clusterCardinality, count &clusterCardinalityPrime,
                                count &elementCardinality, count &elementCardinalityPrime,
                                double symDiff, double &size, double &quality);

    count minCut(index s, index t, std::vector<int> &bestBelongs);

    count gusfield(index &bestS, index &bestT, std::vector<int> &bestBestBelongs);

    count run(const Partition &partition1, const Partition &partition2);
};

} /* namespace NetworKit */
#endif // NETWORKIT_CORRESPONDENCES_CORRESPONDENCES_HPP_
