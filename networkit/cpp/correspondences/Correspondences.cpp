/*
 * Correspondences.cpp
 *
 *  Created on: 22.01.2013
 *      Author: Roland Glantz (roland.glantz@kit.edu)
 */

// networkit-format

#include <functional>
#include <iterator>
#include <numeric>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/correspondences/Correspondences.hpp>

namespace NetworKit {

/**********************************************************************/
/*                              normalize                             */
/**********************************************************************/
void Correspondences::normalize(const Partition &partitionA, const Partition &partitionB,
                                Partition &normalPartitionA, Partition &normalPartitionB) {

    // TRACE("partitionA is ", partitionA.getVector());
    // TRACE("partitionB is ", partitionB.getVector());

    // consecutively renumber the elements of partitionA such that
    // the new number of any element in cluster i is smaller than
    // the new number of any element in cluster j whenever i < j
    std::vector<index> newNumber(partitionA.numberOfElements(), 0);
    std::map<index, count> cardinalityOfClusterA = partitionA.subsetSizeMap();
    std::vector<count> currSlotInCluster(partitionA.numberOfSubsets(), 0);
    for (count c = 1; c < partitionA.numberOfSubsets(); c++) {
        currSlotInCluster[c] = currSlotInCluster[c - 1] + cardinalityOfClusterA[c - 1];
    }
    for (count e = 0; e < partitionA.numberOfElements(); e++) {
        newNumber[e] = currSlotInCluster[partitionA.subsetOf(e)];
        (currSlotInCluster[partitionA.subsetOf(e)])++;
    }

    // normalization of partitionA w.r.t newNumber
    normalPartitionA = Partition(partitionA.numberOfElements());
    normalPartitionA.setUpperBound(partitionA.numberOfSubsets());
    for (count e = 0; e < partitionA.numberOfElements(); e++) {
        normalPartitionA.addToSubset(partitionA.subsetOf(e), newNumber[e]);
    }

    // normalization of partitionB w.r.t newNumber
    normalPartitionB = Partition(partitionB.numberOfElements());
    normalPartitionB.setUpperBound(partitionB.numberOfSubsets());
    for (count e = 0; e < partitionB.numberOfElements(); e++) {
        normalPartitionB.addToSubset(partitionB.subsetOf(e), newNumber[e]);
    }
}

/**********************************************************************/
/*                           getDistributions                         */
/**********************************************************************/
void Correspondences::getDistributions(Partition &partition1, Partition &partition2) {

    numberOfElements = partition1.numberOfElements();
    cardPartition1 = partition1.upperBound();
    cardPartition2 = partition2.upperBound();
    cardinalityOfCluster1 = partition1.subsetSizeMap();
    cardinalityOfCluster2 = partition2.subsetSizeMap();

    // represent clusters by sets, partitions by vectors of sets
    std::vector<std::set<index>> clustersOfPartition1(cardPartition1);
    std::vector<std::set<index>> clustersOfPartition2(cardPartition2);
    for (index i1 = 0; i1 < cardPartition1; i1++) {
        clustersOfPartition1[i1] = partition1.getMembers(i1);
    }
    for (index i2 = 0; i2 < cardPartition2; i2++) {
        clustersOfPartition2[i2] = partition2.getMembers(i2);
    }

    // Get distributions from vectors of sets
    distributions.resize(cardPartition1, std::vector<count>(cardPartition2));
    for (index i1 = 0; i1 < cardPartition1; i1++) {
        for (index i2 = 0; i2 < cardPartition2; i2++) {
            std::vector<index> interset;
            std::set_intersection((clustersOfPartition1[i1]).begin(),
                                  (clustersOfPartition1[i1]).end(),
                                  (clustersOfPartition2[i2]).begin(),
                                  (clustersOfPartition2[i2]).end(), std::back_inserter(interset));
            distributions[i1][i2] = (count)interset.size();
        }
    }
    // TRACE("cardinalityOfCluster1 is ", cardinalityOfCluster1);
    // TRACE("distributions is ", distributions);
}

/**********************************************************************/
/*                                 peak                               */
/**********************************************************************/
count Correspondences::peak(count cardCluster2, count overlap) {
    return (std::min(overlap, cardCluster2 - overlap));
}

/**********************************************************************/
/*                              getBoundPeak                          */
/**********************************************************************/
count Correspondences::getBoundPeak(std::vector<count> &distriB_s, std::vector<count> &distriB_t) {
    count ret = 0;
    for (count i2 = 0; i2 < distriB_s.size(); i2++) {
        ret += std::min(distriB_s[i2], distriB_t[i2]);
    }
    return (ret);
}

/**********************************************************************/
/*                               potFracs                             */
/**********************************************************************/
count Correspondences::potFracs(std::vector<int> &belongs, count &sFracPotNew, count &tFracPotNew,
                                std::vector<count> &distriB_s, std::vector<count> &distriB_t) {

    index bestCluster1 = cardPartition1; // no such cluster
    int64_t maxDiff = 0;                 // absolute difference between sFracPot and tFracPot
    for (index i1 = 0; i1 < cardPartition1; i1++) {
        if (belongs[i1] == 0) { // i1 has not been inserted yet
            index sFracPot = 0;
            index tFracPot = 0;
            for (index i2 = 0; i2 < cardPartition2; i2++) {
                sFracPot += std::min(distriB_s[i2] + distributions[i1][i2], distriB_t[i2]);
                tFracPot += std::min(distriB_s[i2], distriB_t[i2] + distributions[i1][i2]);
            }
            int64_t diff = labs(static_cast<int64_t>(sFracPot) - static_cast<int64_t>(tFracPot));
            if (diff >= maxDiff) {
                maxDiff = diff;
                bestCluster1 = i1;
                sFracPotNew = sFracPot;
                tFracPotNew = tFracPot;
            }
        }
    }
    return (bestCluster1);
}

/**********************************************************************/
/*                            greedyDescent                           */
/**********************************************************************/
index Correspondences::greedyDescent(index s, index t, count &bestFrac, count &currentFrac,
                                     std::vector<index> &position2cluster, count &position,
                                     std::vector<int> &belongs, std::vector<int> &bestBelongs,
                                     std::vector<count> &distriB_s, std::vector<count> &distriB_t) {
    if (currentFrac < bestFrac) {
        count newCluster = cardPartition1; // no such cluster
        while (currentFrac < bestFrac) {
            count sFracPotNew = 0;
            count tFracPotNew = 0;
            newCluster = potFracs(belongs, sFracPotNew, tFracPotNew, distriB_s, distriB_t);
            if (newCluster < cardPartition1) { // newCluster is valid cluster ID for greedy descent
                position2cluster[position] = newCluster;
                position++;
                // check whether cluster with ID newCluster goes to U_{\mathcal{B}_s} or to
                // U_{\mathcal{B}_t} if((s == 4) && (t == 2)) TRACE("sFracPotNew is ", sFracPotNew,
                // " und tFracPotNew is ", tFracPotNew);
                if (sFracPotNew < tFracPotNew) {
                    // if((s == 4) && (t == 2)) TRACE("Cluster ", position2cluster[position - 1], "
                    // goes to s = ", s);
                    belongs[newCluster] = 1;
                    std::transform(distributions[newCluster].begin(),
                                   distributions[newCluster].end(), distriB_s.begin(),
                                   distriB_s.begin(), std::plus<count>());
                } else {
                    // if((s == 4) && (t == 2)) TRACE("Cluster ", position2cluster[position - 1], "
                    // goes to t = ", t);
                    belongs[newCluster] = 2;
                    std::transform(distributions[newCluster].begin(),
                                   distributions[newCluster].end(), distriB_t.begin(),
                                   distriB_t.begin(), std::plus<count>());
                }
                currentFrac = getBoundPeak(distriB_s, distriB_t);
            } else { // newCluster == cardPartition1, i.e., we have a new solution
                if (currentFrac < bestFrac) {
                    bestFrac = currentFrac;
                    std::copy(belongs.begin(), belongs.end(), bestBelongs.begin());
                    // if((s == 4) && (t == 2)) TRACE("Found a solution: ", bestBelongs);
                }
            }
        }
        return (position - 1);
    } else {
        return (position);
    }
}

/**********************************************************************/
/*                             greedyBB                               */
/**********************************************************************/
count Correspondences::greedyBB(index s, index t, count bestFrac, count currentFrac,
                                std::vector<int> &belongs, std::vector<int> &bestBelongs,
                                std::vector<count> &distriB_s, std::vector<count> &distriB_t) {

    count position = 0;                     // initialization
    count maxPosition = cardPartition1 - 2; // s and t are never inserted
    // next vector is for keeping track of when clusters were inserted
    std::vector<index> position2cluster(maxPosition,
                                        cardPartition1); // cardPartition1 means no cluster
    std::vector<bool> doneWith(maxPosition, false);      // needed to avoid infinite loops

    bool done = false;
    while (done == false) {
        position = greedyDescent(s, t, bestFrac, currentFrac, position2cluster, position, belongs,
                                 bestBelongs, distriB_s, distriB_t);

        // if((s == 4) && (t == 2)) TRACE(" position is ", position);
        // if((s == 4) && (t == 2)) TRACE(" currentFrac und bestFrac sind ",  currentFrac, " und ",
        // bestFrac);

        // backtracking
        //...roll back
        doneWith[position] = true;
        index cluster = cardPartition1; // no such cluster;
        while ((doneWith[position] == true) && (position > 0)) {
            cluster = position2cluster[position];
            if (belongs[cluster] == 1) {
                std::transform(distriB_s.begin(), distriB_s.end(), distributions[cluster].begin(),
                               distriB_s.begin(), std::minus<count>());
            }
            if (belongs[cluster] == 2) {
                std::transform(distriB_t.begin(), distriB_t.end(), distributions[cluster].begin(),
                               distriB_t.begin(), std::minus<count>());
            }
            belongs[cluster] = 0;
            doneWith[position] = false;
            position--;
        }

        if (doneWith[position] == false) { // if not done yet
            //...new direction for next greedy descent
            cluster = position2cluster[position];
            if (belongs[cluster] == 1) {
                belongs[cluster] = 2;
                std::transform(distriB_s.begin(), distriB_s.end(), distributions[cluster].begin(),
                               distriB_s.begin(), std::minus<count>());
                std::transform(distributions[cluster].begin(), distributions[cluster].end(),
                               distriB_t.begin(), distriB_t.begin(), std::plus<count>());
            } else {
                if (belongs[cluster] == 2) {
                    belongs[cluster] = 1;
                    std::transform(distriB_t.begin(), distriB_t.end(),
                                   distributions[cluster].begin(), distriB_t.begin(),
                                   std::minus<count>());
                    std::transform(distributions[cluster].begin(), distributions[cluster].end(),
                                   distriB_s.begin(), distriB_s.begin(), std::plus<count>());
                }
            }
            currentFrac = getBoundPeak(distriB_s, distriB_t);
            doneWith[position] = true;

            // if((s == 4) && (t == 2)) TRACE("BACKTRACKING:::currentFrac is ", currentFrac);
            // if((s == 4) && (t == 2)) TRACE("position and cardPartition1 are ", position, " and ",
            // cardPartition1); if((s == 4) && (t == 2)) TRACE("belongs is ", belongs); if((s == 4)
            // && (t == 2)) TRACE("distriB_s is ", distriB_s); if((s == 4) && (t == 2))
            // TRACE("distriB_t is ", distriB_t); if((s == 4) && (t == 2)) TRACE("bestFrac is ",
            // bestFrac); if((s == 4) && (t == 2)) TRACE(" ");
        } else {
            done = true;
        }
    }
    return (bestFrac);
}

/**********************************************************************/
/*                        getBestBelongsPrime                         */
/**********************************************************************/
void Correspondences::getBestBelongsPrime(std::vector<int> &bestBelongs,
                                          std::vector<int> &bestBelongsPrime) {
    for (index i2 = 0; i2 < cardPartition2; i2++) {
        count thresh = (cardinalityOfCluster2[i2]) / 2;
        count fit = 0;
        for (index i1 = 0; i1 < cardPartition1; i1++) {
            if (bestBelongs[i1] == 1) {
                fit += distributions[i1][i2];
            }
        }
        if (fit > thresh) {
            bestBelongsPrime[i2] = 1;
        } else {
            bestBelongsPrime[i2] = 2;
        }
    }
}

/**********************************************************************/
/*                       evaluateCorrespondence                       */
/**********************************************************************/
void Correspondences::evaluateCorrespondence(
    std::vector<int> &bestBelongs, std::vector<int> &bestBelongsPrime, count &clusterCardinality,
    count &clusterCardinalityPrime, count &elementCardinality, count &elementCardinalityPrime,
    double symDiff, double &size, double &quality) {
    // cardinalities w.r.t. bestBelongs
    count clusterCardinalityA = 0;
    count elementCardinalityA = 0;
    for (index i1 = 0; i1 < cardPartition1; i1++) {
        if (bestBelongs[i1] == 1) {
            clusterCardinalityA++;
            elementCardinalityA += cardinalityOfCluster1[i1];
        }
    }

    // cardinalities w.r.t. bestBelongsPrime
    count clusterCardinalityAPrime = 0;
    count elementCardinalityAPrime = 0;
    for (index i2 = 0; i2 < cardPartition2; i2++) {
        if (bestBelongsPrime[i2] == 1) {
            clusterCardinalityAPrime++;
            elementCardinalityAPrime += cardinalityOfCluster2[i2];
        }
    }

    // go with $\mathcal{B}$, as opposed to $\mathcal{C} \setminus \mathcal{B}$, if
    // $U_{\mathcal{B}}$ has fewer elements
    count elementCardinalityB = numberOfElements - elementCardinalityA;
    if (elementCardinalityA <= elementCardinalityB) {
        clusterCardinality = clusterCardinalityA;
        elementCardinality = elementCardinalityA;
        clusterCardinalityPrime = clusterCardinalityAPrime;
        elementCardinalityPrime = elementCardinalityAPrime;
    } else {
        clusterCardinality = cardPartition1 - clusterCardinalityA;
        elementCardinality = elementCardinalityB;
        clusterCardinalityPrime = cardPartition2 - clusterCardinalityAPrime;
        elementCardinalityPrime = numberOfElements - elementCardinalityAPrime;
    }
    double sizeIntersection = (elementCardinality + elementCardinalityPrime - symDiff) / 2.0;
    size = elementCardinality + elementCardinalityPrime - sizeIntersection;
    quality = sizeIntersection / size;
}

/**********************************************************************/
/*                               minCut                               */
/**********************************************************************/
count Correspondences::minCut(index s, index t, std::vector<int> &bestBelongs) {

    // belongs[] expresses membership of a cluster $c$ from partition1 as follows.
    // belongs[c] = 1: $c \in $\mathcal{B}_s$
    // belongs[c] = 2: $c \in $\mathcal{B}_t$
    // belongs[c] = 0: $c \notin $\mathcal{B}_s$, $c \notin $\mathcal{B}_t$
    std::vector<int> belongs(cardPartition1, 0);

    std::vector<count> distriB_s(
        cardPartition2, 0); // distribution of U_{\mathcal{B}_s} over the clusters of partition2
    std::vector<count> distriB_t(
        cardPartition2, 0); // distribution of U_{\mathcal{B}_t} over the clusters of partition2

    belongs[s] = 1;
    std::transform(distributions[s].begin(), distributions[s].end(), distriB_s.begin(),
                   distriB_s.begin(), std::plus<count>());
    belongs[t] = 2;
    std::transform(distributions[t].begin(), distributions[t].end(), distriB_t.begin(),
                   distriB_t.begin(), std::plus<count>());

    if (cardPartition1 > 2) {
        return (greedyBB(s, t, numberOfElements, getBoundPeak(distriB_s, distriB_t), belongs,
                         bestBelongs, distriB_s, distriB_t));
    } else {
        std::copy(belongs.begin(), belongs.end(), bestBelongs.begin());
        return (getBoundPeak(distriB_s, distriB_t));
    }
}

/**********************************************************************/
/*                              gusfield                              */
/**********************************************************************/
count Correspondences::gusfield(index &bestS, index &bestT, std::vector<int> &bestBestBelongs) {

    // Gomory-Hu tree
    std::vector<index> gomoryHuParent(cardPartition1, cardPartition1);
    std::vector<count> cutWithGomoryHuParent(cardPartition1, numberOfElements);

    // build first Gomory-Hu tree (star with vertex 0 in the center
    for (index s = 1; s < cardPartition1; s++) {
        gomoryHuParent[s] = 0;
    }

    // rebuild the Gomory-Hu tree
    count cMin = numberOfElements;
    bestBestBelongs.resize(cardPartition1, 0);
    for (index s = 1; s < cardPartition1; s++) {
        index t = gomoryHuParent[s];
        // TRACE(" s und t sind ", s, " und  ", t);
        std::vector<int> bestBelongs(cardPartition1, 0); // best belongs per s-t cut
        count c = minCut(s, t, bestBelongs);
        INFO("   minCut zwischen ", s, " und  ", t, " ist ", c);

        std::vector<int> bestBelongsPrime(cardPartition2, 0);
        getBestBelongsPrime(bestBelongs, bestBelongsPrime);
        count clusterCardinality = 0;
        count clusterCardinalityPrime = 0;
        count elementCardinality = 0;
        count elementCardinalityPrime = 0;
        double size = 0.0;
        double quality = 0.0;
        evaluateCorrespondence(bestBelongs, bestBelongsPrime, clusterCardinality,
                               clusterCardinalityPrime, elementCardinality, elementCardinalityPrime,
                               (double)c, size, quality);
        INFO("   Cluster cardinalities are ", clusterCardinality, " and  ",
             clusterCardinalityPrime);
        INFO("   Element cardinalities are ", elementCardinality, " and  ",
             elementCardinalityPrime);
        INFO("   Size and quality are ", size, " and ", quality);
        INFO(" ");
        INFO(" ");

        // TRACE("   bestBelongs is ", bestBelongs);
        cutWithGomoryHuParent[s] = c; // label edges of Gomory-Hu tree
        if (c < cMin) {               // update value of minimal Cut
            cMin = c;
            bestS = s;
            bestT = t;
            std::copy(bestBelongs.begin(), bestBelongs.end(), bestBestBelongs.begin());
        }
        // relink edges of Gomory-Hu tree
        for (index i = s + 1; i < cardPartition1; i++) {
            if (gomoryHuParent[i] == t) {  // i has same parent as s, that is t
                if (bestBelongs[i] == 1) { // i is on the same side of the cut as s
                    // TRACE("       gomoryHuParent of ", i, " is set to ", s);
                    gomoryHuParent[i] = s;
                }
            }
        }
    }

    INFO("gomoryHuParent is ", gomoryHuParent);
    INFO("cutWithGomoryHuParent is ", cutWithGomoryHuParent);

    return (cMin);
}

// 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5
// 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5

// 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 4, 4, 5, 5, 5, 4, 4, 4, 5, 5
// 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 4, 4, 5, 5, 5, 4, 4, 4, 5, 5

/**********************************************************************/
/*                                run                                 */
/**********************************************************************/
count Correspondences::run(const Partition &partitionA, const Partition &partitionB) {

    Partition partition1, partition2;
    index bestS, bestT;
    std::vector<int> bestBestBelongs;

    normalize(partitionA, partitionB, partition1, partition2);
    getDistributions(partition1, partition2);
    // TRACE("partition1 is ", partition1.getVector());
    // TRACE("partition2 is ", partition2.getVector());
    count minCut = gusfield(bestS, bestT, bestBestBelongs);
    // TRACE("bestS is ", bestS);
    // TRACE("bestT is ", bestT);

    // std::vector<int> bestBestBelongsPrime(cardPartition2, 0);
    // getBestBelongsPrime(bestBestBelongs, bestBestBelongsPrime);
    // TRACE("bestBestBelongs is ", bestBestBelongs);
    // TRACE("bestBestBelongsPrime is ", bestBestBelongsPrime);
    return (minCut);
}

} /* namespace NetworKit */
