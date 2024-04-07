#ifndef NETWORKIT_MATCHING_B_SUITOR_MATCHER_HPP_
#define NETWORKIT_MATCHING_B_SUITOR_MATCHER_HPP_

#include <algorithm>
#include <map>
#include <set>
#include <networkit/auxiliary/BucketPQ.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/matching/BMatcher.hpp>

namespace NetworKit {

struct DynBNode {
    node id;
    edgeweight weight;

    DynBNode() = default;
    DynBNode(node n, edgeweight w) : id(n), weight(w) {}

    bool operator==(const DynBNode &other) const {
        return id == other.id;
        // return id == other.id && weight == other.weight;
    }
    bool operator!=(const DynBNode &other) const {
        return id != other.id;
        // return id != other.id || weight != other.weight;
    }
};

struct DynBNodeComparator {
    bool operator() (const DynBNode& lhs, const DynBNode& rhs) const{
        return lhs.weight < rhs.weight || (lhs.weight == rhs.weight && lhs.id > rhs.id);
    }
};

struct DynBNodeMatchesInfo {
    // std::vector<DynBNode> partners;
    DynBNode min; // (none, 0) if partners still has free capacity
    count max_size;
    count num_visits;
    node looseEnd;
    node activeLooseEnd;
    std::vector<DynBNode> lastUpdate;
    // std::map<int, DynBNode> partnersHashed;
    std::set<DynBNode, DynBNodeComparator> partnersBucket;

    DynBNodeMatchesInfo() = default;

    DynBNodeMatchesInfo(count b) {
        // partners.reserve(b);
        min = DynBNode(none, 0);
        max_size = b;
        num_visits = 0;
        looseEnd = none;
        activeLooseEnd = none;
    }

    // TODO: fix
    bool hasPartner(node u) {
        bool same = false;
        DynBNode tester{u, none};
        // INFO("Started hasPartner against ", tester.id);
        for(auto it = partnersBucket.begin(); it != partnersBucket.end(); ++it) {
            if(tester == *it) {
                same = true;
            }
        }
        INFO("Has partner: ", u, " True: ", same);
        return same;

        // bool tester = partnersBucket.contains(DynBNode{u,none});
        // INFO("Has partner: ", u, " True: ", tester);
        // return partnersBucket.contains(DynBNode{u,none});
        
        // NOT SO OLD CODE
        // return partnersHashed.contains(u);

        // OLD CODE
        // return std::find_if(partners.begin(), partners.end(),
        //                     [u](const DynBNode &v) { return v.id == u; })
        //        != partners.end();
    }

    DynBNode popMinIfFull() {
        // INFO("Called popMinIfFull for min: ", min.id, " (weight: ", min.weight, ")");
        if (partnersBucket.size() < max_size) {
            return {none, 0};
        } else {
            auto ret = min;
            remove(min.id);
            return ret;
        }

        // NOT SO OLD CODE
        // if (partnersHashed.size() < max_size) {
        //     return {none, 0};
        // } else {
        //     auto ret = min;
        //     remove(min.id);
        //     return ret;
        // }

        // OLD CODE
        // if (partners.size() < max_size) {
        //     return {none, 0};
        // } else {
        //     auto ret = min;
        //     remove(min.id);
        //     return ret;
        // }
    }

    DynBNode insert(const DynBNode &u, bool update = true) {
        // INFO("Want to insert ", u.id);
        if (hasPartner(u.id)) 
            return {none, 0};

        DynBNode prevMin = popMinIfFull();
        // TODO: activate again if working
        // assert(partners.size() < max_size);

        partnersBucket.insert(u);
        
        if (update)
            lastUpdate.emplace_back(u);
        if (partnersBucket.size() >= max_size && !partnersBucket.empty()) {
            min = DynBNode{(*partnersBucket.begin()).id, (*partnersBucket.begin()).weight};
            INFO("Min after insertion: ", min.id, " with weight: ", min.weight);
        }

        return prevMin;


        // NOT SO OLD CODE
        // partnersHashed[u.id] = u;
        // if (update)
        //     lastUpdate.emplace_back(u);
        // if (partnersHashed.size() >= max_size && !partnersHashed.empty()) {
        //     min = partnersHashed.begin()->second;
        // }

        // return prevMin;

        // OLD CODE
        // if (hasPartner(u.id))
        //     return {none, 0};

        // DynBNode prevMin = popMinIfFull();
        // // TODO: activate again if working
        // // assert(partners.size() < max_size);

        // partners.emplace_back(u);
        // if (update)
        //     lastUpdate.emplace_back(u);
        // if (partners.size() >= max_size && !partners.empty()) {
        //     min = *std::min_element(partners.begin(), partners.end(),
        //                             [](const DynBNode &x, const DynBNode &y) {
        //                                 if (x.weight == y.weight) {
        //                                     return x.id > y.id;
        //                                 }
        //                                 return x.weight < y.weight;
        //                             });
        // }

        // return prevMin;
    }

    // TODO: fix
    void remove(node u) {
        looseEnd = u;
        std::set<DynBNode, DynBNodeComparator>::iterator myIt;
        for(auto it = partnersBucket.begin(); it != partnersBucket.end(); ++it) {
            if((*it).id == u) {
                myIt = it;
                break;
            }
        }
        partnersBucket.erase(myIt);
        // partnersBucket.erase(DynBNode{u, none});
        min = DynBNode(none, 0);

        // NOT SO OLD CODE
        // looseEnd = u;
        // partnersHashed.erase(u);
        // min = DynBNode(none, 0);

        // OLD CODE
        // looseEnd = u;
        // partners.erase(std::remove_if(partners.begin(), partners.end(),
        //                               [u](const DynBNode &v) {
        //                                   if (v.id == u) {
        //                                       return true;
        //                                   }
        //                                   return false;
        //                               }),
        //                partners.end());
        // INFO("Removed: ", u);
        // min = DynBNode(none, 0);
    }

    void sort() {
        // OLD CODE
        // std::sort(partners.begin(), partners.end(), [](const DynBNode &u, const DynBNode &v) {
        //     return (u.weight > v.weight || (u.weight == v.weight && u.id < v.id));
        // });
    }

    friend std::ostream &operator<<(std::ostream &out, const DynBNodeMatchesInfo &nmi) {
        out << "[";
        for (auto i = nmi.partnersBucket.begin(); i != nmi.partnersBucket.end(); ++i) {
            out << (*i).id << ' ';
        }
        out << "]";
        return out;

        // NOT OS OLD CODE
        // out << "[";
        // for (auto i = nmi.partnersHashed.begin(); i != nmi.partnersHashed.end(); ++i) {
        //     out << (*i).second.id << ' ';
        // }
        // out << "]";
        // return out;
        
        // OLD CODE
        // out << "[";
        // for (auto i = nmi.partners.begin(); i != nmi.partners.end(); ++i) {
        //     out << (*i).id << ' ';
        // }
        // out << "]";
        // return out;
    }

    bool operator==(const DynBNodeMatchesInfo &other) const {
        bool same = true;
        auto myIt = partnersBucket.begin(); 
        for (auto i = other.partnersBucket.begin(); i != other.partnersBucket.end(); ++i) {
            (*i).id == (*myIt).id ? same = true : same = false;
            std::advance(myIt, 1);
        }
        return same;

        // NOT SO OLD CODE
        // bool same = true;
        // for (auto i = other.partnersHashed.begin(); i != other.partnersHashed.end(); ++i) {
        //     (*i).second == partnersHashed.at((*i).first) ? same = true : same = false;
        // }
        // return same;

        // OLD CODE
        // bool same = true;
        // for (size_t i = 0; i < other.partners.size(); i++) {
        //     partners[i] == other.partners[i] ? same = true : same = false;
        // }
        // return same;
    }

    DynBNode popLastUpdate() {
        // INFO("Started last update");
        if (lastUpdate.size() == 0) {
            return {none, 0};
        } else {
            auto ret = lastUpdate.back();
            // INFO("Last update has item: ", ret.id, " (size: ", lastUpdate.size(), ")");
            lastUpdate.pop_back();
            return ret;
        }
    }

    void reset() {
        // INFO("Started reset.");
        lastUpdate.clear();
        // auto max = lastUpdate.size();
        // for (size_t i = 0; i < max; i++) {
        //     lastUpdate.pop_back();
        // }
        // INFO("Finished reset.");
    };
};

/**
 * @ingroup matching
 * B-Suitor matching finding algorithm.
 */
class BSuitorMatcher : public BMatcher {
public:
    /**
     * Computes a 1/2-approximate maximum weight b-matching of an undirected weighted Graph @c G
     * using the sequential b-Suitor algorithm published by Khan et al. in "Efficient
     * Approximation Algorithms For Weighted B-Matching", SIAM Journal on Scientific Computing,
     * Vol. 38, Iss. 5 (2016).
     *
     * @param G An undirected graph.
     * @param b A vector of @a b values that represents the max number of edges per vertex @a v
     * in the b-Matching (b.at(v)).
     */
    BSuitorMatcher(const Graph &G, const std::vector<count> &b);

    /**
     * @param G An undirected graph.
     * @param b A value @a b that represents the max number of edges per vertex in the
     * b-Matching. Defaults to the ordinary 1-Matching.
     */
    BSuitorMatcher(const Graph &G, count b = 1);

    /**
     * @param G  An undirected graph.
     * @param path  A path to a file containing @a b values that represents the max number of
     * edges per vertex in the b-Matching.
     */
    BSuitorMatcher(const Graph &G, const std::string &path);

    ~BSuitorMatcher() override = default;

    /**
     * Runs the algorithm.
     */
    void run() override;

    void buildBMatching();

protected:
    std::vector<std::unique_ptr<DynBNodeMatchesInfo>> Suitors;
    std::vector<std::unique_ptr<DynBNodeMatchesInfo>> Proposed;
    const std::vector<count> b;

    /**
     * Reads values from a file at @a path into the vector of b-values.
     *
     * @param size
     * @param path
     * @return std::vector<count>
     */
    std::vector<count> readBValuesFromFile(count size, const std::string &path) const;

    /**
     * Iterates up to @a b times over the heaviest neighbors of node @a u and makes
     * them to suitors if eligible.
     *
     * @param u
     */
    void findSuitors(node u);

    /**
     * Finds the heaviest unmatched neighbor that @a u has not yet proposed to
     * if it exists. For equally weighted edges w(u, t), w(u, v) and t < v, w(u, t) is
     * considered smaller than w(u, v) to break ties.
     *
     * @param y
     * @return DynBNode
     */
    DynBNode findPreferred(node u);

    /**
     * Makes @a v a suitor of @a u and recursively calls itself for previous worse
     * suitors of @a u that got replaced with their new best match.
     *
     * @param u
     * @param w
     * @param v
     */
    void makeSuitor(node u, edgeweight w, node v);

    /**
     * Checks the symmetry of pairs of nodes. It must hold that v is in suitors(u) iff u is
     * in suitors(v).
     *
     */
    bool isSymmetrical() const;
};
} // namespace NetworKit

#endif // NETWORKIT_MATCHING_B_SUITOR_MATCHER_HPP_
