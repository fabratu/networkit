#include <networkit/matching/BSuitorMatcher.hpp>

#include <algorithm>
#include <deque>
#include <queue>
#include <unordered_set>
#include <iomanip>
#include <stdexcept>
#include <random>
#include <mutex>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/Parallelism.hpp>

namespace NetworKit {
BSuitorMatcher::BSuitorMatcher(const Graph &G, const std::vector<count> &b) : BMatcher(G, b), b(b) {
    if (G.numberOfSelfLoops() > 0)
        throw std::runtime_error("This algorithm does not support graphs with self-loops.");

    if (G.isDirected())
        throw std::runtime_error("This algorithm does not support directed graphs.");

    if (b.size() != G.numberOfNodes())
        throw std::runtime_error(
            "The number of b values does not match the number of nodes in this graph.");

    const auto n = G.upperNodeIdBound();
    Suitors.reserve(n);
    Proposed.reserve(n);
    for (index i = 0; i < n; i++) {
        Suitors.emplace_back(std::make_unique<DynBNodeMatchesInfo>(b.at(i)));
        Proposed.emplace_back(std::make_unique<DynBNodeMatchesInfo>(b.at(i)));
    }
}

BSuitorMatcher::BSuitorMatcher(const Graph &G, count b)
    : BSuitorMatcher(G, std::vector<count>(G.numberOfNodes(), b)) {}

BSuitorMatcher::BSuitorMatcher(const Graph &G, const std::string &path)
    : BSuitorMatcher(G, readBValuesFromFile(G.numberOfNodes(), path)) {}

std::vector<count> BSuitorMatcher::readBValuesFromFile(count size, const std::string &path) const {
    std::vector<count> b;
    b.reserve(size);
    std::ifstream file(path);
    std::string line;
    int line_number = 1;

    while (std::getline(file, line)) {
        std::istringstream istring(line);
        int val;
        if (!(istring >> val)) {
            throw std::runtime_error("File " + path + " contains an invalid value in line "
                                     + std::to_string(line_number) + ".");
        }
        if (istring >> val) {
            throw std::runtime_error("File " + path + " contains multiple values in line "
                                     + std::to_string(line_number) + ".");
        }
        if (val < 0) {
            throw std::runtime_error("File " + path + " contains a negative value in line "
                                     + std::to_string(line_number) + ".");
        }
        b.emplace_back(val);
        line_number++;
    }
    if (b.size() != size) {
        throw std::runtime_error("The number of values in file " + path
                                 + " does not match the number of nodes in this graph.");
    }
    return b;
}

void BSuitorMatcher::runSequential() {
    G->forNodes([&](node u) { findSuitors(u); });
    hasRun = true;
}

void BSuitorMatcher::findSuitors(node cur) {
    for (index i = 0; i < b.at(cur); i++) {
        auto [pref, heaviest] = findPreferred(cur);
        INFO("Got: ", pref);
        if (pref != none) {
            makeSuitor(cur, heaviest, pref);
        }
    }
}

DynBNode BSuitorMatcher::findPreferred(node u) {
    DynBNode best = DynBNode{none, 0};

    // OLD CODE
    // auto hasProposedTo = [&](node x) -> bool {
    //     return std::any_of(Proposed.at(u)->partners.begin(), Proposed.at(u)->partners.end(),
    //                        [x](const DynBNode &y) { return y.id == x; });
    // };

    for (auto n : G->weightNeighborRange(u)) {
        const DynBNode v = DynBNode(n.first, n.second);
        if (!Proposed.at(u)->hasPartner(v)) {
            // INFO("Checking: ", u , ",", v.id, " weight: ", v.weight);
        // if (!hasProposedTo(v.id)) {
            if (v.weight > best.weight || (v.weight == best.weight && v.id < best.id)) {
                const auto n_suitor_weight = Suitors.at(v.id)->min.weight;

                if (v.weight > n_suitor_weight
                    || (v.weight == n_suitor_weight && u < Suitors.at(v.id)->min.id)) {
                    best = v;
                }
            }
        }
    }
    return best;
}

void BSuitorMatcher::makeSuitor(node u, edgeweight w, node v) {
    // auto smallest = Suitors.at(v)->popMinIfFull();
    auto smallest = Suitors.at(v)->insert(DynBNode(u, w), false);
    Proposed.at(u)->insert(DynBNode(v, w), false);
    INFO("Inserted: ", u, ",", v);

    if (smallest.id != none) {
        Proposed.at(smallest.id)->remove(DynBNode{v, smallest.weight});
        auto [pref, heaviest] = findPreferred(smallest.id);
        if (pref != none) {
            makeSuitor(smallest.id, heaviest, pref);
        }
    }
}

void BSuitorMatcher::runParallel() {
    // Note: this is wrong, since technically we also have to lock cur due to search of proposed


    INFO("Starting in parallel");

    std::vector<std::mutex> nodeLocks(G->numberOfNodes());

    auto numThreads = Aux::getMaxNumberOfThreads();

    std::vector<std::unordered_set<node>> Q(numThreads);
    // Q.reserve(numThreads);

    int chunkSize = std::ceil(static_cast<double>(G->numberOfNodes()) / numThreads);
    for (size_t i = 0; i < numThreads; i++)
    {  
        Q[i].reserve(chunkSize);
    }
    
    G->balancedParallelForNodes([&](node u) {
        int threadNum = omp_get_thread_num();
        Q[threadNum].insert(u);
    });

#pragma omp parallel
    {
        // std::random_device rd;
        // std::mt19937 gen {rd()};
        int threadNum = omp_get_thread_num();
        std::unordered_set<node> *localQ = &Q[threadNum];
        // std::deque<node> *localQ = &Q[threadNum];
        // std::ranges::shuffle(*localQ, gen);
        // std::deque<node> QPrime{};
        std::unordered_set<node> QPrime;

        // size_t localNodes = localQ.size();
        do {
            if (localQ->size() == 0) {
                break;
            }

            for(auto it = localQ->begin(); it != localQ->end(); ++it) {

                node cur = *it;
                // localQ->pop_back();
                size_t i = 1;

                bool isExhausted = false;

                INFO("Thread: ", threadNum, " processing node ", cur);
                while (i <= b.at(cur) && !isExhausted) {
                    nodeLocks[cur].lock();
                    DynBNode p = findPreferred(cur);
                    nodeLocks[cur].unlock();

                    if (p.id != none) {
                        INFO("Thread: ", threadNum, " node ", cur, " found a potential suitor: ", p.id);
                        nodeLocks[p.id].lock();
                        nodeLocks[cur].lock();
                        DynBNode newP = findPreferred(cur);
                        nodeLocks[cur].unlock();

                        if(newP == p) {
                            i++;
                            INFO("Thread: ", threadNum, " node ", cur, " locked and inserted: ", p.id);
                            DynBNode prevP = Suitors.at(p.id)->insert({cur,G->weight(cur, p.id)});
                            Proposed.at(cur)->insert(p, false);
                            if (prevP.id != none) {
                                INFO("Thread: ", threadNum, " node ", cur, " inserted: ", p.id, " and got a loose end ", prevP.id);
                                nodeLocks[prevP.id].lock();
                                if(Proposed.at(prevP.id)->hasPartner(p)) {
                                    INFO("Thread: ", threadNum, " node ", cur, " remove ", p.id, " from loose end ", prevP.id);
                                    Proposed.at(prevP.id)->remove(p);
                                }
                                nodeLocks[prevP.id].unlock();
                                INFO("Thread: ", threadNum, " node ", cur, " unlocked loose end ", prevP.id);
                                QPrime.insert(prevP.id);
                                INFO("Thread: ", threadNum, " node ", cur, " pushed: ", prevP.id, " to Q");
                            }
                        }
                        nodeLocks[p.id].unlock();
                        INFO("Thread: ", threadNum, " node ", cur, " unlocked node ", p.id);
                    } else {
                        isExhausted = true;
                    }

                }
            }

            localQ->clear();

            for( node newNode : QPrime){
                localQ->insert(newNode);
            }
            QPrime.clear();

            // localNodes = localQ.size();
        } while(true);        


    }

     hasRun = true;

}

bool BSuitorMatcher::isSymmetrical() const {
    bool sym = true;
    auto matchedSymmetrical = [&](DynBNode x, DynBNode y) -> bool {
        return Suitors.at(x.id)->hasPartner(y) == Suitors.at(y.id)->hasPartner(x);
    };

    G->forNodes([&](node u) {
        G->forNodes([&](node v) {
            edgeweight weight = G->weight(u,v);
            if (u > v && !matchedSymmetrical(DynBNode{u,weight}, DynBNode{v,weight})) {
                sym = false;
            }
        });
    });
    return sym;
}

void BSuitorMatcher::buildBMatching() {
    // TODO make parallel
    M.reset();
    G->forNodes([&](node x) {
        // assert(Suitors.at(x)->partners.size() <= b.at(x));
        for (auto y : Suitors.at(x)->partnersBucket) {
            // std::cout << x << ": " << y.first << "," << y.second.id << std::endl;
        // for (auto y : Suitors.at(x)->partners) {
            if (y.id != none && x < y.id) {
                // std::cout << "Match: " << x << " and " << y.second.id << std::endl;
                M.match(x, y.id);
            }
        }
    });
}

void BSuitorMatcher::run() {
    if (parallel) {
        runParallel();
    } else {
        runSequential();
    }
}
} // namespace NetworKit
