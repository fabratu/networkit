/*
 * HypergraphLouvain.cpp
 * 
 * Created on: 07.2024
 *     Author: Isaline PLAID
*/

#include <networkit/community/HypergraphLouvain.hpp>
#include <networkit/graph/Hypergraph.hpp>

// Fast binomial coefficient calculation function
// k among n
double BinomialCoef_bis(const double n, const double k) {
  std::vector<double> aSolutions(k);
  aSolutions[0] = n - k + 1;

  for (double i = 1; i < k; ++i) {
    aSolutions[i] = aSolutions[i - 1] * (n - k + 1 + i) / (i + 1);
  }

  return aSolutions[k - 1];
}

namespace NetworKit {
HypergraphLouvain::HypergraphLouvain(const Hypergraph &graph, int iterations, bool randomize, double gamma, int type_contribution, double (*weightFun)(double))
    : CommunityDetectionAlgorithm(graph), gamma(gamma), numberOfIterations(iterations), random(randomize), type_contribution(type_contribution), weightFun (weightFun) {
    this->result = Partition(graph.numberOfNodes());
    this->result.allToSingletons();
}

void HypergraphLouvain::run() {

    handler.assureRunning();
    Hypergraph *currentGraph = const_cast<Hypergraph*>(G);
    
    //For loop to obtain the maxdge
    double s;
    d=0;
    (*currentGraph).forEdges([&](edgeid eid, edgeweight ew) {
        s = ((*currentGraph).order(eid));
        if (s>d){
            d=s;
        }
        if (ew != weightFun(s)){
          currentGraph->setEdgeWeight(eid,weightFun(s));
          DEBUG("weights on hyperedges that are set incorrectly, it's correct with ", currentGraph->getEdgeWeight(eid));
        }
    });

    //Initialization of volumes
    calculateVolumesHypergraph(*currentGraph);

    // communityVolumes_1 corresponds to the volume of our blocks
    // communityVolumes_2 corresponds to the volume of the communities we are building (by greedy move of blocks)
    communityVolumes_2 = communityVolumes_1; // Initialy we have only singleton comm, so communityVolumes_2 = communityVolumes_1
    //unweighted
    communityVolumes_2_unweighted = communityVolumes_1_unweighted;

    // Vector of sum of weight of edges of (i=0 to d)
    EdgeSizeWeight.resize(d+1);
    totalEdgeWeight=0.0;
    (*currentGraph).forEdges([&](edgeid eid, edgeweight ew) {
        EdgeSizeWeight[(*currentGraph).order(eid)] += ew;
        totalEdgeWeight += ew;
    });

    // At initialization each node is in a singleton community
    Partition zeta((*currentGraph).upperNodeIdBound());
    zeta.allToSingletons();
    zeta = result;

    // For each community, we obtain the set of neighboring communities
    NeighborComm.resize(result.upperBound());
    (*currentGraph).forNodes([&](node nid, nodeweight nw) {
      index c = zeta[nid];
      for (edgeid eid: (*currentGraph).edgesOf(nid)){
        for (node nid_prime: (*currentGraph).edgeMembers(eid)){
          NeighborComm[c].insert(zeta[nid_prime]);
        }
      }
    });

    // Main loop
    // greedy move + refinement phase until there are no more community changes, or we exceed the numberOfIterations
    for (int i = 0; i < numberOfIterations; ++i) {
        MoveHypergraph((*currentGraph), zeta);
        communityVolumes_1 = communityVolumes_2;
        communityVolumes_1_unweighted = communityVolumes_2_unweighted;
        
        // Stop the loop if greedy move and refinement phase do no community changes
        if (zeta.numberOfSubsets() == result.numberOfSubsets()) {
          break;
        }

        // Resetting neighbors for a new round of the loop
        zeta = result;
        (*currentGraph).forNodes([&](node nid, nodeweight nw) {
          index c = zeta[nid];
          for (edgeid eid: (*currentGraph).edgesOf(nid)){
            for (node nid_prime: (*currentGraph).edgeMembers(eid)){
              NeighborComm[c].insert(zeta[nid_prime]);
            }
          }
        });
    }

    hasRun = true;
}


// This function computates the volume of each community, as well as the volume of the graph
void HypergraphLouvain::calculateVolumesHypergraph(const Hypergraph &graph){
  auto timer = Aux::Timer();
  timer.start();
  communityVolumes_1.clear();
  communityVolumes_1.resize(result.upperBound());
  communityVolumes_1_unweighted.clear();
  communityVolumes_1_unweighted.resize(result.upperBound());
  GraphVolume=0.0;
  GraphVolume_unweighted=0.0;
  if (true) {
    std::vector<double> threadVolumes(omp_get_max_threads());
    std::vector<double> threadVolumes_unweighted(omp_get_max_threads());
    graph.forNodes([&](node nid, nodeweight nw) {
      {
      for (edgeid eid: graph.edgesOf(nid)){
        edgeweight ew =graph.getEdgeWeight(eid);
#pragma omp atomic
        communityVolumes_1[result[nid]] += ew;
        communityVolumes_1_unweighted[result[nid]] += 1;
        threadVolumes[omp_get_thread_num()] += ew;
        threadVolumes_unweighted[omp_get_thread_num()] += 1;
      }}
    });
    for (const auto vol : threadVolumes) {
      GraphVolume += vol;
    }
    for (const auto vol : threadVolumes_unweighted) {
      GraphVolume_unweighted += vol;
    }
  } 
  TRACE("Calculating Volumes took " + timer.elapsedTag());
}


// This function computates the gain of modularity of moving a set of nodes S from the community c to the community target_c
double HypergraphLouvain::deltaModHypergraph(const Hypergraph &graph, const Partition &zeta, index S, const Partition &p,  index c, index target_c){
  double modularityGain =0.0; // modularityGain = covGain - expCovGain
  double covGain =0.0;
  double expCovGain = 0.0;

  std::set<node> set_S = zeta.getMembers(S);
  std::set<edgeid> border_edge;
  for (node n : set_S){
    for (edgeid eid: graph.edgesOf(n)){
      border_edge.insert(eid);
    }
  }


  if (type_contribution==0 || type_contribution==10) {
    // >>>>> STRICT EDGE CONTRIBUTION

    // Let us compute covGain for strict edge contribution
    // Formulas :
    // covGain = \frac{In_w(S,c')-In_w(S,c)}{|E|_w} 
    // with $In_w(A,B) = \sum_{e \in E : e \nsubseteq A  \wedge  e \nsubseteq B  \wedge  e \subseteq A \cup B} w_e
    double in_c_S = 0.0;
    double in_target_c_S = 0.0;
    bool edgeBelongs_c=true;
    bool edgeBelongs_target_c=true;

    for (edgeid eid: border_edge){
      edgeBelongs_c =false;
      for (node nid: graph.edgeMembers(eid)){
        if (!edgeBelongs_c && zeta[nid]!=S){ // eid is not strictly included in S, so it's possible to take eid into account
          edgeBelongs_c = true;
        }
        if (p[nid] !=c){
          edgeBelongs_c = false;
          break; 
        }
      }
      if (edgeBelongs_c){
        in_c_S +=graph.getEdgeWeight(eid);
      }
                
      edgeBelongs_target_c=false;
      for (node nid: graph.edgeMembers(eid)){
        if (!edgeBelongs_target_c && zeta[nid]!=S){ // eid is not strictly included in S, so it's possible to take eid into account
          edgeBelongs_target_c = true;
        }
        if (p[nid] !=target_c && zeta[nid]!=S){
          edgeBelongs_target_c = false;
          break; 
        }
      }
      if (edgeBelongs_target_c){
        in_target_c_S += graph.getEdgeWeight(eid);
      }
    }
    covGain = (in_target_c_S - in_c_S) / totalEdgeWeight;

    if (type_contribution==0) {
      // FULLY WEIGHTED

      // The volume of S, community of S and target community. 
      double vol_S=communityVolumes_1[S];
      double vol_c = communityVolumes_2[c]; 
      double vol_target_c = communityVolumes_2[target_c];

      // Let us compute expCovGain for strict edge contribution
      // Formulas :
      // expCovGain = \sum_{d\geq 2} \frac{|E_d|_w}{|E|_w vol_w(V)^d}(vol_w(c'+S)^d -vol_w(c')^d -vol_w(c)^d + vol_w(c-S)^d)
      for (double j=0 ; j < d+1; j++){
        expCovGain += (EdgeSizeWeight[j] / (totalEdgeWeight * pow(GraphVolume,j))) * ( pow(vol_target_c + vol_S,j)-pow(vol_target_c,j)  - pow(vol_c,j) + pow(vol_c - vol_S,j));
      }
    }
    else{
      // PARTIALLY WEIGHTED

      // The volume of S, community of S and target community. 
      double vol_S=communityVolumes_1_unweighted[S];
      double vol_c = communityVolumes_2_unweighted[c]; 
      double vol_target_c = communityVolumes_2_unweighted[target_c];

      // Let us compute expCovGain for strict edge contribution
      // Formulas :
      // expCovGain = \sum_{d\geq 2} \frac{|E_d|_w}{|E|_w vol_w(V)^d}(vol_w(c'+S)^d -vol_w(c')^d -vol_w(c)^d + vol_w(c-S)^d)
      for (double j=0 ; j < d+1; j++){
        expCovGain += (EdgeSizeWeight[j] / (totalEdgeWeight * pow(GraphVolume_unweighted,j))) * ( pow(vol_target_c + vol_S,j)-pow(vol_target_c,j)  - pow(vol_c,j) + pow(vol_c - vol_S,j));
      }
    }
  }


  if (type_contribution==1 || type_contribution ==11) {// >>>>> MAJORITY EDGE CONTRIBUTION

    // Let us compute covGain for majority edge contribution
    // Formulas :
    // covGain = \frac{In(S,c')-In(S,c)}{|E|} 
    // with In_w(A,B) = \sum_{e \in E ~:~ e \nsubseteq_{maj} A ~ \wedge ~ e \nsubseteq_{maj} B ~ \wedge ~ e \subseteq_{maj} A \cup B} w_e
    double in_c_S = 0.0;
    double in_target_c_S = 0.0;
    bool edgeBelongs_c=true;
    bool edgeBelongs_target_c=true;

    for (edgeid eid: border_edge){
      edgeBelongs_c =false;
      edgeBelongs_target_c=false;
      std::vector<double> CommEdge(3, 0.0);
      for (node nid: graph.edgeMembers(eid)){
        if (p[nid]==c){CommEdge[0]++;}
        if (p[nid]==target_c){CommEdge[1]++;}
        if (zeta[nid]==S){CommEdge[2]++;}
      }
                
      if (!(CommEdge[1] >= int(graph.order(eid)/ 2. +1)) && !(CommEdge[0]- CommEdge[2] >= int(graph.order(eid)/ 2. +1)) && !(CommEdge[2] >= int(graph.order(eid)/ 2. +1))){
        if (CommEdge[0] >= int(graph.order(eid)/ 2. +1)){
          edgeBelongs_c =true;
        }
        if (CommEdge[1] + CommEdge[2] >= int(graph.order(eid)/ 2. +1)){
          edgeBelongs_target_c =true;
        }
      }
            
      if (edgeBelongs_c){
        in_c_S += graph.getEdgeWeight(eid);
      }

      if (edgeBelongs_target_c){
        in_target_c_S += graph.getEdgeWeight(eid);
      }
    }
    covGain = (in_target_c_S - in_c_S) / totalEdgeWeight;


    if (type_contribution==1) {
      // FULLY WEIGHTED

      // The volume of S, community of S and target community. 
      double vol_S=communityVolumes_1[S];
      double vol_c = communityVolumes_2[c]; 
      double vol_target_c = communityVolumes_2[target_c];

      // Let us compute expCovGain for majority edge contribution
      // Formulas :
      // expCovGain = \sum_{d\geq 2} \frac{|E_d|}{|E| vol(V)^d} \sum_{i=\frac{d}{2}+1}^{d} \binom{d}{i} (vol(c'+S)^i(vol(V)-vol(c'+S))^{d-i} -vol(c')^i(vol(V)-vol(c'))^{d-i} -vol(c)^i(vol(V)-vol(c))^{d-i} + vol(c-S)^i(vol(V)-vol(c-S))^{d-i})
      double sum = 0.0;
      for (double j=0 ; j < d+1; j++){
        sum = 0.0;
        for (double i= int ((j / 2.) + 1); i <= j; i++){
          sum += BinomialCoef_bis(j,i) *(pow((vol_target_c + vol_S),i)*pow((GraphVolume - vol_target_c - vol_S),j-i) -  pow(vol_target_c,i)*pow((GraphVolume - vol_target_c),j-i) - pow(vol_c ,i)*pow((GraphVolume- vol_c),j-i) + pow((vol_c - vol_S),i)*pow((GraphVolume - vol_c + vol_S),j-i));
        }
        expCovGain += (EdgeSizeWeight[j] / (totalEdgeWeight * pow(GraphVolume,j))) * sum;
      }
    } 
    else {
      // PARTIALLY WEIGHTED

      // The volume of S, community of S and target community. 
      double vol_S=communityVolumes_1_unweighted[S];
      double vol_c = communityVolumes_2_unweighted[c]; 
      double vol_target_c = communityVolumes_2_unweighted[target_c];

      // Let us compute expCovGain for majority edge contribution
      // Formulas :
      // expCovGain = \sum_{d\geq 2} \frac{|E_d|}{|E| vol(V)^d} \sum_{i=\frac{d}{2}+1}^{d} \binom{d}{i} (vol(c'+S)^i(vol(V)-vol(c'+S))^{d-i} -vol(c')^i(vol(V)-vol(c'))^{d-i} -vol(c)^i(vol(V)-vol(c))^{d-i} + vol(c-S)^i(vol(V)-vol(c-S))^{d-i})
      double sum = 0.0;
      for (double j=0 ; j < d+1; j++){
        sum = 0.0;
        for (double i= int ((j / 2.) + 1); i <= j; i++){
          sum += BinomialCoef_bis(j,i) *(pow((vol_target_c + vol_S),i)*pow((GraphVolume_unweighted - vol_target_c - vol_S),j-i) -  pow(vol_target_c,i)*pow((GraphVolume_unweighted - vol_target_c),j-i) - pow(vol_c ,i)*pow((GraphVolume_unweighted- vol_c),j-i) + pow((vol_c - vol_S),i)*pow((GraphVolume_unweighted - vol_c + vol_S),j-i));
        }
        expCovGain += (EdgeSizeWeight[j] / (totalEdgeWeight * pow(GraphVolume_unweighted,j))) * sum;
      }
    }
  } 
  
  modularityGain= covGain - gamma * expCovGain; 

  return modularityGain;
}


// Greedy Move Phase
void HypergraphLouvain::MoveHypergraph(const Hypergraph &graph, const Partition &zeta){
  // The zeta partition corresponds to our block of nodes that we will move in this phase of greedy move
  // Initialy zeta = result 

  count moved = 0;
  count totalNodes = 0;

  std::vector<bool> inQueue(zeta.upperBound(), false); // Boolean array allowing you not to add the same node twice in the queue
  std::queue<index> queue; // nodes that remain to be processed

  //uint64_t vectorSize = communityVolumes_1.capacity();
  int upperBound = result.upperBound();


  std::vector<index> currentBlocks; // block to cover during a loop
  std::vector<index> newNodes;
  newNodes.reserve(WORKING_SIZE);

  // We add all node in queue 
  std::unique_ptr<std::atomic<bool>[]> exists(new std::atomic<bool>[zeta.upperBound()] {});
  zeta.parallelForEntries([&](index, index s) {
    if (s != none) {
      exists[s] = true;
    }
  });
  count k = 0; // number of actually existing clusters (singleton comm)
  for (index i = 0; i < static_cast<index>(zeta.upperBound()); ++i) {
    if (exists[i]) {
      k++;
      currentBlocks.push_back(i);
      inQueue[i] = true;
    }
  }

  // We randomize the order of the blocks
  if (random) {
    auto &mt = Aux::Random::getURNG();
    std::shuffle(currentBlocks.begin(), currentBlocks.end(), mt);
  }

  // Main Loop : Greedy Move
  do {
    handler.assureRunning();
    for (index u : currentBlocks) {

      //Check u not empty
      std::set<node> S= zeta.getMembers(u);
      if (S.empty()){
        continue;
      }
      index current_comm = result[*(S.begin())];

      if (!inQueue[u]){
        continue;
      } // index in currentBlocks must be in the queue

      double maxDelta = std::numeric_limits<double>::lowest();
      index bestCommunity = none;
      bool exit_best_com = false; 

      if (NeighborComm[u].empty()){
        continue;
      }

      // We test all the neighboring communities, we look for the best gain in modularity
      for (index neighborCommunity : NeighborComm[u]){
        node no = zeta.giveOne(neighborCommunity);
        if (no == std::numeric_limits<uint64_t>::max()){
          continue;
        }
        index target_comm = result[no];
        if (target_comm != current_comm) {
          double delta = deltaModHypergraph(graph, zeta, u, result, current_comm, target_comm);
          if (delta > maxDelta) {
            maxDelta = delta;
            bestCommunity = target_comm;
            exit_best_com = true;
          }
        }
      }
      

      // We move the block u into the best community found
      if (exit_best_com){
      /*Aux::Log::setLogLevel("DEBUG");
      INFO("move u: ", u); 
      INFO("delta: ", maxDelta);
      INFO("for: ", bestCommunity);*/
      if (maxDelta <= 0){
        bestCommunity=current_comm;
      }
      else {
        moved++;
      }

      
      for (node nid : zeta.getMembers(u)){
        result[nid] = bestCommunity;
      }
      communityVolumes_2[bestCommunity] += communityVolumes_1[u];
      communityVolumes_2[current_comm]-= communityVolumes_1[u];
      communityVolumes_2_unweighted[bestCommunity] += communityVolumes_1_unweighted[u];
      communityVolumes_2_unweighted[current_comm]-= communityVolumes_1_unweighted[u];
      changed = true;
      inQueue[u] = false;

      // We add in the queue the nodes which must be explored again
      if (bestCommunity!=current_comm){
      for (index neighborCommunity : NeighborComm[u]){
        node no = zeta.giveOne(neighborCommunity);
        if (no == std::numeric_limits<uint64_t>::max()){
          continue;
        }
        index target_comm = result[no];
        if (target_comm != bestCommunity &&  neighborCommunity !=u) {
          if (!inQueue[neighborCommunity]) {
            INFO("  ", neighborCommunity);
            newNodes.push_back(neighborCommunity);
            inQueue[neighborCommunity] = true;
            if (newNodes.size() == WORKING_SIZE) {
              for (node v : newNodes) {
                queue.push(v);
              }
              newNodes.clear();
              newNodes.reserve(WORKING_SIZE);
            }
          }
        }
      } } }
    }

    totalNodes += currentBlocks.size();
    if (!newNodes.empty()) {
      currentBlocks = std::move(newNodes);
      newNodes.clear();
    } else if (!queue.empty()) {
      currentBlocks.clear();
      while (!queue.empty()) {
        currentBlocks.push_back(queue.front());
        queue.pop();
      }
    } else {
      break;
    }
    } while (true);

    result.setUpperBound(upperBound);
    if (Aux::Log::isLogLevelEnabled(Aux::Log::LogLevel::DEBUG)) {
        DEBUG("Total worked: ", totalNodes, " Total moved: ", moved);
    }
  
}

} // namespace NetworKit