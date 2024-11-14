/*
 * Modularity.cpp
 *
 *  Created on: 10.12.2012
 *      Author: Christian Staudt
 */

#include <cmath>
#include <stdexcept>

#include <networkit/auxiliary/Log.hpp>
#include <networkit/community/Coverage.hpp>
#include <networkit/community/Modularity.hpp>

namespace NetworKit {

Modularity::Modularity() : QualityMeasure(), gTotalEdgeWeight(0.0) {}

void Modularity::setTotalEdgeWeight(double totalEdgeWeight) {
    gTotalEdgeWeight = totalEdgeWeight;
}

double Modularity::getQuality(const Partition &zeta, const Graph &G) {
    assert(G.numberOfNodes() <= zeta.numberOfElements());

    Coverage coverage;
    // deprecated: intraEdgeWeightSum / gTotalEdgeWeight;
    double cov = coverage.getQuality(zeta, G);
    //// term $\frac{ \sum_{C \in \zeta}( \sum_{v \in C}
    ///\omega(v) )^2 }{4( \sum_{e \in E} \omega(e) )^2 }$
    double expCov;
    double modularity; // mod = coverage - expected coverage
    if (gTotalEdgeWeight == 0.0) {
        gTotalEdgeWeight = G.totalEdgeWeight(); // compute total edge weight in G
        DEBUG("computed total edge weight: ", gTotalEdgeWeight);
    }

    if (gTotalEdgeWeight == 0.0) {
        ERROR("G: m=", G.numberOfEdges(), "n=", G.numberOfNodes());
        throw std::invalid_argument(
            "Modularity is undefined for graphs without edges (including self-loops).");
    }

    //!< cluster -> sum of the weights of incident edges for all nodes
    std::vector<double> incidentWeightSum(zeta.upperBound(), 0.0);

    // compute volume of each cluster
    G.parallelForNodes([&](node v) {
        // add to cluster weight
        index c = zeta[v];
        assert(zeta.lowerBound() <= c);
        assert(c < zeta.upperBound());
        // account for self-loops a second time
#pragma omp atomic
        incidentWeightSum[c] += G.weightedDegree(v) + G.weight(v, v);
    });

    // compute sum of squared cluster volumes and divide by squared graph volume
    expCov = 0.0;

#pragma omp parallel for reduction(+ : expCov)
    for (omp_index c = static_cast<omp_index>(zeta.lowerBound());
         c < static_cast<omp_index>(zeta.upperBound()); ++c) {
        // squared
        expCov +=
            ((incidentWeightSum[c] / gTotalEdgeWeight) * (incidentWeightSum[c] / gTotalEdgeWeight))
            / 4;
    }

    DEBUG("expected coverage: ", expCov);

    // assert ranges of coverage
    assert(cov <= 1.0);
    assert(cov >= 0.0);
    assert(expCov <= 1.0);
    assert(expCov >= 0.0);

    modularity = cov - expCov;
    DEBUG("modularity = ", modularity);

    // reset totalEdgeWeight
    gTotalEdgeWeight = 0.0;

    assert(!std::isnan(modularity)); // do not return NaN
    // do not return anything not in the range of modularity values
    assert(modularity >= -0.5);
    assert(modularity <= 1);
    return modularity;
}

// TODO check if binom exist in networkit
// Fast binomial coefficient calculation function
// k among n
double BinomialCoefficient(const double n, const double k) {
  std::vector<double> aSolutions(k);
  aSolutions[0] = n - k + 1;

  for (double i = 1; i < k; ++i) {
    aSolutions[i] = aSolutions[i - 1] * (n - k + 1 + i) / (i + 1);
  }

  return aSolutions[k - 1];
}

double Modularity::getQualityHypergraph(const Partition &zeta, const Hypergraph &G, double gamma, int type_contribution) {
    assert(G.numberOfNodes() <= zeta.numberOfElements());
    double cov = 0.0;
    double expCov = 0.0;
    double modularity =0.0; // mod = coverage - expected coverage

    //For loop to obtain the max size edge
    double d= 0;
    G.forEdges([&](edgeid eid, edgeweight ew){
        if (G.order(eid)>d){
            d=G.order(eid);
        }
    });

    if (type_contribution==0) {// >>>>> STRICT EDGE CONTRIBUTION // FULLY WEIGHTED

        // >>> Firstly, we compute the coverage : 
        double intraEdgeWeightSum = 0.0; // sum of all intra weight inside each comminuty
        double totalEdgeWeight = 0.0; // add all edge weights
        std::vector<double> EdgeSizeWeight(d+1, 0.0); // vector of sum of weight of edges of size i (i=0 to d)

#ifdef DEBUG
        if (G.order(0) > 0) {
            ERROR("First edge error");
        }
#endif
        // compute intra-cluster edge weights per cluster
        bool isFirst=true;
        bool edgeBelongs=true;
        index c, c_prime;
        G.forEdges([&](edgeid eid, edgeweight ew) {
            isFirst=true;
            edgeBelongs=true;
            // an edge belongs to a community if all nodes in this edge belong to the same community 
            for (auto nid: G.edgeMembers(eid)) {
                if(isFirst){
                    c = zeta[nid];
                    isFirst=false;
                }
                c_prime = zeta[nid];
                if (c!=c_prime){
                    edgeBelongs=false;
                    break;
                }
            }
            if (edgeBelongs){
                intraEdgeWeightSum += ew;
            }
            EdgeSizeWeight[G.order(eid)] += ew;
            totalEdgeWeight+= ew;
        });
        cov = intraEdgeWeightSum / totalEdgeWeight;

        // >>> Secondly, we compute the expected coverage : 
        std::set<std::set<index>> subsets =zeta.getSubsets(); //obtain all partition 
        std::vector<double> intraVolume(zeta.upperBound(), 0.0);
        double volume=0.0;

        // Compute volume of each community 
        double i=0;
        for (std::set c : subsets){
            for (index n : c){
                for (edgeid eid: G.edgesOf(n)){
                    intraVolume[i]+=G.getEdgeWeight(eid);
                }
            }
            i++;
        }

        //Compute volume of all hypergraph
        G.forNodes([&](node nid, nodeweight nw) {
            for (edgeid eid: G.edgesOf(nid)){
                volume+=G.getEdgeWeight(eid);
            }
        });

        double sum_vol =0.0;
        // Compute expCov : \sum_{d} (\frac{|E_d|}{vol(V)^d} \sum_{c} vol(c)^d)
        for (double j=0 ; j < d+1; j++){
            sum_vol=0.0;
            for (i=0; i < zeta.upperBound(); i++){
                sum_vol += pow(intraVolume[i],j);
            }
            expCov += (EdgeSizeWeight[j] / pow(volume,j)) * sum_vol;
        }
        
        expCov = expCov / totalEdgeWeight;
    }// End if for strict edge contribution : type_contribution == 0 



    if (type_contribution==1) {// >>>>> MAJORITY EDGE CONTRIBUTION // FULLY WEIGHTED
        
        //>>> Firstly, we compute the coverage : 

        double intraEdgeWeightSum = 0.0; // sum of all intra weight inside each comminuty
        double totalEdgeWeight = 0.0; // add all edge weights

        std::vector<double> EdgeSizeWeight(d+1, 0.0); // vector of sum of weight of edges of size i (i=0 to d)

        // compute intra-cluster edge weights per cluster
        bool isFirst=true;
        bool edgeBelongs=true;
        index c, c_prime;
        G.forEdges([&](edgeid eid, edgeweight ew) {
            isFirst=true;
            edgeBelongs=false;
            // an edge belongs to a community if the majority of its nodes belongs to the same community 
            std::vector<double> CommEdge(zeta.upperBound(), 0.0);
            for (auto nid: G.edgeMembers(eid)) {
                CommEdge[zeta[nid]]++;
            }
            for (auto v: CommEdge) {
                if (v >= int(G.order(eid)/ 2. +1)){
                    edgeBelongs=true;
                }
            }
            if (edgeBelongs){
                intraEdgeWeightSum += ew;
            }
            EdgeSizeWeight[G.order(eid)] += ew;
            totalEdgeWeight+= ew;
        });

        cov = intraEdgeWeightSum / totalEdgeWeight;


        //>>> Secondly, we compute the expected coverage : 
        std::set<std::set<index>> subsets =zeta.getSubsets(); //obtain all partition 
        std::vector<double> intraVolume(zeta.upperBound(), 0.0);
        double volume=0.0;

        // Compute volume of each community 
        double i=0;
        for (std::set c : subsets){
            for (index n : c){
                for (edgeid eid: G.edgesOf(n)){
                    intraVolume[i]+=G.getEdgeWeight(eid);
                }
            }
            i++;
        }

        //Compute volume of all hypergraph
        G.forNodes([&](node nid, nodeweight nw) {
            for (edgeid eid: G.edgesOf(nid)){
                volume+=G.getEdgeWeight(eid);
            }
        });

        double sum_vol =0.0;
        double sum =0.0;
        // Compute expCov : \sum_d \frac{|E_d|}{vol(V)^d} \sum_{i= \frac{d}{2} +1}^{d}  \binom{d}{i} \sum_{c \in C} (vol(c))^i (vol(V)- vol(c))^{d-i}
        for (double j=0 ; j < d+1; j++){
            sum=0.0;
            for (i= int ((j / 2.) + 1); i <= j; i++){
                sum_vol=0.0;
                for (double c=0; c < zeta.upperBound(); c++){
                    sum_vol += pow(intraVolume[c],i)* pow(volume - intraVolume[c], j-i);
                }
                sum_vol= BinomialCoefficient(j,i)* sum_vol;
                sum += sum_vol;
            }
            expCov += (EdgeSizeWeight[j] / pow(volume,j)) * sum;
        }

        expCov = expCov / totalEdgeWeight;
    }// End if for majority edge contribution, fully weighted: type_contribution == 01

    
    
    if (type_contribution==10) {// >>>>> STRICT EDGE CONTRIBUTION // PARTIALLY WEIGHTED

        // >>> Firstly, we compute the coverage : 
        double intraEdgeWeightSum = 0.0; // sum of all intra weight inside each comminuty
        double totalEdgeWeight = 0.0; // add all edge weights
        std::vector<double> EdgeSizeWeight(d+1, 0.0); // vector of sum of number of edges of size i (i=0 to d)

        // compute intra-cluster edge weights per cluster
        bool isFirst=true;
        bool edgeBelongs=true;
        double nbEdge = 0.0;
        index c, c_prime;
        G.forEdges([&](edgeid eid, edgeweight ew) {
            nbEdge ++;
            isFirst=true;
            edgeBelongs=true;
            // an edge belongs to a community if all nodes in this edge belong to the same community 
            for (auto nid: G.edgeMembers(eid)) {
                if(isFirst){
                    c = zeta[nid];
                    isFirst=false;
                }
                c_prime = zeta[nid];
                if (c!=c_prime){
                    edgeBelongs=false;
                    break;
                }
            }
            if (edgeBelongs){
                intraEdgeWeightSum += ew;
            }
            EdgeSizeWeight[G.edgeMembers(eid).size()] += ew;
            totalEdgeWeight+= ew;
        });
        cov = intraEdgeWeightSum / totalEdgeWeight;

        // >>> Secondly, we compute the expected coverage : 
        std::set<std::set<index>> subsets =zeta.getSubsets(); //obtain all partition 
        std::vector<double> intraVolume(zeta.upperBound(), 0.0);
        double volume=0.0;

        // Compute volume of each community 
        double i=0;
        for (std::set c : subsets){
            for (index n : c){
                for (edgeid eid: G.edgesOf(n)){
                    intraVolume[i]+= 1;
                }
            }
            i++;
        }

        //Compute volume of all hypergraph
        G.forNodes([&](node nid, nodeweight nw) {
            for (edgeid eid: G.edgesOf(nid)){
                volume+=1;
            }
        });

        double sum_vol =0.0;
        // Compute expCov : \sum_{d} (\frac{|E_d|}{vol(V)^d} \sum_{c} vol(c)^d)
        for (double j=0 ; j < d+1; j++){
            sum_vol=0.0;
            for (i=0; i < zeta.upperBound(); i++){
                sum_vol += pow(intraVolume[i],j);
            }
            expCov += (EdgeSizeWeight[j] / pow(volume,j)) * sum_vol;
        }
        
        expCov = expCov / totalEdgeWeight;
    }// End if for strict edge contribution, partially weighted: type_contribution == 10



    if (type_contribution==11) {// >>>>> MAJORITY EDGE CONTRIBUTION // PARTIALLY WEIGHTED
        
        //>>> Firstly, we compute the coverage : 

        double intraEdgeWeightSum = 0.0; // sum of all intra weight inside each comminuty
        double totalEdgeWeight = 0.0; // add all edge weights
        double nbEdge = 0.0;
        std::vector<double> EdgeSizeWeight(d+1, 0.0); // vector of sum of weight of edges of size i (i=0 to d)

        // compute intra-cluster edge weights per cluster
        bool isFirst=true;
        bool edgeBelongs=true;
        index c, c_prime;
        G.forEdges([&](edgeid eid, edgeweight ew) {
            nbEdge ++;
            isFirst=true;
            edgeBelongs=false;
            // an edge belongs to a community if the majority of its nodes belongs to the same community 
            std::vector<double> CommEdge(zeta.upperBound(), 0.0);
            for (auto nid: G.edgeMembers(eid)) {
                CommEdge[zeta[nid]]++;
            }
            for (auto v: CommEdge) {
                if (v >= int(G.order(eid)/ 2. +1)){
                    edgeBelongs=true;
                }
            }
            if (edgeBelongs){
                intraEdgeWeightSum += ew;
            }
            EdgeSizeWeight[G.order(eid)] += ew;
            totalEdgeWeight+= ew;
        });

        cov = intraEdgeWeightSum / totalEdgeWeight;


        //>>> Secondly, we compute the expected coverage : 
        std::set<std::set<index>> subsets =zeta.getSubsets(); //obtain all partition 
        std::vector<double> intraVolume(zeta.upperBound(), 0.0);
        double volume=0.0;

        // Compute volume of each community 
        double i=0;
        for (std::set c : subsets){
            for (index n : c){
                for (edgeid eid: G.edgesOf(n)){
                    intraVolume[i]+=1;
                }
            }
            i++;
        }

        //Compute volume of all hypergraph
        G.forNodes([&](node nid, nodeweight nw) {
            for (edgeid eid: G.edgesOf(nid)){
                volume+=1;
            }
        });

        double sum_vol =0.0;
        double sum =0.0;
        // Compute expCov : \sum_d \frac{|E_d|}{vol(V)^d} \sum_{i= \frac{d}{2} +1}^{d}  \binom{d}{i} \sum_{c \in C} (vol(c))^i (vol(V)- vol(c))^{d-i}
        for (double j=0 ; j < d+1; j++){
            sum=0.0;
            for (i= int ((j / 2.) + 1); i <= j; i++){
                sum_vol=0.0;
                for (double c=0; c < zeta.upperBound(); c++){
                    sum_vol += pow(intraVolume[c],i)* pow(volume - intraVolume[c], j-i);
                }
                sum_vol= BinomialCoefficient(j,i)* sum_vol;
                sum += sum_vol;
            }
            expCov += (EdgeSizeWeight[j] / pow(volume,j)) * sum;
        }

        expCov = expCov / totalEdgeWeight;
    }// End if for majority edge contribution, partially weighted: type_contribution == 11


    modularity =  cov - gamma * expCov;
    return modularity;
}



double Modularity::deltaModularityHypergraph(const Partition &zeta, const Hypergraph &G, std::set<node> S, index target_c,  double gamma, int type_contribution) {
    double modularityGain =0.0; // modularityGain = covGain - expCovGain
    double covGain =0.0;
    double expCovGain = 0.0;

    //Compute volume of all hypergraph
    double volume = 0.0;
    double volume_unweight = 0.0;
    G.forNodes([&](node nid, nodeweight nw) {
        for (edgeid eid: G.edgesOf(nid)){
            volume+=G.getEdgeWeight(eid);
            volume_unweight += 1;
        }
    });

    //For loop to obtain the max size edge
    double d= 0;
    for (edgeid eid = 0; eid < G.upperEdgeIdBound(); ++eid) {
        if (G.order(eid)>d){
            d=G.order(eid);
        }
    }

    std::vector<double> EdgeSizeWeight(d+1, 0.0); // vector of sum of weight of edges of size i (i=0 to d)
    double totalEdgeWeight = 0.0; // add all edge weights
    G.forEdges([&](edgeid eid, edgeweight ew) {
        EdgeSizeWeight[G.order(eid)] += ew;
        totalEdgeWeight += ew;
    });

    //To calculate expCovGain, I need 3 different volumes
    // The volume of the old community c of nodes of S including the nodes of S   (All S nodes belong to the same community)
    index c= zeta[*S.begin()];
    double vol_c=0.0;
    double vol_c_unweight=0.0;
    std::set<node> set_c;
    set_c = zeta.getMembers(c);
    for (node nid : set_c){
        for (edgeid eid: G.edgesOf(nid)){
            vol_c +=G.getEdgeWeight(eid);
            vol_c_unweight +=1;
        }
    }

    //The volume of S
    double vol_S=0.0;
    double vol_S_unweight=0.0;
    for (node nid : S){
        for (edgeid eid: G.edgesOf(nid)){
            vol_S +=G.getEdgeWeight(eid);
            vol_S_unweight+=1;
        }
    }

    //The volume of target community without the nodes of S
    double vol_target_c=0.0;
    double vol_target_c_unweight=0.0;
    std::set<node> set_target_c;
    set_target_c = zeta.getMembers(target_c);
    for (node nid : set_target_c){
        for (edgeid eid: G.edgesOf(nid)){
            vol_target_c +=G.getEdgeWeight(eid);
            vol_target_c_unweight+=1;
        }
    }


    if (type_contribution==0) {// >>>>> STRICT EDGE CONTRIBUTION
        //covGain = \frac{In_w(S,c')-In_w(S,c)}{|E|_w} 
        //with $In_w(A,B) = \sum_{e \in E ~:~ e \nsubseteq A ~ \wedge ~ e \nsubseteq B ~ \wedge ~ e \subseteq A \cup B} w_e$
        double in_c_S = 0.0;
        double in_target_c_S = 0.0;
        bool edgeBelongs_c=true;
        bool edgeBelongs_target_c=true;
        std::set<edgeid> border_edge;
        for (node n : S){
            for (edgeid eid: G.edgesOf(n)){
                border_edge.insert(eid);
            }
        }

        for (edgeid eid: border_edge){
            edgeBelongs_c =false;
            for (node nid: G.edgeMembers(eid)){
                if (!edgeBelongs_c && (S.find(nid) == S.end())){ // eid is not strictly included in S, so it's possible to take eid into account
                    edgeBelongs_c = true;
                }
                if (zeta[nid] !=c){
                    edgeBelongs_c = false;
                    break; 
                }
            }
            if (edgeBelongs_c){
                in_c_S +=G.getEdgeWeight(eid);
            }
                
            edgeBelongs_target_c=false;
            for (node nid: G.edgeMembers(eid)){
                if (!edgeBelongs_target_c && (S.find(nid) == S.end())){ // eid is not strictly included in S, so it's possible to take eid into account
                    edgeBelongs_target_c = true;
                }
                if (zeta[nid] !=target_c && (S.find(nid) == S.end())){
                    edgeBelongs_target_c = false;
                    break; 
                }
            }
            if (edgeBelongs_target_c){
                in_target_c_S +=G.getEdgeWeight(eid);
            }
        }
        covGain = (in_target_c_S - in_c_S) / totalEdgeWeight;

        //expCovGain = \sum_{d\geq 2} \frac{|E_d|_w}{|E|_w vol_w(V)^d}(vol_w(c'+S)^d -vol_w(c')^d -vol_w(c)^d + vol_w(c-S)^d)
        for (double j=0 ; j < d+1; j++){
            expCovGain += (EdgeSizeWeight[j] / (totalEdgeWeight * pow(volume,j))) * ( pow(vol_target_c + vol_S,j)-pow(vol_target_c,j)  - pow(vol_c,j) + pow(vol_c - vol_S,j));
        }
    }


    if (type_contribution==1) {// >>>>> MAJORITY EDGE CONTRIBUTION
        // covGain = \frac{In(S,c')-In(S,c)}{|E|} 
        // with In_w(A,B) = \sum_{e \in E ~:~ e \nsubseteq_{maj} A ~ \wedge ~ e \nsubseteq_{maj} B ~ \wedge ~ e \subseteq_{maj} A \cup B} w_e
        double in_c_S = 0.0;
        double in_target_c_S = 0.0;
        bool edgeBelongs_c=true;
        bool edgeBelongs_target_c=true;
        std::set<edgeid> border_edge;
        for (node n : S){
            for (edgeid eid: G.edgesOf(n)){
                border_edge.insert(eid);
            }
        }

        for (edgeid eid: border_edge){
            edgeBelongs_c =false;
            edgeBelongs_target_c=false;
            std::vector<double> CommEdge(3, 0.0);
            for (node nid: G.edgeMembers(eid)){
                if (zeta[nid]==c){CommEdge[0]++;}
                if (zeta[nid]==target_c){CommEdge[1]++;}
                if (S.find(nid) != S.end()){CommEdge[2]++;}
            }
                
            if (!(CommEdge[1] >= int(G.order(eid)/ 2. +1)) && !(CommEdge[0]- CommEdge[2] >= int(G.order(eid))/ 2. +1)) {
                if (CommEdge[0] >= int(G.order(eid)/ 2. +1)){
                    edgeBelongs_c =true;
                }
                if (CommEdge[1] + CommEdge[2] >= int(G.order(eid)/ 2. +1)){
                    edgeBelongs_target_c =true;
                }
            }
            
            if (edgeBelongs_c){
                in_c_S +=G.getEdgeWeight(eid);
            }

            if (edgeBelongs_target_c){
                in_target_c_S +=G.getEdgeWeight(eid);
            }
        }
        covGain = (in_target_c_S - in_c_S) / totalEdgeWeight;


        // expCovGain = \sum_{d\geq 2} \frac{|E_d|}{|E| vol(V)^d} \sum_{i=\frac{d}{2}+1}^{d} \binom{d}{i} (vol(c'+S)^i(vol(V)-vol(c'+S))^{d-i} -vol(c')^i(vol(V)-vol(c'))^{d-i} -vol(c)^i(vol(V)-vol(c))^{d-i} + vol(c-S)^i(vol(V)-vol(c-S))^{d-i})
        double sum = 0.0;
        for (double j=0 ; j < d+1; j++){
            sum = 0.0;
            for (double i= int ((j / 2.) + 1); i <= j; i++){
                sum += BinomialCoefficient(j,i) *(pow(vol_target_c + vol_S,i)*pow(volume - vol_target_c - vol_S,j-i) -  pow(vol_target_c,i)*pow(volume - vol_target_c,j-i) - pow(vol_c ,i)*pow(volume - vol_c,j-i) + pow(vol_c - vol_S,i)*pow(volume - vol_c + vol_S,j-i));
            }
            expCovGain += (EdgeSizeWeight[j] / (totalEdgeWeight * pow(volume,j))) * sum;
        }
    } 


    if (type_contribution==10) {// >>>>> STRICT EDGE CONTRIBUTION // PARTIALLY WEIGHTED
        //covGain = \frac{In_w(S,c')-In_w(S,c)}{|E|_w} 
        //with $In_w(A,B) = \sum_{e \in E ~:~ e \nsubseteq A ~ \wedge ~ e \nsubseteq B ~ \wedge ~ e \subseteq A \cup B} w_e$
        double in_c_S = 0.0;
        double in_target_c_S = 0.0;
        bool edgeBelongs_c=true;
        bool edgeBelongs_target_c=true;
        std::set<edgeid> border_edge;
        for (node n : S){
            for (edgeid eid: G.edgesOf(n)){
                border_edge.insert(eid);
            }
        }

        for (edgeid eid: border_edge){
            edgeBelongs_c =false;
            for (node nid: G.edgeMembers(eid)){
                if (!edgeBelongs_c && (S.find(nid) == S.end())){ // eid is not strictly included in S, so it's possible to take eid into account
                    edgeBelongs_c = true;
                }
                if (zeta[nid] !=c){
                    edgeBelongs_c = false;
                    break; 
                }
            }
            if (edgeBelongs_c){
                in_c_S +=G.getEdgeWeight(eid);
            }
                
            edgeBelongs_target_c=false;
            for (node nid: G.edgeMembers(eid)){
                if (!edgeBelongs_target_c && (S.find(nid) == S.end())){ // eid is not strictly included in S, so it's possible to take eid into account
                    edgeBelongs_target_c = true;
                }
                if (zeta[nid] !=target_c && (S.find(nid) == S.end())){
                    edgeBelongs_target_c = false;
                    break; 
                }
            }
            if (edgeBelongs_target_c){
                in_target_c_S +=G.getEdgeWeight(eid);
            }
        }
        covGain = (in_target_c_S - in_c_S) / totalEdgeWeight;

        //expCovGain = \sum_{d\geq 2} \frac{|E_d|_w}{|E|_w vol_w(V)^d}(vol_w(c'+S)^d -vol_w(c')^d -vol_w(c)^d + vol_w(c-S)^d)
        for (double j=0 ; j < d+1; j++){
            expCovGain += (EdgeSizeWeight[j] / (totalEdgeWeight * pow(volume_unweight,j))) * ( pow(vol_target_c_unweight + vol_S_unweight,j)-pow(vol_target_c_unweight,j)  - pow(vol_c_unweight,j) + pow(vol_c_unweight - vol_S_unweight,j));
        }
    }
    

    if (type_contribution==11) {// >>>>> MAJORITY EDGE CONTRIBUTION // PARTIALLY WEIGHTED
        // covGain = \frac{In(S,c')-In(S,c)}{|E|} 
        // with In_w(A,B) = \sum_{e \in E ~:~ e \nsubseteq_{maj} A ~ \wedge ~ e \nsubseteq_{maj} B ~ \wedge ~ e \subseteq_{maj} A \cup B} w_e
        double in_c_S = 0.0;
        double in_target_c_S = 0.0;
        bool edgeBelongs_c=true;
        bool edgeBelongs_target_c=true;
        std::set<edgeid> border_edge;
        for (node n : S){
            for (edgeid eid: G.edgesOf(n)){
                border_edge.insert(eid);
            }
        }

        for (edgeid eid: border_edge){
            edgeBelongs_c =false;
            edgeBelongs_target_c=false;
            std::vector<double> CommEdge(3, 0.0);
            for (node nid: G.edgeMembers(eid)){
                if (zeta[nid]==c){CommEdge[0]++;}
                if (zeta[nid]==target_c){CommEdge[1]++;}
                if (S.find(nid) != S.end()){CommEdge[2]++;}
            }
                
            if (!(CommEdge[1] >= int(G.order(eid)/ 2. +1)) && !(CommEdge[0]- CommEdge[2] >= int(G.order(eid)/ 2. +1)) ){
                if (CommEdge[0] >= int(G.order(eid)/ 2. +1)){
                    edgeBelongs_c =true;
                }
                if (CommEdge[1] + CommEdge[2] >= int(G.order(eid)/ 2. +1)){
                    edgeBelongs_target_c =true;
                }
            }
            
            if (edgeBelongs_c){
                in_c_S +=G.getEdgeWeight(eid);
            }

            if (edgeBelongs_target_c){
                in_target_c_S +=G.getEdgeWeight(eid);
            }
        }
        covGain = (in_target_c_S - in_c_S) / totalEdgeWeight;


        // expCovGain = \sum_{d\geq 2} \frac{|E_d|}{|E| vol(V)^d} \sum_{i=\frac{d}{2}+1}^{d} \binom{d}{i} (vol(c'+S)^i(vol(V)-vol(c'+S))^{d-i} -vol(c')^i(vol(V)-vol(c'))^{d-i} -vol(c)^i(vol(V)-vol(c))^{d-i} + vol(c-S)^i(vol(V)-vol(c-S))^{d-i})
        double sum = 0.0;
        for (double j=0 ; j < d+1; j++){
            sum = 0.0;
            for (double i= int ((j / 2.) + 1); i <= j; i++){
                sum += BinomialCoefficient(j,i) *(pow(vol_target_c_unweight + vol_S_unweight,i)*pow(volume_unweight - vol_target_c_unweight - vol_S_unweight,j-i) -  pow(vol_target_c_unweight,i)*pow(volume_unweight - vol_target_c_unweight,j-i) - pow(vol_c_unweight ,i)*pow(volume_unweight - vol_c_unweight,j-i) + pow(vol_c_unweight - vol_S_unweight,i)*pow(volume_unweight - vol_c_unweight + vol_S_unweight,j-i));
            }
            expCovGain += (EdgeSizeWeight[j] / (totalEdgeWeight * pow(volume_unweight,j))) * sum;
        }
    } 

    /*Aux::Log::setLogLevel("DEBUG");
    INFO("true mod comp ", covGain); 
    INFO("minus ", expCovGain);*/
    modularityGain= covGain - gamma * expCovGain; 

    return modularityGain;
}

} /* namespace NetworKit */
