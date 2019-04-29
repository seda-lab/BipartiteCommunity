// File: louvain.cpp
// -- community detection source file
//-----------------------------------------------------------------------------
// Community detection
// Based on the article "Fast unfolding of community hierarchies in large networks"
// Copyright (C) 2008 V. Blondel, J.-L. Guillaume, R. Lambiotte, E. Lefebvre
//
// This file is part of Louvain algorithm.
// 
// Louvain algorithm is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// Louvain algorithm is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public License
// along with Louvain algorithm.  If not, see <http://www.gnu.org/licenses/>.
//-----------------------------------------------------------------------------
// Author   : E. Lefebvre, adapted by J.-L. Guillaume
// Email    : jean-loup.guillaume@lip6.fr
// Location : Paris, France
// Time	    : February 2008
//-----------------------------------------------------------------------------
// see readme.txt for more details

#include "louvain.h"

using namespace std;


Louvain::Louvain(int nbp, long double epsq, Quality* q) {
  qual = q;
  neigh_weight.resize(qual->size,-1);
  neigh_pos.resize(qual->size);

  neigh_last = 0;

  nb_pass = nbp;
  eps_impr = epsq;
  side = -1; 
}



void
Louvain::vector_init_partition(vector<int> &new_n2c) {

	for(unsigned int node =0; node<new_n2c.size(); node++){
		
		int comm = new_n2c[node];
		int old_comm = qual->n2c[node];

		neigh_comm(node);
		qual->remove(node, old_comm, neigh_weight[old_comm]); //pop node from its old community
      
		int i=0;
		for (i=0 ; i<neigh_last ; i++) {
			int best_comm = neigh_pos[i];
			long double best_nblinks = neigh_weight[neigh_pos[i]];
			if (best_comm==comm) {
				qual->insert(node, best_comm, best_nblinks);
				break;
			}
		}
		if (i==neigh_last){qual->insert(node, comm, 0);} //no links from node to comm
    }
	
}



void
Louvain::file_init_partition(char *filename, vector<int> &first_renum, vector<int> &proj_renum) {

  vector<int> fr, pr;
  vector<int> renum; //renum[old_node] = new_node;
  //if we have relabelled the original nodes
  if( first_renum.size() > 0 ){
	  fr.assign( (*max_element(first_renum.begin(), first_renum.end())) + 1 , -1);
	  for(unsigned int i=0; i<first_renum.size(); ++i){ fr[ first_renum[i] ] = i; }
	  renum = fr;
  }
  //if we relabelled after a projection
  if( proj_renum.size() > 0 ){
	  pr.assign( (*max_element(proj_renum.begin(), proj_renum.end())) + 1 , -1);
	  for(unsigned int i=0; i<proj_renum.size(); ++i){ pr[ proj_renum[i] ] = i; }
	  if( first_renum.size() > 0 ){
		  renum.assign( pr.size(), -1 );
		  for(unsigned int i=0; i<pr.size(); ++i){ renum[ i ] = pr[ fr[i] ]; }
	  } else {
		  renum = pr;
	  }
  }

  ifstream finput;
  finput.open(filename,fstream::in);
  if (finput.is_open() != true) {
    cerr << "The file " << filename << " does not exist" << endl;
    exit(EXIT_FAILURE);
  }
  vector<int> n2c, t_n2c;
  unsigned int ct = 0;
  while (!finput.eof()) {
    unsigned int node, comm;
    finput >> node >> comm;
	
    if (finput) {
		//if(node != ct){cerr << "input the partition in order" << endl; exit(1); }
		t_n2c.push_back(comm);
		ct+=1;
    }
  }
  finput.close();

  n2c = t_n2c;  
  vector_init_partition(n2c);
	
}

// compute the set of neighboring communities of node
// for each community, gives the number of links from node to comm
void
Louvain::neigh_comm(int node) {

  for (unsigned int i=0 ; i<neigh_weight.size() ; i++){ neigh_weight[neigh_pos[i]]=-1; } //set vector to -1s
  
  neigh_last = 0;
  
  pair<vector<int>::iterator, vector<long double>::iterator> p = (qual->g).neighbors(node); //neighbours of target node
  int deg = (qual->g).nb_neighbors(node); //number of neighbours of target node
  
  neigh_pos[0] = qual->n2c[node]; //community of target node
  neigh_weight[neigh_pos[0]] = 0; 
  neigh_last = 1;

  for (int i=0 ; i<deg ; i++) {
    int neigh  = *(p.first+i); //neighbour node
    int neigh_comm = qual->n2c[neigh]; //neighbour's community
    long double neigh_w = *(p.second+i);

    if (neigh!=node) { //don't count self loops
      if (neigh_weight[neigh_comm]==-1) { //first link to community 'neigh_comm'
		neigh_weight[neigh_comm] = 0.0L;
		neigh_pos[neigh_last++] = neigh_comm;
      }
      neigh_weight[neigh_comm] += neigh_w;
    }
  }
  //neigh last = number of communities target node is attached to
  //neigh_pos = communities node is attached to
  //neigh_weight = weight of links between target node and community: neigh_weight[node_comm] = \sum_{j in C} A_nj
}

void
Louvain::partition2graph() {
  vector<int> renumber(qual->size, -1);
  for (int node=0 ; node<qual->size ; node++) {
    renumber[qual->n2c[node]]++;
  }

  int end=0;
  for (int i=0 ; i< qual->size ; i++)
    if (renumber[i]!=-1)
      renumber[i]=end++;

  for (int i=0 ; i< qual->size ; i++) {
    pair<vector<int>::iterator, vector<long double>::iterator> p = (qual->g).neighbors(i);

    int deg = (qual->g).nb_neighbors(i);
    for (int j=0 ; j<deg ; j++) {
      int neigh = *(p.first+j);
      cout << renumber[qual->n2c[i]] << " " << renumber[qual->n2c[neigh]] << endl;
    }
  }
}

void
Louvain::display_partition() {
  vector<int> renumber(qual->size, -1);

  for (int node=0 ; node < qual->size ; node++) {
    renumber[qual->n2c[node]]++;
  }

  int end=0;
  for (int i=0 ; i < qual->size ; i++)
    if (renumber[i]!=-1)
      renumber[i] = end++;

  for (int i=0 ; i < qual->size ; i++)
    cout << i << " " << renumber[qual->n2c[i]] << endl;
}

BinaryGraph
Louvain::partition2graph_binary(vector<long double> &qin) {
  // Renumber communities
  int num_comms = qual->renumber_comms();

  // Compute communities
  vector<vector<int> > comm_nodes(num_comms); //map from node to community
  vector<int> comm_weight(num_comms, 0);	  //weight within each community??
  vector<long double> qin_comm(num_comms, 0); //bidegrees of community
    
  for (int node = 0 ; node < (qual->size) ; node++) {
    comm_nodes[ qual->n2c[node] ].push_back(node);
    comm_weight[ qual->n2c[node] ] += (qual->g).nodes_w[node];
    if( qin.size() > 0){
		qin_comm[ qual->n2c[node] ] += qin[node];
	}
  }
  qin = qin_comm;

  // Compute weighted graph
  BinaryGraph g2;
  
  //resize everything
  int nbc = comm_nodes.size();
  g2.nb_nodes = comm_nodes.size();
  g2.degrees.resize(nbc);
  g2.nodes_w.resize(nbc);
  
  //for every community
  for (int comm=0 ; comm<nbc ; comm++) {
    map<int,long double> m;
    map<int,long double>::iterator it;

    int size_c = comm_nodes[comm].size();
    g2.assign_weight(comm, comm_weight[comm]);
    for (int node=0 ; node<size_c ; node++) { //for every node in this communtiy

      pair<vector<int>::iterator, vector<long double>::iterator> p = (qual->g).neighbors(comm_nodes[comm][node]);//find neighbours
      int deg = (qual->g).nb_neighbors(comm_nodes[comm][node]);
      for (int i=0 ; i<deg ; i++) {
		int neigh = *(p.first+i);
		int neigh_comm = qual->n2c[neigh]; //neighbours' community
		long double neigh_weight = *(p.second+i);

		it = m.find(neigh_comm);
		if (it==m.end()){ //link this community with neighbours' community
			m.insert(make_pair(neigh_comm, neigh_weight));
		} else {
			it->second += neigh_weight;
		}
      }
    }

	//update induced graph
    g2.degrees[comm] = (comm==0)?m.size():g2.degrees[comm-1]+m.size();
    g2.nb_links += m.size();
    for (it = m.begin() ; it!=m.end() ; it++) {
      g2.total_weight += it->second;
      g2.links.push_back(it->first);
      g2.weights.push_back(it->second);
    }
    
  }

  return g2;
}


BinaryGraph
Louvain::partition2bipartitegraph_binary(vector<long double> &qin, vector<long double> &din, vector<int> &n2c, vector<int> &node_map) {
  //compute induced graph on one half of bipartite graph
  BinaryGraph g2;
  
  // Renumber communities
  int num_side_comms = qual->renumber_comms();
  
  // Compute communities
  int num_new_nodes = num_side_comms + ((side==0)?qual->Nright:qual->Nleft); //number of nodes in induced graph
  vector<vector<int> > comm_nodes(num_new_nodes);
  vector<int> comm_weight(num_new_nodes, 0);
  vector<long double> qin_comm(  ((side==0)?num_side_comms:qual->Nleft), 0);
  vector<long double> din_comm(  ((side==1)?num_side_comms:qual->Nright), 0);
  node_map.resize(qual->size); //node_map[ old_node_id ] = new_node_id
  
  for (int node = 0 ; node < (qual->size) ; node++) {
	int new_node = -1;
	if( side == 0 ){ //merge left nodes
		if(node < qual->Nleft){
			new_node = qual->n2c[node]; //label is community label
			qin_comm[new_node] += qin[node];
		} else {
			new_node = (node-qual->Nleft)+num_side_comms; //label is idx on right + num comms
			din_comm[(node-qual->Nleft)] += din[(node-qual->Nleft)];
		} 
	} else { //merge right nodes
		if(node < qual->Nleft){
			new_node = node; //label is idx
			qin_comm[new_node] += qin[node];
		} else {
			new_node = qual->Nleft + qual->n2c[node]; //Nleft + comm label
			din_comm[ qual->n2c[node] ] += din[(node-qual->Nleft)];
		} 
	}
	node_map[node] = new_node;
    comm_nodes[new_node].push_back(node);
    comm_weight[new_node] += (qual->g).nodes_w[node];
  }
  qin = qin_comm;
  din = din_comm;

  g2.nb_nodes = num_new_nodes; 
  g2.degrees.resize(num_new_nodes);
  g2.nodes_w.resize(num_new_nodes);
  g2.bipartite = qin.size();

  //for every new node
  for (int comm=0 ; comm<num_new_nodes ; comm++) {
    map<int,long double> m;
    map<int,long double>::iterator it;

	int size_c = comm_nodes[comm].size();
    g2.assign_weight(comm, comm_weight[comm]);
    //for every old node with this new id 
    for (int node=0 ; node<size_c ; node++) {
	
	 	
      pair<vector<int>::iterator, vector<long double>::iterator> p = (qual->g).neighbors(comm_nodes[comm][node]); //find neighbours
      int deg = (qual->g).nb_neighbors(comm_nodes[comm][node]);
      for (int i=0 ; i<deg ; i++) {
		int neigh = *(p.first+i);
		int neigh_comm = node_map[neigh]; //new id of neighbour
		long double neigh_weight = *(p.second+i);
		it = m.find(neigh_comm);
			if (it==m.end()){
			  m.insert(make_pair(neigh_comm, neigh_weight));  //link this node with neighbours' new_node
			} else {
			  it->second += neigh_weight;
			}
		}
    }
    if( comm == 0 || comm == g2.bipartite){
		g2.degrees[comm] = m.size();
	} else {
		g2.degrees[comm] = g2.degrees[comm-1]+m.size();
	}
    g2.nb_links += m.size();

    for (it = m.begin() ; it!=m.end() ; it++) {
      g2.total_weight += it->second;
      g2.links.push_back(it->first);
      g2.weights.push_back(it->second);
    }
  }
  g2.total_weight /= 2;
 
  //assign partition
  n2c.resize(num_new_nodes);
  for (int node=0 ; node<qual->size; node++) {
      n2c[node_map[node]] = qual->n2c[node];
  }
  return g2;
}


void 
Louvain::shuffle (vector<int> &x)
{
  for (unsigned int i=0; i<x.size(); ++i) {
    int rand_pos = rand()%x.size();
    //swap
    int tmp = x[i];
    x[i] = x[rand_pos];
    x[rand_pos] = tmp;
  }
}

bool
Louvain::one_level() {
	
  bool improvement=false ;
  int nb_moves;
  int nb_pass_done = 0;
  long double new_qual = qual->quality();
  long double cur_qual = new_qual;

  vector<int> random_order(qual->size);
  // repeat while 
  //   there is an improvement of quality
  //   or there is an improvement of quality greater than a given epsilon 
  //   or a predefined number of pass have been done
  do {

	  //random order to traverse nodes
      for (int i=0 ; i < qual->size ; i++){random_order[i] = i; }
      shuffle(random_order);
  
    cur_qual = new_qual;
    nb_moves = 0;
    nb_pass_done++;

    // for each node: remove the node from its community and insert it in the best community
    // if bipartite algorithm, only do nodes on one side!
    for (int node_tmp = 0 ; node_tmp < qual->size ; node_tmp++) {
		
      int node = random_order[node_tmp];
      int node_comm = qual->n2c[node]; //community of node
      long double w_degree = (qual->g).weighted_degree(node); //degree of node

      // computation of all neighboring communities of current node
      neigh_comm(node);
      
	  // take a node out of its community
      qual->remove(node, node_comm, neigh_weight[node_comm]);

      // compute the nearest community for node
      // default choice for future insertion is the former community
      int best_comm = node_comm;
      long double best_nblinks  = 0.0L;
      long double best_increase = 0.0L;
      
      //visit neighbours of node in random order
      vector<int> idx_order(neigh_last);
      for (int nn=0 ; nn < neigh_last ; ++nn ){ idx_order[nn] = neigh_pos[nn]; }
      shuffle(idx_order);
  
      for (int i=0 ; i<neigh_last ; i++) {
		  int ii = idx_order[i];
			long double increase = qual->gain(node, ii, neigh_weight[ii], w_degree); //gain from adding to community of neighbour
			if (increase>best_increase) {
			  best_comm = ii;
			  best_nblinks = neigh_weight[ii];
			  best_increase = increase;
			}
      }

      // insert node in the community yielding highest gain
      qual->insert(node, best_comm, best_nblinks);

      if (best_comm!=node_comm){nb_moves++;}
    }

    new_qual = qual->quality();

    if (nb_moves>0){improvement=true;}

  } while (nb_moves>0 && new_qual-cur_qual > eps_impr);

  qual->renumber_comms();
  return improvement;
}

bool
Louvain::one_level_agg() {
	  
  //communities are renumbered from zero and community to node lists are created;
  int ncomms = qual->fill_c2n();
  
  map<int, map<int, double> > R; //R[c][d] = connection strength from community c to community d
  for(int comm1=0; comm1<ncomms; comm1++){
	  for(unsigned int n=0; n<qual->c2n[comm1].size(); ++n){
		  int node = qual->c2n[comm1][n];
		  //neigh last = number of communities target node is attached to
		  //neigh_pos = communities node is attached to
		  //neigh_weight = weight of links between target node and community: neigh_weight[node_comm] = \sum_{j in C} A_nj
		  neigh_comm(node);
		  for(int i=0; i<neigh_last; ++i){
			  if( neigh_pos[i] != comm1 ){
				  if( R[comm1].find(neigh_pos[i]) == R[comm1].end() ){
					R[comm1][neigh_pos[i]] = neigh_weight[ neigh_pos[i] ];
				  } else {
					R[comm1][neigh_pos[i]] += neigh_weight[ neigh_pos[i] ];
				  }
			  }
	      }
	  }
  } 

  map<int, map<int, double> > S; //gains from joining communities
  double max_gain = 0;
  pair<int, int> join_pair; //pair with the most to gain
  for(int comm1=0; comm1<ncomms; comm1++){
	for(map<int, double>::iterator it=R[comm1].begin(); it != R[comm1].end(); ++it){
		int comm2 = it->first;
		double da = qual->agg_gain(comm1, comm2, it->second);
		S[comm1][comm2] = da;
		if( da > max_gain ){
			max_gain = da; 
			join_pair = make_pair<int, int>(comm1, comm2);
		}
	}  
  } 
  
  
  int joined = 0;
  
  while(max_gain > 0){
	  
	  if( joined % 1000 == 0 ){ cerr << "agg:: it " << joined << " gain " << max_gain << endl; }
	  
	  int comm1 = join_pair.first;
	  int comm2 = join_pair.second;
	  
	  qual->agg(comm1, comm2, R[comm1][comm2] );//join comm1 and comm2 (updates sums)
	  //NB we get rid of comm1 and keep comm2!
	  
	  //join rows of R
	  vector<int> col_ids;
	  for(map<int, double>::iterator it=R[comm1].begin(); it != R[comm1].end(); ++it){
		  int c = it->first;
		  if( c != comm2){ //diagonal elements are in qual
			if( R[comm2].find(c) == R[comm2].end() ){
				R[comm2][c] = it->second;
			} else {
				R[comm2][c] += it->second;
			}
			col_ids.push_back(c);
		  }
	  }
	  //join cols of R
	  for(unsigned int ci=0; ci<col_ids.size(); ++ci){ //works because of symmetrical links!
		  int c = col_ids[ci];
			if( R[c].find(comm2) == R[c].end() ){
				R[c][comm2] = R[c][comm1];
			} else {
				R[c][comm2] += R[c][comm1];
			}
	  }
	  //erase row comm1
	  R.erase(comm1);
	  //erase col comm1
	  R[comm2].erase(comm1);
	  for(unsigned int ci=0; ci < col_ids.size(); ++ci){  R[ col_ids[ci] ].erase( comm1 ); }
	  ncomms -= 1;
	  
	  //compute new gains
	  for(map<int, double>::iterator it=R[comm2].begin(); it != R[comm2].end(); ++it){
		int c = it->first;
		double da = qual->agg_gain(comm2, c, it->second);
		S[comm2][c] = da;
		S[c][comm2] = da;
	  }  
	  //find max gain
	  max_gain = 0;
	  for(map< int, map<int, double> >::iterator it1=R.begin(); it1!=R.end(); ++it1){
		int comm1 = it1->first;
		for(map<int, double>::iterator it=R[comm1].begin(); it != R[comm1].end(); ++it){
			int comm2 = it->first;
			if( comm1 < comm2 ){ //symmetry
				double da = S[comm1][comm2];
				if( da > max_gain ){
					max_gain = da; 
					join_pair = make_pair<int, int>(comm1, comm2);
				}
			}
		}  
	  }
	  joined++;

  }
  
  cerr << "agg joined " << joined << endl;
  return false;
}


