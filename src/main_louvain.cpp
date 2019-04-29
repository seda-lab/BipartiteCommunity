// File: main_louvain.cpp
// -- community detection, sample main file
//-----------------------------------------------------------------------------
// Community detection 
// Based on the article "Fast unfolding of community hierarchies in large networks"
// Copyright (C) 2008 V. Blondel, J.-L. Guillaume, R. Lambiotte, E. Lefebvre
//
// And based on the article
// Copyright (C) 2013 R. Campigotto, P. Conde Céspedes, J.-L. Guillaume
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
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Lesser General Public License for more details.
// 
// You should have received a copy of the GNU Lesser General Public License
// along with Louvain algorithm.  If not, see <http://www.gnu.org/licenses/>.
//-----------------------------------------------------------------------------
// Author   : E. Lefebvre, adapted by J.-L. Guillaume and R. Campigotto
// Email    : jean-loup.guillaume@lip6.fr
// Location : Paris, France
// Time	    : July 2013
//-----------------------------------------------------------------------------
// see README.txt for more details


#include "graph.h"
#include "graph_binary.h"
#include "louvain.h"
#include <unistd.h>
#include <set>

#include "modularity.h"
#include "bimodularity.h"
#include "barbermodularity.h"

using namespace std;


vector< vector<int> > levels;
vector< vector<int> > node_maps;
char *filename = NULL;
char *rel = NULL;
char *filename_part = NULL;
char *filename_proj = NULL;

vector<int> first_renum(0);
vector<int> proj_renum(0);


int randseed = -1;
int nb_pass = 0;
long double precision = 0.000001L;
unsigned short id_qual = 0;

Quality *q;

vector<long double> q_degree;
vector<long double> d_degree;
long double F;
unsigned int bipartite = 0;
unsigned int bipartite_original = 0;
int project = -1;
bool weighted = false;
bool agg = false;
bool lou = true;


void bidegrees(Graph &g, unsigned int bip){
  
  q_degree.resize(0);
  d_degree.resize(0);
  F = 0;
  for (unsigned int i=0 ; i<g.links.size() ; i++) {
	  if(i<bip){
		  q_degree.push_back(0);
	  } else {
		  d_degree.push_back(0);
	  }
	  for (unsigned int j=0 ; j<g.links[i].size() ; j++) {
		    if(i<bip){
				q_degree[i] += g.links[i][j].second;
			} else {
				d_degree[i-bip] += g.links[i][j].second;
			}			
	   }
   }
   for(unsigned int i=0; i<q_degree.size(); ++i){ F += q_degree[i]; }

}


void
usage(char *prog_name, const char *more) {
  cerr << more;
  cerr << "usage: " << endl;
  cerr << "-i input_file: file containing the graph to decompose in communities" << endl;

  cerr << "-q id\tthe quality function used to compute partition of the graph (modularity is chosen by default):" << endl << endl;

  cerr << "\tid = 0\t -> the classical Newman-Girvan criterion (also called \"Modularity\")" << endl;
  cerr << "\tid = 10\t -> the Bipartite Modularity criterion" << endl;
  cerr << "\tid = 11\t -> the Barber Modularity criterion" << endl;
  cerr << endl;
  cerr << "-b int\tnumber of left nodes" << endl;
  cerr << "-a try agglomeration" << endl;
  cerr << "-A only agglomeration" << endl;
  cerr << "-x random seed" << endl;
  cerr << "-w weighted" << endl;
  cerr << "-p projection" << endl;
  cerr << "-P projection file" << endl;
  cerr << "-t starting partition" << endl;
  cerr << "-h show this usage message" << endl;

  exit(0);
}

void
parse_args(int argc, char **argv) {
	
  if (argc<2){ usage(argv[0], "Bad arguments number\n"); }

  for (int i = 1; i < argc; i++) {
    if(argv[i][0] == '-') {		
      switch(argv[i][1]) {
	case 'i':
		if (i==argc-1){ usage(argv[0], "Infile missing\n"); }
		filename = argv[i+1];
		i++;
	break;
	case 't':
		filename_part = argv[i+1];
		i++;
	break;
	case 'a':
        agg = true;
	break;
	case 'A':
		lou = false;
        agg = true;
	break;
	case 'r':
        rel = argv[i+1];
		i++;
	break;
	case 'x':
        randseed = atoi(argv[i+1]);
		i++;
	break;
    case 'q':
		id_qual = (unsigned short)atoi(argv[i+1]);
		i++;
	break;
	case 'w':
		weighted = true;
	break;
	case 'b':
		bipartite = atoi(argv[i+1]);
		bipartite_original = bipartite;
	i++;
	break;
	case 'p':
		project = atoi(argv[i+1]);
		i++;
	break;
	case 'P':
		filename_proj = argv[i+1];
		i++;
	break;
    case 'e':
		precision = atof(argv[i+1]);
		i++;
	break;
    case 'h':
		usage(argv[0], "");
	break;
    default:
		usage(argv[0], "Unknown option\n");
      }
    } 
  }
  if (filename == NULL){ usage(argv[0], "No input file has been provided\n"); }
  
}



void
init_quality(BinaryGraph *g, unsigned short nbc, int qual_type) {

  if (nbc > 0) delete q;

  switch (qual_type) {
  case 0:
    q = new Modularity(*g);
    break;
  case 10:
    if(project == 0){ q = new BiModularity(*g, q_degree); }
    else if(project == 1){ q = new BiModularity(*g, d_degree); }
    else{ cerr << "Only use this with a projection!" << endl; exit(1); }
    break;
  case 11:
    q = new BarberModularity(*g, q_degree, d_degree);
    break;
  default:
    q = new Modularity(*g);
    break;
  }
  
}


int
main(int argc, char **argv) {
  
  parse_args(argc, argv);

  if(randseed > 0){
	  srand(randseed); 
  }else{	
	  srand(time(NULL)+getpid()); 
  }
  
 
  //read the input text file
  Graph gin(filename, weighted, bipartite);
  gin.clean();
  //if we renumber and then project we will lose this map, save it for later.
  if( gin.renumber(rel) ){first_renum = gin.renum; }
  
  if( gin.bipartite > 0 ){ bidegrees(gin, gin.bipartite); bipartite_original = gin.bipartite; }
  if( project >= 0 ){ 
	  gin.project(project);
	  if( gin.links.size() <= 1 ){
		  cout << "No projection" << endl;
		  exit(1);
	  }

	  vector<long double> tmp( gin.links.size() );
	  if(project == 0 && proj_renum.size() > 0){
		  for(unsigned int i=0; i<tmp.size() ; i++){ tmp[i] = q_degree[  proj_renum[i]  ]; } 
		  q_degree = tmp;
	  } else if(project == 1 && proj_renum.size() > 0){
		  for(unsigned int i=0; i<tmp.size() ; i++){ tmp[i] = d_degree[  proj_renum[i]-bipartite_original  ]; }
		  d_degree = tmp;
	  }
	  if(filename_proj){
		  gin.display_txt(true, filename_proj);
	  }
  }
  vector<int> final_n2c(gin.links.size());
  for (unsigned int i=0 ; i<final_n2c.size() ; i++){ final_n2c[i] = i; } //initial community assignment
  

  //gin.display(true);
  unsigned short nb_calls = 0;
  if(lou){ //do the louvain method

	  //generate the binary graph that the louvain method wants
	  BinaryGraph g;
	  g.generate_from_links(gin.links, gin.bipartite);
  
      //map from node to communtiy
	  vector<int> n2c;
	  //map from node to new node id
	  vector<int> new_node;

      init_quality(&g, nb_calls, id_qual);
      nb_calls++;
  
      cerr << endl << "Computation of communities with the " << q->name << " quality function" << endl << endl;

	  bool improvement = true;
	  Louvain c(-1, precision, q);
	  if(id_qual == 11){ c.side = (q->Nleft > q->Nright)?0:1; } //start the bilouvain algorithm on the larger side
	  if(filename_part){ c.file_init_partition(filename_part, first_renum, proj_renum); } //inital partition if provided
	  long double quality = (c.qual)->quality();
	  cerr << "starting quality " << quality << endl;
	  long double new_qual;
	  
	  /*if(filename_part){ 
		  //(c.qual)->print_internal(); 
		  exit(1);}*/
	  
  do {

      cerr << "network size: " 
	   << (c.qual)->g.nb_nodes << " nodes, " 
	   << (c.qual)->g.nb_links << " links, "
	   << (c.qual)->g.total_weight << " weight" << endl;
    
    
		//louvain - try adding things into communities
		improvement = c.one_level();
		new_qual = (c.qual)->quality();

		cerr << "new qual " << new_qual << " from " << quality << endl;
		if(id_qual != 11){  levels.push_back( (c.qual)->n2c ); 	} //qual 11 needs care

		//compute the induced graph
		if(id_qual == 11){
			g = c.partition2bipartitegraph_binary(q_degree, d_degree, n2c, new_node);
			levels.push_back( n2c ); //save this community assignment
			node_maps.push_back( new_node ); //save the map from the start graph to the induced graph
		} else {
			if(project == 1){ g = c.partition2graph_binary(d_degree); }
			else{ g = c.partition2graph_binary(q_degree); }
		}

		nb_calls++;		
		init_quality(&g, nb_calls, id_qual); 
		int old_side = c.side;
		c = Louvain(-1, precision, q);
		if(id_qual == 11){ 
			c.side = (old_side+1)%2; 
			c.vector_init_partition(n2c);
		} 

		cerr << "  quality increased from " << quality << " to " << new_qual << endl;
		quality = new_qual;

	} while(improvement);

    cerr << (c.qual)->name << " quality " << (c.qual)->quality()<< endl;

    
      if(id_qual == 11){
		  for (unsigned int l=0 ; l<levels.size() ; l++){
			  for (unsigned int node=0 ; node<final_n2c.size() ; node++){
				final_n2c[node] = node_maps[l][final_n2c[node]];
			  }
		  } //map from original node ids to final induced graph
		  for (unsigned int node=0 ; node<final_n2c.size() ; node++){
				final_n2c[node] = levels[levels.size()-1][final_n2c[node]];
		  } //map from original node ids to communities of induced graph
	  } else {
		  for (unsigned int l=0 ; l<levels.size() ; l++){
			  for (unsigned int node=0 ; node<final_n2c.size() ; node++){
				final_n2c[node] = levels[l][final_n2c[node]];
			  }
		  }
	  }
	  
	  
  }
  
  
  
  if(agg){
	cerr << "agglomeration " << gin.links.size() << " " << bipartite << endl;
	
	//start from scratch for agg only runs
	Graph gagg(filename, weighted, bipartite);
    gagg.clean();
    if( gagg.renumber(rel) && first_renum.size() == 0){first_renum = gagg.renum; }
    if( bipartite > 0 ){ bidegrees(gagg, gagg.bipartite); }
	  if( project >= 0 ){ 
		  gagg.project(project);
		  vector<long double> tmp( gagg.links.size() );
		  if(project == 0 && proj_renum.size() > 0){
			  for(unsigned int i=0; i<tmp.size() ; i++){ tmp[i] = q_degree[  proj_renum[i]  ]; }
			  q_degree = tmp;
		  } else if(project == 1 && proj_renum.size() > 0){
			  for(unsigned int i=0; i<tmp.size() ; i++){ tmp[i] = d_degree[  proj_renum[i]-bipartite_original  ]; }
			  d_degree = tmp;
		  }
	  }
  
    BinaryGraph ga;
    ga.generate_from_links(gagg.links, gagg.bipartite);

    init_quality(&ga, 0, id_qual);
    Louvain c(-1, precision, q);    
    if( !lou && filename_part ){ c.file_init_partition(filename_part, first_renum, proj_renum); }  //read in initial community assignment
    else{ c.vector_init_partition(final_n2c); } //use community assignment computed above

	cerr << "agg starting quality " << (c.qual)->quality() << endl;
	
	c.one_level_agg();
	(c.qual)->renumber_comms();

	for (unsigned int node=0 ; node<final_n2c.size() ; node++){ final_n2c[node] = (c.qual)->n2c[node]; }

	cerr << "agg final quality " << (c.qual)->quality() << endl;
  }
  
  //some extra checking
  //start from scratch again!
  cerr << "CHECKING" << endl;

  Graph gcheck(filename, weighted, bipartite);
  gcheck.clean();
  if( gcheck.renumber(rel) && first_renum.size() == 0){ first_renum = gcheck.renum; }
  if( bipartite > 0 ){ bidegrees(gcheck, gcheck.bipartite); bipartite_original = gcheck.bipartite; }
  if( project >= 0 ){ 
	  gcheck.project(project);
	  vector<long double> tmp( gcheck.links.size() );
	  if(project == 0 && proj_renum.size() > 0){
		  for(unsigned int i=0; i<tmp.size() ; i++){ tmp[i] = q_degree[  proj_renum[i]  ]; }
		  q_degree = tmp;
	  } else if(project == 1 && proj_renum.size() > 0){
		  for(unsigned int i=0; i<tmp.size() ; i++){ tmp[i] = d_degree[  proj_renum[i]-bipartite_original  ]; }
		  d_degree = tmp;
	  }
  }

  BinaryGraph gbcheck;
  gbcheck.generate_from_links(gcheck.links, gcheck.bipartite);

	    
  init_quality(&gbcheck, 0, id_qual);  
  Louvain c_check(-1, precision, q);  

  c_check.vector_init_partition(final_n2c);
  cerr << "quality " << (c_check.qual)->quality() << " " << (c_check.qual)->name << endl;

  if( proj_renum.size() > 0 || first_renum.size() > 0 ){ //return result with original labels

    vector<int> renum;
    if(proj_renum.size() > 0 && first_renum.size() > 0){
		renum.resize(proj_renum.size());		
		for(unsigned int i=0; i<proj_renum.size(); ++i){ renum[i] = first_renum[ proj_renum[ i ] ]; }
	} else if( proj_renum.size() > 0 ) {
		renum = proj_renum;
	} else {
		renum = first_renum;
	}

    vector<int>::iterator it = max_element( renum.begin(),  renum.end());
	int max_idx = (*it);
	
	vector<int> all_comms( (max_idx + 1) + gin.missing.size(), -1); //all_comms[ old_id ] = community or -1
	set<int> comms; //communities
	for (unsigned int node=0 ; node<final_n2c.size() ; node++){ 
		comms.insert(final_n2c[node]);
		all_comms[ renum[node] ] = final_n2c[node];
	} 
	if( gin.missing.size() > 0 ){
		int nc = comms.size(); //number of communities
		//if projection has removed some of the nodes, assign these to their own community
		for (unsigned int node=0 ; node<gin.missing.size() ; node++){ 
			if( first_renum.size() > 0 ){
				all_comms[ first_renum[ gin.missing[node] ] ] = nc++; 
			} else {
				all_comms[ gin.missing[node] ] = nc++; 
			}
		} 
	}
	for (unsigned int node=0 ; node<all_comms.size() ; node++){ 
		if(all_comms[node] > -1){cout << node << " " << all_comms[node] << endl; }
	}
	
  } else{
	for (unsigned int node=0 ; node<final_n2c.size() ; node++){ cout << node << " " << final_n2c[node] << endl; }
  } 
  
  delete q;
}


