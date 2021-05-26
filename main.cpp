#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <set>
#include <string>
#include <algorithm>
#include <stdlib.h>    /* for exit */

#include "networks/network_types.h"
#include "qualities/quality_types.h"
#include "optimisers/optimiser_types.h"

#include <getopt.h>

using namespace std; 

bool bipartite = false;	
bool leftproj = true;
bool val = false;
int seed = 12345;
set<string> projectors = {"weighted", "hyperbolic"};
string projector = "";
//set<string> qualities = {"modularity", "barber_mod", "projected_mod"};
//string quality = "modularity";
set<string> optimisers = {"sync_labelprop", "async_labelprop", "aggregate", "louvain", "bilouvain", "projected_louvain", "dual_projection"};
string optimiser = "louvain";

string stringify(set<string> &x){
	string s;
	for (auto const& i : x)
	{
		s += i;
		s += ',';
	}
	s.pop_back();	
	return s;
}

void usage(){
	cout << "Finds communities in uni and bi-partite networks" << endl;
	cout << "./main filename [options]" << endl;
	cout << "Options:" << endl;
	cout << "\t--bipartite -b: flag to indicate a bipartite network" << endl;
	cout << "\t--project -p: project a bipartite network using: {" << stringify(projectors) << "} projections. Automatically uses projected_mod and projected_louvain." << endl;
	cout << "\t--right -r: flag to project onto right node set" << endl;
	//cout << "\t--quality -q: evaluate community quality using: {" << stringify(qualities) << "} qualities" << endl;
	cout << "\t--optimisers -o: optimise node labels using: {" << stringify(optimisers) << "} optimisers" << endl;
	cout << "\t--val -v: flag to print the quality value" << endl;
	cout << "\t--seed -s: set the random seed" << endl;
	
	exit(1);
}

int main(int argc, char **argv){
	

   static struct option loptions[] =    {
        {"bipartite",0,0,'b'},		
        {"right",0,0,'r'},		
        {"project",1,0,'p'},		
        //{"quality",1,0,'q'},		
        {"optimiser",1,0,'o'},		
        {"val",0,0,'v'},		
        {"seed",1,0,'s'},		
        {0,0,0,0}
    };
    if(argc<2){
		usage();
	}
	
    int c;
	while ((c = getopt_long(argc, argv, "brvs:p:o:",loptions,NULL)) >= 0) 
    {
		switch (c)
		{
		case 'b': 
			bipartite = true;
			break;
		case 'r': 
			leftproj = false;
			break;		
		case 'v': 
			val = true;
			break;		
		case 's': 
			seed = atoi(optarg);
			break;							
		case 'p': 
			projector = (string)optarg;
			if( projectors.find(projector) == projectors.end() ){
				cout << "Allowed arguments for projection: " << stringify(projectors) << endl;
				exit(1);
			} 
			break;
		/*case 'q': 
			quality = (string)optarg;
			if( qualities.find(quality) == qualities.end() ){
				cout << "Allowed arguments for quality: " << stringify(qualities) << endl;
				exit(1);
			} 
			break;*/		
		case 'o': 
			optimiser = (string)optarg;
			if( optimisers.find(optimiser) == optimisers.end() ){
				cout << "Allowed arguments for optimiser: " << stringify(optimisers) << endl;
				exit(1);
			} 
			break;					
		case '?':
			usage();
			break;
        default: 
			if(optarg!=NULL) {cerr << "Unknown argument:"+(string)optarg << endl; exit(1);}
			else {cerr << "Unknown argument:"; exit(1);}
        }
    }
    char* filename = argv[optind];
    
    if(!bipartite && (optimiser == "bilouvain" || optimiser == "projected_louvain" || optimiser == "dual_projection") ){ cerr << "Can only use " << optimiser << " on bipartite networks" << endl; exit(1); }
    if(bipartite && (optimiser == "sync_labelprop" || optimiser == "async_labelprop") ){ cerr << "Can only use " << optimiser << " on unipartite networks" << endl; exit(1); }
    if(optimiser == "louvain" && optimiser == "barber_mod" ){ cerr << "Can't use louvain on bipartite networks. Try bilouvain." << endl; exit(1); }
    if(projector != "" && !bipartite){ cerr << "Can only project bipartite networks. Use -b flag if " << filename << " is bipartite." << endl; exit(1); }
    
    
    srand(seed);
    if(bipartite){
		BipartiteNetwork B(filename);	
		if(projector != ""){
			if(projector == "weighted"){ 
				if(optimiser == "dual_projection"){
					BarberModularity Q(B);	
					DualProjection< BarberModularity, WeightedProjector > O(Q);
					O.optimise();
					O.print_labels();	
					if(val){cout << "BarberModularity = " << O.val << endl; }	
				} else if(optimiser == "projected_louvain"){
					WeightedProjector P(B, leftproj, false);  //never includes loops
					WeightedNetwork W = P.project();
					ProjectedModularity<WeightedProjector> Q(P, W);	
					ProjectedLouvain< ProjectedModularity<WeightedProjector> > O(Q);
					O.optimise();
					O.print_labels();	
					if(val){cout << "ProjectedModularity = " << O.val << endl; }	
				}
			}
			else if(projector == "hyperbolic"){ 
				if(optimiser == "dual_projection"){
					BarberModularity Q(B);	
					DualProjection< BarberModularity, HyperbolicProjector > O(Q);
					O.optimise();
					O.print_labels();	
					if(val){cout << "BarberModularity = " << O.val << endl; }	
				} else if(optimiser == "projected_louvain"){
					HyperbolicProjector P(B, leftproj, false); 
					WeightedNetwork W = P.project();
					ProjectedModularity<HyperbolicProjector> Q(P, W);	
					ProjectedLouvain< ProjectedModularity<HyperbolicProjector> > O(Q);
					O.optimise();
					O.print_labels();	
					if(val){cout << "ProjectedModularity = " << O.val << endl; }	
				}
			}	
		} else {
			BarberModularity Q(B);	
			if(optimiser == "aggregate"){
					Aggregate< BarberModularity > O(Q, true); //always does maximum aggregation
					O.optimise();
					O.print_labels();	
					if(val){cout << "BarberModularity = " << O.val << endl; }
			}
			else if(optimiser == "bilouvain"){
					BiLouvain< BarberModularity > O(Q);
					O.optimise();
					O.print_labels();	
					if(val){cout << "BarberModularity = " << O.val << endl; }
			}
			
		}
		
	} else {
		WeightedNetwork W(filename);
		Modularity<WeightedNetwork> Q(W);	

		if(optimiser == "sync_labelprop"){
				SyncLabelProp< Modularity<WeightedNetwork> > O(Q);
				O.optimise();
				O.print_labels();	
				if(val){cout << "Modularity = " << O.val << endl; }
		}
		else if(optimiser == "async_labelprop"){
				AsyncLabelProp< Modularity<WeightedNetwork> > O(Q);
				O.optimise();
				O.print_labels();	
				if(val){cout << "Modularity = " << O.val << endl; }
		}
		else if(optimiser == "louvain"){
				Louvain< Modularity<WeightedNetwork> > O(Q);
				O.optimise();
				O.print_labels();	
				if(val){cout << "Modularity = " << O.val << endl; }
		}
		else if(optimiser == "aggregate"){
				Aggregate< Modularity<WeightedNetwork> > O(Q, true); //always does maximum aggregation
				O.optimise();
				O.print_labels();	
				if(val){cout << "Modularity = " << O.val << endl; }
		}
	
	}
	
	return 1;
}
