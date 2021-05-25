#include <iostream>
#include "../optimiser.h"
#include "louvain.h"
#include "bilouvain.h"
#include "projected_louvain.h"
#include "aggregate.h"


using namespace std; 

int main(int argc, char **argv){
	
	if(argc==1){ cout << "Options are louvain, aggregate, bilouvain, projected" << endl; }
	else if(string(argv[1]) == "louvain"){
		vector< string > filenames = { "../../test_data/graph1.txt", "../../test_data/graph2.txt", "../../test_data/graph3.txt"};

		for(auto &filename : filenames){
			cout << "\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\" << endl;
			cout << "\\CHECKING " << filename << endl;
			cout << "\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\" << endl;
			WeightedNetwork W(filename.c_str());
			W.print_basic();
			Modularity<WeightedNetwork> Q(W);			
			Q.print_basic();
			
			Louvain< Modularity<WeightedNetwork> > L(Q);
			L.print_basic();
			L.optimise();
			L.print_basic();
			
		}
	} else if(string(argv[1]) == "aggregate"){
		vector< string > filenames = { "../../test_data/graph1.txt", "../../test_data/graph2.txt", "../../test_data/graph3.txt"};

		for(auto &filename : filenames){
			cout << "\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\" << endl;
			cout << "\\CHECKING " << filename << endl;
			cout << "\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\" << endl;
			WeightedNetwork W(filename.c_str());
			W.print_basic();
			Modularity<WeightedNetwork> Q(W);			
			Q.print_basic();
			
			Aggregate< Modularity<WeightedNetwork> > L(Q, true);
			L.print_basic();
			L.optimise();
			L.print_basic();
			
		}
	} else if(string(argv[1]) == "bilouvain"){
		vector< string > filenames = { "../../test_data/bipartite1.txt"};

		for(auto &filename : filenames){
			cout << "\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\" << endl;
			cout << "\\CHECKING " << filename << endl;
			cout << "\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\" << endl;
			BipartiteNetwork W(filename.c_str());
			W.print_basic();
			BarberModularity Q(W);			
			Q.print_basic();
			
			BiLouvain< BarberModularity > L(Q);
			L.print_basic();
			L.optimise();
			L.print_basic();
			
		}
	} else if(string(argv[1]) == "projected"){
		vector< string > filenames = { "../../test_data/bipartite2.txt"};

		for(auto &filename : filenames){
			cout << "\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\" << endl;
			cout << "\\CHECKING " << filename << endl;
			cout << "\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\" << endl;
			for(unsigned loops=0; loops<2; ++loops){
				cout<< "LOOPS === " << loops << endl;
				BipartiteNetwork B(filename.c_str());
				B.print_basic();
				HyperbolicProjector P(B, true, bool(loops) ); //test hyperbolic on bipartite 2
				WeightedNetwork W = P.project();
				W.print_basic();
				ProjectedModularity<HyperbolicProjector> Q(P, W);
				Q.print_basic();
				
				ProjectedLouvain< ProjectedModularity<HyperbolicProjector> > L(Q);
				L.print_basic();
				L.optimise();
				L.print_basic();
			}
			
		}
		
	} else {
		cout << "Options are louvain, bilouvain, projected" << endl;
	}


	return 0;
}
