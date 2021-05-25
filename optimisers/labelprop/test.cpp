#include <iostream>
#include "../optimiser.h"
#include "labelprop.h"
#include "sync_labelprop.h"
#include "async_labelprop.h"


using namespace std; 

int main(int argc, char **argv){
	
	if(argc==1){ cout << "Options are sync, async" << endl; }
	else if(string(argv[1]) == "sync"){
		vector< string > filenames = { "../../test_data/graph1.txt", "../../test_data/graph2.txt", "../../test_data/graph3.txt"};

		for(auto &filename : filenames){
			cout << "\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\" << endl;
			cout << "\\CHECKING " << filename << endl;
			cout << "\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\" << endl;
			WeightedNetwork W(filename.c_str());
			W.print_basic();
			Modularity<WeightedNetwork> Q(W);			
			Q.print_basic();
			
			SyncLabelProp< Modularity<WeightedNetwork> > L( Q , true );
			L.optimise();
			L.print_basic();
			
		}
	} else if(string(argv[1]) == "async"){
		vector< string > filenames = { "../../test_data/graph1.txt", "../../test_data/graph2.txt", "../../test_data/graph3.txt"};

		for(auto &filename : filenames){
			cout << "\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\" << endl;
			cout << "\\CHECKING " << filename << endl;
			cout << "\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\" << endl;
			WeightedNetwork W(filename.c_str());
			W.print_basic();
			Modularity<WeightedNetwork> Q(W);			
			Q.print_basic();
			
			AsyncLabelProp< Modularity<WeightedNetwork> > L( Q , true );
			L.optimise();
			L.print_basic();
			
		}
	} else {
		cout << "Options are sync, async" << endl;
	}


	return 0;
}
