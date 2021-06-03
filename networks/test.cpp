#include <iostream>
#include "network.h"
#include "weighted_network.h"

using namespace std; 

int main(int argc, char **argv){

	if(argc==1){ cout << "Options are induced, join" << endl; }
	else if(string(argv[1]) == "join"){

			string filename = "../test_data/graph3.txt";
			
			cout << "\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\" << endl;
			cout << "\\CHECKING " << filename << endl;
			cout << "\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\" << endl;
			WeightedNetwork W(filename.c_str());
			cout << "BEFORE" << endl;
			W.print_basic();
			cout << "--------------------------------------------------------" << endl;
			W.join(0,2);
			W.join(4,5);
			W.join(1,3); 
			cout << "AFTER JOINED 2->0, 5->4, 3->1" << endl;
			W.print_basic();
			cout << "--------------------------------------------------------" << endl;			

			
			WeightedNetwork W2(filename.c_str());
			cout << "BEFORE" << endl;
			W2.print_basic();
			cout << "--------------------------------------------------------" << endl;
			W2.join(5,0);
			W2.join(2,3);
			cout << "AFTER JOINED 0->5, 3->2" << endl;
			W2.print_basic();
			cout << "--------------------------------------------------------" << endl;			

	} else if(string(argv[1]) == "induced"){
		vector< string > filenames = { 
			"../test_data/graph1.txt", 
			"../test_data/graph2.txt", 
			"../test_data/graph3.txt"};
			
		vector< vector< unordered_map<int,int> > > labels = {
			{//graph1	
				{ {0,0},{1,0},{2,0},{3,1},{4,1},{5,1} },
				{ {0,0},{1,1},{2,0},{3,1},{4,2},{5,2} }
			},
			{//graph2	
				{ {0,0},{1,0},{2,0},{3,1},{4,1},{5,1},{6,2},{7,2} },
				{ {0,0},{1,0},{2,0},{3,0},{4,0},{5,0},{6,1},{7,1} }
			},
			{//graph3	
				{ {0,0},{1,0},{2,0},{3,1},{4,1},{5,1} },
				{ {0,0},{1,1},{2,0},{3,1},{4,2},{5,2} }
			}
		};

		for(unsigned f=0; f<filenames.size(); ++f){
			string filename = filenames[f];
			
			cout << "\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\" << endl;
			cout << "\\CHECKING " << filename << endl;
			cout << "\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\" << endl;
			for(unsigned i=0; i<labels[f].size(); ++i){
				WeightedNetwork W(filename.c_str());
				cout << "BEFORE" << endl;
				W.print_basic();
				cout << "--------------------------------------------------------" << endl;							
				WeightedNetwork WI = W.induced_graph( labels[f][i] );
				cout << "AFTER ( "; for(auto &j : labels[f][i]){ cout << j.first << ":" << j.second << " "; } cout << ")" << endl;		
				WI.print_basic();
				cout << "--------------------------------------------------------" << endl;							

			}
		}

	} else {
		cout << "Options are induced, join" << endl;
	}


	return 0;
}
