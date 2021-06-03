#include <iostream>
#include "quality.h"
/*#include "modularity.h"
#include "barber_modularity.h"
#include "projected_modularity.h"
#include "cp_modularity.h"
#include "inv_modularity.h"*/
#include "block_modularity.h"
#include <unordered_map>


using namespace std; 

int main(int argc, char **argv){
	
	if(argc==1){ cout << "Options are block" << endl; }
	else if(string(argv[1]) == "block"){
		vector< string > filenames = { 
			"../test_data/graph1.txt",
			/*"../test_data/graph4.txt",
			"../test_data/graph5.txt",
			"../test_data/graph6.txt",
			"../test_data/graph7.txt",
			"../test_data/graph8.txt",
			"../test_data/graph9.txt",
			"../test_data/graph10.txt"*/
		};

		vector< vector< unordered_map<int,int> > > labels = {
			 //graph1
			{ 
				{ {0,0},{1,0},{2,0},{3,1},{4,2},{5,2} }//,
				//{ {0,0},{1,0},{2,0},{3,1},{4,1},{5,1} }
				//{0,1,2,3,4,5}, {0,0,0,0,0,0}, 
				//{0,0,0,1,1,1}, {0,1,0,2,3,3},
		   		//{0,1,0,1,2,2}, {0,1,0,2,2,2}
			},
			/*//graph4
			{ 
				//0 1 2 3 4 5 6 7 8 9  10  11		             
				 {0,0,0,0,0,0,1,1,1,1, 1,  1}, //community
				 {0,0,1,1,1,1,2,3,3,2, 3,  3}, //CP
				 {0,1,2,3,4,5,6,7,8,9,10,11},
			},
			//graph5
			{
				//0 1 2 3 4 5 6 7 8 9  10  11  12  13  14  15
				 {0,0,0,0,0,0,1,1,1,1, 1,  2  ,2  ,2  ,2  ,2}, //community
				 {1,0,0,0,0,0,1,0,0,0, 0,  1  ,0  ,0  ,0  ,0}, //bipartite
				 {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15},
			},
			//graph6
			{ 
				{0,0,0,0, 1,1,1,1,  2,2,2,2},
				{0,0,0,0, 0,0,0,0,  1,1,1,1},
				{0,1,2,3,4,5,6,7,8,9,10,11},
				{0,0,0,0,0,0,0,0,0,0,0, 0}
			},
			//graph7
			{ 
				{0,0,0,0, 1,1,1,1,  2,2,2,2},
				{0,0,0,0, 0,0,0,0,  1,1,1,1},
				{0,1,2,3,4,5,6,7,8,9,10,11},
				{0,0,0,0,0,0,0,0,0,0,0, 0}
			},
			//graph8
			{ 
				
				{0,0,0,0, 1,1,1,1,1,1},
				{0,0,0,0, 1,1,1, 2,2,2},
				{0,1,2,3,4,5,6,7,8,9},
				{0,0,0,0,0,0,0,0,0,0}
			},
			//graph9
			{ 
				
				{0,0,0,0, 1,1,1,1, 2,2,2},
				{0,0,0,0, 1,1,1,1, 1,1,1},
				{0,0,0,0, 0,0,0,0, 1,1,1},
				{0,1,2,3,4,5,6,7,8,9,10},
				{0,0,0,0,0,0,0,0,0,0,0}
			},
			//graph9
			{ 
				
				{0,0,0,0,1,1,1,1},
				{0,0,1,1,2,2,2,2},
				{0,0,1,1,2,2,3,3},
				{0,1,2,3,4,5,6,7},
				{0,0,0,0,0,0,0,0}
			}*/
		};

		for(unsigned f=0; f<filenames.size(); ++f){
			string filename = filenames[f];
			
			cout << "\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\" << endl;
			cout << "\\CHECKING " << filename << endl;
			cout << "\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\" << endl;
			for(unsigned i=0; i<labels[f].size(); ++i){
				WeightedNetwork W(filename.c_str());
				BlockModularity<WeightedNetwork> Q(W,labels[f][i]);		
				Q.eval();
				Q.print_basic();
				
				cout << "x = ( "; for(auto &k : labels[f][i]){ cout << k.first << ":" << k.second << " "; } cout << ") Q = " << Q.val << endl;
				
				
				/*Q.remove(3,3,0);
				Q.GI = W.induced_graph( Q.node_to_comm );
				cout << "AFTER REMOVAL" << endl;
				Q.eval();
				Q.print_basic();
				double Qbefore = Q.val;
				
				Q.insert(3,5,0);
				Q.GI = W.induced_graph( Q.node_to_comm );
				cout << "AFTER INSERTION" << endl;
				Q.eval();				
				Q.print_basic();
				double Qafter = Q.val;
				
				cout << "Qbefore " << Qbefore << " Qafter " << Qafter << " diff = " << Qafter - Qbefore << endl;*/
			}
		}
	} else {
		cout << "Options are block" << endl;
	}


	return 0;
}
