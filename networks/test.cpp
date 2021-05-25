#include <iostream>
#include "network.h"
#include "weighted_network.h"
#include "bipartite_network.h"
#include "projector.h"

using namespace std; 

int main(int argc, char **argv){

	if(argc==1){ cout << "Options are induced, bipartite, projector" << endl; }
	else if(string(argv[1]) == "induced"){
		vector< string > filenames = { "../test_data/graph1.txt", "../test_data/graph2.txt", "../test_data/graph3.txt"};

		for(auto &filename : filenames){
			cout << "\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\" << endl;
			cout << "\\CHECKING " << filename << endl;
			cout << "\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\" << endl;
			//Weighted
			WeightedNetwork W(filename.c_str());
			W.print_basic();
			if(W.num_nodes == 6){
				vector<int> x{0,0,0,1,1,1};
				WeightedNetwork WI = W.induced_graph(x);
				WI.print_basic();
				/* graph1 ->    (6)-0--[1]--1-(6)
				 * graph3 ->    (6)-0--[1]--1-(8)
				 */
			} else if(W.num_nodes == 8){
				vector<int> x{0,0,0,1,1,1,2,2};
				WeightedNetwork WI = W.induced_graph(x);
				WI.print_basic();
				/* graph2 ->    (6)-0--[1]--1-(6)
				 * 							|
				 * 							2-(4)
				 */
			 } 
		}
	} else if(string(argv[1]) == "bipartite"){
		vector< string > filenames = {"../test_data/graph1.txt" , "../test_data/bipartite1.txt",  "../test_data/bipartite2.txt"};

		for(auto &filename : filenames){
			cout << "\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\" << endl;
			cout << "\\CHECKING " << filename << endl;
			cout << "\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\" << endl;
			//Weighted
			BipartiteNetwork B(filename.c_str());
			if(!B.is_bipartite()){ cout << "This is not a bipartite network" << endl; continue; }
			B.print_basic();
			if(B.num_nodes == 10){
				vector<int> x{0,0,0,1,1,1};
				BipartiteNetwork BI = B.induced_graph_side(x, true); //left induce
				BI.print_basic();
				vector<int> y{0,0,1,1};						
				BipartiteNetwork BIr = B.induced_graph_side(y, false); //right induce
				BIr.print_basic();
				/* bipartite1 left ->    
				 *    0--[3]--6
				 *     \_____7
				 *    1_____/
				 * 	  |\
				 *    | \__[2]_8
				 *    |
				 * 	  |__[3]__9
				 * 
				 */
				/* bipartite1 right ->    
				 * 0-------10--
				 * 1------/ | |
				 * 2--[2]--/ /
				 *           |
				 * 3---------|
				 *  \___[2]___
				 * 			  11
				 * 4----[2]--/ |
				 * 5-----------/
				 * 
				 */
			} else if(B.num_nodes == 16){
				vector<int> x{0,0,0,0,1,1,1,1,1,1,1};				
				BipartiteNetwork BI = B.induced_graph_side(x, true); //left induce
				BI.print_basic();
				vector<int> y{0,0,1,1,1};										
				BipartiteNetwork BIr = B.induced_graph_side(y, false); //right induce
				BIr.print_basic();
			}
		}
	} else if(string(argv[1]) == "projector"){
		vector< string > filenames = { "../test_data/bipartite1.txt"};

		for(auto &filename : filenames){
			cout << "\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\" << endl;
			cout << "\\CHECKING " << filename << endl;
			cout << "\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\" << endl;
			//Weighted
			BipartiteNetwork B(filename.c_str());
			if(!B.is_bipartite()){ cout << "This is not a bipartite network" << endl; continue; }
			B.print_basic();
			
			
			/* LEFT, NOLOOP
			 * dsum = 6+2+2+6 = 16
			 * 
			 * 0--2---3--[2]--4
			 *  \ |    \      |
			 *   1      5-----
			 *   
			 * LEFT, LOOP
			 * dsum = 9+4+4+9 = 26
			 * Loops : {0:1,1:1,2:2,3:3,4:2,5:1}
			 * 
			 * 
			 * RIGHT, NOLOOP
			 * dsum = 0+0+2+6+2+0 = 10
			 * 6--7--8
			 *    |  |
			 *    | [2]
			 *    |  |
			 *     \ |
			 *       9
			 * RIGHT, LOOP
			 * Loops : {6:3,7:2,8:2,9:3}
			 * dsum = 1+1+4+9+4+1 = 20
			 * 
			 */
			
			cout << "Weighted Projection" << endl;
			for(unsigned left=0; left<2;++left){
				for(unsigned loops=0; loops<2;++loops){
					cout << "=============================================" << endl;
					if(left){ cout << "left, "; }
					else{ cout << "right, "; }
				
					if(loops){ cout << "include loops, "; }
					else{ cout << "don't include loops, "; }
					
					WeightedProjector WP(B, bool(left), bool(loops) );
					cout << " dsum = " << WP.dsum() << endl;
					
					WeightedNetwork GP = WP.project();
					GP.print_basic();					
					cout << "=============================================" << endl;
				}
			}
			
			do{
				cout << "Binary Projection" << endl;
				BinaryProjector BP(B);
				cout << " dsum = " << BP.dsum() << endl;	
				WeightedNetwork GP = BP.project();
				GP.print_basic();
			} while(0);
			
			do{
			/* Left, no loops, weights = {2, 1, 1, 2}
			 * 0-1 = 1/2
			 * 0-2 = 1/2
			 * 1-2 = 1/2
			 * 2-3 = 1/1
			 * 3-4 = 1/1+1/2 = 3/2
			 * 3-5 = 1/2
			 * 4-5 = 1/2
			 * */	
				cout << "Hyperbolic Projection: Left, no loops" << endl;
				HyperbolicProjector HP(B);
				cout << " dsum = " << HP.dsum() << endl;	
				WeightedNetwork GP = HP.project();
				GP.print_basic();
			} while(0);
			
			do{	
			/* Left, with loops, weights = {3, 2, 2, 3}
			 * 0-1 = 1/3
			 * 0-2 = 1/3
			 * 1-2 = 1/3
			 * 2-3 = 1/2
			 * 3-4 = 1/2+1/3 = 5/6
			 * 3-5 = 1/3
			 * 4-5 = 1/3
			 * 0-0 = 1/3
			 * 1-1 = 1/3
			 * 2-2 = 1/3+1/2 = 5/6
			 * 3-3 = 1/2+1+2+1/3 = 4/3
			 * 4-4 = 1/2+1/3 = 5/6
			 * 5-5 = 1/3
			 * */	
				cout << "Hyperbolic Projection: Left, loops" << endl;
				HyperbolicProjector HP(B, true, true);
				cout << " dsum = " << HP.dsum() << endl;	
				WeightedNetwork GP = HP.project();
				GP.print_basic();
			} while(0);
			
			
		}
	} else {
		cout << "Options are induced, bipartite, projector" << endl;
	}


	return 0;
}
