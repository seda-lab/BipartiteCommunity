#include <iostream>
#include "quality.h"
#include "modularity.h"
#include "barber_modularity.h"
#include "projected_modularity.h"


using namespace std; 

int main(int argc, char **argv){
	
	if(argc==1){ cout << "Options are modularity, barber, projected" << endl; }
	else if(string(argv[1]) == "modularity"){
		vector< string > filenames = { "../test_data/graph1.txt", "../test_data/graph2.txt", "../test_data/graph3.txt"};

		for(auto &filename : filenames){
			cout << "\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\" << endl;
			cout << "\\CHECKING " << filename << endl;
			cout << "\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\" << endl;
			
			/* graph1
			 * in[0] = 2
			 * in[1] = 2
			 * in[2] = 2
			 * tot[0] = 4
			 * tot[1] = 6
			 * tot[2] = 4
			 */ 
		    /* graph2
			 * in[0] = 6
			 * in[1] = 6
			 * in[2] = 4
			 * tot[0] = 7
			 * tot[1] = 8
			 * tot[2] = 5
			 */ 
			/* graph3
			 * in[0] = 2
			 * in[1] = 2
			 * in[2] = (edge)+(self) = 2+2=4
			 * tot[0] = 4
			 * tot[1] = 6
			 * tot[2] = 2+4 = 6
			 */ 
			WeightedNetwork W(filename.c_str());
			vector<int> x;
			if(W.num_nodes==6){x={0,1,0,1,2,2};}
			else{x={0,0,0,1,1,1,2,2};}
			W.print_basic();
			Modularity<WeightedNetwork> Q(W,x);			
			Q.print_basic();
			
			WeightedNetwork GI = Q.G.induced_graph( Q.node_to_comm );		
			

			for(int node=0; node<W.num_nodes; ++node){
				vector<double> dc = Q.links_to_comms(node);
				cout << node << " is attached to communities (0,1,2) with weight ( ";
				for(unsigned j=0; j<dc.size(); ++j){
					cout << dc[j] << " ";
				} cout << ")" << endl;
				
				int comm = Q.node_to_comm[node];

				Q.remove(node, comm, dc[comm]);
				double Qbefore = Q.eval() + (Q.G.selfloop(node) - Q.G.degrees[node]*Q.G.degrees[node]/Q.G.num_links )/Q.G.num_links; //add isolated node				
				double dq = Q.gain(node, comm, dc[comm]);
				Q.insert(node, comm, dc[comm]);
				double Qafter = Q.eval();

				cout << "node " << node << " gain calc = " << dq*2/Q.G.num_links << " == " << Qafter-Qbefore << " = Qdiff" << endl;
			}

			
			double Qbefore = Q.eval();
			int comm1 = 0;
			int comm2 = 1;
			double v12 = 0; for(unsigned j=0; j<GI.get_nbrs[comm1].size(); ++j){ if(j==unsigned(comm2)){ v12 = GI.get_nbrs[comm1][comm2].second;} }
			double jg = Q.join_gain(comm1, comm2, v12);
			Q.join(comm1, comm2, v12);
			double Qagg = Q.eval();	
			cout << "join " << comm1 << " and " << comm2 << " gain calc = " << 2*jg/Q.G.num_links << " == " << Qagg-Qbefore << " = Qdiff" << endl;
			
			//works because no connection between 0 and 2, otherwise calc induced again
			Qbefore = Qagg;
			comm1 = 1;
			comm2 = 2;
			v12 = 0; for(unsigned j=0; j<GI.get_nbrs[comm1].size(); ++j){ if(j==unsigned(comm2)){ v12 = GI.get_nbrs[comm1][comm2].second;} }
			jg = Q.join_gain(comm1, comm2, v12);
			Q.join(comm1, comm2, v12);
			Qagg = Q.eval();	
			cout << "join " << comm1 << " and " << comm2 << " gain calc = " << 2*jg/Q.G.num_links << " == " << Qagg-Qbefore << " = Qdiff" << endl;
			
		}
		
		
		
		
	} else if(string(argv[1]) == "barber"){
		
		vector< string > filenames = { "../test_data/bipartite1.txt", "../test_data/bipartite2.txt"};

		for(auto &filename : filenames){
			cout << "\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\" << endl;
			cout << "\\CHECKING " << filename << endl;
			cout << "\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\" << endl;
			
			/* bipartite1
			 * in[0] = 3
			 * in[1] = 1
			 * in[2] = 3
			 * tot[0] = 4
			 * tot[1] = 3
			 * tot[2] = 3
			 * tot[3+0] = 3
			 * tot[3+1] = 2
			 * tot[3+2] = 5
			 */ 
			BipartiteNetwork W(filename.c_str());
			vector<int> x;
			if(W.num_nodes==10){x={0,0,0,1,2,2,0,1,2,2};}
			      //0 1 2 3 4 5 6 7 8 9 10  11 12  13,14,15
			else{x={0,0,0,1,2,2,2,2,2,2,2  ,0, 1,  2,2,2};}
			W.print_basic();
			BarberModularity Q(W,x);			
			Q.print_basic();
		
			for(int node=0; node<W.num_nodes; ++node){
				vector<double> dc = Q.links_to_comms(node);
				cout << node << " is attached to communities (0,1,2) with weight ( ";
				for(unsigned j=0; j<dc.size(); ++j){
					cout << dc[j] << " ";
				} cout << ")" << endl;
				
				int comm = Q.node_to_comm[node];

				Q.remove(node, comm, dc[comm]);
				double Qbefore = Q.eval(); //isolated node contributes 0 to Qb				
				double dq = Q.gain(node, comm, dc[comm]);
				Q.insert(node, comm, dc[comm]);
				double Qafter = Q.eval();

				cout << "node " << node << " gain calc = " << dq*2/Q.G.num_links << " == " << Qafter-Qbefore << " = Qdiff" << endl;
			}
			
			WeightedNetwork GI = Q.G.induced_graph( Q.node_to_comm );		

			double Qbefore = Q.eval();
			int comm1 = 0;
			int comm2 = 1;
			double v12 = 0; for(unsigned j=0; j<GI.get_nbrs[comm1].size(); ++j){ if(j==unsigned(comm2)){ v12 = GI.get_nbrs[comm1][comm2].second;} }
			double jg = Q.join_gain(comm1, comm2, v12);
			Q.join(comm1, comm2, v12);
			double Qagg = Q.eval();	
			cout << "join " << comm1 << " and " << comm2 << " gain calc = " << jg*2/Q.G.num_links << " == " << Qagg-Qbefore << " = Qdiff" << endl;
			
			Qbefore = Qagg;
			comm1 = 1;
			comm2 = 2;
			v12 = 0; for(unsigned j=0; j<GI.get_nbrs[comm1].size(); ++j){ if(j==unsigned(comm2)){ v12 = GI.get_nbrs[comm1][comm2].second;} }
			jg = Q.join_gain(comm1, comm2, v12);
			Q.join(comm1, comm2, v12);
			Qagg = Q.eval();	
			cout << "join " << comm1 << " and " << comm2 << " gain calc = " << jg*2/Q.G.num_links << " == " << Qagg-Qbefore << " = Qdiff" << endl;
			
		}
	} else if(string(argv[1]) == "projected"){
		vector< string > filenames = { "../test_data/bipartite1.txt"};

		for(auto &filename : filenames){
			cout << "\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\" << endl;
			cout << "\\CHECKING " << filename << endl;
			cout << "\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\" << endl;
			
			/* bipartite1 ->
			 * 
			 * 0--2---3--[2]--4
			 *  \ |    \      |
			 *   1      5-----
			 * 
			 * in[0] = 2
			 * in[1] = 2
			 * in[2] = 2
			 * tot[0] = 2
			 * tot[1] = 5
			 * tot[2] = 3
			 */ 
			BipartiteNetwork W(filename.c_str());
			WeightedProjector P(W, true, false);
			WeightedNetwork WP = P.project();
			WP.print_basic(true);
			vector<int> x;
			if(WP.num_nodes==6){x={0,0,1,1,2,2};}
			ProjectedModularity<WeightedProjector> Q(P, WP, x);
			Q.print_basic();

			for(int node=0; node<WP.num_nodes; ++node){
				vector<double> dc = Q.links_to_comms(node);
				cout << node << " is attached to communities (0,1,2) with weight ( ";
				for(unsigned j=0; j<dc.size(); ++j){
					cout << dc[j] << " ";
				} cout << ")" << endl;
				
				int comm = Q.node_to_comm[node];

				Q.remove(node, comm, dc[comm]);
				double Qbefore = Q.eval() + (Q.G.selfloop(node) - Q.P.G.degrees[P.start+node]*Q.P.G.degrees[P.start+node]*Q.dsum )/Q.G.num_links; //add isolated node				
				double dq = Q.gain(node, comm, dc[comm]);
				Q.insert(node, comm, dc[comm]);
				double Qafter = Q.eval();

				cout << "node " << node << " gain calc = " << dq*2/Q.G.num_links << " == " << Qafter-Qbefore << " = Qdiff" << endl;
			}
			
			WeightedNetwork GI = Q.G.induced_graph( Q.node_to_comm );		

			double Qbefore = Q.eval();
			int comm1 = 0;
			int comm2 = 1;
			double v12 = 0; for(unsigned j=0; j<GI.get_nbrs[comm1].size(); ++j){ if(j==unsigned(comm2)){ v12 = GI.get_nbrs[comm1][comm2].second;} }
			double jg = Q.join_gain(comm1, comm2, v12);
			Q.join(comm1, comm2, v12);
			double Qagg = Q.eval();	
			cout << "join " << comm1 << " and " << comm2 << " gain calc = " << jg*2/Q.G.num_links << " == " << Qagg-Qbefore << " = Qdiff" << endl;
			
			Qbefore = Qagg;
			comm1 = 1;
			comm2 = 2;
			v12 = 0; for(unsigned j=0; j<GI.get_nbrs[comm1].size(); ++j){ if(j==unsigned(comm2)){ v12 = GI.get_nbrs[comm1][comm2].second;} }
			jg = Q.join_gain(comm1, comm2, v12);
			Q.join(comm1, comm2, v12);
			Qagg = Q.eval();	
			cout << "join " << comm1 << " and " << comm2 << " gain calc = " << jg*2/Q.G.num_links << " == " << Qagg-Qbefore << " = Qdiff" << endl;
			
			
		}
	} else {
		cout << "Options are modularity, barber, projected" << endl;
	}


	return 0;
}
