// File: modularity.h
// -- quality functions (for Modularity criterion) header file
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


#ifndef BARBERMODULARITY_H
#define BARBERMODULARITY_H

#include "quality.h"

using namespace std;


class BarberModularity: public Quality {
 public:



  //double Nleft, Nright;
  long double F;
  vector<long double> in, tot_left, tot_right; // used to compute the quality participation of each community
  int verbose;
  
  BarberModularity(BinaryGraph & gr, vector<long double> &q_in, vector<long double> &d_in);
  ~BarberModularity();

  inline void remove(int node, int comm, long double dnodecomm);

  inline void insert(int node, int comm, long double dnodecomm);

  inline long double gain(int node, int comm, long double dnodecomm, long double w_degree);

  inline long double agg_gain(int comm1, int comm2, long double v_12);
  inline void agg(int comm1, int comm2, long double v_12);
  
  long double quality();
  
};


inline void
BarberModularity::remove(int node, int comm, long double dnodecomm) {
  assert(node>=0 && node<size);

  in[comm]  -= 2.0L*dnodecomm + g.nb_selfloops(node);
  if( tot_left[comm] > 0 && node < Nleft){tot_left[comm] -= qin[node];} 
  if( tot_right[comm] > 0 && node >= Nleft){tot_right[comm] -= din[node-Nleft];}
  
  n2c[node] = -1;
    
}

inline void
BarberModularity::insert(int node, int comm, long double dnodecomm) {
  assert(node>=0 && node<size);

  in[comm]  += 2.0L*dnodecomm + g.nb_selfloops(node);
  if( node < Nleft){ tot_left[comm] += qin[node]; }
  if( node >= Nleft){tot_right[comm] += din[node-Nleft]; }
  
  n2c[node] = comm;

}

inline long double
BarberModularity::gain(int node, int comm, long double dnc, long double degc) {
  assert(node>=0 && node<size);

  long double m   = g.total_weight;
  if( node < Nleft ){
	return (dnc - (tot_right[comm]*qin[node]/m))/m;
  } else {
	return (dnc - (tot_left[comm]*din[node-Nleft]/m))/m;
  }
  cerr << "logic error!" << endl;
  exit(1);
  return 0;
}


inline long double
BarberModularity::agg_gain(int comm1, int comm2, long double v_12) {
  long double m   = g.total_weight;  
  return (v_12/m) - ((tot_left[comm1]*tot_right[comm2] + tot_left[comm2]*tot_right[comm1] )/(m*m));
}

inline void
BarberModularity::agg(int comm1, int comm2, long double v_12) {
  
  in[comm2]  += in[comm1] + 2*v_12;
  in[comm1] = 0;

  tot_left[comm2] += tot_left[comm1];
  tot_left[comm1] = 0;
 
  tot_right[comm2] += tot_right[comm1];
  tot_right[comm1] = 0;
   
  for(unsigned int i=0; i<c2n[comm1].size(); ++i){ 
	  n2c[ c2n[comm1][i] ] = comm2; 
	  c2n[comm2].push_back(c2n[comm1][i]); 
  }
  c2n[comm1].resize(0);
  
}

#endif // BARBERMODULARITY_H
