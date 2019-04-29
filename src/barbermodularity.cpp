// File: modularity.cpp
// -- quality functions (for Newman-Girvan Modularity criterion) source file
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


#include "barbermodularity.h"

using namespace std;


BarberModularity::BarberModularity(BinaryGraph & gr, vector<long double> &q_in, vector<long double> &d_in):Quality(gr,"Barber Modularity") {
  
  verbose = 1;
  n2c.resize(size);

  in.resize(size);
  tot_left.resize(size);
  tot_right.resize(size);

  qin = q_in;
  din = d_in;
  Nleft = qin.size();
  Nright = din.size();

  if( size != Nleft+Nright ){
	  cerr << "bipartite size error " << size << " " << Nleft << " " << Nright << endl;
	  exit(1);
  }

  // initialization
  F = 0;
  for (int i=0 ; i<size ; i++) {
    n2c[i] = i;
    in[i]  = g.nb_selfloops(i);
    if(i<Nleft){
		tot_left[i] = qin[i]; 
		tot_right[i] = 0;
		F += qin[i];
	} else {
		tot_left[i] = 0; 
		tot_right[i] = din[i-Nleft];
	}
	
  }

  
}

BarberModularity::~BarberModularity() {
  in.clear();
  tot_left.clear();
  tot_right.clear();
}

long double
BarberModularity::quality() {
  long double qp  = 0.0L;
  long double m = g.total_weight;

  for (int i=0 ; i<size ; i++) {
      qp += (1.0/m)*(in[i]/2.0 - (tot_left[i]*tot_right[i]) / m);
  }

  return qp;
  
}

