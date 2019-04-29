# BipartiteCommunity
Modifications of the Louvain Algorithm to detect communities in biaprtite networks.

Examples:
Python script that runs all the algorithms on a bipartite graph
```python3 run_louvain.py bipartite 10 6```
Notes: 
	Assumes file ends in .txt, don't include it!
	10 = number of random restarts of louvain
	6 = number of nodes on the left hand side

unipartite, maximise standard modularity
```./louvain -i graph1.txt -q 0```

unipartite, repeated links transformed into weights, maximise standard modularity
```./louvain -i graph2.txt -q 0```

unipartite, weirdly numbered nodes, maximise standard modularity
```./louvain -i graph3.txt -q 0```

unipartite, self loops, maximise standard modularity
```./louvain -i graph4.txt -q 0```

unipartite, weighted, maximise standard modularity
```./louvain -i weighted.txt -q 0 -w```

bipartite, maximise barber modularity
```./louvain -i bipartite.txt -q 11 -b 6```

bipartite, maximise barber modularity, try agglomerative clustering on induced graph
```./louvain -i bipartite.txt -q 11 -b 6 -a```

left projected, maximise bimodularity
```./louvain -i bipartite.txt -q 10 -b 6 -p 0```


