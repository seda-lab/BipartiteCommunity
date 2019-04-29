import sys
import os
from subprocess import Popen, PIPE
import shlex
	
def run_louvain_command(execpath, filename, ltype, B, num_its, weighted, output=True):

	if ltype == 0: ##Naive - standard louvain, optimising Q
		cmd = execpath  + "louvain -i " + filename + ".txt -q 0 -p 0 -b " + B + " -P " + filename + "_proj0.txt";
		label = "standard"
	elif ltype == 1: 
		cmd = execpath  + "louvain -i " + filename + ".txt -q 10 -p 0 -b " + B;
		label = "proj"
	elif ltype == 2:
		cmd = execpath  + "louvain -i " + filename + ".txt -q 11 -b " + B + " -a";
		label = "barber"
	else:
		print("incorrect ltype");
		sys.exit(1);
	if weighted: cmd += " -w";
	
	best_mod = -1;
	best_out = ""
	for i in range(num_its):
	
		p = Popen( shlex.split(cmd), stdout=PIPE, stderr=PIPE)
		stdout, stderr = p.communicate();
		#print( "stdout", stdout.decode("utf-8").strip().split("\n") );
		#print( "stderr", stderr.decode("utf-8").strip().split("\n") ); 
		txt = stderr.decode("utf-8").strip().split("\n")[-1].split()[1];
		if txt == "size": 
			return None;
		mod = float( txt )
		if mod > best_mod:
			best_mod = mod;
			best_out = stdout;
	
	if output:		
		with open(filename + ".tree." + label, 'w') as outfile:
			txt = best_out.decode("utf-8").strip().split("\n");
			if txt == "No projection":
				return None;
			for line in txt:
				outfile.write("{}\n".format(line.strip()) );
		
	return best_mod;

def run_louvain_command_bi(execpath, filename, ltype, B, num_its, weighted, output=True):

	if ltype == 0: ##Naive - standard louvain, optimising Q
		cmd = execpath  + "louvain -i " + filename + ".txt -q 0";
		label = "naive"
	elif ltype == 1: 
		cmd = execpath  + "louvain -i " + filename + ".txt -q 11 -b " + B + " -a";
		label = "bilou"
	else:
		print("incorrect ltype");
		sys.exit(1);
	if weighted: cmd += " -w";

	best_mod = -1;
	best_out = ""
	for i in range(num_its):
		p = Popen( shlex.split(cmd), stdout=PIPE, stderr=PIPE)
		stdout, stderr = p.communicate();
		#print( "stdout", stdout.decode("utf-8").strip().split("\n") );
		#print( "stderr", stderr.decode("utf-8").strip().split("\n") ); 
		txt = stderr.decode("utf-8").strip().split("\n")[-1].split()[1];
		if txt == "size": 
			return None;
		mod = float( txt )		
		if mod > best_mod:
			best_mod = mod;
			best_out = stdout;
	
	if output:		
		with open(filename + ".tree."+label, 'w') as outfile:
			txt = best_out.decode("utf-8").strip().split("\n");
			if txt == "No projection":
				return None;
			for line in txt:
				outfile.write("{}\n".format(line.strip()) );
		
	return best_mod;
	
def run_louvain_command_dual_proj(execpath, filename, ltype, B, num_its, weighted, output=True):

	###################
	##left projection##
	###################
	cmd = execpath  + "louvain -i " + filename + ".txt -q 10 -b " + B + " -p 0"; 
	if weighted: cmd += " -w";

	best_mod_left = -1;
	best_out_left = ""
	for i in range(num_its):
	
		p = Popen( shlex.split(cmd), stdout=PIPE, stderr=PIPE)
		stdout, stderr = p.communicate();
		#print( "stdout", stdout.decode("utf-8").strip().split("\n") );
		#print( "stderr", stderr.decode("utf-8").strip().split("\n") ); 
		txt = stderr.decode("utf-8").strip().split("\n")[-1].split()[1];
		if txt == "size": 
			return None;
		mod = float( txt )
		if mod > best_mod_left:
			best_mod_left = mod;
			best_out_left = stdout;
	

	
	####################
	##right projection##
	####################
	##let's always use left projection
	left_ids = {};
	right_ids = {};
	with open(filename + ".txt", 'r') as infile:
		for line in infile:
			words = line.split();
			left_ids[ int(words[0]) ] = int(words[0])
			right_ids[ int(words[1]) ] = int(words[1])
	##renumber
	for i,r in enumerate(sorted(right_ids)):
		right_ids[r] = i;
	nr = len(right_ids)
	for i,l in enumerate(sorted(left_ids)):
		left_ids[l] = i+nr;
	##write the file	
	edges = {}
	with open(filename + ".txt", 'r') as infile:
		for line in infile:
			words = line.split();
			l=left_ids[ int(words[0]) ]
			r=right_ids[ int(words[1]) ]
			if min(l,r) in edges:
				edges[min(l,r)].append(max(l,r))
			else:
				edges[min(l,r)] = [max(l,r)]
			
	with open(filename + "_dp.txt", 'w') as outfile:
		for e in edges:
			for ee in sorted(edges[e]):
				outfile.write("{} {}\n".format( e, ee ) )
	
	cmd = execpath  + "louvain -i " + filename + "_dp.txt -q 10 -b " + str(len(right_ids)) + " -p 0"; 
	if weighted: cmd += " -w";

	best_mod_right = -1;
	best_out_right = ""
	for i in range(num_its):
	
		p = Popen( shlex.split(cmd), stdout=PIPE, stderr=PIPE)
		stdout, stderr = p.communicate();
		#print( "stdout", stdout.decode("utf-8").strip().split("\n") );
		#print( "stderr", stderr.decode("utf-8").strip().split("\n") ); 
		txt = stderr.decode("utf-8").strip().split("\n")[-1].split()[1];
		if txt == "size": 
			return None;
		mod = float( txt )		
		if mod > best_mod_right:
			best_mod_right = mod;
			best_out_right = stdout;

			
	comms = set();
	label = "dp1";
	with open(filename + ".tree."+label, 'w') as outfile:
		##communities found on left
		nl = 0;
		for line in best_out_left.decode("utf-8").strip().split("\n"):
			words = line.strip().split();
			comms.add( words[1] )
			outfile.write("{}\n".format(line.strip()) );
			nl += 1

		##communities found on right
		nc = len(comms);
		##old -> new ----> new -> old
		n_right_ids = { right_ids[r]:r for r in right_ids };
		for line in best_out_right.decode("utf-8").strip().split("\n"):
			nl += 1;
			words = line.strip().split();
			outfile.write("{} {}\n".format( n_right_ids[ int(words[0]) ], int(words[1])+nc) );

	###################
	##dual projection##
	###################
	starting_partition = filename + ".tree."+label
	cmd = execpath  + "louvain -i " + filename + ".txt -q 11 -b " + B + " -A -t " + starting_partition; 
	#cmd = execpath  + "louvain -i " + filename + ".txt -q 0 -t " + starting_partition ; 
	if weighted: cmd += " -w";

	best_mod = -1;
	best_out = ""
	
	p = Popen( shlex.split(cmd), stdout=PIPE, stderr=PIPE)
	stdout, stderr = p.communicate();
	#print( "stdout", stdout.decode("utf-8").strip().split("\n") );
	#print( "stderr", stderr.decode("utf-8").strip().split("\n") ); 
	mod = float( stderr.decode("utf-8").strip().split("\n")[-1].split()[1] )
	if mod > best_mod:
		best_mod = mod;
		best_out = stdout;
	
	label = "dualproj";
	if output:		
		with open(filename + ".tree."+label, 'w') as outfile:
			for line in best_out.decode("utf-8").strip().split("\n"):
				outfile.write("{}\n".format(line.strip()) );
						
	return best_mod;
	
				
def find_partition(filename, num_its, B, dostandard=True, doproj=True, dobarber=True,  weighted=False):

	num_its = int(num_its);
	B = str(B);

	execpath=os.getenv("HOME") + "/Dropbox/bipartite/louvain/v7/"

	###find communities
	mod0 = None;
	mod1 = None;
	mod2 = None;
	if dostandard:
		mod0 = run_louvain_command(execpath, filename, 0, B, num_its, weighted);
		print("done standard", mod0);
	if doproj:
		mod1 = run_louvain_command(execpath, filename, 1, B, num_its, weighted);
		print("done proj", mod1);
	if(dobarber): 
		mod2 = run_louvain_command(execpath, filename, 2, B, num_its, weighted);
		print("done barber", mod2);
	
	return (mod0, mod1, mod2);

def find_partition_bi(filename, num_its, B, donaive=True, dobilou=True, dodual=True, weighted=False):

	num_its = int(num_its);
	B = str(B);

	execpath=os.getcwd() + "/"; #os.getenv("HOME") + "louvain/"

	##find communities
	mod0 = None;
	mod1 = None;
	mod2 = None;
	if donaive:
		mod0 = run_louvain_command_bi(execpath, filename, 0, B, num_its, weighted);
		print("done naive", mod0);
	if dobilou:
		mod1 = run_louvain_command_bi(execpath, filename, 1, B, num_its, weighted);
		print("done bilou", mod1);
	if dodual:
		mod2 = run_louvain_command_dual_proj(execpath, filename, 1, B, num_its, weighted);
		print("done dual", mod2);

	return (mod0, mod1, mod2);
		
if __name__ == '__main__':
	
	print("projected")
	print( find_partition(sys.argv[1], sys.argv[2], sys.argv[3]) )
	print("unprojected")
	print( find_partition_bi(sys.argv[1], sys.argv[2], sys.argv[3]) )
	
	
	
	
