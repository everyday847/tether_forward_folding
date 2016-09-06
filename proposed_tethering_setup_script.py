#!/usr/bin/python
import subprocess
import sys
import argparse
import os

#print os.environ[ 'ROSETTA_TOOLS']


"""
Note that if your BPs are both from one PDB, then all you have to do is just
pass the one PDB twice. This is the simplest interface that allows you to pass
two PDBs, though!
"""

out = open( "README_TETHER_SETUP", "w" )

class Helix:
	def __init__( self ):
		self.seq1 = "" # built forwards
		self.seq2 = "" # built backwards
		self.nums1 = []
		self.nums2 = []
		self.len  = 0 

def extract_base_pairs_from_pdbs( pdbs, bps, whole_ribosome ):
	def matches( pdb_line, bp_specification ):
		"""
		A PDB line has chain and resnum to match a base pair with.
		TODO: also match segid, icode.
		"""
		#print bp_specification
		each_base = bp_specification.split('-')
		#print "each_base: %s" %each_base
		#if len(pdb_line) < 21:
		#	print pdb_line.strip()
		ch  = pdb_line[21]
		num = pdb_line[22:26].strip()
		if ch == each_base[0][0] and num == each_base[0][1:]:
			return True
		if ch == each_base[1][0] and num == each_base[1][1:]:
			return True
		return False
		

	bpfile = open( "bps.pdb", "w" ) 
	if whole_ribosome is not None:
		ribfile = open( "rib.pdb", "w" ) 
	if len(pdbs) != 2:
		print "wrong number pdbs"
		exit(1)
	if len(bps) != 2:
		print "wrong number bps"
		exit(1)
	for i in (0,1):
		pdb = open( pdbs[i] ).readlines()
		bp = bps[i]
		# first write BPs to rib.pdb, if needed
		if whole_ribosome is not None:
			for line in pdb:
				if line[0:6] == 'ATOM  ' or line[0:6] == 'HETATM':
					if matches( line, bp ):
						line = line[0:21] + 'L' + line[22:]
						ribfile.write(line)
		# then everything else
		for line in pdb:
			if line[0:6] == 'ATOM  ' or line[0:6] == 'HETATM':
				if matches( line, bp ):
					line = line[0:21] + 'L' + line[22:]
					bpfile.write(line)
				elif whole_ribosome is not None:
					if i == 0:
						ribfile.write(line[0:21]+'A'+line[22:])
					else: #if i == 1:
						ribfile.write(line[0:21]+'B'+line[22:])


def helices_from_secstruct_and_fasta( sec_struct, fasta ):
	"""
	Given a secondary structure and fasta, give all the helices
	"""
	helices = []	
	if len(sec_struct) != len(fasta):
		print "Error: sec struct and fasta of different lengths!"
		exit(1)
	#start from beginning and end of secstruct
	begin_iter = 0
	end_iter = len( sec_struct ) - 1 
	residues_left = True
	helix_temp = Helix()
	inside_a_helix = False
	while begin_iter < end_iter:
		if sec_struct[begin_iter] == '(' and sec_struct[end_iter] == ')':
			inside_a_helix = True
			helix_temp.seq1 = helix_temp.seq1   + fasta[begin_iter]
			helix_temp.seq2 = fasta[end_iter]   +   helix_temp.seq2 
			# second helix numbers are offset due to the single char
			# separator, while we need to offset first helix ourselves
			helix_temp.nums1 = helix_temp.nums1 + [begin_iter+1]
			helix_temp.nums2 = [end_iter]   + helix_temp.nums2
			helix_temp.len += 1
		else:
			if inside_a_helix:
				inside_a_helix = False
				helices.append( helix_temp )
				helix_temp = Helix()
		
		# increment/decrement logic
		if inside_a_helix:
			begin_iter += 1
			end_iter -= 1
		elif sec_struct[begin_iter] != '(' and sec_struct[end_iter] != ')':
			begin_iter += 1
			end_iter -= 1
		elif sec_struct[begin_iter] == '(':
			end_iter -= 1
		elif sec_struct[end_iter] == ')':
			begin_iter += 1

	print fasta
	print sec_struct
	for helix in helices:
		print helix.seq1, helix.seq2			
		print helix.nums1, helix.nums2			
	return helices

def parse_fasta( fasta_file ):
	"""
	Parse a fasta file; at the moment, one that we assume is one line long.
	(That means no name line is provided -- this is fine, because we will
	generate that type of fasta file in setting up the job.
	"""
	for line in fasta_file:
		return line.strip()

parser = argparse.ArgumentParser(description='Obtain initial helices for de novo.')
parser.add_argument('--fasta', dest='fasta', nargs='?', help='target fasta')
parser.add_argument('--fasta_file', dest='fasta_file', nargs='?', help='path to a fasta file')
parser.add_argument('--secstruct', dest='secstruct', nargs='?', help='target secstruct')
parser.add_argument('--secstruct_file', dest='secstruct_file', nargs='?', help='path to a secstruct file')

parser.add_argument('--pdbs', dest='pdbs', help='the two PDB files (possibly the same file!) in which the starting and ending base pair are found, respectively')
parser.add_argument('--bps', dest='bps', help='string describing the base pairs for start and end, in each pdb\'s numbering (A56-B23 C41-D93); add icode and segid as needed after underscores') 
parser.add_argument('--execute', dest='execute', help='pass any value if you want the setup script to run automatically') 
parser.add_argument('--sep', dest='sep', help='sequence separator used in fasta (default space)' )
parser.add_argument('--whole_ribosome', dest='whole_ribosome', help='retain the whole ribosome in the input')
parser.add_argument('--vdw_screen', dest='vdw_screen', help='use a van der Waals screen strategy for the region around the base pairs')

args = parser.parse_args()

pdbs_list = args.pdbs.split()
bps_list = args.bps.split()
sep = args.sep
if sep == None: sep = ' '
fasta_file = None
sec_struct_file = None
fasta = None
sec_struct = None
if args.fasta is not None:
	#print args.fasta
	fasta = args.fasta
if args.fasta_file is not None:
	#print args.fasta_file
	try:
		fasta_file = open(args.fasta_file)
	except:
		print "Error: fasta file provided not found."
		exit
	fasta = parse_fasta( fasta_file )
if args.secstruct is not None:
	#print args.fasta
	sec_struct = args.secstruct
if args.secstruct_file is not None:
	#print args.fasta_file
	try:
		secstruct_file = open(args.secstruct_file)
	except:
		print "Error: secondary structure file provided not found."
		exit(1)
	sec_struct = parse_secstruct( secstruct_file )
	 
assert(sec_struct is not None and fasta is not None)

helices = helices_from_secstruct_and_fasta( sec_struct, fasta )


# cleave base pairs based on CL specification
extract_base_pairs_from_pdbs( pdbs_list, bps_list, args.whole_ribosome )
#extract_surroundings_of_bps_for_vdw( args.pdbs, args.bps, fasta )
# renumber base pairs accordingly
# residue numbering scheme: first bp of pdb1  is 1, sequence length
#                           second bp of pdb2 is before and after seqsep 
out.write( "renumber_pdb_in_place.py bps.pdb L:%s L:%s L:%s L:%s\n\n" % (1, len(fasta)-1, fasta.index(sep), fasta.index(sep)+1) )

input_files = [ "bps.pdb" ]
helix_index = 1
for helix in helices:
	out.write("rna_helix.linuxclangrelease -rna::corrected_geo  -score:rna_torsion_potential RNA11_based_new -chemical::enlarge_H_lj -o helix_%s.pdb -minimize_all -seq %s %s\n" % (helix_index, helix.seq1, helix.seq2) )
	# Now, we have to correct the chain and numbering
	out.write("replace_chain_inplace.py helix_%s.pdb L\n" % ( helix_index ) )
	out.write("renumber_pdb_in_place.py helix_%s.pdb L:%s-%s L:%s-%s\n\n" % ( helix_index, min(helix.nums1), max(helix.nums1), min(helix.nums2), max(helix.nums2) ) )
	input_files.append( "helix_%s.pdb" % helix_index )	
	
	if min(helix.nums1) == 1 and max(helix.nums2) == len(fasta)-1:
		# This is one bps helix
		out.write( "rna_graft.linuxclangrelease -s bps.pdb helix_%s.pdb -o new_helix_%s.pdb\n" % ( helix_index, helix_index ) ) 
		out.write( "mv new_helix_%s.pdb helix_%s.pdb\n" % ( helix_index, helix_index ) )
		out.write("replace_chain_inplace.py helix_%s.pdb L\n" % ( helix_index ) )
	if max(helix.nums1) == fasta.index(sep) and min(helix.nums2) == fasta.index(sep)+1:
		out.write( "rna_graft.linuxclangrelease -s bps.pdb helix_%s.pdb -o new_helix_%s.pdb\n" % ( helix_index, helix_index ) ) 
		out.write( "mv new_helix_%s.pdb helix_%s.pdb\n" % ( helix_index, helix_index ) )
		out.write("replace_chain_inplace.py helix_%s.pdb L\n" % ( helix_index ) )
	
	helix_index += 1



# graft any helix that includes a BP residue (i.e. first or last of fasta) so that rigid body alignment is maintained
helix_index = 1
for helix in helices:
	print helix.nums1, helix.nums2
	helix_index += 1
	
outfasta = open( "input.fasta", "w" )
outfasta.write( "> input L:%s-%s\n%s\n" % ( 1, len(fasta)-1, fasta ) )
outfasta.close()

# Some of these will be added as customizable, but really you can edit the
# flags yourself	#	
flags = open( "auto_flags", "w" )
flags.write( "-nstruct 1000\n" )
flags.write( "-cycles 100\n" )
flags.write( "-secstruct '%s'\n" % sec_struct )
flags.write( "#          '%s'\n" % fasta )
flags.write( "-fasta input.fasta\n" )
flags.write( "-out:file:silent tether.out\n" )
flags.write( "-minimize_rna true\n" )
flags.write( "-s %s\n" % ( " ".join( input_files ) ) )
flags.write( "-output_res_num L:%s-%s\n" % ( 1, len(fasta)-1 ) )
flags.write( "-ignore_zero_occupancy false\n" ) 
# no idea why these are the values
#flags.write( "-cutpoint_closed L:%s L:%s\n" % ( max(helices[0].nums1), fasta.index(sep) ) ) 
# This seems better supported
flags.write( "-chain_connection L:%s L:%s L:%s L:%s\n" % ( 1, fasta.index(sep), fasta.index(sep)+1, len(fasta)-1 ) )

#-input_res L:1-13 L:89-101 L:20-31 L:39-44 L:49-50 L:51-52 L:61-66 L:71-82

if args.whole_ribosome is not None:
	flags.write( "-grid_vdw_weight 0.3\n" )
	flags.write( "-VDW_rep_screen_include_sidechains true\n" )
	# align to first set of bps, IMO
	# make sure rib gets the bps in it, and they're written out first!
	flags.write( "-VDW_rep_screen_info vdw_screen_minimal.pdb 1-2 1-2\n" )

flags.close()
out.close()

if args.execute is not None:
	print "Executing README_TETHER_SETUP"
	subprocess.call( [ "sh", "README_TETHER_SETUP" ] )

READMESETUP = open( "README_SETUP", "w" )
READMESETUP.write( "rosetta_submit.py command.sh OUT 256 48 -save_logs\n" )
READMESETUP.close()
COMMAND = open( "command.sh", "w" )
COMMAND.write( "/biox3/home/amw579/Rosetta/main/source/bin/rna_denovo.cxx11.linuxclangrelease  @ auto_flags\n" )
COMMAND.close()
