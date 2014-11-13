import heapq
import itertools
import sys
import fileinput


# returns (index of beginning of largest overlap, length of largest overlap)
def overlap(s1, s2):
	if (s1 == s2):
		return (0, len(s1))

	# len(s1) must be >= len(s2)
	if len(s2) > len(s1):
		s1, s2 = s2, s1
	

	# check for full overlaps (#shorter read is fully contained within longer read)
	for i in range(0, len(s1)-len(s2)+1):
		if s1[i:i+len(s2)] == s2:
			return (i, len(s2))
	
	largest_overlap = (None, 0) #(index of beginning of overlap, length of overlap)
			
	# check for overlaps on right end of longer read
	for i in range(len(s1)-len(s2)+1, len(s1)):
		if s1[i:] == s2[:len(s1)-i]:
			if (len(s1)-i > largest_overlap[1]):
				largest_overlap = (i, len(s1)-i)

	#check for overlaps on left end of longer read
	for i in range(-len(s2)+1, 0):
		if s1[:len(s2)+i] == s2[-i:]:
			if (len(s2[-i:]) > largest_overlap[1]):
				largest_overlap = (i, len(s2[-i:]))

	return largest_overlap

# overlap_index corresponds to s1
def merge(s1, s2, overlap_index):
	if (s1 == s2):
		return s1

	if overlap_index >= 0 and overlap_index <= len(s1)-len(s2):
		return s1

	if overlap_index > len(s1)-len(s2):
		return s1 + s2[len(s1)-overlap_index:]

	if overlap_index < 0:
		return s2[:-overlap_index] + s1

# short_reads = []
# for line in fifleinput.input():
# 	line = line.strip()
# 	short_reads.append(line)
# short_reads = sys.stdin.read().splitlines()

with open('reads13.txt') as f:
    short_reads = f.read().splitlines()

# short_reads = ["AGTCA", "ACAGT", "GTCAGTCC", "GTCCCA", "TCAGTC"]
short_reads = set(short_reads)
read_pairs = set()
pairs = itertools.combinations(short_reads, 2)
for pair in pairs:
	if len(pair[1]) > len(pair[0]):
		pair = (pair[1], pair[0])
	overlaps = overlap(pair[0], pair[1])
	overlap_length = overlaps[1]
	overlap_index = overlaps[0]
	read_pairs.add((overlap_length, (pair[0], pair[1], overlap_index)))

while(len(short_reads) > 1):
	most_overlapped_pair = max(read_pairs)
	read1 = most_overlapped_pair[1][0]
	read2 = most_overlapped_pair[1][1]
	overlap_index = most_overlapped_pair[1][2]
	short_reads.remove(read1)
	short_reads.remove(read2)
	read_pairs = set([pair for pair in read_pairs if pair[1][0] != read1 and pair[1][0] != read2 and pair[1][1] != read1 and pair[1][1] != read2])
	merged_read = merge(read1, read2, overlap_index)
	for read in short_reads:
		pair = (merged_read, read)
		if len(pair[1]) > len(pair[0]):
			pair = (pair[1], pair[0])
		overlaps = overlap(pair[0], pair[1])
		overlap_length = overlaps[1]
		overlap_index = overlaps[0]
		read_pairs.add((overlap_length, (pair[0], pair[1], overlap_index)))
	short_reads.add(merged_read)

(result,) = short_reads
sys.stdout.write("%s\n" % result)



