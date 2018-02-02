import sys, math
from Bio import SeqIO
import pdb

def split_fasta(f, n):
	"""
	Takes a fasta file, returns list of lists, corresponding to positions that will split it into about n chunks
	The output consists of chrom:start-stop with one chunk per line
	"""
	class myrec:
		def __init__(self, name, start, stop):
			self.name = name
			self.start = start
			self.stop = stop
		def __str__(self):
			return("{}:{}-{}".format(self.name,self.start+1,self.stop))
	recs = []
	total = 0
	for rec in SeqIO.parse(f, "fasta"):
		recs.append(myrec(rec.id, 0, len(rec)))
		total += len(rec)
	
	bins = list(range(0, total, int(math.ceil(total/n))))[1:] + [total]
	current = bins.pop(0)
	# pdb.set_trace()	
	outbins = [[]]
	runningSum = 0
	for rec in recs:
		recLen = rec.stop - rec.start
		spaceLeft = current - runningSum # how much space in the bin
		while recLen > spaceLeft:						
			if bins:
				if spaceLeft > 0:
					outbins[-1].append(myrec(rec.name, rec.start, rec.start + spaceLeft))
				rec.start = rec.start + spaceLeft
				current = bins.pop(0)
				runningSum += spaceLeft
				outbins.append([]) #start new bin
				recLen = rec.stop - rec.start
				spaceLeft = current - runningSum
			else:
				outbins[-1].append(myrec(rec.name, rec.start, rec.stop))
				# at this point there are no more bins and we're done
				break
		else:
			# add rest of record to the current bin
			outbins[-1].append(myrec(rec.name, rec.start, rec.stop))
			runningSum += rec.stop - rec.start
	outdict = dict(zip(map(lambda x: str(x),range(len(outbins))), outbins))			
#	for i in outdict:
#	 	print(" ".join(map(lambda x: str(x), outdict[i])))
	return(outdict)

if __name__ == "__main__":
	split_fasta(sys.argv[1],int(sys.argv[2]))
