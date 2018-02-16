import sys, pdb
blocks = []
blockStart = ""

with open(sys.argv[1]) as infile:
	for line in map(lambda x: x.rstrip().split(), filter(lambda x: x[0] != "#", infile)):
		line[1] = int(line[1])
		if not blockStart:
			blockStart = (line[0], line[1]) # initialize block
			blockLast = (line[0], line[1])	# last element in block
			variants = 1
		for rec in line[9:]:
			dp = int(rec.split(":")[1])
			if dp <= 4 or rec[:3] == "./." or rec[:3] == "0/1" or line[0] != blockLast[0]:
				# there is an ambiguity or unphased genotype, or chromsomes crossed
				if blockStart != blockLast:
					blocks.append((blockStart, blockLast[1], variants))
				blockStart = (line[0], line[1])
				blockLast = (line[0], line[1])
				variants = 1
				break
		else:
			blockLast = (line[0], line[1])
			variants += 1

for block in blocks:
	if block[1] - block[0][1] > 2000:
		print("{}\t{}\t{}\t{}".format(block[0][0], block[0][1], block[1], block[2]))
