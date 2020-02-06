#!/usr/bin/env python

def composition(text, k):
	composition = []
	ks = len(text)
	for i in range(ks-k+1) :
		kmer = text[i:i+k]
		if kmer not in composition :
			composition.append(kmer)
	return sorted(composition)


def GenomePath2String(path) :
	path.reverse()
	text = path[0]

	for elem in path[1:] :
		text = elem[0] + text

	return text
#----------------------------------------------------
def composition_main(fp) :
	k = int(fp.readline().rstrip())
	text = fp.readline().rstrip()
	print k
	comp = composition(text, k)
	for elem in comp:
		print elem

def GenomePath2String_main(fp) :
	path = []
	while True :
		line = fp.readline()
		if not line: break
		path.append(line.rstrip())

	text = GenomePath2String(path)
	print text

#-----------------------------------------------------

fn = raw_input ('Input file : ')
fp = open(fn, 'r')
#composition_main(fp)
GenomePath2String_main(fp)
		
