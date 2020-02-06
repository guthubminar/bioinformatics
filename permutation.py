#!/usr/bin/env python

def printPermutation(perm) :
	if perm[0] > 0 : 
		s = '(+' + str(perm[0])
	else : 
		s = '(' + str(perm[0])

	for i in perm[1:] :
		if i > 0 : 
			s = s + ' +' + str(i)
		else : 
			s = s + ' ' + str(i)
	s = s + ')'
	print s

def PermutationReversals(perm) :
# uses a greedy algorithm to apply reversals to a permutation to get identity permutation
	perm = [0] + perm  # pad with a zero for convenience of index
	for k in range(1, len(perm)) :
		if perm[k] == k : continue
		if perm[k] == -k :
			perm[k] = k
			printPermutation(perm[1:])
			continue
		for j in range(k, len(perm)) :
			if perm[j] == k or perm[j] == -k : 
				subPerm = perm[j:k-1:-1]
				subPerm = [-x for x in subPerm]
				perm = perm[:k] + subPerm + perm[j+1:]
				printPermutation(perm[1:])
				break
		if perm[k] == -k : 
				perm[k] = k
				printPermutation(perm[1:])

def CountBreakpoints(perm) :
	n = len(perm)
	perm = [0] + perm + [n+1]
	bp = 0
	for k in range(n+1) :
		if perm[k+1] - perm[k] != 1 : bp += 1

	print bp
#--------------------------------------------------------------------------------------------
fn = raw_input('File :')
fp = open(fn, 'r')
perm = fp.readline().rstrip()
length = len(perm)
perm = perm[1:length-1].split(' ')
perm = [int(x) for x in perm]
# PermutationReversals(perm)
CountBreakpoints(perm)
