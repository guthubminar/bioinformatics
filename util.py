#!/usr/bin/env python
# Author: Rajendran Panda, Jan. 1, 2015

class BitStore :
	def __init__(self, numBits = 1<<62) :
		self.numBits = numBits
		numIntegers = numBits/30
		reminder = numBits%30
		if reminder > 0 : numIntegers += 1
		self.Store = []
		for i in range(numIntegers) :
			self.Store.append(0)
		print self.Store

	def setBit(self, bitNum, value=1) :
		if bitNum > self.numBits :
			print 'Illegal Bit Index : ', bitNum
			exit() 
		intNum = bitNum/30
		bit = bitNum % 30
		if value == 1: 
			mask = 1<<(30-bit)
			self.Store[intNum] = self.Store[intNum] | mask
		else : 
			val = self.getBit(bitNum)
			if val == value : return
			else : self.Store[intNum] -= 1<<(30-bit)
		print self.Store[intNum]

	def getBit(self, bitNum) :
		if bitNum > self.numBits :
			print 'Illegal Bit Index : ', bitNum
			exit() 
		intNum = bitNum/32
		bit = bitNum % 32
		x = self.Store[intNum] and 1<<(30-bit)
		if x == 1<<(30-bit) : return 1
		else : return 0

#--------Unit test----------------
S = BitStore(247)

S.setBit(122, 1)
S.setBit(102)
S.setBit(133, 0)
S.setBit(65)

print S.getBit(122)
print S.getBit(102)
print S.getBit(133)
print S.getBit(113)
print S.getBit(65)
print S.getBit(300)
