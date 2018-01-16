from Bio import AlignIO
from Bio import Phylo
from Bio import SeqIO
from Bio import pairwise2
from Bio.Align.Applications import ClustalwCommandline
from Bio.Align.Applications import MuscleCommandline
from Bio.Phylo.TreeConstruction import DistanceCalculator
import numpy as np
# from asciitree import LeftAligned
from collections import OrderedDict as OD


class Node(object):
	def __init__(self, name, data):
		self.__data = data
		self.__name = name
		self.__rchild = None
		self.__lchild = None
		self.__parent = None
	
	def isLeaf(self):
		return self.__rchild is None and self.__lchild is None
	
	def setData(self, data):
		self.__data = data
	
	def getData(self):
		return self.__data
	
	def getRchild(self):
		return self.__rchild
	
	def getLchild(self):
		return self.__lchild
	
	def setRchild(self, node):
		self.__rchild = node
		self.__rchild.setParent(self)
	
	def setLchild(self, node):
		self.__lchild = node
		self.__lchild.setParent(self)
	
	def setParent(self, node):
		self.__parent = node
		
	def getParent(self):
		return self.__parent
	def getName(self):
		return self.__name
	
	def __repr__(self):
		return self.__name


## construct tree object - constructor gets path to sequences file
class PhylTree(object):
	def __init__(self, inputSequences):
		self.__dists = {}
		seqs = []
		print("reading fastas")
		for seq in SeqIO.parse(inputSequences, "fasta"):
			seqs.append((seq.id, seq.seq._data))
		
		self.__dm = np.zeros((len(seqs), len(seqs)))
		### pairwise gives shit topos
		# for i in range(len(seqs)):
		# 	for j in range(len(seqs)):
		# 		s = pairwise2.align.globalxx(seqs[i][1], seqs[j][1], score_only=True)
		# 		s = s / min(len(seqs[i][1]), len(seqs[j][1]))
		# 		self.__dists[(seqs[i][0], seqs[j][0])] = s
		# 		# self.__dists[(seqs[j][0], seqs[i][0])] = s
		# 		self.__dm[i,j] = s
		# 		# self.__dm[j,i] = s
		# print(self.__dm)
		print("MSA")
	
		cline = ClustalwCommandline("C:\Program Files (x86)\ClustalW2\clustalw2",infile=inputSequences, outfile="outAlign.aln")
		cline()
		aln = AlignIO.read('outAlign.aln', 'clustal')
		
		print("Alignment output to 'alignment.txt'")
		with open('alignment.txt', 'w') as f:
			for s in aln._records:
				n = len(str(s.id))
				p = str(" " * (10 - n))
				f.write(s.id + p + '\t' + str(s.seq) + '\n')
		
		calculator = DistanceCalculator('identity')
		dm = calculator.get_distance(aln)
		# print(dm)
		
		self.__nodes = dm.names
		for i in range(len(dm.matrix)):
			for j in range(i + 1):
				self.__dists[(dm.names[i], dm.names[j])] = dm.matrix[i][j]
				self.__dists[(dm.names[j], dm.names[i])] = dm.matrix[i][j]
				self.__dm[i, j] = dm.matrix[i][j]
				self.__dm[j, i] = dm.matrix[i][j]
		
		nodes = []
		print("Neighbour joining")
		for leaf in aln._records:
			nodes.append(Node(leaf.id, str(leaf.seq)))
		self.__root = self.neighborJoin(nodes, self.__dm, 0)
		
	def setAllEdgeDist(self, val):
		for key in self.__dists.keys():
			self.__dists[key] = val
		

	def getRoot(self):
		return self.__root
	
	
	def updateDist(self, na, nb, s):
		self.__dists[(na.getName(), nb.getName())] = s
		self.__dists[(nb.getName(), na.getName())] = s
	
	
	def neighborJoin(self, forest, D, idx):
		rptsum = lambda arr: np.repeat(sum(arr), np.size(arr))
		mapvsum = lambda mat: np.matrix(map(rptsum, mat))
		idxmin = lambda mat: np.unravel_index(np.argmin(mat), np.shape(mat))
		wraparr = lambda x: [x]
		
		# Remove rows and columns of matrix with the listed indices
		withoutIndices = lambda m, ids: np.delete(np.delete(m, ids, axis=0), ids, axis=1)
		# Append a vector as both a row and a column
		appendRowCol = lambda m, v: np.hstack((np.vstack((m, [v])), map(wraparr, v + [0])))
		
		# def neighborJoin(D, forest, idx):
		if len(D) == 2:
			root = Node('inner' + str(idx), None)
			root.setLchild(forest[0])
			root.setRchild(forest[1])
			self.__dists[('inner' + str(idx), forest[0].getName())] = 0
			self.__dists[forest[0].getName(), ('inner' + str(idx))] = 0
			self.__dists[('inner' + str(idx), forest[1].getName())] = 0
			self.__dists[(forest[1].getName(), 'inner' + str(idx))] = 0
			
			return root
		SH = mapvsum(D)
		SV = SH.transpose()
		I = np.identity(len(D))
		M = D + (np.multiply(I, SH + SV) - SH - SV) / (len(D) - 2)
		i, j = idxmin(M)
		u = [(D[i, k] + D[j, k] - D[i, j]) / 2 for k in range(len(D))]
		# su  = (D[i,j] + SH[i,j] - SV[i,j]) / 2
		dui = (D[i, j] + SH[i, j] - SV[i, j]) / 2
		duj = (D[j, i] + SH[j, i] - SV[j, i]) / 2
		nodeU = Node('inner' + str(idx), None)
		nodeU.setLchild(forest[i])
		nodeU.setRchild(forest[j])
		self.__dists[('inner' + str(idx), forest[i].getName())] = dui
		self.__dists[(forest[i].getName(), 'inner' + str(idx))] = dui
		self.__dists[('inner' + str(idx), forest[j].getName())] = duj
		self.__dists[(forest[j].getName(), 'inner' + str(idx))] = duj
		forest = np.hstack((forest, [nodeU]))
		D = appendRowCol(D, u)
		D = withoutIndices(D, (i, j))
		forest = np.delete(forest, (i, j))
		return self.neighborJoin(forest, D, idx + 1)
	
	def getDist(self, na, nb):
		return self.__dists[(na.getName(), nb.getName())]

	def getAllDists(self):
		return self.__dists
	def drawTree(self):
		Phylo.draw_ascii(self.__repr__())
	
	def __repr__(self):
		xml = '<?xml version="1.0" encoding="UTF-8"?> ' \
			  ' <phyloxml xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://www.phyloxml.org http://www.phyloxml.org/1.10/phyloxml.xsd" ' \
			  'xmlns="http://www.phyloxml.org">' \
			  ' <phylogeny rooted="true"> ' \
			  ' <name>An example</name> <clade>'
		
		def makeXML(root):
			if root.isLeaf():
				return "<name>%s</name>\n" % root.getName()
			res = '<name>%s</name>\n' % root.getName()
			res += '<clade>\n'
			res += makeXML(root.getLchild())
			res += '</clade>\n'
			res += '<clade>\n'
			res += makeXML(root.getRchild())
			res += '</clade>\n'
			return res
		
		xml += makeXML(self.__root)
		xml += ' </clade> \n</phylogeny> \n</phyloxml>'
		
		with open('tmptree.xml', 'w') as f:
			f.write(xml)
		tree = Phylo.read('tmptree.xml', 'phyloxml')
		return tree  ## creates a tree
# tree = PhylTree('seqs.fasta')
#
# # prints the tree
# Phylo.draw_ascii(tree.__repr__())
