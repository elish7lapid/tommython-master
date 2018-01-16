import numpy as np
from estimate_branch_length import *

import PhylTree

ALPHABET_SIZE = 4
ALPHABET = ['A', 'T', 'C', 'G', '-']
ROOT_BASE_CASE = 'None'
STATIONARY = 0.2  # it is 0.2 but it doesnt affect the probability
NORMALIZE_FACTOR = 1


class TPU:
	def __init__(self, tree):
		self.__T = tree
		self.__root = tree.getRoot()
	
	def postOrder(self, root, stream_arr):
		"""
		:param root: of type Node
		:param stream_arr: an empty array to output the results to
		"""
		if root is None:
			return
		self.postOrder(root.getLchild(), stream_arr)
		self.postOrder(root.getRchild(), stream_arr)
		stream_arr.append(root)
	
	def preOrder(self, root, stream_arr):
		"""
		:param root: of type Node
		:param stream_arr: an empty array to output the results to
		"""
		if root is None:
			return
		stream_arr.append(root)
		self.preOrder(root.getLchild(), stream_arr)
		self.preOrder(root.getRchild(), stream_arr)
	
	def init_U(self, tree_list):
		"""
		initializes a dictionary with entry names corresponding to the
		nodes names in the ordered_tree.
		each entry points to another dictionary with entries corresponding
		to the ALPHABET of the sequences
		:param tree_list: an array containing all the nodes of the tree
		:return: the U dictionary
		"""
		U = {}
		for node in tree_list:
			U[node.getName()] = {}
		
		return U
	
	def init_U_edges(self, tree_list, root):
		"""
		initializes a dictionary with entry names corresponding to edges in
		a tree. a name of an edge between "node1" and "node2" is the
		concatenation of their names: "node1node2"
		each entry points to another dictionary with entries corresponding
		to the ALPHABET of the sequences
		:param root: the tree's root.
		:param tree_list: an array containing all the nodes of the tree
		:return: the U_edges dictionary
		"""
		U = {}
		for node in tree_list:
			if not node.isLeaf():
				name = node.getName()
				for other_node in [node.getLchild(), node.getRchild()]:
					other_name = other_node.getName()
					U[name + other_name] = {}
					U[other_name + name] = {}
		
		# a base case- initializing the edge from the root to a non-real
		# node that comes before it. the entry is: "rootNameNone"
		U[root.getName() + ROOT_BASE_CASE] = {x: 1 for x in ALPHABET}
		return U

	def init_trans_mat(self, tree_list):
		"""
		initializes a dictionary with entry names corresponding to edges in
		a tree. a name of an edge between "node1" and "node2" is the
		concatenation of their names: "node1node2"
		each entry points to another dictionary with entries corresponding
		to the ALPHABET of the sequences
		:param root: the tree's root.
		:param tree_list: an array containing all the nodes of the tree
		:return: the U_edges dictionary
		"""
		U = {}
		for node in tree_list:
			if not node.isLeaf():
				name = node.getName()
				for other_node in [node.getLchild(),node.getRchild()]:
					other_name = other_node.getName()
					U[name + other_name] = 0

		return U

	def singleLetterUp(self, post_tree, letter_idx, U_edges):
		"""
		calculates the Up matrix (up algorithm) for a single station
		in all of the sequences. uses the pseudo code learned in class
		:param post_tree: array of all tree nodes in post-order
		:param letter_idx: the index of the current station in the
		sequences
		:param U_edges: initialized U_edges matrix (see above
		documentation)
		:return: the U matrix of the up algorithm.
		"""
		
		U = self.init_U(post_tree)
		
		for node in post_tree:
			node_sequence = node.getData()
			i = node.getName()
			if node.isLeaf():
				for letter in ALPHABET:
					if node_sequence[letter_idx] == letter:
						U[i][letter] = 1
					else:
						U[i][letter] = 0
			
			else:
				for letter in ALPHABET:
					U[i][letter] = 1
					
					children = [node.getLchild(), node.getRchild()]
					
					for child in children:
						j = child.getName()
						U_edges[i + j][letter] = \
							np.sum([U[j][b] * probability(letter, b,   self.__T.getDist(node, child))for b in ALPHABET])
						U[i][letter] *= U_edges[i + j][letter]
		
		return U
	
	def singleLetterDown(self, pre_tree, U_edges):
		"""
		calculates the Down matrix (down algorithm) for a single station
		in all of the sequences. uses the pseudo code learned in class
		:param pre_tree: all tree nodes in pre-order
		:param U_edges:the U_edges matrix after going through Up algorithm
		 (see above documantation)
		"""
		for node in pre_tree:
			i = node.getName()
			children = [] if node.isLeaf() else [node.getLchild(), node.getRchild()]
			if node.getParent() is None:
				father_name = ROOT_BASE_CASE
			
			else:
				father_name = node.getParent().getName()
			
			for child_id in range(len(children)):
				other_child = 1 if child_id == 0 else 0
				other_name = children[other_child].getName()
				j = children[child_id].getName()
				for letter in ALPHABET:
					U_edges[j + i][letter] = \
						np.sum([probability(b, letter, self.__T.getDist(node, children[child_id])) *
						U_edges[i + other_name][b] * U_edges[i + father_name][b] for b in ALPHABET]) * NORMALIZE_FACTOR
		
		return
	
	def initial_tree_data(self, tree_list):
		for node in tree_list:
			if not node.isLeaf():
				node.setData('')

	def UpDown(self, root, post_tree, pre_tree):
		"""
		calculates up and down for all sequences and returns a tree with
		calculated sequences in all nodes
		:param root:
		:param post_tree: all tree nodes in an array of post-order
		:param pre_tree: "" "" "" "" of pre-order
		:return: the log likelihood of the leafs of the updates tree
		"""
		
		self.initial_tree_data(post_tree)
		r = root.getName()
		trans_mat = self.init_trans_mat(post_tree)
		
		log_likelihood = 0
		
		for letter_id in range(len(post_tree[0].getData())):
			U_edges = self.init_U_edges(post_tree, root)
			U = self.singleLetterUp(post_tree, letter_id, U_edges)
			self.singleLetterDown(pre_tree, U_edges)

			leafs_prob = np.sum([(U[r][a] * STATIONARY) for a in ALPHABET])
			
			for node in post_tree:
				if node.isLeaf():
					continue

				father = ROOT_BASE_CASE

				if node.getParent() is not None:
					father = node.getParent().getName()

				i = node.getName()
				left = node.getLchild()
				right = node.getRchild()

				for child in [left, right]:
					j = child.getName()
					other_name = left.getName() if child.getName() == right.getName() else right.getName()
					trans_mat[i + j] += (np.sum([U[j][a]*
					                            probability(a,a,self.__T.getDist(node, child))*
					                            U_edges[i+other_name][a]*U_edges[i+father][a] for a in ALPHABET])/(5*leafs_prob))
				
				letters_prob = []
				
				if node.getParent() is None:
					father_name = ROOT_BASE_CASE
				else:
					father_name = node.getParent().getName()

				for letter in ALPHABET:
					calc = U[i][letter] * STATIONARY * NORMALIZE_FACTOR*U_edges[i + father_name][letter]
					letters_prob.append(calc / leafs_prob)

				chosen_letter = np.argmax(letters_prob)
				update_data = node.getData() + ALPHABET[chosen_letter]
				node.setData(update_data)

			log_likelihood += np.log(leafs_prob)

		return log_likelihood, trans_mat
	
	def updateBranchLengths(self, post_tree, trans_mat):
		"""
		updates branch lengths using the JC
		:param tree_list: a list of all tree nodes
		"""
		length_seq = len(post_tree[0].getData())
		for node in post_tree:
			i = node.getName()
			if not node.isLeaf():
				for child in [node.getRchild(), node.getLchild()]:
					child_name = child.getName()
					N1 = trans_mat[i + child_name]
					length = estimate_time(N1, length_seq - N1)
					self.__T.updateDist(node, child, length)

	def update_tree_with_EM(self, converge=0.001):
		"""
		uses an EM algorithm to update a tree's branches lengths and
		inner node's sequences
		:param root:
		:param converge: threshold for convergence of the EM algorithm
		"""

		root = self.__root
		post_tree = []
		pre_tree = []
		self.postOrder(root, post_tree)
		self.preOrder(root, pre_tree)

		prev_likely = -1
		cur_likely = 0
		self.__T.setAllEdgeDist(0.5)
		cnt = 0
		while np.abs(cur_likely - prev_likely) > converge:

			if cnt % 2 == 0:
				print("iter: %d, likely: %f" % (cnt, cur_likely))
			prev_likely = cur_likely

			cur_likely, trans_mat = self.UpDown(root, post_tree, pre_tree)
			self.updateBranchLengths(post_tree, trans_mat)
			cnt += 1

		print("iter: %d, likely: %f" % (cnt, cur_likely))
