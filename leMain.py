import sys
from TreeParametersUpdater import TPU
from PhylTree import PhylTree

if __name__ == '__main__':
	fastas = sys.argv[1]
	tree = PhylTree(fastas)
	tree.drawTree()
	updater = TPU(tree)
	updater.update_tree_with_EM()
