from PhylTree import Node
from PhylTree import PhylTree
import estimate_branch_length
import matplotlib.pyplot as plt


def pumping(node, tree, i, alpha):
        if node.getRchild() == None:
            return 1
        else:
            return estimate_branch_length.probability(node[i], node.getRchild[i], tree.getDist(node, node.getRchild)*alpha)* \
                   estimate_branch_length.probability(node[i], node.getLchild[i], tree.getDist(node, node.getLchild)*alpha) * \
                   pumping(node.getRchild, tree, i, alpha) * pumping(node.getLchild, tree, i, alpha)


def makePumping(node, tree):
    alphas = [0.1, 0.3, 0.5, 1, 2, 5, 20, 50, 100]
    vals = []
    for i in range(0, len(node.getData()), 10): # todo: to think - maybe to do it for every position instead of windows of 10
        the_max = 0
        alp = 0
        for a in alphas:
            val = 0
            for j in range(i, i+10):
                val = val*pumping(node, tree, j, a)
            if the_max < val :
                the_max = val
                alp = a
        vals.append(alp)

    plt.scatter(vals, alphas)
    plt.show() #todo to think of how we want thr graph
    return vals
