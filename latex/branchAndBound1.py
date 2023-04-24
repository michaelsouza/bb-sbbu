import numpy as np

class Solution:
    def __init__(self):        
        self.val = np.inf
        self.sol = None

    
def evaluate(perm):
  w = [2, 5, 7, 5, 3, 1, 4, 9, 6, 8]
  cost = 0
  for i in sorted(range(len(perm)), reverse=True):
    # remainder of division by len(w) is the index of the weight
    j = perm[i] % len(w)
    cost += 10**w[j]
  return cost


def swap(perm, i, j):
  perm[i], perm[j] = perm[j], perm[i]


def calculateLowerBound(perm, level):
  return evaluate(perm[:(level + 1)])


# Recursive search function to construct and evaluate permutations
def search(perm, level, xopt:Solution):
    n = len(perm)
    # If the permutation is complete, evaluate and update bestSolution and bestValue
    if level == n:
        currentValue = evaluate(perm)
        if currentValue < xopt.val:
            print('New best solution: {} with value: {}'.format(perm, currentValue))
            xopt.val = currentValue
            xopt.sol[:] = perm
    # Otherwise, construct and search through different permutations
    else:
        for i in range(level, n):
            # Swap elements to construct a new permutation
            swap(perm, level, i)

            # Calculate lower bound for pruning
            lowerBound = calculateLowerBound(perm, level)

            # If lower bound is better than bestValue, continue the search recursively
            if lowerBound < xsol.val:
                search(perm, level + 1, xopt)
            else:
                print('Pruned permutation: {} with lower bound: {}'.format(perm[:level], lowerBound))

            # Swap elements back to restore the original permutation
            swap(perm, level, i)


def branch_and_bound(n, xsol):
    # Create the initial permutation and start the search from level 0
    permutation = [i for i in range(n)]
    search(permutation, 0, xsol)    


if __name__ == "__main__":
    n = 7 # Number of elements in the permutation
    # Create an initial solution
    xsol = Solution()
    xsol.sol = [i for i in range(n)]
    xsol.val = evaluate(xsol.sol)
    print("Initial solution: {}".format(xsol.sol))
    print("Initial value: {}".format(xsol.val))

    branch_and_bound(n, xsol)
    print("Best solution: {}".format(xsol.sol))
    print("Best value: {}".format(xsol.val))
