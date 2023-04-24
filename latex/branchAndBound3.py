from typing import List
import random


SEGMENTS = [ (),
    set([5,6]),                 # s_1
    set([7,8,9,10]),            # s_2
    set([11,12,13,14,15,16]),   # s_3
    set([17]),                  # s_4
    set([18]),                  # s_5
    set([19]),                  # s_6
    set([20]),                  # s_7
]

EDGES = [[],
    [1,2],          # e_1
    [2,3,4],        # e_2
    [3,4,5,6,7],    # e_3
    [4,5,6],        # e_4
    [7]             # e_5
    ]

def cost_function(permutation: List[int]) -> int:
    cost = 1
    S = set() # set all covered segments
    for i in permutation:
        edge = EDGES[i]
        for j in edge:
            if j in S:
                continue
            S = S | j
            s = SEGMENTS[j]
            cost += 2**len(s)
    return 0 if cost == 1 else cost

def relaxed_cost_function(partial_permutation: List[int], remaining: List[int]) -> int:    
    # Calculate the cost of the current partial permutation
    partial_cost = cost_function(partial_permutation)

    # Add the minimum cost for the remaining numbers
    S = set() # set all covered segments
    for i in partial_permutation:
        edge = EDGES[i]
        for j in edge:
            if j in S:
                continue
            S = S | j

    for j in range(1, len(SEGMENTS)):
        if j in S:
            continue
        partial_cost += 2**len(SEGMENTS[j])
    return partial_cost

def branch_and_bound(n: int) -> List[int]:
    def dfs(partial_permutation: List[int], remaining: List[int], relaxed_cost:int, depth: int):
        nonlocal best_permutation, min_cost

        print(f"{'  ' * depth}Current partial_permutation: {partial_permutation}")
        print(f"{'  ' * depth}Current remaining: {remaining}")
        print(f"{'  ' * depth}Current relaxed_cost: {relaxed_cost}")


        if not remaining:
            cost = cost_function(partial_permutation)
            if cost < min_cost:
                min_cost = cost
                best_permutation = partial_permutation.copy()
                print(f"{'  ' * depth}New best_permutation: {best_permutation}")
                print(f"{'  ' * depth}New min_cost: {min_cost}")                
            return

        for i, edge in enumerate(remaining):
            new_partial_permutation = partial_permutation + [edge]
            relaxed_cost = relaxed_cost_function(new_partial_permutation, remaining)

            if relaxed_cost < min_cost:
                remaining.pop(i)
                dfs(new_partial_permutation, remaining, relaxed_cost, depth + 1)
                remaining.insert(i, edge)
            else:
                print(f"{'  ' * depth}Pruning branch with relaxed_cost: {relaxed_cost}")
                print(f"{'  ' * depth}Current new_partial_permutation: {new_partial_permutation}")

    best_permutation = None
    min_cost = float('inf')    
    edges = [1,2,3,4,5]    

    dfs([], edges, 0, 0)
    return best_permutation

if __name__ == "__main__":
    n = 4
    minimal_cost_permutation = branch_and_bound(n)
    print("Minimal cost permutation:", minimal_cost_permutation)
    print("Cost:", cost_function(minimal_cost_permutation))
