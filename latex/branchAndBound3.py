from typing import List
import random

def cost_function(permutation: List[int]) -> int:
    return sum(abs(permutation[i] - permutation[i + 1]) for i in range(len(permutation) - 1))

def relaxed_cost_function(partial_permutation: List[int], remaining: List[int]) -> int:    
    # Calculate the cost of the current partial permutation
    partial_cost = cost_function(partial_permutation)

    # Add the minimum cost for the remaining numbers
    partial_cost += (len(remaining) - 1) * 1

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

        for i, num in enumerate(remaining):
            new_partial_permutation = partial_permutation + [num]
            relaxed_cost = relaxed_cost_function(new_partial_permutation, remaining)

            if relaxed_cost < min_cost:
                remaining.pop(i)
                dfs(new_partial_permutation, remaining, relaxed_cost, depth + 1)
                remaining.insert(i, num)
            else:
                print(f"{'  ' * depth}Pruning branch with relaxed_cost: {relaxed_cost}")
                print(f"{'  ' * depth}Current new_partial_permutation: {new_partial_permutation}")

    best_permutation = None
    min_cost = float('inf')
    # set the seed to get the same random numbers every time
    random.seed(0)
    # create a list of five random numbers
    numbers = [random.randint(1, 100) for _ in range(5)]    

    dfs([], numbers, 0, 0)
    return best_permutation

if __name__ == "__main__":
    n = 4
    minimal_cost_permutation = branch_and_bound(n)
    print("Minimal cost permutation:", minimal_cost_permutation)
    print("Cost:", cost_function(minimal_cost_permutation))
