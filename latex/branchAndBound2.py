from typing import List

def cost_function(permutation: List[int]) -> int:
    return sum(abs(permutation[i] - permutation[i + 1]) for i in range(len(permutation) - 1))

def branch_and_bound(n: int) -> List[int]:
    def dfs(partial_permutation: List[int], remaining: List[int], partial_cost: int, depth: int):
        nonlocal best_permutation, min_cost

        print(f"{'  ' * depth}Current partial_permutation: {partial_permutation}")
        print(f"{'  ' * depth}Current remaining: {remaining}")
        print(f"{'  ' * depth}Current partial_cost: {partial_cost}")

        if not remaining:
            if partial_cost < min_cost:
                min_cost = partial_cost
                best_permutation = partial_permutation.copy()
                print(f"{'  ' * depth}New best_permutation: {best_permutation}")
                print(f"{'  ' * depth}New min_cost: {min_cost}")
            return

        for i, num in enumerate(remaining):
            new_partial_cost = partial_cost + (abs(num - partial_permutation[-1]) if partial_permutation else 0)
            if new_partial_cost < min_cost:
                remaining.pop(i)
                partial_permutation.append(num)
                dfs(partial_permutation, remaining, new_partial_cost, depth + 1)
                partial_permutation.pop()
                remaining.insert(i, num)

    best_permutation = None
    min_cost = float('inf')
    numbers = list(range(1, n + 1))

    dfs([], numbers, 0, 0)
    return best_permutation

if __name__ == "__main__":
    n = 4
    minimal_cost_permutation = branch_and_bound(n)
    print("\nFinal minimal cost permutation:", minimal_cost_permutation)
    print("Final cost:", cost_function(minimal_cost_permutation))
