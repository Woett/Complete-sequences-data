import math
import numpy as np
import matplotlib.pyplot as plt
from typing import List, Tuple, Dict, Optional
from fractions import Fraction

# --- CONFIGURATION ---
t_min_initial = Fraction(10, 10)
t_max_initial = Fraction(30, 10)
alpha_min_initial = Fraction(13, 10)
alpha_max_initial = Fraction(14, 10)
max_iterations = 20
N_terms = 20
plot_interval = 1
plot_aspect_ratio = 1.0
fixed_denom = 100*(2**max_iterations)
# --- END CONFIGURATION ---


# --- Helper functions ---
def format_with_denominator(frac: Fraction, fixed_denominator: int) -> str:
    """
    Returns a string representation of a fraction with a specified denominator.
    """
    if frac.denominator == fixed_denominator:
        return f"{frac.numerator}/{fixed_denominator}"
    
    # Calculate the scaling factor
    if fixed_denominator % frac.denominator != 0:
        # Handle cases where the denominator can't be exactly converted
        # You might want to raise an error or round
        print(f"Warning: Cannot exactly convert {frac} to denominator {fixed_denominator}")
        
    scaling_factor = fixed_denominator // frac.denominator
    new_numerator = frac.numerator * scaling_factor
    
    return f"{new_numerator}/{fixed_denominator}"

def lemma_check_py(svec: List[int]) -> Tuple:
    """
    Return (True, r, X, prefix) if lemma condition holds for some r (1..N-1),
    using prefix s1..s_{r+1}. Otherwise return (False,).
    """
    N = len(svec)
    for r in range(1, N):
        S = sum(svec[:r])
        if svec[r] > S + 1:
            continue

        bs = np.zeros(S + 1, dtype=bool)
        bs[0] = True
        for si in svec[:r]:
            for k in range(S - si, -1, -1):
                if bs[k]:
                    bs[k + si] = True

        block = svec[r]  # s_{r+1}
        limitX = S - block + 1

        if limitX >= 1:
            for X in range(1, limitX + 1):
                ok = True
                for m in range(X, X + block):
                    if m > S or (not bs[m]):
                        ok = False
                        break
                if ok:
                    return (True, r, X, svec[:r + 1])
    return (False,)

def gen_svec(alpha: Fraction, t: Fraction, N: int) -> List[int]:
    """Generates the sequence s_n = floor(t * alpha^n)."""
    return [int(math.floor(t * (alpha ** n))) for n in range(1, N + 1)]

# --- Main simulation logic ---

class Rectangle:
    def __init__(self, t_min: Fraction, alpha_min: Fraction, t_max: Fraction, alpha_max: Fraction):
        self.t_min = t_min
        self.alpha_min = alpha_min
        self.t_max = t_max
        self.alpha_max = alpha_max

def get_overlap_with_dummy(svec1: List[int], svec2: List[int]) -> Tuple[List[int], Optional[Tuple[int,int,int]]]:
    """
    Compute overlap, inserting dummy x at the first mismatch.
    Return (sequence_with_x, dummy_info)
    where dummy_info = (pos, low, high).
    """
    overlap: List[int] = []
    min_len = min(len(svec1), len(svec2))
    mismatch_found = False
    dummy_info = None

    for i in range(min_len):
        if not mismatch_found:
            if svec1[i] == svec2[i]:
                overlap.append(svec1[i])
            else:
                # first mismatch → insert dummy
                low, high = min(svec1[i], svec2[i]), max(svec1[i], svec2[i])
                overlap.append("x")  # placeholder
                dummy_info = (i, low, high)
                mismatch_found = True
        else:
            # after mismatch: keep if difference ≤ 1
            if abs(svec1[i] - svec2[i]) > 1:
                break
            if svec1[i] == svec2[i]:
                overlap.append(svec1[i])

    # if no mismatch at all, sequences identical
    if not mismatch_found:
        last_index = min_len - 1
        val = svec1[last_index]
        overlap = svec1[:]  # whole thing
        dummy_info = (last_index, val, val)

    return overlap, dummy_info

def plot_rectangles(colored_rects: List[Rectangle], uncolored_rects: List[Rectangle],
                    iteration: int, t_min_limit: Fraction, t_max_limit: Fraction,
                    alpha_min_limit: Fraction, alpha_max_limit: Fraction, plot_aspect_ratio: Fraction):
    """Plot rectangles: green = colored, red = to be subdivided."""
    fig, ax = plt.subplots(figsize=(8, 8 * float(plot_aspect_ratio)))

    for rect in colored_rects:
        rect_patch = plt.Rectangle(
            (float(rect.t_min), float(rect.alpha_min)),
            float(rect.t_max - rect.t_min),
            float(rect.alpha_max - rect.alpha_min),
            facecolor='green',
            edgecolor='black',
            linewidth=0.5
        )
        ax.add_patch(rect_patch)

    for rect in uncolored_rects:
        rect_patch = plt.Rectangle(
            (float(rect.t_min), float(rect.alpha_min)),
            float(rect.t_max - rect.t_min),
            float(rect.alpha_max - rect.alpha_min),
            fill=False,
            edgecolor='red',
            linewidth=1.0,
            linestyle='--'
        )
        ax.add_patch(rect_patch)

    plt.ylim(float(alpha_min_limit), float(alpha_max_limit))
    plt.xlim(float(t_min_limit), float(t_max_limit))

    ax.set_title(f'Rectangle Subdivision: Iteration {iteration}')
    ax.set_xlabel('t')
    ax.set_ylabel('alpha')
    plt.grid(True)
    plt.show()

def main():
    initial_rect = Rectangle(t_min_initial, alpha_min_initial, t_max_initial, alpha_max_initial)
    uncolored_rects = [initial_rect]
    colored_rects: List[Rectangle] = []

    # Caches
    svec_cache: Dict[Tuple[Fraction, Fraction], List[int]] = {}
    lemma_cache: Dict[Tuple[int, ...], Tuple] = {}
    prefix_good: Dict[Tuple[int, ...], bool] = {}
    prefix_range_good: Dict[Tuple[Tuple[int, ...], int, int], bool] = {}

    # Stats
    stats = {
        "lemma_checks_total": 0,
        "lemma_checks_done": 0,
        "lemma_prefix_hits": 0,
        "lemma_range_hits": 0,
    }

    def lemma_check_cached(svec: List[int]) -> Tuple:
        key = tuple(svec)
        stats["lemma_checks_total"] += 1
        if key in lemma_cache:
            return lemma_cache[key]
        res = lemma_check_py(svec)
        lemma_cache[key] = res
        stats["lemma_checks_done"] += 1
        if res[0]:
            prefix_good[key] = True
        return res

    for iteration in range(max_iterations+1):
        rects_to_process_this_iter = list(uncolored_rects)
        uncolored_rects = []

        for rect in rects_to_process_this_iter:
            t_min, alpha_min = rect.t_min, rect.alpha_min
            t_max, alpha_max = rect.t_max, rect.alpha_max

            if (alpha_min, t_min) not in svec_cache:
                svec_cache[(alpha_min, t_min)] = gen_svec(alpha_min, t_min, N_terms)
            svec_bl = svec_cache[(alpha_min, t_min)]

            if (alpha_max, t_max) not in svec_cache:
                svec_cache[(alpha_max, t_max)] = gen_svec(alpha_max, t_max, N_terms)
            svec_tr = svec_cache[(alpha_max, t_max)]

            overlap_svec, dummy_info = get_overlap_with_dummy(svec_bl, svec_tr)

            lemma_success = False
            if len(overlap_svec) > 3:
                # Case 0: already known good prefix
                prefix_tuple = tuple(x for x in overlap_svec if isinstance(x,int))
                if prefix_tuple in prefix_good:
                    lemma_success = True
                    stats["lemma_prefix_hits"] += 1
                else:
                    # Case 1: plain overlap (replace "x" with nothing)
                    pure_overlap = [x for x in overlap_svec if isinstance(x,int)]
                    if lemma_check_cached(pure_overlap)[0]:
                        prefix_good[prefix_tuple] = True
                        lemma_success = True
                    # Case 2: try variable dummy
                    elif dummy_info is not None:
                        pos, low, high = dummy_info
                        range_key = (prefix_tuple, low, high)
                        if range_key in prefix_range_good:
                            lemma_success = True
                            stats["lemma_range_hits"] += 1
                        else:
                            all_good = True
                            for x in range(low, high + 1):
                                candidate = [x if v=="x" else v for v in overlap_svec]
                                if not lemma_check_cached(candidate)[0]:
                                    all_good = False
                                    break
                            if all_good:
                                prefix_range_good[range_key] = True
                                lemma_success = True

            if lemma_success:
                print(f" ({format_with_denominator(t_min,fixed_denom)},{format_with_denominator(t_max, fixed_denom)},{format_with_denominator(alpha_min,fixed_denom)},{format_with_denominator(alpha_max,fixed_denom)}),")

                #print(f"\nCoordinates for the bottom left point of the rectangle: ({format_with_denominator(t_min,fixed_denom)},{format_with_denominator(alpha_min,fixed_denom)})")
                #print(f"Coordinates for the top right point of the rectangle: ({format_with_denominator(t_max, fixed_denom)},{format_with_denominator(alpha_max,fixed_denom)})")
                #print(f"Sequence for the bottom left point of the rectangle = {svec_bl}")
                #print(f"Sequence for the top right point of the rectangle = {svec_tr}")
                #print(f"Overlap between sequences with one variable = {overlap_svec}, with x between {dummy_info[1]} and {dummy_info[2]}")
                colored_rects.append(rect)
            else:
                t_mid = (t_min + t_max) / Fraction(2)
                alpha_mid = (alpha_min + alpha_max) / Fraction(2)
                uncolored_rects.append(Rectangle(t_min, alpha_mid, t_mid, alpha_max))
                uncolored_rects.append(Rectangle(t_mid, alpha_mid, t_max, alpha_max))
                uncolored_rects.append(Rectangle(t_min, alpha_min, t_mid, alpha_mid))
                uncolored_rects.append(Rectangle(t_mid, alpha_min, t_max, alpha_mid))

        #if (iteration + 1) % plot_interval == 0 or (iteration + 1) == max_iterations or not uncolored_rects:
            #plot_rectangles(colored_rects, uncolored_rects, iteration + 1,
                            #t_min_initial, t_max_initial,
                            #alpha_min_initial, alpha_max_initial,
                            #plot_aspect_ratio)

        if not uncolored_rects:
            print(f"\nAll rectangles have been colored. Calculation finished after {iteration} iterations.")
            break

    if uncolored_rects:
        print(f"\nNot all rectangles have been colored. Remaining uncolored rectangles after {iteration} iterations: {len(uncolored_rects)}")
    print(f"\nTotal colored rectangles: {len(colored_rects)}")
    #print("\n--- Stats ---")
    #print(f"Total lemma checks requested: {stats['lemma_checks_total']}")
    #print(f"Actual lemma_check_py runs:   {stats['lemma_checks_done']}")
    #print(f"Prefix cache hits:            {stats['lemma_prefix_hits']}")
    #print(f"Range cache hits:             {stats['lemma_range_hits']}")

if __name__ == "__main__":
    main()