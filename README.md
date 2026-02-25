In this file you can find data and code for my paper 'Completeness of exponentially increasing sequences'. So for more information it makes most sense to just read that paper, but I will try to quickly provide some context on these files here anyway.

Every line in the five .txt-files contains a fourtuple of rationals $(t_1, t_2, \alpha_1, \alpha_2)$, which represent $x$ and $y$ coordinates of the bottom left and the top right corner of a rectangle in the plane. The variable that I call $t$ runs along the $x$-axis and the variable $\alpha$ runs along the $y$-axis. For example, the $x$ and $y$ coordinates from the first file titled 'Rectangles for t between 1 and 3 and alpha between 1.3 and 1.4', partition the region of length $3-1 = 2$ and height $1.4 - 1.3 = 0.1$ into $27373$ little rectangles.

Now, why are we subdividing the plane into tiny rectangles?

Well, given any point $(t, \alpha)$ in the plane, we are interested in whether the corresponding sequence

$S_t(\alpha) = (\lfloor t \alpha^1 \rfloor, \lfloor t \alpha^2 \rfloor, \ldots)$

is complete. That is, whether every large enough integer can be written as the sum of distinct terms of this sequence $S_t(\alpha)$. This is the content of Erdős Problem #349 (https://www.erdosproblems.com/349).

Lemma $5$ in my paper shows how a finite calculation can suffice in order to conclude that a sequence $S_t(\alpha)$ is complete. In fact, as explained in the proof of Proposition $8$, if these finite calculations are essentially the same for the pairs $(t_1, \alpha_1)$ and $(t_2, \alpha_2)$, then one can further conclude the completeness of all sequences $S_t(\alpha)$ with $t_1 \le t \le t_2$ and $\alpha_1 \le \alpha \le \alpha_2$. That is, partitioning a region into rectangles and then performing the necessary computations for the bottom left and the top right corners of all rectangles, potentially proves the completeness for all pairs $(t, \alpha)$ in the region. And this GitHub file provides exactly such partitions. 

I have also added the Python code I used to find these partitions, where one can add in the currently commented out print statements under 'if lemma_success' in order to get more information on the sequences $S_t(\alpha)$ corresponding to the bottom left and the top right corners of the rectangles. For a visual view of what's happening when the region is divided into smaller and smaller rectangles until it works, one can also add the currently commented out line 'plot_rectangles(colored_rects, uncolored_rects, iteration + 1, #t_min_initial, t_max_initial, #alpha_min_initial, alpha_max_initial, #plot_aspect_ratio)' back into the code.
