# Monte Carlo Estimation of π (Pi)

It demonstrates how to estimate the value of π using a "Monte Carlo simulation". Random points are generated inside a square, and the fraction that falls inside the unit circle is used to approximate π.

Algorithm:
1. Generate random (x, y) points uniformly in [-1,1] × [-1,1].
2. Count the proportion of points inside the circle (x^2 + y^2 <= 1).
3. Estimate π as:
   pi_hat = 4 *(points inside/total points)
4. Use a while-loop that continues sampling until the estimate reaches the desired number of significant figures(using a confidence interval half-width criterion).

How to use:
-Using MATLAB and run: for example, if requiring 3 significant figures, edit "pi_est = user_pi_precision(3);", or multiple targets (2, 3, and 4 s.f.), edit "pi_est = user_pi_precision([2 3 4]);". The other parts in the function are optional to edit and will be set to default if the user skips them.
