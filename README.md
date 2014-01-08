Mathematica Markov Chain Monte Carlo
====================================
Mathematica package containing a [Markov chain Monte Carlo] [1] routine I wrote for the purpose of fitting models to data. Contains various examples and documentation. MCMC allows the posterior probability distribution of a model's parameters to be stochastically sampled based on the parameter likelihood, thus maximizing the sampling efficiency. Features:
 * Assuming Gaussian data errors, yields statistically rigorous parameter errors (using chi square statistic to evaluate parameter likelihood)
 * Handles both real-valued and discrete-valued model parameters; allows for vector-valued independent and dependent variables
 * Uses Metropolis algorithm with decaying exponential proposal distribution
 * Progress monitor; support for auto save/resume

  [1]: http://en.wikipedia.org/wiki/Markov_chain_Monte_Carlo        "Markov chain Monte Carlo"
