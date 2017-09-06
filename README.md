# Voorn-Overbeek
Construct phase diagram for polyelectrolyte mixtures using Voorn-Overbeek theory [1]

For help: python main.py --help

Optional arguments: 
  --N N; polymer length (assume N1 = N2) (default N=1000)
  --sigma sigma; polymer charge fraction (assume sigma1 = sigma2) (default sigma = 0.44)
  --alpha alpha; reduced temperature (see eqn. (4) of reference [1]) (default alpha = 3.655)
  --outfile filename.dat; name of output file to write data
  
  additionally one can include a Flory-Huggins interaction parameter acting between polymer types (chi12) using
  --chi12 chi12; Flory Huggins paramter between polymer1 and polymer2 (default chi12 = 0.0)
  
  
  References:
  [1] J.T.G. Overbeek and M.J. Voorn, "Phase Separation in Polyelectrolyte Solutions. Theory of Complex Coacervation", Journal of Cellular Physiology, Vol 49, Issue S1, Pg. 7-26 (1957).
  
