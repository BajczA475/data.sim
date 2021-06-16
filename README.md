# data.sim
Repository for a package of functions for simulating random data sets

The functions in this package help you generate random data sets for various purposes, most notably for testing model performance for generalized linear (mixed) modeling within a frequentist or Bayestian framework. 
sim.pred.data() generates simulated X variable data given certain instructions and specifications.
sim.obs.data() generated simulated Y variable data given inputted X data and certain instructions and specifications. 
The typical workflow will involve using sim.pred.data to generate the data to be fed into sim.obs.data. 

See the manuals for the functions for more details as well as worked examples. 
