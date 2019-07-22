#!/usr/binenv python 
#
# Tutorial on how to use dadi
# this file is an edited version of the example file supplied with the dadi distribution
# 
# More explanations will be added, soon, 
#           but you will need a copy of the manual
#           and the README (see tutorial folder 
#           to follow
# 2015 Peter Beerli
import numpy 
import dadi

# In demographic_models.py, we've defined a custom model for this problem
#
import demographic_models

# Load the original data
data = dadi.Spectrum.from_file('YRI_CEU.fs')
ns = data.sample_sizes
print 'DADI analysis ================================================'
print 'Raw data: file YRY_CEU.fs'
print 'samples sizes:', ns
print 
S = data.S()
print "Segregating sites:",S
Fst = data.Fst()
print "FST:",Fst
print "Marginalized Spectrum populations aka combining YRI + CEU"
margdata = data.marginalize([1])
theta_w =  margdata.Watterson_theta()
print "Watterson's theta:", theta_w
# These are the grid point settings will use for extrapolation.
# 
pts_l = [40,50,60]
#pts_l = [80,90,100]

# The Demographics1D and Demographics2D modules contain a few simple models,
# mostly as examples. We could use one of those.
#func = dadi.Demographics2D.split_mig
# Instead, we'll work with our custom model
func = demographic_models.prior_onegrow_mig
paramlist = ["AncSize after growth","bottleneck Popsize2","Popsize2","migration m","time for ancestor to grow","Population split"]


# Now let's optimize parameters for this model.

# The upper_bound and lower_bound lists are for use in optimization.
# Occasionally the optimizer will try wacky parameter values. We in particular
# want to exclude values with very long times, very small population sizes, or
# very high migration rates, as they will take a long time to evaluate.

# For the onegrow_mig model:
# Parameters are: (nu1F, nu2B, nu2F, m, Tp, T)
#                  popsize
upper_bound = [100, 100, 100, 10, 3, 3]
lower_bound = [1e-2, 1e-2, 1e-2, 0, 0, 0]

# This is our initial guess for the parameters, which is somewhat arbitrary.
p0 = [2,0.1,2,1,0.2,0.2]
# Make the extrapolating version of our demographic model function.
func_ex = dadi.Numerics.make_extrap_log_func(func)

# Perturb our parameters before optimization. This does so by taking each
# parameter a up to a factor of two up or down.
p0 = dadi.Misc.perturb_params(p0, fold=1, upper_bound=upper_bound,
                              lower_bound=lower_bound)
# Do the optimization. By default we assume that theta is a free parameter,
# since it's trivial to find given the other parameters. If you want to fix
# theta, add a multinom=False to the call.
# The maxiter argument restricts how long the optimizer will run. For real 
# runs, you will want to set this value higher (at least 10), to encourage
# better convergence. You will also want to run optimization several times
# using multiple sets of intial parameters, to be confident you've actually
# found the true maximum likelihood parameters.
print('Beginning optimization ************************************************')
popt = dadi.Inference.optimize_log(p0, data, func_ex, pts_l, 
                                   lower_bound=lower_bound,
                                   upper_bound=upper_bound,
                                   verbose=len(p0), maxiter=5)
# The verbose argument controls how often progress of the optimizer should be
# printed. It's useful to keep track of optimization process.
print('Finshed optimization **************************************************')

# These are the actual best-fit model parameters, which we found through
# longer optimizations and confirmed by running multiple optimizations.
# We'll work with them through the rest of this script.
popt_lonrun = [1.881, 0.0710, 1.845, 0.911, 0.355, 0.111]
print('Best-fit parameters: {0}'.format(popt))
print('Best-fit parameters: {0}'.format(popt_lonrun),'(precooked values)')
print "Table"
print "-------------------------------------------------------------------------"
print '{0:25s}  {1:15s}  {2:8s}  {3:8s}  {4:10s}'.format('Parameter','  Estimate','Best','Lower','Upper')
for i in range(len(paramlist)):
    print '{0:25s}  {1:10.5f}  {2:10.5f}  {3:10.5f}  {4:10.5f}'.format(paramlist[i],popt[i],popt_lonrun[i],lower_bound[i],upper_bound[i])
print "-------------------------------------------------------------------------"
# Calculate the best-fit model AFS.
model = func_ex(popt, ns, pts_l)
# Likelihood of the data given the model AFS.
ll_model = dadi.Inference.ll_multinom(model, data)
print('Maximum log composite likelihood: {0}'.format(ll_model))
# The optimal value of theta given the model.
theta = dadi.Inference.optimal_sfs_scaling(model, data)
print('Optimal value of theta: {0}'.format(theta))

# Plot a comparison of the resulting fs with the data.
import pylab
pylab.figure(1)
dadi.Plotting.plot_2d_comp_multinom(model, data, vmin=1, resid_range=3,pop_ids =('YRI','CEU'))
# This ensures that the figure pops up. It may be unecessary if you are using
# ipython.
#pylab.show()
# Save the figure
#pylab.savefig('YRI_CEU_test.png', dpi=50)

# Let's generate some data using ms, if you have it installed.
mscore = demographic_models.prior_onegrow_mig_mscore(popt)
# I find that it's most efficient to simulate with theta=1, average over many
# iterations, and then scale up.
#print ns
#print mscore
mscommand = dadi.Misc.ms_command(1., ns, mscore, int(1e5))
# We use Python's os module to call this command from within the script.
import os
scaledtheta=1.0
datasize = ns
reps = 100000
rec = 0
chrL=10*scaledtheta
mscommand = dadi.Misc.ms_command(scaledtheta, datasize, mscore, 1, recomb=rec, rsites=chrL)
print "Example bootstrap sample"
msfs=dadi.Spectrum.from_ms_file(os.popen(mscommand))
pylab.figure(2)
pylab.title('Example bootstrap sample')
dadi.Plotting.plot_2d_comp_multinom(model, theta*msfs, vmin=1,
                                        pop_ids=('YRI','CEU'))
#    pylab.show()
#pylab.savefig('YRI_CEU_boot.png', dpi=50)

for ii in range(100):
    msfs=dadi.Spectrum.from_ms_file(os.popen(mscommand))
    ## fold simulated data because original data is also folded
    #    msfold = msfs.fold()
    msfs.to_file('bootstraps/{0:02d}.fs'.format(ii)) 

# Estimate parameter uncertainties using the Godambe Information Matrix, to
# account for linkage in the data. To use the GIM approach, we need to have
# spectra from bootstrapping our data.  Let's load the ones we've provided for
# the example.  
# (We're using Python list comprehension syntax to do this in one line.)
all_boot = [dadi.Spectrum.from_file('bootstraps/{0:02d}.fs'.format(ii)) 
            for ii in range(100)]
uncerts = dadi.Godambe.GIM_uncert(func_ex, pts_l, all_boot, popt, data, 
                                  multinom=True)
# uncert contains the estimated standard deviations of each parameter, with
# theta as the final entry in the list.
print('Estimated parameter standard deviations from GIM: {0}'.format(uncerts))

# For comparison, we can estimate uncertainties with the Fisher Information
# Matrix, which doesn't account for linkage in the data and thus underestimates
# uncertainty. (Although it's a fine approach if you think your data is truly
# unlinked.)
uncerts_fim = dadi.Godambe.FIM_uncert(func_ex, pts_l, popt, data, multinom=True)
print('Estimated parameter standard deviations from FIM: {0}'.format(uncerts_fim))

print('Factors by which FIM underestimates parameter uncertainties: {0}'.format(uncerts/uncerts_fim))

# What if we fold the data?
# These are the optimal parameters when the spectrum is folded. They can be
# found simply by passing data.fold() to the above call to optimize_log. 
##popt_fold =  array([1.907,  0.073,  1.830,  0.899,  0.425,  0.113])
##uncerts_folded = dadi.Godambe.GIM_uncert(func_ex, pts_l, all_boot, popt_fold, 
#                                         data.fold(), multinom=True)
##print('Folding increases parameter uncertainties by factors of: {0}'.format(uncerts_folded/uncerts))

# Let's do a likelihood-ratio test comparing models with and without migration.
# The no migration model is implemented as 
# demographic_models.prior_onegrow_nomig
func_nomig = demographic_models.prior_onegrow_nomig
func_ex_nomig = dadi.Numerics.make_extrap_log_func(func_nomig)
upper_bound_nomig = [100, 100, 100, 3, 3]
lower_bound_nomig = [1e-2, 1e-2, 1e-2, 0, 0]
# This is our initial guess for the parameters, which is somewhat arbitrary.
p0_nomig = [2,0.1,2,0.2,0.2]
p0_nomig = dadi.Misc.perturb_params(p0_nomig, fold=1, upper_bound=upper_bound_nomig,
                              lower_bound=lower_bound_nomig)
print('Beginning optimization ************************************************')
popt_nomig = dadi.Inference.optimize_log(p0_nomig, data, func_ex_nomig, pts_l, 
                                   lower_bound=lower_bound_nomig,
                                   upper_bound=upper_bound_nomig,
                                   verbose=len(p0_nomig), maxiter=3)
# The verbose argument controls how often progress of the optimizer should be
# printed. It's useful to keep track of optimization process.
print('Finished optimization **************************************************')
nparamlist = ["AncSize after growth","bottleneck Popsize2","Popsize2","time for ancestor to grow","Population split"]
popt2 = popt.copy()
popt2 = numpy.delete(popt2,3)
print "Table"
print "-------------------------------------------------------------------------"
print '{0:25s}  {1:15s}  {2:8s}  {3:8s}  {4:10s}'.format('Parameter','  Estimate','Full','Lower','Upper')
for i in range(len(nparamlist)):
    print '{0:25s}  {1:10.5f}  {2:10.5f}  {3:10.5f}  {4:10.5f}'.format(nparamlist[i],popt_nomig[i],popt2[i],lower_bound_nomig[i],upper_bound_nomig[i])
print "-------------------------------------------------------------------------"

# These are the best-fit parameters, which we found by multiple optimizations
#popt_nomig = array([ 1.897,  0.0388,  9.677,  0.395,  0.070])
model_nomig = func_ex_nomig(popt_nomig, ns, pts_l)
ll_nomig = dadi.Inference.ll_multinom(model_nomig, data)

# Since LRT evaluates the complex model using the best-fit parameters from the
# simple model, we need to create list of parameters for the complex model
# using the simple (no-mig) best-fit params.  Since evalution is done with more
# complex model, need to insert zero migration value at corresponding migration
# parameter index in complex model. And we need to tell the LRT adjust function
# that the 3rd parameter (counting from 0) is the nested one.
#p_lrt = [1.897,  0.0388,  9.677, 0, 0.395,  0.070]
p_lrt = popt_nomig
print p_lrt
p_lrt = numpy.insert(p_lrt,3,0)
print p_lrt

adj = dadi.Godambe.LRT_adjust(func_ex, pts_l, all_boot, p_lrt, data, 
                              nested_indices=[3], multinom=True)
D_adj = adj*2*(ll_model - ll_nomig)
print('Adjusted D statistic: {0:.4f}'.format(D_adj))

# Because this is test of a parameter on the boundary of parameter space 
# (m cannot be less than zero), our null distribution is an even proportion 
# of chi^2 distributions with 0 and 1 d.o.f. To evaluate the p-value, we use the
# point percent function for a weighted sum of chi^2 dists.
pval = dadi.Godambe.sum_chi2_ppf(D_adj, weights=(0.5,0.5))
print('p-value for rejecting no-migration model: {0:.4f}'.format(pval))
