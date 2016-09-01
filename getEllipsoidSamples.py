import os
import select
import sys
import stat
from functools import partial
from optparse import OptionParser, OptionGroup

import numpy as np
import time

import lal
import lalsimulation as lalsim

import glue.lal
from glue.ligolw import utils, ligolw, lsctables, table, ilwd
lsctables.use_in(ligolw.LIGOLWContentHandler)
from glue.ligolw.utils import process
from glue import pipeline

#from pylal import series
from lal import series


from lalinference.rapid_pe import lalsimutils as lsu
import effectiveFisher as eff
from lalinference.rapid_pe import common_cl

'''
This function gives the samples from the 3D ambiguity ellipsoid around the triggered mass1,
mass2 and chi1 value.

Example of the function call:
samples = getSamples(1135654616, 10.0, 1.4, -0.5, 1000, {'H1=../psds_2016.xml.gz'}, {'L1=../psds_2016.xml.gz'}, saveData=True, plot=True, path=/path/where/you/want/the/output/to/go)

'''



def getSamples(graceid, mass1, mass2, chi1, samples, h_PSD, l_PSD, saveData=False, plot=False, path=False, show=False):
    m1_SI = mass1 * lal.MSUN_SI
    m2_SI = mass2 * lal.MSUN_SI
    min_mc_factor, max_mc_factor = 0.9, 1.1
    min_eta, max_eta = 0.05, 0.25
    min_chi1, max_chi1 = -0.99, 0.99
    # Control evaluation of the effective Fisher grid
    NMcs = 5
    NEtas = 5
    NChis = 5
    match_cntr = 0.9 # Fill an ellipsoid of match = 0.9
    wide_match = 1 - (1 - match_cntr)**(2/3.0)
    fit_cntr = match_cntr # Do the effective Fisher fit with pts above this match
    Nrandpts = samples # Requested number of pts to put inside the ellipsoid

    template_min_freq = 40.
    ip_min_freq = 40.
    lambda1, lambda2 = 0, 0

    #
    # Setup signal and IP class
    #
    param_names = ['Mc', 'eta', 'spin1z']
    McSIG = lsu.mchirp(m1_SI, m2_SI)
    etaSIG = lsu.symRatio(m1_SI, m2_SI)
    chiSIG = chi1

    PSIG = lsu.ChooseWaveformParams(
            m1=m1_SI, m2=m2_SI, spin1z=chi1,
            lambda1=lambda1, lambda2=lambda2,
            fmin=template_min_freq,
            approx=lalsim.GetApproximantFromString('SpinTaylorT4')
            )

    # Find a deltaF sufficient for entire range to be explored
    PTEST = PSIG.copy()

    # Check the waveform generated in the corners for the
    # longest possible waveform
    PTEST.m1, PTEST.m2 = lsu.m1m2(McSIG*min_mc_factor, min_eta)
    deltaF_1 = lsu.findDeltaF(PTEST)
    PTEST.m1, PTEST.m2 = lsu.m1m2(McSIG*min_mc_factor, max_eta)
    deltaF_2 = lsu.findDeltaF(PTEST)
    # set deltaF accordingly
    PSIG.deltaF = min(deltaF_1, deltaF_2)

    PTMPLT = PSIG.copy()

    psd_map = common_cl.parse_cl_key_value(h_PSD)
    for inst, psdfile in psd_map.items():
        if psd_map.has_key(psdfile):
            psd_map[psdfile].add(inst)
        else:
            psd_map[psdfile] = set([inst])
        del psd_map[inst]

    for psdf, insts in psd_map.iteritems():
        #xmldoc = utils.load_filename(psdf, contenthandler=series.LIGOLWContentHandler) ## CHECK: This is for old series in pylal
        xmldoc = utils.load_filename(psdf, contenthandler=series.PSDContentHandler)
        # FIXME: How to handle multiple PSDs
        for inst in insts:
#            psd = series.read_psd_xmldoc(xmldoc)[inst] ## CHECK: This is for old series in pylal
            psd = series.read_psd_xmldoc(xmldoc, root_name=None)[inst]
#            psd_f_high = len(psd.data)*psd.deltaF ## CHECK: This for old series in pylal
            psd_f_high = len(psd.data.data)*psd.deltaF
            f = np.arange(0, psd_f_high, psd.deltaF)
            fvals = np.arange(0, psd_f_high, PSIG.deltaF)

            def anon_interp(newf):
#                return np.interp(newf, f, psd.data) ## CHECK: This for old series in pylal
                return np.interp(newf, f, psd.data.data)
            eff_fisher_psd = np.array(map(anon_interp, fvals))

    analyticPSD_Q = False

    IP = lsu.Overlap(fLow = ip_min_freq,
        deltaF = PSIG.deltaF,
        psd = eff_fisher_psd,
        analyticPSD_Q = analyticPSD_Q
        )

    hfSIG = lsu.norm_hoff(PSIG, IP)

    # Find appropriate parameter ranges
    min_mc = McSIG * min_mc_factor
    max_mc = McSIG * max_mc_factor
    param_ranges = eff.find_effective_Fisher_region(PSIG, IP, wide_match,
            param_names, [[min_mc, max_mc],[min_eta, max_eta], [min_chi1, max_chi1]])


    # setup uniform parameter grid for effective Fisher
    pts_per_dim = [NMcs, NEtas, NChis]

    Mcpts, etapts, chipts = eff.make_regular_1d_grids(param_ranges, pts_per_dim)
    etapts = map(lsu.sanitize_eta, etapts)
    McMESH, etaMESH, chiMESH = eff.multi_dim_meshgrid(Mcpts, etapts, chipts)
    McFLAT, etaFLAT, chiFLAT = eff.multi_dim_flatgrid(Mcpts, etapts, chipts)
    dMcMESH = McMESH - McSIG
    detaMESH = etaMESH - etaSIG
    dchiMESH = chiMESH - chiSIG
    dMcFLAT = McFLAT - McSIG
    detaFLAT = etaFLAT - etaSIG
    dchiFLAT = chiFLAT - chiSIG

    grid = eff.multi_dim_grid(Mcpts, etapts, chipts)

    # Change units on Mc
    dMcFLAT_MSUN = dMcFLAT / lal.MSUN_SI
    dMcMESH_MSUN = dMcMESH / lal.MSUN_SI
    McMESH_MSUN = McMESH / lal.MSUN_SI
    McSIG_MSUN = McSIG / lal.MSUN_SI

    # Evaluate ambiguity function on the grid
    rhos = np.array(eff.evaluate_ip_on_grid(hfSIG, PTMPLT, IP, param_names, grid))
    rhogrid = rhos.reshape(NMcs, NEtas, NChis)

    # Fit to determine effective Fisher matrix
    # Adapt the match value to make sure all the Evals are positive

    evals = np.array([-1, -1, -1])
    count = 0
    start = time.time()
    match_cntrs = np.array([0.9, 0.99, 0.999])
    while np.any( np.array( [np.real(evals[0]), np.real(evals[1]), np.real(evals[2])] ) < 0 ):
        if count>0: print 'At least one of the eval is negative: switching to match of ' + str(match_cntrs[count])
        wide_match = 1 - (1 - match_cntrs[count])**(2/3.0)
        fit_cntr = match_cntrs[count] # Do the effective Fisher fit with pts above this match
        cut = rhos > fit_cntr
        if np.sum(cut) >= 6:
            fitgamma = eff.effectiveFisher(eff.residuals3d, rhos[cut], dMcFLAT_MSUN[cut], detaFLAT[cut], dchiFLAT[cut])
            # Find the eigenvalues/vectors of the effective Fisher matrix
            gam = eff.array_to_symmetric_matrix(fitgamma)
            evals, evecs, rot = eff.eigensystem(gam)
            count += 1
            if (count >= 3) and np.any( np.array( [np.real(evals[0]), np.real(evals[1]), np.real(evals[2])] ) < 0 ):
                return adapt_failure()
                sys.exit()
        else:
            return adapt_failure()
            sys.exit()

    #
    # Distribute points inside predicted ellipsoid of certain level of overlap
    #
    r1 = np.sqrt(2.*(1.-match_cntr)/np.real(evals[0])) # ellipse radii ...
    r2 = np.sqrt(2.*(1.-match_cntr)/np.real(evals[1])) # ... along eigen-directions
    r3 = np.sqrt(2.*(1.-match_cntr)/np.real(evals[2])) # ... along eigen-directions

    NN = 0
    NN_total = 0
    cart_grid = [[0., 0., 0.]]
    sph_grid = [[0., 0., 0.]]
    while NN < Nrandpts:
        NN_total += 1
        r = np.random.rand()
        ph = np.random.rand() * 2.*np.pi
        costh = np.random.rand()*2. - 1.
        sinth = np.sqrt(1. - costh * costh)
        th = np.arccos(costh)
        rrt = r**(1./3.)
        x1 = r1 * rrt * sinth * np.cos(ph)
        x2 = r2 * rrt * sinth * np.sin(ph)
        x3 = r3 * rrt * costh
        rand_Mc = x1 * lal.MSUN_SI + McSIG # Mc (kg)
        rand_eta = x2 + etaSIG # eta
        rand_chi = x3 + chiSIG
        condition1 = rand_eta > 0
        condition2 = rand_eta <= 0.25
        condition3 = np.abs(rand_chi) < 1.0
        joint_condition = condition1 * condition2 * condition3
        if joint_condition:
            cart_grid.append([x1, x2, x3])
            sph_grid.append([rrt, th, ph])
            NN += 1
    cart_grid = np.array(cart_grid)
    sph_grid = np.array(sph_grid)
    print 'Selected ' + str(NN) + ' points from ' + str(NN_total) + ' random samples within the ellipsoid'

    # Rotate to get coordinates in parameter basis
    cart_grid = np.array([ np.real( np.dot(rot, cart_grid[i]))
        for i in xrange(len(cart_grid)) ])
    # Put in convenient units,
    # change from parameter differential (i.e. dtheta)
    # to absolute parameter value (i.e. theta = theta_true + dtheta)
    rand_dMcs_MSUN, rand_detas, rand_dChis = tuple(np.transpose(cart_grid)) # dMc, deta, dchi
    rand_Mcs = rand_dMcs_MSUN * lal.MSUN_SI + McSIG # Mc (kg)
    rand_etas = rand_detas + etaSIG # eta
    rand_chis = rand_dChis + chiSIG

    # Prune points with unphysical values of eta from cart_grid
    rand_etas = np.array(map(partial(lsu.sanitize_eta, exception=np.NAN), rand_etas))
    cart_grid = np.transpose((rand_Mcs,rand_etas,rand_chis)) #### CHECK ####
    phys_cut = ~np.isnan(cart_grid).any(1) # cut to remove unphysical pts
    cart_grid = cart_grid[phys_cut]
    keep_phys_spins = np.abs(cart_grid[:,2]) < 1.0 #### CHECK ####
    cart_grid = cart_grid[keep_phys_spins] #### CHECK ####

    # Output Cartesian and spherical coordinates of intrinsic grid
    indices = np.arange(len(cart_grid))
    Mcs_MSUN, etas, chis = np.transpose(cart_grid) #### CHECK ####
    Mcs_MSUN = Mcs_MSUN / lal.MSUN_SI
    radii, thetas, phis = np.transpose(sph_grid[phys_cut][keep_phys_spins]) #### CHECK ####
    outgrid = np.transpose((Mcs_MSUN,etas,chis)) #### CHECK ####

    if saveData:
        if path: data_dir = path + '/intrinsic_grids'
        else: data_dir = 'intrinsic_grids'
        os.system('mkdir -p ' + data_dir)
        np.savetxt(data_dir + '/intrinsic_grid_' + str(graceid) + '_samples.dat', outgrid, fmt='%f\t%f\t%f')

    if plot:
        import pylab as pl
        from mpl_toolkits.mplot3d import Axes3D
        mc = outgrid[:,0]
        eta = outgrid[:,1]
        sz1 = outgrid[:,2]
        fig = pl.figure(figsize=(12, 10))
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(mc[1:], eta[1:], sz1[1:], zdir='z', s=3, c='b', alpha=0.2, depthshade=True)

        xlim = ax.get_xlim3d()
        ylim = ax.get_ylim3d()
        zlim = ax.get_zlim3d()
        print xlim
        print ylim
        print zlim

        ax.scatter(mc, eta, zlim[0]*np.ones(len(outgrid)), zdir='z', s=3, c='k', alpha=0.04, depthshade=True)
        ax.scatter(mc, ylim[1]*np.ones(len(outgrid)), sz1, zdir='z', s=3, c='k', alpha=0.04, depthshade=True)
        ax.scatter(xlim[0]*np.ones(len(outgrid)), eta, sz1, zdir='z', s=3, c='k', alpha=0.04, depthshade=True)
        ax.scatter(mc[0], eta[0], sz1[0], s=100, c='r', depthshade=True)

        ax.scatter(mc[0], eta[0], zlim[0], s=100, c='r', alpha=0.1, depthshade=True)
        ax.scatter(mc[0], ylim[1], sz1[0], s=100, c='r', alpha=0.1, depthshade=True)
        ax.scatter(xlim[0], eta[0], sz1[0], s=100, c='r', alpha=0.1, depthshade=True)
        ax.set_xlim3d(xlim)
        ax.set_ylim3d(ylim)
        ax.set_zlim3d(zlim)

        ax.set_xlabel('$\\mathcal{M}$')
        ax.set_ylabel('$\\eta$')
        ax.set_zlabel('$\\chi_1$')
        pl.title('$m_1$ = ' + str(mass1) + '$M_{\\odot}$: $m_2$ = ' + str(mass2) + '$M_{\\odot}$: $\\chi_1$ = ' + str(chi1), y=1.03)
        pl.grid()
        pl.rcParams.update({'font.size': 13})

        if path: plot_dir = path + '/ellipsoid_sample_plots'
        else: plot_dir = 'ellipsoid_sample_plots'
        os.system('mkdir -p ' + plot_dir)
        pl.savefig(plot_dir + '/ellipsoid_sample_' + str(graceid) + '_plot.png')
        if show:
            pl.show()

    return outgrid



def adapt_failure():
    print 'Could not find an all positive set Evals in three attempts... Quitting program'
    outgrid = np.array([np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan])
    np.savetxt('intrinsic_grid_Failed.dat', outgrid, newline="\t")
    return np.array([np.nan, np.nan, np.nan])




