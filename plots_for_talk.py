"""
Run like:
> python plots_for_talk.py -r 3545-3549

Sorry, Ann.
"""
import argparse
import array
import copy
import math
import os
import sys
import ROOT
import numpy as np 
ROOT.gROOT.SetBatch()

bigpdf = False

def main():

    ops = options()
    if not ops.r:
        fatal("Please provide run number with -r")
    if not ops.r in ["3513", "3518", "3519", "3520", "3518-3520", "3522", "3523", "3526", "3527", "3528", "3530", "3539",
                     "3545", "3546", "3547", "3548", "3545-3547", "3545-3548", "3549", "3545-3549"]:
        fatal("Script is not configured with this run (%s)" % (ops.r))

    dotrigger = ops.r in ["3522", "3523", "3539"]
    dopacman  = ops.r in ["3522"]
    doonlyres = ops.r in ["3513"]

    rfile = ROOT.TFile.Open(filename(ops.r, ops.a))
    pdf   = "plots_for_talk_%s%s.pdf" % (ops.r, "" if not ops.a else "_"+ops.a)

    pacman          = rfile.Get("histograms/strip_position_vs_board")
    q_vs_ch         = [rfile.Get("histograms/strip_q_vs_ch_%s" % (bo)) for bo in xrange(8)]
    pdo_vs_ch       = [rfile.Get("histograms/strip_pdo_vs_ch_%s" % (bo)) for bo in xrange(8)]
    tdo_vs_ch       = [rfile.Get("histograms/strip_tdo_vs_ch_%s" % (bo)) for bo in xrange(8)]
    tdoc_vs_ch      = [rfile.Get("histograms/strip_tdoc_vs_ch_%s" % (bo)) for bo in xrange(8)]
    zdri_vs_ch      = [rfile.Get("histograms/strip_zdri_vs_ch_%s" % (bo)) for bo in xrange(8)]
    zpos_vs_ch      = [rfile.Get("histograms/strip_zpos_vs_ch_%s" % (bo)) for bo in xrange(8)]
    ch_vs_time      = [rfile.Get("histograms/strip_event_vs_ch_%s"  % (bo)) for bo in xrange(8)]
    dups_vs_ch      = rfile.Get("histograms/dups_vs_ch")
    clus_per_board  = rfile.Get("histograms/clus_vs_board")
    strips_per_clus = rfile.Get("histograms/hits_per_clus_vs_board")
    delta_bcid      = rfile.Get("histograms/timediff_vs_board")
    track_angle_6   = rfile.Get("histograms/track_angle_6")
    track_angle_7   = rfile.Get("histograms/track_angle_7")
    track_angle_8   = rfile.Get("histograms/track_angle_8")
    track_angle_N   = rfile.Get("histograms/track_angle_N_denom")
    track_angle_N2  = rfile.Get("histograms/track_angle_N_numer")
    residuals       = rfile.Get("histograms/track_N1_board_vs_residual")
    residuals_norm  = rfile.Get("histograms/track_N1_board_vs_residual_norm")
    residuals_tpc   = rfile.Get("histograms/track_N1_board_vs_utpc")
    residuals_art   = rfile.Get("histograms/track_N1_board_vs_residual_art")
    residuals10     = rfile.Get("histograms/track_N1_board_vs_residual_10deg")
    res_vs_theta_0  = rfile.Get("histograms/track_N1_theta_x_vs_residual_0")
    res_vs_theta_1  = rfile.Get("histograms/track_N1_theta_x_vs_residual_1")
    res_vs_theta_2  = rfile.Get("histograms/track_N1_theta_x_vs_residual_2")
    res_vs_theta_3  = rfile.Get("histograms/track_N1_theta_x_vs_residual_3")
    res_vs_theta_4  = rfile.Get("histograms/track_N1_theta_x_vs_residual_4")
    res_vs_theta_5  = rfile.Get("histograms/track_N1_theta_x_vs_residual_5")
    res_vs_theta_6  = rfile.Get("histograms/track_N1_theta_x_vs_residual_6")
    res_vs_theta_7  = rfile.Get("histograms/track_N1_theta_x_vs_residual_7")
    res_vs_x_0      = rfile.Get("histograms/track_N1_x_vs_residual_0")
    res_vs_x_1      = rfile.Get("histograms/track_N1_x_vs_residual_1")
    res_vs_x_2      = rfile.Get("histograms/track_N1_x_vs_residual_2")
    res_vs_x_3      = rfile.Get("histograms/track_N1_x_vs_residual_3")
    res_vs_x_4      = rfile.Get("histograms/track_N1_x_vs_residual_4")
    res_vs_x_5      = rfile.Get("histograms/track_N1_x_vs_residual_5")
    res_vs_x_6      = rfile.Get("histograms/track_N1_x_vs_residual_6")
    res_vs_x_7      = rfile.Get("histograms/track_N1_x_vs_residual_7")
    res_vs_x_10d_0  = rfile.Get("histograms/track_N1_x_vs_residual_10deg_0")
    res_vs_x_10d_1  = rfile.Get("histograms/track_N1_x_vs_residual_10deg_1")
    res_vs_x_10d_2  = rfile.Get("histograms/track_N1_x_vs_residual_10deg_2")
    res_vs_x_10d_3  = rfile.Get("histograms/track_N1_x_vs_residual_10deg_3")
    res_vs_x_10d_4  = rfile.Get("histograms/track_N1_x_vs_residual_10deg_4")
    res_vs_x_10d_5  = rfile.Get("histograms/track_N1_x_vs_residual_10deg_5")
    res_vs_x_10d_6  = rfile.Get("histograms/track_N1_x_vs_residual_10deg_6")
    res_vs_x_10d_7  = rfile.Get("histograms/track_N1_x_vs_residual_10deg_7")
    tpc_vs_theta_0  = rfile.Get("histograms/track_N1_theta_x_vs_utpc_0")
    tpc_vs_theta_1  = rfile.Get("histograms/track_N1_theta_x_vs_utpc_1")
    tpc_vs_theta_2  = rfile.Get("histograms/track_N1_theta_x_vs_utpc_2")
    tpc_vs_theta_3  = rfile.Get("histograms/track_N1_theta_x_vs_utpc_3")
    tpc_vs_theta_4  = rfile.Get("histograms/track_N1_theta_x_vs_utpc_4")
    tpc_vs_theta_5  = rfile.Get("histograms/track_N1_theta_x_vs_utpc_5")
    tpc_vs_theta_6  = rfile.Get("histograms/track_N1_theta_x_vs_utpc_6")
    tpc_vs_theta_7  = rfile.Get("histograms/track_N1_theta_x_vs_utpc_7")
    xy_ss           = rfile.Get("histograms/track_x_vs_y_0123_vs_4567")
    xy_sd           = rfile.Get("histograms/track_x_vs_y_0236_vs_1457")
    xy_dd           = rfile.Get("histograms/track_x_vs_y_0256_vs_1347")
    mult_vs_theta   = rfile.Get("histograms/track_clusmult_vs_theta")
    clus_on_track   = rfile.Get("histograms/track_hits_vs_time_fid")
    eff_numer       = rfile.Get("histograms/track_obs_hit")
    eff_denom       = rfile.Get("histograms/track_exp_hit")
    scint_0         = rfile.Get("histograms/track_scint_0_vs_time")
    scint_1         = rfile.Get("histograms/track_scint_1_vs_time")
    scint_2         = rfile.Get("histograms/track_scint_2_vs_time")
    scint_3         = rfile.Get("histograms/track_scint_3_vs_time")
    scint_4         = rfile.Get("histograms/track_scint_4_vs_time")
    scint_5         = rfile.Get("histograms/track_scint_5_vs_time")

    track_unc_0     = rfile.Get("histograms/track_unc_0")
    track_unc_1     = rfile.Get("histograms/track_unc_1")
    track_unc_6     = rfile.Get("histograms/track_unc_6")
    track_unc_7     = rfile.Get("histograms/track_unc_7")

    zres_vs_ch_0    = [rfile.Get("histograms/strip_zres_vs_ch_BC0_%s"     % (bo)) for bo in xrange(8)]
    zres_vs_ch_1    = [rfile.Get("histograms/strip_zres_vs_ch_BC1_%s"     % (bo)) for bo in xrange(8)]
    zres_vs_ch_2    = [rfile.Get("histograms/strip_zres_vs_ch_BC2_%s"     % (bo)) for bo in xrange(8)]
    zres_vs_ch_3    = [rfile.Get("histograms/strip_zres_vs_ch_BC3_%s"     % (bo)) for bo in xrange(8)]
    zres_vs_ch      = [rfile.Get("histograms/strip_zres_vs_ch_%s"         % (bo)) for bo in xrange(8)]
    zpos_vs_ztrack  = [rfile.Get("histograms/strip_zpos_vs_ztrack_%s"     % (bo)) for bo in xrange(8)]
    tpc_phi_vs_phi  = [rfile.Get("histograms/tpc_phi_vs_phi_%s"           % (bo)) for bo in xrange(8)]
    tpc_vs_bary     = [rfile.Get("histograms/track_N1_theta_x_vs_comp_%s" % (bo)) for bo in xrange(8)]
    dx_01_bary      = rfile.Get("histograms/track_diff01_bary_theta")
    dx_67_bary      = rfile.Get("histograms/track_diff67_bary_theta")
    dx_01_btpc      = rfile.Get("histograms/track_diff01_btpc_theta")
    dx_67_btpc      = rfile.Get("histograms/track_diff67_btpc_theta")
    dx_01_utpc      = rfile.Get("histograms/track_diff01_utpc_theta")
    dx_67_utpc      = rfile.Get("histograms/track_diff67_utpc_theta")
    dx_01_comb      = rfile.Get("histograms/track_diff01_comb_theta")
    dx_67_comb      = rfile.Get("histograms/track_diff67_comb_theta")
    tpc_fitprobs    = rfile.Get("histograms/track_N1_board_vs_prob")
    tpc_unc_0       = rfile.Get("histograms/tpc_uncertainty_0")
    tpc_unc_1       = rfile.Get("histograms/tpc_uncertainty_1")
    tpc_unc_6       = rfile.Get("histograms/tpc_uncertainty_6")
    tpc_unc_7       = rfile.Get("histograms/tpc_uncertainty_7")
    res01_barybary  = rfile.Get("histograms/track_diff01_btpc_norm")
    res67_barybary  = rfile.Get("histograms/track_diff67_btpc_norm")
    res06_barybary  = rfile.Get("histograms/track_diff06_btpc_norm")
    res17_barybary  = rfile.Get("histograms/track_diff17_btpc_norm")
    res01_utpcutpc  = rfile.Get("histograms/track_diff01_utpc_norm")
    res67_utpcutpc  = rfile.Get("histograms/track_diff67_utpc_norm")
    res06_utpcutpc  = rfile.Get("histograms/track_diff06_utpc_norm")
    res17_utpcutpc  = rfile.Get("histograms/track_diff17_utpc_norm")
    res00_utpcbary  = rfile.Get("histograms/track_N1_theta_x_vs_comp_01")
    res11_utpcbary  = rfile.Get("histograms/track_N1_theta_x_vs_comp_11")
    res66_utpcbary  = rfile.Get("histograms/track_N1_theta_x_vs_comp_61")
    res77_utpcbary  = rfile.Get("histograms/track_N1_theta_x_vs_comp_71")

    art_th      = rfile.Get("histograms/trig_theta")
    artn        = rfile.Get("histograms/trig_art")
    clun        = rfile.Get("histograms/trig_mm")
    pairs       = rfile.Get("histograms/trig_dbc_pairs")
    bcwin       = rfile.Get("histograms/trig_dbc_vs_N")
    dtheta      = rfile.Get("histograms/trig_dtheta_vsNX")
    dthetaok    = rfile.Get("histograms/trig_dtheta_vsNX_ok")
    bcwin_theta = rfile.Get("histograms/trig_dbc_vs_theta")
    trig_art_vs_event = rfile.Get("histograms/trig_art_vs_evt")

    dx_vs_theta_0 = rfile.Get("histograms/trig_dx_vs_theta_0")
    dx_vs_theta_1 = rfile.Get("histograms/trig_dx_vs_theta_1")
    dx_vs_theta_6 = rfile.Get("histograms/trig_dx_vs_theta_6")
    dx_vs_theta_7 = rfile.Get("histograms/trig_dx_vs_theta_7")
    dx_vs_theta_N = copy.copy(dx_vs_theta_0)
    if dotrigger:
        dx_vs_theta_N.Add(dx_vs_theta_1)
        dx_vs_theta_N.Add(dx_vs_theta_6)
        dx_vs_theta_N.Add(dx_vs_theta_7)
        dx_vs_theta_N.SetName(dx_vs_theta_N.GetName()+"N")

    dtheta_vs_theta_all = rfile.Get("histograms/trig_dtheta_vs_theta_all")
    dtheta_vs_theta_n24 = rfile.Get("histograms/trig_dtheta_vs_theta_near24")
    dtheta_vs_theta_n12 = rfile.Get("histograms/trig_dtheta_vs_theta_near12")

    artpos_by_hit  = rfile.Get("histograms/trig_artpos_hit_vs_clussize")
    artpos_by_bcid = rfile.Get("histograms/trig_artpos_bci_vs_clussize")
    artpos_by_pdo  = rfile.Get("histograms/trig_artpos_pdo_vs_clussize")

    trig_xres = rfile.Get("histograms/trig_vs_mmfe_x_vs_theta")

    # draw stuff
    if bigpdf:
        open_pdf(pdf)

    # time machine to 2016
    if doonlyres:
        residuals_N1(residuals,   0, 1, pdf, tag=None)
        residuals_N1(residuals,   6, 7, pdf, tag=None)
        if bigpdf:
            close_pdf(pdf)
        return

    # MMFE
    if dopacman:
        pacman_example(pacman, 6897, 2, pdf)
    #for bo in xrange(8):
    #    strip_pdo_vs_channel(pdo_vs_ch[bo], bo, pdf)
    for bo in xrange(8):
        strip_charge_vs_channel(q_vs_ch[bo], bo, pdf)
    for bo in xrange(8):
        strip_tdo_vs_channel(tdo_vs_ch[bo], bo, pdf)
    for bo in xrange(8):
        strip_tdoc_vs_channel(tdoc_vs_ch[bo], bo, pdf)
    for bo in xrange(8):
        if bo in [0, 1, 6, 7]:
            strip_zpos_vs_channel(zpos_vs_ch[bo], bo, pdf)
    for bo in xrange(8):
        if bo in [0, 1, 6, 7]:
            strip_zpos_vs_channel(zdri_vs_ch[bo], bo, pdf, ontrack=True)
    #for bo in xrange(8):
    #    strip_duplicates_vs_channel(dups_vs_ch, bo, pdf)
    #for bo in xrange(8):
    #    strip_lifetime(ch_vs_time[bo], bo, event_range(ops.r), pdf)
            
    #close_pdf(pdf)
    #sys.exit(1)

    for bo in xrange(8):
        if bo in [0, 1, 6, 7]:
            strip_zres(zres_vs_ch[bo],   bo, pdf)
            strip_zres(zres_vs_ch_0[bo], bo, pdf, "0")
            strip_zres(zres_vs_ch_1[bo], bo, pdf, "1")
            strip_zres(zres_vs_ch_2[bo], bo, pdf, "2")
            strip_zres(zres_vs_ch_3[bo], bo, pdf, "3")
    for bo in xrange(8):
        if bo in [0, 1, 6, 7]:
            strip_zpos_vs_ztrack(zpos_vs_ztrack[bo], bo, pdf)
    for bo in xrange(8):
        if bo in [0, 1, 6, 7]:
            phi_vs_phi(tpc_phi_vs_phi[bo], bo, pdf)
    tpc_fitprob(tpc_fitprobs,  pdf)
    tpc_unc(tpc_unc_0, 0, pdf)
    tpc_unc(tpc_unc_1, 1, pdf)
    tpc_unc(tpc_unc_6, 6, pdf)
    tpc_unc(tpc_unc_7, 7, pdf)

    track_unc(track_unc_0, 0, pdf)
    track_unc(track_unc_1, 1, pdf)
    track_unc(track_unc_6, 6, pdf)
    track_unc(track_unc_7, 7, pdf)

    #for bo in xrange(8):
    #    if bo in [0, 1, 6, 7]:
    #        residuals_vs_theta(tpc_vs_bary[bo], pdf, tpc=False, comp=True)

    residuals_1D(res01_barybary, "barybary01", pdf)
    residuals_1D(res67_barybary, "barybary67", pdf)
    residuals_1D(res06_barybary, "barybary06", pdf)
    residuals_1D(res17_barybary, "barybary17", pdf)
    residuals_1D(res01_utpcutpc, "utpcutpc01", pdf, rms=True)
    residuals_1D(res67_utpcutpc, "utpcutpc67", pdf, rms=True)
    residuals_1D(res06_utpcutpc, "utpcutpc06", pdf, rms=True)
    residuals_1D(res17_utpcutpc, "utpcutpc17", pdf, rms=True)
    residuals_1D(res00_utpcbary, "utpcbary00", pdf, rms=True)
    residuals_1D(res11_utpcbary, "utpcbary11", pdf, rms=True)
    residuals_1D(res66_utpcbary, "utpcbary66", pdf, rms=True)
    residuals_1D(res77_utpcbary, "utpcbary77", pdf, rms=True)

    dxrange = 25
    dx_ij_vs_theta(dx_01_bary, False, pdf, dxrange)
    dx_ij_vs_theta(dx_67_bary, False, pdf, dxrange)
    dx_ij_vs_theta(dx_01_btpc, False, pdf, dxrange, tag="#theta > 10 deg.")
    dx_ij_vs_theta(dx_67_btpc, False, pdf, dxrange, tag="#theta > 10 deg.")
    dx_ij_vs_theta(dx_01_utpc, False, pdf, dxrange, tag="#theta > 10 deg.")
    dx_ij_vs_theta(dx_67_utpc, False, pdf, dxrange, tag="#theta > 10 deg.")
    dx_ij_vs_theta(dx_01_comb, False, pdf, dxrange, tag="#theta > 10 deg.")
    dx_ij_vs_theta(dx_67_comb, False, pdf, dxrange, tag="#theta > 10 deg.")

    dxrange = 5
    dx_ij_vs_theta(dx_01_bary, False, pdf, dxrange)
    dx_ij_vs_theta(dx_67_bary, False, pdf, dxrange)
    dx_ij_vs_theta(dx_01_btpc, False, pdf, dxrange, tag="#theta > 10 deg.")
    dx_ij_vs_theta(dx_67_btpc, False, pdf, dxrange, tag="#theta > 10 deg.")
    dx_ij_vs_theta(dx_01_utpc, False, pdf, dxrange, tag="#theta > 10 deg.")
    dx_ij_vs_theta(dx_67_utpc, False, pdf, dxrange, tag="#theta > 10 deg.")
    dx_ij_vs_theta(dx_01_comb, False, pdf, dxrange, tag="#theta > 10 deg.")
    dx_ij_vs_theta(dx_67_comb, False, pdf, dxrange, tag="#theta > 10 deg.")

    dxrange = 2
    dx_ij_vs_theta(dx_01_bary, False, pdf, dxrange)
    dx_ij_vs_theta(dx_67_bary, False, pdf, dxrange)
    dx_ij_vs_theta(dx_01_btpc, False, pdf, dxrange, tag="#theta > 10 deg.")
    dx_ij_vs_theta(dx_67_btpc, False, pdf, dxrange, tag="#theta > 10 deg.")
    dx_ij_vs_theta(dx_01_utpc, False, pdf, dxrange, tag="#theta > 10 deg.")
    dx_ij_vs_theta(dx_67_utpc, False, pdf, dxrange, tag="#theta > 10 deg.")
    dx_ij_vs_theta(dx_01_comb, False, pdf, dxrange, tag="#theta > 10 deg.")
    dx_ij_vs_theta(dx_67_comb, False, pdf, dxrange, tag="#theta > 10 deg.")

    clusters_per_board(clus_per_board, pdf)
    strips_per_cluster(strips_per_clus, pdf)
    strips_per_cluster_vs_theta(mult_vs_theta, pdf)
    dbcid_mmfe_scint(delta_bcid, pdf)
    dbcid_mmfe_scint(delta_bcid, pdf, "0")
    dbcid_mmfe_scint(delta_bcid, pdf, "7")
    clusters_on_track(clus_on_track, pdf)
    track_angle(track_angle_7, 7, pdf)
    track_angle(track_angle_8, 8, pdf)
    # track_angle(track_angle_N, 0, pdf)

    if ops.r != "3539":
        residuals_N1(residuals,      0, 1, pdf, art=False, tag=None, tpc=False)
        residuals_N1(residuals,      6, 7, pdf, art=False, tag=None, tpc=False)
        residuals_N1(residuals_norm, 0, 1, pdf, art=False, tag=None, tpc=False, norm=True)
        residuals_N1(residuals_norm, 6, 7, pdf, art=False, tag=None, tpc=False, norm=True)
        #residuals_N1(residuals10, 0, 1, pdf, art=False, tag="10deg")
        #residuals_N1(residuals10, 6, 7, pdf, art=False, tag="10deg")
        residuals_vs_theta(res_vs_theta_0, pdf)
        residuals_vs_theta(res_vs_theta_1, pdf)
        residuals_vs_theta(res_vs_theta_2, pdf)
        residuals_vs_theta(res_vs_theta_3, pdf)
        residuals_vs_theta(res_vs_theta_4, pdf)
        residuals_vs_theta(res_vs_theta_5, pdf)
        residuals_vs_theta(res_vs_theta_6, pdf)
        residuals_vs_theta(res_vs_theta_7, pdf)
        residuals_N1(residuals_tpc, 0, 1, pdf, art=False, tag=">10deg", tpc=True)
        residuals_N1(residuals_tpc, 6, 7, pdf, art=False, tag=">10deg", tpc=True)
        residuals_vs_theta(tpc_vs_theta_0, pdf, tpc=True)
        residuals_vs_theta(tpc_vs_theta_1, pdf, tpc=True)
        #residuals_vs_theta(tpc_vs_theta_2, pdf, tpc=True)
        #residuals_vs_theta(tpc_vs_theta_3, pdf, tpc=True)
        #residuals_vs_theta(tpc_vs_theta_4, pdf, tpc=True)
        #residuals_vs_theta(tpc_vs_theta_5, pdf, tpc=True)
        residuals_vs_theta(tpc_vs_theta_6, pdf, tpc=True)
        residuals_vs_theta(tpc_vs_theta_7, pdf, tpc=True)
        #residuals_vs_x(res_vs_x_0, pdf)
        #residuals_vs_x(res_vs_x_1, pdf)
        #residuals_vs_x(res_vs_x_2, pdf)
        #residuals_vs_x(res_vs_x_3, pdf)
        #residuals_vs_x(res_vs_x_4, pdf)
        #residuals_vs_x(res_vs_x_5, pdf)
        #residuals_vs_x(res_vs_x_6, pdf)
        #residuals_vs_x(res_vs_x_7, pdf)
        #residuals_vs_x(res_vs_x_10d_0, pdf, angled=True)
        #residuals_vs_x(res_vs_x_10d_1, pdf, angled=True)
        #residuals_vs_x(res_vs_x_10d_2, pdf, angled=True)
        #residuals_vs_x(res_vs_x_10d_3, pdf, angled=True)
        #residuals_vs_x(res_vs_x_10d_4, pdf, angled=True)
        #residuals_vs_x(res_vs_x_10d_5, pdf, angled=True)
        #residuals_vs_x(res_vs_x_10d_6, pdf, angled=True)
        #residuals_vs_x(res_vs_x_10d_7, pdf, angled=True)
        # xy_difference(xy_ss, "sameX_sameUV", pdf)
        # xy_difference(xy_sd, "sameX_diffUV", pdf)
        # xy_difference(xy_dd, "diffX_diffUV", pdf)

    scintillator_slopes(rfile, pdf)
    for (isc, scint) in enumerate([scint_0, scint_1, scint_2,
                                   scint_3, scint_4, scint_5,
                                   ]):
        scintillator_slope_vs_time(isc, scint, pdf)
    vmm_efficiency(eff_numer, eff_denom, pdf, ops.r)

    # ART
    if dotrigger:

        if bigpdf:
            divide_pdf(pdf)

        art_n(artn, None, pdf)
        art_n(artn, clun, pdf)
        trigger_art_vs_event(trig_art_vs_event, pdf)
        art_theta(art_th, pdf)
        
        art_window(bcwin, pdf)
        art_window_vs_theta(bcwin_theta, pdf)
        bcwindow_vs_N(bcwin, pdf)
        art_pairs(pairs, pdf)
        
        art_dx_vs_theta(dx_vs_theta_0, pdf)
        art_dx_vs_theta(dx_vs_theta_1, pdf)
        art_dx_vs_theta(dx_vs_theta_6, pdf)
        art_dx_vs_theta(dx_vs_theta_7, pdf)
        art_dx_vs_theta(dx_vs_theta_N, pdf)
        
        theta_resolution_all(dtheta, "3VMM", pdf)
        theta_resolution_all(dtheta_vs_theta_n24.ProjectionY(), "pm24", pdf)
        art_dtheta_vs_theta(dtheta_vs_theta_all, "3VMM", pdf)
        art_dtheta_vs_theta(dtheta_vs_theta_n24, "pm24", pdf)
        if int(ops.r) < 3535:
            art_dtheta_rms(rfile, pdf, ops.r)
            art_dtheta_efficiency(rfile, 15, pdf, ops.r)
        #art_dtheta(rfile, pdf, road="08", nx="3X")
        #art_dtheta(rfile, pdf, road="08", nx="4X")
        
        x_resolution(trig_xres, pdf)

        residuals_N1(residuals_art, 0, 1, pdf, art=True, tag=None)
        residuals_N1(residuals_art, 6, 7, pdf, art=True, tag=None)

        if int(ops.r) != 3539:
            art_position_within_cluster(artpos_by_hit,  pdf)
            art_position_within_cluster(artpos_by_bcid, pdf)
            art_position_within_cluster(artpos_by_pdo,  pdf)
        else:
            trigger_efficiency_vs_angle(track_angle_N2.ProjectionX(), track_angle_N.ProjectionX(), pdf)
        
    if bigpdf:
        close_pdf(pdf)

def filename(run, align=""):
    if   run == "3518-3520" : return "test_3518_3519_3520.root"
    elif run == "3518"      : return "test_3518.root"
    elif run == "3519"      : return "test_3519.root"
    elif run == "3520"      : return "test_3520.root"
    elif run == "3522"      : return "test.root"
    elif run == "3523"      : return "test_3523.root"
    elif run == "3526"      : return "test_3526.root"
    elif run == "3527"      : return "test_3527.root"
    elif run == "3528"      : return "test_3528.root"
    elif run == "3530"      : return "test_3530.root"
    elif run == "3513"      : return "test_3513.root"
    elif run == "3539"      : return "test_3539.root"
    elif run == "3545"      : return "test_3545.root"
    elif run == "3546"      : return "test_3546.root"
    elif run == "3547"      : return "test_3547.root"
    elif run == "3548"      : return "test_3548.root"
    elif run == "3549"      : return "test_3549.root"
    elif run == "3545-3547" : return "test_3545-3547.root" if not align else "test_3545-3547_%s.root" % (align)
    elif run == "3545-3548" : return "test_3545-3548.root" if not align else "test_3545-3548_%s.root" % (align)
    # elif run == "3545-3549" : return "test_3545-3549.root" if not align else "test_3545-3549_%s.root" % (align)
    #elif run == "3545-3549" : return "tpc_Run3545-3549_no_offset.root"
    elif run == "3545-3549" : return "tpc_Run3545-3549_BC0p5_fixedoffsets.root"
    #elif run == "3545-3549" : return "tpc_Run3545-3549_vetosuspBC.root"
    #elif run == "3545-3549" : return "tpc_Run3545-3549_BC0p5.root"
    #elif run == "3545-3549" : return "tpc_Run3545-3549_updated_trackuncfix.root"
    #elif run == "3545-3549" : return "tpc_Run3545-3549.root"
    #elif run == "3545-3549" : return "../../utpc/octuplet_analysis/tpc_Run3545-3549_TDOonly.root"
    # elif run == "3545-3549" : return "test_3545-3549.root"
    fatal("Dont have filename for %s" % (run))

def event_range(run):
    if   run == "3518-3520" : return 600
    elif run == "3518"      : return 250
    elif run == "3519"      : return 250
    elif run == "3520"      : return 250
    elif run == "3522"      : return 310
    elif run == "3523"      : return 700
    elif run == "3526"      : return 280
    elif run == "3527"      : return 280
    elif run == "3528"      : return 280
    elif run == "3530"      : return 280
    elif run == "3539"      : return 110
    elif run == "3545"      : return 250
    elif run == "3546"      : return 100
    elif run == "3547"      : return 300
    elif run == "3548"      : return 300
    elif run == "3549"      : return 300
    fatal("Dont have event range for %s" % (run))

def time_range(run):
    if   run == "3545-3547" : return 13
    elif run == "3545-3548" : return 20
    elif run == "3545-3549" : return 27
    elif run == "3545"      : return 5
    elif run == "3546"      : return 2.2
    elif run == "3547"      : return 5
    elif run == "3548"      : return 6
    elif run == "3549"      : return 8
    fatal("Dont have event range for %s" % (run))

def pacman_example(h2, event, board, pdf):

    rootlogon()
    pacman = h2.ProjectionX("pacman", h2.GetYaxis().FindBin(board), h2.GetYaxis().FindBin(board))
    style(pacman)
    pacman.SetLineColor(ROOT.kBlack)
    pacman.SetFillColor(ROOT.kRed)
    pacman.SetMaximum(90)
    pacman.GetXaxis().SetRangeUser(170, 490)
    pacman.GetYaxis().SetTitle("Strip charge [fC]")
    pacman.GetYaxis().SetTitleOffset(1.4)

    examp = ROOT.TLatex(0.49, 0.75, "Pacman clustering")
    blurb = ROOT.TLatex(0.49, 0.70, "Event %i, Board %i" % (event, board))
    clus0 = ROOT.TLatex(0.25, 0.80, "cluster")
    clus1 = ROOT.TLatex(0.60, 0.43, "cluster")
    clus2 = ROOT.TLatex(0.78, 0.43, "cluster")

    canv = ROOT.TCanvas("pacman", "pacman", 800, 800)
    canv.Draw()
    pacman.Draw("hist")
    for tex in [examp, blurb, clus0, clus1, clus2]:
        style(tex)
        tex.Draw()
    save(canv, pdf, "pacman_example.pdf")

def strip_charge_vs_channel(h2, board, pdf):

    rootlogon()
    ROOT.gStyle.SetPadRightMargin(0.20)
    canv = ROOT.TCanvas("charge_%s" % (board), "charge_%s" % (board), 800, 800)
    canv.Draw()

    h2.RebinY(4)
    h2.GetXaxis().SetTitle("Strip number")
    h2.GetYaxis().SetTitle("Charge [fC]")
    h2.GetZaxis().SetTitle("Events")
    style(h2)
    h2.GetZaxis().SetTitleOffset(1.6)
    h2.Draw("colz")

    latex = ROOT.TLatex(0.20, 0.95, "Board %s" % (board))
    latex.SetTextSize(0.040)
    latex.SetTextFont(42)
    latex.SetNDC()
    latex.Draw()

    save(canv, pdf, "strip_charge_mmfe%s.pdf" % (board))

def strip_pdo_vs_channel(h2, board, pdf):

    rootlogon()
    ROOT.gStyle.SetPadRightMargin(0.20)
    canv = ROOT.TCanvas("pdo_%s" % (board), "pdo_%s" % (board), 800, 800)
    canv.Draw()

    h2.RebinY(4)
    h2.GetXaxis().SetTitle("Strip number")
    h2.GetYaxis().SetTitle("PDO [counts]")
    h2.GetZaxis().SetTitle("Events")
    h2.GetYaxis().SetRangeUser(0, 1200)
    style(h2)
    h2.GetZaxis().SetTitleOffset(1.6)
    h2.Draw("colz")

    latex = ROOT.TLatex(0.20, 0.95, "Board %s" % (board))
    latex.SetTextSize(0.040)
    latex.SetTextFont(42)
    latex.SetNDC()
    latex.Draw()

    save(canv, pdf, "strip_pdo_mmfe%s.pdf" % (board))

def strip_tdo_vs_channel(h2, board, pdf):

    rootlogon()
    ROOT.gStyle.SetPadRightMargin(0.20)
    canv = ROOT.TCanvas("tdo_%s" % (board), "tdo_%s" % (board), 800, 800)
    canv.Draw()

    h2.RebinY(4)
    h2.GetXaxis().SetTitle("Strip number")
    h2.GetYaxis().SetTitle("TDO [counts]")
    h2.GetZaxis().SetTitle("Events")
    style(h2)
    h2.GetZaxis().SetTitleOffset(1.6)
    h2.Draw("colz")

    latex = ROOT.TLatex(0.20, 0.95, "Board %s" % (board))
    latex.SetTextSize(0.040)
    latex.SetTextFont(42)
    latex.SetNDC()
    latex.Draw()

    save(canv, pdf, "strip_tdo_mmfe%s.pdf" % (board))

def strip_tdoc_vs_channel(h2, board, pdf):

    rootlogon()
    ROOT.gStyle.SetPadRightMargin(0.20)
    canv = ROOT.TCanvas("tdoc_%s" % (board), "tdoc_%s" % (board), 800, 800)
    canv.Draw()

    # h2.RebinY(4)
    h2.GetXaxis().SetTitle("Strip number")
    h2.GetYaxis().SetTitle("TDO [ns]")
    h2.GetZaxis().SetTitle("Events")
    style(h2)
    h2.GetYaxis().SetRangeUser(-10.,60.)
    h2.GetZaxis().SetTitleOffset(1.6)
    h2.Draw("colz")

    latex = ROOT.TLatex(0.20, 0.95, "Board %s" % (board))
    latex.SetTextSize(0.040)
    latex.SetTextFont(42)
    latex.SetNDC()
    latex.Draw()

    save(canv, pdf, "strip_tdoc_mmfe%s.pdf" % (board))

def strip_zpos_vs_channel(h2, board, pdf, ontrack=False):

    rootlogon()
    ROOT.gStyle.SetPadRightMargin(0.20)
    canv = ROOT.TCanvas("%s_%s" % ("zdri" if ontrack else "zpos", board), 
                        "%s_%s" % ("zdri" if ontrack else "zpos", board), 
                        800, 800)
    canv.Draw()

    # h2.RebinY(4)
    h2.GetXaxis().SetTitle("Strip number")
    h2.GetYaxis().SetTitle("Drift position [mm]")
    h2.GetZaxis().SetTitle("Events")
    style(h2)
    h2.GetZaxis().SetTitleOffset(1.6)
    h2.Draw("colz")

    latex = ROOT.TLatex(0.20, 0.95, "Board %s" % (board))
    latex.SetTextSize(0.040)
    latex.SetTextFont(42)
    latex.SetNDC()
    latex.Draw()

    latex2 = ROOT.TLatex(0.60, 0.95, "For #muTPC" if ontrack else "All strips")
    latex2.SetTextSize(0.040)
    latex2.SetTextFont(42)
    latex2.SetNDC()
    latex2.Draw()

    save(canv, pdf, "strip_%s_mmfe%s.pdf" % ("zdri" if ontrack else "zpos", board))

    # -------- projection -----------
    rootlogon()
    ROOT.gStyle.SetPadLeftMargin(0.18)
    proj = h2.ProjectionY()
    style(proj)
    proj.GetXaxis().SetRangeUser(-8, 16)
    proj.GetYaxis().SetTitleOffset(2.0)
    proj.SetLineWidth(2)
    proj.SetLineColor(ROOT.kBlack)
    proj.SetFillColor(ROOT.kYellow)
    proj.GetYaxis().SetTitle("Strips")
    canv = ROOT.TCanvas("%s_%s_proj" % ("zdri" if ontrack else "zpos", board), 
                        "%s_%s_proj" % ("zdri" if ontrack else "zpos", board), 
                        800, 800)
    canv.Draw()
    proj.Draw("histsame")
    latex.SetX(latex.GetX()+0.1)
    latex2.SetX(latex2.GetX()+0.1)
    latex.Draw()
    latex2.Draw()
    save(canv, pdf, "strip_%s_mmfe%s_proj.pdf" % ("zdri" if ontrack else "zpos", board))

    # -------- profile -----------
    rootlogon()
    ROOT.gStyle.SetPadLeftMargin(0.18)
    prof = h2.ProfileX()
    style(prof)
    #prof.GetXaxis().SetRangeUser(-8, 16)
    prof.GetYaxis().SetRangeUser(-1.2, 6.6)
    #prof.GetYaxis().SetRangeUser(2., 2.7+.7)
    prof.GetYaxis().SetTitleOffset(1.6)
    prof.SetLineWidth(2)
    prof.SetLineColor(ROOT.kBlack)
    #prof.SetFillColor(ROOT.kYellow)
    prof.GetYaxis().SetTitle("Drift average [mm]")
    prof.GetXaxis().SetTitle("Strip number")
    prof.GetZaxis().SetTitle("Events")
    canv = ROOT.TCanvas("%s_%s_prof" % ("zdri" if ontrack else "zpos", board), 
                        "%s_%s_prof" % ("zdri" if ontrack else "zpos", board), 
                        800, 800)
    canv.Draw()
    prof.Draw()
    #prof.Draw("histsame")
    latex.SetX(latex.GetX()+0.1)
    latex2.SetX(latex2.GetX()+0.1)
    latex.Draw()
    latex2.Draw()
    save(canv, pdf, "strip_%s_mmfe%s_prof.pdf" % ("zdri" if ontrack else "zpos", board))

def strip_zres(h2, board, pdf, bc=None):

    rootlogon()
    zres = h2.ProjectionY()
    zres.GetXaxis().SetRangeUser(-10, 10)
    zres.GetXaxis().SetTitle("#Deltaz(#muTPC, track) [mm]")
    zres.GetYaxis().SetTitle("Events")
    style(zres)
    zres.SetLineWidth(2)
    zres.SetLineColor(ROOT.kBlack)
    zres.SetFillColor(ROOT.kGreen-7)

    canv = ROOT.TCanvas("canv", "canv", 800, 800)
    canv.Draw()
    zres.Draw("histsame")
    zres.GetYaxis().SetTitleOffset(1.9)
    if not bc:
        saveas = "zres_mmfe%s.pdf" % (board)
        text = "Board %i" % (board)
    else:
        saveas = "zres_mmfe%s_BC%s.pdf" % (board, bc)
        text = "Board %i and BC%%4 = %s" % (board, bc)
    latex = ROOT.TLatex(0.20, 0.95, text)
    style(latex)
    latex.Draw()
    save(canv, pdf, saveas)

def strip_zpos_vs_ztrack(h2, board, pdf):

    rootlogon()
    ROOT.gStyle.SetPadLeftMargin(0.12)
    ROOT.gStyle.SetPadRightMargin(0.20)
    canv = ROOT.TCanvas("zz_%s" % (board), "zz_%s" % (board), 800, 800)
    canv.Draw()

    h2.GetXaxis().SetTitle("z_{track} [mm]")
    h2.GetYaxis().SetTitle("z_{#muTPC} [mm]")
    h2.GetZaxis().SetTitle("Strips")
    style(h2)
    h2.GetYaxis().SetTitleOffset(1.2)
    h2.GetZaxis().SetTitleOffset(1.6)
    h2.Draw("colz")

    latex = ROOT.TLatex(0.20, 0.95, "Board %s" % (board))
    latex.SetTextSize(0.040)
    latex.SetTextFont(42)
    latex.SetNDC()
    latex.Draw()

    save(canv, pdf, "strip_zpos_vs_ztrack_mmfe%s.pdf" % (board))

def phi_vs_phi(h2, board, pdf):

    rootlogon()
    ROOT.gStyle.SetPadLeftMargin(0.12)
    ROOT.gStyle.SetPadRightMargin(0.20)
    canv = ROOT.TCanvas("phi_vs_phi_%s" % (board), "phi_vs_phi_%s" % (board), 800, 800)
    canv.Draw()

    h2.GetXaxis().SetRangeUser(-0.7, 0.7)
    h2.GetYaxis().SetRangeUser(-2, 2)
    h2.GetXaxis().SetTitle("Track slope")
    h2.GetYaxis().SetTitle("#muTPC slope")
    h2.GetZaxis().SetTitle("Events")
    style(h2)
    h2.GetYaxis().SetTitleOffset(1.2)
    h2.GetZaxis().SetTitleOffset(1.6)
    h2.Draw("colz")

    latex = ROOT.TLatex(0.20, 0.95, "Board %s" % (board))
    latex.SetTextSize(0.040)
    latex.SetTextFont(42)
    latex.SetNDC()
    latex.Draw()

    save(canv, pdf, "theta_vs_theta_mmfe%s.pdf" % (board))

def tpc_fitprob(h2, pdf):

    rootlogon()
    for bo in [0, 1, 6, 7]:

        canv = ROOT.TCanvas("fitprob_%s" % (bo), "fitprob_%s" % (bo), 800, 800)
        canv.Draw()

        proj = h2.ProjectionY(str(bo), h2.GetXaxis().FindBin(bo), h2.GetXaxis().FindBin(bo))
        style(proj)
        proj.SetLineWidth(2)
        proj.SetLineColor(ROOT.kBlack)
        proj.SetFillColor(ROOT.kMagenta-10)
        proj.GetXaxis().SetTitle("#muTPC fit probability")
        proj.GetYaxis().SetTitle("Events")
        proj.Draw("histsame")

        latex = ROOT.TLatex(0.20, 0.95, "Board %s" % (bo))
        latex.SetTextSize(0.040)
        latex.SetTextFont(42)
        latex.SetNDC()
        latex.Draw()

        save(canv, pdf, "tpc_fitprob_mmfe%s.pdf" % (bo))

def tpc_unc(h1, bo, pdf):

    rootlogon()

    canv = ROOT.TCanvas("fitunc_%s" % (bo), "fitunc_%s" % (bo), 800, 800)
    canv.Draw()

    style(h1)
    h1.SetLineWidth(2)
    h1.SetLineColor(ROOT.kBlack)
    h1.SetFillColor(ROOT.kGreen+1)
    h1.Rebin(2)
    h1.GetXaxis().SetRangeUser(-2, 30)
    h1.GetXaxis().SetTitle("#muTPC uncertainty [mm]")
    h1.GetYaxis().SetTitle("Events")
    h1.GetYaxis().SetTitleOffset(1.4)
    h1.Draw("histsame")

    latex = ROOT.TLatex(0.20, 0.95, "Board %s" % (bo))
    latex.SetTextSize(0.040)
    latex.SetTextFont(42)
    latex.SetNDC()
    latex.Draw()

    canv.SetLogy(1)
    save(canv, pdf, "tpc_uncertainty_mmfe%s.pdf" % (bo))
    canv.SetLogy(0)

def track_unc(h1, bo, pdf):

    rootlogon()

    canv = ROOT.TCanvas("trackunc_%s" % (bo), "trackunc_%s" % (bo), 800, 800)
    canv.Draw()

    style(h1)
    h1.SetLineWidth(2)
    h1.SetLineColor(ROOT.kBlack)
    h1.SetFillColor(ROOT.kBlue-4)
    h1.Rebin(2)
    h1.GetXaxis().SetRangeUser(0, 1.3)
    h1.GetXaxis().SetTitle("Track uncertainty [mm]")
    h1.GetYaxis().SetTitle("Events")
    h1.GetYaxis().SetTitleOffset(1.9)
    h1.Draw("histsame")

    latex = ROOT.TLatex(0.20, 0.95, "Board %s" % (bo))
    latex.SetTextSize(0.040)
    latex.SetTextFont(42)
    latex.SetNDC()
    latex.Draw()

    #canv.SetLogy(1)
    save(canv, pdf, "track_uncertainty_mmfe%s.pdf" % (bo))
    #canv.SetLogy(0)

def dx_ij_vs_theta(h2, fabs, pdf, dxrange, tag=None, onesided=True):

    if not fabs and -dxrange < h2.GetYaxis().GetBinCenter(1):
        fatal("This dxrange doesnt work (lo) :: %s " % (dxrange))
    if  dxrange > h2.GetYaxis().GetBinCenter(h2.GetNbinsY()):
        fatal("This dxrange doesnt work (hi) :: %s " % (dxrange))

    bary = "bary" in h2.GetName()
    btpc = "btpc" in h2.GetName()
    utpc = "utpc" in h2.GetName()
    comb = "comb" in h2.GetName()
    d01  = "01"   in h2.GetName()
    d67  = "67"   in h2.GetName()

    meas, boards = None, None
    if bary: meas = "bary"
    if btpc: meas = "bary"
    if utpc: meas = "uTPC"
    if comb: meas = "comb"
    if d01: boards = "01"
    if d67: boards = "67"

    rootlogon()
    ROOT.gStyle.SetPadLeftMargin(0.12)
    ROOT.gStyle.SetPadRightMargin(0.20)
    name = "diff_%s_%s_%s_%s" % (meas.replace("comb", "combined"), boards, dxrange if dxrange else "all", "10deg" if tag else "all")
    canv = ROOT.TCanvas(name, name, 800, 800)
    canv.Draw()

    h2.GetXaxis().SetTitle("#theta_{MM} [deg.]" if not fabs else "#left|#theta_{MM}#right| [deg.]")
    h2.GetYaxis().SetTitle("x_{%s,i} #minus x_{%s,j} [mm]" % (meas.replace("uTPC", "#muTPC"), meas.replace("uTPC", "#muTPC")))
    h2.GetZaxis().SetTitle("Events")
    #h2.GetZaxis().SetTitle("Clusters")

    if d01:
        h2.GetYaxis().SetTitle(h2.GetYaxis().GetTitle().replace(",i", ",0").replace(",j", ",1"))
    if d67:
        h2.GetYaxis().SetTitle(h2.GetYaxis().GetTitle().replace(",i", ",6").replace(",j", ",7"))
    if fabs:
        if tag:
            h2.GetXaxis().SetRangeUser(10, 25)
    else:
        if tag:
            h2.GetXaxis().SetRangeUser(-25, -10)

    h2.GetYaxis().SetRangeUser(-dxrange, dxrange)
    
    style(h2)
    h2.GetYaxis().SetTitleOffset(1.2)
    h2.GetZaxis().SetTitleOffset(1.6)
    h2.Draw("colz")

    latex = ROOT.TLatex(0.20, 0.95, "Board %s vs. %s" % (boards[0], boards[1]))
    latex.SetTextSize(0.040)
    latex.SetTextFont(42)
    latex.SetNDC()
    # latex.Draw()

    if tag:
        latex2 = ROOT.TLatex(0.60, 0.95, tag)
        latex2.SetTextSize(0.040)
        latex2.SetTextFont(42)
        latex2.SetNDC()
        latex2.Draw()

    save(canv, pdf, "%s.pdf" % (name))

    # ---- profile ----
    rootlogon()
    ROOT.gStyle.SetPadLeftMargin(0.14)
    if onesided:
        h2_ = h2.RebinX(8, h2.GetName()+"_rebinx") # AW
        prof = h2_.ProfileX("_pfx", h2.GetYaxis().FindBin(-dxrange), h2.GetYaxis().FindBin(dxrange), "s")
    else:
        prof = h2.ProfileX("_pfx", h2.GetYaxis().FindBin(-dxrange), h2.GetYaxis().FindBin(dxrange), "s")
    prof1 = None

    # find the rms devs + errors
    nbinsx = h2.GetXaxis().GetNbins()
    sigs = []
    sigerrs = []
    if onesided:
        for i in np.arange(0,nbinsx/8):
            h2_ = h2.ProjectionY("proj_%d"%(i*8+1), i*8+1, i*8+8) # AW
            h2_.GetXaxis().SetRangeUser(-dxrange, dxrange)
            sigs.append(h2_.GetRMS())
            sigerrs.append(h2_.GetRMSError())

    if onesided:
        name1 = prof.GetName() + "_onesided"
        title = ";%s;%s" % (prof.GetXaxis().GetTitle(), prof.GetYaxis().GetTitle())
        nbins = prof.GetNbinsX()
        xmin  = prof.GetXaxis().GetBinLowEdge(1)
        xmax  = prof.GetXaxis().GetBinLowEdge(nbins+1)
        prof1 = ROOT.TH1F(name1, title, nbins, xmin, xmax)
        for bin in xrange(1, nbins+1):
            #if prof.GetBinError(bin) > 0.01:
            #    prof1.SetBinContent(bin, prof.GetBinError(bin))
            #prof1.SetBinError(  bin, 0.000001)
            prof1.SetBinContent(bin, sigs[bin-1])
            prof1.SetBinError(bin, sigerrs[bin-1])
        prof = prof1

    style(prof)
    prof.SetMarkerStyle(20)
    prof.SetMarkerSize(1)
    prof.SetMarkerColor(ROOT.kBlack)
    prof.SetLineColor(ROOT.kBlack)
    if not onesided:
        if dxrange == 10:
            prof.GetYaxis().SetRangeUser(-2, 2)
        elif dxrange == 2:
            prof.GetYaxis().SetRangeUser(-2, 2)
        elif dxrange == 1:
            prof.GetYaxis().SetRangeUser(-2, 2)
        else:
            prof.GetYaxis().SetRangeUser(-2, 2)
    else:
        if dxrange > 10:
            prof.SetMaximum(2)
        else:
            prof.SetMaximum(1.2)
    prof.GetYaxis().SetTitleOffset(1.4)
    prof.GetYaxis().SetTitle("RMS of %s" % (h2.GetYaxis().GetTitle()))
    if onesided:
        prof.SetMinimum(0)
        if not tag:
            # prof.GetXaxis().SetRangeUser(-25, 15)
            prof.GetXaxis().SetRangeUser(-30, 25)
    for bin in xrange(1, prof.GetNbinsX()+1):
        if not onesided:
            prof.SetBinContent(bin, 0)
    canv = ROOT.TCanvas(name+"prof", name+"prof", 800, 800)
    canv.Draw()
    canv.SetGrid()

    if onesided:
        tmp = ROOT.gStyle.GetErrorX()
        ROOT.gStyle.SetErrorX(0.0001)
        #ROOT.gStyle.SetErrorX(1)
        prof.Draw("psame")
        ROOT.gStyle.SetErrorX(tmp)
        ROOT.gStyle.SetErrorX(0.0001)
    else:
        prof.Draw("same")
    # latex.Draw()
    if tag:
        latex2.SetX(latex2.GetX()+0.15)
        latex2.Draw()

    save(canv, pdf, "%s_%s.pdf" % (name, "prof"))

    # ---- projection ----
    rootlogon()
    ROOT.gStyle.SetPadLeftMargin(0.18)
    proj = h2.ProjectionY()
    style(proj)
    proj.SetLineWidth(2)
    proj.SetLineColor(ROOT.kBlack)
    proj.SetFillColor(ROOT.kCyan-7)
    proj.GetXaxis().SetRangeUser(-dxrange, dxrange)
    #rms = proj.GetRMS()
    #proj.GetXaxis().SetRangeUser(-5, 5)
    proj.Rebin(2)
    proj.GetYaxis().SetTitle("Events / %2.2f mm"%(proj.GetBinWidth(1)))
    #proj.GetYaxis().SetTitle("Clusters")
    proj.GetYaxis().SetTitleOffset(2.0)
    canv = ROOT.TCanvas(name+"proj", name+"proj", 800, 800)
    canv.Draw()
    proj.Draw("histsame")
    # latex.Draw()
    if tag:
        latex2.Draw()

    latexRMS1 = ROOT.TLatex(0.23, 0.70, "RMS = %.2f mm" % (proj.GetRMS()))
    #latexRMS2 = ROOT.TLatex(0.23, 0.65, "from %i#minus%i" % (-dxrange, dxrange))
    for latexRMS in [latexRMS1]:
        latexRMS.SetTextSize(0.040)
        latexRMS.SetTextFont(42)
        latexRMS.SetNDC()
        latexRMS.Draw()

    save(canv, pdf, "%s_%s.pdf" % (name, "proj"))

def residuals_1D(h1, tag, pdf, rms=False):

    if tag == "barybary01": xaxis = "(x_{bary.,0} #minus x_{bary.,1}) / #sqrt{#sigma^{2}_{bary.,0} + #sigma^{2}_{bary.,1}}"
    if tag == "barybary67": xaxis = "(x_{bary.,6} #minus x_{bary.,7}) / #sqrt{#sigma^{2}_{bary.,6} + #sigma^{2}_{bary.,7}}"
    if tag == "barybary06": xaxis = "(x_{bary.,0} #minus x_{bary.,6}) / #sqrt{#sigma^{2}_{bary.,0} + #sigma^{2}_{bary.,6}}"
    if tag == "barybary17": xaxis = "(x_{bary.,1} #minus x_{bary.,7}) / #sqrt{#sigma^{2}_{bary.,1} + #sigma^{2}_{bary.,7}}"
    if tag == "utpcutpc01": xaxis = "(x_{#muTPC,0} #minus x_{#muTPC,1}) / #sqrt{#sigma^{2}_{#muTPC,0} + #sigma^{2}_{#muTPC,1}}"
    if tag == "utpcutpc67": xaxis = "(x_{#muTPC,6} #minus x_{#muTPC,7}) / #sqrt{#sigma^{2}_{#muTPC,6} + #sigma^{2}_{#muTPC,7}}"
    if tag == "utpcutpc06": xaxis = "(x_{#muTPC,0} #minus x_{#muTPC,6}) / #sqrt{#sigma^{2}_{#muTPC,0} + #sigma^{2}_{#muTPC,6}}"
    if tag == "utpcutpc17": xaxis = "(x_{#muTPC,1} #minus x_{#muTPC,7}) / #sqrt{#sigma^{2}_{#muTPC,1} + #sigma^{2}_{#muTPC,7}}"
    if tag == "utpcbary00": xaxis = "(x_{#muTPC,0} #minus x_{bary.,0}) / #sqrt{#sigma^{2}_{#muTPC,0} + #sigma^{2}_{bary.,0}}"
    if tag == "utpcbary11": xaxis = "(x_{#muTPC,1} #minus x_{bary.,1}) / #sqrt{#sigma^{2}_{#muTPC,1} + #sigma^{2}_{bary.,1}}"
    if tag == "utpcbary66": xaxis = "(x_{#muTPC,6} #minus x_{bary.,6}) / #sqrt{#sigma^{2}_{#muTPC,6} + #sigma^{2}_{bary.,6}}"
    if tag == "utpcbary77": xaxis = "(x_{#muTPC,7} #minus x_{bary.,7}) / #sqrt{#sigma^{2}_{#muTPC,7} + #sigma^{2}_{bary.,7}}"
        
    rootlogon()
    style(h1)
    h1.SetLineWidth(2)
    h1.SetLineColor(ROOT.kBlack)
    h1.SetFillColor(ROOT.kRed-7 if not "bary" in tag else ROOT.kViolet-9)
    # rms = h1.GetRMS()
    ROOT.gStyle.SetPadBottomMargin(0.14)
    h1.Rebin(4)
    h1.GetXaxis().SetTitle(xaxis)
    h1.GetYaxis().SetTitle("Events / %2.2f" %(h1.GetBinWidth(1)))
    h1.GetYaxis().SetTitleOffset(1.8)
    h1.GetXaxis().SetTitleOffset(1.3)
    canv = ROOT.TCanvas(tag, tag, 800, 800)
    canv.Draw()
    h1.Draw("histsame")

    latexRMS = ROOT.TLatex(0.23, 0.95, "RMS = %.2f mm" % (h1.GetRMS()))
    latexRMS.SetTextSize(0.040)
    latexRMS.SetTextFont(42)
    latexRMS.SetNDC()
    if rms:
        latexRMS.Draw()

    save(canv, pdf, "diff_%s.pdf" % (tag))
    
def strip_duplicates_vs_channel(h2, board, pdf):

    rootlogon()
    dups = h2.ProjectionX("dupsx", h2.GetYaxis().FindBin(board), h2.GetYaxis().FindBin(board))
    dups.GetXaxis().SetTitle("Strip number")
    dups.GetYaxis().SetTitle("Number of duplicate strips")
    style(dups)
    dups.SetLineColor(ROOT.kBlack)
    dups.SetFillColor(ROOT.kYellow)

    canv = ROOT.TCanvas("canv", "canv", 800, 800)
    canv.Draw()
    dups.Draw("hist")
    dups.GetYaxis().SetTitleOffset(1.4)
    latex = ROOT.TLatex(0.20, 0.95, "Board %i" % (board))
    style(latex)
    latex.Draw()
    canv.SetLogy(1)
    save(canv, pdf, "duplicates_mmfe%s.pdf" % (board))
    canv.SetLogy(0)

def clusters_per_board(h2, pdf):

    rootlogon()
    canv = ROOT.TCanvas("canv", "canv", 800, 800)
    canv.Draw()

    clusmult = h2.ProjectionY()
    clusmult.GetXaxis().SetTitle("Clusters per board")
    clusmult.GetYaxis().SetTitle("Boards")
    style(clusmult)
    clusmult.SetLineWidth(2)
    clusmult.SetLineColor(ROOT.kBlack)
    clusmult.SetFillColor(210)
    clusmult.GetXaxis().SetRangeUser(-0.5, 5.5)
    clusmult.Draw("hist")

    latex = ROOT.TLatex(0.60, 0.66, "All boards")
    style(latex)
    latex.SetTextSize(0.050)
    latex.Draw()
    save(canv, pdf, "clusters_per_board_lin.pdf")
    
def strips_per_cluster(h2, pdf):

    rootlogon()
    canv = ROOT.TCanvas("canv", "canv", 800, 800)
    canv.Draw()

    stripmult = h2.ProjectionY()
    stripmult.GetXaxis().SetTitle("Strips per cluster")
    stripmult.GetYaxis().SetTitle("Clusters")
    style(stripmult)
    stripmult.SetLineWidth(2)
    stripmult.SetLineColor(ROOT.kBlack)
    stripmult.SetFillColor(ROOT.kViolet)
    stripmult.GetXaxis().SetRangeUser(-0.5, 12.5)
    stripmult.Draw("hist")

    latex = ROOT.TLatex(0.60, 0.66, "All boards")
    style(latex)
    latex.SetTextSize(0.050)
    latex.Draw()
    save(canv, pdf, "strips_per_cluster_lin.pdf")

def strips_per_cluster_vs_theta(h2, pdf):

    rootlogon()
    ROOT.gStyle.SetPadRightMargin(0.24)
    canv = ROOT.TCanvas("canv", "canv", 800, 800)
    canv.Draw()
    h2.GetXaxis().SetTitle("#theta#lower[0.5]{#scale[0.7]{track}} [degrees]")
    h2.GetYaxis().SetTitle("Strips per cluster")
    h2.GetZaxis().SetTitle("Clusters")
    style(h2)
    h2.GetYaxis().SetTitleOffset(1.5)
    h2.GetZaxis().SetTitleOffset(2.0)
    h2.GetXaxis().SetRangeUser(-30, 25)
    #h2.GetYaxis().SetRangeUser(-0.5, 10.5)
    h2.Draw("colzsame")
    #save(canv, pdf, "strips_per_cluster_vs_theta.pdf")

    rootlogon()
    h2.RebinX(4)
    pfx = h2.ProfileX("_pfx", 0, h2.GetNbinsY()+1, "")
    style(pfx)
    pfx.GetYaxis().SetLabelOffset(0.01)
    pfx.GetYaxis().SetTitle("Mean strips per cluster")
    pfx.GetYaxis().SetTitleOffset(1.5)
    pfx.GetXaxis().SetRangeUser(-30, 25)
    pfx.SetMarkerStyle(20)
    pfx.SetMarkerSize(1.3)
    pfx.SetMarkerColor(ROOT.kBlue)
    pfx.SetLineColor(ROOT.kBlue)
    pfx.SetLineWidth(2)
    pfx.SetMaximum(6.0)
    pfx.SetMinimum(0)
    canv = ROOT.TCanvas("canvpfx", "canvpfx", 800, 800)
    canv.SetGrid()
    canv.Draw()
    pfx.Draw("same")
    ROOT.gPad.RedrawAxis()
    save(canv, pdf, "strips_per_cluster_vs_theta_profile.pdf")

def clusters_on_track(h2, pdf):

    ops = options()
    rootlogon()
    ROOT.gStyle.SetPadRightMargin(0.24)
    canv = ROOT.TCanvas("canv", "canv", 800, 800)
    canv.Draw()
    # h2.RebinX(4)
    h2.GetXaxis().SetTitle("Time [days]")
    h2.GetYaxis().SetTitle("Clusters per track")
    h2.GetZaxis().SetTitle("Tracks")
    h2.GetYaxis().SetLabelOffset(0.01)
    style(h2)
    h2.GetYaxis().SetTitleOffset(1.5)
    h2.GetZaxis().SetTitleOffset(2.0)
    h2.GetXaxis().SetNdivisions(505)
    h2.GetXaxis().SetRangeUser(0, time_range(ops.r))
    h2.GetYaxis().SetRangeUser(3.5, 9.5)
    h2.Draw("colzsame")
    save(canv, pdf, "clusters_per_track_vs_event.pdf")

    rootlogon()
    ROOT.gStyle.SetPadLeftMargin(0.20)
    proj = h2.ProjectionY()
    style(proj)
    proj.GetYaxis().SetTitleOffset(2.0)
    proj.GetXaxis().SetTitle("Clusters per track")
    proj.GetYaxis().SetTitle("Tracks")
    proj.SetLineColor(ROOT.kBlack)
    proj.SetFillColor(ROOT.kAzure+7)
    proj.SetLineWidth(2)
    canv = ROOT.TCanvas("canvproj", "canvproj", 800, 800)
    canv.Draw()
    proj.Draw("same")
    ROOT.gPad.RedrawAxis()
    save(canv, pdf, "clusters_per_track.pdf")

def dbcid_mmfe_scint(h2, pdf, board=None):

    if not board:
        timediff_vs_board = h2.ProjectionY()
    else:
        timediff_vs_board = h2.ProjectionY(str(board), h2.GetXaxis().FindBin(board), h2.GetXaxis().FindBin(board))
    timediff_vs_board.GetXaxis().SetTitle("Trigger BCID #minus strip BCID")
    timediff_vs_board.GetYaxis().SetTitle("Strips")
    style(timediff_vs_board)
    timediff_vs_board.SetLineWidth(2)
    timediff_vs_board.SetLineColor(ROOT.kBlack)
    timediff_vs_board.SetFillColor(ROOT.kOrange+1)
    timediff_vs_board.GetXaxis().SetRangeUser(-0.5, 60.5)

    canv = ROOT.TCanvas("canv", "canv", 800, 800)
    canv.Draw()
    timediff_vs_board.Draw("hist")

    latex = ROOT.TLatex(0.22, 0.75, "All boards" if not board else "Board %s" % (board))
    style(latex)
    latex.SetTextSize(0.050)
    latex.Draw()

    if not board:
        save(canv, pdf, "bcid_strip_v_trigger_lin.pdf")
        canv.SetLogy(1)
        save(canv, pdf, "bcid_strip_v_trigger_log.pdf")
        canv.SetLogy(0)
    else:
        save(canv, pdf, "bcid_strip_v_trigger_lin_%s.pdf" % (board))
        canv.SetLogy(1)
        save(canv, pdf, "bcid_strip_v_trigger_log_%s.pdf" % (board))
        canv.SetLogy(0)

def track_angle(h2, nclus, pdf):

    rootlogon()
    canv = ROOT.TCanvas("canv", "canv", 800, 800)
    canv.Draw()
    angle = h2.ProjectionX("track_param_x_1D")
    angle.GetYaxis().SetTitle("Tracks")
    angle.GetXaxis().SetTitle("#theta#lower[0.5]{#scale[0.7]{track}} [degrees]")
    style(angle)
    angle.Rebin(2)
    angle.SetLineWidth(2)
    angle.SetLineColor(ROOT.kBlack)
    angle.SetFillColor(ROOT.kAzure+1 if nclus==7 else ROOT.kOrange+1)
    angle.GetXaxis().SetRangeUser(-30, 25)
    angle.Draw("hist")
    latex = ROOT.TLatex(0.20, 0.95, "%i tracks with %i clusters" % (angle.Integral(), nclus))
    style(latex)
    latex.Draw()
    save(canv, pdf, "track_angle_%shits.pdf" % (nclus))

def residuals_N1(h2, board_a, board_b, pdf, art=True, tag=None, tpc=False, norm=False):

    rootlogon()
    residual_a = h2.ProjectionY("residual_a", h2.GetXaxis().FindBin(board_a), h2.GetXaxis().FindBin(board_a))
    residual_b = h2.ProjectionY("residual_b", h2.GetXaxis().FindBin(board_b), h2.GetXaxis().FindBin(board_b))
    residuals = [residual_a, residual_b]
    colors = {0: ROOT.kBlack, 1: ROOT.kRed, 6: ROOT.kBlue, 7: 210}
    latexs = []

    gaus_a = ROOT.TF1("gaus_a", "gaus(0)")
    gaus_b = ROOT.TF1("gaus_b", "gaus(0)")
    gausses = [gaus_a, gaus_b]

    for (residual, gaus) in zip(residuals, gausses):
        residual.Rebin(2)
        board  = board_a if residual==residual_a else board_b
        offset = 0       if residual==residual_a else 1
        style(residual)
        if art:
            det = "x_{ART}"
        elif tpc:
            det = "x_{#muTPC}"
        else:
            det = "x_{cluster}"
        residual.GetYaxis().SetTitle("Events / %2.2f" % (residual.GetXaxis().GetBinWidth(1)) )
        #residual.GetYaxis().SetTitle("Tracks with #geq6 clusters")
        residual.GetXaxis().SetTitle("%s #minus x_{track, proj.} [mm]" % (det))
        if norm:
            title = residual.GetXaxis().GetTitle()
            title = title.replace("cluster", "bary.")
            title = title.replace(" [mm]", "")
            title = "(%s)" % (title)
            title = title + " / #sigma_{bary.}"
            residual.GetXaxis().SetTitle(title)
        residual.SetLineWidth(3)
        residual.SetMarkerStyle(8)
        residual.SetMarkerSize(0.8)
        residual.SetMarkerColor(colors[board])
        #residual.SetLineWidth(2)
        residual.SetLineColor(ROOT.kBlack)
        mean = residual.GetMean()
        rms  = residual.GetRMS()
        show_overflow(residual)
        residual.SetLineColor(colors[board])
        if norm:
            residual.GetXaxis().SetRangeUser(-5, 5)
        else:
            residual.GetXaxis().SetRangeUser(-2.4, 2.4)
        gaus.SetRange(mean-0.3, mean+0.3)
        residual.Fit(gaus, "QRNM")
        gaus.SetLineWidth(4)
        gaus.SetLineColor(colors[board])
        sigma = gaus.GetParameter(2)
        bo   = ROOT.TLatex(0.24, 0.70-0.06*offset, "Board %i" % (board))
        # mean = ROOT.TLatex(0.65, 0.70-0.06*offset, "%6.2f"   % (mean))
        # rms  = ROOT.TLatex(0.80, 0.70-0.06*offset, "%6.2f"   % (rms))
        rms   = ROOT.TLatex(0.645, 0.70-0.06*offset, "%6.2f"   % (rms))
        sigma = ROOT.TLatex(0.80,  0.70-0.06*offset, "%6.2f"   % (sigma))
        for latex in [bo, rms, sigma]:
            if norm:
                if latex==sigma:
                    continue
                if latex==rms:
                    latex.SetX(latex.GetX()+0.09)
            latex.SetTextColor(colors[board])
            latexs.append(latex)
        
    # latexs.append( ROOT.TLatex(0.66, 0.86, "Mean      RMS") )
    if not norm:
        latexs.append( ROOT.TLatex(0.66,  0.76, "RMS      #sigma#lower[0.5]{#scale[0.7]{gaus.}}") )
        latexs.append( ROOT.TLatex(0.735, 0.58, "[mm]") )
    else:
        latexs.append( ROOT.TLatex(0.75,  0.76, "RMS") )

    if tag=="10deg":
        latexs.append(ROOT.TLatex(0.18, 0.95, "#theta < 10 deg."))
    if tag==">10deg":
        latexs.append(ROOT.TLatex(0.18, 0.95, "#theta > 10 deg."))

    canv = ROOT.TCanvas("canv", "canv", 800, 800)
    canv.Draw()
    for residual in sorted(residuals, key=lambda h: h.GetMaximum(), reverse=True):
        residual.Draw("psame")
    for gaus in gausses:
        if norm:
            continue
        gaus.Draw("same")
    for latex in latexs:
        style(latex)
        latex.Draw()

    if art:
        det = "art" 
    elif tpc:
        det = "utpc"
    else:
        det = "clus"

    if norm:  name = "track_residuals_boards%s%s_%s_%s.pdf" % (board_a, board_b, det, "norm")
    elif tag: name = "track_residuals_boards%s%s_%s_%s.pdf" % (board_a, board_b, det, tag.replace(">", ""))
    else:     name = "track_residuals_boards%s%s_%s.pdf"    % (board_a, board_b, det)
        
    save(canv, pdf, name)

def residuals_vs_theta(h2, pdf, tpc=False, comp=False):

    rootlogon()
    obj = "cluster" if not tpc else "uTPC"

    # th2
    ROOT.gStyle.SetPadRightMargin(0.18)
    style(h2)
    h2.GetXaxis().SetRangeUser(-27, 16)
    # h2.GetYaxis().SetRangeUser(-11, 11)
    h2.GetXaxis().SetTitle("#theta#lower[0.5]{#scale[0.7]{track}} [degrees]")
    h2.GetYaxis().SetTitle("x#lower[0.5]{#scale[0.7]{%s}} #minus x#lower[0.5]{#scale[0.7]{track, proj.}} [mm]" % (obj))
    h2.GetZaxis().SetTitle("Events")
    h2.GetXaxis().SetNdivisions(505)
    h2.GetYaxis().SetTitleOffset(1.6)
    board_number = h2.GetName()[-1]
    text_board = "All boards" if board_number=="N" else "Board %s" % (board_number)
    texBo = ROOT.TLatex(0.20, 0.95, text_board)
    style(texBo)
    canv = ROOT.TCanvas("canv", "canv", 800, 800)
    canv.Draw()
    h2.Draw("colzsame")
    texBo.Draw()
    ROOT.gPad.RedrawAxis()
    if not comp:
        save(canv, pdf, "residuals_vs_theta_%s_%s.pdf" % (obj, board_number))

    # profile
    rootlogon()
    colors = {0: ROOT.kBlack,
              1: ROOT.kRed, 
              2: ROOT.kOrange,
              3: ROOT.kViolet,
              4: ROOT.kCyan,
              5: ROOT.kGray+1,
              6: ROOT.kBlue, 
              7: 210}
    # pfx = h2.ProfileX("_pfx", 0, h2.GetNbinsY()+1)
    pfx = h2.ProfileX("_pfx")
    style(pfx)
    (ymax, ymin) = (1.5, -1.5)
    pfx.GetYaxis().SetRangeUser(ymin, ymax)
    pfx.GetXaxis().SetRangeUser(-27, 16)
    pfx.GetYaxis().SetTitle("Mean of x#lower[0.5]{#scale[0.7]{%s}} #minus x#lower[0.5]{#scale[0.7]{track, proj.}} [mm]" % (obj))
    # pfx.GetXaxis().SetTitle("track slope (xz)")
    pfx.GetXaxis().SetNdivisions(505)
    pfx.GetYaxis().SetTitleOffset(1.6)
    pfx.SetMarkerStyle(20)
    pfx.SetMarkerSize(1.3)
    pfx.SetMarkerColor(colors[int(board_number)])
    pfx.SetLineColor(colors[int(board_number)])
    pfx.SetLineWidth(2)
    canv = ROOT.TCanvas("canvpfx", "canvpfx", 800, 800)
    canv.SetGrid()
    canv.Draw()
    pfx.Draw("same")
    texBo.Draw()
    ROOT.gPad.RedrawAxis()
    if not comp:
        save(canv, pdf, "residuals_vs_theta_%s_%s_profile.pdf" % (obj, board_number))

    # proj
    if comp:
        name = "residuals_vs_theta_comp_proj_%s" % (board_number)
        proj = h2.ProjectionY()
        style(proj)
        proj.GetXaxis().SetTitle("(x_{#muTPC} #minus x_{bary.}) / #sqrt{#sigma^{2}_{#muTPC} + #sigma^{2}_{bary.}}")
        proj.GetXaxis().SetRangeUser(-5, 5)
        proj.GetXaxis().SetTitleOffset(1.3)
        ROOT.gStyle.SetPadBottomMargin(0.14)
        #proj.GetXaxis().SetRangeUser(-8, 16)
        proj.SetLineWidth(2)
        proj.SetLineColor(ROOT.kBlack)
        proj.SetFillColor(ROOT.kBlue-10)
        proj.GetYaxis().SetTitle("Clusters")
        canv = ROOT.TCanvas(name, name, 800, 800)
        canv.Draw()
        proj.Draw("histsame")
        texBo.Draw()
        save(canv, pdf, "%s.pdf" % (name))

def residuals_vs_x(h2, pdf, angled=False):

    rootlogon()

    # th2
    ROOT.gStyle.SetPadRightMargin(0.18)
    style(h2)
    # h2.GetXaxis().SetRangeUser(-27, 16)
    h2.GetYaxis().SetRangeUser(-11, 11)
    h2.GetXaxis().SetTitle("x#lower[0.5]{#scale[0.7]{cluster}} [mm]")
    h2.GetYaxis().SetTitle("x#lower[0.5]{#scale[0.7]{cluster}} #minus x#lower[0.5]{#scale[0.7]{track, proj.}} [mm]")
    h2.GetZaxis().SetTitle("Events")
    h2.GetXaxis().SetNdivisions(505)
    h2.GetYaxis().SetTitleOffset(1.6)
    board_number = h2.GetName()[-1]
    text_board = "All boards" if board_number=="N" else "Board %s" % (board_number)
    if angled:
        text_board = "%s, #theta > 10^{o}" % (text_board)
    texBo = ROOT.TLatex(0.20, 0.95, text_board)
    style(texBo)
    canv = ROOT.TCanvas("canv", "canv", 800, 800)
    canv.Draw()
    h2.Draw("colzsame")
    texBo.Draw()
    ROOT.gPad.RedrawAxis()
    save(canv, pdf, "residuals_vs_x_%s.pdf" % (board_number if not angled else "10deg_%s" % (board_number)))

    # profile
    rootlogon()
    colors = {0: ROOT.kBlack,
              1: ROOT.kRed, 
              2: ROOT.kOrange,
              3: ROOT.kViolet,
              4: ROOT.kCyan,
              5: ROOT.kGray+1,
              6: ROOT.kBlue, 
              7: 210}
    pfx = h2.ProfileX("_pfx", 0, h2.GetNbinsY()+1)
    style(pfx)
    (ymax, ymin) = (1.5, -1.5)
    pfx.GetYaxis().SetRangeUser(ymin, ymax)
    # pfx.GetXaxis().SetRangeUser(-27, 16)
    pfx.GetYaxis().SetTitle("Mean of x#lower[0.5]{#scale[0.7]{cluster}} #minus x#lower[0.5]{#scale[0.7]{track, proj.}} [mm]")
    # pfx.GetXaxis().SetTitle("track slope (xz)")
    pfx.GetXaxis().SetNdivisions(505)
    pfx.GetYaxis().SetTitleOffset(1.6)
    pfx.SetMarkerStyle(20)
    pfx.SetMarkerSize(1.3)
    pfx.SetMarkerColor(colors[int(board_number)])
    pfx.SetLineColor(colors[int(board_number)])
    pfx.SetLineWidth(2)
    canv = ROOT.TCanvas("canvpfx", "canvpfx", 800, 800)
    canv.SetGrid()
    canv.Draw()
    pfx.Draw("same")
    texBo.Draw()
    ROOT.gPad.RedrawAxis()
    save(canv, pdf, "residuals_vs_x_%s_profile.pdf" % (board_number if not angled else "10deg_%s" % (board_number)))

def xy_difference(h2, tag, pdf):

    rootlogon()
    ROOT.gStyle.SetPadRightMargin(0.21)
    style(h2)
    h2.GetXaxis().SetTitle("#Deltax [mm]")
    h2.GetYaxis().SetTitle("#Deltay [mm]")
    h2.GetZaxis().SetTitle("Events with 8-cluster tracks")
    h2.GetZaxis().SetTitleOffset(1.85)
    xlatex = ROOT.TLatex(0.18, 0.96, "RMS, x = %.2f" % (h2.GetRMS(1)))
    ylatex = ROOT.TLatex(0.55, 0.96, "RMS, y = %i"   % (h2.GetRMS(2)))

    canv = ROOT.TCanvas("canv", "canv", 800, 800)
    canv.Draw()
    h2.Draw("colz")
    for latex in [xlatex, ylatex]:
        style(latex)
        latex.Draw()
    save(canv, pdf, "dxdy_%s.pdf" % (tag))

def scintillator_slopes(rfile, pdf):

    rootlogon()
    slopes_20 = rfile.Get("histograms/track_angle_x_scint_20")
    slopes_21 = rfile.Get("histograms/track_angle_x_scint_21")
    slopes2D = copy.copy(slopes_20)
    slopes2D.Add(slopes_21)

    colors = [ROOT.kBlack, ROOT.kBlue, ROOT.kRed, 210, ROOT.kViolet, ROOT.kOrange+7]
    latexs = []
    slopes = []
    for isl in xrange(6):
        slopes.append( slopes2D.ProjectionX(str(isl), isl+2, isl+2) )
        x = slopes[-1].GetMean()
        y = slopes[-1].GetMaximum()*1.1
        # if isl == 0: y *= 1.5
        if isl == 0: x -= 4
        if isl == 1: x -= 2
        if isl == 2: x -= 2
        if isl == 3: x += 0
        if isl == 4: x += 0
        if isl == 5: x += 2
        slopes[-1].GetXaxis().SetTitle("#theta#lower[0.5]{#scale[0.7]{track}} [degrees]")
        slopes[-1].GetYaxis().SetTitle("Tracks")
        slopes[-1].GetYaxis().SetTitleOffset(1.2)
        slopes[-1].SetLineWidth(2)
        slopes[-1].GetXaxis().SetRangeUser(-28, 12)
        slopes[-1].SetMaximum(slopes[-1].GetMaximum()*1.2)
        latexs.append( ROOT.TLatex(x, y, "scint. %i" % (isl)) )
        slopes[-1].SetLineColor(colors[isl])
        latexs[-1].SetTextColor(colors[isl])
        style(slopes[-1])
        style(latexs[-1])
        latexs[-1].SetTextAlign(22)
        latexs[-1].SetNDC(0)

    canv = ROOT.TCanvas("track_angle_x_scint", "track_angle_x_scint", 800, 800)
    canv.Draw()
    for slope in sorted(slopes, key=lambda h: h.GetMaximum(), reverse=True):
        slope.Draw("histsame")
    for latex in latexs:
        latex.Draw()
    save(canv, pdf, "scint_angle.pdf")

def scintillator_slope_vs_time(isc, scint, pdf):

    ops = options()
    rootlogon()
    ROOT.gStyle.SetPadRightMargin(0.16)
    ROOT.gStyle.SetPadLeftMargin(0.15)

    style(scint)
    scint.GetXaxis().SetRangeUser(0, time_range(ops.r))
    scint.GetYaxis().SetRangeUser(-30, 20)
    scint.GetXaxis().SetTitle("Time [days]")
    scint.GetYaxis().SetTitle("#theta#lower[0.5]{#scale[0.7]{track}} [degrees]")
    scint.GetZaxis().SetTitle("Tracks")
    scint.GetYaxis().SetTitleOffset(1.3)
    scint.GetZaxis().SetTitleOffset(1.2)
    scint.GetXaxis().SetNdivisions(505)

    tex1 = ROOT.TLatex(0.20, 0.95, "Run %s" % (ops.r))
    tex2 = ROOT.TLatex(0.58, 0.95, "Bottom scint. %i" % (isc))
    texs = [tex1, tex2]

    canv = ROOT.TCanvas("scint_vs_time_%i" % (isc), "scint_vs_time_%i" % (isc), 800, 800)
    canv.Draw()
    scint.Draw("colzsame")
    for tex in texs:
        style(tex)
        tex.Draw()
    save(canv, pdf, "scint_vs_time_%i.pdf" % (isc))

def vmm_efficiency(numer, denom, pdf, run=None):

    rootlogon()
    ROOT.gStyle.SetPadRightMargin(0.18)
    ROOT.gStyle.SetPadLeftMargin(0.15)

    eff = copy.copy(numer)
    style(eff)
    eff.Divide(numer, denom, 1.0, 1.0, "B")
    eff.Scale(100)
    eff.GetXaxis().SetRangeUser(0, 7)
    eff.GetXaxis().SetTitle("VMM number")
    eff.GetYaxis().SetTitle("Board number")
    eff.GetZaxis().SetTitle("Efficiency [%]")
    eff.GetYaxis().SetTitleOffset(1.3)
    eff.GetYaxis().SetLabelOffset(0.01)
    eff.SetMaximum(100)
    eff.SetMinimum(0)
    eff.SetMarkerSize(1.3)
    ROOT.gStyle.SetPaintTextFormat(".0f")

    ncontours = 200
    stops = array.array("d", [0.00, 0.50, 0.75, 1.00])
    red   = array.array("d", [1.00, 1.00, 1.00, 0.00])
    blue  = array.array("d", [0.00, 0.00, 0.00, 0.00])
    green = array.array("d", [0.00, 0.66, 1.00, 1.00])
    ROOT.TColor.CreateGradientColorTable(len(stops), stops, red, green, blue, ncontours)
    ROOT.gStyle.SetNumberContours(ncontours)

    tex = ROOT.TLatex(0.20, 0.95, "Run %s" % (run)) if run else None

    canv = ROOT.TCanvas("canv", "canv", 800, 800)
    canv.Draw()
    eff.Draw("colztextsame")
    if tex:
        style(tex)
        tex.Draw()
    save(canv, pdf, "vmm_efficiency.pdf")
    ROOT.gStyle.SetPalette(1)

def strip_lifetime(h2, board, ymax, pdf):

    rootlogon()
    ROOT.gStyle.SetPadRightMargin(0.18)
    tex = ROOT.TLatex(0.20, 0.95, "Board %i" % (board))
    style(tex)
    style(h2)
    h2.RebinY(8)
    h2.GetXaxis().SetTitle("Strip number")
    h2.GetYaxis().SetTitle("Event number / 1000")
    h2.GetZaxis().SetTitle("Hits recorded")
    h2.GetYaxis().SetRangeUser(0, ymax) # 310

    canv = ROOT.TCanvas("canv", "canv", 800, 800)
    canv.Draw()
    h2.Draw("colzsame")
    tex.Draw()
    save(canv, pdf, "strip_hits_vs_event_%s.pdf" % (board))

def trigger_art_vs_event(h2, pdf):

    ops = options()
    rootlogon()
    ROOT.gStyle.SetPadRightMargin(0.20)
    style(h2)
    h2.RebinX(4)
    h2.GetXaxis().SetTitle("Event number / 10^{3}")
    h2.GetYaxis().SetTitle("Number of hits in trigger")
    h2.GetZaxis().SetTitle("Triggers")
    h2.GetYaxis().SetTitleOffset(1.3)
    h2.GetZaxis().SetTitleOffset(1.6)
    h2.GetYaxis().SetLabelOffset(0.01)
    h2.GetXaxis().SetRangeUser(0, event_range(ops.r))
    h2.GetYaxis().SetRangeUser(2.5, 9.5)
    canv = ROOT.TCanvas("canv", "canv", 800, 800)
    canv.Draw()
    h2.Draw("colz")
    save(canv, pdf, "trigger_hits_vs_event.pdf")

def art_position_within_cluster(h2, pdf):

    rootlogon()
    def metric(name):
        if "_hit_" in name: return "channel"
        if "_bci_" in name: return "BCID"
        if "_pdo_" in name: return "PDO"

    # COLZ
    ncontours = 200
    stops = array.array("d", [0.00, 0.33, 0.66, 1.00])
    red   = array.array("d", [1.00, 1.00, 1.00, 0.80])
    blue  = array.array("d", [1.00, 0.00, 0.00, 0.00])
    green = array.array("d", [1.00, 1.00, 0.66, 0.00])
    ROOT.TColor.CreateGradientColorTable(len(stops), stops, red, green, blue, ncontours)
    ROOT.gStyle.SetNumberContours(ncontours)

    # 2D
    style(h2)
    ROOT.gStyle.SetPadRightMargin(0.18)
    ROOT.gStyle.SetPadLeftMargin(0.15)
    h2.GetYaxis().SetTitleOffset(1.3)
    h2.GetYaxis().SetLabelOffset(0.01)
    h2.GetXaxis().SetTitle("ART position in MMFE cluster, by %s" % (metric(h2.GetName())))
    h2.GetYaxis().SetTitle("Cluster size [strips]")
    ROOT.gStyle.SetPaintTextFormat("g");
    canv = ROOT.TCanvas("canv", "canv", 800, 800)
    canv.Draw()
    h2.Draw("colzsametext35")
    ROOT.gPad.RedrawAxis()
    save(canv, pdf, "trigger_index_in_cluster_%s.pdf" % (metric(h2.GetName())))

    # pretty colors
    colors = [ROOT.kBlue+1, ROOT.kAzure-3, ROOT.kAzure-9, ROOT.kRed-9, ROOT.kRed-4, ROOT.kRed+2]

    # 1D overlays
    rootlogon()
    stack = ROOT.THStack("stack", ";%s;Events" % (h2.GetXaxis().GetTitle()))
    clus1 = h2.ProjectionX("clus1", h2.GetYaxis().FindBin(1), h2.GetYaxis().FindBin(1))
    clus2 = h2.ProjectionX("clus2", h2.GetYaxis().FindBin(2), h2.GetYaxis().FindBin(2))
    clus3 = h2.ProjectionX("clus3", h2.GetYaxis().FindBin(3), h2.GetYaxis().FindBin(3))
    clus4 = h2.ProjectionX("clus4", h2.GetYaxis().FindBin(4), h2.GetYaxis().FindBin(4))
    clus5 = h2.ProjectionX("clus5", h2.GetYaxis().FindBin(5), h2.GetYaxis().FindBin(5))
    clus6 = h2.ProjectionX("clus6", h2.GetYaxis().FindBin(6), h2.GetNbinsY())
    clusters = [clus1, clus2, clus3, clus4, clus5, clus6]
    for (icl, clus) in enumerate(clusters):
        style(clus)
        clus.SetLineWidth(2)
        clus.SetLineColor(ROOT.kBlack)
        clus.SetFillColor(colors[icl])
        stack.Add(clus)

    stack.Draw()
    stack.GetXaxis().SetRangeUser(-0.5, 5.5)
    style(stack)
    ROOT.gStyle.SetPadLeftMargin(0.20)
    stack.GetYaxis().SetTitleOffset(2.2)

    legend = ROOT.TLegend(0.53, 0.55, 0.73, 0.82)
    style(legend)
    legend.AddEntry(clus1, " 1 strip", "f")
    legend.AddEntry(clus2, " 2 strips", "f")
    legend.AddEntry(clus3, " 3 strips", "f")
    legend.AddEntry(clus4, " 4 strips", "f")
    legend.AddEntry(clus5, " 5 strips", "f")
    legend.AddEntry(clus6, " #geq6 strips", "f")

    canv = ROOT.TCanvas("canv2", "canv2", 800, 800)
    canv.Draw()
    #for clus in sorted(clusters, key=lambda h: h.Integral(), reverse=True):
    #    clus.Draw("histsame")
    #for clus in sorted(clusters, key=lambda h: h.Integral()):
    #    clus.Draw("histsame")
    stack.Draw()
    legend.Draw()
    ROOT.gPad.RedrawAxis()
    save(canv, pdf, "trigger_index_in_cluster_%s_stack.pdf" % (metric(h2.GetName())))
    ROOT.gStyle.SetPalette(1)

    # keep it simple: 3 hits in cluster
    rootlogon()
    canv = ROOT.TCanvas("canv3", "canv3", 800, 800)
    canv.Draw()
    clus3.GetXaxis().SetRangeUser(-0.5, 4.5)
    clus3.SetNdivisions(505)
    clus3.Draw("histsame")
    save(canv, pdf, "trigger_index_in_cluster_%s_clus3.pdf" % (metric(h2.GetName())))

    # keep it simple: 4 hits in cluster
    rootlogon()
    canv = ROOT.TCanvas("canv4", "canv4", 800, 800)
    canv.Draw()
    clus4.GetXaxis().SetRangeUser(-0.5, 4.5)
    clus4.SetNdivisions(505)
    clus4.Draw("histsame")
    save(canv, pdf, "trigger_index_in_cluster_%s_clus4.pdf" % (metric(h2.GetName())))

def art_dx_vs_theta(h2, pdf):

    rootlogon()

    # th2
    ROOT.gStyle.SetPadRightMargin(0.18)
    style(h2)
    h2.RebinX(2)
    h2.RebinY(2)
    h2.GetXaxis().SetRangeUser(-27, 16)
    h2.GetYaxis().SetRangeUser(-11, 11)
    h2.GetXaxis().SetTitle("#theta#lower[0.5]{#scale[0.7]{MMFE track}} [deg.]")
    h2.GetYaxis().SetTitle("x#lower[0.5]{#scale[0.7]{ART}} #minus x#lower[0.5]{#scale[0.7]{MMFE track}} [mm]")
    h2.GetZaxis().SetTitle("Events")
    h2.GetXaxis().SetNdivisions(505)
    h2.GetYaxis().SetTitleOffset(1.6)
    board_number = h2.GetName()[-1]
    text_board = "All boards" if board_number=="N" else "Board %s" % (board_number)
    texRo = ROOT.TLatex(0.20, 0.95, "Road size: #pm24 strips")
    texBo = ROOT.TLatex(0.68, 0.95, text_board)
    style(texRo)
    style(texBo)
    canv = ROOT.TCanvas("canv", "canv", 800, 800)
    canv.Draw()
    h2.Draw("colzsame")
    texRo.Draw()
    texBo.Draw()
    ROOT.gPad.RedrawAxis()
    save(canv, pdf, "trigger_art_dx_vs_theta_%s.pdf" % (board_number))

    # profile
    rootlogon()
    pfx = h2.ProfileX("_pfx", 0, h2.GetNbinsY()+1, "s")
    style(pfx)
    (ymax, ymin) = (3, -3)
    pfx.GetYaxis().SetRangeUser(ymin, ymax)
    pfx.GetXaxis().SetRangeUser(-27, 16)
    pfx.GetYaxis().SetTitle("Mean, RMS of x#lower[0.5]{#scale[0.7]{ART}} #minus x#lower[0.5]{#scale[0.7]{MMFE track}} [mm]")
    pfx.GetXaxis().SetNdivisions(505)
    pfx.GetYaxis().SetTitleOffset(1.6)
    pfx.SetMarkerStyle(20)
    pfx.SetMarkerSize(1.3)
    pfx.SetMarkerColor(ROOT.kBlue)
    pfx.SetLineColor(ROOT.kBlue)
    pfx.SetLineWidth(2)
    canv = ROOT.TCanvas("canvpfx", "canvpfx", 800, 800)
    canv.SetGrid()
    canv.Draw()
    pfx.Draw("same")
    texRo.Draw()
    texBo.Draw()
    ROOT.gPad.RedrawAxis()
    save(canv, pdf, "trigger_art_dx_vs_theta_%s_profile.pdf" % (board_number))
   
def art_dtheta_vs_theta(h2, tag, pdf):

    rootlogon()

    # th2
    ROOT.gStyle.SetPadRightMargin(0.18)
    # ROOT.gStyle.SetPalette(-1)
    style(h2)
    h2.RebinX(2)
    h2.RebinY(2)
    # h2.GetYaxis().SetRangeUser(-30, 30)
    h2.GetXaxis().SetRangeUser(-27, 16)
    h2.GetXaxis().SetTitle("#theta#lower[0.5]{#scale[0.7]{MMFE}} [deg.]")
    h2.GetYaxis().SetTitle("#theta#lower[0.5]{#scale[0.7]{MMTP}} #minus #theta#lower[0.5]{#scale[0.7]{MMFE}} [mrad]")
    h2.GetZaxis().SetTitle("Events")
    h2.GetXaxis().SetNdivisions(505)
    h2.GetYaxis().SetTitleOffset(1.6)
    if h2.GetName().endswith("_all"):
        road = "default (3 VMM)"
    else:
        road = "#pm%s strips" % (h2.GetName().split("near")[1])
    tex = ROOT.TLatex(0.20, 0.95, "Road size: %s" % (road))
    style(tex)
    canv = ROOT.TCanvas("canv", "canv", 800, 800)
    canv.Draw()
    h2.Draw("colzsame")
    tex.Draw()
    ROOT.gPad.RedrawAxis()
    save(canv, pdf, "trigger_dtheta_vs_theta_%s.pdf" % (tag))

    # profile
    rootlogon()
    pfx = h2.ProfileX("_pfx", 0, h2.GetNbinsY()+1, "s")
    style(pfx)
    (ymax, ymin) = (35, -35) if h2.GetName().endswith("_all") else (15, -15)
    pfx.GetYaxis().SetRangeUser(ymin, ymax)
    pfx.GetXaxis().SetRangeUser(-27, 16)
    pfx.GetYaxis().SetTitle("Mean, RMS of #theta#lower[0.5]{#scale[0.7]{MMTP}} #minus #theta#lower[0.5]{#scale[0.7]{MMFE}} [mrad]")
    pfx.GetXaxis().SetNdivisions(505)
    pfx.GetYaxis().SetTitleOffset(1.6)
    pfx.SetMarkerStyle(20)
    pfx.SetMarkerSize(1.3)
    pfx.SetMarkerColor(210)
    pfx.SetLineColor(210)
    pfx.SetLineWidth(2)
    canv = ROOT.TCanvas("canvpfx", "canvpfx", 800, 800)
    canv.SetGrid()
    canv.Draw()
    pfx.Draw("same")
    tex.Draw()
    ROOT.gPad.RedrawAxis()
    save(canv, pdf, "trigger_dtheta_vs_theta_%s_profile.pdf" % (tag))
   
def art_dtheta(rfile, pdf, road=None, nx=None):

    road_tag = "all" if road==None else "near"+road
    nxxx_tag = "NX"  if   nx==None else nx
    histpath = "histograms/trig_dtheta_%s_%s" % (road_tag, nxxx_tag)

    rootlogon()
    ROOT.gStyle.SetOptStat(111111)
    hist = rfile.Get(histpath)
    if not hist:
        fatal("Couldnt find %s" % (histpath))
    style(hist)
    hist.SetLineColor(ROOT.kBlue)
    hist.SetLineWidth(2)
    hist.GetYaxis().SetTitle("Events")
    hist.GetXaxis().SetTitle("#theta#lower[0.5]{#scale[0.7]{MMTP}} #minus #theta#lower[0.5]{#scale[0.7]{MMFE}} [mrad]")
    hist.GetXaxis().SetRangeUser(-80, 80)

    sel = "All triggers" if nx==None else "%s triggers" % (nx)
    texNX = ROOT.TLatex(0.21, 0.7, sel)
    style(texNX)
    texNX.SetTextColor(ROOT.kBlue)
    texNX.SetTextSize(0.040)

    sel = "No road req." if road==None else "Road: #pm%s strips" % (road)
    texRo = ROOT.TLatex(0.21, 0.64, sel)
    style(texRo)
    texRo.SetTextColor(ROOT.kBlue)
    texRo.SetTextSize(0.040)

    canv = ROOT.TCanvas("canv", "canv", 800, 800)
    canv.Draw()

    hist.Draw("histsame")
    texNX.Draw()
    texRo.Draw()

    ROOT.gPad.RedrawAxis()
    save(canv, pdf, "trigger_dtheta.pdf")

def art_dtheta_rms(rfile, pdf, run=None):

    rootlogon()
    ROOT.gStyle.SetLabelSize(0.05, 'xyz')
    ROOT.gStyle.SetTitleSize(0.05, 'xyz')
    ROOT.gStyle.SetPadTopMargin(0.06)
    ROOT.gStyle.SetPadRightMargin(0.05)
    ROOT.gStyle.SetPadBottomMargin(0.16)
    ROOT.gStyle.SetPadLeftMargin(0.16)
    ROOT.gStyle.SetTitleOffset(1.2, 'x')
    ROOT.gStyle.SetTitleOffset(1.2, 'y')

    nears = [64, 48, 36, 24, 16, 12, 8, 4]
    rmsNX = []
    rms3X = []
    rms4X = []
    for near in nears:
        h = rfile.Get("histograms/trig_dtheta_near%02i_NX" % (near))
        rmsNX.append( h.GetRMS() )
        h = rfile.Get("histograms/trig_dtheta_near%02i_3X" % (near))
        rms3X.append( h.GetRMS() )
        h = rfile.Get("histograms/trig_dtheta_near%02i_4X" % (near))
        rms4X.append( h.GetRMS() )

    import array
    xs = array.array("d", [near*2 for near in nears])
    titlex = "Effective road size [strips]"
    titley = "RMS of #theta#lower[0.5]{#scale[0.7]{local}} #minus #theta#lower[0.5]{#scale[0.7]{MM}} [mrad]"

    # 3X
    ys = array.array("d", rms3X)
    gr3X = ROOT.TGraph(len(xs), xs, ys)
    gr3X.SetTitle(";%s;%s" % (titlex, titley))
    gr3X.SetMarkerColor(ROOT.kRed)
    gr3X.SetLineColor(ROOT.kRed)
    gr3X.SetLineWidth(2)
    gr3X.SetMarkerStyle(20)
    gr3X.SetMarkerSize(1.3)
    # 4X
    ys = array.array("d", rms4X)
    gr4X = ROOT.TGraph(len(xs), xs, ys)
    gr4X.SetTitle(";%s;%s" % (titlex, titley))
    gr4X.SetMarkerColor(210)
    gr4X.SetLineColor(210)
    gr4X.SetLineWidth(2)
    gr4X.SetMarkerStyle(20)
    gr4X.SetMarkerSize(1.3)
    # NX
    ys = array.array("d", rmsNX)
    grNX = ROOT.TGraph(len(xs), xs, ys)
    grNX.SetTitle(";%s;%s" % (titlex, titley))
    grNX.SetMarkerColor(210)
    grNX.SetLineColor(210)
    grNX.SetLineWidth(3)
    grNX.SetMarkerStyle(20)
    grNX.SetMarkerSize(1.3)

    texs = []
    for iw,word in enumerate(["All triggers", "3X triggers", "4X triggers"]):
        texs.append(ROOT.TLatex(0.25, 0.8-0.06*iw, word))
        style(texs[-1])
        if "All" in word: texs[-1].SetTextColor(ROOT.kBlack)
        if "3X"  in word: texs[-1].SetTextColor(ROOT.kRed)
        if "4X"  in word: texs[-1].SetTextColor(210)
    texs = []
    lines = []

    texs.append( ROOT.TLatex(0.20, 0.95, "Harvard CRTS") )
    texs.append( ROOT.TLatex(0.75, 0.95, "Run %s" % (run)) )

    texs.append( ROOT.TLatex(0.40, 0.75, "1/4 VMM") )
    texs.append( ROOT.TLatex(0.41, 0.72, "#pm1 road") )
    texs.append( ROOT.TLatex(0.25, 0.75, "1/8 VMM") )
    texs.append( ROOT.TLatex(0.26, 0.72, "#pm1 road") )
    lines.append( ROOT.TLine(0.444, 0.69, 0.444, 0.55) )
    lines.append( ROOT.TLine(0.305, 0.69, 0.305, 0.43) )

    canv = ROOT.TCanvas("canv", "canv", 800, 800)
    canv.Draw()

    multi = ROOT.TMultiGraph()
    multi.Add(grNX, "PL")
    # multi.Add(gr3X, "PL")
    # multi.Add(gr4X, "PL")
    multi.SetTitle(grNX.GetTitle())
    multi.SetMinimum(0)
    multi.SetMaximum(13)
    multi.Draw("A")
    multi.GetXaxis().SetLimits(-1, 135)
    multi.Draw("A")
    for tex in texs:
        style(tex)
        if tex in texs[-4:]:
            tex.SetTextSize(0.030)
        tex.Draw()
    for line in lines:
        style(line)
        line.SetLineWidth(1)
        line.Draw()

    ROOT.gPad.RedrawAxis()
    save(canv, pdf, "trigger_dtheta_rms_vs_roadsize.pdf")

def art_dtheta_efficiency(rfile, cut, pdf, run=None):

    rootlogon()
    ROOT.gStyle.SetLabelSize(0.05, 'xyz')
    ROOT.gStyle.SetTitleSize(0.05, 'xyz')
    ROOT.gStyle.SetPadTopMargin(0.06)
    ROOT.gStyle.SetPadRightMargin(0.05)
    ROOT.gStyle.SetPadBottomMargin(0.16)
    ROOT.gStyle.SetPadLeftMargin(0.16)
    ROOT.gStyle.SetTitleOffset(1.2, 'x')
    ROOT.gStyle.SetTitleOffset(1.2, 'y')

    nears  = [64, 48, 36, 24, 16, 12, 8, 4]
    effsNX = []
    effs3X = []
    effs4X = []

    for near in nears:
        # NX
        h = rfile.Get("histograms/trig_dtheta_near%02i_NX" % (near))
        numer = h.Integral(h.GetXaxis().FindBin(-cut), h.GetXaxis().FindBin(cut))
        denom = h.Integral(0, h.GetNbinsX()+1)
        effsNX.append(100*numer/denom)
        # 3X
        h = rfile.Get("histograms/trig_dtheta_near%02i_3X" % (near))
        numer = h.Integral(h.GetXaxis().FindBin(-cut), h.GetXaxis().FindBin(cut))
        denom = h.Integral(0, h.GetNbinsX()+1)
        effs3X.append(100*numer/denom)
        # 4X
        h = rfile.Get("histograms/trig_dtheta_near%02i_4X" % (near))
        numer = h.Integral(h.GetXaxis().FindBin(-cut), h.GetXaxis().FindBin(cut))
        denom = h.Integral(0, h.GetNbinsX()+1)
        effs4X.append(100*numer/denom)

    import array
    xs = array.array("d", [near*2 for near in nears])
    dtheta = "#theta#lower[0.5]{#scale[0.7]{local}} #minus #theta#lower[0.5]{#scale[0.7]{MM}}"
    titlex = "Effective road size [strips]"
    titley = "Efficiency of %s < 15 mrad [%%]" % (dtheta)

    # 3X
    ys = array.array("d", effs3X)
    gr3X = ROOT.TGraph(len(xs), xs, ys)
    gr3X.SetTitle(";%s;%s" % (titlex, titley))
    gr3X.SetMarkerColor(ROOT.kRed)
    gr3X.SetLineColor(ROOT.kRed)
    gr3X.SetLineWidth(2)
    gr3X.SetMarkerStyle(20)
    gr3X.SetMarkerSize(1.3)
    # 4X
    ys = array.array("d", effs4X)
    gr4X = ROOT.TGraph(len(xs), xs, ys)
    gr4X.SetTitle(";%s;%s" % (titlex, titley))
    gr4X.SetMarkerColor(210)
    gr4X.SetLineColor(210)
    gr4X.SetLineWidth(2)
    gr4X.SetMarkerStyle(20)
    gr4X.SetMarkerSize(1.3)
    # NX
    ys = array.array("d", effsNX)
    grNX = ROOT.TGraph(len(xs), xs, ys)
    grNX.SetTitle(";%s;%s" % (titlex, titley))
    grNX.SetMarkerColor(ROOT.kBlue)
    grNX.SetLineColor(ROOT.kBlue)
    grNX.SetLineWidth(3)
    grNX.SetMarkerStyle(20)
    grNX.SetMarkerSize(1.3)

    texs = []
    for iw,word in enumerate(["All triggers", "3X triggers", "4X triggers"]):
        texs.append(ROOT.TLatex(0.65, 0.8-0.06*iw, word))
        style(texs[-1])
        if "All" in word: texs[-1].SetTextColor(ROOT.kBlack)
        if "3X"  in word: texs[-1].SetTextColor(ROOT.kRed)
        if "4X"  in word: texs[-1].SetTextColor(210)
    texs = []
    lines = []

    texs.append( ROOT.TLatex(0.20, 0.95, "Harvard CRTS") )
    texs.append( ROOT.TLatex(0.75, 0.95, "Run %s" % (run)) )

    texs.append( ROOT.TLatex(0.40, 0.25, "1/4 VMM") )
    texs.append( ROOT.TLatex(0.41, 0.22, "#pm1 road") )
    texs.append( ROOT.TLatex(0.25, 0.25, "1/8 VMM") )
    texs.append( ROOT.TLatex(0.26, 0.22, "#pm1 road") )
    lines.append( ROOT.TLine(0.444, 0.31, 0.444, 0.49) )
    lines.append( ROOT.TLine(0.305, 0.31, 0.305, 0.69) )

    canv = ROOT.TCanvas("canv", "canv", 800, 800)
    canv.Draw()

    multi = ROOT.TMultiGraph()
    # multi.Add(gr3X, "PL")
    # multi.Add(gr4X, "PL")
    multi.Add(grNX, "PL")
    multi.SetTitle(grNX.GetTitle())
    multi.SetMinimum(93)
    multi.SetMaximum(100.9)
    multi.Draw("A")
    multi.GetXaxis().SetLimits(-1, 135)
    multi.Draw("A")
    for tex in texs:
        style(tex)
        if tex in texs[-4:]:
            tex.SetTextSize(0.030)
        tex.Draw()
    for line in lines:
        style(line)
        line.SetLineWidth(1)
        line.Draw()

    ROOT.gPad.RedrawAxis()
    save(canv, pdf, "trigger_dtheta_eff_vs_roadsize.pdf")

def art_theta(th, pdf):
    rootlogon()
    style(th)
    th.GetXaxis().SetRangeUser(-30, 20)
    th.GetXaxis().SetTitle("#theta#lower[0.5]{#scale[0.7]{trigger}} [deg.]")
    th.GetYaxis().SetTitle("Events")
    th.SetLineColor(ROOT.kBlack)
    th.SetFillColor(ROOT.kYellow)
    th.SetLineWidth(3)
    canv = ROOT.TCanvas("art_theta", "art_theta", 800, 800)
    canv.Draw()
    ROOT.gPad.SetLogy(0)
    th.Draw("histsame")
    ROOT.gPad.RedrawAxis()
    save(canv, pdf, "trigger_theta.pdf")

def art_n(th1, th2=None, pdf=None):
    rootlogon()
    ths = [th1, th2] if th2 else [th1]
    for th in ths:
        style(th)
        th.GetXaxis().SetTitle("Number of hits in track")
        th.GetYaxis().SetTitle("Events")
        th.SetLineWidth(3)
        th.GetXaxis().SetRangeUser(1.5, 9.5)
    th1.SetLineColor(ROOT.kCyan-6)
    th1.SetFillColor(0)
    th1.SetMaximum(th1.GetMaximum()*1.07)
    if th2:
        th2.SetLineColor(ROOT.kRed-9)

    legend = ROOT.TLegend(0.27, 0.69, 0.47, 0.82)
    style(legend)
    legend.AddEntry(th2, " MM", "f")
    legend.AddEntry(th1, " TP", "f")

    canv = ROOT.TCanvas("art_n", "art_n", 800, 800)
    canv.Draw()
    ROOT.gPad.SetLogy(0)
    for th in ths:
        th.Draw("histsame")
    legend.Draw()
    ROOT.gPad.RedrawAxis()
    save(canv, pdf, "trigger_nart.pdf")

def art_pairs(th, pdf):
    rootlogon()
    style(th)
    th.GetXaxis().SetTitle("#DeltaBC(ART pairs)")
    th.GetYaxis().SetTitle("Events")
    th.SetLineColor(ROOT.kBlack)
    th.SetFillColor(ROOT.kAzure-4)
    th.SetLineWidth(3)

    gaus = ROOT.TF1("gaus", "gaus(0)")
    gaus.SetRange(-3.7, 3.7)
    th.Fit(gaus, "QRNM")
    gaus.SetLineWidth(3)
    gaus.SetLineColor(ROOT.kRed)

    tex = ROOT.TLatex(0.67, 0.8, "#sigma = %.2f BC" % (gaus.GetParameter(2)))
    style(tex)
    tex.SetTextColor(ROOT.kRed)
    tex.SetTextSize(0.05)

    canv = ROOT.TCanvas("art_pairs", "art_pairs", 800, 800)
    canv.Draw()
    ROOT.gPad.SetLogy(0)
    th.Draw("histsame")
    gaus.Draw("same")
    tex.Draw()
    ROOT.gPad.RedrawAxis()
    save(canv, pdf, "trigger_dbc_pairs.pdf")

def art_window_vs_theta(th2, pdf):
    rootlogon()
    th = th2.ProjectionX("x", 0, th2.GetYaxis().FindBin(-10))
    style(th)
    # th.GetXaxis().SetRangeUser(-30, 20)
    th.GetXaxis().SetTitle("BC Window of Trigger ARTs")
    th.GetYaxis().SetTitle("Events")
    th.SetLineColor(ROOT.kBlack)
    th.SetFillColor(ROOT.kYellow)
    th.SetLineWidth(3)
    tex = ROOT.TLatex(0.22, 0.81, "#theta#lower[0.5]{#scale[0.7]{MMTP}} > 10 deg.")
    style(tex)
    canv = ROOT.TCanvas("art_window", "art_window", 800, 800)
    canv.Draw()
    ROOT.gPad.SetLogy(0)
    th.Draw("histsame")
    tex.Draw()
    ROOT.gPad.RedrawAxis()
    save(canv, pdf, "trigger_dbc_vs_theta.pdf")

def art_window(th2, pdf):
    rootlogon()
    th = th2.ProjectionX()
    style(th)
    # th.GetXaxis().SetRangeUser(-30, 20)
    th.GetXaxis().SetTitle("BC Window of Trigger ARTs")
    th.GetYaxis().SetTitle("Events")
    th.SetLineColor(ROOT.kBlack)
    th.SetFillColor(ROOT.kOrange-4)
    th.SetLineWidth(3)
    canv = ROOT.TCanvas("art_window", "art_window", 800, 800)
    canv.Draw()
    ROOT.gPad.SetLogy(0)
    th.Draw("histsame")
    ROOT.gPad.RedrawAxis()
    save(canv, pdf, "trigger_dbc.pdf")

def bcwindow_vs_N(th2, pdf, norm=False):
    rootlogon()
    proj = {}
    texs = {}
    for art in xrange(4, 9):
        bin = th2.GetYaxis().FindBin(art)
        proj[art] = th2.ProjectionX(str(art), bin, bin)
        proj[art].SetLineColor(color(art-3))
        proj[art].GetXaxis().SetTitle("BC Window")
        proj[art].GetYaxis().SetTitle("Events, normalized" if norm else "Events")
        show_overflow(proj[art])
        if norm:
            print art, proj[art].Integral(0, proj[art].GetNbinsX()+1)
            proj[art].Scale(1/proj[art].Integral(0, proj[art].GetNbinsX()+1))
            proj[art].SetMaximum(0.4)
        #else:
        #    proj[art].SetMaximum(7000)
        proj[art].SetLineWidth(3)
        style(proj[art])
        texs[art] = ROOT.TLatex(0.25, 1.1-art*0.06, "N(ART) = %i" % (art))
        texs[art].SetTextColor(color(art-3))
        style(texs[art])
    
    canv = ROOT.TCanvas("bcwindow_stacked", "bcwindow_stacked", 800, 800)
    canv.Draw()
    ROOT.gPad.SetLogy(0)
    for art in sorted(proj.values(), key=lambda h: h.GetMaximum(), reverse=True):
        art.Draw("histsame")
    for tex in sorted(texs.keys()): 
        texs[tex].Draw()

    ROOT.gPad.RedrawAxis()
    save(canv, pdf, "trigger_dbc_vs_nart.pdf")

def x_resolution(h2, pdf):

    # projection
    rootlogon()
    proj = h2.ProjectionY()
    show_overflow(proj)
    style(proj)
    proj.SetLineColor(ROOT.kBlack)
    proj.SetFillColor(ROOT.kViolet-2)
    proj.GetXaxis().SetTitle("#LTx#GT#lower[0.5]{#scale[0.7]{MMTP}} #minus #LTx#GT#lower[0.5]{#scale[0.7]{MMFE}} [mm]")
    proj.GetYaxis().SetTitle("Events")
    proj.SetLineWidth(3)
    proj.GetXaxis().SetRangeUser(-4, 4)
    tex0 = ROOT.TLatex(0.24, 0.75, "RMS: %.2f mm" % (proj.GetRMS()))
    tex1 = ROOT.TLatex(0.20, 0.95, "Road size: #pm24 strips")
    texs = [tex0, tex1]

    canv = ROOT.TCanvas("canv", "canv", 800, 800)
    canv.Draw()
    proj.Draw("histsame")
    for tex in texs:
        style(tex)
        tex.Draw()
    ROOT.gPad.SetLogy(1)
    ROOT.gPad.RedrawAxis()
    save(canv, pdf, "trigger_x_residual_log.pdf")
    ROOT.gPad.SetLogy(0)
    ROOT.gPad.RedrawAxis()
    save(canv, pdf, "trigger_x_residual_lin.pdf")

    # profile
    rootlogon()
    h2.RebinX(10)
    pfx = h2.ProfileX("_pfx", 0, h2.GetNbinsY()+1, "s")
    style(pfx)
    pfx.GetYaxis().SetRangeUser(-1.5, 1.5)
    pfx.GetXaxis().SetRangeUser(-27, 16)
    pfx.GetYaxis().SetTitle("Mean, RMS of %s" % (proj.GetXaxis().GetTitle()))
    pfx.GetXaxis().SetTitle("#theta#lower[0.5]{#scale[0.7]{MMTP}} [degrees]")
    pfx.GetXaxis().SetNdivisions(505)
    pfx.GetYaxis().SetTitleOffset(1.6)
    pfx.SetMarkerStyle(20)
    pfx.SetMarkerSize(1.3)
    pfx.SetMarkerColor(ROOT.kViolet-2)
    pfx.SetLineColor(ROOT.kViolet-2)
    pfx.SetLineWidth(2)
    canv = ROOT.TCanvas("canvpfx", "canvpfx", 800, 800)
    canv.SetGrid()
    canv.Draw()
    pfx.Draw("same")
    tex.Draw()
    ROOT.gPad.RedrawAxis()
    save(canv, pdf, "trigger_x_residual_vs_theta.pdf")

def theta_resolution_all(th2, tag, pdf, norm=False):
    rootlogon()
    proj = th2.ProjectionX() if isinstance(th2, ROOT.TH2) else th2
    show_overflow(proj)
    style(proj)
    proj.SetLineColor(ROOT.kBlack)
    proj.SetFillColor(8)
    proj.GetXaxis().SetTitle("#theta#lower[0.5]{#scale[0.7]{MMTP}} #minus #theta#lower[0.5]{#scale[0.7]{MMFE}} [mrad]")
    proj.GetYaxis().SetTitle("Events")
    proj.Rebin(2)
    proj.SetLineWidth(3)
    if isinstance(th2, ROOT.TH2):
        proj.GetXaxis().SetRangeUser(-150, 150)
    else:
        proj.GetXaxis().SetRangeUser(-60, 60)

    texs = ROOT.TLatex(0.23, 0.78, "RMS: %.1f mrad" % (proj.GetRMS()))
    texs.SetTextColor(ROOT.kBlack)
    style(texs)
    accept15 = proj.Integral(proj.FindBin(-15), proj.FindBin(15)) / proj.Integral(0, proj.GetNbinsX()+1)
    a15s = ROOT.TLatex(0.65, 0.78, "%4.1f%%" % (accept15*100))
    a15s.SetTextColor(ROOT.kBlack)
    style(a15s)
    accept99 = proj.Integral(proj.FindBin(-99), proj.FindBin(99)) / proj.Integral(0, proj.GetNbinsX()+1)
    a99s = ROOT.TLatex(0.80, 0.78, "%4.1f%%" % (accept99*100))
    a99s.SetTextColor(ROOT.kBlack)
    style(a99s)
    roadsize = "default (3 VMM)" if isinstance(th2, ROOT.TH2) else "#pm24 strips"
    road = ROOT.TLatex(0.20, 0.95, "Road size: %s" % (roadsize))
    road.SetTextColor(ROOT.kBlack)
    style(road)

    legs = [ROOT.TLatex(0.63, 0.84, "#Delta#theta < 15"),
            ROOT.TLatex(0.78, 0.84, "#Delta#theta < 100")]
    for leg in legs:
        style(leg)

    lines = [ROOT.TLine(0.62, 0.83, 0.76, 0.83),
             ROOT.TLine(0.78, 0.83, 0.92, 0.83)]
    for line in lines:
        line.SetLineWidth(2)
        style(line)

    canv = ROOT.TCanvas("theta_resolution1D_lin", "theta_resolution1D_lin", 800, 800)
    canv.Draw()
    ROOT.gPad.SetLogy(0)
    proj.Draw("histsame")
    texs.Draw()
    a15s.Draw()
    a99s.Draw()
    road.Draw()
    for obj in legs+lines: 
        obj.Draw()

    ROOT.gPad.SetLogy(0)
    save(canv, pdf, "trigger_theta_resolution_%s_lin.pdf" % (tag))
    ROOT.gPad.SetLogy(1)
    save(canv, pdf, "trigger_theta_resolution_%s_log.pdf" % (tag))

def theta_resolution(th2, pdf, norm=False):
    rootlogon()
    proj = {}
    texs = {}
    a15s = {}
    a99s = {}
    for art in xrange(2, 5):
        bin = th2.GetYaxis().FindBin(art)
        proj[art] = th2.ProjectionX(str(art), bin, bin)
        proj[art].SetLineColor(color(art-1))
        proj[art].GetXaxis().SetTitle("#theta#lower[0.5]{#scale[0.7]{MMTP}} #minus #theta#lower[0.5]{#scale[0.7]{MMFE}} [mrad]")
        proj[art].GetYaxis().SetTitle("Events, normalized" if norm else "Events")
        show_overflow(proj[art])
        if norm:
            proj[art].Scale(1/proj[art].Integral(0, proj[art].GetNbinsX()+1))
            proj[art].SetMaximum(0.17)
        proj[art].Rebin(2)
        proj[art].SetLineWidth(3)
        style(proj[art])
        if norm: texs[art] = ROOT.TLatex(0.25, 0.9-art*0.06, "%iX triggers"      % (art))
        else:    texs[art] = ROOT.TLatex(0.21, 0.9-art*0.06, "%iX triggers (%i)" % (art, proj[art].Integral(0, proj[art].GetNbinsX()+1)))
        texs[art].SetTextColor(color(art-1))
        style(texs[art])
        accept15 = proj[art].Integral(proj[art].FindBin(-15), proj[art].FindBin(15)) / proj[art].Integral(0, proj[art].GetNbinsX()+1)
        a15s[art] = ROOT.TLatex(0.65, 0.9-art*0.06, "%4.1f%%" % (accept15*100))
        a15s[art].SetTextColor(color(art-1))
        style(a15s[art])
        accept99 = proj[art].Integral(proj[art].FindBin(-99), proj[art].FindBin(99)) / proj[art].Integral(0, proj[art].GetNbinsX()+1)
        a99s[art] = ROOT.TLatex(0.80, 0.9-art*0.06, "%4.1f%%" % (accept99*100))
        a99s[art].SetTextColor(color(art-1))
        style(a99s[art])

    legs = [ROOT.TLatex(0.63, 0.84, "#Delta#theta < 15"),
            ROOT.TLatex(0.78, 0.84, "#Delta#theta < 100"),
            ]
    for leg in legs:
        style(leg)

    lines = [ROOT.TLine(0.62, 0.83, 0.76, 0.83),
             ROOT.TLine(0.78, 0.83, 0.92, 0.83),
             ]
    for line in lines:
        line.SetLineWidth(2)
        style(line)
    for art in proj.values():
        art.SetMaximum(1.2*art.GetBinContent(art.GetMaximumBin()))
        art.GetXaxis().SetRangeUser(-290, 290)

    canv = ROOT.TCanvas("theta_resolution_lin", "theta_resolution_lin", 800, 800)
    canv.Draw()
    ROOT.gPad.SetLogy(0)

    for art in sorted(proj.values(), key=lambda h: h.GetMaximum(), reverse=True):
        art.Draw("histsame")
    for tex in reversed(sorted(texs.keys())): texs[tex].Draw()
    for acc in reversed(sorted(a15s.keys())): a15s[acc].Draw()
    for acc in reversed(sorted(a99s.keys())): a99s[acc].Draw()
    for leg in legs:                          leg.Draw()
    for lin in lines:                         lin.Draw()
    save(canv, pdf, "trigger_theta_resolution_vs_n_lin.pdf")

    canv = ROOT.TCanvas("theta_resolution_log", "theta_resolution_log", 800, 800)
    canv.Draw()
    ROOT.gPad.SetLogy(1)
    for art in proj: proj[art].SetMinimum(2e-4 if norm else 0.5)
    for art in proj: proj[art].SetMaximum(0.4  if norm else 3*proj[art].GetMaximum())
    for art in reversed(sorted(proj.keys())): proj[art].Draw("histsame")
    for tex in reversed(sorted(texs.keys())): texs[tex].Draw()
    for acc in reversed(sorted(a15s.keys())): a15s[acc].Draw()
    for acc in reversed(sorted(a99s.keys())): a99s[acc].Draw()
    for leg in legs:                          leg.Draw()
    for lin in lines:                         lin.Draw()
    save(canv, pdf, "trigger_theta_resolution_vs_n_log.pdf")

def trigger_efficiency_vs_angle(numer, denom, pdf):

    rootlogon()
    canv = ROOT.TCanvas("trigeff_angle", "trigeff_angle", 800, 800)
    canv.Draw()

    bins = [-35, -20, -18, -16, -14, -12, -10]
    bins = bins + range(-10, -5)
    bins = bins + [i/10.0 for i in range(-50, 40, 2)]
    bins = bins + [4, 5, 6, 7, 8, 35]
    bins = array.array("d", bins)
    numer = numer.Rebin(len(bins)-1, numer.GetName(), bins)
    denom = denom.Rebin(len(bins)-1, denom.GetName(), bins)

    line1 = ROOT.TLine(-4, 90, -4, 109)
    line2 = ROOT.TLine( 4, 90,  4, 109)
    lines = [line1, line2]
    for line in lines:
        line.SetLineStyle(7)
        line.SetLineColor(ROOT.kGray)
        
    ratio = copy.copy(numer)
    ratio.Divide(numer, denom, 1.0, 1.0, "B")
    ratio.Scale(100)
    ratio.SetMinimum(0)
    # ratio.SetMinimum(90)
    ratio.SetMaximum(109)
    ratio.GetYaxis().SetTitle("Efficiency [%]")
    ratio.SetFillColor(ROOT.kYellow)
    ratio.SetLineColor(ROOT.kBlack)
    ratio.SetLineWidth(2)
    style(ratio)
    ratio.GetYaxis().SetTitleOffset(1.65)
    ratio.GetXaxis().SetTitle("#theta_{MM} [degrees]")
    ratio.Draw("histsame")
    #for line in lines:
    #    line.Draw()

    save(canv, pdf, "trigger_efficiency_vs_angle.pdf")

def color(it):
    if it == 0: return ROOT.kBlack
    if it == 1: return ROOT.kBlue
    if it == 2: return ROOT.kRed
    if it == 3: return 210
    if it == 4: return ROOT.kOrange+1
    if it == 5: return ROOT.kViolet
    if it == 6: return ROOT.kPink+7
    if it == 7: return ROOT.kCyan-7
    if it == 8: return ROOT.kGray

def fatal(msg):
    import sys
    sys.exit("Fatal error %s" % (msg))

def rootlogon():
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetPadTopMargin(0.06)
    ROOT.gStyle.SetPadRightMargin(0.05)
    ROOT.gStyle.SetPadBottomMargin(0.12)
    ROOT.gStyle.SetPadLeftMargin(0.16)
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)
    ROOT.gStyle.SetPaintTextFormat(".2f")
    ROOT.gStyle.SetTextFont(42)
    ROOT.gStyle.SetFillColor(10)

def style(hist):
    if isinstance(hist, ROOT.TLegend):
        hist.SetBorderSize(0)
        hist.SetFillColor(0)
        hist.SetFillStyle(0)
        hist.SetTextSize(0.04)
        return
    if isinstance(hist, ROOT.TLatex):
        hist.SetTextSize(0.040)
        hist.SetTextFont(42)
        hist.SetNDC()
        return
    if isinstance(hist, ROOT.TLine):
        hist.SetNDC()
        hist.SetLineWidth(2)
        return
    size = 0.045
    hist.GetXaxis().SetTitleSize(size)
    hist.GetXaxis().SetLabelSize(size)
    hist.GetYaxis().SetTitleSize(size)
    hist.GetYaxis().SetLabelSize(size)
    if isinstance(hist, ROOT.TH2):
        hist.GetXaxis().SetTitleOffset(1.1)
        hist.GetYaxis().SetTitleOffset(1.8)
        hist.GetZaxis().SetTitleOffset(1.4)
        hist.GetZaxis().SetTitleSize(size)
        hist.GetZaxis().SetLabelSize(size)
    elif isinstance(hist, ROOT.TH1) or isinstance(hist, ROOT.TProfile):
        hist.GetXaxis().SetTitleOffset(1.1)
        hist.GetYaxis().SetTitleOffset(1.85)
    
def show_overflow(hist, show_underflow=True, show_overflow=True):
    """ h/t Josh """
    nbins          = hist.GetNbinsX()
    underflow      = hist.GetBinContent(   0   )
    underflowerror = hist.GetBinError  (   0   )
    overflow       = hist.GetBinContent(nbins+1)
    overflowerror  = hist.GetBinError  (nbins+1)
    firstbin       = hist.GetBinContent(   1   )
    firstbinerror  = hist.GetBinError  (   1   )
    lastbin        = hist.GetBinContent( nbins )
    lastbinerror   = hist.GetBinError  ( nbins )
    if show_underflow and underflow != 0:
        newcontent = underflow + firstbin
        if firstbin == 0 :
            newerror = underflowerror
        else:
            newerror = math.sqrt( underflowerror * underflowerror + firstbinerror * firstbinerror )
        hist.SetBinContent(1, newcontent)
        hist.SetBinError  (1, newerror)
        hist.SetBinContent(0, 0)
        hist.SetBinError  (0, 0)
    if show_overflow and overflow != 0:
        newcontent = overflow + lastbin
        if lastbin == 0 :
            newerror = overflowerror
        else:
            newerror = math.sqrt( overflowerror * overflowerror + lastbinerror * lastbinerror )
        hist.SetBinContent(nbins,   newcontent)
        hist.SetBinError  (nbins,   newerror)
        hist.SetBinContent(nbins+1, 0)
        hist.SetBinError  (nbins+1, 0)

def options():
    parser = argparse.ArgumentParser(usage=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-r", default=None, help="Run number")
    parser.add_argument("-a", default="",   help="Alignment used")
    return parser.parse_args()

def save(canvas, pdf, name):
    if bigpdf:
        canvas.Print(pdf, "pdf")
    else:
        dir = pdf.rstrip(".pdf")
        if not os.path.isdir(dir):
            os.makedirs(dir)
        canvas.SaveAs(os.path.join(dir, name))

def open_pdf(pdf):
    rootlogon()
    canv = ROOT.TCanvas("open", "open", 800, 800)
    canv.Draw()
    tex = ROOT.TLatex(0.5, 0.5, "Beginning of plots PDF")
    tex.SetNDC()
    tex.SetTextAlign(22)
    tex.Draw()
    canv.Print(pdf+"(", "pdf")

def divide_pdf(pdf):
    rootlogon()
    canv = ROOT.TCanvas("open", "open", 800, 800)
    canv.Draw()
    tex = ROOT.TLatex(0.5, 0.5, "End of MMFE8 plots. Start of trigger plots.")
    tex.SetNDC()
    tex.SetTextAlign(22)
    tex.Draw()
    canv.Print(pdf, "pdf")

def close_pdf(pdf):
    rootlogon()
    canv = ROOT.TCanvas("close", "close", 800, 800)
    canv.Draw()
    tex = ROOT.TLatex(0.5, 0.5, "End of plots PDF")
    tex.SetNDC()
    tex.SetTextAlign(22)
    tex.Draw()
    canv.Print(pdf+")", "pdf")


if __name__ == "__main__":
    main()
