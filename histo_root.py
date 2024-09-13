from __future__ import division
import ROOT
import math
from math import *
import sys
import csv
import numpy as np
import awkward as ak
import uproot as uproot
import matplotlib.pyplot as plt
ROOT.gROOT.SetBatch(True)

# defining trucation function                                                 
def truncate(number, decimals=0):                                             
    """                                                                       
    Returns a value truncated to a specific number of decimal places.         
    """                                                                       
    if not isinstance(decimals, int):                                         
        raise TypeError("decimal places must be an integer.")                 
    elif decimals < 0:                                                        
        raise ValueError("decimal places has to be 0 or more.")               
    elif decimals == 0:                                                       
        return math.trunc(number)                                             
    factor = 10.0 ** decimals                                                 
    return math.trunc(number * factor) / factor 

def find_index(vector, value):
    for index, item in enumerate(vector):
        if item == value:
            return index
    return 0


filename = 'histo.root'
tf=ROOT.TFile.Open(filename)
tree_c3=tf.Get('ticlDumper/trackstersMerged')
tree_lc=tf.Get('ticlDumper/clusters')
tree_rh=tf.Get('ticlDumper/rechits')
tree_tpe=tf.Get('ticlDumper/simtrackstersCP')
c3_hits = ROOT.TH1D('c3hits', 'c3_hits',150,0,150)
all_hits = ROOT.TH1D('all_hits', 'all_hits',100,0,20)
z_hits = ROOT.TH1D('z_hits', 'z_hits',100,0,700)
tpe = ROOT.TH1D('tpe', 'tpe',150,0,150)
reso1 = ROOT.TH1D('reso_hist', 'reso_hist',50,0,2)
all_rh_in_ts = []
counter=0
energy=0
rhInLC_idx_ev=None
rhInLC_idx_ev1=None
all_e=0

all_hists_sum=0
for it,t in enumerate(tree_c3):
    if it==285 or it ==446:
        continue
    all_e=0
    energy=0
    lcInTs_idx_ev=t.vertices_indexes
    tree_lc.GetEntry(it)
    rhInLC_idx_ev=tree_lc.rechits
    tree_rh.GetEntry(it)
    tree_tpe.GetEntry(it)
    tpe.Fill(tree_tpe.regressed_energy[0])
#    all_hits_sum=tree_rh.energy
    for i in range(len(tree_rh.noise)):
        if tree_rh.noise[i]<2.5 and tree_rh.position_z[i]>0:
            z_hits.Fill(tree_rh.energy[i])
            all_e+=tree_rh.energy[i]
    for i in range(len(tree_rh.energy)):
        all_hits.Fill(tree_rh.energy[i])
    # print("tpe:", tree_tpe.regressed_energy[0])
    # print("z+ rechits summed: ",all_e)
    # print("all rechits summed: ",all_hits_sum)
    for i in range(len(rhInLC_idx_ev)):
        for j in range(len(rhInLC_idx_ev[i])):
            rh_idx=find_index(tree_rh.ID,rhInLC_idx_ev[i][j])
            #            if tree_rh.position_z[rh_idx]<=0:
            #               print(tree_rh.position_z[rh_energy])
            c3_hits.Fill(tree_rh.energy[rh_idx])            
            energy+=tree_rh.energy[rh_idx]
            counter+=1
#    c3_hits.Fill(energy)            
    reso1.Fill(energy/tree_tpe.regressed_energy[0])
#    print("c3 energy: ",energy)
#    break
    
# tree_rh.GetEntry(5)
# all_e=0
# tree_tpe.GetEntry(5)
# print("true particle energy: ",tree_tpe.regressed_energy[0])
# print("number of rh in rechits: ",len(tree_rh.energy))
# cnt=0
# #for i in range(len(tree_rh.energy)):
# #    print(tree_rh.ID[i])
# #    print(tree_rh.energy[i])
#     # cnt+=1
#     # all_e+=tree_rh.energy[i]
#     # all_hits.Fill(tree_rh.energy[i])


all_hits.SetLineColor(ROOT.kRed)
c3_hits.SetLineColor(ROOT.kBlue)
z_hits.SetLineColor(ROOT.kGreen)
tpe.SetLineColor(ROOT.kBlack)
all_hits.SetLineWidth(2)
c3_hits.SetLineWidth(2)
z_hits.SetLineWidth(2)
tpe.SetLineWidth(2)
reso1.SetLineWidth(2)
reso1.GetXaxis().SetRangeUser(0,3)
all_hits.GetXaxis().SetTitle("Energy [ GeV]")
all_hits.GetYaxis().SetTitle("Count")
reso1.GetXaxis().SetTitle("summed E/True Particle E [ GeV]")
reso1.GetYaxis().SetTitle("Count")
leg1=ROOT.TLegend(0.50,.70,.68,.90)
leg1.SetBorderSize(0)
leg1.SetFillColor(0)
leg1.SetFillStyle(0)
leg1.SetTextSize(0.035)
leg1.AddEntry(tpe,'True Particle Energy','L')
leg1.AddEntry(c3_hits,'summed Clue3d hits','L')
leg1.AddEntry(all_hits,'all summed rechits','L')
leg1.AddEntry(z_hits,'+z hits only','L')


c3=ROOT.TCanvas("1can1va2345s11","can11va1324521s1",800,600)
all_hits.Draw("HIST")
z_hits.Draw("SAME")
c3_hits.Draw("SAME")
#tpe.Draw("SAME")
leg1.Draw("SAME")
c3.SetLogy()
c3.Draw()
c3.SaveAs("hits_distribution_summed.png")


c4=ROOT.TCanvas("canvas","canvas",800,600)
reso1.Draw("HIST")
c4.SetLogy()
c4.Draw()
c4.SaveAs("reso.png")
