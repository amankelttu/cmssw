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
ROOT.gStyle.SetOptStat(0000000)

# defining trucation function
def function1(graph,it,max_val=0):
    colors = [ROOT.kBlack, ROOT.kCyan,ROOT.kBlue,ROOT.kRed, ROOT.kGreen, ROOT.kMagenta, ROOT.kOrange, ROOT.kYellow, ROOT.kSpring]
    #    max_val=int(max_val)+100
    #graph.GetYaxis().SetRangeUser(0,max_val)
    graph.SetLineColor(colors[it])
    graph.SetLineWidth(2)

def find_max(g1,g2,g3,g4,g5,g6,g7):
    input=[]
    input.append(max(g1))
    input.append(max(g2))
    input.append(max(g3))
    input.append(max(g4))
    input.append(max(g5))
    input.append(max(g6))
    input.append(max(g7))
    return max(input)
def find_min_index(vector):
    # Initialize variables to store the minimum value and its index
    if not vector.size():  # Check if the vector is empty
        return -1  # Return -1 if the vector is empty
    min_value = float('inf')
    min_index = -1
    
    # Iterate over the vector to find the minimum value
    for i in range(vector.size()):
        if vector[i] < min_value:
            min_value = vector[i]
            min_index = i

    return min_index

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

def overlay_ratio_plot(h1,h2):
    # Create a canvas
    iter=sys.argv[2]
    c = ROOT.TCanvas("c", "Overlay + Ratio", 800, 800)    
    # Create two pads, one for the histograms and one for the ratio
    pad1 = ROOT.TPad("pad1", "Overlay", 0, 0.3, 1, 1.0)
    pad1.SetBottomMargin(0)  # Upper and lower margin
    pad1.SetLogy(1)
    pad1.Draw()
    pad1.cd()  # Move to this pad
    # Draw the first histogram
    h1.Draw("HIST")
    # Draw the second histogram on the same canvas
    h2.Draw("HIST SAME")
    # Now go back to the main canvas to create the ratio pad
    c.cd()
    pad2 = ROOT.TPad("pad2", "Ratio", 0, 0.05, 1, 0.3)
    pad2.SetTopMargin(0)
    pad2.SetBottomMargin(0.2)
    pad2.Draw()
    pad2.cd()  # Move to this pad

    # Create the ratio histogram
    h_ratio = h1.Clone("h_ratio")
    h_ratio.SetLineColor(ROOT.kBlack)
    h_ratio.SetTitle("")
    h_ratio.Divide(h2)  # Ratio h1 / h2

    # Adjust y-axis range for better visibility
    h_ratio.SetMinimum(0.5)
    h_ratio.SetMaximum(1.5)

    # Draw the ratio plot
    h_ratio.Draw("HIST")

    # Add some labels to the ratio plot
    h_ratio.GetYaxis().SetTitle("Ratio h1/h2")
    h_ratio.GetXaxis().SetTitle("X-axis title")
    h_ratio.GetYaxis().SetNdivisions(505)
    h_ratio.GetYaxis().SetTitleSize(20)
    h_ratio.GetYaxis().SetTitleFont(43)
    h_ratio.GetYaxis().SetTitleOffset(1.55)
    h_ratio.GetXaxis().SetTitleSize(20)
    h_ratio.GetXaxis().SetTitleFont(43)
    h_ratio.GetXaxis().SetTitleOffset(3.0)

    # Go back to the main pad for further work, if needed
    c.cd()
    # Save the canvas as an image or display it
    c.SaveAs(f"plots/overlay_ratio_plot_{iter}.png")
    c.Draw()

def fill_th2d_with_rz_energy(data,event):
    iter=sys.argv[2]
    # Create a TH2D histogram with (r, z) axes
    h2 = ROOT.TH2D("h2", "r vs z with energy", 300, 0, 300, 300, 300, 600)  # Binning and ranges can be adjusted
    # Loop over the data, calculate r, and fill the histogram
    for x, y, z, energy in data:
        r = math.sqrt(x**2 + y**2)  # Calculate radial distance r
        h2.Fill(r, z, energy)  # Fill the histogram with r, z and weight by energy
    # Create a canvas to draw the histogram
    c = ROOT.TCanvas("c", "r vs z with energy", 800, 600)
    h2.Draw("COLZ")  # Use COLZ option to show color mapping of the bin contents

    # Save the canvas as an image or display it
    c.SaveAs(f"plots/rz_energy_plot_{iter}.png")
    c.Draw()

iter=sys.argv[2]
print(iter)
filename = sys.argv[1]
tf=ROOT.TFile.Open(filename)

tree_c3=tf.Get('ticlDumper/trackstersTiclCandidate')
tree_lc=tf.Get('ticlDumper/clusters')
tree_rh=tf.Get('ticlDumper/rechits')
tree_tpe=tf.Get('ticlDumper/simtrackstersCP')
tree_can=tf.Get('ticlDumper/candidates')
assEv=tf.Get('ticlDumper/associations')

c3_hits = ROOT.TH1D('c3hits', 'c3_hits',200,0,250)
all_hits = ROOT.TH1D('all_hits', 'all_hits',200,0,250)
noise_hits = ROOT.TH1D('all123_hits', 'all_hi12ts',200,0,250)
z_hits = ROOT.TH1D('all_hi123ts', 'all_123hits',200,0,250)
tpe = ROOT.TH1D('tpe', 'tpe',200,0,250)
ind_rechit_c3 = ROOT.TH1D('individual hitsc3', 'individual hitsc3',200,0,10)
num_hits_c3 = ROOT.TH1D('#of hitsc3', '#of hitsc3',200,0,2300)
ind_rechit_rh = ROOT.TH1D('individual hits', 'individual hits',200,0,10)
num_hits_rh = ROOT.TH1D('#of hits', '#of hits',200,0,2300)
reso = ROOT.TH1D('reso_hist', 'reso_hist',50,0,3)
reso1 = ROOT.TH1D('reso_hist1', 'reso_hist1',50,0,3)
reso2 = ROOT.TH1D('reso_hist2', 'reso_hist2',50,0,3)

all_rh_in_ts = []
energy=0
rhInLC_idx_ev=None
rhInLC_idx_ev1=None
all_e=0
energy1=0
all_hists_sum=0
event_scale=0
idt=tree_c3.GetEntries()
ratio_c3=np.zeros(idt)
ratio_rh=np.zeros(idt)
ratio_cand=np.zeros(idt)
ratio_raw=np.zeros(idt)
data=[]

for it,t in enumerate(tree_c3):
    print(it)
    if it==285 or it ==446 or it==100:
        break
        continue 
    all_z=0
    energy=0
    hit_count_c3=0
    hit_count_rh=0
    #set trees to match
    assEv.GetEntry(it)
    tree_lc.GetEntry(it)
    tree_rh.GetEntry(it)
    tree_tpe.GetEntry(it)
    trpe=tree_tpe.regressed_energy[0]
    if trpe<=1: #incase an event fails to reco
        event_scale+=1
        continue
    #lines needed to  find their regressed energy outputs
    simCand_idx=0
    simToReco = assEv.ticlCandidate_simToReco_CP[simCand_idx]
    sharedE = assEv.ticlCandidate_simToReco_CP_sharedE[simCand_idx]
    score = assEv.ticlCandidate_simToReco_CP_score[simCand_idx]
    index=find_min_index(score)
    recoenergy=-1
    rawenergy=-1
    if index!=-1:
        tid = simToReco[index]
        if sharedE[index]/trpe>0.50:
            recoenergy = t.raw_energy[tid]
            print("raw_energy: ",recoenergy)
            rawenergy = t.regressed_energy[tid]
            print("regressed energy: ",rawenergy)
    lcInTs_idx_ev=t.vertices_indexes
    rhInLC_idx_ev=tree_lc.rechits
    for i in range(len(tree_rh.energy)):
        if tree_rh.noise[i]>3:
            hit_count_rh+=1
            ind_rechit_rh.Fill(tree_rh.energy[i])
            all_z+=tree_rh.energy[i]    
    for i in range(len(rhInLC_idx_ev)):
        for j in range(len(rhInLC_idx_ev[i])):
            rh_idx=find_index(tree_rh.ID,rhInLC_idx_ev[i][j])
            if tree_rh.noise[rh_idx]>3:
                hit_count_c3+=1
                ind_rechit_c3.Fill(tree_rh.energy[rh_idx])
                energy+=tree_rh.energy[rh_idx]
                data.append((tree_rh.position_x[rh_idx],tree_rh.position_y[rh_idx],tree_rh.position_z[rh_idx],tree_rh.energy[rh_idx]))    
    ratio_c3[it]=energy/trpe
    ratio_rh[it]=all_z/trpe
    ratio_cand[it]=recoenergy/trpe
    ratio_raw[it]=rawenergy/trpe
    c3_hits.Fill(energy)
    tpe.Fill(trpe)
    noise_hits.Fill(recoenergy)
    all_hits.Fill(all_z)
    reso.Fill(energy/trpe)
    reso1.Fill(all_z/trpe)
    reso2.Fill(recoenergy/trpe)
    num_hits_rh.Fill(hit_count_rh)
    num_hits_c3.Fill(hit_count_c3)


# Run the function

std_c3=np.std(ratio_c3)
mean_c3=np.mean(ratio_c3)
resolution_c3=std_c3/mean_c3

std_rh=np.std(ratio_rh)
mean_rh=np.mean(ratio_rh)
resolution_rh=std_rh/mean_rh

std_cand=np.std(ratio_cand)
mean_cand=np.mean(ratio_cand)
resolution_cand=std_cand/mean_cand

std_raw=np.std(ratio_raw)
mean_raw=np.mean(ratio_raw)
resolution_raw=std_raw/mean_raw

plt.figure(2)
fig1,ax1=plt.subplots()
ax1.grid(True)
ax1.xaxis.set_tick_params(labelsize=10)
ax1.yaxis.set_tick_params(labelsize=10)
ax1.scatter(149,resolution_c3,label='C3 resolution',marker='x',color='blue')
ax1.scatter(151,resolution_rh,label='Rechit resoltion',color='red')
ax1.scatter(152,resolution_cand,label='TICL cand regressed energy resolution',marker='x',color='purple')
ax1.scatter(153,resolution_raw,label='TICL cand raw energy resolution',marker='x',color='green')
ax1.set_xlabel('Ebeam [GeV]',fontsize=12)
ax1.set_ylabel(r'$\frac{RMS}{\langle E_{reco} \rangle}$',fontsize=12)
ax1.set_title('Resoluion {iter}',fontsize=14)
ax1.set_ylim(0,1)
pos = ax1.get_position()
ax1.set_position([pos.x0, pos.y0, pos.width * 0.9, pos.height])
ax1.legend(loc='upper right',fontsize=10)
plt.tight_layout()
plt.savefig(f'plots/resoltuion_{iter}.png',dpi=300)
plt.close()                       

function1(all_hits,3)
function1(noise_hits,5)
function1(c3_hits,1)
function1(tpe,0)

function1(reso,1)
function1(reso1,3)
function1(reso2,5)
function1(ind_rechit_c3,1)
function1(ind_rechit_rh,3)

function1(num_hits_c3,1)
function1(num_hits_rh,3)

# ind_rechit_c3.Scale(1/(1000-event_scale))
# ind_rechit_rh.Scale(1/(1000-event_scale))
# num_hits_c3.Scale(1/(1000-event_scale))
# num_hits_rh.Scale(1/(1000-event_scale))

overlay_ratio_plot(ind_rechit_c3,ind_rechit_rh)
fill_th2d_with_rz_energy(data,1000-event_scale)


num_hits_rh.GetXaxis().SetTitle("# of Hits")
num_hits_rh.GetYaxis().SetTitle("Count")
ind_rechit_rh.GetXaxis().SetTitle("Individual rechits[GeV]")
ind_rechit_rh.GetYaxis().SetTitle("Count")
all_hits.GetXaxis().SetTitle("Energy [ GeV]")
all_hits.GetYaxis().SetTitle("Count")
reso.GetXaxis().SetTitle("summed E/True Particle E [ GeV]")
reso.GetYaxis().SetTitle("Count")

leg1=ROOT.TLegend(0.70,.70,.98,.90)
leg1.SetBorderSize(0)
leg1.SetFillColor(0)
leg1.SetFillStyle(0)
leg1.SetTextSize(0.035)
leg1.AddEntry(tpe,'True Particle Energy','L')
leg1.AddEntry(c3_hits,'summed Clue3d hits','L')
leg1.AddEntry(all_hits,'summed rechits Noise, -z removed','L')
leg1.AddEntry(noise_hits,'TICL candidate','L')

#all_hits.GetYaxis().SetRangeUser(0,25)
c3=ROOT.TCanvas("1can1va2345s11","can11va1324521s1",800,600)
all_hits.Draw("HIST")
noise_hits.Draw("SAMEHIST")
c3_hits.Draw("SAMEHIST")
tpe.Draw("SAMEHIST")
leg1.Draw("SAME")
c3.SetLogy()
c3.Draw()
c3.SaveAs(f"plots/hits_{iter}.png")


c4=ROOT.TCanvas("canvas","canvas",800,600)
reso.Draw("HIST")
reso1.Draw("SAMEHIST")
reso2.Draw("SAMEHIST")
leg1.Draw("SAME")
c4.SetLogy()
c4.Draw()
c4.SaveAs(f"plots/reso_{iter}.png")


leg2=ROOT.TLegend(0.70,.70,.98,.90)
leg2.SetBorderSize(0)
leg2.SetFillColor(0)
leg2.SetFillStyle(0)
leg2.SetTextSize(0.035)
leg2.AddEntry(c3_hits,'Clue3d Rechits','L')
leg2.AddEntry(all_hits,'all Rechits','L')


c5=ROOT.TCanvas("c5anvas","can5vas",800,600)
num_hits_rh.Draw("HIST")
num_hits_c3.Draw("SAMEHIST")

leg2.Draw("SAME")
c5.SetLogy()
c5.Draw()
c5.SaveAs(f"plots/num_hits_{iter}.png")


# result = num_hits_rh.Clone("result")
# result.Add(num_hits_c3, -1.0);

# c6=ROOT.TCanvas("c5an6vas","can5va6s",800,600)
# result.Draw("HIST")
# c6.SetLogy()
# c6.Draw()
# c6.SaveAs(f"plots/diff_of_hits_{iter}.png")

