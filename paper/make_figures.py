#!/usr/bin/env python3
"""Figures: two deep models straddle the simple interpretable baselines on native proteins."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":9.5,
    "axes.spines.top":False,"axes.spines.right":False,
    "axes.edgecolor":"#444444","axes.linewidth":0.8})
DEEP="#D95F0E"; DEEP2="#E8943F"; SIMPLE="#2C7FB8"; SIMPLE2="#5AB4AC"; CALIB="#7B5EA7"

def labels_above(ax, bars, vals, his, fs=8.2):
    for bar,v,h in zip(bars,vals,his):
        ax.text(bar.get_x()+bar.get_width()/2, h+0.009, f"{v:.3f}",
                ha="center", va="bottom", fontsize=fs, fontweight="bold")

def bracket(ax, x1, x2, y, text, fs=7.8, color="#222"):
    ax.plot([x1,x1,x2,x2],[y,y+0.006,y+0.006,y], c=color, lw=1.0)
    ax.text((x1+x2)/2, y+0.010, text, ha="center", va="bottom", fontsize=fs, color=color)

fig, axes = plt.subplots(1, 2, figsize=(9.0, 4.6), gridspec_kw={"width_ratios":[4,3.5]})

# ---- eSOL: NetSolP / SWI / naive / RP3Net ----
ax=axes[0]
lab=["NetSolP\n(deep)","SWI\n(simple)","Heuristic\n(naive)","RP3Net\n(deep)"]
vals=[0.792,0.745,0.725,0.709]; lo=[0.771,0.723,0.702,0.686]; hi=[0.812,0.767,0.746,0.731]
cols=[DEEP,SIMPLE2,SIMPLE,DEEP]
err=[[v-l for v,l in zip(vals,lo)],[h-v for v,h in zip(vals,hi)]]
bars=ax.bar(range(4),vals,0.62,yerr=err,capsize=3.5,color=cols,error_kw=dict(ecolor="#222",lw=1.0))
labels_above(ax,bars,vals,hi)
ax.axhline(0.5,ls=":",c="#999",lw=1)
bracket(ax,0,1,0.875,"p < 0.001")
bracket(ax,1,3,0.815,"SWI > RP3Net (p = 0.006)")
ax.set_xticks(range(4)); ax.set_xticklabels(lab,fontsize=8.2)
ax.set_ylabel("AUROC (95% CI)"); ax.set_ylim(0.45,0.94)
ax.set_title("eSOL cytoplasmic (native, n = 2,154)",fontsize=10.5,fontweight="bold")

# ---- SoluProt: NetSolP / RP3Net / SWI / naive / ProteinSol ----
ax=axes[1]
l2=["NetSolP","RP3Net","SWI","Heuristic\n(naive)","ProteinSol\n(eSOL-calib.)"]
v2=[0.633,0.620,0.598,0.581,0.542]; lo2=[0.612,0.599,0.578,0.562,0.521]; hi2=[0.652,0.640,0.617,0.601,0.562]
c2=[DEEP,DEEP2,SIMPLE2,SIMPLE,CALIB]
e2=[[v-l for v,l in zip(v2,lo2)],[h-v for v,h in zip(v2,hi2)]]
b2=ax.bar(range(5),v2,0.64,yerr=e2,capsize=2.8,color=c2,error_kw=dict(ecolor="#222",lw=1.0))
labels_above(ax,b2,v2,hi2,fs=7.6)
ax.axhline(0.5,ls=":",c="#999",lw=1)
bracket(ax,0,1,0.690,"n.s.",fs=7.4,color="#666")
ax.set_xticks(range(5)); ax.set_xticklabels(l2,fontsize=7.2)
ax.set_ylim(0.45,0.94); ax.set_yticklabels([])
ax.set_title("SoluProt (heterolog., leak-free, n = 3,060)",fontsize=10,fontweight="bold")

fig.suptitle("Two SOTA deep models disagree by more than the gap to a simple baseline on native proteins",
             fontsize=10,fontweight="bold",y=1.005)
fig.tight_layout(rect=[0,0,1,0.99])
fig.savefig("fig1_auroc.png",dpi=200,bbox_inches="tight"); plt.close(fig)

# ---- Fig 2: threshold robustness, all four, leak-free ----
fig,ax=plt.subplots(figsize=(6.6,3.6))
thr=[60,70,80]
ax.plot(thr,[0.781,0.792,0.799],"-o",c=DEEP,lw=2,ms=6,label="NetSolP (deep)")
ax.plot(thr,[0.747,0.745,0.750],"-D",c=SIMPLE2,lw=2,ms=6,label="SWI (simple)")
ax.plot(thr,[0.719,0.725,0.727],"-^",c=SIMPLE,lw=2,ms=6,label="Heuristic (naive)")
ax.plot(thr,[0.716,0.709,0.711],"-s",c=DEEP2,lw=2,ms=6,label="RP3Net (deep)")
ax.set_xticks(thr); ax.set_xlabel("eSOL solubility threshold (%)")
ax.set_ylabel("AUROC (cytoplasmic)"); ax.set_ylim(0.66,0.82)
ax.legend(frameon=False,fontsize=8.3,loc="center left",bbox_to_anchor=(1.01,0.5),ncol=1)
ax.set_title("Threshold robustness (leak-free, n = 2,154)",fontsize=10.5,fontweight="bold")
fig.tight_layout(); fig.savefig("fig2_threshold.png",dpi=200,bbox_inches="tight"); plt.close(fig)
print("figures rebuilt")
