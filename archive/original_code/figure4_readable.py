from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
HERE=Path(__file__).resolve().parent
out=HERE/'figures';out.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none'})
for name in ['Baseline','BBB_pass','ADME']:
    a=np.load(HERE.parent/f'recalculation/repaired/global_models/seed42/{name}/shap_values.npz');sv=a['shap'][:,:-1];x=a['values'];features=a['features'];ranks=np.argsort(np.abs(sv).mean(0))[-20:]
    fig,ax=plt.subplots(figsize=(7.1,8),layout='constrained');rng=np.random.default_rng(42)
    for row,j in enumerate(ranks):
        lo,hi=np.quantile(x[:,j],[.02,.98]);color=np.clip((x[:,j]-lo)/max(hi-lo,1e-12),0,1)
        ax.scatter(sv[:,j],row+rng.uniform(-.24,.24,len(x)),c=color,cmap='coolwarm',vmin=0,vmax=1,s=8,alpha=.65,rasterized=True)
    ax.axvline(0,color='grey',lw=.7);ax.set(yticks=range(20),yticklabels=features[ranks],xlabel='SHAP contribution',title=name.replace('_',' '));ax.tick_params(axis='y',labelsize=11)
    cb=fig.colorbar(plt.cm.ScalarMappable(norm=plt.Normalize(0,1),cmap='coolwarm'),ax=ax,shrink=.5,pad=.025);cb.set_ticks([0,1],labels=['Low','High']);cb.set_label('Feature value')
    for ext in ['png','svg']:fig.savefig(out/f'figure4_{name}.{ext}',dpi=300)
    plt.close(fig)
