from paths import ROOT, SOURCE, PROCESSED, RESULTS, REPAIR_OUTPUT
"""Updated global comparison and SHAP figures, using saved repaired-model predictions."""
import os,json,pickle
from pathlib import Path
import numpy as np,pandas as pd
from catboost import CatBoostRegressor,Pool
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
HERE=Path(__file__).resolve().parent
os.environ.setdefault('REANALYSIS_DATA',str(PROCESSED/'audited_parent_matched.csv.gz'));os.environ.setdefault('REANALYSIS_OUT',str(RESULTS))
from analysis import OUT,load,PROTEINS
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none'})
figdir=OUT/'figures';figdir.mkdir(exist_ok=True);td=OUT/'tables';td.mkdir(exist_ok=True)
def save(fig,name):
    for ext in ['png','svg']:fig.savefig(figdir/(name+'.'+ext),dpi=220,bbox_inches='tight')
    plt.close(fig)

def main():
    df=load().set_index('ligand_id');out=OUT/'global_models';allgroups=[]
    for name in ['Baseline','PCA','ADME','Plain','BBB pass']:
        rd=out/'seed42'/name.replace(' ','_');features=pd.read_csv(rd/'features.csv').feature.tolist();med=pd.read_csv(rd/'imputation_medians.csv',index_col=0).iloc[:,0]
        ids=pd.read_csv(rd/'shap_sample.csv').ligand_id;raw=df.loc[ids,med.index].apply(pd.to_numeric,errors='coerce').fillna(med)
        if name=='PCA':
            with (rd/'pca.pkl').open('rb') as f:pca=pickle.load(f)
            pc=pca.transform(raw[PROTEINS]);raw=raw.drop(columns=PROTEINS)
            for j in range(3):raw[f'PC{j+1}']=pc[:,j]
        x=raw[features];model=CatBoostRegressor();model.load_model(str(rd/'model.cbm'));sv=model.get_feature_importance(Pool(x),type='ShapValues',thread_count=2)
        assert np.allclose(sv[:,:-1].sum(axis=1)+sv[:,-1],model.predict(x),atol=1e-8)
        np.savez_compressed(rd/'shap_values.npz',ligand_id=ids.to_numpy(),shap=sv,features=np.array(features),values=x.to_numpy())
        gr=np.array(['Fingerprints' if c.startswith('FP_') else 'PhysChem' if c in ['MW, g/mol','MW','logP'] else 'Proteins' if c in PROTEINS or c.startswith('PC') else 'ADME' for c in features])
        for group in ['Fingerprints','Proteins','PhysChem','ADME']:
            ix=np.flatnonzero(gr==group);value=float(np.abs(sv[:,ix].sum(axis=1)).mean()) if len(ix) else np.nan
            allgroups.append({'model':name,'group':group,'mean_abs_grouped_shap':value,'n_features':len(ix),'n_sample':len(x)})
        if name not in ['Baseline','BBB pass','ADME']:continue
        ranks=np.argsort(np.abs(sv[:,:-1]).mean(axis=0))[-20:];fig,ax=plt.subplots(figsize=(9,7),layout='constrained');rng=np.random.default_rng(42)
        for row,j in enumerate(ranks):
            vals=x.iloc[:,j].to_numpy();lo,hi=np.quantile(vals,[.02,.98]);color=np.clip((vals-lo)/max(hi-lo,1e-12),0,1)
            ax.scatter(sv[:,j],row+rng.uniform(-.25,.25,len(x)),c=color,cmap='coolwarm',vmin=0,vmax=1,s=7,alpha=.65,rasterized=True)
        ax.axvline(0,color='grey',lw=.6);ax.set(yticks=range(20),yticklabels=[features[j] for j in ranks],xlabel='SHAP contribution to predicted -log10 LD50',title=f'{name}: 512 held-out molecules, split seed 42')
        cbar=fig.colorbar(plt.cm.ScalarMappable(norm=plt.Normalize(0,1),cmap='coolwarm'),ax=ax,shrink=.7);cbar.set_ticks([0,1],labels=['Low','High']);cbar.set_label('Feature value (scaled within each feature)')
        save(fig,'figure4_'+name.replace(' ','_'))
    g=pd.DataFrame(allgroups);g.to_csv(td/'figure3_grouped_shap.csv',index=False)
    order=['Baseline','Plain','BBB pass','ADME'];roworder=['Fingerprints','Proteins','PhysChem','ADME'];z=g.pivot(index='group',columns='model',values='mean_abs_grouped_shap').loc[roworder,order]
    fig,ax=plt.subplots(figsize=(8.2,5),layout='constrained');cmap=plt.get_cmap('Blues').copy();cmap.set_bad('#eeeeee');im=ax.imshow(z,cmap=cmap,vmin=0)
    for i in range(4):
        for j in range(4):v=z.iloc[i,j];ax.text(j,i,'N/A' if pd.isna(v) else f'{v:.3f}',ha='center',va='center',color='black')
    ax.set(xticks=range(4),xticklabels=order,yticks=range(4),yticklabels=roworder,title='Grouped SHAP on repaired data, split seed 42');fig.colorbar(im,ax=ax,label='Mean |sum of SHAP within group|',shrink=.8);save(fig,'figure3_revised')
    metrics=pd.read_csv(out/'metrics.csv')
    for ix,row in metrics.iterrows():
        p=pd.read_csv(out/f'seed{row.seed}'/row.model.replace(' ','_')/'test_predictions.csv');yt=p.y.to_numpy();yp=p.prediction.to_numpy()
        metrics.loc[ix,'mse']=np.mean((yt-yp)**2);metrics.loc[ix,'ccc']=2*np.mean((yt-yt.mean())*(yp-yp.mean()))/(yt.var()+yp.var()+(yt.mean()-yp.mean())**2)
    metrics.to_csv(out/'metrics.csv',index=False);metrics[metrics.seed.eq(42)].to_csv(td/'table1_primary_global_models.csv',index=False)
    metrics.groupby('model')[['rmse','mae','mse','r2','spearman','ccc','coverage_0.9','width_0.9','interval_score_90']].agg(['mean','std','min','max']).to_csv(out/'five_split_summary.csv')
    order=['Baseline','PCA','ADME','Plain','BBB pass'];fig,axs=plt.subplots(1,2,figsize=(11.5,4.8),layout='constrained')
    for ax,metric in zip(axs,['rmse','mae']):
        for seed in range(42,47):
            r=metrics[metrics.seed.eq(seed)].set_index('model').loc[order];ax.plot(range(5),r[metric],marker='o',lw=1,alpha=.7,label=str(seed))
        ax.set(xticks=range(5),xticklabels=['Baseline','PCA','ADME','Plain','BBB subset'],ylabel=metric.upper())
    axs[0].set_title('Point prediction across five structural partitions');axs[1].legend(title='Split seed',frameon=False);save(fig,'figure1_revised')
    fig,axs=plt.subplots(1,2,figsize=(11.5,4.8),layout='constrained')
    for j,level in enumerate([.8,.9,.95]):
        r=metrics.groupby('model').mean(numeric_only=True).loc[order];axs[0].plot(range(5),r[f'width_{level}'],marker='o',label=f'{level:.0%}');axs[1].plot(range(5),r[f'coverage_{level}']*100,marker='o',label=f'{level:.0%}');axs[1].axhline(level*100,ls=':',lw=.7,color='grey')
    for ax in axs:ax.set_xticks(range(5),['Baseline','PCA','ADME','Plain','BBB subset']);ax.legend(title='Nominal level',frameon=False)
    axs[0].set(ylabel='Mean full interval width',title='Global residual intervals');axs[1].set(ylabel='Observed coverage (%)',title='Coverage must accompany interval width');save(fig,'figure2_revised')
    print('GLOBAL FIGURES COMPLETE: SHAP additivity checked for five models.',flush=True)
if __name__=='__main__':main()
