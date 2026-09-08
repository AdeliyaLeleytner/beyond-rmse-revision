from paths import ROOT, SOURCE, PROCESSED, RESULTS, REPAIR_OUTPUT
from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from analysis import OUT,REPO,load

plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none'})
FIG=OUT/'figures';FIG.mkdir(exist_ok=True)
TABLE=OUT/'tables';TABLE.mkdir(exist_ok=True)
BLUE='#356c9a';ORANGE='#bd6b37';GREY='#747e86'

def save(fig,name):
    fig.savefig(FIG/(name+'.png'),dpi=220,bbox_inches='tight')
    fig.savefig(FIG/(name+'.svg'),bbox_inches='tight')
    plt.close(fig)

def main():
    overview=pd.read_csv(OUT/'overview.csv');all_metrics=[];all_quad=[];selected=[]
    for p in sorted((OUT/'runs').glob('*/complete.json')):
        cfg=json.loads(p.read_text());m=pd.read_csv(p.parent/'interval_metrics.csv');m.insert(0,'run',cfg['run']);all_metrics.append(m)
        q=pd.read_csv(p.parent/'quadrants.csv');q.insert(0,'run',cfg['run']);all_quad.append(q)
        selected.append({'run':cfg['run'],'selected_markers':';'.join(cfg['selected_markers'])})
    metrics=pd.concat(all_metrics,ignore_index=True);quad=pd.concat(all_quad,ignore_index=True)
    metrics.to_csv(TABLE/'all_interval_metrics.csv',index=False);quad.to_csv(TABLE/'all_quadrants.csv',index=False);pd.DataFrame(selected).to_csv(TABLE/'selected_markers.csv',index=False)
    primary=overview[overview.run.str.startswith('butina04_')].sort_values('run');random=overview[overview.run.str.startswith('random_size_')].sort_values('run')
    summary=[]
    for name,group in [('cluster_disjoint',primary),('molecule_random_matched_sizes',random)]:
        for col in ['rmse','global_coverage','joint_coverage','global_width','joint_width','global_score','joint_score']:
            summary.append({'regime':name,'metric':col,'mean':float(group[col].mean()),'sd_between_splits':float(group[col].std()),'min':float(group[col].min()),'max':float(group[col].max())})
    pd.DataFrame(summary).to_csv(TABLE/'five_split_summary.csv',index=False)
    for source,dest in [('runs/butina04_seed42/interval_metrics.csv','table2_primary_intervals.csv'),('runs/butina04_seed42/quadrants.csv','table3_primary_quadrants.csv')]:
        d=pd.read_csv(OUT/source);d.to_csv(TABLE/dest,index=False)
        (TABLE/Path(dest).with_suffix('.tex')).write_text('% Numerical table; see CSV for full precision.\n'+'\n'.join(' & '.join(str(v).replace('_',r'\_') for v in row)+r' \\' for row in [list(d.columns)]+d.round(4).values.tolist()))
    fig,axs=plt.subplots(2,2,figsize=(11.6,8.2),layout='constrained')
    ax=axs[0,0]
    for i in range(5):ax.plot([0,1],[random.iloc[i].rmse,primary.iloc[i].rmse],color=GREY,alpha=.5,lw=1)
    ax.scatter(np.zeros(5),random.rmse,color=BLUE,s=35);ax.scatter(np.ones(5),primary.rmse,color=ORANGE,s=35)
    ax.set(xticks=[0,1],xticklabels=['Molecule-random','Separate clusters'],ylabel='Test RMSE',title='a  Point prediction across splits',xlim=(-.35,1.35))
    ax=axs[0,1]
    for j,method in enumerate(['global','joint']):
        m=metrics[(metrics.run.str.startswith('butina04_'))&metrics.method.eq(method)].sort_values('run');y=m.coverage.to_numpy()*100;x=np.arange(5)+(j-.5)*.12
        ax.errorbar(x,y,yerr=np.vstack([y-m.coverage_ci_low.to_numpy()*100,m.coverage_ci_high.to_numpy()*100-y]),fmt='o',color=[GREY,BLUE][j],capsize=3,label=['Global','Local biological'][j])
    ax.axhline(90,color='black',ls='--',lw=1);ax.set(xticks=np.arange(5),xticklabels=range(42,47),xlabel='Split seed',ylabel='Test coverage (%)',title='b  Coverage varies across splits');ax.legend(frameon=False,fontsize=9)
    ax=axs[1,0]
    for i in range(5):ax.plot([0,1],[primary.iloc[i].global_width,primary.iloc[i].joint_width],color=GREY,alpha=.5)
    ax.scatter(np.zeros(5),primary.global_width,color=GREY);ax.scatter(np.ones(5),primary.joint_width,color=BLUE)
    ax.set(xticks=[0,1],xticklabels=['Global','Local biological'],ylabel='Mean full interval width',title='c  Interval width across splits',xlim=(-.35,1.35))
    ax=axs[1,1];ab=pd.read_csv(OUT/'analogue_ablation/summary.csv');full=ab.loc[ab.arm.eq('all_training'),'rmse'].iloc[0];removed=ab.loc[ab.arm.eq('remove_test_cluster_members'),'rmse'].iloc[0];rand=ab[ab.arm.str.startswith('random_')].rmse
    ax.scatter([0],[full],color=GREY,s=55);ax.scatter(np.ones(5)+np.linspace(-.07,.07,5),rand,color=BLUE,s=30);ax.scatter([2],[removed],color=ORANGE,s=55)
    ax.set(xticks=[0,1,2],xticklabels=['Full train','Random removal\n(same size)','Remove test-\ncluster members'],ylabel='RMSE on the same test molecules',title='d  Fixed-test analogue-removal control',xlim=(-.4,2.4))
    save(fig,'recalculation_comparison')

    # Replacement Fig.5: all four outcomes smoothed on held-out test molecules only.
    name='butina04_seed42';pred=pd.read_csv(OUT/'runs'/name/'test_predictions.csv');cfg=json.loads((OUT/'runs'/name/'complete.json').read_text());df=load();man=pd.read_csv(OUT/'runs'/name/'split_manifest.csv');tr=man.loc[man.role.eq('train'),'row_index'].to_numpy()
    cols=json.loads((REPO/'results/tables/conformal_run_config.json').read_text())['proteins']
    trainx=np.nanmean(df.iloc[tr][cols].to_numpy(float),axis=1);trainy=df.iloc[tr][cfg['selected_markers']].to_numpy(float).sum(axis=1)
    xgrid=np.linspace(*np.quantile(trainx,[.01,.99]),80);ygrid=np.linspace(*np.quantile(trainy,[.01,.99]),80);xx,yy=np.meshgrid(xgrid,ygrid)
    sx,sy=trainx.std(),trainy.std();grid=np.column_stack([xx.ravel(),yy.ravel()]);pcoords=pred[['mean_E_44','phenotype_burden']].to_numpy();scale=np.array([sx,sy])
    maps={};support={}
    for h in [.45,.75]:
        dist=((grid[:,None,:]/scale-pcoords[None,:,:]/scale)**2).sum(axis=2);w=np.exp(-.5*dist/h**2);sw=w.sum(axis=1);ne=sw**2/np.maximum((w*w).sum(axis=1),1e-300);wn=w/np.maximum(sw[:,None],1e-300)
        if h==.45:
            mean=wn@pred.y.to_numpy();var=wn@(pred.y.to_numpy()**2)-mean**2;maps['mean']=mean;maps['variance']=np.maximum(var,0);maps['width']=wn@pred.width_joint.to_numpy()
            for k in ['mean','variance','width']:support[k]=ne.copy()
        else:maps['coverage']=wn@pred.covered_joint.to_numpy();support['coverage']=ne.copy()
    fig,axs=plt.subplots(2,2,figsize=(11.5,8.5),layout='constrained');export={'mean_E':grid[:,0],'phenotype_burden':grid[:,1]}
    for ax,key,title,label in zip(axs.flat,['mean','variance','width','coverage'],['a  Local mean endpoint','b  Local endpoint variance','c  Mean local interval width','d  Observed local coverage'],['-log10 LD50 (mol/kg)','Endpoint variance','Full interval width','Coverage']):
        values=maps[key].copy();values[support[key]<120]=np.nan;export[key]=values;export[key+'_neff']=support[key]
        cmap=plt.get_cmap('viridis' if key!='coverage' else 'cividis').copy();cmap.set_bad('#ededed');z=values.reshape(xx.shape)
        im=ax.pcolormesh(xx,yy,z,shading='auto',cmap=cmap);fig.colorbar(im,ax=ax,label=label,shrink=.9)
        ax.axvline(cfg['quadrant_thresholds']['mean_E'],lw=.7,ls=':',color='black');ax.axhline(cfg['quadrant_thresholds']['phenotype'],lw=.7,ls=':',color='black')
        ax.set(title=title,xlabel='Mean docking score (kcal/mol)',ylabel='Phenotype burden')
    save(fig,'figure5_revised_seed42');pd.DataFrame(export).to_csv(TABLE/'figure5_grid.csv',index=False)
    q=quad[quad.run.str.startswith('butina04_')];order=['Pheno+ MIE+','Pheno+ MIE-','Pheno- MIE+','Pheno- MIE-']
    fig,axs=plt.subplots(1,2,figsize=(11.6,4.5),layout='constrained')
    for seed in range(42,47):
        t=q[q.run.eq(f'butina04_seed{seed}')].set_index('quadrant').loc[order]
        axs[0].plot(range(4),t.joint_width,marker='o',lw=1,label=str(seed),alpha=.8)
        axs[1].plot(range(4),t.coverage*100,marker='o',lw=1,alpha=.8)
    for ax in axs:ax.set_xticks(range(4),[s.replace(' ','\n') for s in order])
    axs[0].set(ylabel='Mean full interval width',title='a  Quadrant interval widths');axs[0].legend(title='Split seed',frameon=False,ncol=2,fontsize=8)
    axs[1].axhline(90,color='black',lw=1,ls='--');axs[1].set(ylabel='Observed coverage (%)',title='b  Local coverage is uneven')
    save(fig,'quadrant_robustness')
    dumpfile={'primary_seed':42,'primary_mean_rmse':float(primary.rmse.mean()),'primary_mean_global_coverage':float(primary.global_coverage.mean()),'primary_mean_joint_coverage':float(primary.joint_coverage.mean()),'primary_mean_global_width':float(primary.global_width.mean()),'primary_mean_joint_width':float(primary.joint_width.mean()),'primary_mean_global_score':float(primary.global_score.mean()),'primary_mean_joint_score':float(primary.joint_score.mean()),'random_mean_rmse':float(random.rmse.mean()),'mean_relative_width_reduction':float(((primary.global_width-primary.joint_width)/primary.global_width).mean()),'mean_relative_score_reduction':float(((primary.global_score-primary.joint_score)/primary.global_score).mean())}
    (OUT/'report_numbers.json').write_text(json.dumps(dumpfile,indent=2)+'\n')
    print(json.dumps(dumpfile,indent=2))
if __name__=='__main__':main()
