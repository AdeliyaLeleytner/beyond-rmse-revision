"""Execute notebooks in fresh kernels and export standalone HTML previews."""
from pathlib import Path
import argparse,json,os,sys,time
import nbformat
from nbclient import NotebookClient
from nbconvert import HTMLExporter
ROOT=Path(__file__).resolve().parents[1]
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--pattern',default='*.ipynb');args=ap.parse_args()
    kernel=ROOT/'.runtime/kernels/revision';kernel.mkdir(parents=True,exist_ok=True)
    (kernel/'kernel.json').write_text(json.dumps({'argv':[sys.executable,'-m','ipykernel_launcher','-f','{connection_file}'],'display_name':'Revision Python','language':'python'}))
    os.environ['JUPYTER_PATH']=str(ROOT/'.runtime')
    html=ROOT/'.runtime/html';html.mkdir(parents=True,exist_ok=True)
    record=[]
    for path in sorted((ROOT/'notebooks').glob(args.pattern)):
        started=time.time();print('EXECUTING',path.name,flush=True)
        nb=nbformat.read(path,as_version=4)
        client=NotebookClient(nb,timeout=900,kernel_name='revision',resources={'metadata':{'path':str(ROOT)}})
        client.execute()
        for cell in nb.cells:
            assert all(o.output_type!='error' for o in cell.get('outputs',[]))
            # Keep portable, user-neutral execution metadata.
            cell.metadata.pop('execution',None)
        nb.metadata.kernelspec={'display_name':'Python 3 (revision)','language':'python','name':'python3'}
        nbformat.validate(nb);nbformat.write(nb,path)
        body,_=HTMLExporter(template_name='lab').from_notebook_node(nb)
        (html/(path.stem+'.html')).write_text(body)
        record.append({'notebook':path.name,'seconds':round(time.time()-started,2),'code_cells':sum(c.cell_type=='code' for c in nb.cells),'embedded_figures':sum('image/png' in o.get('data',{}) for c in nb.cells for o in c.get('outputs',[])),'errors':0})
        print('PASSED',record[-1],flush=True)
    report=ROOT/'provenance/notebook_execution.json'
    if report.exists() and args.pattern!='*.ipynb':
        old=json.loads(report.read_text());names={x['notebook'] for x in record};record=[x for x in old if x['notebook'] not in names]+record
    report.write_text(json.dumps(sorted(record,key=lambda x:x['notebook']),indent=2)+'\n')
if __name__=='__main__':main()
