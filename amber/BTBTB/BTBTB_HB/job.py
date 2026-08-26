##tetracene層内計算
import os
os.environ['HOME'] ='/vol0303/data/hp260444'
import pandas as pd
import argparse
import subprocess
import numpy as np

def init_process(args):
    auto_dir = f'/vol0303/data/hp260444/Working/fugaku/amber/BTBTB/{args.auto_dir}'
    monomer_name=args.monomer_name
    df_init=pd.read_csv(os.path.join(auto_dir,'step1_init_params.csv'))
    z_list=[np.round(z,1) for z in np.linspace(-2.0,2.0,21)]
    for z in z_list:
        dir_name = f'{z}'
        os.makedirs(os.path.join(auto_dir,f'{dir_name}'), exist_ok=True)
        df_init_=df_init[df_init['z'] == z]
        df_init_.to_csv(os.path.join(auto_dir,f'{dir_name}/step1_init_params.csv'),index=False)
        os.chdir(os.path.join(auto_dir,f'{dir_name}'))
        job_lines=[
        '#!/bin/bash \n',
        '#PJM -L "rscgrp=small"\n',
        '#PJM -L "node=1"\n',
        '#PJM -L "elapse=2:20:00"\n',
        '#PJM -L "freq=2200,eco_state=2"\n',
        '#PJM -g hp260444\n',
        '#PJM -x PJM_LLIO_GFSCACHE=/vol0004:/vol0003\n',
        '#PJM --llio localtmp-size=20Gi\n',
        '#PJM -S\n',
        '#PJM "--norestart"\n',
        '\n',
        'source /vol0303/data/hp260444/venv/bin/activate\n',
        '\n',
        f'python /vol0303/data/hp260444/Working/fugaku/amber/BTBTB/src/step1_8_xyz.py --auto-dir {args.auto_dir}/{dir_name} --monomer-name {monomer_name} --num-nodes 5\n',
        '\n',
        '#sleep 12 \n'
            ]
        with open(os.path.join(auto_dir,f'{dir_name}/job.sh'),'w')as f:
            f.writelines(job_lines)
        subprocess.run(['pjsub',os.path.join(auto_dir,f'{dir_name}/job.sh')])

def update_value_in_df(df,index,key,value):
    df.loc[index,key]=value
    return df

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--isTest',action='store_true')
    parser.add_argument('--auto-dir',type=str,help='path to dir which includes gaussian, gaussview and csv')
    parser.add_argument('--monomer-name',type=str,help='name of monomer to be calculated')
    ##maxnum-machine2 がない
    args = parser.parse_args()

    print("----main process----")
    init_process(args)
    print("----finish process----")    