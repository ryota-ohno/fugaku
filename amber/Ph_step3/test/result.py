##tetracene層内計算
import os
os.environ['HOME'] ='/data/group1/z40145w'
import pandas as pd
import argparse
import subprocess
import numpy as np

def result_process(args):
    subprocess.run(['rm','*.sh.*'])
    auto_dir = f'/vol0303/data/hp260444/Working/fugaku/amber/Ph_step3/{args.auto_dir}'
    df_tot=[]
    phi_list=[int(phi) for phi in np.linspace(0,170,18)]
    for phi in phi_list:
        dir_name = f'{phi}'
        path_dir=os.path.join(auto_dir,f'{dir_name}')
        df=pd.read_csv(os.path.join(path_dir,'step3.csv'))
        df_tot.append(df)
    df_=pd.concat(df_tot, ignore_index=True)
    df_.to_csv(os.path.join(auto_dir,'step3.csv'),index=False)
    
def update_value_in_df(df,index,key,value):
    df.loc[index,key]=value
    return df

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--isTest',action='store_true')
    parser.add_argument('--auto-dir',type=str,help='path to dir which includes gaussian, gaussview and csv')
    ##maxnum-machine2 がない
    args = parser.parse_args()

    print("----main process----")
    result_process(args)
    print("----finish process----")
    