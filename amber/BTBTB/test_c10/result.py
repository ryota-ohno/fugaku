##tetracene層内計算
import os
os.environ['HOME'] ='/vol0303/data/hp260444'
import pandas as pd
import argparse
import subprocess
import numpy as np

def main_process(args):
    auto_dir = f'/vol0303/data/hp260444/Working/fugaku/amber/BTBTB/{args.auto_dir}'
    df_init=pd.read_csv(os.path.join(auto_dir,'step1_init_params.csv'))
    while True:
        theta_list=[10.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0, 80.0]
        for theta in theta_list:
            dir_name = f'{theta}'
            path_dir=os.path.join(auto_dir,f'{dir_name}')
            df_init_=pd.read_csv(os.path.join(path_dir,'step1_init_params.csv'))
            if len(df_init_)==len(df_init_[df_init_['status']=='Done']):
                df_init.loc[df_init['theta'] == theta, 'status'] = 'Done'
                df_init.to_csv(os.path.join(auto_dir,'step1_init_params.csv'),index=False)
        df_init_done=df_init[df_init['status']=='Done']
        if len(df_init)==len(df_init_done):
            break
    print('Done')    

def result_process(args):
    subprocess.run(['rm','*.sh.*'])
    auto_dir = f'/vol0303/data/hp260444/Working/fugaku/amber/BTBTB/{args.auto_dir}'
    df_init=pd.read_csv(os.path.join(auto_dir,'step1_init_params.csv'))
    df_tot=[]
    theta_list=[10.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0, 80.0]
    for theta in theta_list:
        dir_name = f'{theta}'
        path_dir=os.path.join(auto_dir,f'{dir_name}')
        df=pd.read_csv(os.path.join(path_dir,'step1.csv'))
        df_tot.append(df)
    df_=pd.concat(df_tot, ignore_index=True)
    df_.to_csv(os.path.join(auto_dir,'step1.csv'),index=False)
    
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
    main_process(args)
    result_process(args)
    print("----finish process----")
    