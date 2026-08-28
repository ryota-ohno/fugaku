##tetracene層内計算
import os ##ok
import numpy as np ##ok
import time ##ok
from make_step3 import exec_gjf##計算した点のxyzfileを出す
from utils import get_E
import argparse ##ok 
import shutil ##ok
import csv ##ok
import pandas as pd ##ok

def main_process(args):
    auto_dir = f'/vol0303/data/hp260444/Working/fugaku/amber/Ph_step3/{args.auto_dir}'
    os.makedirs(auto_dir, exist_ok=True)
    os.makedirs(os.path.join(auto_dir,'amber'), exist_ok=True)
    os.makedirs(os.path.join(auto_dir,'gaussview'), exist_ok=True)
    amber_path=os.path.join(auto_dir,'amber')
    shutil.copy(f'/vol0303/data/hp260444/Working/fugaku/amber/Ph_step3/src/FF_calc.in',amber_path)
    shutil.copy(f'/vol0303/data/hp260444/Working/fugaku/amber/Ph_step3/monomer/{args.monomer_name}.frcmod',amber_path)
    
    auto_csv_path = os.path.join(auto_dir,'step1.csv')
    if not os.path.exists(auto_csv_path):
        df=pd.DataFrame(columns=['cx','cy','cz','theta','a','b','z','A2','phi','E','E1','E2','E3','E4','E5','E6','E7','E8','E9','status'])
        df.to_csv(auto_csv_path,index=False)

    for i in range(1,8):
        path = os.path.join(auto_dir,f'step1_{i}.csv')
        if not os.path.exists(path):
            df_=pd.DataFrame(columns=['cx','cy','cz','theta','a','b','z','A2','phi',f'E{i}','status','file_name'])
            df_.to_csv(path,index=False)
            
    file_mono=f'/vol0303/data/hp260444/Working/fugaku/amber/Ph_step3/monomer/{args.monomer_name}.out'
    E_mono=get_E(file_mono)[0]
    os.chdir(os.path.join(auto_dir,'amber'))
    isOver = False
    while not(isOver):
        #check
        isOver = listen(auto_dir,args.monomer_name,E_mono,args.num_nodes)##argsの中身を取る
        time.sleep(1)

def listen(auto_dir,monomer_name,E_mono,num_nodes):##args自体を引数に取るか中身をばらして取るかの違い

    fixed_param_keys = ['theta','a','b','A2','z','phi'];opt_param_keys = ['cx','cy','cz']

    for i in range(1,8):
        auto_csv_=os.path.join(auto_dir,f'step1_{i}.csv');df_E_ = pd.read_csv(auto_csv_)
        df_prg_ = df_E_.loc[df_E_['status']=='InProgress',fixed_param_keys+opt_param_keys+['file_name']]
        for idx, row in df_prg_.iterrows():
            params_dict_ = row[fixed_param_keys + opt_param_keys + ['file_name']].to_dict()
            file_name_=params_dict_['file_name']
            log_filepath_ = os.path.join(*[auto_dir,'amber',file_name_])
            if not(os.path.exists(log_filepath_)):
                continue
            E_list_=get_E(log_filepath_)
            if len(E_list_)!=1 :##get Eの長さは計算した分子の数
                continue
            else:
                E = round(float(E_list_[0]) - 2 * E_mono, 4)##8分子に向けてep1,ep2作成　ep1:b ep2:a
                df_E_.loc[idx, [f'E{i}','status']] = [E,'Done']
        df_E_.to_csv(auto_csv_,index=False)
        
    auto_csv = os.path.join(auto_dir,'step1.csv')
    df_E = pd.read_csv(auto_csv)
    df_prg = df_E.loc[df_E['status']=='InProgress',fixed_param_keys+opt_param_keys]
    for idx,row in df_prg.iterrows():
        params_dict = {**row[fixed_param_keys + opt_param_keys].to_dict(), 'status': 'Done'}
        E_list=[]
        for i in range(1,8):
            auto_csv_=os.path.join(auto_dir,f'step1_{i}.csv');df_E_ = pd.read_csv(auto_csv_)
            s_= filter_df(df_E_, params_dict)[f'E{i}']
            if len(s_)>0:
                E_list.append(s_.values[0])
        if len(E_list)<7:
            continue
        E1,E2,E3,E4,E5,E6,E7=E_list
        E=E1+E2+E3+E4+E5+E6*2+E7*2
        df_E.loc[idx, ['E','E1','E2','E3','E4','E5','E6','E7','E8','E9','status']] = [E,E1,E2,E3,E4,E5,E6,E7,E6,E7,'Done']
    df_E.to_csv(auto_csv,index=False)
#####実質的にはここで一回切るくらいのイメージ
    dict_matrix = get_params_dict(auto_dir,num_nodes)##更新分を流す a1/HOME/HASEGAWALABz2まで取得
    new=0
    if len(dict_matrix)!=0:#終わりがまだ見えないなら
        df_E= pd.read_csv(os.path.join(auto_dir,'step1.csv'))
        for params_dict in dict_matrix:
            alreadyCalculated = check_calc_status(df_E,params_dict)
            if not(alreadyCalculated):
                new+=1
                df_E_filtered = filter_df(df_E, params_dict)
                if len(df_E_filtered) == 0:
                    df_newline = pd.Series({**params_dict,'E':0.,'E1':0.,'E2':0.,'E3':0.,'E4':0.,'E5':0.,'E6':0.,'E7':0.,'E8':0.,'E9':0.,'status':'InProgress'})
                    df_E=pd.concat([df_E,df_newline.to_frame().T],axis=0,ignore_index=True)        

        df_E_list=[];file_path_list=[]
        for i in range(1,8):
            auto_csv_=os.path.join(auto_dir,f'step1_{i}.csv');df_E_ = pd.read_csv(auto_csv_)
            for params_dict in dict_matrix:
                alreadyCalculated_ = check_calc_status(df_E_,params_dict)
                if not(alreadyCalculated_):    
                    file_name = exec_gjf(auto_dir, monomer_name, {**params_dict}, structure_type=i)
                    df_newline_ = pd.Series({**params_dict,f'E{i}':0.,'status':'InProgress','file_name':file_name})
                    df_E_=pd.concat([df_E_,df_newline_.to_frame().T],axis=0,ignore_index=True)
            df_E_list.append(df_E_);file_path_list.append(auto_csv_)
                        
    if new>0:
        df_E.to_csv(auto_csv,index=False)
        for i in range(7):
            df_= df_E_list[i];file_path=file_path_list[i]
            df_.to_csv(file_path,index=False)
        
    init_params_csv=os.path.join(auto_dir, 'step1_init_params.csv')
    df_init_params = pd.read_csv(init_params_csv)
    df_init_params_done = filter_df(df_init_params,{'status':'Done'})
    isOver = True if len(df_init_params_done)==len(df_init_params) else False
    return isOver

def check_calc_status(df_E,params_dict):
    if len(df_E)==0:
        return False
    df_E_filtered = filter_df(df_E, params_dict)
    df_E_filtered = df_E_filtered.reset_index(drop=True)
    try:
        status = get_values_from_df(df_E_filtered,0,'status')
        return status=='Done'
    except KeyError:
        return False

def get_params_dict(auto_dir, num_nodes):
    fixed_param_keys = ['theta','a','b','A2','z','phi'];opt_param_keys = ['cx','cy','cz']
    init_params_csv=os.path.join(auto_dir, 'step1_init_params.csv')
    df_init_params = pd.read_csv(init_params_csv)
    df_cur = pd.read_csv(os.path.join(auto_dir, 'step1.csv'))
    df_init_params_inprogress = df_init_params[df_init_params['status']=='InProgress']
    df_init_params_notyet = df_init_params[df_init_params['status']=='NotYet']
    
    if len(df_init_params_notyet)>0.1:#最初の立ち上がり時
        dict_matrix_init=[]
        if len(df_init_params_inprogress) < num_nodes:
            df_init_params_notyet = df_init_params[df_init_params['status']=='NotYet']
            for index in df_init_params_notyet.index:
                df_init_params = update_value_in_df(df_init_params,index,'status','InProgress')
                df_init_params.to_csv(init_params_csv,index=False)
                params_dict = df_init_params.loc[index,fixed_param_keys+opt_param_keys].to_dict()
                dict_matrix_init.append(params_dict)
                if len(df_init_params_inprogress) + len(dict_matrix_init) >= num_nodes:
                    return dict_matrix_init
            return dict_matrix_init
    
    dict_matrix=[]
    for index in df_init_params_inprogress.index:##こちら側はinit_params内のある業に関する探索が終わった際の新しい行での探索を開始するもの ###ここを改良すればよさそう
        df_init_params = pd.read_csv(init_params_csv)
        init_params_dict = df_init_params.loc[index,fixed_param_keys+opt_param_keys].to_dict()
        fixed_params_dict = df_init_params.loc[index,fixed_param_keys].to_dict()
        isDone, opt_params_matrix = get_opt_params_dict(df_cur, init_params_dict,fixed_params_dict)
        if isDone:
            opt_params_dict={'cx':opt_params_matrix[0][0],'cy':opt_params_matrix[0][1],'cz':opt_params_matrix[0][2]}
            df_init_params = update_value_in_df(df_init_params,index,'status','Done')
            df_init_params.to_csv(init_params_csv,index=False)
        else:
            for j in range(len(opt_params_matrix)):
                opt_params_dict={'cx':opt_params_matrix[j][0],'cy':opt_params_matrix[j][1],'cz':opt_params_matrix[j][2]}
                df_inprogress = filter_df(df_cur, {**fixed_params_dict,**opt_params_dict,'status':'InProgress'})
                if len(df_inprogress)>=1:
                    continue
                else:
                    d={**fixed_params_dict,**opt_params_dict}
                    dict_matrix.append(d)
    return dict_matrix
        
def get_opt_params_dict(df_cur, init_params_dict,fixed_params_dict):
    df_val = filter_df(df_cur, fixed_params_dict)
    cx_init_prev = init_params_dict['cx']; cy_init_prev = init_params_dict['cy']; cz_init_prev = init_params_dict['cz']
    a = init_params_dict['a']; b = init_params_dict['b']; theta = init_params_dict['theta']
    z = init_params_dict['z']; A2 = init_params_dict['A2']; phi = init_params_dict['phi']
    while True:
        E_list=[];heri_list=[]
        para_list=[]
        for cx in [cx_init_prev]:
            for cy in [cy_init_prev-0.1,cy_init_prev,cy_init_prev+0.1]:
                for cz in [cz_init_prev-0.1,cz_init_prev,cz_init_prev+0.1]:
                    cx = np.round(cx,1);cy = np.round(cz,1);cz = np.round(cz,1)
                    df_val_ab = df_val[
                        (df_val['cz']==cz)&(df_val['cy']==cy)&(df_val['cz']==cz)&
                        (df_val['a']==a)&(df_val['b']==b)&(df_val['theta']==theta)&
                        (df_val['z']==z)&(df_val['A2']==A2)&(df_val['phi']==phi)&
                        (df_val['status']=='Done')]
                    if len(df_val_ab)==0:
                        para_list.append([cx,cy,cz])
                        continue
                    heri_list.append([cx,cy,cz]);E_list.append(df_val_ab['E'].values[0])
        if len(para_list) != 0:
            return False,para_list
        cx_init,cy_init,cz_init = heri_list[np.argmin(np.array(E_list))]
        if cx_init==cx_init_prev and cy_init==cy_init_prev and cz_init==cz_init_prev:
            return True,[[cx_init,cy_init,cz_init]]
        else:
            cx_init_prev=cx_init;cy_init_prev=cy_init;cz_init_prev=cz_init

def get_values_from_df(df,index,key):
    return df.loc[index,key]

def update_value_in_df(df,index,key,value):
    df.loc[index,key]=value
    return df

def filter_df(df, dict_filter):
    for k, v in dict_filter.items():
        if type(v)==str:
            df=df[df[k]==v]
        else:
            df=df[df[k]==v]
    df_filtered=df
    return df_filtered

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    
    parser.add_argument('--isTest',action='store_true')
    parser.add_argument('--auto-dir',type=str,help='path to dir which includes amber, gaussview and csv')
    parser.add_argument('--monomer-name',type=str,help='monomer name')
    parser.add_argument('--num-nodes',type=int,help='num nodes')
    args = parser.parse_args()

    print("----main process----")
    main_process(args)
    print("----finish process----")
    