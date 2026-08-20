##tetracene層内計算
import os ##ok
import numpy as np ##ok
import time ##ok
from make_8_xyz_new import exec_gjf##計算した点のxyzfileを出す
from utils import get_E
import argparse ##ok 
import shutil ##ok
import csv ##ok
import pandas as pd ##ok

def main_process(args):
    auto_dir = f'/vol0303/data/hp260444/Working/fugaku/amber/BTBTB/{args.auto_dir}'
    os.makedirs(auto_dir, exist_ok=True)
    os.makedirs(os.path.join(auto_dir,'amber'), exist_ok=True)
    os.makedirs(os.path.join(auto_dir,'gaussview'), exist_ok=True)
    amber_path=os.path.join(auto_dir,'amber')
    shutil.copy(f'/vol0303/data/hp260444/Working/fugaku/amber/BTBTB/src/FF_calc.in',amber_path)
    shutil.copy(f'/vol0303/data/hp260444/Working/fugaku/amber/BTBTB/monomer/{args.monomer_name}.frcmod',amber_path)
    
    auto_csv_path = os.path.join(auto_dir,'step1.csv')
    if not os.path.exists(auto_csv_path):
        df=pd.DataFrame(columns=['theta','A2','phi','a','b','z','E','E1','E2','E3','status'])
        df.to_csv(auto_csv_path,index=False)
        
    auto_csv_path1 = os.path.join(auto_dir,'step1_1.csv')
    if not os.path.exists(auto_csv_path1):
        df1=pd.DataFrame(columns=['theta','A2','phi','a','z','E1','status','file_name'])
        df1.to_csv(auto_csv_path1,index=False)
        
    auto_csv_path2 = os.path.join(auto_dir,'step1_2.csv')
    if not os.path.exists(auto_csv_path2):
        df2=pd.DataFrame(columns=['theta','A2','phi','b','z','E2','status','file_name'])
        df2.to_csv(auto_csv_path2,index=False)
        
    auto_csv_path3 = os.path.join(auto_dir,'step1_3.csv')
    if not os.path.exists(auto_csv_path3):
        df3=pd.DataFrame(columns=['theta','A2','phi','a','b','z','E3','status','file_name'])
        df3.to_csv(auto_csv_path3,index=False)                

    auto_csv_path4 = os.path.join(auto_dir,'step1_4.csv')
    if not os.path.exists(auto_csv_path4):
        df4=pd.DataFrame(columns=['theta','A2','phi','a','b','z','E3','status','file_name'])
        df4.to_csv(auto_csv_path3,index=False)                
            
    df_mono=pd.read_csv(f'/vol0303/data/hp260444/Working/fugaku/amber/BTBTB/monomer/{args.monomer_name}_mono.csv')
    os.chdir(os.path.join(auto_dir,'amber'))
    isOver = False
    while not(isOver):
        #check
        isOver = listen(auto_dir,args.monomer_name,df_mono,args.num_nodes)##argsの中身を取る
        time.sleep(5)

def listen(auto_dir,monomer_name,df_mono,num_nodes):##args自体を引数に取るか中身をばらして取るかの違い

    fixed_param_keys = ['theta','A2','phi','z'];opt_param_keys_1 = ['a'];opt_param_keys_2 = ['b']

    auto_csv_1 = os.path.join(auto_dir,'step1_1.csv');df_E_1 = pd.read_csv(auto_csv_1)
    df_prg_1 = df_E_1.loc[df_E_1['status']=='InProgress',fixed_param_keys+opt_param_keys_1+['file_name']]
    for idx, row in df_prg_1.iterrows():
        params_dict1_ = row[fixed_param_keys + opt_param_keys_1 + ['file_name']].to_dict()
        file_name1=params_dict1_['file_name']
        log_filepath1 = os.path.join(*[auto_dir,'amber',file_name1])
        phi1=params_dict1_['phi']
        E_mono1 = df_mono.loc[df_mono['phi'] == phi1, 'E'].iloc[0]
        if not(os.path.exists(log_filepath1)):
            continue
        E_list1=get_E(log_filepath1)
        if len(E_list1)!=1 :##get Eの長さは計算した分子の数
            continue
        else:
            E1 = round(float(E_list1[0]) - 2 * E_mono1, 4)##8分子に向けてep1,ep2作成　ep1:b ep2:a
            df_E_1.loc[idx, ['E1','status']] = [E1,'Done']
    df_E_1.to_csv(auto_csv_1,index=False)
    
    auto_csv_2 = os.path.join(auto_dir,'step1_2.csv');df_E_2 = pd.read_csv(auto_csv_2)
    df_prg_2 = df_E_2.loc[df_E_2['status']=='InProgress',fixed_param_keys+opt_param_keys_2+['file_name']]
    for idx, row in df_prg_2.iterrows():
        params_dict2_ = row[fixed_param_keys + opt_param_keys_2 + ['file_name']].to_dict()
        file_name2=params_dict2_['file_name']
        log_filepath2 = os.path.join(*[auto_dir,'amber',file_name2])
        phi2=params_dict2_['phi']
        E_mono2 = df_mono.loc[df_mono['phi'] == phi2, 'E'].iloc[0]
        if not(os.path.exists(log_filepath2)):
            continue
        E_list2=get_E(log_filepath2)
        if len(E_list2)!=1 :##get Eの長さは計算した分子の数
            continue
        else:
            E2 = round(float(E_list2[0]) - 2 * E_mono2, 4)##8分子に向けてep1,ep2作成　ep1:b ep2:a
            df_E_2.loc[idx, ['E2','status']] = [E2,'Done']
    df_E_2.to_csv(auto_csv_2,index=False)
    
    auto_csv_3 = os.path.join(auto_dir,'step1_3.csv');df_E_3 = pd.read_csv(auto_csv_3)
    df_prg_3 = df_E_3.loc[df_E_3['status']=='InProgress',fixed_param_keys+opt_param_keys_1+opt_param_keys_2+['file_name']]
    for idx, row in df_prg_3.iterrows():
        params_dict3_ = row[fixed_param_keys + opt_param_keys_1 + opt_param_keys_2 + ['file_name']].to_dict()
        file_name3=params_dict3_['file_name']
        log_filepath3 = os.path.join(*[auto_dir,'amber',file_name3])
        phi3=params_dict3_['phi']
        E_mono3 = df_mono.loc[df_mono['phi'] == phi3, 'E'].iloc[0]
        if not(os.path.exists(log_filepath3)):
            continue
        E_list3=get_E(log_filepath3)
        if len(E_list3)!=1 :##get Eの長さは計算した分子の数
            continue
        else:
            E3 = round(float(E_list3[0]) - 2 * E_mono3, 4)##8分子に向けてep1,ep2作成　ep1:b ep2:a
            df_E_3.loc[idx, ['E3','status']] = [E3,'Done']
    df_E_3.to_csv(auto_csv_3,index=False)

    auto_csv_4 = os.path.join(auto_dir,'step1_4.csv');df_E_4 = pd.read_csv(auto_csv_4)
    df_prg_4 = df_E_4.loc[df_E_4['status']=='InProgress',fixed_param_keys+opt_param_keys_1+opt_param_keys_2+['file_name']]
    for idx, row in df_prg_4.iterrows():
        params_dict4_ = row[fixed_param_keys + opt_param_keys_1 + opt_param_keys_2 + ['file_name']].to_dict()
        file_name4=params_dict4_['file_name']
        log_filepath4 = os.path.join(*[auto_dir,'amber',file_name4])
        phi4=params_dict4_['phi']
        E_mono4 = df_mono.loc[df_mono['phi'] == phi4, 'E'].iloc[0]
        if not(os.path.exists(log_filepath4)):
            continue
        E_list4=get_E(log_filepath4)
        if len(E_list4)!=1 :##get Eの長さは計算した分子の数
            continue
        else:
            E4 = round(float(E_list4[0]) - 2 * E_mono4, 4)##8分子に向けてep1,ep2作成　ep1:b ep2:a
            df_E_4.loc[idx, ['E4','status']] = [E4,'Done']
    df_E_4.to_csv(auto_csv_4,index=False)
        
    auto_csv = os.path.join(auto_dir,'step1.csv')
    df_E = pd.read_csv(auto_csv)
    df_prg = df_E.loc[df_E['status']=='InProgress',fixed_param_keys+opt_param_keys_1+opt_param_keys_2]
    
    for idx,row in df_prg.iterrows():
        params_dict1_ = {**row[fixed_param_keys + opt_param_keys_1].to_dict(), 'status': 'Done'}
        params_dict2_ = {**row[fixed_param_keys + opt_param_keys_2].to_dict(), 'status': 'Done'}
        params_dict3_ = {**row[fixed_param_keys + opt_param_keys_1 + opt_param_keys_2].to_dict(), 'status': 'Done'}
        params_dict4_ = {**row[fixed_param_keys + opt_param_keys_1 + opt_param_keys_2].to_dict(), 'status': 'Done'}
        s1=filter_df(df_E_1, params_dict1_)['E1'];s2=filter_df(df_E_2, params_dict2_)['E2'];s3=filter_df(df_E_3, params_dict3_)['E3'];s4=filter_df(df_E_4, params_dict4_)['E4']
        if (len(s1) == 0) or (len(s2) == 0) or (len(s3) == 0) or (len(s4) == 0):
            continue
        E1 = float(s1.values[0]);E2 = float(s2.values[0]);E3 = float(s3.values[0]);E4 = float(s4.values[0])
        E=(E1+E2+E3+E4)
        df_E.loc[idx, ['E','E1','E2','E3','E4','status']] = [E,E1,E2,E3,E4,'Done']
    df_E.to_csv(auto_csv,index=False)
#####実質的にはここで一回切るくらいのイメージ
    dict_matrix = get_params_dict(auto_dir,num_nodes)##更新分を流す a1/HOME/HASEGAWALABz2まで取得
    new=0
    if len(dict_matrix)!=0:#終わりがまだ見えないなら
        df_E= pd.read_csv(os.path.join(auto_dir,'step1.csv'))
        df_E_1 = pd.read_csv(auto_csv_1);df_E_2 = pd.read_csv(auto_csv_2);df_E_3 = pd.read_csv(auto_csv_3);df_E_4 = pd.read_csv(auto_csv_4)
        for i in range(len(dict_matrix)):
            params_dict=dict_matrix[i]
            params_dict1 = {k: v for k, v in params_dict.items() if (k in fixed_param_keys) or (k in opt_param_keys_1)}
            params_dict2 = {k: v for k, v in params_dict.items() if (k in fixed_param_keys) or (k in opt_param_keys_2)}
            params_dict3 = params_dict;params_dict4 = params_dict
            alreadyCalculated = check_calc_status(df_E,params_dict)
            if not(alreadyCalculated):
                new+=1
                df_E_filtered = filter_df(df_E, params_dict)
                if len(df_E_filtered) == 0:
                    df_newline = pd.Series({**params_dict,'E':0.,'E1':0.,'E2':0.,'E3':0.,'status':'InProgress'})
                    df_E=pd.concat([df_E,df_newline.to_frame().T],axis=0,ignore_index=True)
            alreadyCalculated1 = check_calc_status(df_E_1,params_dict1)
            if not(alreadyCalculated1):    
                file_name = exec_gjf(auto_dir, monomer_name, {**params_dict1}, structure_type=1)
                df_newline_1 = pd.Series({**params_dict1,'E1':0.,'status':'InProgress','file_name':file_name})
                df_E_1=pd.concat([df_E_1,df_newline_1.to_frame().T],axis=0,ignore_index=True)
            alreadyCalculated2 = check_calc_status(df_E_2,params_dict2)
            if not(alreadyCalculated2):    
                file_name = exec_gjf(auto_dir, monomer_name, {**params_dict2}, structure_type=2)
                df_newline_2 = pd.Series({**params_dict2,'E2':0.,'status':'InProgress','file_name':file_name})
                df_E_2=pd.concat([df_E_2,df_newline_2.to_frame().T],axis=0,ignore_index=True)
            alreadyCalculated3 = check_calc_status(df_E_3,params_dict3)
            if not(alreadyCalculated3):    
                file_name = exec_gjf(auto_dir, monomer_name, {**params_dict3}, structure_type=3)
                df_newline_3 = pd.Series({**params_dict3,'E3':0.,'status':'InProgress','file_name':file_name})
                df_E_3=pd.concat([df_E_3,df_newline_3.to_frame().T],axis=0,ignore_index=True)
            alreadyCalculated4 = check_calc_status(df_E_4,params_dict4)
            if not(alreadyCalculated4):    
                file_name = exec_gjf(auto_dir, monomer_name, {**params_dict4}, structure_type=4)
                df_newline_4 = pd.Series({**params_dict4,'E4':0.,'status':'InProgress','file_name':file_name})
                df_E_4=pd.concat([df_E_4,df_newline_4.to_frame().T],axis=0,ignore_index=True)
                        
    if new>0:
        df_E.to_csv(auto_csv,index=False);df_E_1.to_csv(auto_csv_1,index=False);df_E_2.to_csv(auto_csv_2,index=False);df_E_3.to_csv(auto_csv_3,index=False);df_E_4.to_csv(auto_csv_4,index=False)

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
    fixed_param_keys = ['theta','A2','phi','z'];opt_param_keys_1 = ['a'];opt_param_keys_2 = ['b']
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
                params_dict = df_init_params.loc[index,fixed_param_keys+opt_param_keys_1+opt_param_keys_2].to_dict()
                dict_matrix_init.append(params_dict)
                if len(df_init_params_inprogress) + len(dict_matrix_init) >= num_nodes:
                    return dict_matrix_init
            return dict_matrix_init
    
    dict_matrix=[]
    for index in df_init_params_inprogress.index:##こちら側はinit_params内のある業に関する探索が終わった際の新しい行での探索を開始するもの ###ここを改良すればよさそう
        df_init_params = pd.read_csv(init_params_csv)
        init_params_dict = df_init_params.loc[index,fixed_param_keys+opt_param_keys_1+opt_param_keys_2].to_dict()
        fixed_params_dict = df_init_params.loc[index,fixed_param_keys].to_dict()
        isDone, opt_params_matrix = get_opt_params_dict(df_cur, init_params_dict,fixed_params_dict)
        if isDone:
            opt_params_dict={'a':opt_params_matrix[0][0],'b':opt_params_matrix[0][1]}
            df_init_params = update_value_in_df(df_init_params,index,'status','Done')
            df_init_params.to_csv(init_params_csv,index=False)
        else:
            for i in range(len(opt_params_matrix)):
                opt_params_dict={'a':opt_params_matrix[i][0],'b':opt_params_matrix[i][1]}
                df_inprogress = filter_df(df_cur, {**fixed_params_dict,**opt_params_dict,'status':'InProgress'})
                if len(df_inprogress)>=1:
                    continue
                else:
                    d={**fixed_params_dict,**opt_params_dict}
                    dict_matrix.append(d)
    return dict_matrix
        
def get_opt_params_dict(df_cur, init_params_dict,fixed_params_dict):
    df_val = filter_df(df_cur, fixed_params_dict)
    a_init_prev = init_params_dict['a']; b_init_prev = init_params_dict['b']; theta = init_params_dict['theta']
    z = init_params_dict['z']; A2 = init_params_dict['A2']; phi = init_params_dict['phi']
    while True:
        E_list=[];heri_list=[]
        para_list=[]
        for a in [a_init_prev-0.1,a_init_prev,a_init_prev+0.1]:
            for b in [b_init_prev-0.1,b_init_prev,b_init_prev+0.1]:
                a = np.round(a,1);b = np.round(b,1)
                df_val_ab = df_val[
                    (df_val['a']==a)&(df_val['b']==b)&(df_val['theta']==theta)&
                    (df_val['z']==z)&(df_val['A2']==A2)&(df_val['phi']==phi)&
                    (df_val['status']=='Done')]
                if len(df_val_ab)==0:
                    para_list.append([a,b])
                    continue
                heri_list.append([a,b]);E_list.append(df_val_ab['E'].values[0])
        if len(para_list) != 0:
            return False,para_list
        a_init,b_init = heri_list[np.argmin(np.array(E_list))]
        if a_init==a_init_prev and b_init==b_init_prev:
            return True,[[a_init,b_init]]
        else:
            a_init_prev=a_init;b_init_prev=b_init

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
    