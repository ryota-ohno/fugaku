import os ##ok
import subprocess ##ok
from utils import Rod
import numpy as np ##ok
import pandas as pd ##ok

def get_monomer_xyzR(monomer_name,Ta,Tb,Tc,A2,A3,phi):  
    T_vec = np.array([Ta,Tb,Tc])
    df_mono=pd.read_csv(f'/vol0303/data/hp260444/Working/fugaku/amber/BXBX/monomer/{monomer_name}.csv')
    atoms_array_xyzR=df_mono[['atom','X','Y','Z']].values
    xyz_array = atoms_array_xyzR[:,1:];atom_array = atoms_array_xyzR[:,0].reshape((-1,1))

    ex = np.array([1.,0.,0.]); ez = np.array([0.,0.,1.])
    xyz_array = np.matmul(xyz_array,Rod(-ex,A2).T)#
    xyz_array = np.matmul(xyz_array,Rod(ez,A3).T)#
    xyz_array = xyz_array + T_vec
    
    C0_index = 19;C1_index = 29####
    C0=xyz_array[C0_index];C1=xyz_array[C1_index];n1=C1-C0;n1/=np.linalg.norm(n1)
    
    xyz_array[C1_index:] = np.matmul((xyz_array[C1_index:]-C0),Rod(n1,phi).T) + C0
    return np.concatenate([xyz_array,atom_array],axis=1)
        
def get_xyzR_lines(monomer_name,xyzr_array):
    df_mono=pd.read_csv(f'/vol0303/data/hp260444/Working/fugaku/amber/BXBX/monomer/{monomer_name}.csv')
    mol=len(df_mono)

    bond_lines=[]
    with open(f'/vol0303/data/hp260444/Working/fugaku/amber/BXBX/monomer/{monomer_name}_bond.txt')as f:
        for line in f:
            num, a1, a2, btype = line.strip().split()              # ['1','1','2','ar']
            bond_lines.append([int(num), int(a1), int(a2), btype])
    bond=len(bond_lines)

    line1=f'@<TRIPOS>MOLECULE\nBXBX_dimer \n   {int(2*mol)}  {int(2*bond)}       2     0     0\nSMALL\nrbcc\n\n\n@<TRIPOS>ATOM\n'
    line2='@<TRIPOS>BOND\n'
    line3=f'@<TRIPOS>SUBSTRUCTURE\n     1 RES1        1 GROUP             0 ****  ****    0  \n     2 RES2       {int(mol+1)} GROUP             0 ****  ****    0 \n\n'

    para_list=[]
    with open(f'/vol0303/data/hp260444/Working/fugaku/amber/BXBX/monomer/{monomer_name}.mol2')as f:
        for line in f:
            #print(line)
            s=line.split()
            if len(s)==9:
                para_list.append([s[5],float(s[8])])
            if (line.find('BOND')>-1):
                break

    lines=[]
    lines.append(line1)
    for i in range(mol):
        x,y,z,atom=xyzr_array[i]
        atom_type,charge=para_list[i]
        lines.append(f'  {i+1} {atom} {x} {y} {z} {atom_type} 1 RES1 {charge}\n')
    for i in range(mol):
        x,y,z,atom=xyzr_array[i+mol]
        atom_type,charge=para_list[i]
        lines.append(f'  {i+1+mol} {atom} {x} {y} {z} {atom_type} 2 RES2 {charge}\n')   
    lines.append(line2)
    for bond,atom1,atom2,type in bond_lines:
        line=f'{bond} {atom1} {atom2} {type}\n'
        lines.append(line)
    for bond,atom1,atom2,type in bond_lines:
        line=f'{bond+len(bond_lines)} {atom1+mol} {atom2+mol} {type}\n'
        lines.append(line)
    lines.append(line3)
    return lines

# 実行ファイル作成
def get_one_exe(auto_dir,file_name,monomer_name):
    file_basename = os.path.splitext(file_name)[0]
    lines_job=[
'#!/bin/bash\n','\n',
'module use /vol0004/apps/isv/Amber20/modulefiles \n',####
'module load Amber20 \n','\n',#####
#f'parmchk2 -i {file_basename}.mol2 -f mol2 -o {file_basename}.frcmod\n',
f'tleap -f {file_basename}_tleap.in\n',
f'sander -O -i FF_calc.in -o {file_basename}.out -p {file_basename}.prmtop -c {file_basename}.inpcrd -r min.rst -ref {file_basename}.inpcrd\n'
f'rm {file_basename}.inpcrd\n',
f'rm {file_basename}.prmtop\n',
]
    
    lines_tleap=['source /vol0004/apps/isv/Amber20/Amber2021_20251202/dat/leap/cmd/leaprc.gaff\n',
f'MOL = loadmol2 {file_basename}.mol2\n',
f'loadamberparams {monomer_name}.frcmod\n',
f'saveamberparm MOL {file_basename}.prmtop {file_basename}.inpcrd\n',
'quit\n']
    file_job = os.path.join(auto_dir,f'amber/job_{file_basename}.sh')
    file_tleap = os.path.join(auto_dir,f'amber/{file_basename}_tleap.in')
    
    with open(file_job,'w')as f:
        f.writelines(lines_job)
    with open(file_tleap,'w')as f:
        f.writelines(lines_tleap)

    return file_job,f'{file_basename}.out'

######################################## 特化関数 ########################################

##################gaussview##################
def make_xyzfile(monomer_name,params_dict,structure_type):
    a = float(params_dict.get('a',0.0));b = float(params_dict.get('b',0.0)); z = float(params_dict.get('z',0.0))
    A2 = float(params_dict.get('A2',0.0)); A3 = float(params_dict.get('theta',0.0))
    phi = float(params_dict.get('phi',0.0))

    monomer_array_i = get_monomer_xyzR(monomer_name,0,0,0,A2,A3,phi)
    
    monomer_array_p1 = get_monomer_xyzR(monomer_name,a,0,0,A2,A3,phi)##1,2がb方向
    monomer_array_p2 = get_monomer_xyzR(monomer_name,0,b,2*z,A2,A3,phi)##1,2がb方向
    monomer_array_t1 = get_monomer_xyzR(monomer_name,a/2,b/2,z,A2,-A3,-phi)##1,2がb方向
    monomer_array_t2 = get_monomer_xyzR(monomer_name,-a/2,-b/2,-z,A2,-A3,-phi)##1,2がb方向
    
    xyz_list=['400 \n','polyacene9 \n']##4分子のxyzファイルを作成
    
    if structure_type == 1:##隣接8分子について対称性より3分子でエネルギー計算
        monomers_array_4 = np. concatenate([monomer_array_i,monomer_array_p1])
    elif structure_type == 2:##隣接8分子について対称性より3分子でエネルギー計算
        monomers_array_4 = np.concatenate([monomer_array_i,monomer_array_p2])
    elif structure_type == 3 :##隣接8分子について対称性より3分子でエネルギー計算
        monomers_array_4 = np.concatenate([monomer_array_i,monomer_array_t1])
    elif structure_type == 4 :##隣接8分子について対称性より3分子でエネルギー計算
        monomers_array_4 = np.concatenate([monomer_array_i,monomer_array_t2])
    
    for x,y,z,atom in monomers_array_4:
        line = '{} {} {} {}\n'.format(atom,x,y,z)     
        xyz_list.append(line)
    
    return xyz_list

def make_xyz(monomer_name,params_dict,structure_type):
    xyzfile_name = ''
    xyzfile_name += monomer_name
    for key,val in params_dict.items():
        val=float(val)
        xyzfile_name += '_{}_{}'.format(key,val)
    return xyzfile_name + f'_{structure_type}.xyz'

def make_gjf_xyz(auto_dir,monomer_name,params_dict,structure_type):
    a = float(params_dict.get('a',0.0));b = float(params_dict.get('b',0.0)); z = float(params_dict.get('z',0.0))
    A2 = float(params_dict.get('A2',0.0)); A3 = float(params_dict.get('theta',0.0))
    phi = float(params_dict.get('phi',0.0))

    monomer_array_i = get_monomer_xyzR(monomer_name,0,0,0,A2,A3,phi)
    
    monomer_array_p1 = get_monomer_xyzR(monomer_name,a,0,0,A2,A3,phi)##1,2がb方向
    monomer_array_p2 = get_monomer_xyzR(monomer_name,0,b,2*z,A2,A3,phi)##1,2がb方向
    monomer_array_t1 = get_monomer_xyzR(monomer_name,a/2,b/2,z,A2,-A3,-phi)##1,2がb方向
    monomer_array_t2 = get_monomer_xyzR(monomer_name,-a/2,-b/2,-z,A2,-A3,-phi)##1,2がb方向
    
    dimer_array_p1 = np.concatenate([monomer_array_i,monomer_array_p1]);dimer_array_p2 = np.concatenate([monomer_array_i,monomer_array_p2])
    dimer_array_t1 = np.concatenate([monomer_array_i,monomer_array_t1]);dimer_array_t2 = np.concatenate([monomer_array_i,monomer_array_t2])
    
    line_list_dimer_p1 = get_xyzR_lines(monomer_name,dimer_array_p1);line_list_dimer_p2 = get_xyzR_lines(monomer_name,dimer_array_p2)
    line_list_dimer_t1 = get_xyzR_lines(monomer_name,dimer_array_t1);line_list_dimer_t2 = get_xyzR_lines(monomer_name,dimer_array_t2)
    
    if structure_type == 1:##隣接8分子について対称性より3分子でエネルギー計算
        gij_xyz_lines = line_list_dimer_p1 
    elif structure_type == 2:##隣接8分子について対称性より3分子でエネルギー計算
        gij_xyz_lines = line_list_dimer_p2 
    elif structure_type == 3:##隣接8分子について対称性より3分子でエネルギー計算
        gij_xyz_lines = line_list_dimer_t1 
    elif structure_type == 4:##隣接8分子について対称性より3分子でエネルギー計算
        gij_xyz_lines = line_list_dimer_t2 
    
    file_name = get_file_name_from_dict(monomer_name,params_dict,structure_type)
    os.makedirs(os.path.join(auto_dir,'amber'),exist_ok=True)
    gij_xyz_path = os.path.join(auto_dir,'amber',file_name)
    with open(gij_xyz_path,'w') as f:
        f.writelines(gij_xyz_lines)
    
    return file_name

def get_file_name_from_dict(monomer_name,params_dict,structure_type):
    file_name = ''
    file_name += monomer_name
    for key,val in params_dict.items():
        val=float(val)
        file_name += '_{}_{}'.format(key,val)
    return file_name + f'_{structure_type}.mol2'
    
def exec_gjf(auto_dir, monomer_name, params_dict,structure_type):
    xyz_dir = os.path.join(auto_dir,'gaussview')
    xyzfile_name = make_xyz(monomer_name, params_dict,structure_type)
    xyz_path = os.path.join(xyz_dir,xyzfile_name)
    xyz_list = make_xyzfile(monomer_name,params_dict,structure_type)
    with open(xyz_path,'w') as f:
        f.writelines(xyz_list)
    
    file_name = make_gjf_xyz(auto_dir, monomer_name, params_dict,structure_type)
    file_job,log_file_name = get_one_exe(auto_dir,file_name,monomer_name)
    subprocess.run(['chmod','+x',file_job])
    subprocess.run([file_job])
    return log_file_name
    
############################################################################################