import os.path

from Scripts.Utils import *
from Scripts.Classes import *
from Scripts.Global_Value import *
import Scripts.Global_Value
from Scripts.Prof import Compute_B_Factor
from Scripts.Run_DisEMBL import Run_DisEMBL,Generate_Res_DisEMBL
# from Caps import Compute_Co_Evo
from Scripts.Run_Sift import Run_Sift
from Scripts.Pymol import *
from Scripts.Run_Ring3 import Run_Ring,Judge_Bond_of_Ring,Devide_Res_of_Ring_by_Layers
from Scripts.MSA import find_pssm_score
import Scripts.AAindex
from Scripts.AAindex import Get_Mutation_Index_List_from_Matrix,Get_Mutation_Index_List_from_Index
import multiprocessing
import signal
from Scripts.Log import Log

def Signal_Handler(sig, frame):
    print("Received signal to terminate.")
    Log("Received signal to terminate.")
    os.kill(os.getpid(), signal.SIGTERM)


data_list=[]
def Feature_Extraction(table_path, table_name, features_obj_list:list, process_num:int):
    Log('Reading task table for Features Extraction')
    print('Reading task table for Features Extraction')
    with open(table_path+table_name,'r') as table:
        lines=table.readlines()
        if lines[0]!='id,wt_aa_short,mut_aa_short,loc,t_loc,wt_pdb_name,wt_pdb_path,mut_pdb_name,mut_pdb_path,wt_fasta_path,mut_fasta_path,wt_pssm_path,mut_pssm_path,wt_psi_blast_path,mut_psi_blast_path,wt_blastp_path,mut_blastp_path,pH,temperature,ddg,is_beta\n':
            error_obj.Something_Wrong(Feature_Extraction.__name__)
            exit(1)
        for line in lines[1:]:
            if line!='' and line!='\n':
                data_list.append(line.replace('\n',''))

    Log('Checking if all data have different WT and MUT AA')
    print('Checking if all data have different WT and MUT AA')

    for i in range(len(data_list)):
        item_list = str(data_list[i]).split(',')
        ID = item_list[0]
        postfix=ID.split('_')[-2]
        if postfix[0]==postfix[-1]:
            Log(f'{ID} has same WT and MUT AA, has been removed')
            print(f'{ID} has same WT and MUT AA, has been removed')
            data_list[i]='remove'

    for i in range(len(data_list) - 1, -1, -1):
        if data_list[i]=='remove':
            data_list.pop(i)

    temp_list=[]
    for i in range(len(data_list)):
        item_list=str(data_list[i]).split(',')
        ID=item_list[0]
        if ID in temp_list:
            Log(f'{ID} is repeated, has been removed')
            print(f'{ID} is repeated, has been removed')
            data_list[i]='remove'
        else:
            temp_list.append(ID)

    for i in range(len(data_list) - 1, -1, -1):
        if data_list[i]=='remove':
            data_list.pop(i)

    if len(data_list)!=len(temp_list):
        error_obj.Something_Wrong(Feature_Extraction.__name__)
        exit(1)


    Log('Aligning PDB with Pymol')
    print('Aligning PDB with Pymol')
    for data in data_list:
        item_list=str(data).split(',')
        if len(item_list)!=21:
            error_obj.Something_Wrong(Feature_Extraction.__name__)
            exit(1)
        wt_pdb_path=item_list[6]
        mut_pdb_path=item_list[8]
        Pymol_Clean_Align_PDB_Pair(wt_pdb_path, mut_pdb_path, wt_pdb_path, mut_pdb_path)

        Is_Beta=item_list[20]
        if Is_Beta=='1':
            Change_TER(wt_pdb_path)
            Change_TER(mut_pdb_path)

    Log('Begin with multiple process')
    print('Begin with multiple process')

    task_count=0
    pool = multiprocessing.Pool(process_num)
    process_res_list = []
    signal.signal(signal.SIGINT, Signal_Handler)
    for data in data_list:
        item_list=str(data).split(',')
        if len(item_list)!=21:
            error_obj.Something_Wrong(Feature_Extraction.__name__)
            exit(1)
        task_count+=1
        obj=Feature_Object()

        # if not Detail_Extraction(obj,item_list):
        #     error_obj.Something_Wrong(Feature_Extraction.__name__)
        #     exit(1)
        #features_obj_list.append(obj)

        arg=(obj,item_list,task_count)
        res=pool.apply_async(Detail_Extraction,arg)
        process_res_list.append(res)

    pool.close()
    pool.join()

    Log('End with multiple process and check multiple process results')
    print('End with multiple process and check multiple process results')

    if len(process_res_list)!=task_count:
        Log('Multiple process results wrong, program will abort')
        print('Multiple process results wrong, program will abort')
        exit(1)

    res_count=0
    failed_count=0
    for process_res in process_res_list:
        item_list = str(data_list[res_count]).split(',')
        ID = item_list[0]
        res_count+=1

        if not process_res.successful():
            Log(f'task {res_count}, ID {ID} has false return, maybe has exception in there')
            print(f'task {res_count}, ID {ID} has false return, maybe has exception in there')
            failed_count += 1
        else:
            if process_res.get() is False:
                Log(f'task {res_count}, ID {ID} has failed and been filtered')
                print(f'task {res_count}, ID {ID} has failed and been filtered')
                failed_count += 1
            else:
                obj = process_res.get()
                if isinstance(obj, Feature_Object):
                    features_obj_list.append(obj)
                else:
                    Log(f'task {res_count}, ID {ID} has unknown return and been filtered')
                    print(f'task {res_count}, ID {ID} has unknown return and been filtered')
                    failed_count += 1

    Log(f'Has failed {failed_count} data')
    print(f'Has failed {failed_count} data')









def Detail_Extraction(obj:Feature_Object,basic_list:list,task_count:int):
    Log(f'Processing task {task_count}')
    print(f'Processing task {task_count}')
    try:
        Log(f'Task {task_count}, ID (waiting to get): Features Extraction 1: Extracting task data')
        print(f'Task {task_count}, ID (waiting to get): Features Extraction 1: Extracting task data')
        obj.ID = basic_list[0]
        obj.WT_Amino_Acid_short = basic_list[1]
        obj.MUT_Amino_Acid_short = basic_list[2]
        obj.Loc_of_Mutation = int(basic_list[3])
        obj.True_Loc_of_Mutation = int(basic_list[4])
        obj.WT_Structure.PDB_Name = basic_list[5]
        obj.WT_Structure.PDB_path = basic_list[6]
        obj.MUT_Structure.PDB_Name = basic_list[7]
        obj.MUT_Structure.PDB_path = basic_list[8]
        obj.WT_Sequence_path = basic_list[9]
        obj.MUT_Sequence_path = basic_list[10]
        obj.WT_PSSM_Path = basic_list[11]
        obj.MUT_PSSM_Path = basic_list[12]
        obj.WT_PSI_BLAST_Path = basic_list[13]
        obj.MUT_PSI_BLAST_Path = basic_list[14]
        obj.WT_BLASTP_Path = basic_list[15]
        obj.MUT_BLASTP_Path = basic_list[16]
        obj.pH = float(basic_list[17])
        obj.Temperature = float(basic_list[18])
        obj.Experimental_DDG = float(basic_list[19])
        obj.Is_Beta=basic_list[20]
        if obj.Experimental_DDG>0.5:
            obj.Experimental_DDG_Classification=1
        elif obj.Experimental_DDG<-0.5:
            obj.Experimental_DDG_Classification=-1
        else:
            obj.Experimental_DDG_Classification=0
    except:
        error_obj.Something_Wrong(Detail_Extraction.__name__)
        return False


    #clean pdb
    # print('Features Extraction 2: Aligning PDB with Pymol')
    # res_list=Pymol_Clean_Align_PDB_Pair(obj.WT_Structure.PDB_path,obj.MUT_Structure.PDB_path,obj.WT_Structure.PDB_path,obj.MUT_Structure.PDB_path)
    # obj.RMSD_WT_MUT=res_list[0]


    #####
    # WT_Amino_Acid
    Log(f'Task {task_count}, ID {obj.ID}: Features Extraction 2: Generating basic info')
    print(f'Task {task_count}, ID {obj.ID}: Features Extraction 2: Generating basic info')
    if not Get_Reasearched_Amino_Acid(obj.WT_Amino_Acid, obj.WT_Structure.PDB_Name, obj.WT_Structure.PDB_path, obj.True_Loc_of_Mutation, obj.WT_Amino_Acid_short):
        error_obj.Something_Wrong(Detail_Extraction.__name__)
        return False

    # MUT_Amino_Acid
    if not Get_Reasearched_Amino_Acid(obj.MUT_Amino_Acid, obj.MUT_Structure.PDB_Name, obj.MUT_Structure.PDB_path, obj.True_Loc_of_Mutation, obj.MUT_Amino_Acid_short):
        error_obj.Something_Wrong(Detail_Extraction.__name__)
        return False

    # WT_Amino_Acid_List
    if not Get_All_Amino_Acid(obj.WT_Amino_Acid_List, obj.WT_Structure.PDB_Name, obj.WT_Structure.PDB_path):
        error_obj.Something_Wrong(Detail_Extraction.__name__)
        return False
    obj.WT_Amino_Acid_List_Layer1=Get_Surrounding_AA(obj.WT_Amino_Acid,obj.WT_Amino_Acid_List,obj.Cutoff1)
    obj.WT_Amino_Acid_List_Layer2 = Get_Surrounding_AA(obj.WT_Amino_Acid, obj.WT_Amino_Acid_List, obj.Cutoff2)
    obj.WT_Amino_Acid_List_Layer3 = Get_Surrounding_AA(obj.WT_Amino_Acid, obj.WT_Amino_Acid_List, obj.Cutoff3)

    # MUT_Amino_Acid_List
    if not Get_All_Amino_Acid(obj.MUT_Amino_Acid_List, obj.MUT_Structure.PDB_Name, obj.MUT_Structure.PDB_path):
        error_obj.Something_Wrong(Detail_Extraction.__name__)
        return False
    obj.MUT_Amino_Acid_List_Layer1=Get_Surrounding_AA(obj.MUT_Amino_Acid,obj.MUT_Amino_Acid_List,obj.Cutoff1)
    obj.MUT_Amino_Acid_List_Layer2 = Get_Surrounding_AA(obj.MUT_Amino_Acid, obj.MUT_Amino_Acid_List, obj.Cutoff2)
    obj.MUT_Amino_Acid_List_Layer3 = Get_Surrounding_AA(obj.MUT_Amino_Acid, obj.MUT_Amino_Acid_List, obj.Cutoff3)

    # WT_Seq
    if not Read_Seq_from_AA_List(obj.WT_Seq,obj.WT_Amino_Acid_List):
        error_obj.Something_Wrong(Detail_Extraction.__name__)
        return False

    # MUT_Seq
    if not Read_Seq_from_AA_List(obj.MUT_Seq,obj.MUT_Amino_Acid_List):
        error_obj.Something_Wrong(Detail_Extraction.__name__)
        return False

    diff_count = 0
    for key in obj.WT_Seq.keys():
        for i in range(len(obj.WT_Seq[key])):
            if obj.WT_Seq[key][i] != obj.MUT_Seq[key][i]:
                diff_count += 1

    if diff_count != 1:
        error_obj.Something_Wrong(__name__, 'diff_count')
        return False

    res=Fetch_Chain_ID_from_Seq(obj.True_Loc_of_Mutation,obj.WT_Seq,obj.WT_Amino_Acid_short)
    if res is False:
        error_obj.Something_Wrong(Detail_Extraction.__name__)
        return False
    else:
        obj.Chain_ID_of_Mut=res

    ######



    #Ring_Bond_List, Num_HBOND_Ring, Num_SSBOND_Ring, Num_IONIC_Ring, Num_VDW_Ring, Num_PICATION_Ring, Num_PIPISTACK_Ring, Num_IAC_Ring,
    Log(f'Task {task_count}, ID {obj.ID}: Features Extraction 3: Running Ring3')
    print(f'Task {task_count}, ID {obj.ID}: Features Extraction 3: Running Ring3')
    if os.path.isfile(Ring_Path+'ring'):
        res_dict=Run_Ring(obj.WT_Structure.PDB_path,Ring_Path,obj.WT_Ring_Bond_List,TMP_Path,f'ring3_res_{obj.ID}_WT')
        if res_dict is False:
            error_obj.Something_Wrong(Detail_Extraction.__name__)
        else:
            obj.WT_Num_HBOND_Ring=res_dict['HBOND']
            obj.WT_Num_SSBOND_Ring=res_dict['SSBOND']
            obj.WT_Num_IONIC_Ring=res_dict['IONIC']
            obj.WT_Num_VDW_Ring=res_dict['VDW']
            obj.WT_Num_PICATION_Ring=res_dict['PICATION']
            obj.WT_Num_PIPISTACK_Ring=res_dict['PIPISTACK']


            res_dict=Devide_Res_of_Ring_by_Layers(obj.WT_Ring_Bond_List,obj.WT_Amino_Acid_List_Layer1)
            obj.WT_Num_HBOND_Ring_Layer1 = res_dict['HBOND']
            obj.WT_Num_SSBOND_Ring_Layer1 = res_dict['SSBOND']
            obj.WT_Num_IONIC_Ring_Layer1 = res_dict['IONIC']
            obj.WT_Num_VDW_Ring_Layer1 = res_dict['VDW']
            obj.WT_Num_PICATION_Ring_Layer1 = res_dict['PICATION']
            obj.WT_Num_PIPISTACK_Ring_Layer1 = res_dict['PIPISTACK']

            res_dict=Devide_Res_of_Ring_by_Layers(obj.WT_Ring_Bond_List,obj.WT_Amino_Acid_List_Layer2)
            obj.WT_Num_HBOND_Ring_Layer2 = res_dict['HBOND']
            obj.WT_Num_SSBOND_Ring_Layer2 = res_dict['SSBOND']
            obj.WT_Num_IONIC_Ring_Layer2 = res_dict['IONIC']
            obj.WT_Num_VDW_Ring_Layer2 = res_dict['VDW']
            obj.WT_Num_PICATION_Ring_Layer2 = res_dict['PICATION']
            obj.WT_Num_PIPISTACK_Ring_Layer2 = res_dict['PIPISTACK']

            res_dict=Devide_Res_of_Ring_by_Layers(obj.WT_Ring_Bond_List,obj.WT_Amino_Acid_List_Layer3)
            obj.WT_Num_HBOND_Ring_Layer3 = res_dict['HBOND']
            obj.WT_Num_SSBOND_Ring_Layer3 = res_dict['SSBOND']
            obj.WT_Num_IONIC_Ring_Layer3 = res_dict['IONIC']
            obj.WT_Num_VDW_Ring_Layer3 = res_dict['VDW']
            obj.WT_Num_PICATION_Ring_Layer3 = res_dict['PICATION']
            obj.WT_Num_PIPISTACK_Ring_Layer3 = res_dict['PIPISTACK']

            res_dict = Run_Ring(obj.MUT_Structure.PDB_path, Ring_Path, obj.MUT_Ring_Bond_List,TMP_Path,f'ring3_res_{obj.ID}_MUT')
            if res_dict is False:
                error_obj.Something_Wrong(Detail_Extraction.__name__)
            else:
                obj.MUT_Num_HBOND_Ring = res_dict['HBOND']
                obj.MUT_Num_SSBOND_Ring = res_dict['SSBOND']
                obj.MUT_Num_IONIC_Ring = res_dict['IONIC']
                obj.MUT_Num_VDW_Ring = res_dict['VDW']
                obj.MUT_Num_PICATION_Ring = res_dict['PICATION']
                obj.MUT_Num_PIPISTACK_Ring = res_dict['PIPISTACK']

                res_dict = Devide_Res_of_Ring_by_Layers(obj.MUT_Ring_Bond_List, obj.MUT_Amino_Acid_List_Layer1)
                obj.MUT_Num_HBOND_Ring_Layer1 = res_dict['HBOND']
                obj.MUT_Num_SSBOND_Ring_Layer1 = res_dict['SSBOND']
                obj.MUT_Num_IONIC_Ring_Layer1 = res_dict['IONIC']
                obj.MUT_Num_VDW_Ring_Layer1 = res_dict['VDW']
                obj.MUT_Num_PICATION_Ring_Layer1 = res_dict['PICATION']
                obj.MUT_Num_PIPISTACK_Ring_Layer1 = res_dict['PIPISTACK']

                res_dict = Devide_Res_of_Ring_by_Layers(obj.MUT_Ring_Bond_List, obj.MUT_Amino_Acid_List_Layer2)
                obj.MUT_Num_HBOND_Ring_Layer2 = res_dict['HBOND']
                obj.MUT_Num_SSBOND_Ring_Layer2 = res_dict['SSBOND']
                obj.MUT_Num_IONIC_Ring_Layer2 = res_dict['IONIC']
                obj.MUT_Num_VDW_Ring_Layer2 = res_dict['VDW']
                obj.MUT_Num_PICATION_Ring_Layer2 = res_dict['PICATION']
                obj.MUT_Num_PIPISTACK_Ring_Layer2 = res_dict['PIPISTACK']

                res_dict = Devide_Res_of_Ring_by_Layers(obj.MUT_Ring_Bond_List, obj.MUT_Amino_Acid_List_Layer3)
                obj.MUT_Num_HBOND_Ring_Layer3 = res_dict['HBOND']
                obj.MUT_Num_SSBOND_Ring_Layer3 = res_dict['SSBOND']
                obj.MUT_Num_IONIC_Ring_Layer3 = res_dict['IONIC']
                obj.MUT_Num_VDW_Ring_Layer3 = res_dict['VDW']
                obj.MUT_Num_PICATION_Ring_Layer3 = res_dict['PICATION']
                obj.MUT_Num_PIPISTACK_Ring_Layer3 = res_dict['PIPISTACK']



                obj.Diff_Num_HBOND_Ring = obj.MUT_Num_HBOND_Ring - obj.WT_Num_HBOND_Ring
                obj.Diff_Num_SSBOND_Ring = obj.MUT_Num_SSBOND_Ring - obj.WT_Num_SSBOND_Ring
                obj.Diff_Num_IONIC_Ring = obj.MUT_Num_IONIC_Ring -obj.WT_Num_IONIC_Ring
                obj.Diff_Num_VDW_Ring = obj.MUT_Num_VDW_Ring -obj.WT_Num_VDW_Ring
                obj.Diff_Num_PICATION_Ring = obj.MUT_Num_PICATION_Ring - obj.WT_Num_PICATION_Ring
                obj.Diff_Num_PIPISTACK_Ring = obj.MUT_Num_PIPISTACK_Ring - obj.WT_Num_PIPISTACK_Ring

                obj.Diff_Num_HBOND_Ring_Layer1 = obj.MUT_Num_HBOND_Ring_Layer1 - obj.WT_Num_HBOND_Ring_Layer1
                obj.Diff_Num_SSBOND_Ring_Layer1 = obj.MUT_Num_SSBOND_Ring_Layer1 - obj.WT_Num_SSBOND_Ring_Layer1
                obj.Diff_Num_IONIC_Ring_Layer1 = obj.MUT_Num_IONIC_Ring_Layer1 -obj.WT_Num_IONIC_Ring_Layer1
                obj.Diff_Num_VDW_Ring_Layer1 = obj.MUT_Num_VDW_Ring_Layer1 -obj.WT_Num_VDW_Ring_Layer1
                obj.Diff_Num_PICATION_Ring_Layer1 = obj.MUT_Num_PICATION_Ring_Layer1 - obj.WT_Num_PICATION_Ring_Layer1
                obj.Diff_Num_PIPISTACK_Ring_Layer1 = obj.MUT_Num_PIPISTACK_Ring_Layer1 - obj.WT_Num_PIPISTACK_Ring_Layer1

                obj.Diff_Num_HBOND_Ring_Layer2 = obj.MUT_Num_HBOND_Ring_Layer2 - obj.WT_Num_HBOND_Ring_Layer2
                obj.Diff_Num_SSBOND_Ring_Layer2 = obj.MUT_Num_SSBOND_Ring_Layer2 - obj.WT_Num_SSBOND_Ring_Layer2
                obj.Diff_Num_IONIC_Ring_Layer2 = obj.MUT_Num_IONIC_Ring_Layer2 -obj.WT_Num_IONIC_Ring_Layer2
                obj.Diff_Num_VDW_Ring_Layer2 = obj.MUT_Num_VDW_Ring_Layer2 -obj.WT_Num_VDW_Ring_Layer2
                obj.Diff_Num_PICATION_Ring_Layer2 = obj.MUT_Num_PICATION_Ring_Layer2 - obj.WT_Num_PICATION_Ring_Layer2
                obj.Diff_Num_PIPISTACK_Ring_Layer2 = obj.MUT_Num_PIPISTACK_Ring_Layer2 - obj.WT_Num_PIPISTACK_Ring_Layer2

                obj.Diff_Num_HBOND_Ring_Layer3 = obj.MUT_Num_HBOND_Ring_Layer3 - obj.WT_Num_HBOND_Ring_Layer3
                obj.Diff_Num_SSBOND_Ring_Layer3 = obj.MUT_Num_SSBOND_Ring_Layer3 - obj.WT_Num_SSBOND_Ring_Layer3
                obj.Diff_Num_IONIC_Ring_Layer3 = obj.MUT_Num_IONIC_Ring_Layer3 -obj.WT_Num_IONIC_Ring_Layer3
                obj.Diff_Num_VDW_Ring_Layer3 = obj.MUT_Num_VDW_Ring_Layer3 -obj.WT_Num_VDW_Ring_Layer3
                obj.Diff_Num_PICATION_Ring_Layer3 = obj.MUT_Num_PICATION_Ring_Layer3 - obj.WT_Num_PICATION_Ring_Layer3
                obj.Diff_Num_PIPISTACK_Ring_Layer3 = obj.MUT_Num_PIPISTACK_Ring_Layer3 - obj.WT_Num_PIPISTACK_Ring_Layer3


                res_dict=Judge_Bond_of_Ring(obj.WT_Ring_Bond_List, obj.WT_Amino_Acid)
                obj.Is_WT_HBOND = res_dict['HBOND']
                obj.Is_WT_SSBOND = res_dict['SSBOND']
                obj.Is_WT_IONIC = res_dict['IONIC']
                obj.Is_WT_VDW = res_dict['VDW']
                obj.Is_WT_PICATION = res_dict['PICATION']
                obj.Is_WT_PIPISTACK = res_dict['PIPISTACK']

                res_dict = Judge_Bond_of_Ring(obj.MUT_Ring_Bond_List, obj.MUT_Amino_Acid)
                obj.Is_MUT_HBOND = res_dict['HBOND']
                obj.Is_MUT_SSBOND = res_dict['SSBOND']
                obj.Is_MUT_IONIC = res_dict['IONIC']
                obj.Is_MUT_VDW = res_dict['VDW']
                obj.Is_MUT_PICATION = res_dict['PICATION']
                obj.Is_MUT_PIPISTACK = res_dict['PIPISTACK']
    else:
        Log(f'Skipped running Ring3. All feature values from Ring3 have been set to initial values.')
        print(f'Skipped running Ring3. All feature values from Ring3 have been set to initial values.')


    #HD_Cluster_List, Num_HD_Cluster_Protlego
    Log(f'Task {task_count}, ID {obj.ID}: Features Extraction 4: Running Protlego')
    print(f'Task {task_count}, ID {obj.ID}: Features Extraction 4: Running Protlego')
    obj.WT_Num_HD_Cluster_Protlego=Run_Prolego(obj.WT_Structure.PDB_path,obj.WT_HD_Cluster_List,Main_Location)
    obj.WT_Num_HD_Cluster_Protlego_Layer1=  Devide_Res_of_HD_Cluster_by_Layers(obj.WT_HD_Cluster_List,obj.WT_Amino_Acid_List_Layer1)
    obj.WT_Num_HD_Cluster_Protlego_Layer2 = Devide_Res_of_HD_Cluster_by_Layers(obj.WT_HD_Cluster_List,obj.WT_Amino_Acid_List_Layer2)
    obj.WT_Num_HD_Cluster_Protlego_Layer3 = Devide_Res_of_HD_Cluster_by_Layers(obj.WT_HD_Cluster_List,obj.WT_Amino_Acid_List_Layer3)
    obj.WT_Max_HD_Cluster_Area=Get_Max_Area(obj.WT_HD_Cluster_List,obj.WT_Amino_Acid_List)


    obj.MUT_Num_HD_Cluster_Protlego=Run_Prolego(obj.MUT_Structure.PDB_path,obj.MUT_HD_Cluster_List,Main_Location)
    obj.MUT_Num_HD_Cluster_Protlego_Layer1=  Devide_Res_of_HD_Cluster_by_Layers(obj.MUT_HD_Cluster_List,obj.MUT_Amino_Acid_List_Layer1)
    obj.MUT_Num_HD_Cluster_Protlego_Layer2 = Devide_Res_of_HD_Cluster_by_Layers(obj.MUT_HD_Cluster_List,obj.MUT_Amino_Acid_List_Layer2)
    obj.MUT_Num_HD_Cluster_Protlego_Layer3 = Devide_Res_of_HD_Cluster_by_Layers(obj.MUT_HD_Cluster_List,obj.MUT_Amino_Acid_List_Layer3)
    obj.MUT_Max_HD_Cluster_Area=Get_Max_Area(obj.MUT_HD_Cluster_List,obj.MUT_Amino_Acid_List)


    obj.Diff_Num_HD_Cluster_Protlego=obj.MUT_Num_HD_Cluster_Protlego-obj.WT_Num_HD_Cluster_Protlego
    obj.Diff_Num_HD_Cluster_Protlego_Layer1 = obj.MUT_Num_HD_Cluster_Protlego_Layer1-obj.WT_Num_HD_Cluster_Protlego_Layer1
    obj.Diff_Num_HD_Cluster_Protlego_Layer2 = obj.MUT_Num_HD_Cluster_Protlego_Layer2-obj.WT_Num_HD_Cluster_Protlego_Layer2
    obj.Diff_Num_HD_Cluster_Protlego_Layer3 = obj.MUT_Num_HD_Cluster_Protlego_Layer3-obj.WT_Num_HD_Cluster_Protlego_Layer3
    obj.Diff_Max_HD_Cluster_Area=obj.MUT_Max_HD_Cluster_Area-obj.WT_Max_HD_Cluster_Area



    res=Judge_If_in_Cluster(obj.WT_Amino_Acid,obj.WT_HD_Cluster_List)
    obj.Is_WT_HD_Cluster=res[0]
    obj.WT_HD_Cluster_Area=res[1]
    res=Judge_If_in_Cluster(obj.MUT_Amino_Acid,obj.MUT_HD_Cluster_List)
    obj.Is_MUT_HD_Cluster=res[0]
    obj.MUT_HD_Cluster_Area=res[1]
    obj.Diff_HD_Cluster_Area=obj.MUT_HD_Cluster_Area-obj.WT_HD_Cluster_Area




    # Amino_Acid_Categories
    Log(f'Task {task_count}, ID {obj.ID}: Features Extraction 5: Calculating AA categories and Running DSSP to get RSA')
    print(f'Task {task_count}, ID {obj.ID}: Features Extraction 5: Calculating AA categories and Running DSSP to get RSA')
    Compute_AA_Categories(obj.WT_Amino_Acid_List,obj.WT_Pct_Amino_Acid_Categories,obj.WT_Num_Amino_Acid_Categories)
    Compute_AA_Categories(obj.WT_Amino_Acid_List_Layer1,obj.WT_Pct_Amino_Acid_Categories_Layer1,obj.WT_Num_Amino_Acid_Categories_Layer1)
    Compute_AA_Categories(obj.WT_Amino_Acid_List_Layer2,obj.WT_Pct_Amino_Acid_Categories_Layer2,obj.WT_Num_Amino_Acid_Categories_Layer2)
    Compute_AA_Categories(obj.WT_Amino_Acid_List_Layer3,obj.WT_Pct_Amino_Acid_Categories_Layer3,obj.WT_Num_Amino_Acid_Categories_Layer3)

    if shutil.which("dssp"):
        res_list=Run_Dssp(obj.Is_Beta,obj.Dssp_List,obj.WT_Structure.PDB_Name, obj.WT_Structure.PDB_path,obj.WT_Seq)
        if res_list is False:
            error_obj.Something_Wrong(Detail_Extraction.__name__)
        else:
            obj.WT_Pct_Buried_Residue=res_list[0]
            obj.WT_Pct_Exposed_Residue=res_list[1]
            obj.WT_Pct_Secondary_Structure=res_list[2]

            res_list=Devide_Res_of_DSSP_by_Layers(obj.Dssp_List,obj.WT_Amino_Acid_List_Layer1,obj.WT_Pct_Secondary_Structure_Layer1)
            obj.WT_Pct_Buried_Residue_Layer1=res_list[0]
            obj.WT_Pct_Exposed_Residue_Layer1 = res_list[1]

            res_list=Devide_Res_of_DSSP_by_Layers(obj.Dssp_List,obj.WT_Amino_Acid_List_Layer2,obj.WT_Pct_Secondary_Structure_Layer2)
            obj.WT_Pct_Buried_Residue_Layer2=res_list[0]
            obj.WT_Pct_Exposed_Residue_Layer2 = res_list[1]

            res_list=Devide_Res_of_DSSP_by_Layers(obj.Dssp_List,obj.WT_Amino_Acid_List_Layer3,obj.WT_Pct_Secondary_Structure_Layer3)
            obj.WT_Pct_Buried_Residue_Layer3=res_list[0]
            obj.WT_Pct_Exposed_Residue_Layer3 = res_list[1]


            res_list=Get_Res_of_DSSP(obj.Is_Beta,obj.WT_Structure.PDB_Name,obj.WT_Structure.PDB_path,obj.WT_Seq,obj.WT_Amino_Acid)
            obj.WT_RSA=res_list[0]
            obj.WT_Is_Buried_or_Exposed=res_list[1]
            obj.WT_Secondary_Structure=res_list[2]
            obj.WT_Secondary_Structure_Char=res_list[3]
            obj.WT_Psi=res_list[4]
            obj.WT_Phi=res_list[5]


            res_list = Get_Res_of_DSSP(obj.Is_Beta,obj.MUT_Structure.PDB_Name, obj.MUT_Structure.PDB_path, obj.MUT_Seq, obj.MUT_Amino_Acid)
            obj.MUT_RSA = res_list[0]
            obj.MUT_Is_Buried_or_Exposed = res_list[1]
            obj.MUT_Secondary_Structure = res_list[2]
            obj.MUT_Secondary_Structure_Char=res_list[3]
            obj.MUT_Psi=res_list[4]
            obj.MUT_Phi=res_list[5]



            obj.Diff_RSA=obj.MUT_RSA-obj.WT_RSA
            obj.Diff_Psi=obj.MUT_Psi-obj.WT_Psi
            obj.Diff_Phi=obj.MUT_Phi-obj.WT_Phi
    else:
        Log(f'Skipped running DSSP. All feature values from DSSP have been set to initial values.')
        print(f'Skipped running DSSP. All feature values from DSSP have been set to initial values.')




    # Pharmacophore
    Log(f'Task {task_count}, ID {obj.ID}: Features Extraction 6: Running Rdkit to get Pharmacophore info')
    print(f'Task {task_count}, ID {obj.ID}: Features Extraction 6: Running Rdkit to get Pharmacophore info')
    is_bonding_WT=Check_Available_PDB_with_Rdkit(obj.WT_Structure.PDB_path)
    is_bonding_MUT = Check_Available_PDB_with_Rdkit(obj.MUT_Structure.PDB_path)
    if is_bonding_WT and is_bonding_MUT:
        is_bonding=True
    else:
        is_bonding=False

    if not Run_Rdikit(obj.WT_Structure.PDB_path, Rdkit_Path, Rdkit_Fdef_Name, obj.WT_Num_Pharmacophore_Categories,obj.WT_Amino_Acid,0.0,is_bonding):
        error_obj.Something_Wrong(Detail_Extraction.__name__,'rkit can not read pdb')
    else:
        if not Run_Rdikit(obj.MUT_Structure.PDB_path, Rdkit_Path, Rdkit_Fdef_Name, obj.MUT_Num_Pharmacophore_Categories,obj.MUT_Amino_Acid,0.0,is_bonding):
            error_obj.Something_Wrong(Detail_Extraction.__name__,'rkit can not read pdb')
        else:
            if not Subtract_Dict(obj.WT_Num_Pharmacophore_Categories,obj.MUT_Num_Pharmacophore_Categories,obj.Diff_Num_Pharmacophore_Categories):
                error_obj.Something_Wrong(Detail_Extraction.__name__)


    if not Run_Rdikit(obj.WT_Structure.PDB_path, Rdkit_Path, Rdkit_Fdef_Name, obj.WT_Num_Pharmacophore_Categories_Layer1,obj.WT_Amino_Acid,obj.Cutoff1,is_bonding):
        error_obj.Something_Wrong(Detail_Extraction.__name__,'rkit can not read pdb')
    else:
        if not Run_Rdikit(obj.MUT_Structure.PDB_path, Rdkit_Path, Rdkit_Fdef_Name, obj.MUT_Num_Pharmacophore_Categories_Layer1,obj.MUT_Amino_Acid,obj.Cutoff1,is_bonding):
            error_obj.Something_Wrong(Detail_Extraction.__name__,'rkit can not read pdb')
        else:
            if not Subtract_Dict(obj.WT_Num_Pharmacophore_Categories_Layer1,obj.MUT_Num_Pharmacophore_Categories_Layer1,obj.Diff_Num_Pharmacophore_Categories_Layer1):
                error_obj.Something_Wrong(Detail_Extraction.__name__)


    if not Run_Rdikit(obj.WT_Structure.PDB_path, Rdkit_Path, Rdkit_Fdef_Name, obj.WT_Num_Pharmacophore_Categories_Layer2,obj.WT_Amino_Acid,obj.Cutoff2,is_bonding):
        error_obj.Something_Wrong(Detail_Extraction.__name__,'rkit can not read pdb')
    else:
        if not Run_Rdikit(obj.MUT_Structure.PDB_path, Rdkit_Path, Rdkit_Fdef_Name, obj.MUT_Num_Pharmacophore_Categories_Layer2,obj.MUT_Amino_Acid,obj.Cutoff2,is_bonding):
            error_obj.Something_Wrong(Detail_Extraction.__name__,'rkit can not read pdb')
        else:
            if not Subtract_Dict(obj.WT_Num_Pharmacophore_Categories_Layer2,obj.MUT_Num_Pharmacophore_Categories_Layer2,obj.Diff_Num_Pharmacophore_Categories_Layer2):
                error_obj.Something_Wrong(Detail_Extraction.__name__)


    if not Run_Rdikit(obj.WT_Structure.PDB_path, Rdkit_Path, Rdkit_Fdef_Name, obj.WT_Num_Pharmacophore_Categories_Layer3,obj.WT_Amino_Acid,obj.Cutoff3,is_bonding):
        error_obj.Something_Wrong(Detail_Extraction.__name__,'rkit can not read pdb')
    else:
        if not Run_Rdikit(obj.MUT_Structure.PDB_path, Rdkit_Path, Rdkit_Fdef_Name, obj.MUT_Num_Pharmacophore_Categories_Layer3,obj.MUT_Amino_Acid,obj.Cutoff3,is_bonding):
            error_obj.Something_Wrong(Detail_Extraction.__name__,'rkit can not read pdb')
        else:
            if not Subtract_Dict(obj.WT_Num_Pharmacophore_Categories_Layer3,obj.MUT_Num_Pharmacophore_Categories_Layer3,obj.Diff_Num_Pharmacophore_Categories_Layer3):
                error_obj.Something_Wrong(Detail_Extraction.__name__)







    #B Factor
    Log(f'Task {task_count}, ID {obj.ID}: Features Extraction 7: Running PROFbval in container to get B-factor')
    print(f'Task {task_count}, ID {obj.ID}: Features Extraction 7: Running PROFbval in container to get B-factor')
    if Scripts.Global_Value.D_or_S != '-':
        if (Scripts.Global_Value.D_or_S=='S' and os.path.isfile(Singularity_Container_Path)) or (Scripts.Global_Value.D_or_S=='D' and os.path.isfile('./src/Prof_Source/myprof.tar')):
            res=Compute_B_Factor(obj.WT_Seq,obj.Chain_ID_of_Mut,TMP_Path,f'prof_res_{obj.ID}_WT',obj.WT_PSI_BLAST_Path,Main_Location,obj.True_Loc_of_Mutation,obj.WT_Amino_Acid_short)
            if res is False:
                error_obj.Something_Wrong(Detail_Extraction.__name__)
            else:
                obj.WT_B_Factor=res

                res = Compute_B_Factor(obj.MUT_Seq, obj.Chain_ID_of_Mut,TMP_Path,f'prof_res_{obj.ID}_MUT', obj.MUT_PSI_BLAST_Path,
                                   Main_Location, obj.True_Loc_of_Mutation, obj.MUT_Amino_Acid_short)
                if res is False:
                    error_obj.Something_Wrong(Detail_Extraction.__name__)
                else:
                    obj.MUT_B_Factor = res
                    obj.Diff_B_Factor=obj.MUT_B_Factor-obj.WT_B_Factor
        else:
            Log(f'Skipped running PROFbval. All feature values from PROFbval have been set to initial values.')
            print(f'Skipped running PROFbval. All feature values from PROFbval have been set to initial values.')
    else:
        Log(f'Skipped running PROFbval. All feature values from PROFbval have been set to initial values.')
        print(f'Skipped running PROFbval. All feature values from PROFbval have been set to initial values.')


    #FoldX
    Log(f'Task {task_count}, ID {obj.ID}: Features Extraction 8: Running FoldX')
    print(f'Task {task_count}, ID {obj.ID}: Features Extraction 8: Running FoldX')
    if os.path.isfile(FoldX_Path+FoldX_Name):
        if not Run_FoldX(FoldX_Path,FoldX_Name,obj.WT_Structure.PDB_path,obj.WT_Amino_Acid_short,obj.MUT_Amino_Acid_short,obj.True_Loc_of_Mutation,obj.Chain_ID_of_Mut,obj.WT_FoldX_Energy_Term_Dict,obj.Diff_FoldX_Energy_Term_Dict,TMP_Path,f'foldx_res_{obj.ID}'):
            error_obj.Something_Wrong(Detail_Extraction.__name__)
    else:
        Log(f'Skipped running FoldX. All feature values from FoldX have been set to initial values.')
        print(f'Skipped running FoldX. All feature values from FoldX have been set to initial values.')

    #NMA
    Log(f'Task {task_count}, ID {obj.ID}: Features Extraction 9: Running Bio3D to get NMA')
    print(f'Task {task_count}, ID {obj.ID}: Features Extraction 9: Running Bio3D to get NMA')
    if shutil.which("Rscript"):
        res=Run_NMA(obj.WT_Structure.PDB_path,obj.MUT_Structure.PDB_path,obj.True_Loc_of_Mutation,obj.Chain_ID_of_Mut,obj.WT_Seq,R_NMA_Path,R_NMA_App_Name,TMP_Path,f'nma_res_{obj.ID}')
        if res is False:
            error_obj.Something_Wrong(Detail_Extraction.__name__,'NMA failed, maybe PDB is too big and need more memories')
        else:
            obj.WT_NMA_Fluctuation=res['wt_fluctuation_loc']
            obj.MUT_NMA_Fluctuation=res['mut_fluctuation_loc']
            obj.Overall_Rmsip=res['rmsip']
            obj.Diff_NMA_Fluctuation=obj.MUT_NMA_Fluctuation-obj.WT_NMA_Fluctuation
    else:
        Log(f'Skipped running Bio3D. All feature values from Bio3D have been set to initial values.')
        print(f'Skipped running Bio3D. All feature values from Bio3D have been set to initial values.')

    #Length
    Log(f'Task {task_count}, ID {obj.ID}: Features Extraction 10: Running DisEMBL')
    print(f'Task {task_count}, ID {obj.ID}: Features Extraction 10: Running DisEMBL')
    if os.path.isfile(DisEMBL_Path+'DisEMBL.py'):
        res=Run_DisEMBL(obj.WT_Seq,obj.Chain_ID_of_Mut,obj.WT_Structure.PDB_Name,DisEMBL_Path,TMP_Path,f'disembl_res_{obj.ID}')
        if res is False:
            error_obj.Something_Wrong(Detail_Extraction.__name__)
        else:
            obj.COILS_line=res[0]
            obj.REM465_line=res[1]
            obj.HOTLOOPS_line=res[2]
            res=Generate_Res_DisEMBL(obj.COILS_line,obj.REM465_line,obj.HOTLOOPS_line,obj.Chain_ID_of_Mut,obj.WT_Amino_Acid_List)
            if res is False:
                error_obj.Something_Wrong(Detail_Extraction.__name__)
            else:
                obj.WT_Pct_coils=res['COILS_Pct']
                obj.WT_Whole_Length_coils = res['COILS_Length']
                obj.WT_Pct_rem465 = res['REM465_Pct']
                obj.WT_Whole_Length_rem465 = res['REM465_Length']
                obj.WT_Pct_hotloop = res['HOTLOOPS_Pct']
                obj.WT_Whole_Length_hotloop = res['HOTLOOPS_Length']
    else:
        Log(f'Skipped running DisEMBL. All feature values from DisEMBL have been set to initial values.')
        print(f'Skipped running DisEMBL. All feature values from DisEMBL have been set to initial values.')



    #SIFT
    Log(f'Task {task_count}, ID {obj.ID}: Features Extraction 11: Running SIFT')
    print(f'Task {task_count}, ID {obj.ID}: Features Extraction 11: Running SIFT')
    if os.path.exists(SIFT_Path + 'blimps/') and os.path.exists(SIFT_Path+'bin/'):
        res=Run_Sift(obj.WT_Structure.PDB_Name,obj.WT_Amino_Acid_short,obj.MUT_Amino_Acid_short,obj.True_Loc_of_Mutation,SIFT_Path,WT_MSA_Path,obj.WT_Seq,obj.Chain_ID_of_Mut,obj.WT_BLASTP_Path,TMP_Path,f'sift_res_{obj.ID}')
        if res is False:
            error_obj.Something_Wrong(Detail_Extraction.__name__)
        else:
            obj.SIFT_Score=res
    else:
        Log(f'Skipped running SIFT. All feature values from SIFT have been set to initial values.')
        print(f'Skipped running SIFT. All feature values from SIFT have been set to initial values.')



    #
    Log(f'Task {task_count}, ID {obj.ID}: Features Extraction 12: Calculating features on AA site')
    print(f'Task {task_count}, ID {obj.ID}: Features Extraction 12: Calculating features on AA site')
    try:
        res_list=Get_Mutation_Description(obj.WT_Amino_Acid,obj.MUT_Amino_Acid,obj.WT_Secondary_Structure_Char,obj.MUT_Secondary_Structure_Char)
        obj.WT_AA_Type=res_list[0]
        obj.MUT_AA_Type=res_list[1]
        obj.Descri_AA=res_list[2]
        obj.Descri_SS=res_list[3]

        res_dict=Judge_AA_Categories(obj.WT_Amino_Acid)
        obj.Is_WT_Uncharged_Polar=res_dict['uncharged_polar']
        obj.Is_WT_Positively_Charged_Polar=res_dict['positively_charged_polar']
        obj.Is_WT_Negatively_Charged_Polar=res_dict['negatively_charged_polar']
        obj.Is_WT_Nonpolar=res_dict['nonpolar']
        obj.Is_WT_Aliphatic=res_dict['aliphatic']
        obj.Is_WT_Aromatic=res_dict['aromatic']
        obj.Is_WT_Heterocyclic=res_dict['heterocyclic']
        obj.Is_WT_Sulfur_Containing=res_dict['sulfur_containing']

        res_dict = Judge_AA_Categories(obj.MUT_Amino_Acid)
        obj.Is_MUT_Uncharged_Polar = res_dict['uncharged_polar']
        obj.Is_MUT_Positively_Charged_Polar = res_dict['positively_charged_polar']
        obj.Is_MUT_Negatively_Charged_Polar = res_dict['negatively_charged_polar']
        obj.Is_MUT_Nonpolar = res_dict['nonpolar']
        obj.Is_MUT_Aliphatic = res_dict['aliphatic']
        obj.Is_MUT_Aromatic = res_dict['aromatic']
        obj.Is_MUT_Heterocyclic = res_dict['heterocyclic']
        obj.Is_MUT_Sulfur_Containing = res_dict['sulfur_containing']

        obj.Descri_HBOND = Return_4_type(obj.Is_WT_HBOND,obj.Is_MUT_HBOND)
        obj.Descri_SSBOND = Return_4_type(obj.Is_WT_SSBOND,obj.Is_MUT_SSBOND)
        obj.Descri_IONIC = Return_4_type(obj.Is_WT_IONIC,obj.Is_MUT_IONIC)
        obj.Descri_VDW = Return_4_type(obj.Is_WT_VDW,obj.Is_MUT_VDW)
        obj.Descri_PICATION = Return_4_type(obj.Is_WT_PICATION,obj.Is_MUT_PICATION)
        obj.Descri_PIPISTACK = Return_4_type(obj.Is_WT_PIPISTACK,obj.Is_WT_PIPISTACK)

        obj.Descri_HD_Cluster = Return_4_type(obj.Is_WT_HD_Cluster,obj.Is_MUT_HD_Cluster)

        obj.Descri_Buried_or_Exposed = Return_4_type(obj.WT_Is_Buried_or_Exposed,obj.MUT_Is_Buried_or_Exposed)

        obj.Descri_Uncharged_Polar = Return_4_type(obj.Is_WT_Uncharged_Polar,obj.Is_MUT_Uncharged_Polar)
        obj.Descri_Positively_Charged_Polar = Return_4_type(obj.Is_WT_Positively_Charged_Polar,obj.Is_MUT_Positively_Charged_Polar)
        obj.Descri_Negatively_Charged_Polar = Return_4_type(obj.Is_WT_Negatively_Charged_Polar,obj.Is_WT_Negatively_Charged_Polar)
        obj.Descri_Nonpolar = Return_4_type(obj.Is_WT_Nonpolar,obj.Is_MUT_Nonpolar)
        obj.Descri_Aliphatic = Return_4_type(obj.Is_WT_Aliphatic,obj.Is_MUT_Aliphatic)
        obj.Descri_Aromatic = Return_4_type(obj.Is_WT_Aromatic,obj.Is_MUT_Aromatic)
        obj.Descri_Heterocyclic = Return_4_type(obj.Is_WT_Heterocyclic,obj.Is_MUT_Heterocyclic)
        obj.Descri_Sulfur_Containing = Return_4_type(obj.Is_WT_Sulfur_Containing,obj.Is_MUT_Sulfur_Containing)
    except:
        pass


    res_list=find_pssm_score(obj.WT_PSSM_Path,obj.WT_Amino_Acid_List,obj.WT_Amino_Acid,obj.WT_Seq,obj.Chain_ID_of_Mut,5)
    if res_list is False:
        error_obj.Something_Wrong(Detail_Extraction.__name__)
    else:
        obj.WT_PSSM_Score=res_list[5]
        obj.WT_PSSM_Score_F1=res_list[0]
        obj.WT_PSSM_Score_F2 = res_list[1]
        obj.WT_PSSM_Score_F3 = res_list[2]
        obj.WT_PSSM_Score_F4 = res_list[3]
        obj.WT_PSSM_Score_F5 = res_list[4]
        obj.WT_PSSM_Score_B1 = res_list[6]
        obj.WT_PSSM_Score_B2 = res_list[7]
        obj.WT_PSSM_Score_B3 = res_list[8]
        obj.WT_PSSM_Score_B4 = res_list[9]
        obj.WT_PSSM_Score_B5 = res_list[10]
        obj.WT_PSSM_Score_Aver = res_list[11]

        res_list = find_pssm_score(obj.MUT_PSSM_Path, obj.MUT_Amino_Acid_List, obj.MUT_Amino_Acid,obj.MUT_Seq,obj.Chain_ID_of_Mut,5)
        if res_list is False:
            error_obj.Something_Wrong(Detail_Extraction.__name__)
        else:
            obj.MUT_PSSM_Score = res_list[5]
            obj.MUT_PSSM_Score_F1 = res_list[0]
            obj.MUT_PSSM_Score_F2 = res_list[1]
            obj.MUT_PSSM_Score_F3 = res_list[2]
            obj.MUT_PSSM_Score_F4 = res_list[3]
            obj.MUT_PSSM_Score_F5 = res_list[4]
            obj.MUT_PSSM_Score_B1 = res_list[6]
            obj.MUT_PSSM_Score_B2 = res_list[7]
            obj.MUT_PSSM_Score_B3 = res_list[8]
            obj.MUT_PSSM_Score_B4 = res_list[9]
            obj.MUT_PSSM_Score_B5 = res_list[10]
            obj.MUT_PSSM_Score_Aver = res_list[11]


            obj.Diff_PSSM_Score=obj.MUT_PSSM_Score-obj.WT_PSSM_Score
            obj.Diff_PSSM_Score_Aver=obj.MUT_PSSM_Score_Aver-obj.WT_PSSM_Score_Aver


    Log(f'Task {task_count}, ID {obj.ID}: Features Extraction 13: Calculating AAindex features')
    print(f'Task {task_count}, ID {obj.ID}: Features Extraction 13: Calculating AAindex features')
    obj.WT_AAindex1=Get_Mutation_Index_List_from_Index(obj.WT_Amino_Acid_short,obj.MUT_Amino_Acid_short,Scripts.AAindex.aaindex1_list,obj.Diff_AAindex1)
    Get_Mutation_Index_List_from_Matrix(f'{obj.WT_Amino_Acid_short}{obj.MUT_Amino_Acid_short}',Scripts.AAindex.aaindex2_list,obj.Descri_AAindex2)
    Get_Mutation_Index_List_from_Matrix(f'{obj.WT_Amino_Acid_short}{obj.MUT_Amino_Acid_short}', Scripts.AAindex.aaindex3_list,obj.Descri_AAindex3)


    Log(f'Task {task_count}, ID {obj.ID}: Have finished features extraction, and normally return')
    print(f'Task {task_count}, ID {obj.ID}: Have finished features extraction, and normally return')
    #return True
    return obj










