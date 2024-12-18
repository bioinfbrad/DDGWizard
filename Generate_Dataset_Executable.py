import os
import warnings
warnings.filterwarnings("ignore")
from Scripts.Log import Log, Init_Log

from Scripts.Feature_Extracting import *
from Scripts.Run_Modeller import *
from Scripts.MSA import *
from Scripts.Record import Record_Feature_Table
from Scripts.Global_Value import *
import Scripts.Global_Value
from Scripts.Init import Init
from Scripts.Utils import *
import argparse
from Scripts.Docker import Docker_Init_Container,Docker_Remove_Container


if __name__ == '__main__':
    Init_Log()
    parser = argparse.ArgumentParser(description='Input arguments')

    parser.add_argument('--raw_dataset_path', type=str, default='')
    parser.add_argument('--db_folder_path', type=str, default='')
    parser.add_argument('--db_name', type=str, default='')
    parser.add_argument('--if_reversed_data', type=int, default=1)
    parser.add_argument('--blast_process_num', type=int, default=1)
    parser.add_argument('--container_type', type=str, default='-')
    parser.add_argument('--mode',type=str,default='whole')
    parser.add_argument('--process_num', type=int, default=1)

    current_directory = os.getcwd()
    script_directory = os.path.dirname(os.path.abspath(__file__))
    if current_directory!=script_directory:
        error_obj.Something_Wrong(__name__,"Please cd to top folder of program!!!")
        exit(1)

    Log('Processing input arguments')
    print('Processing input arguments')
    args = parser.parse_args()
    if args.raw_dataset_path=='' or args.db_folder_path=='' or args.db_name=='' or args.if_reversed_data not in [0,1] or args.blast_process_num<1 or args.blast_process_num>200 or args.process_num>200 or args.process_num<1:
        error_obj.Something_Wrong(__name__)
        exit(1)
    if str(args.raw_dataset_path).split('.')[len(str(args.raw_dataset_path).split('.'))-1]!='csv':
        error_obj.Something_Wrong(__name__)
        exit(1)
    if not os.path.exists(args.raw_dataset_path):
        error_obj.Something_Wrong(__name__)
        exit(1)
    if not os.path.isdir(args.db_folder_path):
        error_obj.Something_Wrong(__name__)
        exit(1)
    if str(args.container_type) not in ['D','S','-']:
        error_obj.Something_Wrong(__name__)
        exit(1)
    if str(args.mode) not in ['blast_only','model_only','whole']:
        error_obj.Something_Wrong(__name__)
        exit(1)

    Scripts.Global_Value.Raw_Dataset_file=args.raw_dataset_path
    Scripts.Global_Value.MSA_DB_Path=args.db_folder_path
    Scripts.Global_Value.MSA_DB_Name=args.db_name
    Scripts.Global_Value.Is_Use_Reverse_Data=args.if_reversed_data
    Scripts.Global_Value.BLAST_Process_Num=args.blast_process_num
    Scripts.Global_Value.D_or_S = args.container_type
    Scripts.Global_Value.Mode = args.mode
    Scripts.Global_Value.Process_Num = args.process_num

    Log(f'Your input arguments:\n--raw_dataset_path:{Scripts.Global_Value.Raw_Dataset_file}\n--db_folder_path:{Scripts.Global_Value.MSA_DB_Path}\n--db_name:{Scripts.Global_Value.MSA_DB_Name}\n--if_reversed_data:{Scripts.Global_Value.Is_Use_Reverse_Data}\n--blast_process_num:{Scripts.Global_Value.BLAST_Process_Num}\n--container_type:{Scripts.Global_Value.D_or_S}\n--mode:{Scripts.Global_Value.Mode}\n--process_num:{Scripts.Global_Value.Process_Num}\n')
    print(f'Your input arguments:\n--raw_dataset_path:{Scripts.Global_Value.Raw_Dataset_file}\n--db_folder_path:{Scripts.Global_Value.MSA_DB_Path}\n--db_name:{Scripts.Global_Value.MSA_DB_Name}\n--if_reversed_data:{Scripts.Global_Value.Is_Use_Reverse_Data}\n--blast_process_num:{Scripts.Global_Value.BLAST_Process_Num}\n--container_type:{Scripts.Global_Value.D_or_S}\n--mode:{Scripts.Global_Value.Mode}\n--process_num:{Scripts.Global_Value.Process_Num}\n')

    if Scripts.Global_Value.D_or_S=='D':
        Log('Initing Docker')
        print('Initing Docker')
        Docker_Init_Container(Docker_Container_Name,Docker_Image_ID)

    try:
        Log('Initing configuration')
        print('Initing configuration')
        Init()

        Log('Reading raw dataset ')
        print('Reading raw dataset ')
        Raw_Data_List = Read_CSV(Scripts.Global_Value.Raw_Dataset_file)

        Log('Clearing folders')
        print('Clearing folders')
        Clean_All_Res_Folder(Table_Path,Features_Table_Path,Raw_PDB_Path,WT_PDB_Path,MUT_PDB_Path,WT_Fasta_Path,MUT_Fasta_Path,WT_PSSM_Data_Path,MUT_PSSM_Data_Path,WT_PSI_BLAST_Data_Path,MUT_PSI_BLAST_Data_Path,WT_BLASTP_Data_Path,MUT_BLASTP_Data_Path,[1,1,0,0,0,1,1,0,0,0,0,0,0])

        Log('Preparing task table')
        print('Preparing task table')
        Prepare_Table(Raw_Data_List,Table_Path,Res_Table_Name,Clean_Path,Raw_PDB_Path,WT_PDB_Path,MUT_PDB_Path,WT_Fasta_Path,MUT_Fasta_Path,WT_PSSM_Data_Path,MUT_PSSM_Data_Path,WT_PSI_BLAST_Data_Path,MUT_PSI_BLAST_Data_Path)


        if Scripts.Global_Value.Mode=='whole' or Scripts.Global_Value.Mode=='model_only':
            Log('Modelling MUT models')
            print('Modelling MUT models')
            Prepare_MUT_Models(Table_Path,Res_Table_Name,MUT_PDB_Path,Scripts.Global_Value.Process_Num)

        if Scripts.Global_Value.Mode=='whole' or Scripts.Global_Value.Mode=='blast_only':
            Log('Preparing Blast files')
            print('Preparing Blast files')
            Prepare_Blast_Files(Table_Path,Res_Table_Name,WT_PSSM_Data_Path,MUT_PSSM_Data_Path,WT_PSI_BLAST_Data_Path,MUT_PSI_BLAST_Data_Path,WT_BLASTP_Data_Path,MUT_BLASTP_Data_Path,Scripts.Global_Value.MSA_DB_Path,Scripts.Global_Value.MSA_DB_Name)

        if Scripts.Global_Value.Mode=='whole' and Scripts.Global_Value.Is_Use_Reverse_Data:
            Log('Adding reverse task')
            print('Adding reverse task')
            Add_Reverse_Data(Table_Path,Res_Table_Name)

        if Scripts.Global_Value.Mode=='whole':
            Feature_Object_List = []

            Log('Beginning features extraction')
            print('Beginning features extraction')
            Feature_Extraction(Table_Path,Res_Table_Name,Feature_Object_List,Scripts.Global_Value.Process_Num)

            Log('Recording features results')
            print('Recording features results')
            if not Record_Feature_Table(Feature_Object_List,Features_Table_Path):
                error_obj.Something_Wrong(__name__)
                exit(1)

    except:
        error_obj.Something_Wrong(__name__)
        Clean_with_Error(Docker_Container_Name)
        exit(1)

    if Scripts.Global_Value.D_or_S=='D':
        Log('Removing Docker container')
        print('Removing Docker container')
        Docker_Remove_Container(Docker_Container_Name)

    Remove_FoldX_Resource()
    Clean_Main_Directory()


    Log('Cleaning temporary folder in ./src/TMP/')
    print('Cleaning temporary folder in ./src/TMP/')
    import shutil
    shutil.rmtree(TMP_Path)
    os.mkdir(TMP_Path)

    Log('Have finished')
    print('Have finished')
    exit(0)













