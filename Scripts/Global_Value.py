import os

Is_Use_Reverse_Data=True

BLAST_Process_Num=4

Process_Num=1

Singularity_Container_Path='./src/Prof_Source/myprof.sif'

D_or_S='-'

Mode=''

Docker_Container_Name='myprof'
Docker_Image_ID='7d9fe6723898'

Raw_Dataset_file=''


Raw_PDB_Path='./src/Raw_PDB/'
WT_PDB_Path='./src/WT_PDB/'
MUT_PDB_Path='./src/Mut_PDB/'
WT_Fasta_Path='./src/WT_Fasta/'
MUT_Fasta_Path='./src/Mut_Fasta/'
WT_PSSM_Data_Path='./src/PSSM_Data/WT/'
MUT_PSSM_Data_Path='./src/PSSM_Data/Mut/'
WT_PSI_BLAST_Data_Path='./src/PSI_BLAST_Data/WT/'
MUT_PSI_BLAST_Data_Path='./src/PSI_BLAST_Data/Mut/'
WT_BLASTP_Data_Path='./src/BLASTP_Data/WT/'
MUT_BLASTP_Data_Path='./src/BLASTP_Data/Mut/'

TMP_Path='./src/TMP/'

Table_Path='./src/Data_Table/'
Res_Table_Name='data_table.txt'

Pred_Table_Path='./src/Pred_Table/'
Pred_Table_Name='pred_table.txt'
Pred_Res_Path='./src/Pred_Res/'

Is_Pred=0



Ring_Path='./bin/ring-3.0.0/ring/bin/'
Ring_Expired_Date='2033-05-04'


FoldX_Path='./bin/FoldX_5.0/'

FoldX_Name='foldx_20241231'
FoldX_Expired_Date='2024-12-31'
Rotabase_Path='./src/Others/rotabase.txt'

Rdkit_Path='./bin/rdkit_2023_3_1/'

Rdkit_Fdef_Name='BaseFeatures.fdef'


Features_Table_Path='./src/Feature_Res/'
Features_Table_Name='features_table.csv'

Model_Path='./Models/'





MSA_DB_Path=''
MSA_DB_Name=''




Main_Location=os.path.abspath('./')+'/'
Home_Location=os.path.expandvars('$HOME')

R_NMA_Path='./bin/R_NMA/'
R_NMA_App_Name='NMA.R'

DisEMBL_Path='./bin/DisEMBL_1_4/'

BLAST_Path='./bin/ncbi_blast_2_13_0+/bin/'


# WT_MSA_Path='./src/WT_MSA/'
WT_MSA_Path=''
SIFT_Path='./bin/sift6_2_1/'




Clean_Path='./bin/clean/'

AAIndex1_Path='./src/AAindex/aaindex1'
AAIndex2_Path='./src/AAindex/aaindex2'
AAIndex3_Path='./src/AAindex/aaindex3'

Log_Path='./src/log.txt'
