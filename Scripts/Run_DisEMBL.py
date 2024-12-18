try:
    from bin.DisEMBL_1_4.DisEMBL import runDisEMBLpipeline
except:
    pass
# from Utils import Clean_Main_Directory
from Scripts.Error import error_obj
from Scripts.Classes import Researched_Amino_Acid

def Run_DisEMBL(seq_dict:dict,chain_id_mut,name,path,temp_path,o_folder_name):
    '''
    :purpose: Call DisEMBL to compute COILS, REM465 and  HOTLOOPS info
    :param seq_dict: Input sequence dict
    :param chain_id_mut: Input chain id
    :param name: name
    :param path:Input path of DisEMBL
    :param temp_path: TMP Path
    :param o_folder_name: output folder name
    :outpath: like TMP/disembl_res_ID
    :return: A list of info/False
    :process: 1. make outpath
              2. make fasta file
              3. Run DisEMBL and output in outpath
              4. Read result in outpath
    '''
    outpath=temp_path+o_folder_name+'/'
    import os,shutil
    if os.path.exists(outpath):
        shutil.rmtree(outpath)
    os.mkdir(outpath)
    try:
        with open(f'{outpath}/temp.fasta','w') as temp:
            temp.write(f'>{name}_{chain_id_mut}\n')
            temp.write(seq_dict[chain_id_mut])
        runDisEMBLpipeline(8, 8, 4, 1.2, 1.4, 1.2, f'{outpath}/temp.fasta', path,outpath)
        with open(f'{outpath}/out.txt','r') as out:
            lines=out.readlines()
            count=0
            COILS_line=''
            REM465_line=''
            HOTLOOPS_line=''
            for line in lines:
                if line.find('>')!=-1 and line.find('COILS')!=-1:
                    count+=1
                    COILS_line=line
                if line.find('>')!=-1 and line.find('REM465')!=-1:
                    count+=1
                    REM465_line=line
                if line.find('>')!=-1 and line.find('HOTLOOPS')!=-1:
                    count+=1
                    HOTLOOPS_line=line
            if count!=3 or COILS_line=='' or REM465_line=='' or HOTLOOPS_line=='':
                error_obj.Something_Wrong(Run_DisEMBL.__name__)
                return False
            # Clean_Main_Directory()
            return [COILS_line,REM465_line,HOTLOOPS_line]
    except:
        error_obj.Something_Wrong(Run_DisEMBL.__name__)
        # Clean_Main_Directory()
        return False


def Get_Length_by_Line(line:str,aa_list:list[Researched_Amino_Acid]):
    try:
        a_l=[]
        for aa in aa_list:
            a_l.append(aa.Num)
        l=line.split()
        l_=[]
        for i in l[2:]:
            l_.append(i.replace('\n','').replace(',','').replace(' ',''))
        if len(l_)==0:
            return 0
        whole_length=0
        for i in l_:
            div=i.split('-')
            begin=int(div[0])
            end=int(div[1])
            l__=range(begin,end+1)
            for j in l__:
                if j in a_l:
                    whole_length+=1
        return whole_length
    except:
        error_obj.Something_Wrong(Get_Length_by_Line.__name__)
        return False

def Generate_Res_DisEMBL(COILS_line,REM465_line,HOTLOOPS_line,chain_id_mut,aa_list:list[Researched_Amino_Acid]):
    COILS_length = Get_Length_by_Line(COILS_line, aa_list)
    REM465_length = Get_Length_by_Line(REM465_line, aa_list)
    HOTLOOPS_length = Get_Length_by_Line(HOTLOOPS_line, aa_list)
    whole_length = 0
    for aa in aa_list:
        if aa.Chain_ID == chain_id_mut:
            whole_length += 1
    length_dict = {'COILS_Length': COILS_length, 'COILS_Pct': float(COILS_length / whole_length),
                   'REM465_Length': REM465_length, 'REM465_Pct': float(REM465_length / whole_length),
                   'HOTLOOPS_Length': HOTLOOPS_length, 'HOTLOOPS_Pct': float(HOTLOOPS_length / whole_length)}
    return length_dict