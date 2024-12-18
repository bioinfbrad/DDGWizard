#!/usr/bin/env python

import os

from aa_additional import full_aa_names
from aa_additional import modified_aa

import argparse

bad_inserted_res = False
bad_alternate_pos = False
bad_modified_res = False
bad_missing_dns = False 

pdb_file = ""

def Check_and_Record_PDB(count, res_buffer):
    global pdb_file
    Is_CA = False
    Is_N = False
    Is_C = False
    for line in res_buffer:
        atom_name = line[12:16]
        occupancy = float(line[55:60])
        if atom_name == " CA " and occupancy > 0.0:
            Is_CA = True
        if atom_name == " N  " and occupancy > 0.0:
            Is_N = True
        if atom_name == " C  " and occupancy > 0.0:
            Is_C = True

    if Is_CA and Is_N and Is_C:
        for line in res_buffer:
            new_num = '%4d ' % count
            line_edit = line[0:22] + new_num + line[27:]
            pdb_file = pdb_file + line_edit
        return True
    return False

def Open_PDB(filename):
    if os.path.exists(filename):
        print(f"pdb file is at {filename}, ready for cleaning")

    file_prefix = os.path.basename(filename)
    pdb_num = file_prefix[:-4]
    lines = open(filename, 'r').readlines()

    return lines, pdb_num

parser = argparse.ArgumentParser()
parser.add_argument('--pdb', type=str, default='')
args = parser.parse_args()

lines, filename_pdb = Open_PDB(args.pdb)

old_res_num = '   '
count = 1

res_buffer = []

for line in lines:

    if line.startswith('ENDMDL'): break  
    if len(line) > 21 :
        if line[0:4] != "ATOM" and line[0:6] != 'HETATM':
            continue

        line_edit = line
        res_name = line[17:20]

        if res_name in modified_aa:
            ori_res_name = res_name
            res_name = modified_aa[res_name]
            line_edit = 'ATOM  '+line[6:17]+ res_name + line[20:]

            if ori_res_name == "MSE":
                if (line_edit[12:14] == 'SE'):
                    line_edit = line_edit[0:12]+' S'+line_edit[14:]
                if len(line_edit) > 75:
                    if (line_edit[76:78] == 'SE'):
                        line_edit = line_edit[0:76]+' S'+line_edit[78:]
            else:
                bad_modified_res = True

        if res_name not in full_aa_names:
            continue

        res_num = line_edit[22:27]

        if not res_num == old_res_num:
            if res_buffer != []:
                if not Check_and_Record_PDB(count, res_buffer):
                    bad_missing_dns = True
                else:
                    count = count + 1

            res_buffer = []

        old_res_num = res_num

        inserted_res = line[26]
        if inserted_res != ' ':
            bad_inserted_res = True

        alternate_pos = line[16]
        if alternate_pos != ' ':
            bad_alternate_pos = True
            if alternate_pos == 'A':
                line_edit = line_edit[:16]+' '+line_edit[17:]
            else:
                continue

        res_buffer.append(line_edit)

if res_buffer != []: 
    if not Check_and_Record_PDB(count, res_buffer):
        bad_missing_dns = True
    else:
        count = count + 1

print(filename_pdb, bad_alternate_pos,  bad_inserted_res,  bad_modified_res,  bad_missing_dns)

outfile = filename_pdb + ".pdb"

outid = open(outfile, 'w')
outid.write(pdb_file)
outid.write("TER\n")
outid.close()


