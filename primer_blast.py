import os
from Bio import SeqIO 
import json
os.chdir('/home/ahmed/Downloads')
with open('sxtF.fasta','w') as f:
    f.write('>SXTF\n')
    f.write('TTATCGTTTCGATGGC')

with open('sxtB.fasta','w') as f:
    f.write('>SXTB\n')
    f.write('GCTCTTCTTGTCCGTTC')
    
#%%
#list of files    
A=os.listdir('/home/ahmed/Lab_HD/fluvialis_genomics/Vibrio_fluvialis/03.assembly_annotation/test_folder/Plasmid/ice_removed')

os.chdir('/home/ahmed/Lab_HD/fluvialis_genomics/Vibrio_fluvialis/03.assembly_annotation/test_folder/plasmid_database')    

for db in A:
    temp_node1='/home/ahmed/Downloads/sxtF.fasta'
    json_file1='/home/ahmed/Downloads/output1.json'
    os.system('blastn -task blastn-short -db %s -query %s -outfmt 15  -out %s '%(db,temp_node1,json_file1))
    D1 = open(json_file1, 'r').read()
    D1 = json.loads(D1)
    ref_rem=[]
    
    temp_node2='/home/ahmed/Downloads/sxtB.fasta'
    json_file2='/home/ahmed/Downloads/output2.json'
    os.system('blastn -task blastn-short -db %s -query %s -outfmt 15  -out %s '%(db,temp_node2,json_file2))
    D2 = open(json_file2, 'r').read()
    D2 = json.loads(D2)
    ref_rem=[]    
    lenD1=len(D1['BlastOutput2'][0]['report']['results']['search']['hits'])
    lenD2=len(D2['BlastOutput2'][0]['report']['results']['search']['hits'])
    a=min([lenD1,lenD2])
    for i in range(a):
        if D1['BlastOutput2'][0]['report']['results']['search']['hits']:
            identity1=D1['BlastOutput2'][0]['report']['results']['search']['hits'][i]['hsps'][0]['identity']
            alig_len1=D1['BlastOutput2'][0]['report']['results']['search']['hits'][i]['hsps'][0]['align_len']
            query_len1=D1['BlastOutput2'][0]['report']['results']['search']['query_len']
            ref_id1=D1['BlastOutput2'][0]['report']['results']['search']['hits'][i]['description'][0]['id']
            Perc_idn1=identity1/alig_len1*100
            Coverage1=alig_len1/query_len1*100
            
            if D2['BlastOutput2'][0]['report']['results']['search']['hits']:
                identity2=D2['BlastOutput2'][0]['report']['results']['search']['hits'][i]['hsps'][0]['identity']
                alig_len2=D2['BlastOutput2'][0]['report']['results']['search']['hits'][i]['hsps'][0]['align_len']
                query_len2=D2['BlastOutput2'][0]['report']['results']['search']['query_len']
                ref_id2=D2['BlastOutput2'][0]['report']['results']['search']['hits'][i]['description'][0]['id']
                Perc_idn2=identity2/alig_len2*100
                Coverage2=alig_len2/query_len2*100
            
                if Perc_idn1==100 and Coverage1==100 and Perc_idn2==100 and Coverage2==100:
                    
                    print('')
                    print('===========================')
                    print(db)
    
    
                    print(ref_id1)
                    print('Query length1 : ' +  str(query_len1))
                    print('alig_len length1 : ' +  str(alig_len1))
                    
                    print('coverage1: ' + str(Coverage1))
                    print('percent identity1: ' + str(Perc_idn1))
                    print(D1['BlastOutput2'][0]['report']['results']['search']['hits'][i]['hsps'][0]['evalue'])
                    print(D1['BlastOutput2'][0]['report']['results']['search']['hits'][i]['hsps'][0])
                    print('2')
                    print(ref_id2)
                    print('Query length2 : ' +  str(query_len2))
                    print('alig_len length2 : ' +  str(alig_len2))
                    
                    print('coverage2: ' + str(Coverage2))
                    print('percent identity2: ' + str(Perc_idn2))
                    print(D2['BlastOutput2'][0]['report']['results']['search']['hits'][i]['hsps'][0]['evalue'])
                    print(D2['BlastOutput2'][0]['report']['results']['search']['hits'][i]['hsps'][0])
