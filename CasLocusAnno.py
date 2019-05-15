#!/usr/bin/env python
# -*- coding: utf-8 -*-

# notice
'''
CasLocusAnno: A tool for annotate Cas proteins, Cas loci, and their corresponding (sub)type
We developed this tool refered to Makarova et al's method
we also introduced MCL algorithm to handle the overlapping segments of PSI-BLAST results
Therefore, if you used this script, please cite the following work:
1. An updated evolutionary classification of CRISPR–Cas systems
2. Stijn van Dongen, Graph Clustering by Flow Simulation. PhD thesis, University of Utrecht
3. CasLocusAnno: a web-based server for annotating Cas loci and their corresponding (sub)types

Version 1.0
Author: Chuand Dong (chuand@cefg.cn)
'''

import sys, os, Bio, getopt
import logging
from Bio import SeqIO
from multiprocessing import Pool
from multiprocessing.dummy import Pool as ThreadPool
from operator import itemgetter
from itertools import groupby

help_info = \
'''
Usage:
CasLocusAnno.py -i <inputFile> -o <OutFile>

Options:
-h, --help          Check help of CasLocusAnno.py
-v, --version       Check the currtent release version
-i, --input         Input file name 
-o, --out           Out file name
-i, -o are required parameters.
-t is optional parameter.

The code can be copiable, distributable, modifiable, and usable without any restrictions.
Report bugs to <chuand@cefg.cn> or <fbguo@uestc.edu.cn>
'''

script_path = os.path.dirname(os.path.realpath(__file__))
mcl_command = script_path + "/bin/mcl"
makedb_command = script_path + "/bin/makeblastdb"
psiblast_command = script_path + "/bin/psiblast"
zcurve_command = script_path+"/bin/zcurve3.0"
profile_dir = script_path + "/profiles/"
temp_dir = script_path+"/temp/"
profile_ini = script_path+"/configure/profile.ini"

class additionalException(Exception):
	def __init__(self, errorInfo):
		self.errorInfo = errorInfo

def options(argv):
	#Parameters checking step
	try:
		longOpt_list = ['help', 'version', 'input=', 'out=']
		opts, args = getopt.getopt(argv, '-h-v-i:-o:', longOpt_list)
	except getopt.error, msg:
		print("Wrong parameters")
		print(help_info)
		exit()
	for option, value in opts:
		if option in ("-h", "--help"):
			print(help_info)
			exit()
		if option in ("-v", "--version"):
			print("CasLocusAnno.py release v_1.0")
			exit()
		if option in ("-i", "--input"):
			inputFile = value
		if option in ("-o", "--out"):
			outFile = value
	return inputFile, outFile

def makeTemp():
	#This method used to built temporary file to store data
	if os.path.exists(temp_dir):
		os.system("rm -rf {0}*".format(temp_dir))
	else:
		os.mkdir(temp_dir)

def filtration_and_completion(potential_cas_locus):
	if potential_cas_locus:
		continuous=[]
		potential_cas_locus=list(set(potential_cas_locus))
		potential_cas_locus.sort()
		for k, g in groupby(enumerate(potential_cas_locus),lambda (i,x):i-x):
			continuous.append(list(map(itemgetter(1), g)))
		continuous_num = []
		for i in continuous:
			continuous_num.append(len(i))
		max_continuous_num = max(continuous_num)
		num_mccs = continuous_num.count(max_continuous_num)
		if max_continuous_num>=3:
			loop_start = continuous_num.index(max_continuous_num)
			new_potential_cas_locus = continuous[loop_start]
			if loop_start == 0:
				loop_index = range(1, len(continuous))
				for i1 in loop_index:
					gap_num = continuous[i1][0]-new_potential_cas_locus[-1]
					if gap_num<=2 and len(continuous[i1])>=2:
						new_potential_cas_locus.append(new_potential_cas_locus[-1]+1)
						new_potential_cas_locus.extend(continuous[i1])
				continuous = new_potential_cas_locus
				return continuous
			elif loop_start == len(continuous)-1:
				loop_index = range(0, len(continuous)-1)
				loop_index.reverse()
				for i1 in loop_index:
					gap_num = new_potential_cas_locus[0]-continuous[i1][-1]
					if gap_num<=2 and len(continuous[i1])>=2:
						new_potential_cas_locus.extend(continuous[i1])
						new_potential_cas_locus.append(continuous[i1][-1]+1)
						new_potential_cas_locus.sort()
				continuous = new_potential_cas_locus
				return continuous
			else:
				up_loop_index = range(0, loop_start)
				up_loop_index.reverse()
				down_loop_index = range(loop_start+1, len(continuous))
				for i1 in up_loop_index:
					gap_num = new_potential_cas_locus[0]-continuous[i1][-1]
					if gap_num<=2 and len(continuous[i1])>=2:
						new_potential_cas_locus.extend(continuous[i1])
						new_potential_cas_locus.append(continuous[i1][-1]+1)
						new_potential_cas_locus.sort()
				for i1 in down_loop_index:
					gap_num = continuous[i1][0]-new_potential_cas_locus[-1]
					if gap_num<=2 and len(continuous[i1])>=2:
						new_potential_cas_locus.append(new_potential_cas_locus[-1]+1)
						new_potential_cas_locus.extend(continuous[i1])
						new_potential_cas_locus.sort()
				continuous = new_potential_cas_locus
				return continuous
		else:
			continuous=[]
			return continuous

def seq_index(file):
	fasta_f = SeqIO.parse(file, "fasta")
	id2index={}
	index2id={}
	index2seq={}
	temp_index=1
	for seq_record in fasta_f:
		id2index[str(seq_record.id)]=temp_index
		index2id[temp_index]= str(seq_record.id)
		index2seq[temp_index] = str(seq_record.seq)
		temp_index=temp_index+1
		
	return id2index,index2id,index2seq

def make_db(inputFile, inputFileName):
	dbmake = "{0} -in {1} -dbtype prot -out {2} ".format(makedb_command, inputFile, temp_dir+inputFileName)
	os.popen(dbmake)

def search_list(e_value, inBaseName):
	psiblast_list = []
	profile_list = os.listdir(profile_dir)
	for i in profile_list:
		profile = profile_dir+i
		db = temp_dir + inBaseName
		out = temp_dir+i+"_"+inBaseName+".blast"
		psi_command = "{0} -in_msa {1} -db {2} -evalue {3} -outfmt\
		'6 sseqid sstart send evalue' -out {4}  -comp_based_stats 0 > {5}log 2>&1".\
		 format(psiblast_command, profile, db, e_value, out, temp_dir)
		psiblast_list.append(psi_command)
	return psiblast_list

def highest_scoring_profile(sr_path,InfileName,result_file):
	out = open(result_file,"w")
	profile_list = os.listdir(sr_path)
	for i in profile_list:
		file_name = temp_dir+i+"_"+InfileName+".blast"
		if os.path.exists(file_name):
			highest_scoring_subject = open(file_name).read().split("\n")[0]
			if highest_scoring_subject!='':
				out.write(highest_scoring_subject+"\t"+i+"\n")
	out.close()

def deal_file(file, out_anno, profile2cas):
	annotation =open(out_anno,"w")
	f = open(file).read().split("\n")[0:-1]
	cas_dict = {}
	sr2seq={}
	sr2evalue={}
	for i in f:
		line_info = i.split("\t")
		seq_id = line_info[0]
		cas_dict[seq_id] = []
		sr2seq[line_info[-1]] = seq_id
		sr2evalue[line_info[-1]] = float(line_info[3])
	for i in f:
		line_info = i.split("\t")
		seq_id = line_info[0]
		seq_id_info = line_info[1:5]
		cas_dict[seq_id].append(seq_id_info)
	for i in cas_dict.keys():
		single_sr = []
		
		if len(cas_dict[i])==1:
			annotation.write(sr2seq[cas_dict[i][0][-1]]+"\t"+str(sr2evalue[cas_dict[i][0][-1]])+"\t"+cas_dict[i][0][-1]+"\t"+profile2cas[cas_dict[i][0][-1]]+"\n")
		else:
			mcl=open(file+"sr.mcl","w")
			for j in cas_dict[i]:
				
				j_range = range(int(j[0]),int(j[1]))
				for k in cas_dict[i]:
					if j !=k:
						k_range = range(int(k[0]),int(k[1]))
						j_k_or = float(len(list(set(j_range).intersection(set(k_range)))))/min(len(j_range),len(k_range))
						if j_k_or>=0.9:
							mcl.write(j[-1]+" "+k[-1]+" "+"1"+"\n")
							
						else:
							mcl.write(j[-1]+" "+k[-1]+" "+"0"+"\n")
							
			mcl.close()
			os.system("{0} {1}sr.mcl --abc -o {1}.out.sr.mcl > {2}log 2>&1".format(mcl_command, file, temp_dir))
			f_mcl = open(file+".out.sr.mcl").read().split("\n")[0:-1]
			for i_1 in f_mcl:
				info_i_1 = i_1.split("\t")
				temp_e = float(sr2evalue[info_i_1[0]])
				best_sr = info_i_1[0]
				for i_2 in info_i_1:
					
					if float(sr2evalue[i_2]) <temp_e:
						best_sr = i_2
						temp_e=float(sr2evalue[i_2])
				annotation.write(sr2seq[best_sr]+"\t"+str(sr2evalue[best_sr])+"\t"+best_sr+"\t"+profile2cas[best_sr]+"\n")
	annotation.close()

def deal_anno(file,id2index,index2id,index2seq,vicility):
	anno_file = open(file).read().split("\n")[0:-1]
	first_round_hit =[]
	seq_num=len(id2index.keys())
	all_index_list=[]
	first_round_index=[]
	neighbor_sequence_file=open(vicility,"w")
	for i in anno_file:
		sequence_id=i.split("\t")[0]
		first_round_hit.append(sequence_id)
		first_round_index.append(id2index[sequence_id])
	for j in first_round_hit:
		if id2index[j]>20:
			up_index= range(id2index[j]-20,id2index[j])
			if id2index[j]+20<=seq_num:
				down_index=range(id2index[j]+1,id2index[j]+21)
			else:
				down_index=range(id2index[j]+1,seq_num)
		else:
			up_index= range(1,id2index[j])
			max_key=max(id2index.values())
			if max_key >= id2index[j]+21:
				down_index=range(id2index[j]+1,id2index[j]+21)
			else:
				down_index=range(id2index[j]+1,max_key)
		all_index_list.extend(up_index)
		all_index_list.extend(down_index)
		all_index_list=list(set(all_index_list))
		all_index_list.sort()
	for k in all_index_list:
		if k not in first_round_index:
			neighbor_sequence_file.write(">"+index2id[k]+"\n")
			neighbor_sequence_file.write(index2seq[k]+"\n")
	neighbor_sequence_file.close()

def profileMapName():
	profile_contents = open(profile_ini).read().split("\n")
	profile2cas_gene={}
	for i in profile_contents:
		i_info = i.split("\t")
		profile2cas_gene[i_info[0]+".sr"]=i_info[1]
	return profile2cas_gene

def merge_trim(anno1,anno2,id2index,index2id,index2seq,out_anno):
	cas_locus_info=open(out_anno,"w")
	profile2cas_gene={}
	profile2cas_type={}
	profile_contents = open(profile_ini).read().split("\n")
	effector_module_cas_list=[]
	all_cas_type=[]
	cas_core_list=[]
	for i in profile_contents:
		i_info = i.split("\t")
		profile2cas_gene[i_info[0]+".sr"]=i_info[1]
		profile2cas_type[i_info[0]+".sr"]=i_info[3]
		cas_gene_type=i_info[3]
		all_cas_type.extend(cas_gene_type.split(","))
		if i_info[4] == "effector":
			effector_module_cas_list.append(i_info[1])
		if i_info[2]=="core":
			cas_core_list.append(i_info[1])
	all_cas_type=list(set(all_cas_type))
	seq_num=len(id2index.keys())
	f1=open(anno1).read().split("\n")[0:-1]
	cas_index=[]
	real_cas={}
	merge_index=[]
	psi_2=[]
	all_blast_index=[]
	index2record={}
	for i1 in f1:
		seq_id=i1.split("\t")[0]
		cas_index.append(id2index[seq_id])
		index2record[id2index[seq_id]]=[]
	for i1 in f1:
		seq_id=i1.split("\t")[0]
		index2record[id2index[seq_id]].append(i1.split("\t"))
	cas_index=list(set(cas_index))
	cas_index.sort()
	f2=open(anno2).read().split("\n")[0:-1]
	locus_dict={}
	for i2 in f2:
		seq_id=i2.split("\t")[0]
		psi_2.append(id2index[seq_id])
		index2record[id2index[seq_id]]=[]
	for i2 in f2:
		seq_id=i2.split("\t")[0]
		index2record[id2index[seq_id]].append(i2.split("\t"))
	for i in cas_index:
		if i>=5:
			up_locus=range(i-5,i)
			if i+5<=seq_num:
				down_locus=range(i+1,i+6)
			else:
				down_locus=range(i+1,seq_num)
		else:
			up_locus= range(1,i)
			down_locus=range(i+1,i+6)
		for locus_up in up_locus:
			pass
		for locus_down in down_locus:
			pass
		real_cas[i]=[up_locus, down_locus]
	j = 0
	str1 = ''
	for i, item in enumerate(cas_index):
		if i > 0:
			if cas_index[i]  >= cas_index[i-1] + 10:
				tmp = cas_index[j:i]
				if len(tmp) == 1:
					str1 += str(tmp[0]) + ','
					merge_index.append(tmp[0])
				else:
					str1 += str(tmp[0]) + "~" + str(tmp[-1]) + ','
					merge_index.append([tmp[0],tmp[-1]])
				j = i
	tmp2 = cas_index[j:]
	if len(tmp2) == 1:
		str1 += str(tmp2[0]) + ','
		merge_index.append(tmp2[0])
	else:
		str1 += str(tmp2[0]) + "~" + str(tmp2[-1]) + ','
		merge_index.append([tmp2[0],tmp2[-1]])
	res = str1[:-1].split(",")
	all_blast_index.extend(cas_index)
	all_blast_index.extend(psi_2)
	for i in merge_index:
		if isinstance(i,int):
			potential_index=range(i-5,i+6)
			real_key=i
			locus_dict[real_key]=[]
		else:
			potential_index=range(i[0]-5,i[1]+6)
			real_key=i[0]
			locus_dict[real_key]=[]
		for j in potential_index:
			if j in all_blast_index:
				locus_dict[real_key].append(j)
		if len(locus_dict[real_key])>2:
			psi_blast_cas_list=[]
			psi_blast_cas_type_list=[]
			continuous = filtration_and_completion(locus_dict[real_key])
			locus_dict[real_key] = continuous
			for k in locus_dict[real_key]:
				if k in index2record.keys():
					if len(index2record[k])==1:
						profile_name=index2record[k][0][-2]
						cas_locus_info.write(str(k)+"\t"+"\t".join(index2record[k][0])+"\t"+profile2cas_gene[profile_name]+"\t"+profile2cas_type[profile_name]+"\n")
						psi_blast_cas_type_list.extend(profile2cas_type[profile_name].split(","))
						psi_blast_cas_list.append(profile2cas_gene[profile_name])
					else:
						for m1 in index2record[k]:
							profile_name=m1[-2]
							cas_locus_info.write(str(k)+"\t"+"\t".join(m1)+"\t"+profile2cas_gene[profile_name]+"\t"+profile2cas_type[profile_name]+"\n")
							psi_blast_cas_type_list.extend(profile2cas_type[profile_name].split(","))
							psi_blast_cas_list.append(profile2cas_gene[profile_name])
				else:
					cas_locus_info.write(str(k)+"\t"+index2id[k]+"\t"+"-"+"\t"+"-"+"\t"+"-"+"\t-"+"\n")
			temp_cutoff=0
			if set(psi_blast_cas_list) & set(cas_core_list):
				temp_cutoff=1
			if temp_cutoff==0:
				cas_locus_info.write ("Cas locus type: Not a valid cas locus"+"\n")
			else:
				if set(effector_module_cas_list) & set(psi_blast_cas_list):
					type2cas={}
					type2cas['I-D']=['cas10d', 'csc1gr5', 'csc2gr7']
					type2cas['I-A']=['cas5a', 'cas8a1', 'cas8a2', 'cas8a3', 'cas8a4', 'cas8a5', 'cas8a6', 'cas8a7', 'cas8a8']
					type2cas['I-F']=['cas5f', 'cas7f', 'cas8f']
					type2cas['I-U']=['cas5u', 'cas8u1', 'cas8u2', 'csb1gr7', 'csb2gr5']
					type2cas['I-B']=['cas8b1', 'cas8b10', 'cas8b2', 'cas8b3', 'cas8b4', 'cas8b5', 'cas8b6', 'cas8b7', 'cas8b8', 'cas8b9']
					type2cas['I-C']=['cas8c']
					type2cas['I-E']=['cas8e']
					type2cas['III-B']=['cmr8gr7']
					type2cas['V-A']=['cpf1']
					type2cas['IV-A']=['csf1gr8', 'csf2gr7', 'csf3gr5']
					type2cas['III-A']=['csm5gr7']
					type2cas['III-D']=['csx10gr5']
					type2cas['V-B']=['c2c1']
					type2cas['VI-A']=['c2c2']
					type2cas['VI-C']=['c2c3']
					type2cas['V-D']=['Cas12d']
					type2cas['V-E']=['Cas12e']
					type2cas['V-F']=['Cas14e','Cas14f','Cas14a','Cas14b','Cas14d','Cas14c','Cas14u','Cas14h']
					type2cas['VI-B']=['cas13b']
					type2cas['VI-C']=['cas13c']
					type2cas['VI-D']=['Cas13d']
					one_cas2one_type=[]
					for keys_in_type2cas in type2cas.keys():
						one_cas2one_type.extend(type2cas[keys_in_type2cas])
					if set(psi_blast_cas_list) & set(one_cas2one_type):
						subtype=[]
						for keys_in_type2cas in type2cas.keys():
							if set(psi_blast_cas_list) & set(type2cas[keys_in_type2cas]):
								subtype.append(keys_in_type2cas)
						if len(subtype)==1:
							cas_locus_info.write ("Cas locus type: CAS-"+subtype[0]+"\n")
						else:
							cas_type2vote_num={}
							for m_3 in all_cas_type:
								cas_type2vote_num[m_3]=psi_blast_cas_type_list.count(m_3)
							max_num=0
							temp_cas=""
							votes_number_list=[]
							new_cas2type={}
							for m_4 in cas_type2vote_num.keys():
								if cas_type2vote_num[m_4]!=0:
									if m_4.count("-") ==2:
										new_cas2type[m_4]= psi_blast_cas_type_list.count("-".join(m_4.split("-")[0:2]))+psi_blast_cas_type_list.count(m_4)
									else:
										new_cas2type[m_4]=psi_blast_cas_type_list.count(m_4)
							final_cas=[]
							max_value=max(new_cas2type.values())
							for each_key in new_cas2type.keys():
								if new_cas2type[each_key]==max_value:
									final_cas.append(each_key)
							if len(final_cas)==1:
								cas_locus_info.write ("Cas locus type: "+",".join(final_cas)+"\n")
							else:
								cas_locus_info.write ("Cas locus type: Ambiguously classified "+",".join(final_cas)+"\n")
					else:
						cas_type2vote_num={}
						for m_3 in all_cas_type:
							cas_type2vote_num[m_3]=psi_blast_cas_type_list.count(m_3)
						max_num=0
						temp_cas=""
						votes_number_list=[]
						new_cas2type={}
						for m_4 in cas_type2vote_num.keys():
							if cas_type2vote_num[m_4]!=0:
								if m_4.count("-") ==2:
									new_cas2type[m_4]= psi_blast_cas_type_list.count("-".join(m_4.split("-")[0:2]))+psi_blast_cas_type_list.count(m_4)
								else:
									new_cas2type[m_4]=psi_blast_cas_type_list.count(m_4)
						final_cas=[]
						max_value=max(new_cas2type.values())
						for each_key in new_cas2type.keys():
							if new_cas2type[each_key]==max_value:
								final_cas.append(each_key)
						if len(final_cas)==1:
							cas_locus_info.write ("Cas locus type: "+",".join(final_cas)+"\n")
						else:
							cas_locus_info.write ("Cas locus type: Ambiguously classified "+",".join(final_cas)+"\n")
				else:
					cas_type2vote_num={}
					for m_3 in all_cas_type:
						cas_type2vote_num[m_3]=psi_blast_cas_type_list.count(m_3)
					votes_number_list=[]
					new_cas2type={}
					for m_4 in cas_type2vote_num.keys():
						if cas_type2vote_num[m_4]!=0:
							if m_4.count("-") ==2:
								new_cas2type[m_4]= psi_blast_cas_type_list.count("-".join(m_4.split("-")[0:2]))+psi_blast_cas_type_list.count(m_4)
							else:
								new_cas2type[m_4]=psi_blast_cas_type_list.count(m_4)
					final_cas=[]
					max_value=max(new_cas2type.values())
					for each_key in new_cas2type.keys():
						if new_cas2type[each_key]==max_value:
							final_cas.append(each_key)
					cas_locus_info.write ("Cas locus type: Partial Cas locus "+",".join(final_cas)+"\n")
			cas_locus_info.write ("====="+"\n\n")
		else:
			cas_locus_info.write ("Cas locus type: No Cas locus"+"\n")
			cas_locus_info.write ("====="+"\n\n")
	cas_locus_info.close()

def formatOut(fileIn):
	pass
	'''
	明天搞这个
	'''

def geneAnno(inFile):
	inFileProtein = os.path.basename(inFile)+".protein"
	gene_anno_cmd = "{0} {1} {2} -p".format(zcurve_command, inFile, inFileProtein)
	print(gene_anno_cmd)
	os.system(gene_anno_cmd)
	return inFileProtein

def main(argv=sys.argv[1:]):
	inputFile, outFile = options(argv)
	inputFilePath = os.path.dirname(inputFile)
	inputFileName = os.path.basename(inputFile)
	makeTemp()
	make_db(inputFile, inputFileName)
	pool = ThreadPool(18)
	psi_command_list = search_list("1e-06", inputFileName)
	pool.map(os.popen, psi_command_list)
	blast1 = outFile + ".hit1"
	highest_scoring_profile(profile_dir, inputFileName, blast1)
	anno1 = outFile+".anno1"
	profile2cas = profileMapName()
	deal_file(blast1,anno1, profile2cas)
	id2index,index2id,index2seq = seq_index(inputFile)
	deal_anno(anno1,id2index,index2id,index2seq, "{0}{1}.neighbor".format(temp_dir,inputFileName))
	make_db("{0}{1}.neighbor".format(temp_dir,inputFileName), inputFileName+".neighbor")
	psi_command_list = search_list("0.01", "{0}.neighbor".format(inputFileName))
	pool.map(os.popen, psi_command_list)
	blast2 = outFile + ".hit2"
	highest_scoring_profile(profile_dir, "{0}.neighbor".format(inputFileName), blast2)
	pool.close()
	anno2 = outFile+".anno2"
	deal_file(blast2, anno2, profile2cas)
	merge_trim(anno1, anno2, id2index, index2id, index2seq, outFile)
	if blast1:
		os.remove(blast1)
	if blast2:
		os.remove(blast2)
	if outFile+".hit2sr.mcl":
		os.remove(outFile+".hit2sr.mcl")
	if outFile+".hit1sr.mcl":
		os.remove(outFile+".hit1sr.mcl")
	if outFile+".hit2.out.sr.mcl":
		os.remove(outFile+".hit2.out.sr.mcl")
	if outFile+".hit1.out.sr.mcl":
		os.remove(outFile+".hit1.out.sr.mcl")
	os.system("rm -rf {0}*".format(temp_dir))

if __name__ == "__main__":
	main()
