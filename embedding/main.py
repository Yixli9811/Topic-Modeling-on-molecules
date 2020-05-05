import matplotlib.pyplot as plt
from matplotlib import pyplot as mp
import numpy as np 
import sys 
import os
import csv
import random
import math
from rdkit import Chem
from rdkit.Chem import PandasTools
from rdkit.Chem.Draw import IPythonConsole
from mol2vec.features import mol2alt_sentence, mol2sentence, MolSentence, DfVec, sentences2vec
from gensim.models import word2vec


#filter failure mols 
def readmols(suppl):
	ok=[]
	failures=[]
	sio = sys.stderr = StringIO()
	for i,m in enumerate(suppl):
		if m is None:
			failures.append(i)
			sio = sys.stderr = StringIO() # reset the error logger
		else:
			ok.append(m)
	return ok, failures


def clean_file(file_in, cleaned_fn):

	mpt_file = file_in
	df = pd.read_csv(mpt_file)

	df.sort_values("smiles", inplace=True)
	df.drop_duplicates(subset = 'smiles', keep=False, inplace=True, ignore_index=True)

	smiles_str = df['smiles'].values
	
	#mpc = data['mpC'].values

	#df = pd.DataFrame(data=smiles_str, columns=["smiles"])
	#df['mpC'] = mpc

	mols = [] 
	fail_idx = []
	errors = [] 
	for i, smi in enumerate(smiles_str):
		try:
			mol = Chem.MolFromSmiles(smi)
			mols.append(mol)

		except TypeError:
			errors.append(i)
	


	read_mols, fails = readmols(mols)
	fails.extend(errors)

	print("# of total smiles: ", len(smiles_str))
	print("# of fails smiles: ", len(fails))
	print("# of valid smiles: ", len(read_mols))
	
	#backup dropped smi string index 
	with open('./dropped.txt', 'w') as f:
		f.write(str(fails))

	#rewrite clean file 
	new_df = df.drop(fails)

	new_df.reset_index(drop=True, inplace=True)

	new_df.to_csv(cleaned_fn)
	print("write to file done")



def mol2vec(fin_name, fout_name, clean=False):
	
	#clean_data, removing smiles string can't convert to molecules 
	#We may improve this latter. Only do once 

	if clean:
		print('cleaning data...')
		clean_file(fin_name, fin_name)

	clean_data = pd.read_csv(fin_name)

	#Load pre-trained model 
	model = word2vec.Word2Vec.load('./models/model_300dim.pkl')

	print('making vec data...')
	#convert to sentences 
	mols = [Chem.MolFromSmiles(smi) for smi in clean_data['smiles'].values]
	sentences = [MolSentence(mol2alt_sentence(mol, 1)) for mol in mols]

	#convert to vectors 
	vecs = [DfVec(x) for x in sentences2vec(sentences, model, unseen='UNK')]
	vec_values = np.array([v.vec for v in vecs])

	# Form dataframe 
	cols = ['vec_'+str(i) for i in range(300)]
	df = pd.DataFrame(vec_values, columns=cols)
	df.insert(0, "smiles", clean_data['smiles'].values, True) 

	df.to_csv(fout_name)


	return vec_values


#This bolck code is going to conver molecule names to SMILES strings 
#We may find better way to do this in the future. 

def names_to_smiles():
	bpt_fn = './bpt_data.xls'
	data = pd.read_excel(bpt_fn)

	names = data['Substance Name'].values[:10]
	
	for ids in names :
		print(ids, CIRconvert(ids))

#names_to_smiles()

def main():

	### Make mol2vec file ####
	'''
	#mpt_file = './data/BradleyMeltingPointDataset.xlsx'
	mpt_cleaned = './data/cleaned_MeltingPt.csv'
	vec_file = './data/cleaned_MeltingPt_vecs.csv'
  
  mol2vec(mpt_cleaned, fout_name=vec_file, clean=True)

main()
  
  
