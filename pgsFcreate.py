import sys
import numpy as np
import pandas as pd
from bed_reader import open_bed,sample_file
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import Lasso
from sklearn.linear_model import Ridge
from sklearn.metrics import r2_score
from scipy.stats import pearsonr
# inPath=sys.argv[1] # This is the genome file path
# annotPath=sys.argv[2] # This is the annotation file path
savePath=sys.argv[1]
Type=sys.argv[2]
# savePath="/u/scratch/b/boyang19/tmp/u/flashscratch/b/boyang19/GxG/code/robustness/LD_pgs/results/annot/"
# Type="Ridge"

# famsort=pd.read_csv('/u/scratch/b/boyang19/Angela/data/geno/filter4_wes.fam',delim_whitespace=True,header=None)
# annot_df = pd.read_csv(annotPath,sep=' ',header=None)
# famdf = pd.read_csv(f'/u/project/sgss/UKBB/data/wes/plink-23155/ukb23155_c4_b0_v1.fam',delim_whitespace=True,header=None)
# # file_name = sample_file(f'{inPath}.bed')
# columns=['FID','IID','feature','y','z','q']
# famsort.columns=['FID','IID','x','y','z','q']
# famdf.columns=columns
# assert famdf.shape[0] == len(np.unique(famdf.FID))
def converter(X,scheme='unweighted',traitPath=None):
    miss_mean_imputer = SimpleImputer(strategy='mean')
    miss_mean_imputer = miss_mean_imputer.fit(X)
    data = data=miss_mean_imputer.transform(X)
    scaler = StandardScaler()
    data = scaler.fit_transform(data)
    imputed_df = pd.DataFrame(data=data)
    top_values = []
    print(f'LD region shape is {X.shape}')
    if scheme == "unweighted":
        aifeature = np.mean(imputed_df.values,axis=1)
        # print(f'artificial feature: {len(aifeature)}')
        scaler = StandardScaler()
        aifeature=scaler.fit_transform(aifeature.reshape(-1,1))
        return aifeature

    if scheme == 'Ridge':
        clf = Ridge()
        pheno = pd.read_csv(traitPath,delim_whitespace=True)
        miss_mean_imputer = SimpleImputer(strategy='mean')
        miss_mean_imputer = miss_mean_imputer.fit(pheno.pheno.values.reshape(-1,1))
        pheno_val = miss_mean_imputer.transform(pheno.pheno.values.reshape(-1,1))
        clf.fit(imputed_df,pheno_val)
        aifeature = clf.predict(imputed_df)
        # aifeature = np.mean(imputed_df.values,axis=1)
        # print(f'artificial feature: {len(aifeature)}')
        scaler = StandardScaler()
        aifeature=scaler.fit_transform(aifeature.reshape(-1,1))
        print(f'correlation is {pearsonr(pheno_val.flatten(),aifeature.flatten())}')
        return aifeature
        

        # return aifeature
    else:
        raise NotImplementedError(f'scheme: {scheme} is not implemented')
    
    
#genPath='/u/project/sgss/UKBB/data/cal/filter4'
#bimfile = pd.read_csv(f'{genPath}.bim',header=None,delimiter='\t')
genPath='/u/home/a/aanand2/rGxG/GENIE_gxg/example/small'
bimfile = pd.read_csv(f'{genPath}.bim',header=None,delimiter='\t')

covPath="/u/scratch/a/aanand2/no_epistasis_pheno/sim_geno/cau_ratio_0.01_h2_0.25_pheno/"

# bimfile=pd.read_csv('/u/home/b/boyang19/project-sriram/GxG/imputed/all_filter.bim',sep='\t',header=None)
# bimfile.columns = ['chr', 'chrpos', 'MAF', 'pos','MAJ','MIN']

columns = ['chr', 'ID', 'unkown', 'loci', 'Ma', 'Mi']
bimfile.columns = columns
bimfile['chr'] = bimfile['chr'].astype(str)
bimfile['loci'] = bimfile['loci'].astype(int)
N = bimfile.shape[0]
    
    
bed = open_bed(f'{genPath}.bed')


sigpath="/u/scratch/a/aanand2/PGS_test/sig_cand_mkup/"
cand_df = pd.read_csv(f'{sigpath}/1_sig_test.txt',header=None,delimiter=' ')

annot_local=np.zeros(N)
annot_dist = np.zeros(N)
annot_LD = np.zeros(N)
data = []
# bimfile.pos = bimfile.pos.astype(int)
for Index,row in cand_df.iterrows():
    print(row)
    trait,rIndex,CHR1 = row # rIndex is the relative index in the bim file
    #covdf = pd.read_csv(covPath+trait+".covar",delim_whitespace=True)
    covdf = pd.read_csv("/u/home/a/aanand2/rGxG/GENIE_gxg/example/small.cov",delim_whitespace=True)

    # for rIndex in rIndices:
        # print(rIndex)
    inforow = bimfile.iloc[rIndex]
    CHR, snpindex  = inforow['chr'],inforow['loci']
    print(inforow)
        # print(inforow)
        # if CHR==CHR1:
        #     break
    assert str(CHR) == str(CHR1)
    # print(trait, pIndex, rIndex, CHR1)
    LDpath='/u/home/a/aanand2/LDblock/EUR/'
    LDname=f'{LDpath}fourier_ls-chr{CHR1}.bed'
    LD_pd=pd.read_csv(LDname,delim_whitespace=True)
    start_index = np.where(LD_pd.start<=snpindex)[0][-1]
    stop_index = np.where(LD_pd.stop>=snpindex)[0][0]
    start_pindex = LD_pd.start[start_index]
    stop_pindex = LD_pd.stop[stop_index]
    loci = bimfile.iloc[rIndex,3]
    Chr2 = bimfile.iloc[rIndex,0]


    val = bed.read(index=np.s_[:,(bimfile.chr==Chr2)&(bimfile.loci>=start_pindex)&(bimfile.loci<=stop_pindex)])
    #aifeature=converter(val,scheme=Type,traitPath=covPath+trait+".pheno")
    aifeature=converter(val,scheme=Type,traitPath=covPath+"1"+".pheno")
    covdf['f_LD'] = aifeature
    covdf.to_csv(f"{savePath}/{trait}-{rIndex}-{CHR}.{Type}.covar",sep=" ",index=False)

    subannot_LD = annot_LD.copy()
    subannot_LD[(bimfile.chr==Chr2)&(bimfile.loci>=start_pindex)&(bimfile.loci<=stop_pindex)] = 1	
    subannot_dist=annot_dist.copy()
    rm1 = bimfile.chr!=Chr2
    rm2 = (bimfile.chr==Chr2)&((bimfile.loci<start_pindex)|(bimfile.loci>stop_pindex))
    subannot_dist[rm1|rm2]=1
    print(f'Index {rIndex}, 1s are {np.sum(subannot_dist)} in different chromosome')
    subannot = np.concatenate((subannot_LD.reshape(-1,1),subannot_dist.reshape(-1,1)),axis=1)
    np.savetxt(f'{savePath}/{trait}-{rIndex}-{CHR}.{Type}.annot',subannot_dist,fmt='%i')
    np.savetxt(f'{savePath}/{trait}-{rIndex}-{CHR}.{Type}.target',aifeature)

LDsummary = pd.DataFrame(data=data,columns=['trait','Index','Chr','Pindex','Start_index','Stop_index'])
LDsummary.to_csv('./LDsummary_mkup.csv')