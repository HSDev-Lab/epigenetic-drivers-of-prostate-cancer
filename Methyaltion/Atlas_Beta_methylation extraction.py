from os import listdir
from os.path import isfile, join


mypath = "/Users/tobymilne-clark/Documents/Data/Atlas/Prostate_Immune"

onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]

os.chdir(mypath)

beta_files = load_methylation_frequencies(onlyfiles,'beta')

betas = beta_files[:, :, 0] / beta_files[:, :, 1]

regions = pd.read_csv('/Users/tobymilne-clark/Documents/Data/Atlas/markers.csv')



DMRS = []
for j in range(len(regions)):
    DMRS.append(np.mean(betas[range(regions.iloc[j]['startCpG'],regions.iloc[j]['endCpG']),],axis=0))

pd.DataFrame(data=DMRS,columns=onlyfiles).to_csv('/Users/tobymilne-clark/Documents/Data/Atlas/DMRS_cellatlas_betavalues.csv')