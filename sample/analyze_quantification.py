import pandas as pd
import numpy as np
import sys, getopt
from scipy.stats import spearmanr

true = None
estimated = None
kupid = None

args = sys.argv[1:]
options = "tek"

try:
  opts, args = getopt.getopt(args, "t:e:k:")
except getopt.error as err:
  print(str(err))
  
for opt, arg in opts:
  if opt == '-t':
    true = arg
  elif opt == '-e':
    estimated = arg
  elif opt == 'k':
    kupid = arg

true_tpm = pd.read_csv(true)
estim_tpm = pd.read_csv(estimated)
KuPID_tpm = pd.read_csv(KuPID)

both = pd.merge(true_tpm, estim_tpm, on='Transcript', how='outer')
both_KuPID = pd.merge(true_tpm, KuPID_tpm, on='Transcript', how='outer')

both.fillna(0, inplace=True)
both_KuPID.fillna(0, inplace=True)
                                                                                       
corr, p_value = spearmanr(both['TPM'], both['True TPM'])                                
corr_KuPID, p_value = spearmanr(both_KuPID['True TPM'], both_KuPID['Scaled TPM'])

print("Readset\tSpearman Coeffiecient")
print("Non-Processed Reads\t" + str(corr))
print("KuPID Reads\t" + str(corr_KuPID))


