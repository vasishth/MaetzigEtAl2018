import re
import datetime
import pandas as pd
import numpy as np
from itertools import compress

######################################################################
######################   PREPROCESSING   #############################
######################################################################

## read paramsearch.txt and experiment-data.txt
paramSearch = pd.read_table('./paramsearch.txt', delim_whitespace=True, quotechar='"', header=None)
paramSearch.columns = ['exp', 'simset', 'file', 'params']
expData = pd.read_table('./experiment-data.txt', delim_whitespace=True, quotechar='"')

## get unique conditions from expData
conditions = set(expData['cond'])
conditions = [entry for entry in conditions]
conditions = [re.sub(r"\*", "", entry) for entry in conditions]
conditions = [re.sub(r"C15", "", entry) for entry in conditions]
conditions.sort()  ## this modifies the object itself

## create results data frame; number of lines is (len(simset) * len(conditions))
## because for each simulation set with a given parameter vector, there are 
## len(conditions) observations
resultsDF = pd.DataFrame(np.repeat(np.asarray(paramSearch['simset']), len(conditions)))
resultsDF['condition'] = list(conditions) * len(paramSearch['simset'])
resultsDF.columns = ['simset', 'condition']

## extract parameter names and values from paramSearch
params = [re.sub(r"\(|\)", "", entry) for entry in paramSearch['params']]
params = np.asarray([entry.split() for entry in params])
paramNames = params[:, 0::2][0]
paramValues = np.round(params[:, 1::2].astype(float), decimals=2)
paramValues = pd.DataFrame(np.repeat(paramValues, len(conditions), axis=0))
paramValues.columns = paramNames

## combine the data frames
resultsDF = pd.concat([resultsDF, paramValues], axis=1)
resultsDF['sentCount'] = np.zeros(resultsDF.shape[0])

######################################################################
#######################   PREPROCESSING END   ########################
######################################################################

## initialising data structures for storing correctCounts
sentCounts = np.zeros(int(resultsDF.shape[0]/len(conditions)))
correctCountsVerb = {key: 0 for key in conditions}  ## initialises dict with conditions as keys, counts go in here
correctCountsAntecedent = {key: 0 for key in list(compress(conditions, [len(entry) > 2 for entry in conditions]))}

correctCountsSR = np.zeros(int(resultsDF.shape[0]/len(conditions)))
correctCountsOR = np.zeros(int(resultsDF.shape[0]/len(conditions)))
correctCountsVerbSRREFL = np.zeros(int(resultsDF.shape[0]/len(conditions)))
correctCountsVerbORREFL = np.zeros(int(resultsDF.shape[0]/len(conditions)))
correctCountsVerbSRPRON = np.zeros(int(resultsDF.shape[0]/len(conditions)))
correctCountsVerbORPRON = np.zeros(int(resultsDF.shape[0]/len(conditions)))
correctCountsAntecedentSRREFL = np.zeros(int(resultsDF.shape[0]/len(conditions)))
correctCountsAntecedentORREFL = np.zeros(int(resultsDF.shape[0]/len(conditions)))
correctCountsAntecedentSRPRON = np.zeros(int(resultsDF.shape[0]/len(conditions)))
correctCountsAntecedentORPRON = np.zeros(int(resultsDF.shape[0]/len(conditions)))

## define identifiers and correct check strings
identifierSR = 'who hugged'
identifierREFL = 'herself'
identifierPRON = 'her'
correctVerbSimpleSR = r'"BOY" is "subject" of "HUGGED"'
correctVerbSR = r'"WOMAN" is "subject" of "HUGGED"'
correctVerbOR = r'"GIRL" is "subject" of "HUGGED"'
correctAntecedentREFL = r'"WOMAN" is "antecedent" of "HERSELF"'
correctAntecedentPRON = r'"GIRL" is "antecedent" of "HER"'

## loops through simset, reads the "i-relations.txt" subsequently
for i in paramSearch['simset']:
    fileName = "{}-relations.txt".format(i)
    with open(fileName, "r") as f:
        relations = f.readlines()
    
    for j in range(len(relations)):
        lineNr = j
        if re.search('SENTENCE', relations[lineNr]) is not None:
            sentCounts[i-1] += 1  ## index is decremented because simset starts at 1, not 0
            lineNr += 1
            try:
                if re.search(identifierSR, relations[j]) is not None:
                    if re.search(identifierREFL, relations[j]) is not None:
                        if re.search(correctVerbSR, relations[lineNr]) is not None:
                            correctCountsVerbSRREFL[i-1] += 1
                        if re.search(correctAntecedentREFL, relations[lineNr+2]) is not None:
                            correctCountsAntecedentSRREFL[i-1] += 1
                    elif re.search(identifierPRON, relations[j]) is not None:
                        if re.search(correctVerbSR, relations[lineNr]) is not None:
                            correctCountsVerbSRPRON[i-1] += 1
                        if re.search(correctAntecedentPRON, relations[lineNr+2]) is not None:
                            correctCountsAntecedentSRPRON[i-1] += 1
                    else:
                        if re.search(correctVerbSimpleSR, relations[lineNr]) is not None:
                            correctCountsSR[i-1] += 1
                else:
                    if re.search(identifierREFL, relations[j]) is not None:
                        if re.search(correctVerbOR, relations[lineNr]) is not None:
                            correctCountsVerbORREFL[i-1] += 1
                        if re.search(correctAntecedentREFL, relations[lineNr+2]) is not None:
                            correctCountsAntecedentORREFL[i-1] += 1
                    elif re.search(identifierPRON, relations[j]) is not None:
                        if re.search(correctVerbOR, relations[lineNr]) is not None:
                            correctCountsVerbORPRON[i-1] += 1
                        if re.search(correctAntecedentPRON, relations[lineNr+2]) is not None:
                            correctCountsAntecedentORPRON[i-1] += 1
                    else:
                        if re.search(correctVerbOR, relations[lineNr]) is not None:
                            correctCountsOR[i-1] += 1
            except IndexError:
                print("End of file!")

resultsDF['sentCount'] = np.repeat(sentCounts, len(conditions))
resultsDF['correctCountEmbV'] = np.column_stack((correctCountsOR, correctCountsVerbORPRON, correctCountsVerbORREFL, correctCountsSR, correctCountsVerbSRPRON, correctCountsVerbSRREFL)).flatten() 
resultsDF['correctCountAntecedent'] = np.column_stack(([None] * int(len(correctCountsOR)), correctCountsAntecedentORPRON, correctCountsAntecedentORREFL, [None] * int(len(correctCountsSR)), correctCountsAntecedentSRPRON, correctCountsAntecedentSRREFL)).flatten() 

## divide sentCounts by len(conditions) for accuracy calculation
resultsDF['accuracyEmbV'] = resultsDF['correctCountEmbV'] / (resultsDF['sentCount'] / len(conditions))
resultsDF['accuracyAntecedent'] = resultsDF['correctCountAntecedent'] / (resultsDF['sentCount'] / len(conditions))

print(resultsDF)
resultsFileName = 'accuracies-' + expData['exp'][0] + '-' + datetime.datetime.now().strftime('%Y-%m-%d-%H-%M') + '.csv'
resultsDF.to_csv(resultsFileName)
