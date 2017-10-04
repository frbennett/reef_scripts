

# 

def compare(resultsSet1,resultsSet2):
    print('Comparing %s %s'%(resultsSet1, resultsSet2),end='')
    files1 = findFiles(resultsSet1)
    files2 = [f.replace(resultsSet1,resultsSet2) for f in files1]
    pairs = zip(files1,files2)
   
    tmp = {lhs.split('/')[-1]:fileComparison(lhs,rhs) for lhs,rhs in pairs}
    print(' DONE')
    return {key:val for key,val in tmp.items() if not val is None}

    