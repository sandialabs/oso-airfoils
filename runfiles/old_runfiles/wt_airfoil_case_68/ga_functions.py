import copy
import random
import numpy as np
import multiprocessing
import math
from kulfan import Kulfan

def numberCoding(val, mode, base = 'binary', nbe=8, nbm=23, encodeType=None, encodeLength=None):
    if base == 'binary':
        if mode == 'encode':
            if encodeType is None or encodeType==float:
                if val == 0:
                    ex = ''.rjust(nbe,'0')
                    man = ''.rjust(nbm,'0')
                    return '0' + ex + man
                if val == float('inf'):
                    ex = ''.rjust(nbe,'1')
                    man = ''.rjust(nbm,'0')
                    return '0' + ex + man
                if val == -1*float('inf'):
                    ex = ''.rjust(nbe,'1')
                    man = ''.rjust(nbm,'0')
                    return '1' + ex + man
                if np.isnan(val):
                    ex = ''.rjust(nbe,'1')
                    man = '1'+''.rjust(nbm-1,'0')
                    return '0' + ex + man
                
                val = float(val)
                sgn = abs(val)/val

                res = ''
                if sgn == 1:
                    res += '0'
                else:
                    res += '1'

                if val < 0:
                    wholeNumber = abs(math.ceil(val))
                    fracNumber = abs(val)-wholeNumber
                else:
                    wholeNumber = math.floor(val)
                    fracNumber = val-wholeNumber

                binWholeNumber = bin(wholeNumber)[2:]

                if fracNumber == 0:
                    binFracNumber = ''.rjust(nbm,'0')
                else:
                    binFracNumber = ''
                    itrVal = fracNumber
                    maxIter = 200
                    ctr = 0
                    flg = 0
                    for i in range(0,maxIter):
                        itrVal *= 2
                        if itrVal >= 1:
                            binFracNumber += '1'
                            itrVal -=1
                            flg = 1
                        else:
                            binFracNumber += '0'
                        if flg == 1:
                            ctr += 1
                        if ctr > nbm:
                            break
                    if ctr < nbm:
                        binFracNumber = binFracNumber+ ''.rjust(nbm-ctr+1,'0')

                joinedNumber = binWholeNumber + '.' + binFracNumber
                joinedNoDecimal = binWholeNumber + binFracNumber
                if float(joinedNumber) == 0:
                    ex = ''.rjust(nbe,'0')
                    man = ''.rjust(nbm,'0')
                else:
                    ex_us = math.floor(math.log10(float(joinedNumber)))
                    sft = int(2**nbe/2 - 1)
                    if ex_us <= sft and ex_us >= -sft:
                        ex = bin(ex_us + sft)[2:].rjust(nbe,'0')
                        ix = joinedNoDecimal.find('1')
                        man_long = joinedNoDecimal[ix+1:]
                        man = man_long[0:nbm]
                    elif ex_us > sft:
                        # Infinity
                        ex = ''.rjust(nbe,'1')
                        man = ''.rjust(nbm,'0')
                    else: #ex_us < sft
                        # zero
                        ex = ''.rjust(nbe,'0')
                        man = ''.rjust(nbm,'0')
                res += ex + man
                return res

            elif encodeType == int:
                if val.is_integer():
                    val = int(val)
                try:
                    rawString = "{0:b}".format(val)
                except:
                    if val.is_integer():
                        val = int(val)
                    rawString = "{0:b}".format(val)
                    
                if len(rawString)>encodeLength:
                    raise ValueError('Value %d exceeds encoding length of %d bits'%(val, encodeLength))
                else:
                    newString = rawString.zfill(encodeLength)
                return newString

            else:
                raise ValueError('Invalid Encoding Type')

        elif mode == 'decode':
            if encodeType is None or encodeType==float:
                sgnBit = val[0]
                expBits = val[1:1+nbe]
                manBits = val[1+nbe:]

                zeroExp = ''.rjust(nbe,'0')
                zeroMan = ''.rjust(nbm,'0')
                infExp  = ''.rjust(nbe,'1')
                if expBits == zeroExp:
                    if manBits == zeroMan:
                        return  0.0
                if expBits == infExp:
                    if manBits == zeroMan:
                        if sgnBit == '1':
                            return  -1 * float("inf")
                        elif sgnBit == 0:
                            return float("inf")
                    else:
                        return float('nan')
                    
                man = '1'+manBits
                ex  = 0
                for i in range(0,len(expBits)):
                    ex += 2**i * float(expBits[-(i+1)])
                sft = int(2**nbe/2 - 1)
                ex_us = int(ex - sft)
                
                if ex_us >= 0:
                    man += ''.rjust(ex_us,'0')
                    binWholeNumber = man[0:ex_us+1] 
                    wholeNumber = 0
                    for i in range(0,len(binWholeNumber)):
                        wholeNumber += 2**i * float(binWholeNumber[-(i+1)])
                    binFracNumber = man[ex_us+1:]
                else:
                    wholeNumber = 0.0
                    binFracNumber = ''.rjust(abs(ex_us)-1,'0') + man

                fracNumber = 0.0
                for i in range(0,len(binFracNumber)):
                    fracNumber += float(binFracNumber[i]) * 0.5**(1+i)

                if sgnBit == '1':
                    res =  -1 * (wholeNumber + fracNumber)
                elif sgnBit == '0':
                    res = (wholeNumber + fracNumber)

                return res

            elif encodeType == int:
                res = int(val,2)
                return res

            else:
                raise ValueError('Invalid Encoding Type')
    
        else:
            raise ValueError('Invalid Mode')

    else:
        raise ValueError('Invalid Encoding Base')

def encodeChromosome(x, ets, lowerBounds, upperBounds):
    # vrs = formulation.variables_only
    # ets = [v.type_ for v in vrs]
    chromosome = ''
    for i in range(0,len(x)):
        xv = x[i]
        et = ets[i]
        # v = vrs[i]
        if et == float:
            en = None
        # elif et == bool:
        #     en = 1
        else: # int
            ub = math.ceil(upperBounds[i])
            en = len(bin(ub))-1
        gene = numberCoding(xv, 'encode', 
                            nbe = 6, 
                            nbm = 8, 
                            encodeType=et, encodeLength = en)
        chromosome += gene
    return chromosome

def decodeChromosome(chromosome, ets, lowerBounds, upperBounds):
    # vrs = formulation.variables_only
    # ets = [v.type_ for v in vrs]
    # lnGeneCon = 1 + formulation.solverOptions.nbe + formulation.solverOptions.nbm
    lnGeneCon = 15
    lenGenes = [0]
    x = []
    for i in range(0,len(ets)):
        fstix = sum(lenGenes)
        # v = vrs[i]
        et = ets[i]
        if et == float:
            en = None
        # elif et == bool:
        #     en = 1
        else: # int
            ub = math.ceil(upperBounds[i])
            en = len(bin(ub))-1

        if et == float:
            lenGenes.append(lnGeneCon)
        elif et == int:
            lenGenes.append(en)
        # elif et == bool:
            # lenGenes.append(1)
        else:
            raise ValueError('Invalid encoding type')

        lstix = fstix + lenGenes[-1]
        gene = chromosome[fstix:lstix]
        xv = numberCoding(gene, 'decode', 
                          nbe = 6, 
                          nbm = 8, 
                          encodeType=et, encodeLength=en)
        x.append(xv)

    return np.array(x)
            
def mutateChromosome(chromosome, probOfMutation):
    if probOfMutation >= random.uniform(0,1):
        cl = list(chromosome)
        bitFlip = random.randrange(0,len(cl))
        if cl[bitFlip] == '1':
            cl[bitFlip] = '0'
        else:
            cl[bitFlip] = '1'
        return "".join(cl)
    else:
        return chromosome

def crossoverChromosomes(parent1, parent2, Ncrossovers=1):
    pl1 = list(parent1)
    pl2 = list(parent2)
    for i in range(0,Ncrossovers):
        insertionPoint = random.randrange(0,len(pl1))
        c1 = pl1[0:insertionPoint] + pl2[insertionPoint:]
        c2 = pl2[0:insertionPoint] + pl1[insertionPoint:]
        pl1 = c1
        pl2 = c2
    return "".join(c1), "".join(c2)
    
def breedDesignVectors(x1, x2, normalizationVector, encodingTypes, lowerBounds, upperBounds, Ncrossovers=3, probabilityOfMutation=0.10):
    # normalizationVector = [v.guess.to(v.units).magnitude for v in formulation.variables_only]
    ec1 = np.array(x1)/np.array(normalizationVector)
    ec2 = np.array(x2)/np.array(normalizationVector)

    chm1 = encodeChromosome(ec1, encodingTypes, lowerBounds, upperBounds)
    chm2 = encodeChromosome(ec2, encodingTypes, lowerBounds, upperBounds)

    [child1, child2] = crossoverChromosomes(chm1, chm2, Ncrossovers)
    child1M = mutateChromosome(child1, probabilityOfMutation)
    child2M = mutateChromosome(child2, probabilityOfMutation)

    xres1 = decodeChromosome(child1M, encodingTypes, lowerBounds, upperBounds)
    xres2 = decodeChromosome(child2M, encodingTypes, lowerBounds, upperBounds)

    dc1 = np.array(xres1) * np.array(normalizationVector)
    dc2 = np.array(xres2) * np.array(normalizationVector)
   
    return dc1, dc2

def breedDesignVectorsParallel(ipts):
    pid = ipts['pid']
    x1 = ipts['parent1']
    x2 = ipts['parent2']
    normalizationVector = ipts['normalizationVector']
    encodingTypes = ipts['encodingTypes']
    upperBounds = ipts['upperBounds'] 
    lowerBounds = ipts['lowerBounds'] 
    # formulation = ipts['formulation']
    children = breedDesignVectors(x1,x2,normalizationVector, encodingTypes, lowerBounds, upperBounds, Ncrossovers=3, probabilityOfMutation=0.3)
    return children

def newGeneration(fitnessFunction, population, normalizationVector, encodingTypes, lowerBounds, upperBounds, tau, processType='series', initalize=False):
    population = np.array(population)
    Nvars = len(normalizationVector)
    
    if len(population)%4!=0:
        raise ValueError('Population length must be evenly divisible by 4')

    if initalize:
        ipList = []
        for i in range(0,len(population)):
            ins = {}
            ins['pid'] = i
            ins['individual'] = population[i]
            ins['tau'] = tau
            # ins['formulation'] = formulation
            # ins['bounds'] = bounds
            # ins['constraints'] = constraints
            # ins['norms'] = normVector
            ipList.append(ins)
    
        if processType == 'parallel':
            pool = multiprocessing.Pool()  
            result = pool.map(fitnessFunction, ipList)
            pool.close()
        
        if processType == 'series':
            result = []
            for i in range(0,len(population)):
                result.append(fitnessFunction(ipList[i]))

        fitness = np.array([result[i][0] for i in range(0,len(result))])
        # constraintStatus = np.array([result[i][1] for i in range(0,len(result))]).astype(float)
        
        # print(fitness)
        data = np.append(population,np.array([fitness]).T,axis=1)
        # print(data)
        # data = np.append(data,np.array([constraintStatus]).T,axis=1)
        
        # data = np.append(resultantPop,np.array([fitness]).T,axis=1)
    
        # constraintStatus = np.array([result[i][1] for i in range(0,len(result))]).astype(float)
        # cl = np.array([result[i][2] for i in range(0,len(result))])

        for ii in range(1,len(result[0])):
            da = np.array([result[i][ii] for i in range(0,len(result))]).astype(float)
            data = np.append(data,np.array([da]).T,axis=1)
    
        sortedData = data[data[:,Nvars].argsort()]
        dataSave = copy.deepcopy(sortedData)

        return sortedData, dataSave

    # -------------------------------
    # Tournament Selection, no elites
    # -------------------------------
    sortedData = population[population[:,Nvars].argsort()]
    # sortedData = sortedData[0,0:-2]
    remainingMembers = list(range(0,len(sortedData)))
    breedingPop = []
    for i in range(0,int(len(remainingMembers)/2)):
        ix1 = random.randrange(0,len(remainingMembers))
        ix2 = ix1 - random.randrange(1,len(remainingMembers))
        pi1 = remainingMembers[ix1]
        pi2 = remainingMembers[ix2]
        parent1 = sortedData[pi1]
        parent2 = sortedData[pi2]
        if sortedData[pi1][Nvars] <= sortedData[pi2][Nvars]:
            breedingPop.append(parent1)
        else:
            breedingPop.append(parent2)

        remainingMembers.remove(pi1)
        remainingMembers.remove(pi2)
          
    ipList = []
    resultantPop = []
    for ii in [0,1]: #each parent must produce 4 offspring
        remainingMembers = list(range(0,len(breedingPop)))
        for i in range(0,int(len(remainingMembers)/2)):
            ins = {}
            ins['pid'] = (ii+1) * i
            ix1 = random.randrange(0,len(remainingMembers))
            ix2 = ix1 - random.randrange(1,len(remainingMembers))
            pi1 = remainingMembers[ix1]
            pi2 = remainingMembers[ix2]
            parent1 = breedingPop[pi1][0:Nvars]
            parent2 = breedingPop[pi2][0:Nvars]
            ins['parent1'] = copy.deepcopy(parent1)
            ins['parent2'] = copy.deepcopy(parent2)
            ins['normalizationVector']=normalizationVector
            ins['encodingTypes'] = encodingTypes
            ins['upperBounds'] = upperBounds
            ins['lowerBounds'] = lowerBounds
            # ins['formulation'] = copy.deepcopy(formulation)
            ipList.append(ins)
            remainingMembers.remove(pi1)
            remainingMembers.remove(pi2)

    if processType == 'parallel':
        pool = multiprocessing.Pool()  
        childrenList = pool.map(breedDesignVectorsParallel, ipList)
        pool.close()
        
    if processType == 'series':
        childrenList = []
        for i in range(0,len(ipList)):
            children = breedDesignVectorsParallel(ipList[i])
            childrenList.append(children)
    
    for i in range(0,len(childrenList)):
        resultantPop.append(childrenList[i][0])
        resultantPop.append(childrenList[i][1])

    # -------------------------------
    # Can include other types of selection in the future
    # -------------------------------
    ipList = []
    for i in range(0,len(resultantPop)):
        ins = {}
        ins['pid'] = i
        
        te_gap_lookup = {
            '15':  0.00196,
            '18':  0.00230,
            '21':  0.00262,
            '24':  0.00751,
            '27':  0.01012,
            '30':  0.01828,
            '33':  0.02644,
            '36':  0.02896,
        }

        afl = Kulfan(TE_gap = te_gap_lookup[str(int(100*tau))])
        K = resultantPop[i]
        Ku = K[0:int(len(K)/2)]
        Kl = K[int(len(K)/2):]
        afl.upperCoefficients = Ku
        afl.lowerCoefficients = Kl
        afl.scaleThickness(tau)
        resultantPop[i] = afl.upperCoefficients.magnitude.tolist() + afl.lowerCoefficients.magnitude.tolist()
        
        ins['individual'] = resultantPop[i]
        ins['tau'] = tau
        # ins['formulation'] = formulation
        # ins['bounds'] = bounds
        # ins['constraints'] = constraints
        # ins['norms'] = normVector
        ipList.append(ins)
    
    if processType == 'parallel':
        pool = multiprocessing.Pool()  
        result = pool.map(fitnessFunction, ipList)
        pool.close()
    
    if processType == 'series':
        result = []
        for i in range(0,len(resultantPop)):
            result.append(fitnessFunction(ipList[i]))


    assert(len(result) == len(population))
            
    fitness = np.array([result[i][0] for i in range(0,len(result))])
    data = np.append(resultantPop,np.array([fitness]).T,axis=1)
    
    # constraintStatus = np.array([result[i][1] for i in range(0,len(result))]).astype(float)
    # cl = np.array([result[i][2] for i in range(0,len(result))])

    for ii in range(1,len(result[0])):
        da = np.array([result[i][ii] for i in range(0,len(result))]).astype(float)
        data = np.append(data,np.array([da]).T,axis=1)
    

    # data = np.append(data,np.array([constraintStatus]).T,axis=1)
    # data = np.append(data,np.array([cl]).T,axis=1)
    sortedData = data[data[:,Nvars].argsort()]
    dataSave = copy.deepcopy(sortedData)

    return sortedData, dataSave

# from temp import ff

# pop = [[random.uniform(-10,10)] for i in range(0,10)]
# print(pop)
# pop,tsh = newGeneration(ff, pop, normalizationVector = [1], encodingTypes=[float], lowerBounds=[-10], upperBounds=[10], processType='parallel', initalize=True)
# print(pop)
# for i in range(0,15):
#     pop,tsh = newGeneration(ff, pop, normalizationVector = [1], encodingTypes=[float], lowerBounds=[-10], upperBounds=[10], processType='series', initalize=False)
#     print(pop)