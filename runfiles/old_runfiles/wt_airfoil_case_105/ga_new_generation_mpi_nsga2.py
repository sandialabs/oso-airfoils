import copy
import random
import numpy as np
import math
from kulfan import Kulfan

from newMember import newMember

from ga_functions import ( numberCoding, 
                           encodeChromosome, 
                           decodeChromosome, 
                           mutateChromosome, 
                           crossoverChromosomes, 
                           breedDesignVectors, 
                           breedDesignVectorsParallel,
                           NSGA_sort,
                           crowding_sort,
                           )

def newGeneration(  fitnessFunction, 
                    population, 
                    normalizationVector, 
                    encodingTypes, 
                    lowerBounds, 
                    upperBounds, 
                    tau, 
                    initalize=False,
                    comm = None  ,
                    CL = None,
                    re = None):

    population = np.array(population)
    if len(population)%4!=0:
        raise ValueError('Population length must be evenly divisible by 4')
    Nvars = len(normalizationVector)

    size = comm.Get_size()
    rank = comm.Get_rank()
    # print(rank, size)

    if initalize:
        # ipList = []
        result = []
        for i in range(0,len(population)):
            if i%size == rank:
                ins = {}
                ins['pid'] = i
                ins['individual'] = population[i]
                ins['tau'] = tau
                ins['CL'] = CL
                # ins['rLD'] = rLD
                ins['re'] = re
                # ins['file_system'] = file_system
                # ins['N_retries'] = N_retries
                result.append(fitnessFunction(ins))

        result = comm.gather(result, root=0)

        if rank == 0:
            result_temp = [None]*len(population)
            for i in range(0,len(result)):
                for j in range(0,len(result[i])): 
                    result_temp[int(result[i][j][0])] = result[i][j]
                    # result_temp.append(result[i][j])
            result = result_temp
            fitness = np.array([result[i][1] for i in range(0,len(result))])
            data = np.append(population,np.array([fitness]).T,axis=1)
            for ii in range(2,len(result[0])):
                da = np.array([result[i][ii] for i in range(0,len(result))]).astype(float)
                data = np.append(data,np.array([da]).T,axis=1)

            sortedData = NSGA_sort(data, Nvars, 2)
            #add a column of zeros for the crowding distance since this will show up later
            # crowding_distances = np.zeros((sortedData.shape[0], 1))
            # sortedData = np.hstack((sortedData, crowding_distances))

            # sortedData = data[data[:,Nvars].argsort()]
            return np.array(sortedData)
        else:
            return None



    # -------------------------------
    # Tournament Selection, no elites
    # -------------------------------
    if rank == 0:
        sortedData = NSGA_sort(population, Nvars, 2)
        #add a column of zeros for the crowding distance since this will show up later
        # crowding_distances = np.zeros((sortedData.shape[0], 1))
        # sortedData = np.hstack((sortedData, crowding_distances))
        # sortedData = population[population[:,Nvars].argsort()]
        remainingMembers = list(range(0,len(sortedData)))
        breedingPop = []
        for i in range(0,int(len(remainingMembers)/2)):
            ix1 = random.randrange(0,len(remainingMembers))
            ix2 = ix1 - random.randrange(1,len(remainingMembers))
            if ix2 < 0:
                ix2 += len(remainingMembers)
            pi1 = remainingMembers[ix1]
            pi2 = remainingMembers[ix2]
            parent1 = sortedData[pi1]
            parent2 = sortedData[pi2]
            if ix1 <= ix2:
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
                ipList.append(ins)
                remainingMembers.remove(pi1)
                remainingMembers.remove(pi2)

        childrenList = []
        for i in range(0,len(ipList)):
            children = breedDesignVectorsParallel(ipList[i])
            childrenList.append(children)
        
        for i in range(0,len(childrenList)):
            resultantPop.append(childrenList[i][0])
            resultantPop.append(childrenList[i][1])

        ############
        ############
        ############
        # This bit adjusts so that all candidates have correct thickness
        # Necessary for searches to be feasible
        # Can in theory comment out
        ############
        for i in range(0,len(resultantPop)):
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
            if any([math.isnan(kv) for kv in K]):
                N_k = len(K)
                K = newMember(int(N_k/2),tau,1)[0]#[random.uniform(0.1,0.8)] + [random.uniform(-0.1,0.8) for j in range(0,int(N_k/2)-1)]   +   [random.uniform(-0.8,-0.1)] + [random.uniform(-0.8,0.4) for j in range(0,int(N_k/2)-1)]  #newMember(N_k,tau)
            Ku = K[0:int(len(K)/2)]
            Kl = K[int(len(K)/2):]
            afl.upperCoefficients = Ku
            afl.lowerCoefficients = Kl

            try:
                afl.scaleThickness(tau)
                resultantPop[i] = afl.upperCoefficients.magnitude.tolist() + afl.lowerCoefficients.magnitude.tolist()
                if any([math.isnan(rpv) for rpv in resultantPop[i]]):
                    N_k = len(resultantPop[i])
                    resultantPop[i] = newMember(int(N_k/2),tau,1)[0]#[random.uniform(0.1,0.8)] + [random.uniform(-0.1,0.8) for j in range(0,int(N_k/2)-1)]   +   [random.uniform(-0.8,-0.1)] + [random.uniform(-0.8,0.4) for j in range(0,int(N_k/2)-1)]  #newMember(N_k,tau)
            except:
                N_k = len(K)
                resultantPop[i] = newMember(int(N_k/2),tau,1)[0]#[random.uniform(0.1,0.8)] + [random.uniform(-0.1,0.8) for j in range(0,int(N_k/2)-1)]   +   [random.uniform(-0.8,-0.1)] + [random.uniform(-0.8,0.4) for j in range(0,int(N_k/2)-1)]  #newMember(N_k,tau)    

            assert(not any([math.isnan(rpv) for rpv in resultantPop[i]]))


        ############
        ############
        ############

    else:
        resultantPop = None

    resultantPop = comm.bcast(resultantPop,root=0)  

    # -------------------------------
    # Can include other types of selection in the future
    # -------------------------------
    # ipList = []
    result = []
    for i in range(0,len(resultantPop)):
        if i%size == rank:
            ins = {}
            ins['pid'] = i
            ins['individual'] = resultantPop[i]
            ins['tau'] = tau
            ins['CL'] = CL
            ins['rLD'] = rLD
            ins['re'] = re
            ins['file_system'] = file_system
            ins['N_retries'] = N_retries
            result.append(fitnessFunction(ins))

    result = comm.gather(result, root=0)

    if rank == 0:
        result_temp = [None]*len(resultantPop)
        for i in range(0,len(result)):
            for j in range(0,len(result[i])):
                result_temp[int(result[i][j][0])] = result[i][j]
                # result_temp.append(result[i][j])
        result = result_temp
        fitness = np.array([result[i][1] for i in range(0,len(result))])
        data = np.append(resultantPop,np.array([fitness]).T,axis=1)
        for ii in range(2,len(result[0])):
            da = np.array([result[i][ii] for i in range(0,len(result))]).astype(float)
            data = np.append(data,np.array([da]).T,axis=1)

        data_combined = np.append(population[:,0:-1], data, axis=0)


        # do the nsga thing...
        # for i in range(0,len(population)):
        #     resultantPop.append(population[i][0:Nvars])
        
        resultantPop = np.array( NSGA_sort(data_combined,Nvars,2) )
        rpop_length = len(resultantPop)

        # get frequency of numbers in last column of resultantPop
        front_numbers = resultantPop[:, -1]
        unique, counts = np.unique(front_numbers, return_counts=True)
        front_count_dict = dict(zip(unique, counts))
        total_num = 0
        finalPop = []
        last_front = None
        for i in range(1,max(front_numbers.astype(int))+1):
            if total_num + front_count_dict[i] >= rpop_length:
                last_front = i
                break
            else:
                finalPop.append( resultantPop[front_numbers == i] )
                total_num += front_count_dict[i]

        remaining_slots = rpop_length - total_num
        # get all the individuals in the last front
        last_front_inds = resultantPop[front_numbers == last_front]
        last_front_inds = crowding_sort(last_front_inds, Nvars, 2)
        # take the top of the list and add to finalPop
        for i in range(0,remaining_slots):
            finalPop.append( last_front_inds[i] )

        return finalPop

        # sortedData = data[data[:,Nvars].argsort()]
        # dataSave = copy.deepcopy(sortedData)
        # return sortedData#, dataSave
    else:
        return None

