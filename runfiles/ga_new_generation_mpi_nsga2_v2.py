"""
NSGA-II generation step with two fixes applied (versus
``ga_new_generation_mpi_nsga2.py``):

  Fix #1 — Design-vector de-duplication of the combined (parent + child)
           population BEFORE non-dominated sorting. Identical genomes (which
           NSGA_sort cannot distinguish, since it only looks at objective
           columns) were piling up in front 1 every generation. This also
           obviates fix #3 (clone tie-breaking in ``crowding_sort``): once
           the combined pool contains no duplicate design vectors,
           ``crowding_sort`` never sees clones at identical objective
           coordinates, so its arbitrary tie-break ordering has no
           opportunity to fix duplicates into the next generation.

  Fix #2 — Hard cap on front-1 occupancy in the next population. Once the
           Pareto front grows past ``front1_cap_fraction * N_pop`` (default
           0.5), the multi-front structure that drives NSGA-II selection
           collapses and the algorithm degenerates to crowding-distance
           selection on a single front (which, in the presence of even a
           few near-clones, runs away into the diversity loss observed in
           the collapsed run). The cap forces excess front-1 members to be
           displaced by lower fronts; if not enough lower-front material
           exists, the worst (by crowding) front-1 members fill the
           remaining slots so the population size is preserved.

Drop-in usage: change ``from ga_new_generation_mpi_nsga2 import newGeneration``
to ``from ga_new_generation_mpi_nsga2_v2 import newGeneration``.

Optional new parameters in ``params``:
    front1_cap_fraction : float in (0, 1], default 0.5
        Maximum fraction of the next population allowed to come from front 1.
"""

import copy
import random
import sys
import numpy as np
import math
from kulfan import Kulfan

def cprint(x):
    sys.stdout.flush()
    print(x)

from geometry_functions import TE_gap_function
from newMember import newMember

from ga_functions import (
    numberCoding,
    encodeChromosome,
    decodeChromosome,
    mutateChromosome,
    crossoverChromosomes,
    breedDesignVectors,
    breedDesignVectorsParallel,
    NSGA_sort,
    crowding_sort,
)


def _dedup_by_design_vector(data, Nvars, decimals=12):
    """Keep only the FIRST occurrence of each unique design vector.

    ``data`` rows are assumed to have layout::

        [ design vars (Nvars) | objectives | extras | front_number? ]

    Order is preserved (so if ``data`` was already pareto-sorted, the
    representative kept for each duplicate group is the best-ranked one).
    """
    if len(data) == 0:
        return data
    arr = np.asarray(data, dtype=float)
    seen = set()
    keep = []
    for i in range(arr.shape[0]):
        key = tuple(np.round(arr[i, 0:Nvars], decimals).tolist())
        if key in seen:
            continue
        seen.add(key)
        keep.append(i)
    return arr[keep]


def newGeneration(  fitnessFunction,
                    population,
                    normalizationVector,
                    encodingTypes,
                    lowerBounds,
                    upperBounds,
                    initalize=False,
                    comm=None,
                    params=None,
                    ):

    tau = params['tau']
    if 'TE_gap' in params:
        TE_gap = params['TE_gap']
    else:
        TE_gap = TE_gap_function(tau)

    if 'N_crossovers' in params:
        N_crossovers = max([int(params['N_crossovers']), 1])
    else:
        N_crossovers = 3

    if 'probability_of_mutation' in params:
        probability_of_mutation = max([min([1.0, params['probability_of_mutation']]), 0.0])
    else:
        probability_of_mutation = 0.3

    if 'maximum_parent_fraction' in params:
        maximum_parent_fraction = max([min([1.0, params['maximum_parent_fraction']]), 0.0])
    else:
        maximum_parent_fraction = 0.7

    # ---- Fix #2 parameter ----
    if 'front1_cap_fraction' in params and params['front1_cap_fraction'] is not None:
        front1_cap_fraction = max([min([1.0, float(params['front1_cap_fraction'])]), 1e-6])
    else:
        front1_cap_fraction = 0.5

    if 'N_mutations' in params and params['N_mutations'] is not None:
        N_mutations = max([int(params['N_mutations']), 1])
    else:
        N_mutations = 1

    population = np.array(population)
    if len(population) % 4 != 0:
        raise ValueError('Population length must be evenly divisible by 4')
    Nvars = len(normalizationVector)

    size = comm.Get_size()
    rank = comm.Get_rank()

    if initalize:
        result = []
        for i in range(0, len(population)):
            if i % size == rank:
                ins = {}
                ins['pid'] = i
                ins['individual'] = population[i]
                ins['params'] = params
                result.append(fitnessFunction(ins))

        result = comm.gather(result, root=0)

        if rank == 0:
            result_temp = [None] * len(population)
            for i in range(0, len(result)):
                for j in range(0, len(result[i])):
                    result_temp[int(result[i][j][0])] = result[i][j]
            result = result_temp
            fitness = np.array([result[i][1] for i in range(0, len(result))])
            data = np.append(population, np.array([fitness]).T, axis=1)

            for ii in range(2, len(result[0])):
                da = np.array([result[i][ii] for i in range(0, len(result))]).astype(float)
                data = np.append(data, np.array([da]).T, axis=1)

            sortedData = NSGA_sort(data, Nvars, 2)
            return np.array(sortedData)
        else:
            return None

    # -------------------------------
    # Tournament Selection, no elites
    # -------------------------------
    if rank == 0:
        sortedData = np.array(NSGA_sort(population, Nvars, 2))
        remainingMembers = list(range(0, len(sortedData)))
        breedingPop = []
        for i in range(0, int(len(remainingMembers) / 2)):
            ix1 = random.randrange(0, len(remainingMembers))
            ix2 = ix1 - random.randrange(1, len(remainingMembers))
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
        for ii in [0, 1]:  # each parent must produce 4 offspring
            remainingMembers = list(range(0, len(breedingPop)))
            for i in range(0, int(len(remainingMembers) / 2)):
                ins = {}
                ins['pid'] = (ii + 1) * i
                ix1 = random.randrange(0, len(remainingMembers))
                ix2 = ix1 - random.randrange(1, len(remainingMembers))
                pi1 = remainingMembers[ix1]
                pi2 = remainingMembers[ix2]
                parent1 = breedingPop[pi1][0:Nvars]
                parent2 = breedingPop[pi2][0:Nvars]
                ins['parent1'] = copy.deepcopy(parent1)
                ins['parent2'] = copy.deepcopy(parent2)
                ins['normalizationVector'] = normalizationVector
                ins['encodingTypes'] = encodingTypes
                ins['upperBounds'] = upperBounds
                ins['lowerBounds'] = lowerBounds
                ins['Ncrossovers'] = N_crossovers
                ins['probabilityOfMutation'] = probability_of_mutation
                ins['N_mutations'] = N_mutations
                ipList.append(ins)
                remainingMembers.remove(pi1)
                remainingMembers.remove(pi2)

        childrenList = []
        for i in range(0, len(ipList)):
            children = breedDesignVectorsParallel(ipList[i])
            childrenList.append(children)

        for i in range(0, len(childrenList)):
            resultantPop.append(childrenList[i][0])
            resultantPop.append(childrenList[i][1])

        ############
        # Thickness re-scaling of children
        ############
        for i in range(0, len(resultantPop)):

            afl = Kulfan(TE_gap=TE_gap)
            K = resultantPop[i]
            if any([math.isnan(kv) for kv in K]):
                N_k = len(K)
                K = newMember(int(N_k / 2), tau, 1)[0]
            Ku = K[0:int(len(K) / 2)]
            Kl = K[int(len(K) / 2):]
            afl.upperCoefficients = Ku
            afl.lowerCoefficients = Kl

            try:
                afl.scaleThickness(tau)
                resultantPop[i] = afl.upperCoefficients.magnitude.tolist() + afl.lowerCoefficients.magnitude.tolist()
                if any([math.isnan(rpv) for rpv in resultantPop[i]]):
                    N_k = len(resultantPop[i])
                    resultantPop[i] = newMember(int(N_k / 2), tau, 1)[0]
            except Exception:
                N_k = len(K)
                resultantPop[i] = newMember(int(N_k / 2), tau, 1)[0]

            assert (not any([math.isnan(rpv) for rpv in resultantPop[i]]))
    else:
        resultantPop = None

    resultantPop = comm.bcast(resultantPop, root=0)

    # -------------------------------
    # Fitness evaluation of children (distributed)
    # -------------------------------
    result = []
    for i in range(0, len(resultantPop)):
        if i % size == rank:
            ins = {}
            ins['pid'] = i
            ins['individual'] = resultantPop[i]
            ins['params'] = params
            result.append(fitnessFunction(ins))

    result = comm.gather(result, root=0)

    if rank == 0:
        result_temp = [None] * len(resultantPop)
        for i in range(0, len(result)):
            for j in range(0, len(result[i])):
                result_temp[int(result[i][j][0])] = result[i][j]
        result = result_temp
        fitness = np.array([result[i][1] for i in range(0, len(result))])
        data = np.append(resultantPop, np.array([fitness]).T, axis=1)
        for ii in range(2, len(result[0])):
            da = np.array([result[i][ii] for i in range(0, len(result))]).astype(float)
            data = np.append(data, np.array([da]).T, axis=1)

        # Combine top fraction of parents with all children
        last_parent_index = int(maximum_parent_fraction * len(population))
        data_combined = np.append(sortedData[0:last_parent_index, 0:-2], data, axis=0)

        # -------------------------------------------------------------
        # Fix #1: dedup by design vector before non-dominated sorting.
        # sortedData (parents) is already pareto-ordered, and children
        # follow; np.append preserves that order, so keeping the first
        # occurrence of each genome retains the best-ranked copy.
        # -------------------------------------------------------------
        n_before = len(data_combined)
        data_combined = _dedup_by_design_vector(data_combined, Nvars)
        n_after = len(data_combined)
        if n_after < n_before:
            cprint(f"[v2] dedup removed {n_before - n_after} duplicate genomes "
                   f"({n_after} unique going into NSGA_sort)")

        resultantPop = np.array(NSGA_sort(data_combined, Nvars, 2))
        rpop_length = len(population)

        # Front-number bookkeeping
        front_numbers = resultantPop[:, -1]
        unique, counts = np.unique(front_numbers, return_counts=True)
        front_count_dict = dict(zip(unique.astype(int).tolist(), counts.astype(int).tolist()))
        max_front = int(front_numbers.max())

        # -------------------------------------------------------------
        # Fix #2: cap front-1 occupancy at front1_cap_fraction * N_pop.
        # Excess front-1 members are displaced by subsequent fronts;
        # if those don't supply enough rows, the worst (by crowding)
        # front-1 members refill the remainder so |finalPop| == N_pop.
        # -------------------------------------------------------------
        front1_cap = max(1, int(math.ceil(front1_cap_fraction * rpop_length)))
        front1_mask = (front_numbers == 1)
        front1_inds = resultantPop[front1_mask]

        finalPop = []
        front1_excess_sorted = None  # crowding-sorted leftovers if cap fires

        if len(front1_inds) > front1_cap:
            cprint(f"[v2] front-1 cap engaged: |F1|={len(front1_inds)} > cap={front1_cap}")
            f1_sorted = np.array(crowding_sort(front1_inds, Nvars, 2))
            # crowding_sort appends Nobj (=2) columns; drop them to restore
            # the [.. , front_number] layout used by the rest of the loop.
            f1_kept = f1_sorted[:front1_cap, 0:-2]
            front1_excess_sorted = f1_sorted[front1_cap:, 0:-2]
            for row in f1_kept:
                finalPop.append(row.tolist())
            total_num = front1_cap
            next_front_start = 2
        else:
            for row in front1_inds:
                finalPop.append(row.tolist())
            total_num = len(front1_inds)
            next_front_start = 2

        # Pack subsequent fronts whole until adding the next one would
        # overflow rpop_length.
        last_front = None
        for i in range(next_front_start, max_front + 1):
            cnt = front_count_dict.get(i, 0)
            if cnt == 0:
                continue
            if total_num + cnt >= rpop_length:
                last_front = i
                break
            this_front = resultantPop[front_numbers == i]
            for tf in this_front:
                finalPop.append(tf.tolist())
            total_num += cnt

        remaining_slots = rpop_length - total_num
        if remaining_slots > 0:
            if last_front is not None:
                # Standard NSGA-II: crowding-sort the splitting front.
                last_front_inds = resultantPop[front_numbers == last_front]
                last_front_inds = np.array(crowding_sort(last_front_inds, Nvars, 2))
                for i in range(0, remaining_slots):
                    finalPop.append(last_front_inds[i][0:-2].tolist())
            elif front1_excess_sorted is not None:
                # All other fronts were absorbed but rpop_length not met
                # because the front-1 cap fired. Refill from the displaced
                # front-1 excess (already crowding-sorted, best first).
                take = min(remaining_slots, len(front1_excess_sorted))
                for i in range(take):
                    finalPop.append(front1_excess_sorted[i].tolist())
                if take < remaining_slots:
                    # Defensive: should not occur because total excess >=
                    # |F1| - cap, and we only reach here if combined pool
                    # was entirely front 1, meaning |F1| >= rpop_length.
                    cprint(f"[v2] WARNING: short by {remaining_slots - take} rows "
                           f"after refill; padding with front-1 best.")
                    pad_src = np.array(crowding_sort(front1_inds, Nvars, 2))[:, 0:-2]
                    for i in range(remaining_slots - take):
                        finalPop.append(pad_src[i % len(pad_src)].tolist())
            else:
                # No splitting front, no excess: combined pool was smaller
                # than rpop_length (only possible if dedup removed enough
                # rows that the pool itself is too small). Pad by repeating
                # the best individuals.
                cprint(f"[v2] WARNING: combined pool smaller than N_pop after dedup; "
                       f"padding {remaining_slots} rows from best ranked.")
                for i in range(remaining_slots):
                    finalPop.append(resultantPop[i % len(resultantPop)].tolist())

        assert len(finalPop) == rpop_length, (
            f"finalPop size {len(finalPop)} != N_pop {rpop_length}"
        )

        return finalPop

    else:
        return None
