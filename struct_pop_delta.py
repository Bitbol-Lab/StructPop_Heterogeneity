# -*- coding: utf-8 -*-
'''
Code for simulating the serial dilution model.
We consider a population with D demes and N_types types of individuals. The population composition is represented as an array e.g. [[1, 99], [0,100]] has 2 demes ; the first one has 1 mutant and 99 wild-types, the second has 100 wild-types.
The graph structure is implemented through the migration matrix.
'''

import numpy as np
from numba import njit
from numba import set_num_threads
set_num_threads(1)


#%%______________________________________Migration matrices________________________________________________


def define_clique(D,m):
    '''
    Function to define the clique migration matrix.

    Parameters
    ----------
    D : int
        Number of demes.
    m : float
        Per-capita migration rate.

    Returns
    -------
    migration_matrix : array
        Migration matrix for the clique.
    '''
    migration_matrix=np.zeros((D,D))
    for i in range(D):
        for j in range(D):
            if not i==j :
                migration_matrix[i,j]=m
        migration_matrix[i,i]=1-(D-1)*m
    return migration_matrix


def define_cycle(D,mA,mC):
    '''
    Function to define the cycle migration matrix.

    Parameters
    ----------
    D : int
        Number of demes.
    mA : float
        Migration rate, anticlockwise direction.
    mC : float
        Migration rate, clockwise direction.

    Returns
    -------
    migration_matrix : array
        Migration matrix for the cycle.
    '''
    migration_matrix=np.zeros((D,D))
    for i in range(D-1):
        migration_matrix[i,i+1]=mA
        migration_matrix[i+1,i]=mC
    migration_matrix[D-1,0]=mA
    migration_matrix[0,D-1]=mC
    for j in range(D):
        migration_matrix[j,j]=1-mA-mC
    return migration_matrix


def define_star(D,mI,mO,equal_contribution=True):
    '''
    Function to define the star migration matrix.

    Parameters
    ----------
    D : int
        Number of demes.
    mI : float
        Inward migration rate.
    mO : float
        Outward migration rate.
    equal_contribution : bool, optional
        If True: sum on lines = 1 (all demes send the same number of individuals on average).
        If False: sum on columns = 1 (all demes receive the same number on average).
        The default is True.

    Returns
    -------
    migration_matrix : array
        Migration matrix for the star.

    '''
    migration_matrix=np.zeros((D,D))
    for i in range(1,D):
        migration_matrix[i,0]=mI
        migration_matrix[0,i]=mO
        if equal_contribution:
            migration_matrix[i,i]=1-mI
        else:
            migration_matrix[i,i]=1-mO
    if equal_contribution:
        migration_matrix[0,0]=1-(D-1)*mO 
    else:
        migration_matrix[0,0]=1-(D-1)*mI 
    return migration_matrix


def random_graph(adj_matrix,convention,alpha_values=None):
    '''
    Generates a random migration matrix from a Dirichlet distribution, with alpha values given as parameters.
    Convention is either 'eqsize' (demes have size K on average) or 'eqcon' (demes contribute by K on average).

    Parameters:
    ----------
    adj_matrix : array
        Adjacency matrix representing the graph structure.
    convention : str
        Either 'eqsize' or 'eqcon', specifying the type of equilibrium to consider.
    alpha_values : array, optional
        Alpha values for sampling the Dirichlet distribution. If None, uniform values are used.

    Returns
    -------
    M : array
        Random migration matrix.
    '''
    D=np.shape(adj_matrix)[0]
    M=np.zeros((D,D))
    if alpha_values is None:
        alphas=np.ones((D,D))
    else:
        alphas=alpha_values

    neighbour_out=np.array([np.argwhere(adj_matrix[i,:]==1).flatten() for i in range(D)],dtype=object)
    neighbour_in=np.array([np.argwhere(adj_matrix[:,i]==1).flatten() for i in range(D)],dtype=object)
    
    if convention=='eqsize':
        for i in range(D):
            alphas_in=np.array([alphas[j,i] for j in neighbour_in[i]])
            sample_rates=np.random.dirichlet(alphas_in)
            M[neighbour_in[i].astype(np.int64),i]=sample_rates

    elif convention=='eqcon':
        for i in range(D):
            alphas_out=np.array([alphas[i,j] for j in neighbour_out[i]])
            sample_rates=np.random.dirichlet(alphas_out)
     
            M[i,neighbour_out[i]]=sample_rates
                
    return M


def define_line(D,mR,mL, equal_contribution = True):
    '''
    Function to define the line migration matrix.

    Parameters
    ----------
    D : int
        Number of demes.
    mR : float
        Rightward migration rate.
    mL : float
        Leftward migration rate.
    equal_contribution : bool, optional
        If True: sum on lines = 1 (all demes send the same number of individuals on average).
        If False: sum on columns = 1 (all demes receive the same number on average).
        The default is True.

    Returns
    -------
    migration_matrix : array
        Migration matrix for the line.
    '''
    migration_matrix=np.zeros((D,D))
    for i in range(D-1):
        migration_matrix[i,i+1]=mR
        migration_matrix[i+1,i]=mL
    for j in range(1, D-1):
        migration_matrix[j,j]=1-mR-mL
        
    if equal_contribution:
        migration_matrix[0,0] = 1-mR
        migration_matrix[D-1,D-1] = 1-mL
    else:
        migration_matrix[0,0] = 1-mL
        migration_matrix[D-1,D-1] = 1-mR

    return migration_matrix



#%%_________________________________________Growth and dilution/migration_____________________________________________

@njit
def growth_event(in_numbers,fitnesses,t):
    '''
    Describe the deterministic growth event for time t.

    Parameters
    ----------
    in_numbers : List or array
        Initial demes' compositions before the growth phase.
    fitnesses : list or array
        list with the mutant and wild-type fitnesses, i.e., [1, 1+s]. It can also be a list of lists (in gradient case).
    t : int
        Time duration of the growth event.

    Returns
    -------
    Array
        Returns the number of individuals (same shape as in_numbers) after the growth time t.
    '''
    return (in_numbers * np.exp(fitnesses * t)).astype(np.float64)


@njit
def dilution_migration_event(in_numbers,migration_matrix,K):
    '''
    Dilution and migration event.

    Parameters
    ----------
    in_numbers : array or list
        Composition of the demes at the end of the growth phase.
    migration_matrix : array
        Migration matrix, one of those defined above with parameters specific for my case.
    K : int
        Bottleneck number.

    Returns
    -------
    new_numbers : array 
        Composition of the demes after migration and dilution steps.
    '''
    D, N_types= np.shape(in_numbers)

    new_numbers=np.zeros((D, N_types), dtype=np.int64)

    migrants_ij=np.empty(2, dtype=np.int64)
    
    for i in range(D):
        Ni=np.sum(in_numbers[i,:], dtype = np.int64)
        if Ni<1:
            print('extinct deme', i)
            
        p=in_numbers[i,0]/Ni
        for j in np.arange(D):
            mij=migration_matrix[i,j]

            p0=np.float64(max(min(K*p*mij/Ni,1),0))
            p1=np.float64(max(min(K*(1-p)*mij/Ni,1),0))
            migrants_ij[0]=np.random.binomial(Ni, p0,1)[0]
            migrants_ij[1]=np.random.binomial(Ni, p1,1)[0]

            new_numbers[j, 0]+=migrants_ij[0]
            new_numbers[j, 1]+=migrants_ij[1]

    return new_numbers

#%%_________________________________________Extinction or fixation____________________________________________________________

@njit
def extinct_mutant(numbers):
    ''' Check if the mutant population is extinct. '''
    D = numbers.shape[0]
    for i in range(D):
        if numbers[i,0]>0:
            return False
    return True

@njit
def extinct_wild(numbers):
    ''' Check if the wild-type population is extinct. '''
    D, N_types= numbers.shape
    for i in range(D):
        for j in range(N_types -1):
            if numbers[i,j+1]>0:
                return False
    return True

#%%____________________________________Complete simulation for one trajectory_____________________________________________

@njit
def cycle(in_numbers, migration_matrix, fitnesses, nb_cycles, growth_factor, K, start_follow_numbers, size_follow_numbers, start_cycle, print_frequency):
    ''' Simulate one cycle of serial dilutions.

    Parameters
    ----------
    in_numbers : array
        Conditions of the system before beginning cycle.
    migration_matrix : array
        Matrix containing the migration rates for the spatial structure.
    fitnesses : array
        Fitnesses of mutants and wild-types in the demes.
    nb_cycles : int
        A priori number of how many cycles we do before observing fixation or extinction.
    growth_factor : int
        t, the deterministic growth duration.
    K : inte
        Bottleneck size.
    start_follow_numbers : array
        The array of initial start numbers to follow. Can be also 'None'.
    size_follow_numbers : int
        Size of the output array we want to save.
    start_cycle : int
        Number of the cycle at which we want to start saving the dynamics.
    print_frequency : int
        The frequency at which we print the state.
    save_dynamics : bool, optional
        If True: saves the trajectories
        If False: does not save the trajectories.
        The default is False.

    Returns
    -------
    follow_numbers : array
        The array containing the numbers tracked over the cycles.
    end_cycle : int
        The cycle at which the simulation ended.
    fixation : bool
        Whether fixation occurred (True) or not (False).
    '''
    D, N_types= np.shape(in_numbers)

    if start_follow_numbers is None:
        follow_numbers=np.zeros((size_follow_numbers,D,N_types), dtype=np.int64)
    else:
        follow_numbers=start_follow_numbers.copy()

    fixation=True
    numbers=in_numbers.copy()

    #Booleans that checks if mutants or wild-types are extinct
    keep_going=True

    for i in range(nb_cycles):
        end_cycle = nb_cycles
        numbers1=growth_event(numbers,fitnesses,growth_factor)
        numbers=dilution_migration_event(numbers1,migration_matrix,K)

        if (start_cycle+i)%print_frequency==0 and ((i+start_cycle)/print_frequency)<size_follow_numbers:
            follow_numbers[int(i+start_cycle), :, :]=numbers

        if extinct_mutant(numbers):
            keep_going=False
            fixation=False
            end_cycle=i+start_cycle
            follow_numbers[int((i+start_cycle))+1:]=numbers
            break


        if extinct_wild(numbers):
            keep_going=False
            fixation=True
            end_cycle=i+start_cycle
            follow_numbers[int((i+start_cycle))+1:]=numbers
            break
    
    #If mutants are not extinct or fixed at the end of nb_cycles cycles, we keep going
    if keep_going:
        follow_numbers, end_cycle, fixation= cycle(numbers, migration_matrix, fitnesses, nb_cycles, growth_factor, K, follow_numbers, size_follow_numbers, start_cycle+end_cycle, print_frequency)
    return follow_numbers, end_cycle, fixation


#%%________________________________________Fixation probability computed on several simulations_______________________________________________

@njit
def fixation_probability(in_numbers, migration_matrix, fitnesses, nb_sim, nb_cycles, growth_factor, K, size_follow_numbers=10000, print_frequency=1, save_dynamics=False):
    ''' Computes the fixation probability over several simulations, using the 'cycle' function as of above, starting from in_numbers, on nb_sim simulations.
    Parameters (see above).

    Returns:
    - average_extinction_cycle : float
        The average time at which extinction occurred.
    - ci95_ec : float
        The 95% confidence interval for the average extinction time.
    - average_fixation_cycle : float
        The average time at which fixation occurred.
    - ci95_fc : float
        The 95% confidence interval for the average fixation time.
    - proba : float
        The fixation probabilitys.
    '''
    fix_count=0

    fix_cycle=np.zeros(nb_sim)
    ex_cycle=np.zeros(nb_sim)

    for i in range(nb_sim):
        start_cycle=0
        start_follow_numbers=None
        follow_numbers, end_cycle, fixation = cycle(in_numbers, migration_matrix, fitnesses, nb_cycles, growth_factor, K, start_follow_numbers, size_follow_numbers, start_cycle, print_frequency)
        
        if fixation :
            fix_count+=1
            fix_cycle[i] = end_cycle

        else:
            ex_cycle[i] = end_cycle

    #Number of extinctions
    ex_count=nb_sim-fix_count

    if fix_count>0:
        average_fixation_cycle = np.sum(fix_cycle)/fix_count
        variance_fixation_cycle = np.var(fix_cycle)
        ci95_fc = 1.96 * np.sqrt(variance_fixation_cycle) / np.sqrt(fix_count)
    else :
        average_fixation_cycle = 0
    if ex_count>0:
        average_extinction_cycle = np.sum(ex_cycle)/ex_count
        variance_extinction_cycle = np.var(ex_cycle)
        ci95_ec = 1.96 * np.sqrt(variance_extinction_cycle) / np.sqrt(ex_count)
    else:
        average_extinction_cycle=0

    proba = fix_count/nb_sim

    return average_extinction_cycle, ci95_ec, average_fixation_cycle, ci95_fc, proba


