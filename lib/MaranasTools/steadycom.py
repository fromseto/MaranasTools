#/usr/bin/python

# author: Lin Wang, Joshua Chan
"""
 SteadyCom 

"""
import pulp
import json
from fba_tools.fba_toolsClient import fba_tools
import os
import pandas as pd
import re

def parse_reactant(reactant, sign):
    """
    sign should be -1 or 1
    returns {'stoich': int, 'cpd': string, 'compartment': string}
    """
    m = re.match('\((?P<stoich>.+)\)\s*(?P<cpd>[a-zA-Z0-9_]+)\[(?P<compartment>[a-zA-Z0-9_]+)\]', reactant)
    if not m:
        raise ValueError("can't parse {}".format(reactant))
    ret_val = m.groupdict()
    ret_val['stoich'] = float(ret_val['stoich']) * sign
    return ret_val

def parse_equation(equation):
    left_side, right_side = re.split('\s*<?=>?\s*', equation)
    
    reactants = list()
    if left_side:
        left_cpds = re.split('\s+\+\s+', left_side)
        reactants = reactants + [parse_reactant(r, -1) for r in left_cpds]
    if right_side:
        right_cpds = re.split('\s+\+\s+', right_side)
        reactants = reactants + [parse_reactant(r, 1) for r in right_cpds]
    return reactants

def build_s_matrix(df):
    s_matrix = dict()
    reactions = []
    for i in range(len(df)):
        rxn = df.id[i]
        reactions.append(rxn)
        try:
            reactants = parse_equation(df.equation[i])
        except:
            raise ValueError("can't parse equation {} - {}".format(i, df.equation[i]))
        
        for cpd in reactants:
            if cpd['cpd'] not in s_matrix:
                s_matrix[cpd['cpd']] = dict()
            if rxn not in s_matrix[cpd['cpd']]:
                s_matrix[cpd['cpd']][rxn] = dict()
            s_matrix[cpd['cpd']][rxn] = cpd['stoich']
    return s_matrix,reactions

def loop_for_steadycom(param):
    mu = 0.5
    X0 = 0.5 ## what is the value for x_o?

    # mu_bounds = {}
    # mu_bounds['LB'] = None
    # mu_bounds['UB'] = None
    LB = None
    UB = None
    lp_prob = simulate_steadycom(param,mu)
    # solve the model
    pulp_solver = pulp.solvers.GLPK_CMD(path=None, keepFiles=0, mip=1, msg=1, options=[])

    # solve for growth rate
    while LB is not None 
            and UB is not None 
            and abs(LB-UB) < 0.00001:

        lp_prob.solve(pulp_solver)
        obj_val = pulp.value(lp_prob.objective)

        X0 = 1 # total biomass (scale to 1 so that each X_k is relative abundance)
        if obj_val >= X0:
            # mu_bounds['LB'] = mu
            LB = mu
            mu = max(obj_val/X0,1.01)
        else:
            # mu_bounds['UB'] = mu
            UB = mu
            mu = max(obj_val/X0,0.99)

        # change the constraint including growth rate (mu)
        # lp_prob[""] = ...

    # solve FVA
    # change the objective function to each flux
    # lp_prob[""] = ...


def construct_steadycom(param,mu):

    model_inputs = param['model_inputs']

    # fetch information of S matrix for each organism k
    # k: index of organism
    # i: index of metabolites
    # j: index of reactions

    fba_client = fba_tools(self.callback_url)  # or however your callback url is set
                                               # then when you need the files

    S = {} # build S matrix for each FBA model k

    # define sets for steadycom
    reactions = {} # get reaction info for each FBA model k
    metabolites = {} # get metaboite info for each FBA model k
    organisms = []

    for model_input in model_inputs:
        model_upa = model_input['model_upa']
        files = fba_client.model_to_tsv_file({
            'workspace_name': workspace_name,  # from params
            'model_name': model_upa                     # also params
        })

        # files will have two "File" objects. you should be able to get to them as:
        # files['compounds_file']['path']
        # and
        # files['reactions_file']['path']
        
        # model_file = os.path.join(os.sep, "Users", "wjriehl", "Desktop", "iMR1_799.TSV", "iMR1_799-reactions.tsv")
        
        model_file = files['reactions_file']['path']
        model_df = pd.read_table(model_file)
        Sij,reactions = build_s_matrix(model_df)
        k = model_upa['id']
        organisms.append(k)
        S[k] = Sij
        metabolites[k] = Sij.keys()
        reactions[k] = reactions
        reactions_EX[k] = 
        reactions_biomass[k] = 

        # for model_input in model_inputs:
        #     GSM = model_input['model_upa']

        #     organism_id = GSM['id']
        #     modelreactions = GSM['modelreactions']
        #     modelcompounds = GSM['modelcompounds']

        #     S_matrix_ji = []
        #     for modelreaction in modelcompounds: # j
        #         S_matrix_ji[modelreaction] = {}
        #         reagents = modelreaction['modelReactionReagents']
        #         for reagent in reagents: # i
        #             met = reagent['modelcompound_ref']
        #             S_matrix_ji[modelreaction][met] = reagent['coefficient']

        #     S_matrix_ij = S_matrix_ji.transpose()


    #------- define variables
    X = pulp.LpVariable.dicts("X", organisms,
                              lowBound=0, upBound=1, cat='Continuous')

    for k in organisms:
        v = pulp.LpVariable.dicts("v", (k,reactions[k]),
                              lowBound=-M, upBound=M, cat='Continuous')

    #------- define LP problem
    lp_prob = pulp.LpProblem("SteadyCom", pulp.LpMaximize)

    #------- define objective function
        lp_prob += pulp.lpSum([X[k] for k in organisms])

    # define flux balance constraints
    for k in organisms:
        for i in S[k].keys():
            dot_S_v = pulp.lpSum([S[k][i][j] * v[k][j]
                                  for j in reactions[k]])
            condition = dot_S_v == 0
            lp_prob += condition#, label  

            for j in reactions[k]:
                lp_prob += v[k][j] <= UB[k][j] * X[k]
                lp_prob += v[k][j] >= LB[k][j] * X[k]

        lp_prob += v['bio1'][k] - X[k]*mu

    # constraints for medium (joshua: please add it here)
    

    return lp_prob