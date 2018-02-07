#/usr/bin/python

# author: Lin Wang, Joshua Chan
"""
 SteadyCom 

to do:
1. get LB and UB for each reactions in GSM
2. define the id of biomass reaction as input parameter
"""
import pulp
import json
from fba_tools.fba_toolsClient import fba_tools
import os
import pandas as pd
import re
from Workspace.WorkspaceClient import Workspace as workspaceService
from DataFileUtil.DataFileUtilClient import DataFileUtil

def parse_reactant(reactant, sign):
    """
    sign should be -1 or 1
    returns {'stoich': int, 'cpd': string, 'compartment': string}
    """
    # m = re.match('\((?P<stoich>.+)\)\s*(?P<cpd>[a-zA-Z0-9_]+)\[(?P<compartment>[a-zA-Z0-9_]+)\]', reactant)
    m = re.match('\((?P<stoich>\d*\.\d+|\d+)\)\s*(?P<cpd>[a-zA-Z0-9_]+)\[(?P<compartment>[a-zA-Z0-9_]+)\]', reactant)

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
    metabolites_EX = []
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

            if 'e' in cpd['compartment']: metabolites_EX.append(cpd['cpd'])
    return s_matrix,reactions,set(metabolites_EX)

def loop_for_steadycom(param,config,callback_url):
    mu = 0.1
    LB = None
    UB = None

    lp_prob,v,X,k,reactions_biomass = construct_steadycom(param,mu,config,callback_url)
    # solve the model
    pulp_solver = pulp.solvers.GLPK_CMD(path=None, keepFiles=0, mip=1, msg=1, options=[])

    # solve for growth rate
    while (LB == None) or (UB == None) or (abs(LB-UB) > 0.00001):
        # add constraints based on growth rate
        for k,bio_id in reactions_biomass.iteritems():
            lp_prob += v[k][bio_id] - X[k]*mu == 0                            
        
        lp_prob.solve(pulp_solver)
        obj_val = pulp.value(lp_prob.objective)

        X0 = 1 # total biomass (scale to 1 so that each X_k is relative abundance)
        if obj_val >= X0:
            # mu_bounds['LB'] = mu
            LB = mu
            mu = max(obj_val/X0,1.01)*mu
        else:
            # mu_bounds['UB'] = mu
            UB = mu
            mu = max(obj_val/X0,0.99)*mu

        print "--------------------------------"
        print "LB: ",LB
        print "UB: ",UB
        print "mu: ",mu
    return {'growth rate': mu}

def construct_steadycom(param,mu,config,callback_url):

    model_inputs = param['model_inputs']
    media = param['medium_upa']
    # print media
    # print param.keys()
    workspace_name = param['workspace_name']
    # print "--------------------"
    # print workspace_name
    # print "--------------------"
    ws = workspaceService(config['workspace-url'])
    # dataUtil = DataFileUtil(url=config['workspace-url'])
    # ws_id = dataUtil.ws_name_to_id(workspace_name)

    meida_object = ws.get_objects2({'objects': [{'name': media, 
                                                'workspace': workspace_name}]})
    
    print meida_object
    media_metabolites = meida_object['data'][0]['data']['mediacompounds']
    print media_metabolites
    media_dict = {}
    for met_info in media_metabolites:
        compound_ref = met_info["compound_ref"]
        compound_id = compound_ref.split('/')[-1]
        # print compound_id
        media_dict[compound_id] = met_info["maxFlux"]


    # fetch information of S matrix for each organism k
    # k: index of organism
    # i: index of metabolites
    # j: index of reactions

    fba_client = fba_tools(callback_url)  # or however your callback url is set
                                               # then when you need the files


    S = {} # build S matrix for each FBA model k

    # define sets for steadycom
    reactions = {} # get reaction info for each FBA model k
    reactions_EX = {}
    reactions_biomass = {}

    metabolites = {} # get metaboite info for each FBA model k
    metabolites_EX = {}
    organisms = []

    for model_input in model_inputs:
        model_upa_ref = model_input['model_upa']
        model_id = model_upa_ref.split('/')[1]
        model_upa = ws.get_objects2({'objects': [{'objid': model_id, 
                            'workspace': workspace_name}]})['data'][0]['data']
        files = fba_client.model_to_tsv_file({
            'workspace_name': workspace_name,  # from params
            'model_name': model_upa['name']                     # also params
        })

        # files will have two "File" objects. you should be able to get to them as:
        # files['compounds_file']['path']
        # and
        # files['reactions_file']['path']
        
        # model_file = os.path.join(os.sep, "Users", "wjriehl", "Desktop", "iMR1_799.TSV", "iMR1_799-reactions.tsv")
        
        model_file = files['reactions_file']['path']
        model_df = pd.read_table(model_file)

        Sij,rxns,mets_EX = build_s_matrix(model_df)

        k = model_upa['name']
        organisms.append(k)
        S[k] = Sij
        metabolites[k] = Sij.keys()
        reactions[k] = rxns
        reactions_biomass[k] = 'BIOMASS-Ecoli-core-w-GAM_c0'#model_upa['biomasses'][0].id
        metabolites_EX[k] = list(mets_EX)

    metabolites_com = []
    for k,met_list in metabolites_EX.iteritems():
        metabolites_com = list(set(metabolites_com + met_list))
    # metabolites_com = metabolites_com
    
    reactions_all = []
    for k,rxn_list in reactions.iteritems():
        reactions_all = list(set(reactions_all + rxn_list))
    # reactions_all = set(reactions_all)

    #------- define variables
    X = pulp.LpVariable.dicts("X", organisms,
                              lowBound=0, upBound=1, cat='Continuous')

    # for k in organisms:
    M = 1000
    print "***************"
    for r in reactions_all:
        if "BIOMASS" in r:
            print "found:", r
    print "***************"

    v = pulp.LpVariable.dicts("v", (organisms,reactions_all),
                          lowBound=-M, upBound=M, cat='Continuous')
    vex = pulp.LpVariable.dicts("vex", (organisms,metabolites_com),
                          lowBound=-M, upBound=M, cat='Continuous')
    e = pulp.LpVariable.dicts("e", metabolites_com,
                              lowBound=0, upBound=M, cat='Continuous')

    #------- define LP problem
    lp_prob = pulp.LpProblem("SteadyCom", pulp.LpMaximize)

    #------- define objective function
    lp_prob += pulp.lpSum([X[k] for k in organisms])

    biomass_labels = []
    # define flux balance constraints
    UB = 1000
    LB = -1000
    for k in organisms:
        for i in S[k].keys():
            dot_S_v = pulp.lpSum([S[k][i][j] * v[k][j]
                                  for j in S[k][i].keys()])
            if i in metabolites_EX[k]:
                condition = dot_S_v == vex[k][i]
            else:
                condition = dot_S_v == 0
            lp_prob += condition#, label  

            for j in reactions[k]:
                if j == 'ATPM_c0':
                    lp_prob += v[k][j] <= UB * X[k] #UB[k][j] * X[k]
                    lp_prob += v[k][j] >= 10 * X[k] #LB[k][j] * X[k]
                lp_prob += v[k][j] <= UB * X[k] #UB[k][j] * X[k]
                lp_prob += v[k][j] >= LB * X[k] #LB[k][j] * X[k]

    # constraints for medium (joshua: please add it here)
    for i in metabolites_com:
        # some exchange does not actually exist
        for k in organisms:
            if i not in metabolites_EX[k]:
                lp_prob += vex[k][i] == 0

        community_constraint = pulp.lpSum([vex[k][i] for k in organisms])
        if i in media_dict.keys():
            condition_comm = community_constraint - e[i] == -media_dict[i]
        else:
            condition_comm = community_constraint - e[i] == 0
        lp_prob += condition_comm
    return lp_prob,v,X,k,reactions_biomass

if __name__ == '__main__':
    loop_for_steadycom(None)