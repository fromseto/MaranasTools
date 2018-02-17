#/usr/bin/python

# author: Lin Wang
"""
 MinRxnFlux program 

"""
import pulp
import json
import re
import pandas as pd
import os
from fba_tools.fba_toolsClient import fba_tools
from Workspace.WorkspaceClient import Workspace as workspaceService


# GUROBI_CMD_OPTIONS = [('Threads', 2), ('TimeLimit', 1200), \
#                 ('MIPGapAbs', 1e-6), ('MIPGap', 1e-6), ('CliqueCuts', 2)]
# pulp_solver = pulp.solvers.GUROBI_CMD(path=None, keepFiles=0, mip=1, msg=1,
#                 options=GUROBI_CMD_OPTIONS)
pulp_solver = pulp.solvers.GLPK_CMD(path=None, keepFiles=0, mip=1, msg=1, options=[])

def construct_metS(list_of_mets,params,config):
    ws = workspaceService(config['workspace-url'])
    callback_url = os.environ['SDK_CALLBACK_URL']
    fba_client = fba_tools(callback_url, service_ver="beta")

    model_upa_ref = params['model_upa']
    workspace_name = params['workspace_name']

    print model_upa_ref
    model_name = model_upa_ref.split('/')[1]
    # model_upa = ws.get_objects2({'objects': [{'objid': model_id, 
    #                     'workspace': workspace_name}]})['data'][0]['data']
    # model_name = model_upa['name']

    # model_name = params['model_upa']
    print model_name

    use_fulldb = 0 # default use just a GSM
    if params['use_heterologous_steps'] == True:
        use_fulldb = 1

    model_files = fba_client.model_to_tsv_file({
    'model_name': model_name, #'iMR1_799',
    'workspace_name': workspace_name, #'lqw5322:narrative_1515706033382'
    'fulldb': use_fulldb
    })
    # pprint(model_files)

    mets_tsv = model_files['compounds_file']['path']
    model_df = pd.read_table(mets_tsv,index_col='id')

    met_S = {}


    for met_id in list_of_mets:
        met_info = {'C': 0,'H': 0,'N': 0,'O': 0,'P': 0,'S': 0,'F': 0,'Cl': 0,\
            'Mg': 0,'Fe': 0,'Se': 0,'Co': 0,'As': 0,'Br': 0,'I': 0,'R': 0,\
            'charge':None,'dGf':None}
        # get elements
        formula = model_df.loc[met_id]['formula']

        elements = re.findall(r'([A-Z][a-z]*)(\d*)', formula)
        # for element in elements:
        #     for ele,value in element.iteritems():
        #         met_info[ele] = float(value)
        elements = dict(elements)
        for key, value in elements.iteritems():
            if value == '': 
                met_info[key] = 1
            else:
                met_info[key] = float(value)

        # get charge and dG
        met_info['charge'] = model_df.loc[met_id]['charge']
        met_info['dGf'] = model_df.loc[met_id]['deltag']
        met_S[met_id] = met_info
    return met_S, model_files

def simulate_optStoic(params,config):
    # met_S = json.load(open("/kb/module/data/met_details_dict_v2.json"))

    # substrate_metabolite = params['start_compound']
    # target_metabolite = params['target_compound']

    # list_of_mets = [substrate_metabolite, target_metabolite]

    list_of_mets = []
    for data in params['reactant_stoichs']:
        list_of_mets.append(data['start_compound_id'])
    for data in params['product_stoichs']:
        list_of_mets.append(data['target_compound_id'])

    # add proton to the metabolit list because many time it is required
    if 'cpd00067_c0' not in list_of_mets:
        list_of_mets.append('cpd00067_c0')
    
    met_S,model_files = construct_metS(list_of_mets,params,config)
    # print met_S
    # print met_S['C00001']
    metabolites = met_S.keys()
    # substrate_metabolite = "C00267"
    # target_metabolite = "C00033"
    dGmax = 5
    M = 100
    # objective = 'MaxTargetYield'
    EPS = 1e-5
    # allow_list_metabolite = ['C00007',#    /*o2*/
    #                         'C00267',#    /*glc*/
    #                         # 'C00011',#    /*co2*/
    #                         'C00080',#    /*h+*/
    #                         'C00001',#    /*h2o*/
    #                         # 'C02457',#    /*13pdo*/
    #                         'C00033',#    /*acetate*/
    #                         ]

    #------- define variables
    # zp      objective function
    # s(i)    stoichiometric coefficient of metabolite i
    # ;
    s = pulp.LpVariable.dicts("s", metabolites, lowBound=-M, upBound=M, \
                                cat='Continuous')


    #------- define LP problem
    lp_prob = pulp.LpProblem("OptStoic", pulp.LpMaximize)
    

    #------- define objective function
    # obj..                     zp =e= sum(i$pdt(i),s(i))/(-s('%substrate%'));
    objective_stoic = params['objective']
    lp_prob += s[objective_stoic], "MaxTargetYield"

    # if objective == 'MaxTargetYield':
        # lp_prob += s[target_metabolite], "MaxTargetYield"

    #------- define constraints
    # stoic(j)$(elem(j))..   sum(i,s(i)*m(i,j)) =e= 0;
    # bal(j)$(dG(j))..      sum(i,s(i)*m(i,j)) =l= dGmax*(-s('%substrate%'));

    elements = ['C','H','N','O','P','S','F','Cl','Mg','Fe','Se','Co','As','Br',
                'I','R','charge']
    for j in elements:
        lp_prob += pulp.lpSum(s[i]*met_S[i][j] for i in metabolites) == 0

    lp_prob += pulp.lpSum(s[i]*met_S[i]['dGf'] for i in metabolites) <= \
                dGmax,'element_'+ j
                # dGmax*s[substrate_metabolite],'element_'+ j

    for data in params['reactant_stoichs']:
        if data['fixed_stoich'] != None:
            stoic = -data['fixed_stoich']
            lp_prob += s[data['start_compound_id']] == stoic
    for data in params['product_stoichs']:
        if data['fixed_stoich'] != None:
            stoic = -data['fixed_stoich2']
            lp_prob += s[data['target_compound_id']] == stoic

    # lp_prob += s[substrate_metabolite] == -1

    # if objective == 'MaxTargetYield':
    #     for i in metabolites:
    #         if met_S[i]['C'] > 0 and i != substrate_metabolite: 
    #             lp_prob += s[i] >= 0

    # lp_prob += pulp.lpSum(s[i] for i in metabolites) == 2
    # for i in metabolites:
    #     if i not in allow_list_metabolite:
    #         lp_prob += s[i] == 0

    #------- solve the problem
    lp_prob.solve(pulp_solver)

    Ex_stoic = {}
    for i in metabolites:
        if s[i].varValue is not None:
            if s[i].varValue > EPS or s[i].varValue < -EPS:
                print i, s[i].varValue
                Ex_stoic[i] = s[i].varValue

    # # LW: selct one of the solution from optStoic
    # # this is my implementation as a way to selct stoichiometry from a number
    # # of solutions. (different from the paper)
    # # redefine the objective as minimize sum of all s except proton
    # z = pulp.LpVariable.dicts("z", metabolites, lowBound=0, upBound=M, \
    #                             cat='Continuous')

    # metabolites_no_proton = []
    # for i in metabolites:
    #     if i != "C00080":
    #         metabolites_no_proton.append(i)

    # lp_prob.setObjective(pulp.lpSum(z[i] for i in metabolites_no_proton))
    # lp_prob.sense = pulp.LpMinimize
    # for i in metabolites:
    #     lp_prob += s[i] <= z[i]
    #     lp_prob += -s[i] <= z[i]

    # lp_prob += s[target_metabolite] == s[target_metabolite].varValue

    # lp_prob.solve(pulp_solver)

    # print "-----------------------"
    # Ex_stoic = {}
    # for i in metabolites:
    #     if s[i].varValue is not None:
    #         if s[i].varValue > EPS or s[i].varValue < -EPS:
    #             print i, s[i].varValue
    #             Ex_stoic[i] = s[i].varValue
    return Ex_stoic,model_files

if __name__ == '__main__':
    parameters = ''
    simulate_optStoic(None)