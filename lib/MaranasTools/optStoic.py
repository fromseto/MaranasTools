#/usr/bin/python

# author: Lin Wang
"""
 MinRxnFlux program 

"""
import pulp
import json

# GUROBI_CMD_OPTIONS = [('Threads', 2), ('TimeLimit', 1200), \
#                 ('MIPGapAbs', 1e-6), ('MIPGap', 1e-6), ('CliqueCuts', 2)]
# pulp_solver = pulp.solvers.GUROBI_CMD(path=None, keepFiles=0, mip=1, msg=1,
#                 options=GUROBI_CMD_OPTIONS)
pulp_solver = pulp.solvers.GLPK_CMD(path=None, keepFiles=0, mip=1, msg=1, options=[])

def simulate_optStoic(substrate_metabolite, target_metabolite):
    met_S = json.load(open("/kb/module/data/met_details_dict_v2.json"))
    # print met_S['C00001']
    metabolites = met_S.keys()
    # substrate_metabolite = "C00267"
    # target_metabolite = "C00033"
    dGmax = 5
    M = 100
    objective = 'MaxTargetYield'
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
    if objective == 'MaxTargetYield':
        lp_prob += s[target_metabolite], "MaxTargetYield"

    #------- define constraints
    # stoic(j)$(elem(j))..   sum(i,s(i)*m(i,j)) =e= 0;
    # bal(j)$(dG(j))..      sum(i,s(i)*m(i,j)) =l= dGmax*(-s('%substrate%'));

    elements = ['C','H','N','O','P','S','F','Cl','Mg','Fe','Se','Co','As','Br',
                'I','R','charge']    
    for j in elements:
        lp_prob += pulp.lpSum(s[i]*met_S[i][j] for i in metabolites) == 0

    lp_prob += pulp.lpSum(s[i]*met_S[i]['dGf'] for i in metabolites) <= \
                dGmax*s[substrate_metabolite],'element_'+ j

    lp_prob += s[substrate_metabolite] == -1

    if objective == 'MaxTargetYield':
        for i in metabolites:
            if met_S[i]['C'] > 0 and i != substrate_metabolite: 
                lp_prob += s[i] >= 0

    # lp_prob += pulp.lpSum(s[i] for i in metabolites) == 2
    # for i in metabolites:
    #     if i not in allow_list_metabolite:
    #         lp_prob += s[i] == 0

    #------- solve the problem
    lp_prob.solve(pulp_solver)

    # for i in metabolites:
    #     if s[i].varValue is not None:
    #         if s[i].varValue > EPS or s[i].varValue < -EPS:
    #             print i, s[i].varValue

    # LW: selct one of the solution from optStoic
    # this is my implementation as a way to selct stoichiometry from a number
    # of solutions. (different from the paper)
    # redefine the objective as minimize sum of all s except proton
    z = pulp.LpVariable.dicts("z", metabolites, lowBound=0, upBound=M, \
                                cat='Continuous')

    metabolites_no_proton = []
    for i in metabolites:
        if i != "C00080":
            metabolites_no_proton.append(i)

    lp_prob.setObjective(pulp.lpSum(z[i] for i in metabolites_no_proton))
    lp_prob.sense = pulp.LpMinimize
    for i in metabolites:
        lp_prob += s[i] <= z[i]
        lp_prob += -s[i] <= z[i]

    lp_prob += s[target_metabolite] == s[target_metabolite].varValue

    lp_prob.solve(pulp_solver)

    print "-----------------------"
    Ex_stoic = {}
    for i in metabolites:
        if s[i].varValue is not None:
            if s[i].varValue > EPS or s[i].varValue < -EPS:
                print i, s[i].varValue
                Ex_stoic[i] = s[i].varValue
    return Ex_stoic

if __name__ == '__main__':
    parameters = ''
    simulate_optStoic(None)