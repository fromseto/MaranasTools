#/usr/bin/python

# author: Lin Wang, Joshua Chan
"""
 SteadyCom 

"""
import pulp
import json

pulp_solver = pulp.solvers.GLPK_CMD(path=None, keepFiles=0, mip=1, msg=1, options=[])

def loop_for_steadycom(param):
    mu = 0.5
    X0 = 0.5 ## what is the value for x_o?

    mu_bounds = {}
    mu_bounds['LB'] = None
    mu_bounds['UB'] = None
    while mu_bounds['LB'] is not None and mu_bounds['UB'] is not None:
        obj_val = simulate_steadycom(param,mu)
        if obj_val >= X0:
            mu_bounds['LB'] = mu
            mu = max(obj_val/X0,1.01)
        else:
            mu_bounds['UB'] = mu
            mu = max(obj_val/X0,0.09)

    # root finding algorithm? what does it mean
    # why not just do a bisection search to find mu


def simulate_steadycom(param,mu):

    model_inputs = param['model_inputs']

    list_of_S_matrix = []
    for model_input in model_inputs:
        GSM = model_input['model_upa']

        organism_id = GSM['id']
        modelcompounds = GSM['modelreactions']
        modelcompounds = GSM['modelcompounds']

        S_matrix_ji = []
        for modelreaction in modelcompounds: # j
            S_matrix_ji[modelreaction] = {}
            reagents = modelreaction['modelReactionReagents']
            for reagent in reagents: # i
                met = reagent['modelcompound_ref']
                S_matrix_ji[modelreaction][met] = reagent['coefficient']

        S_matrix = S_matrix_ji.transpose()
        list_of_S_matrix.append(S_matrix)


    #------- define variables
    X = pulp.LpVariable.dicts("v", organisms,
                              lowBound=0, upBound=1, cat='Continuous')
    v = pulp.LpVariable.dicts("v", (reactions,organisms),
                              lowBound=-M, upBound=M, cat='Continuous')
    # mu = pulp.LpVariable("mu", lowBound=0, upBound=1, cat='Continuous')

    #------- define LP problem
    lp_prob = pulp.LpProblem("SteadyCom", pulp.LpMaximize)

    #------- define objective function
        lp_prob += pulp.lpSum([X[o] for o in organisms])

    # define flux balance constraints
    for o in organisms:
        for i in S_matrix.keys():
            dot_S_v = pulp.lpSum([S_matrix[i][j] * v[j][o]
                                  for j in S_matrix_ij[i].keys()])
            condition = dot_S_v == 0
            lp_prob += condition#, label  

            for j in S_matrix_ij[i].keys():
                lp_prob += v[j][o] <= UB[j] * X[o]
                lp_prob += v[j][o] >= LB[j] * X[o]

            lp_prob += v['bio1'][o] - X[o]*mu

    # constraints for medium (joshua: please add it here)
    

    # solve the model
    lp_prob.solve(pulp_solver)
    objective_val = pulp.value(lp_prob.objective)
    return objective_val

if __name__ == '__main__':
    parameters = ''
    simulate_optStoic(None)