#/usr/bin/python

# author: Lin Wang
"""
 MinRxnFlux program 

"""
import pulp
import json

GUROBI_CMD_OPTIONS = [('Threads', 2), ('TimeLimit', 1200), \
                ('MIPGapAbs', 1e-6), ('MIPGap', 1e-6), ('CliqueCuts', 2)]
pulp_solver = pulp.solvers.GUROBI_CMD(path=None, keepFiles=0, mip=1, msg=1,
                options=GUROBI_CMD_OPTIONS)

def run_optStoic(parameters):
    met_S = json.load(open("../../data/met_details_dict_v2.json"))
    # print met_S['C00001']
    metabolites = met_S.keys()
    substrate_metabolite = "C00267"
    target_metabolite = "C00033"
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
    for i in metabolites:
        if s[i].varValue is not None:
            if s[i].varValue > EPS or s[i].varValue < -EPS:
                print i, s[i].varValue


def run_minRxnFlux(parameters):
    db = load_db_v3_phosphoketolase(data_dir)
    logging.info("Database loaded!")

    df_bounds = pd.read_csv('specific_bounds.csv', index_col=0)
    specific_bounds = df_bounds.to_dict(orient='index')

    pulp_solver = load_pulp_solver()

    MAX_ITERATION = 10
    USE_LOOPLESS = False
    CLEANUP = True #set to true if you want to delete all .sol and .lp files

    test = MinRxnFlux(db, objective='MinFlux',
                    specific_bounds=specific_bounds,
                    use_loopless=USE_LOOPLESS,
                    max_iteration=MAX_ITERATION,
                    pulp_solver=pulp_solver,
                    data_filepath=data_dir,
                    result_filepath=res_dir,
                    M=1000)

    lp_prob, pathways = test.solve(outputfile='test_optstoic_mac.txt')

    # Creating kegg model and drawing pathways
    f = open(os.path.join(res_dir, 'test_KeggModel.txt'), 'w+')
    pathway_objects = []

    for ind, res in sorted(test.pathways.iteritems()):
        p = Pathway(id=ind, name='OptStoic', reaction_ids=res['reaction_id'], fluxes=res['flux'])
        p.rearrange_reaction_order()
        pathway_objects.append(p)
        generate_kegg_model(p, filehandle=f)
        graph_title = "{0}_P{1}".format(p.name, p.id)
        draw_pathway(p, imageFileName=os.path.join(res_dir+'/pathway_{0:03d}'.format(p.id)),
                    imageFormat='png', graphTitle=graph_title, darkBackgroundMode=False, debug=True)
    f.close()
    print "Generate kegg_model and draw pathway: Pass!"



def load_db_v3_phosphoketolase(data_dir)

    dbdict = {
        'Sji': 'optstoic_v3_Sji_dict.json',
        'reaction': 'optstoic_v3_reactions.txt',
        'metabolite': 'optstoic_v3_metabolites.txt',
        'reactiontype': 'optstoic_v3_reactiontype.txt',
    }

    DB = Database(data_filepath=data_dir, dbdict=dbdict)
    
    DB.load()

    user_defined_export_rxns_Sji = {
    'EX_glucose': {'C00031': -1.0}, #use C00031 for more general glucose
    # 'EX_D_Xylose': {'C00181': -1.0},
    'EX_hplus': {'C00080': -1.0}, #pulp or gurobi has issue with "h+"
    # 'EX_ac': {'C00033': -1.0},
    'EX_h2o': {'C00001': -1.0}
    }

    user_defined_export_rxns_Sij = Database.transpose_S(
        user_defined_export_rxns_Sji
    )

    DB.set_database_export_reaction(user_defined_export_rxns_Sij)

    return DB

if __name__ == '__main__':
    parameters = ''
    run_optStoic(None)