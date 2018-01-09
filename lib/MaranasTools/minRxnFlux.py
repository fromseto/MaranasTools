#/usr/bin/python

# author: Chiam Yu Ng
"""
 MinRxnFlux program 

"""

import pulp
import os, time, sys
import copy
import random
import string  # to generate random hex code
import json
import pdb
import logging
import pandas as pd
import gams_parser

# Global variables/solver options
EPS = 1e-5
GUROBI_OPTIONS = 'Threads=2 TimeLimit=1800 MIPGapAbs=1e-6 MIPGap=1e-6 CliqueCuts=2'

current_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.normpath(os.path.join(
    current_dir, '../../data/'))
res_dir = os.path.normpath(os.path.join(
    current_dir, '../../data/','out'))

# "linux2" is the platform for lxcluster/hammer
# if sys.platform == 'cygwin' and sys.platform =='linux2':

def load_pulp_solver():
    GUROBI_CMD_OPTIONS = [('Threads', 2), ('TimeLimit', 1800),
                          ('MIPGapAbs', 1e-6), ('MIPGap', 1e-6), ('CliqueCuts', 2)]
    CPLEX_CMD_OPTIONS = ['mip tolerances mipgap 1e-6','mip tolerances absmipgap 1e-6']
    GLPK_CMD_OPTIONS  = ['--clique', '--pcost', '--gomory', '--mipgap', '1e-6']
    #Load solvers in the order of preferences
    pulp_solvers = {'solvers':[],
                    'solver_name':['GUROBI', 'GUROBI_CMD', 'CPLEX_CMD', 'GLPK_CMD']}

    pulp_solvers['solvers'].append(
        pulp.solvers.GUROBI(mip=True, msg=True,
                            timeLimit=1800, MIPGapAbs=1e-6))
    pulp_solvers['solvers'].append(
        pulp.solvers.GUROBI_CMD(path=None, keepFiles=0, mip=1, msg=1,
                                options=GUROBI_CMD_OPTIONS))
    pulp_solvers['solvers'].append(
        pulp.solvers.CPLEX_CMD(path=None, keepFiles=0, mip=1, msg=1,
                               options=CPLEX_CMD_OPTIONS,
                               timelimit=1800))
    pulp_solvers['solvers'].append(
        pulp.solvers.GLPK_CMD(msg=1, mip=1, options=GLPK_CMD_OPTIONS))

    for i, p in enumerate(pulp_solvers['solvers']):
        if p.available():
            pulp_solver = p
            logging.warning("Pulp solver set to %s."%pulp_solvers['solver_name'][i])

            if hasattr(pulp_solver,'tmpDir'):
                pulp_solver.tmpDir = './'

            if 'solver_name' == 'GLPK_CMD':
                logging.warning("GLPK takes a significantly longer time to solve "
                                "OptStoic. Please be patient.")
            break
    else:
        logging.warning("No solver is available! Program will be terminated.")
        return None
    return pulp_solver


class Database(object):
    """optstoic Database class: loading database from GAMS input file"""

    def __init__(self, version='', data_filepath='../data/', dbdict={}):
        self.data_filepath = data_filepath
        self.dbdict = dbdict
        self.reactions = []
        self.metabolites = []
        self.S = {}
        self.Sji = {}
        self.rxntype = []
        self.loops = []
        self.Ninternal = {}
        self.all_excluded_reactions = []
        self.excluded_reactions = dbdict.get('excluded_reactions_list') or []
        self.user_defined_export_rxns = []
        self.blocked_rxns = []

    def load(self):

        # Load S matrix dictionary (S(i,j)) from json
        if 'S' not in self.dbdict:
            self.Sji = json.load(open(os.path.join(self.data_filepath,
                                                   self.dbdict['Sji']), 'r+'))

            self.S = self.transpose_S(self.Sji)
        else:
            try:
                self.S = json.load(open(os.path.join(self.data_filepath,
                                                     self.dbdict['S']), 'r+'))

            except:
                # TODO add function to differentiate json/txt input
                self.S = gams_parser.convert_parameter_table_to_dict(
                    os.path.join(self.data_filepath,
                                 '20160616_optstoic_Sij.txt')
                )
            self.Sji = self.transpose_S(self.S)

        # Load reactions
        logging.debug('Reading reaction file...')
        self.reactions = gams_parser.convert_set_to_list(
            os.path.join(self.data_filepath, self.dbdict['reaction'])
        )

        self.internal_rxns = copy.deepcopy(self.reactions)

        logging.debug('Reading metabolite file...')
        self.metabolites = gams_parser.convert_set_to_list(
            os.path.join(self.data_filepath, self.dbdict['metabolite'])
        )

        logging.debug('Reading blocked reactions file...')
        if 'blocked_rxns' in self.dbdict:
            self.blocked_rxns = gams_parser.convert_set_to_list(
                os.path.join(self.data_filepath, self.dbdict['blocked_rxns'])
            )

        self.all_excluded_reactions = list(
            set(self.excluded_reactions + self.blocked_rxns)
        )

        logging.debug('Reading reaction type file...')
        self.rxntype = gams_parser.convert_parameter_list_to_dict(
            os.path.join(self.data_filepath, self.dbdict['reactiontype']),
            datadict=None
        )

    @staticmethod
    def transpose_S(Sji):
        """Tranpose Sji into Sij and also Sij to Sji dictionary."""
        # Update to pandas 0.19 (using sparse dataframe)
        df_Sji = pd.DataFrame(Sji).T
        Sij = dict(
            (k, v.dropna().to_dict()) for k, v in pd.compat.iteritems(df_Sji)
        )
        return Sij

    @staticmethod
    def to_json(Sdict, filepath):
        with open(filepath, 'w+') as fp:
            json.dump(Sdict, fp, sort_keys=True, indent=4)

    def to_mat_file():
        """write to matlab file"""
        pass

    def get_reaction_type(self, rid, verbose=True):
        try:
            self.rxntype[rid]
        except:
            print "Reaction %s not in database!" % rid
            return None
        else:
            if verbose:
                print "Reaction: {0} is ({1}) {2}".format(
                    rid, self.rxntype[rid], REACTION_TYPE[self.rxntype[rid]]
                )
            return self.rxntype[rid]

    def extend_S_from_file(self, filename='Sij_extension_for_glycolysis.txt'):
        self.S = gams_parser.convert_parameter_table_to_dict(
            os.path.join(self.data_filepath, filename), Sdict=self.S)

    def update_S(self, extension_dict):
        temp_rxn = []
        for met, entries in extension_dict.iteritems():
            if met not in self.S:
                self.S[met] = {}
                self.metabolites.append(met)
            for rxn, coeff in entries.iteritems():
                self.S[met][rxn] = float(coeff)
                if rxn not in self.reactions:
                    self.reactions.append(rxn)
                    temp_rxn.append(rxn)
        return self.S, temp_rxn

    def set_database_export_reaction(self, export_reactions_Sij_dict):
        _, temp_rxn = self.update_S(export_reactions_Sij_dict)
        if len(self.user_defined_export_rxns) != 0:
            logging.warning("Warning: The current list of export reactions\
                will be replaced! %s" % str(self.user_defined_export_rxns))
        self.user_defined_export_rxns = temp_rxn
        for rxn in self.user_defined_export_rxns:
            self.rxntype[rxn] = 4
        return 1

    def update_rxntype(self, new_reaction_type_dict):
        for (r, rtype) in new_reaction_type_dict.iteritems():
            try:
                t0 = self.rxntype[r]
                self.rxntype[r] = rtype
            except KeyError:
                logging.warning('Reaction %s not in database!' % r)
            else:
                logging.info('Reaction %s has been updated from %s to %s.'
                             % (r, REACTION_TYPE[t0], REACTION_TYPE[rtype])
                             )

        return self.rxntype

    def __str__(self):
        return "OptStoic Database(Version='%s')" % self.version



class MinRxnFlux(object):
    """An MinRxnFlux problem Class"""
    def __init__(self, database, objective='MinFlux',
                 specific_bounds={},
                 max_iteration=2,
                 pulp_solver=None,
                 data_filepath=data_dir,
                 result_filepath='result/',
                 M=1000):
        """
        Keyword Arguments:
        database -- an optStoic Database object (equivalent to GSM model)
        objective -- change the objective to MinFlux or MinRxn
        specific_bounds -- LB and UB for exchange reactions. E.g. {'Ex_glc': {'LB': -1, 'UB':-1}}
        pulp_solver -- a pulp.solvers object (load any the user-defined solver)
        max_iteration -- set the default maximum number of iteration
        data_filepath -- filepath for data
        M -- maximum flux bound (default 1000)
        """
        self.objective = objective
        if len(specific_bounds) == 0 or specific_bounds is None:
            raise Exception("specific_bounds must be specified!")
        self.specific_bounds = specific_bounds
        self.M = M
        self.varCat = 'Integer'
        self.max_iteration = max_iteration
        self.data_filepath = data_filepath
        self.result_filepath = result_filepath
        self.database = database
        self.pathways = {}
        self.iteration = 1
        self.lp_prob = None
        self.pulp_solver = pulp_solver
        self.lp_prob_fname = "OptStoic_{0}".format(self.generate_random_string(6))

    @staticmethod
    def generate_random_string(N):
        """LP file is appended with random string when using command line mode
            to prevent overwrite/read issues.
        """
        return ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(N))

    def change_objective(self, new_objective):
        if new_objective not in ['MinFlux', 'MinRxn']:
            raise ValueError("The objective for OptStoic is "
                             "not correctly defined. "
                             "Please use either 'MinFlux' or 'MinRxn'.")
        self.objective = new_objective

    def create_minflux_problem(self):
        """
        Create minflux/minRxn LP problem (a pulp LpProblem object)
        """
        logging.info("Formulating problem...")

        M = self.M

        # Initialize variables
        v = pulp.LpVariable.dicts("v", self.database.reactions,
                                  lowBound=-M, upBound=M, cat=self.varCat)
        vf = pulp.LpVariable.dicts("vf", self.database.reactions,
                                   lowBound=0, upBound=M, cat=self.varCat)
        vb = pulp.LpVariable.dicts("vb", self.database.reactions,
                                   lowBound=0, upBound=M, cat=self.varCat)
        yf = pulp.LpVariable.dicts("yf", self.database.reactions,
                                   lowBound=0, upBound=1, cat='Binary')
        yb = pulp.LpVariable.dicts("yb", self.database.reactions,
                                   lowBound=0, upBound=1, cat='Binary')
        a = pulp.LpVariable.dicts("a", self.database.reactions,
                                  lowBound=0, upBound=1, cat='Binary')
        G = pulp.LpVariable.dicts("G", self.database.reactions,
                                  lowBound=-M, upBound=M, cat='Continuous')

        for j in self.database.reactions:

            if j in self.database.all_excluded_reactions:
                v[j].lowBound = 0
                v[j].upBound = 0
                vf[j].lowBound = 0
                vf[j].upBound = 0
                yf[j].lowBound = 0
                yf[j].upBound = 0
                yb[j].lowBound = 0
                yb[j].upBound = 0

            if self.database.rxntype[j] == 0:
                # Forward irreversible
                v[j].lowBound = 0
                v[j].upBound = M
                yb[j].upBound = 0
                vb[j].upBound = 0

            elif self.database.rxntype[j] == 1:
                # Reversible
                v[j].lowBound = -M
                v[j].upBound = M

            elif self.database.rxntype[j] == 2:
                # Reverse irreversible
                v[j].lowBound = -M
                v[j].upBound = 0
                vf[j].upBound = 0
                yf[j].upBound = 0

            elif self.database.rxntype[j] == 4:
                v[j].lowBound = 0
                v[j].upBound = 0

        # Fix stoichiometry of source/sink metabolites
        # Allow user to change this for generalization
        for rxn, bounds in self.specific_bounds.iteritems():
            v[rxn].lowBound = bounds['LB']
            v[rxn].upBound = bounds['UB']

        LB = {}
        UB = {}
        for j in self.database.reactions:
            LB[j] = v[j].lowBound
            UB[j] = v[j].upBound

        lp_prob = pulp.LpProblem("OptStoic", pulp.LpMinimize)

        
        # Min-Rxn objective
        if self.objective == 'MinRxn':
            condition = pulp.lpSum([yf[j] + yb[j]
                                    for j in self.database.reactions
                                    if self.database.rxntype[j] != 4])
            lp_prob += condition, "MinRxn"

        # Min-Flux objective
        elif self.objective == 'MinFlux':
            condition = pulp.lpSum([vf[j] + vb[j]
                                    for j in self.database.reactions
                                    if self.database.rxntype[j] != 4])
            lp_prob += condition, "MinFlux"

        # Constraints
        # Mass_balance
        for i in self.database.metabolites:
            # If metabolites not involve in any reactions
            if i not in self.database.S:
                continue
            label = "mass_balance_%s" % i
            dot_S_v = pulp.lpSum([self.database.S[i][j] * v[j]
                                  for j in self.database.S[i].keys()])
            condition = dot_S_v == 0
            lp_prob += condition, label

        if self.objective == 'MinRxn':
            for j in self.database.reactions:
                lp_prob += v[j] >= y[j]*LB[j], "cons1_%s"%j
                lp_prob += v[j] <= y[j]*UB[j], "cons2_%s"%j

        if self.objective == 'MinFlux':
            for j in self.database.reactions:
                lp_prob += (v[j] == vf[j] - vb[j]), "flux_%s" % j

                # These constraints ensure that when yf=0 and yb=0 ,
                # no flux goes through the reaction
                lp_prob += vf[j] >= yf[j] * 0.5, "cons1_%s" % j
                lp_prob += vf[j] <= yf[j] * M, "cons2_%s" % j
                lp_prob += vb[j] >= yb[j] * 0.5, "cons3_%s" % j
                lp_prob += vb[j] <= yb[j] * M, "cons4_%s" % j
                # Ensure that either yf or yb can be 1, not both
                lp_prob += yf[j] + yb[j] <= 1, 'cons5_%s' % j
        else:
            a = None
            G = None

        return lp_prob, v, vf, vb, yf, yb, a, G

    def solve(self, exclude_existing_solution=False, outputfile="OptStoic_pulp_result.txt", max_iteration=None):
        """
        Solve OptStoic problem using pulp.solvers interface

        Keyword Arguments:
            outputfile: name of outpufile
            max_iteration: Externally specified maximum number of pathway to be found using OpStoic.
                            If not specified, it will set to the internal max iterations.
        """
        if self.objective not in ['MinFlux', 'MinRxn']:
            raise ValueError("The objective for OptStoic is not correctly defined.Please use either 'MinFlux' or 'MinRxn'.")

        if max_iteration is None:
            max_iteration = self.max_iteration

        logging.info("Finding multiple pathways using Optstoic %s...", self.objective)
        lp_prob, v, vf, vb, yf, yb, a, G = self.create_minflux_problem()

        # Create integer cut for existing pathways
        if exclude_existing_solution and bool(self.pathways):
            self.iteration = max(self.pathways.keys()) + 1
            if self.iteration > max_iteration:
                raise ValueError('Max iteration is less than current '
                                 'iteration. Increase max_iteration '
                                 'before solving!')

            for ind, entry in self.pathways.iteritems():
                rxnlist = list(set(entry['reaction_id']) -
                               set(self.database.user_defined_export_rxns))
                condition = pulp.lpSum(
                    [(1 - yf[j] - yb[j]) for j in rxnlist]) >= 1
                lp_prob += condition, "IntegerCut_%d" % ind


        logging.info("Solving problem...")
        if self.iteration == 1:
            result_output = open(os.path.join(self.result_filepath, outputfile),"w+")
        else:
            result_output = open(os.path.join(self.result_filepath, outputfile),"a+")

        while True and self.iteration <= max_iteration:
            logging.info("Iteration %s", self.iteration)
            #lp_prob.writeLP("OptStoic.lp", mip=1)  #optional
            e1 = time.time()
            lp_prob.solve(self.pulp_solver)
            e2 = time.time()
            logging.info("This iteration solved in %.3f seconds.", (e2-e1))

            # The solution is printed if it was deemed "optimal
            if pulp.LpStatus[lp_prob.status] == "Optimal":
                logging.info("Writing result to output file...")
                result_output.write("\nIteration no.: %d\n" %self.iteration)
                result_output.write("\nModelstat: %s\n" %pulp.LpStatus[lp_prob.status])

                res = {}
                res['reaction_id'] = []
                res['flux'] = []
                res['iteration'] = self.iteration
                res['time'] = (e2-e1)
                res['modelstat'] = "Optimal"

                for j in self.database.reactions:
                    if v[j].varValue is not None:
                        if v[j].varValue > EPS or v[j].varValue < -EPS:
                            res['reaction_id'].append(j)
                            res['flux'].append(v[j].varValue)
                            result_output.write("%s %.8f\n" %(v[j].name, v[j].varValue))

                result_output.write("%s = %.8f\n" % (self.objective, pulp.value(lp_prob.objective)))
                result_output.write("----------------------------------\n\n")

                integer_cut_reactions = list(set(res['reaction_id']) - set(self.database.user_defined_export_rxns))

                self.pathways[self.iteration] = res
                json.dump(self.pathways, open(os.path.join(self.result_filepath,'temp_pathways.json'),'w+'), sort_keys=True, indent=4)

                # Integer cut constraint is added so that
                # the same solution cannot be returned again
                condition = pulp.lpSum([(1 - yf[j] - yb[j])
                                        for j in integer_cut_reactions]) >= 1
                lp_prob += condition, "IntegerCut_%d" % self.iteration
                self.iteration += 1

            # If a new optimal solution cannot be found, end the program
            else:
                break
        result_output.close()

        self.lp_prob = lp_prob

        return self.lp_prob, self.pathways

    def add_existing_pathways(self, user_defined_pathways):
        """
        Add list of existing solutions (pathways) to be
        excluded from being identified.

        Keyword Arguments:
        user_defined_pathways -- pathways output from solve_gurobi_cl()
                                 or solve() or pathway in dictionary format
                                 e.g. {1: 'reaction_id': []}
        """
        if (isinstance(user_defined_pathways, dict) and
           ('reaction_id' in user_defined_pathways.values()[0])):
            self.pathways = copy.deepcopy(user_defined_pathways)
        else:
            raise ValueError("user_defined_pathways must be a "
                             "pathways dictionary "
                             "{1: 'reaction_id': ['R00001', 'R00002']}")

    def reset_pathways(self):
        """
        Reset self.pathways to empty dictionary
        """
        self.pathways = {}

    def __repr__(self):
        return "<OptStoic(objective='%s')>" % (self.objective)

def run_minRxnFlux(optSotic_result_dict):
    db = load_db_rxns(data_dir,optSotic_result_dict)
    logging.info("Database loaded!")

    # df_bounds = pd.read_csv('specific_bounds.csv', index_col=0)
    # specific_bounds = df_bounds.to_dict(orient='index')
    specific_bounds = {}
    for met,stoich in optSotic_result_dict.iteritems():
        if met == "C00080":
            specific_bounds['EX_'+met] = {'LB': -5, 'UB': 5}
        else:
            specific_bounds['EX_'+met] = {'LB': stoich, 'UB': stoich}

    pulp_solver = load_pulp_solver()

    MAX_ITERATION = 10
    USE_LOOPLESS = False
    CLEANUP = True #set to true if you want to delete all .sol and .lp files

    test = MinRxnFlux(db, objective='MinFlux',
                    specific_bounds=specific_bounds,
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

def load_db_rxns(data_dir,optSotic_result_dict):

    dbdict = {
        'Sji': 'optstoic_v3_Sji_dict.json',
        'reaction': 'optstoic_v3_reactions.txt',
        'metabolite': 'optstoic_v3_metabolites.txt',
        'reactiontype': 'optstoic_v3_reactiontype.txt',
    }

    DB = Database(data_filepath=data_dir, dbdict=dbdict)
    
    DB.load()

    # user_defined_export_rxns_Sji = {
    # 'EX_glucose': {'C00031': -1.0}, #use C00031 for more general glucose
    # # 'EX_D_Xylose': {'C00181': -1.0},
    # 'EX_hplus': {'C00080': -1.0}, #pulp or gurobi has issue with "h+"
    # # 'EX_ac': {'C00033': -1.0},
    # 'EX_h2o': {'C00001': -1.0}
    # }
    user_defined_export_rxns_Sji = {}
    for met,stoich in optSotic_result_dict.iteritems():
        user_defined_export_rxns_Sji['EX_'+met] = {met: stoich}


    user_defined_export_rxns_Sij = Database.transpose_S(
        user_defined_export_rxns_Sji
    )

    DB.set_database_export_reaction(user_defined_export_rxns_Sij)

    return DB

if __name__ == '__main__':
    parameters = {'C00033': 3.0,
                'C00267': -1.0,
                'C00080': 3.0}
    run_minRxnFlux(parameters)