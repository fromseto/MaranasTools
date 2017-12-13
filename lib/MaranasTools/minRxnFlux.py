#/usr/bin/python

# author: Chiam Yu Ng
"""
 MinRxnFlux program 
Currently, it has been tested with GLPK, Gurobi and CPLEX solver.

"""

import pulp
import os, time, sys
import copy
import random
import string  # to generate random hex code
import gams_parser
import json
import pdb
import logging
from ..core import database


# Global variables/solver options
EPS = 1e-5
GUROBI_OPTIONS = 'Threads=2 TimeLimit=1800 MIPGapAbs=1e-6 MIPGap=1e-6 CliqueCuts=2'

current_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.normpath(os.path.join(
    current_dir, '../data/', 'optstoic_db_v3'))

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