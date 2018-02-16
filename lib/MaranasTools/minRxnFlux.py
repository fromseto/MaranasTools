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
# import sys
# sys.path.append(os.getcwd())
import gams_parser
from drawpath import *
# from CreateReport import CreateReport
from pulp_scip import SCIP_CMD
import re
from fba_tools.fba_toolsClient import fba_tools

# Global variables/solver options
EPS = 1e-5
GUROBI_OPTIONS = 'Threads=2 TimeLimit=1800 MIPGapAbs=1e-6 MIPGap=1e-6 CliqueCuts=2'

current_dir = os.path.dirname(os.path.abspath(__file__))
# data_dir = os.path.normpath(os.path.join(
#     current_dir, '../../data/'))
# data_dir = os.path.normpath(os.path.join(
#     current_dir, '../../data/modelSEED'))
# res_dir = os.path.normpath(os.path.join(
#     current_dir, '../../data/modelSEED/','out'))

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

        # # Load S matrix dictionary (S(i,j)) from json
        # if 'S' not in self.dbdict:
        #     self.Sji = json.load(open(os.path.join(self.data_filepath,
        #                                            self.dbdict['Sji']), 'r+'))
        #     self.S = self.transpose_S(self.Sji)
        # else:
        #     try:
        #         self.S = json.load(open(os.path.join(self.data_filepath,
        #                                              self.dbdict['S']), 'r+'))

        #     except:
        #         # TODO add function to differentiate json/txt input
        #         self.S = gams_parser.convert_parameter_table_to_dict(
        #             os.path.join(self.data_filepath,
        #                          '20160616_optstoic_Sij.txt')
        #         )
        #     self.Sji = self.transpose_S(self.S)
        self.S = json.load(open(os.path.join(self.data_filepath,
                                                   self.dbdict['Sij']), 'r+'))
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
            # pdb.set_trace()
            if met not in self.S.keys():
                self.S[met] = {}
                self.metabolites.append(met)
            for rxn, coeff in entries.iteritems():
                self.S[met][rxn] = float(coeff)
                if rxn not in self.reactions:
                    self.reactions.append(rxn)
                    temp_rxn.append(rxn)
        return self.S, temp_rxn

    def set_database_export_reaction(self, export_reactions_Sij_dict):
        # pdb.set_trace()
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
                 data_filepath=None,
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

    def create_minflux_problem(self,params,model_files):
        """
        Create minflux/minRxn LP problem (a pulp LpProblem object)
        """
        logging.info("Formulating problem...")

        M = self.M

        # Initialize variables
        v = pulp.LpVariable.dicts("v", self.database.reactions,
                                  lowBound=-M, upBound=M, cat='Continuous')
        vf = pulp.LpVariable.dicts("vf", self.database.reactions,
                                   lowBound=0, upBound=M, cat='Continuous')
        vb = pulp.LpVariable.dicts("vb", self.database.reactions,
                                   lowBound=0, upBound=M, cat='Continuous')
        yf = pulp.LpVariable.dicts("yf", self.database.reactions,
                                   lowBound=0, upBound=1, cat='Binary')
        yb = pulp.LpVariable.dicts("yb", self.database.reactions,
                                   lowBound=0, upBound=1, cat='Binary')
        a = pulp.LpVariable.dicts("a", self.database.reactions,
                                  lowBound=0, upBound=1, cat='Binary')
        G = pulp.LpVariable.dicts("G", self.database.reactions,
                                  lowBound=-M, upBound=M, cat='Continuous')

        # pdb.set_trace()
        # get the set of reactions belong to GSM and fulldb
        if params['use_heterologous_steps'] == True:
            rxns_tsv = model_files['reactions_file']['path']
            rxns_df = pd.read_table(rxns_tsv)
            rxns_df.head()
            native_reactions = rxns_df[rxns_df["in model"]==1]['id'].tolist()
            hetere_reactions = rxns_df[rxns_df["in model"]==0]['id'].tolist()

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

                continue

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

        
        if params['use_heterologous_steps'] == True:
            # Min-Rxn objective
            if self.objective == 'MinRxn':
                condition = pulp.lpSum([yf[j] + yb[j]
                                        for j in native_reactions
                                        if self.database.rxntype[j] != 4]) + \
                            pulp.lpSum([10*yf[j] + 10*yb[j]
                                        for j in hetere_reactions
                                        if self.database.rxntype[j] != 4])
                lp_prob += condition, "MinRxn"

            # Min-Flux objective
            elif self.objective == 'MinFlux':
                condition = pulp.lpSum([vf[j] + vb[j]
                                        for j in native_reactions
                                        if self.database.rxntype[j] != 4]) + \
                            pulp.lpSum([10*vf[j] + 10*vb[j]
                                        for j in hetere_reactions
                                        if self.database.rxntype[j] != 4])
                lp_prob += condition, "MinFlux"
        else:
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
        # for i in self.database.metabolites:
        for i in self.database.S.keys():
            # # If metabolites not involve in any reactions
            # if i not in self.database.S:
            #     continue
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
        # lp_prob.writeLP("./test_OptStoic.lp")
        return lp_prob, v, vf, vb, yf, yb, a, G

    def solve(self, params,model_files,exclude_existing_solution=False, outputfile="OptStoic_pulp_result.txt", max_iteration=None):
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
        lp_prob, v, vf, vb, yf, yb, a, G = self.create_minflux_problem(params,model_files)

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
            lp_prob.writeLP("OptStoic_SEED.lp", mip=1)  #optional
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

                # result_output.write("%s = %.8f\n" % (self.objective, pulp.value(lp_prob.objective)))
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


def cleanup(directory=None):
    """
    Clean up after the job.  At the moment this just means removing the working
    directory, but later could mean other things.
    """

    try:
        shutil.rmtree(directory, ignore_errors=True)  # it would not delete if fold is not empty
        # need to iterate each entry
    except IOError, e:
        log("Unable to remove working directory {0}".format(directory))
        raise

def setupWorkingDir(directory=None):
    """
    Clean up an existing workingdir and create a new one
    """
    try:
        if os.path.exists(directory):
            _cleanup(directory)
        os.mkdir(directory)
    except IOError:
        log("Unable to setup working dir {0}".format(directory))
#         raise

def test_write_fba_model(pathways,db,model_files,workspace_name,res_dir):
    # rxn_df = pd.read_table(data_dir+'/iMR1_799-reactions.tsv')
    # met_df = pd.read_table(data_dir+'/iMR1_799-compounds.tsv')
    rxn_df = pd.read_table(model_files['reactions_file']['path'])
    met_df = pd.read_table(model_files['compounds_file']['path'])

    Sji = db.Sji     

    all_rxns_df = pd.DataFrame()
    for ind, res in sorted(pathways.iteritems()):
        reaction_ids=res['reaction_id']
        rxn_pathway = rxn_df[rxn_df['id'].isin(reaction_ids)]
        rxn_pathway['reference'] = ind
        all_rxns_df = all_rxns_df.append(rxn_pathway)

    rxn_file = os.path.join(res_dir, "rxn_file.tsv")            
    # rxn_file = data_dir+ "/rxn_file.tsv"
    # with open(rxn_file, 'w') as rf:
    #         rf.write(all_rxns_df)
    all_rxns_df.to_csv(rxn_file,index=False,sep='\t')

    all_mets_list = []
    for r in list(set(all_rxns_df['id'].tolist())):
        all_mets_list.extend(Sji[r])
    all_mets_list= list(set(all_mets_list))
    met_pathway_df = met_df[met_df['id'].isin(all_mets_list)]
    
    cpd_file = os.path.join(res_dir, "cpd_file.tsv")
    # cpd_file = data_dir+"/cpd_file.tsv"   
    # with open(cpd_file, 'w') as cf:
    #     cf.write(met_pathway_df)
    met_pathway_df.to_csv(cpd_file,index=False,sep='\t')

    # upload those files as a model and get the reference back.
    # see here for details:
    # https://github.com/kbaseapps/fba_tools/blob/master/fba_tools.spec#L524
    
    callback_url = os.environ['SDK_CALLBACK_URL']
    fba_client = fba_tools(callback_url)
    model_upa = fba_client.tsv_file_to_model({
        'model_file': {'path': rxn_file},
        'compounds_file': {'path': cpd_file},
        'workspace_name': workspace_name,#self.getWsName(),
        'model_name': 'pathways_from_optStoic',
        'biomass': []
    })

    print('UPLOAD A MODEL - GOT UPA')
    print(model_upa['ref'])

    # pprint(self.getWsClient().get_objects2({'objects': [{'ref': model_upa['ref']}]}))

def run_minRxnFlux(optSotic_result_dict, config, params, model_files):

    res_dir = os.path.join(config['scratch'], "optstoic_out")
    data_dir = os.path.join(config['scratch'], "optstoic_input")

    # log("optStoic output dir is {0}".format(res_dir))
    setupWorkingDir(data_dir)
    setupWorkingDir(res_dir)

    db = load_db_rxns(optSotic_result_dict, model_files,data_dir)
    logging.info("Database loaded!")

    # df_bounds = pd.read_csv('specific_bounds.csv', index_col=0)
    # specific_bounds = df_bounds.to_dict(orient='index')
    specific_bounds = {}
    for met,stoich in optSotic_result_dict.iteritems():
        print met
        if met == "C00080":
            specific_bounds['EX_'+met] = {'LB': -10, 'UB': 10}
        else:
            specific_bounds['EX_'+met] = {'LB': stoich, 'UB': stoich}

    # pulp_solver = load_pulp_solver()
    pulp_solver = SCIP_CMD(path=None, keepFiles=0, mip=1, msg=1, options=[])


    MAX_ITERATION = 1
    USE_LOOPLESS = False
    CLEANUP = True #set to true if you want to delete all .sol and .lp files

    test = MinRxnFlux(db, objective='MinFlux',
                    specific_bounds=specific_bounds,
                    max_iteration=MAX_ITERATION,
                    pulp_solver=pulp_solver,
                    data_filepath=data_dir,
                    result_filepath=res_dir,
                    M=1000)

    #### solve and draw pathway to figures
    lp_prob, pathways = test.solve(params,model_files,outputfile='test_optstoic_mac.txt')

    workspace_name = params['workspace_name']
    test_write_fba_model(pathways,db,model_files,workspace_name,res_dir)
    #----to do: draw pathways----
    # # Creating kegg model and drawing pathways
    # f = open(os.path.join(res_dir, 'test_KeggModel.txt'), 'w+')
    # pathway_objects = []

    # for ind, res in sorted(test.pathways.iteritems()):
    #     p = Pathway(id=ind, name='OptStoic', reaction_ids=res['reaction_id'], fluxes=res['flux'])
    #     p.rearrange_reaction_order()
    #     pathway_objects.append(p)
    #     generate_kegg_model(p, filehandle=f)
    #     graph_title = "{0}_P{1}".format(p.name, p.id)
    #     # draw_pathway(p, imageFileName=os.path.join(res_dir+'/pathway_{0:03d}'.format(p.id)),
    #     #             imageFormat='png', graphTitle=graph_title, debug=True)
    #     draw_pathway(p, imageFileName=os.path.join(res_dir+'/pathway_{0:03d}'.format(p.id)),
    #                 imageFormat='png', graphTitle=graph_title, debug=True)
    # f.close()
    # print "Generate kegg_model and draw pathway: Pass!"

    callback_url = os.environ['SDK_CALLBACK_URL']
    report_maker = CreateReport(callback_url, config['scratch'])
    result = report_maker.run(params)
    return result

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
        if 'BIOMASS' in rxn: continue
        reactions.append(rxn)
        try:
            # pdb.set_trace()
            reactants = parse_equation(df.equation[i])
        except:
            raise ValueError("can't parse equation {} - {}".format(i, df.equation[i]))
        
        for cpd in reactants:
            if 'cpd' in cpd['cpd']: 
                cpd_name = cpd['cpd'] + '_' + cpd['compartment']
            else: 
                cpd_name = cpd['cpd']
            # if cpd['cpd'] not in s_matrix:
            #     s_matrix[cpd['cpd']] = dict()
            # if rxn not in s_matrix[cpd['cpd']]:
            #     s_matrix[cpd['cpd']][rxn] = dict()
            # s_matrix[cpd['cpd']][rxn] = cpd['stoich']
            if cpd_name not in s_matrix:
                s_matrix[cpd_name] = dict()
            if rxn not in s_matrix[cpd_name]:
                s_matrix[cpd_name][rxn] = dict()
            s_matrix[cpd_name][rxn] = cpd['stoich']
    return s_matrix

def write_to_required_format(datafiles,data_dir):
    rxns_tsv = datafiles['rxns']
    mets_tsv = datafiles['mets']

    rxns_df = pd.read_table(rxns_tsv)
    # print rxns_df.head(5)
    rxns_df.to_csv(os.path.join(data_dir,'optstoic_v3_reactions.txt'), header=False, index=False,columns=['id'])

    S_matrix = build_s_matrix(rxns_df)
    with open(os.path.join(data_dir,'optstoic_v3_Sij_dict.json'), 'w') as fp:
        json.dump(S_matrix, fp)
    # reaction type (0 = forward irreversible;
    # 1 = reversible; 2 = backward irreversible; 4 = exchange)
    rxn_type_df = rxns_df[['id','direction']]
    rxn_type_df['direction'] = rxn_type_df['direction'].map({'>': 0, '=': 1,'<':'2'})
    rxn_type_df.to_csv(os.path.join(data_dir,'optstoic_v3_reactiontype.txt'), header=False, index=False,sep=' ')

    mets_df = pd.read_table(mets_tsv)
    mets_df.to_csv(os.path.join(data_dir,'optstoic_v3_metabolites.txt'), header=False, index=False,columns=['id'])


def load_db_rxns(optSotic_result_dict, model_files, data_dir):
    # conver to the original file format
    datafiles ={}
    datafiles['mets'] = model_files['compounds_file']['path']# data_dir+'/iMR1_799-reactions.tsv'
    datafiles['rxns'] = model_files['reactions_file']['path']# data_dir+'/iMR1_799-compounds.tsv'

    write_to_required_format(datafiles, data_dir)
    
    dbdict = {
        'Sij': 'optstoic_v3_Sij_dict.json',
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
        user_defined_export_rxns_Sji['EX_'+met] = {met: -1}


    user_defined_export_rxns_Sij = Database.transpose_S(
        user_defined_export_rxns_Sji
    )
    # print user_defined_export_rxns_Sij
    # pdb.set_trace()
    DB.set_database_export_reaction(user_defined_export_rxns_Sij)

    return DB

if __name__ == '__main__':
    # parameters = {'C00033': 3.0,
    #             'C00267': -1.0,
    #             'C00080': 3.0}
    # run_minRxnFlux(parameters)
    optSotic_result_dict = {'cpd00067_c0': 3.0,
                            'cpd00027_c0': -1.0,
                            'cpd00029_c0': 3.0
                            }
    # load_db_rxns(data_dir,optSotic_result_dict)
    run_minRxnFlux(optSotic_result_dict)
