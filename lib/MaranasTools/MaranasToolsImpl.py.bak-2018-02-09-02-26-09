# -*- coding: utf-8 -*-
#BEGIN_HEADER
# The header block is where all import statments should live
import os
from Bio import SeqIO
from pprint import pprint, pformat
from AssemblyUtil.AssemblyUtilClient import AssemblyUtil
from KBaseReport.KBaseReportClient import KBaseReport
from optStoic import *
from minRxnFlux import *
from steadycom import loop_for_steadycom
#END_HEADER


class MaranasTools:
    '''
    Module Name:
    MaranasTools

    Module Description:
    A KBase module: MaranasTools
This sample module contains one small method - filter_contigs.
    '''

    ######## WARNING FOR GEVENT USERS ####### noqa
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    ######################################### noqa
    VERSION = "0.0.1"
    GIT_URL = "https://github.com/fromseto/MaranasTools.git"
    GIT_COMMIT_HASH = "5fd28712673259afabc69e4029fb4e0a3de60fe3"

    #BEGIN_CLASS_HEADER
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.callback_url = os.environ['SDK_CALLBACK_URL']
        self.config = config
        self.shared_folder = config['scratch']
        #END_CONSTRUCTOR
        pass


    def run_optstoic(self, ctx, params):
        """
        :param params: instance of type "OptStoicParams" (model - the FBA
           model to use as a basis for modification start_compound - the
           initial compound to be used as a source for the pathway
           target_compound - the target compound to maximize yield for in the
           pathway max_steps - the maximum number of steps to allow in the
           optimized pathway - any pathway created that has more than this
           number of steps is disqualified use_heterologous_steps - allows
           adding dG_threshold - a threshold free energy value to further
           constrain the path optimization) -> structure: parameter "model"
           of type "model_upa" (An X/Y/Z style reference to an FBA model.),
           parameter "start_compound" of type "compound_id" (The id of a
           compound that exists either in the model or in the biochemistry.),
           parameter "target_compound" of type "compound_id" (The id of a
           compound that exists either in the model or in the biochemistry.),
           parameter "max_steps" of Long, parameter "use_heterologous_steps"
           of type "boolean" (A boolean - 0=false, 1=true @range (0, 1)),
           parameter "dG_threshold" of Double, parameter "workspace_name" of
           String
        :returns: instance of type "OptStoicOutput" (report_name - name of
           the report object that gets generated. report_ref - UPA of the
           report object that gets generated.) -> structure: parameter
           "report_name" of String, parameter "report_ref" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN run_optstoic
        # workspace_name = params['workspace_name']
        print params.keys()
        stoic_dic,model_files = simulate_optStoic(params,self.config)
        # stoic_dic = simulate_optStoic(start_compound, target_compound)

        # stoic_dic = {'cpd00067_c0': 3.0,
        #             'cpd00027_c0': -1.0,
        #             'cpd00029_c0': 3.0
        #             }
        output = run_minRxnFlux(stoic_dic, self.config, params, model_files)
        #END run_optstoic

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method run_optstoic return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]

    def run_steadycom(self, ctx, params):
        """
        :param params: instance of type "SteadyComParams" -> structure:
           parameter "model_inputs" of list of type "ModelInput" ->
           structure: parameter "model_upa" of String, parameter "fixed_gr"
           of Double, parameter "medium_upa" of String, parameter
           "flux_output" of String, parameter "workspace_name" of String
        :returns: instance of type "SteadyComOutput" -> structure: parameter
           "report_name" of String, parameter "report_ref" of String,
           parameter "flux_output" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN run_steadycom
        for key, value in params.iteritems():
          print key,':',value
        output = loop_for_steadycom(params,self.config,self.callback_url)
        #END run_steadycom

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method run_steadycom return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]
    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK",
                     'message': "",
                     'version': self.VERSION,
                     'git_url': self.GIT_URL,
                     'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]
