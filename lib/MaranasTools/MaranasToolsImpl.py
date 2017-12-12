# -*- coding: utf-8 -*-
#BEGIN_HEADER
# The header block is where all import statments should live
import os
from Bio import SeqIO
from pprint import pprint, pformat
from AssemblyUtil.AssemblyUtilClient import AssemblyUtil
from KBaseReport.KBaseReportClient import KBaseReport
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
    GIT_URL = ""
    GIT_COMMIT_HASH = ""

    #BEGIN_CLASS_HEADER
    # Class variables and functions can be defined in this block
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        
        # Any configuration parameters that are important should be parsed and
        # saved in the constructor.
        self.callback_url = os.environ['SDK_CALLBACK_URL']
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
           parameter "dG_threshold" of Double
        :returns: instance of type "OptStoicOutput" (report_name - name of
           the report object that gets generated. report_ref - UPA of the
           report object that gets generated.) -> structure: parameter
           "report_name" of String, parameter "report_ref" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN run_optstoic
        print "it is working"
        output = "awesome"
        #END run_optstoic

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method run_optstoic return value ' +
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
