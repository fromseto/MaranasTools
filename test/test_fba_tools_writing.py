# -*- coding: utf-8 -*-
import unittest
import os  # noqa: F401
import json  # noqa: F401
import time
import requests

from os import environ
try:
    from ConfigParser import ConfigParser  # py2
except:
    from configparser import ConfigParser  # py3

from pprint import pprint  # noqa: F401

from Workspace.WorkspaceClient import Workspace as workspaceService
from MaranasTools.MaranasToolsImpl import MaranasTools
from MaranasTools.MaranasToolsServer import MethodContext
from MaranasTools.authclient import KBaseAuth as _KBaseAuth
from fba_tools.fba_toolsClient import fba_tools


class FbaToolsWritingTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        token = environ.get('KB_AUTH_TOKEN', None)
        config_file = environ.get('KB_DEPLOYMENT_CONFIG', None)
        cls.cfg = {}
        config = ConfigParser()
        config.read(config_file)
        for nameval in config.items('MaranasTools'):
            cls.cfg[nameval[0]] = nameval[1]
        # Getting username from Auth profile for token
        authServiceUrl = cls.cfg['auth-service-url']
        auth_client = _KBaseAuth(authServiceUrl)
        user_id = auth_client.get_user(token)
        # WARNING: don't call any logging methods on the context object,
        # it'll result in a NoneType error
        cls.ctx = MethodContext(None)
        cls.ctx.update({'token': token,
                        'user_id': user_id,
                        'provenance': [
                            {'service': 'MaranasTools',
                             'method': 'please_never_use_it_in_production',
                             'method_params': []
                             }],
                        'authenticated': 1})
        cls.wsURL = cls.cfg['workspace-url']
        cls.wsClient = workspaceService(cls.wsURL)
        cls.serviceImpl = MaranasTools(cls.cfg)
        cls.scratch = cls.cfg['scratch']
        cls.callback_url = os.environ['SDK_CALLBACK_URL']

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.wsClient.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')

    def getWsClient(self):
        return self.__class__.wsClient

    def getWsName(self):
        if hasattr(self.__class__, 'wsName'):
            return self.__class__.wsName
        suffix = int(time.time() * 1000)
        wsName = "test_MaranasTools_" + str(suffix)
        ret = self.getWsClient().create_workspace({'workspace': wsName})  # noqa
        self.__class__.wsName = wsName
        return wsName

    def getImpl(self):
        return self.__class__.serviceImpl

    def getContext(self):
        return self.__class__.ctx

    @unittest.skip('skipping test')
    def test_write_fba_model(self):
        cpd_tsv = """id\tname\tformula\tcharge\taliases
cf00001\tcf00001\tnone\t0\tnone
cf00002\tcf00002\tnone\t0\tnone
cf00003\tcf00003\tnone\t0\tnone
cpd00001\tH2O\tH2O\t0\tnone
cpd00002\tATP\tC10H13N5O13P3\t-3\tnone
cpd00003\tNAD\tC21H26N7O14P2\t-1\tnone
cpd00004\tNADH\tC21H27N7O14P2\t-2\tnone
cpd00005\tNADPH\tC21H27N7O17P3\t-3\tnone
cpd00006\tNADP\tC21H26N7O17P3\t-2\tnone
cpd00007\tO2\tO2\t0\tnone
cpd00008\tADP\tC10H13N5O10P2\t-2\tnone
cpd00009\tPhosphate\tHO4P\t-2\tnone"""

        rxn_tsv = """id\tdirection\tcompartment\tgpr\tname\tenzyme\tpathway\treference\tequation
pkr0000001\t>\tc0\tnone\tpkr0000001\tnone\tnone\tnone\t(1) cpd00002[c0] => (1) cpd00008[c0] + (1) cpd00009[c0]"""
        # write out to scratch file system
        cpd_file = os.path.join(self.scratch, "cpd_file.tsv")
        with open(cpd_file, 'w') as cf:
            cf.write(cpd_tsv)

        rxn_file = os.path.join(self.scratch, "rxn_file.tsv")
        with open(rxn_file, 'w') as rf:
            rf.write(rxn_tsv)

        # upload those files as a model and get the reference back.
        # see here for details:
        # https://github.com/kbaseapps/fba_tools/blob/master/fba_tools.spec#L524
        fba_client = fba_tools(self.callback_url)
        model_upa = fba_client.tsv_file_to_model({
            'model_file': {'path': rxn_file},
            'compounds_file': {'path': cpd_file},
            'workspace_name': self.getWsName(),
            'model_name': 'my_test_fbamodel',
            'biomass': []
        })

        print('UPLOAD A MODEL - GOT UPA')
        print(model_upa['ref'])

        pprint(self.getWsClient().get_objects2({'objects': [{'ref': model_upa['ref']}]}))

    def test_fetch_model_file(self):
        fba_client = fba_tools(self.callback_url, service_ver="beta")
        model_files = fba_client.model_to_tsv_file({
            'model_name': 'iMR1_799',
            'workspace_name': 'lqw5322:narrative_1515706033382',#'lqw5322:narrative_1515702128212',
            'fulldb': 1
        })
        pprint(model_files)

    @unittest.skip('skipping test')    
    def test_input_metDB_from_tsv(self):
        fba_client = fba_tools(self.callback_url)
        model_files = fba_client.model_to_tsv_file({
            'model_name': 'iMR1_799',
            'workspace_name': 'lqw5322:narrative_1515706033382',#'lqw5322:narrative_1515702128212',
            'fulldb': 1
        })
        # pprint(model_files)

        mets_tsv = model_files['compounds_file']['path']
        import pandas as pd
        model_df = pd.read_table(mets_tsv,index_col='id')
        # print model_df.head()

        # met_info = {}

        # get elements
        glucose_formula = model_df.loc['cpd00027_c0']['formula']

        # glucose_formula = 'C6H12O'
        import re
        elements = re.findall(r'([A-Z][a-z]*)(\d*)', glucose_formula)
        # for element in elements:
        #     for ele,value in element.iteritems():
        #         met_info[ele] = float(value)
        met_info = dict(elements)
        for key, value in met_info.iteritems():
            if value == '': 
                met_info[key] = 1
        print met_info
        # get charge and dG
        met_info['charge'] = model_df.loc['cpd00027_c0']['charge']
        met_info['dGf'] = model_df.loc['cpd00027_c0']['deltag']

        print met_info
        # met_db = {}
        # for met in list_of_mets:
        #     met_db[met] = met_info

    @unittest.skip('skipping test')     
    def test_run_optstoic_local(self):
        # inputs = { ... define inputs here ... }
        # inputs = {"start_compound":"C00267","target_compound":"C00033","workspace_name":self.getWsName()}
        print "test_run_optstoic"
        inputs = {}
        inputs['mets'] = ['cpd00027_c0','cpd00029_c0','cpd00001_c0'] # 'cpd00001_c0': H2O
        inputs['start_compound'] = 'cpd00027_c0' # glucose
        inputs['target_compound'] = 'cpd00029_c0' # acetate
        # expected_outputs = { ... defined expected results here ... }
        outputs = self.getImpl().run_optstoic(self.getContext(), inputs)[0]
        # insert some assertion that outputs = expected_outputs below.
        # self.assertEqual('awesome', 'aweome')
        print(outputs)
        self.assertIn('report_name', outputs)
        self.assertIn('report_ref', outputs)

    @unittest.skip('skipping test')
    def test_run_optstoic_online(self):
        # inputs = { ... define inputs here ... }
        # inputs = {"start_compound":"C00267","target_compound":"C00033","workspace_name":self.getWsName()}
        print "test_run_optstoic"
        inputs = {}
        inputs['model'] = '12219/11/1' # 'cpd00001_c0': H2O
        inputs['start_compound'] = 'cpd00027_c0' # glucose
        inputs['target_compound'] = 'cpd00029_c0' # acetate
        inputs['workspace_name'] = 'lqw5322:narrative_1515706033382'
        # expected_outputs = { ... defined expected results here ... }
        outputs = self.getImpl().run_optstoic(self.getContext(), inputs)[0]
        # insert some assertion that outputs = expected_outputs below.
        # self.assertEqual('awesome', 'aweome')
        print(outputs)
        self.assertIn('report_name', outputs)
        self.assertIn('report_ref', outputs)
