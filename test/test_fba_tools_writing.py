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
