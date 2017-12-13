# Author: Chiam Yu Ng, Lin Wang

from __future__ import division
from .config import cofactors, default_params, rxnSji, cofactorsList, kegg_compound
import os
from collections import OrderedDict
import graphviz as gv
import logging
import math

class Reaction(object):
    """Reaction class"""

    def __init__(self, rid=None, flux=1, metabolites={},
                 equation='', reversible=True):
        self.rid = rid
        self.flux = flux
        self.metabolites = metabolites
        self.equation = equation
        self.reversible = reversible

    def get_rid(self):
        return self.rid

    def get_metabolites(self):
        return self.metabolites

    def get_equation(self):
        return self.equation

    def autoset_metabolites(self):
        if len(self.metabolites) > 0:
            print "WARNING: Metabolites exists!"
        else:
            print "Retrieving metabolites from default database"
            self.metabolites = rxnSji[self.rid]
        return self.metabolites

    def flux(self):
        return self.flux

    @property
    def reactants(self):
        if self.flux > 0:
            return [k for k, v in self.metabolites.items() if v < 0]
        else:
            return [k for k, v in self.metabolites.items() if v > 0]

    @property
    def products(self):
        if self.flux > 0:
            return [k for k, v in self.metabolites.items() if v > 0]
        else:
            return [k for k, v in self.metabolites.items() if v < 0]

    def __str__(self):
        return "Reaction('%s')" % self.rid

    def __repr__(self):
        return "Reaction('%s')" % self.rid

    def set_equation(self):
        """Write equation in the direction of the flux"""
        if len(self.equation) != 0:
            print "WARNING: Equation exists!"
        else:
            if len(self.metabolites) == 0:
                print "Metabolites are not available!"
                print "Auto-updating metabolites..."
                self.autoset_metabolites()

            temp_list = []
            for cpd in sorted(self.reactants):
                coeff = abs(self.metabolites[cpd])
                if coeff == 1:
                    temp_list.append(cpd)
                else:
                    temp_list.append('%1.0f %s' % (coeff, cpd))

            eqStr = ' + '.join(temp_list)
            eqStr += ' <=> '

            temp_list = []
            for cpd in sorted(self.products):
                coeff = abs(self.metabolites[cpd])
                if coeff == 1:
                    temp_list.append(cpd)
                else:
                    temp_list.append('%1.0f %s' % (coeff, cpd))
            eqStr += ' + '.join(temp_list)
            self.equation = eqStr
        return self.equation

    @classmethod
    def create_Reaction_list_from_dict(cls, dataDict, excludeExchangeRxn=True):
        """
        Make a list of Reaction object from dataDict,
        excluding exchange reaction

        E.g.
        dataDict = {'reaction_id': ['R00001', 'R00002'], 'flux': [-1, 1]}
        output = [Reaction('R00001'), Reaction('R00002')]

        Keyword arguments:
        dataDict -- dictionary with reaction_id and flux
        excludeExchangeRxn -- Exclude all exchange reactions in the list
                              (default true)
        """
        RxnObjList = []
        for i in range(len(dataDict['reaction_id'])):
            if excludeExchangeRxn:
                if 'EX_' in dataDict['reaction_id'][i]:
                    continue
            tempRxn = cls(dataDict['reaction_id'][i], dataDict['flux'][i])
            # Get the metabolites dictionary {'C00001': -1, ...} for each
            # reaction
            tempRxn.metabolites = rxnSji[tempRxn.rid]
            RxnObjList.append(tempRxn)
        return RxnObjList

class Pathway(object):
    """OptStoic Pathway class"""

    def __init__(self, id=None, name=None,
                 reaction_ids=[], fluxes=[], reactions=None,
                 sourceSubstrateID='C00031', endSubstrateID='C00022',
                 total_flux_no_exchange=None, note={}):
        """
        A Pathway instance can be initialized by either
            (a) List of reaction_ids and fluxes (Let reactions = None)
            (b) List of reaction instances (reaction_ids and fluxes
                will be populated)

        Keyword arguments:
        id -- pathway id
        name -- pathway name
        reaction_ids -- list of reaction IDs (kegg_id) in the pathway
        fluxes -- list of reaction fluxes corresponding to the reaction_ids
        reactions  -- list of reaction object that form the pathway
        sourceSubstrateID -- Kegg compound ID of the source metabolite of
                             the pathway
        endSubstrateID -- Kegg compound ID of the end metabolite of the pathway
        note -- (For debugging purpose) modelstat and solvestat can be added
        """
        self.id = id
        self.name = name
        self.note = note
        self.total_flux_no_exchange = total_flux_no_exchange

        # Iniatilize pathway object using list of reaction_ids and fluxes
        if reactions is None:
            # Check if both reaction_ids and fluxes are list and
            # contain the same number of item
            assert (isinstance(reaction_ids, list) == 1)
            assert (isinstance(fluxes, list) == 1)
            assert len(reaction_ids) == len(fluxes), "number of reactions must equal number of fluxes!"

            # Change EX_h+ to EX_hplus as optstoic pulp fail to read "+" symbol
            self.reaction_ids = ["EX_hplus" if x == "EX_h+" else x for x in reaction_ids]

            self.fluxes = fluxes
            # Create list of reaction objects
            self.reactions = Reaction.create_Reaction_list_from_dict(
                {'reaction_id': self.reaction_ids, 'flux': self.fluxes},
                excludeExchangeRxn=True)

        # Iniatilize pathway object using list of reaction objects
        else:
            self.reactions = reactions
            self.fluxes = [r.flux for r in reactions]
            self.reaction_ids = [r.rid for r in reactions]

        self.reaction_ids_no_exchange = [r for r in reaction_ids if 'EX_' not in r]

        if not self.total_flux_no_exchange:
            self.total_flux_no_exchange = sum(map(
                abs, [r.flux for r in self.reactions]))

        self.rxn_flux_dict = dict(zip(self.reaction_ids, self.fluxes))

        try:
            self.nATP = self.rxn_flux_dict['EX_atp']
        except:
            self.nATP = None
        self.sourceSubstrateID = sourceSubstrateID
        self.endSubstrateID = endSubstrateID

    def get_pathway_dict(self):
        """
        return a dictionary of the {reaction:flux}
        """
        return dict(zip(self.reaction_ids, self.fluxes))

    def update_nATP(self):

        try:
            self.nATP = self.rxn_flux_dict['EX_atp']
        except:
            self.nATP = None

    def get_total_flux(self):
        """
        return total flux through the pathway
        """
        return sum(map(abs, self.fluxes))

    def get_total_flux_no_exchange(self):
        """
        return total flux through the pathway excluding exchange reactions
        """
        return self.total_flux_no_exchange

    def get_modelstat(self):
        if 'modelstat' in self.note:
            return self.note['modelstat']
        else:
            return None

    def get_solvestat(self):
        if 'solvestat' in self.note:
            return self.note['solvestat']
        else:
            return None

    def get_reaction_involving_reactant(self, substrate_ID):
        """get reactions that involve reactant substrate_ID"""
        output = []
        for rxn in self.reactions:
            if substrate_ID in rxn.reactants:
                output.append(rxn)
        return output

    def rearrange_reaction_order(self):
        """
        Try to implement a topological sorting of the pathway
        (Todo: use a different algorithm)

        """
        # Exclude exchange reaction from being rearranged
        sortedRxn = []
        flag = 1

        next_substrate = [self.sourceSubstrateID]
        # Find substrateID in list of reactants
        while flag == 1:
            rstore = []
            for subs in next_substrate:
                if subs == self.endSubstrateID:
                    flag = 0
                    break
                r1 = self.get_reaction_involving_reactant(subs)

                if len(r1) == 0:
                    # print "\nWarning: No reaction found using %s" %subs
                    if len(next_substrate) == 1:
                        flag = 0
                        break
                else:
                    sortedRxn.extend([r.rid for r in r1])
                rstore += r1
            next_substrate = []
            for rxn in rstore:
                next_substrate.extend(rxn.products)
            next_substrate = list(set(next_substrate) - cofactors)

        # Add all reactions to the sorted reactions
        # (as the duplicates will be removed in the following command)
        sortedRxn += self.reaction_ids
        # Make unique ordered list
        sortedRxnUnique = list(OrderedDict.fromkeys(sortedRxn))

        # Raise error if the number of reactions changes after processing
        assert len(sortedRxnUnique) == len(self.reaction_ids), "Error: \
        the number of reactions does not match after processing"

        # Sort all the flux according to the order of the reaction ID
        sortedRxnFlux = [self.rxn_flux_dict[rxn] for rxn in sortedRxnUnique]

        self.reaction_ids = sortedRxnUnique
        self.fluxes = sortedRxnFlux
        self.reactions = Reaction.create_Reaction_list_from_dict(
            {'reaction_id': sortedRxnUnique, 'flux': sortedRxnFlux})
        return self

    def __repr__(self):
        return "<OptStoicPathway(id='%s', numRxn='%s', nATP='%s')>" % (
            self.id, len(self.reaction_ids), self.nATP)

    # @staticmethod
    def get_pathway_similarity_index(self, pathway2):
        """Calculate the jaccard index of two pathways"""
        a = set(self.reaction_ids)
        b = set(pathway2.reaction_ids)
        idscore = len(a & b) / len(a | b)
        return idscore

    def get_pathway_similarity_index_no_exchange(self, pathway2):
        """Calculate the jaccard index of two pathways"""
        a = set(self.reaction_ids_no_exchange)
        b = set(pathway2.reaction_ids_no_exchange)
        idscore = len(a & b) / len(a | b)
        return idscore


    def is_same_pathway_with(self, another_pathway):
        idscore = self.get_pathway_similarity_index(another_pathway)
        if idscore == 1:
            return 1
        else:
            return 0

    @staticmethod
    def pathways_to_dict(list_of_pathways):
        """
        Output list of Pathway instances as dictionary
        """
        output = {}
        for p in list_of_pathways:
            output[p.id] = dict(pathway=p.get_pathway_dict(),
                                num_reaction=len(p.reaction_ids),
                                total_flux_no_exchange=p.get_total_flux_no_exchange(),
                                modelstat=p.get_modelstat(),
                                solvestat=p.get_solvestat())

        return output

class Draw(object):
    kegg_compound['C00009'] = 'Pi'
    kegg_compound['C00013'] = 'PPi'
    kegg_compound['C00236'] = '1,3-Bisphospho-D-glycerate'
    kegg_compound['C00111'] = 'dihydroxyacetone phosphate'

    REACTION_FONT_SIZE = '28'
    color_configs = {}
    color_configs['light'] = dict(COFACTOR_SHAPE="ellipse",  # "diamond"
                                  OTHER_COFACTOR_COLOR='#B6B6B6',
                                  NONCOFACTOR_SHAPE="plaintext",  # "box"
                                  NONCOFACTOR_COLOR='transparent',  # "#FFFFFF"
                                  REACTION_COLOR="#512DA8",
                                  RXN_NODE_COLOR="#323232",
                                  EDGE_COLOR="#323232",  # "#505050"
                                  BACKGROUND_COLOR="transparent",
                                  ALL_FONT_COLOR="black")


    color_configs['light']['colorMapping'] = {
        'C00002': '#F05456', 'C00008': '#FFC000',
        'C00003': '#149B76', 'C00004': '#149B76',
        'C00005': '#2393CB', 'C00006': '#2393CB'
    }

    def generate_kegg_model(pathway, params=default_params,
        filehandle=None, add_ratio_constraints=False):
        """
        Convert the pathway to KEGG model format
        (as input for Component Contribution/MDF)

        Keyword arguments:
        pathway -- pathway object
        params -- Kegg model parameters (default parameters are given)
        filehandle -- If a text file handle is provided,
                      it write the model text to the file (default None)

        """
        params['ENTRY'] = "{0}_{1}ATP_P{2}".format(pathway.name,
                                                   pathway.nATP, pathway.id)
        params['NAME'] = "{0}_{1}ATP_P{2}".format(pathway.name,
                                                  pathway.nATP, pathway.id)

        modeltext = """\
        ENTRY\t\t{ENTRY}
        SKIP\t\t{SKIP}
        NAME\t\t{NAME}
        PH\t\t\t{PH}
        I\t\t\t{I}
        T\t\t\t{T}
        C_RANGE\t\t{C_RANGE[0]:.0e} {C_RANGE[1]:.0e}\n""".format(**params)

        all_bounds = params['BOUND']

        if add_ratio_constraints:
            all_bounds = params['RATIO_BOUND']

            # write the ratios
            for i, (cids, ratios) in enumerate(sorted(params['RATIO'].iteritems())):
                if i == 0:
                    modeltext += "RATIO\t\t{C[0]} {C[1]} {B[0]:e} {B[1]:e}\n".format(C=cids, B=ratios)
                else:
                    modeltext += "\t\t\t{C[0]} {C[1]} {B[0]:e} {B[1]:e}\n".format(C=cids, B=ratios)

        # write the bounds
        for i, (cid, bounds) in enumerate(sorted(all_bounds.iteritems())):
            if i == 0:
                modeltext += "BOUND\t\t{0} {1[0]:e} {1[1]:e}\n".format(cid, bounds)
            else:
                modeltext += "\t\t\t{0} {1[0]:e} {1[1]:e}\n".format(cid, bounds)

        # write the reactions (in the direction of the flux)
        for i, rxn in enumerate(pathway.reactions):
            rxn.set_equation()
            if i == 0:
                modeltext += "REACTION\t{0} {1} (x{2:1.2f})\n".format(
                    rxn.rid, rxn.equation, abs(rxn.flux))
            else:
                modeltext += "\t\t\t{0} {1} (x{2:1.2f})\n".format(
                    rxn.rid, rxn.equation, abs(rxn.flux))

        modeltext += "///\n"
        if filehandle:
            filehandle.write(modeltext)

        return modeltext

    # ########################################

    def load_global_styles(colorConfig):
        """
        Create a global styles dictionary for all Graphviz graph
        """
        colorMapping = colorConfig['colorMapping']
        for c in cofactorsList:
            if c not in colorMapping:
                colorMapping[c] = colorConfig['OTHER_COFACTOR_COLOR']

        global_styles = {
            'graph': {
                'fontsize': '20',
                'fontname': 'Helvetica',
                'bgcolor': colorConfig['BACKGROUND_COLOR'],
                # 'rankdir': 'BT',
            },
            'nodes': {
                'fontname': 'Helvetica',
                'fontsize': '30',
                # 'shape': 'hexagon',
                'fontcolor': colorConfig['ALL_FONT_COLOR'],
                # 'color': 'white',
                # 'style': 'filled',
                # 'fillcolor': '#006699',
            },
            'edges': {
                # 'style': 'dashed',
                # 'color': 'white',
                # 'arrowhead': 'open',
                'fontname': 'Helvetica',
                'fontsize': '24',
                # 'fontcolor': 'white',
            }
        }
        return global_styles, colorMapping


    def apply_styles(graph, styles):
        """
        Apply the styles to a Graphviz graph.
        """
        graph.graph_attr.update(
            ('graph' in styles and styles['graph']) or {}
        )
        graph.node_attr.update(
            ('nodes' in styles and styles['nodes']) or {}
        )
        graph.edge_attr.update(
            ('edges' in styles and styles['edges']) or {}
        )
        return graph


    def draw_pathway(Pathway, imageFileName=None, imageFormat='png',
        graphTitle='', scaleLineWidth=False, scalingFactor=200.0,
        cleanup=True, engine='dot', debug=False):
        """
        Draw a digraph for a Pathway objects and render it as
        the given imageFormat using Graphviz.

        Keyword arguments:
        Pathway -- A Pathway object (pathway.py)
        imageFileName -- Name of the output file (default Pathway.name)
        imageFormat -- Any format that Graphviz can support (default 'png')
        graphTitle -- Title of the output graph
        scaleLineWidth -- If true, scale the penwidth of an edge
                          to a value between 1 and 10. This is useful when
                          fluxes are too large.
                          Else, the penwidth of an edge is absolute value
                          of the flux value.  (default False)
        scalingFactor -- If scaleLineWidth is true,
                         penwidth = (abs(flux)/scalingFactor) * 10 + 1.
                         (E.g. Use the maximum flux values of a
                         pathway as the scaling Factor)
        cleanup -- delete the ".dot" file after drawing
        engine -- Graphviz layout engine used to render the graph.
                  Layout engines = {'circo', 'dot', 'fdp', 'neato', 'nop1', 'nop2',
                                    'osage', 'patchwork', 'sfdp', 'twopi'}
        """
        if debug:
            logging.warning('Debug mode: Drawing pathway.')

        colorConfig = color_configs['light']

        global_styles, colorMapping = load_global_styles(colorConfig)
        metabolite_fontname = 'helvetica bold'
        #bypass issue with transparent color for vector image in AI
        if imageFormat.lower() in ['svg', 'eps']:
            colorConfig['NONCOFACTOR_COLOR'] = '#FFFFFF'

        g = gv.Digraph('G', format=imageFormat, engine=engine)

        # g.graph_attr['ratio']='fill'
        g.graph_attr['size'] = "10, 10"
        if imageFormat == 'png':
            g.graph_attr['dpi'] = '300'
        elif imageFormat == 'svg':
            g.graph_attr['dpi'] = '72'
        g.graph_attr['forcelabels'] = 'true'
        g.graph_attr['labelloc'] = 't'  # top or 'b'
        g.graph_attr['label'] = graphTitle
        g = apply_styles(g, global_styles)
        r_counter = 1

        # Automatically use scaling factor if flux > 10
        # (scaling factor is set to the nearest 100)
        all_f = [abs(f) for f in Pathway.fluxes]

        if max(all_f) > 10:
            scaleLineWidth = True
            scalingFactor = 10**(math.ceil(math.log(max(all_f), 10)))

        for rxn in Pathway.reactions:

            g.node(rxn.rid, shape='point',
                   color=colorConfig['RXN_NODE_COLOR'],
                   xlabel=rxn.rid + '; ' + str(abs(rxn.flux)),
                   fontsize=REACTION_FONT_SIZE,
                   fontcolor=colorConfig['REACTION_COLOR'])

            if scaleLineWidth:
                lineW = '%i' % (10 * abs(rxn.flux) / scalingFactor + 1)

            else:
                if rxn.flux >= 1 or rxn.flux <= -1:
                    lineW = '%s' % (abs(rxn.flux) * 2)

                # penwidth for any flux between -1 < v < 1 is 1.
                else:
                    lineW = '1'

            for met in rxn.reactants:
                if met in cofactorsList:
                    cf = abs(rxn.metabolites[met])
                    if cf > 1:
                        # show stoichiometric coefficient only when it is >= 2
                        clabel = '%i %s' % (cf, kegg_compound[met])
                    else:
                        clabel = kegg_compound[met]
                    g.node(met + '_' + str(r_counter),
                           shape=colorConfig['COFACTOR_SHAPE'],
                           color=colorMapping[met], style="filled", label=clabel)
                    g.edge(met + '_' + str(r_counter), rxn.rid, penwidth=lineW,
                           weight='1', arrowhead="none",
                           color=colorConfig['EDGE_COLOR'])
                else:
                    g.node(met, shape=colorConfig['NONCOFACTOR_SHAPE'],
                           label=kegg_compound[met], fontname=metabolite_fontname,
                           style="filled", color=colorConfig['NONCOFACTOR_COLOR'])
                    g.edge(met, rxn.rid, penwidth=lineW, weight='2',
                           arrowhead="none", color=colorConfig['EDGE_COLOR'])

            for met in rxn.products:
                if met in cofactorsList:
                    cf = abs(rxn.metabolites[met])
                    if cf > 1:
                        clabel = '%i %s' % (cf, kegg_compound[met])
                    else:
                        clabel = kegg_compound[met]
                    g.node(met + '_' + str(r_counter),
                           shape=colorConfig['COFACTOR_SHAPE'],
                           color=colorMapping[met], style="filled", label=clabel)
                    g.edge(rxn.rid, met + '_' + str(r_counter), weight='1',
                           penwidth=lineW, color=colorConfig['EDGE_COLOR'])
                else:
                    g.node(met, shape=colorConfig['NONCOFACTOR_SHAPE'],
                           label=kegg_compound[met], fontname=metabolite_fontname,
                           style="filled", color=colorConfig['NONCOFACTOR_COLOR'])
                    g.edge(rxn.rid, met, penwidth=lineW, weight='2',
                           color=colorConfig['EDGE_COLOR'])
            r_counter += 1

        if imageFileName is None:
            imageFileName = Pathway.name
        g.render(imageFileName, cleanup=cleanup)

        return 1  # g.source


    def test_drawpathway():




if __name__ == "__main__":

    logging.basicConfig(format='%(asctime)s %(levelname)s %(message)s',
                        datefmt='%m/%d/%Y %I:%M:%S %p',
                        level=logging.INFO)
    
    # Test for drawing pathways. Create two version of the pathways.
    testpathway = {'flux': [-1.0, 1.0, 1.0, 1.0, 1.0, -1.0, -1.0, 1.0,
                            1.0, 1.0, -1.0, -1.0, -1.0, -1.0, 2.0, 1.0,
                            1.0, 1.0, -1.0, 1.0],
                   'iteration': 1,
                   'reaction_id': ['R00200', 'R00300', 'R00658', 'R01059',
                                   'R01063', 'R01512', 'R01518', 'R01519',
                                   'R01538', 'R08570', 'EX_glc', 'EX_nad',
                                   'EX_adp', 'EX_phosphate', 'EX_pyruvate',
                                   'EX_nadh', 'EX_atp', 'EX_h2o', 'EX_nadp',
                                   'EX_nadph']}

    p1 = Pathway(id=1, name='OptStoic',
                 reaction_ids=testpathway['reaction_id'],
                 fluxes=testpathway['flux'])
    # outputFilename = 'OptStoic'
    # current_dir = dirname(abspath(__file__))
    # outputFilepath = normpath(join(current_dir,'../result', outputFilename))

    logging.info("Creating 'res' folder in current directory if not exist...")
    outputFilepath = 'res'
    try:
        os.makedirs(outputFilepath)
    except OSError:
        if not os.path.isdir(outputFilepath):
            raise

    draw_pathway(p1, os.path.join(outputFilepath, 'test_drawpath_light'),
                 imageFormat='png', graphTitle='test_pathway_light',
                 scaleLineWidth=False, scalingFactor=200.0)

    logging.info("Testing drawpathway.py: Pass\n")
    return None
