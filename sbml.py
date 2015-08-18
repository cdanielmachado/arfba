
from collections import OrderedDict
from libsbml import SBMLReader, SBMLDocument, SBMLWriter
from framed.core.models import GPRConstrainedModel
from framed.io_utils.sbml import _load_compartments, _load_metabolites, _load_reactions, _load_stoichiometry, _load_cb_parameters, _load_gpr


INHIBITOR_TAG = 'inhibitor'
ACTIVATOR_TAG = 'activator'


class AllostericModel(GPRConstrainedModel):

    def __init__(self, model_id):
        """
        Arguments:
            model_id : String -- a valid unique identifier
        """
        GPRConstrainedModel.__init__(self, model_id)
        self.regulation = OrderedDict()
        self._m_r_reg_lookup = None
        self._r_m_reg_lookup = None

    def _clear_temp(self):
        GPRConstrainedModel._clear_temp(self)
        self._m_r_reg_lookup = None
        self._r_m_reg_lookup = None


    def add_regulators(self, regulators):
        for m_id, r_id, kind in regulators:
            if m_id in self.metabolites and r_id in self.reactions:
                self.regulation[(m_id, r_id)] = kind


    def remove_reactions(self, id_list):
        GPRConstrainedModel.remove_reactions(self, id_list)
        for (m_id, r_id) in self.regulation:
            if r_id in id_list:
                del self.regulation[(m_id, r_id)]
        self._clear_temp()


    def metabolite_reaction_regulatory_lookup_table(self):

        if not self._m_r_reg_lookup:
            self._m_r_reg_lookup = OrderedDict([(m_id, OrderedDict()) for m_id in self.metabolites])

            for (m_id, r_id), kind in self.regulation.items():
                self._m_r_reg_lookup[m_id][r_id] = kind

        return self._m_r_reg_lookup


    def reaction_metabolite_regulatory_lookup_table(self):

        if not self._r_m_reg_lookup:
            self._r_m_reg_lookup = OrderedDict([(r_id, OrderedDict()) for r_id in self.reactions])

            for (m_id, r_id), kind in self.regulation.items():
                self._r_m_reg_lookup[r_id][m_id] = kind

        return self._r_m_reg_lookup


def load_allosteric_model(filename):

    reader = SBMLReader()
    document = reader.readSBML(filename)
    sbml_model = document.getModel()

    if sbml_model is None:
        raise IOError('Failed to load model.')

    model = AllostericModel(sbml_model.getId())
    model.add_compartments(_load_compartments(sbml_model))
    model.add_metabolites(_load_metabolites(sbml_model))
    model.add_reactions(_load_reactions(sbml_model))
    model.add_stoichiometry(_load_stoichiometry(sbml_model))
    model.add_regulators(_load_regulators(sbml_model))
    bounds, coefficients = _load_cb_parameters(sbml_model)
    model.set_bounds(bounds)
    model.set_objective_coefficients(coefficients)
    genes, rules = _load_gpr(sbml_model)
    model.add_genes(genes)
    model.set_rules(rules)

    return model


def _load_regulators(model):

    modifiers = [(modifier.getSpecies(), reaction.getId(), modifier.getNotes().getChild(0).toString())
                 for reaction in model.getListOfReactions()
                 for modifier in reaction.getListOfModifiers()]

    regulators = [(m_id, r_id, 1 if tag == ACTIVATOR_TAG else -1)
                  for m_id, r_id, tag in modifiers
                  if tag == ACTIVATOR_TAG or tag == INHIBITOR_TAG]

    return regulators
