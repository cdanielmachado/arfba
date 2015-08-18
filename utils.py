from copy import deepcopy
from collections import OrderedDict
from framed.core.models import Reaction


def build_perturbed_model(model_ref, perturbation=None):
    model = deepcopy(model_ref)
    if perturbation:
        for r_id, (lb, ub) in perturbation.items():
            model.set_flux_bounds(r_id, lb, ub)
    mapping = make_irreversible(model)
    return model, mapping


def make_irreversible(model):

    mapping = dict()
    to_remove = []
    table = model.reaction_metabolite_lookup_table()
    reg_table = model.reaction_metabolite_regulatory_lookup_table()

    for r_id, reaction in model.reactions.items():
        if reaction.reversible:
            fwd_id = reaction.id + '_f'
            bwd_id = reaction.id + '_b'
            mapping[r_id] = (fwd_id, bwd_id)

            model.add_reaction(Reaction(fwd_id, reaction.name, False))
            model.add_reaction(Reaction(bwd_id, reaction.name, False))

            for m_id, coeff in table[r_id].items():
                model.stoichiometry[(m_id, fwd_id)] = coeff
                model.stoichiometry[(m_id, bwd_id)] = -coeff

            lb, ub = model.bounds[r_id]
            lb_fwd = max(0, lb) if lb is not None else 0
            ub_fwd = max(0, ub) if ub is not None else None
            lb_bwd = max(-ub, 0) if ub is not None else 0
            ub_bwd = max(-lb, 0) if lb is not None else None
            model.set_flux_bounds(fwd_id, lb_fwd, ub_fwd)
            model.set_flux_bounds(bwd_id, lb_bwd, ub_bwd)
            obj = model.objective[r_id]
            model.set_reaction_objective(fwd_id, obj if obj >= 0 else 0)
            model.set_reaction_objective(bwd_id, -obj if obj < 0 else 0)

            model.set_rule(fwd_id, model.rules[r_id])
            model.set_rule(bwd_id, model.rules[r_id])

            for m_id, kind in reg_table[r_id].items():
                model.regulation[(m_id, fwd_id)] = kind
                model.regulation[(m_id, bwd_id)] = kind

            to_remove.append(r_id)

    model.remove_reactions(to_remove)

    return mapping


def compute_turnover(model, v):
    m_r_table = model.metabolite_reaction_lookup_table()
    t = {m_id: 0.5*sum([abs(coeff * v[r_id]) for r_id, coeff in neighbours.items()])
         for m_id, neighbours in m_r_table.items()}
    return t


def merge_fluxes(model, mapping, v):
    v_new = OrderedDict()

    for r_id in model.reactions:
        if r_id in mapping:
            fwd_id, bwd_id = mapping[r_id]
            v_new[r_id] = v[fwd_id] - v[bwd_id]
        else:
            v_new[r_id] = v[r_id]

    return v_new

