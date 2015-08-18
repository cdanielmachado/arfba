
from framed.analysis.simulation import pFBA
from framed.solvers import solver_instance
from framed.solvers.solver import Status, VarType
from .utils import build_perturbed_model, compute_turnover, merge_fluxes


def arFBA(model, reference_constraints, perturbed_constraints, weights):

    TOL = 1e-6
    M = 1e3

    solution_ref = pFBA(model, constraints=reference_constraints)
    v0 = solution_ref.values
    t0 = compute_turnover(model, v0)

    model_irrev, mapping = build_perturbed_model(model, perturbed_constraints)
    solver = solver_instance()
    solver.build_problem(model_irrev)

    m_r_lookup = model_irrev.metabolite_reaction_lookup_table()
    reg_m_r_lookup = model_irrev.metabolite_reaction_regulatory_lookup_table()
    regulators = [m_id for m_id, targets in reg_m_r_lookup.items() if len(targets) > 0]

    for m_id in regulators:
        solver.add_variable('t_' + m_id, 0, None, persistent=False, update_problem=False)

    for fwd_id, bwd_id in mapping.values():
        solver.add_variable('y_' + fwd_id, vartype=VarType.BINARY, persistent=False, update_problem=False)
        solver.add_variable('y_' + bwd_id, vartype=VarType.BINARY, persistent=False, update_problem=False)

    for (m_id, r_id), kind in model.regulation.items():
        if v0[r_id] > TOL and t0[m_id] > TOL:
            diff_pos = 'd+_{}_{}'.format(m_id, r_id)
            diff_neg = 'd-_{}_{}'.format(m_id, r_id)
            solver.add_variable(diff_pos, 0, None, persistent=False, update_problem=False)
            solver.add_variable(diff_neg, 0, None, persistent=False, update_problem=False)

    solver.update()

    for m_id in regulators:
        lhs = {r_id: coeff for r_id, coeff in m_r_lookup[m_id].items() if coeff > 0}
        lhs['t_' + m_id] = -1
        solver.add_constraint('ct_' + m_id, lhs.items(), persistent=False, update_problem=False)

    for r_id, (fwd_id, bwd_id) in mapping.items():
        solver.add_constraint('c_' + fwd_id, {fwd_id: 1, 'y_' + fwd_id: -M}.items(), '<', 0, persistent=False, update_problem=False)
        solver.add_constraint('c_' + bwd_id, {bwd_id: 1, 'y_' + bwd_id: -M}.items(), '<', 0, persistent=False, update_problem=False)
        solver.add_constraint('rev_' + r_id, {'y_' + fwd_id: 1, 'y_' + bwd_id: 1}.items(), '<', 1, persistent=False, update_problem=False)

    for (m_id, r_id), kind in model.regulation.items():
        if v0[r_id] > TOL and t0[m_id] > TOL:
            diff_pos = 'd+_{}_{}'.format(m_id, r_id)
            diff_neg = 'd-_{}_{}'.format(m_id, r_id)
            if r_id in mapping:
                fwd_id, bwd_id = mapping[r_id]
                if kind == 1:
                    lhs = {diff_pos: 1, fwd_id: -1/v0[r_id], bwd_id: -1/v0[r_id], 't_' + m_id: 1/t0[m_id]}
                    solver.add_constraint('c' + diff_pos, lhs.items(), '>', 0, persistent=False, update_problem=False)
                    lhs = {diff_neg: 1, fwd_id: 1/v0[r_id], bwd_id: 1/v0[r_id], 't_' + m_id: -1/t0[m_id]}
                    solver.add_constraint('c' + diff_neg, lhs.items(), '>', 0, persistent=False, update_problem=False)
                else:
                    lhs = {diff_pos: 1, fwd_id: -1/v0[r_id], bwd_id: -1/v0[r_id], 't_' + m_id: -1/t0[m_id]}
                    solver.add_constraint('c' + diff_pos, lhs.items(), '>', -2, persistent=False, update_problem=False)
                    lhs = {diff_neg: 1, fwd_id: 1/v0[r_id], bwd_id: 1/v0[r_id], 't_' + m_id: 1/t0[m_id]}
                    solver.add_constraint('c' + diff_neg, lhs.items(), '>', 2, persistent=False, update_problem=False)
            else:
                if kind == 1:
                    lhs = {diff_pos: 1, r_id: -1/v0[r_id], 't_' + m_id: 1/t0[m_id]}
                    solver.add_constraint('c' + diff_pos, lhs.items(), '>', 0, persistent=False, update_problem=False)
                    lhs = {diff_neg: 1, r_id: 1/v0[r_id], 't_' + m_id: -1/t0[m_id]}
                    solver.add_constraint('c' + diff_neg, lhs.items(), '>', 0, persistent=False, update_problem=False)
                else:
                    lhs = {diff_pos: 1, r_id: -1/v0[r_id], 't_' + m_id: -1/t0[m_id]}
                    solver.add_constraint('c' + diff_pos, lhs.items(), '>', -2, persistent=False, update_problem=False)
                    lhs = {diff_neg: 1, r_id: 1/v0[r_id], 't_' + m_id: 1/t0[m_id]}
                    solver.add_constraint('c' + diff_neg, lhs.items(), '>', 2, persistent=False, update_problem=False)

    solver.update()

    objective = {r_id: -1 for r_id in model_irrev.reactions}

    for (m_id, r_id) in model.regulation.keys():
        if (m_id, r_id) in weights and v0[r_id] > TOL and t0[m_id] > TOL:
            diff_pos = 'd+_{}_{}'.format(m_id, r_id)
            diff_neg = 'd-_{}_{}'.format(m_id, r_id)
            objective[diff_pos] = -weights[(m_id, r_id)]
            objective[diff_neg] = -weights[(m_id, r_id)]

    solution = solver.solve_lp(objective)

    if solution.status == Status.OPTIMAL:
        v = merge_fluxes(model, mapping, solution.values)
    else:
        v = None

    return v, solution