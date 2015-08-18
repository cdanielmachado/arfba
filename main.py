from scipy import randn
from .sbml import load_allosteric_model
from .simulation import arFBA

def main():
	model = load_allosteric_model('arFBA/model.xml')

	reference_condition = {'R_EX_glc_e': (-2.93, -2.93),
	                       'R_Biomass_Ecoli_core_w_GAM': (0.2, 0.2)}

	perturbed_condition = {'R_EX_glc_e': (-6-63, -6.63),
	                       'R_Biomass_Ecoli_core_w_GAM': (0.5, 0.5)}

	weights = {key: 10**randn() for key in model.regulation.keys()}

	v, sol = arFBA(model, reference_condition, perturbed_condition, weights)

	for r_id, val in v.items():
		if val:
			print r_id.ljust(12), val

if __name__ == "__main__":
	main()