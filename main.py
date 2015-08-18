""" allosteric regulation FBA (arFBA)

@author: Daniel Machado

   Copyright 2015 Centre of Biological Engineering,
   University of Minho.

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

"""

from scipy import randn
from .sbml import load_allosteric_model
from .simulation import arFBA


def main():
	""" To run a default test simulation, run the following command from outside the source folder:

	> python -m arfba.main

	"""

	model = load_allosteric_model('arfba/model.xml')

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