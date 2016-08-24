'''
Example comparing API for running CPD and POLIMI

use this comparison for defining common robust API
'''

from __future__ import division, absolute_import
from __future__ import print_function, unicode_literals

import pkp.polimi
import pkp.cpd
from test_cpd import ua, pa


operating_conditions = [[0, 400],
                        [0.01, 1300],
                        [0.02, 1300]]
p = 101325
# polimi
print('Run Polimi ...')
polimi = pkp.polimi.Polimi(ultimate_analysis=ua,
                           proximate_analysis=pa,
                           pressure=p,
                           name='Polimi test',
                           mechanism='COAL.xml')
polimi.operating_conditions = operating_conditions
polimi_res = polimi.run()

# CPD
print('Run CPD...')
cpd_fname = 'CPD_test'
cpd = pkp.cpd.CPD(ultimate_analysis=ua,
                  proximate_analysis=pa,
                  pressure=p,
                  name='CPD test')
cpd.path = './test'  # add path to set property
cpd.operating_conditions = operating_conditions
cpd.set_parameters(dt=1e-5, increment=2, dt_max=1e-5,
                   basename='test', nmr_parameters=None)
cpd_res = cpd.run()
