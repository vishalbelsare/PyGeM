import matplotlib
matplotlib.use('agg')

import nose

test_defaults = [
    'tests/test_freeform.py',
    'tests/test_idwparams.py',
    'tests/test_affine.py',
    'tests/test_idw.py',
    'tests/test_khandler.py',
    'tests/test_mdpahandler.py',
    'tests/test_openfhandler.py',
    'tests/test_package.py',
    'tests/test_radial.py',
    'tests/test_rbfparams.py',
    'tests/test_stlhandler.py',
    'tests/test_unvhandler.py',
    'tests/test_vtkhandler.py',
]


test_cad = [
    'tests/test_igeshandler.py',
    'tests/test_nurbshandler.py',
    'tests/test_stephandler.py',
]

default_argv = ['--tests'] + test_defaults
cad_argv = ['--tests'] + test_cad

try:
    import pygem.cad
    nose.run(argv=cad_argv)
except:
    print('module OCC not found, skip tests for CAD')

nose.run(argv=default_argv)





