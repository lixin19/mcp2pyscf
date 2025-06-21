from pyscf import gto
from pyscf import scf

mol = gto.M(atom='N 0 0 0; O 0 0 1.2', basis='6-311+gss', spin=1,charge=0)
mol.build()
mf = scf.ROHF(mol).run(conv_tol=.5)
mf = scf.newton(mf).set(conv_tol=1e-9)
mf.verbose = 5
mf.kernel()

