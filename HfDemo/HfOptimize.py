import numpy as np
import matplotlib as plt
import mayavi.mlab as mv
import gto


MAX_HF_ITERS = 10



class Atom:
    def __init__(this, position, charge):
        this.position = position
        this.charge = charge
        
def make_initial_guess(atoms):
    return []
    
def build_fock_matrix(atoms, current_wavefunctions):
    return []
    
def get_fock_eigenfunctions(matrix):
    return []
    
def hf_converged(w1, w2):
    return True

def hf_optimize(atoms):
    new_wavefunctions = make_initial_guess(atoms)
    
    for i in range(1, MAX_HF_ITERS):
        current_wavefunctions = new_wavefunctions
        fock_matrix = build_fock_matrix(atoms, current_wavefunctions)
        new_wavefunctions = get_fock_eigenfunctions(fock_matrix)
        
        if hf_converged(current_wavefunctions, new_wavefunctions):
            break

    
    
def gto_val(alpha, center, momentum, pos):
    delta = pos - center
    return (delta[0]**momentum[0])*(delta[1]**momentum[1])*(delta[2]**momentum[2])*np.exp(-alpha*np.dot(delta, delta))
    

X, Y, Z = np.mgrid[-2:2:100j, -2:2:20j, -2:2:20j]

c0 = np.array([0., 0., 0.])
m0 = np.array([1., 0., 0.])
m1 = np.array([0., 1., 0.])

cfunc = lambda x, y, z: (gto_val(1, c0, m0, np.array([x, y, z])) + gto_val(1, c0, m1, np.array([x, y, z])))

mv.clf()
mv.contour3d(X, Y, Z, np.vectorize(cfunc), contours=4)

