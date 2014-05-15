# -*- coding: utf-8 -*-
# Released under BSD license

import numpy as np
import itertools
import scipy.special

def boys_function(x, n):
    return scipy.special.gamma(n + 0.5) * scipy.special.gammainc(n, x) / (2 * x**(n + 0.5))

class SimpleGTO:
    
    def __init__(this, center, momentum, falloff):
        this.center = np.array(center)
        this.momentum = np.array(momentum)
        this.falloff = falloff

    # Evaluate Hermitian expansion
    @staticmethod
    def get_Sij(i, j, total_falloff, cog, xa, xb):
        
        if i < 0 or j < 0:
            return 0

        if i == 0 and j == 0:
            return np.sqrt(np.pi/total_falloff)

        if i <= j:
            return (cog - xb)*(SimpleGTO.get_Sij(i, j-1, total_falloff, cog, xa, xb)) + \
                   (i*SimpleGTO.get_Sij(i-1, j-1, total_falloff, cog, xa, xb)) + \
                    (j - 1)*(SimpleGTO.get_Sij(i, j-2, total_falloff, cog, xa, xb)) / (2*total_falloff)
                   
        if j < i:
            return (cog - xa)*SimpleGTO.get_Sij(i-1, j, total_falloff, cog, xa, xb) + \
                   ((i - 1)*SimpleGTO.get_Sij(i-2, j, total_falloff, cog, xa, xb) + \
                    j*SimpleGTO.get_Sij(i-1, j-1, total_falloff, cog, xa, xb)) / (2*total_falloff)
    

    @staticmethod
    def get_Dij(e, i, j, total_falloff, cog, xa, xb):
        return SimpleGTO.get_Dij_rec(e, i, j, total_falloff, cog, xa, xb, i, j)        
        
    @staticmethod
    def get_Dij_rec(e, i, j, total_falloff, cog, xa, xb, a, b):

        if i < 0 or j < 0:
            return 0
            
        if (i == 0 and j == 0) and e != 0:
            return 0
            return (SimpleGTO.get_Dij_rec(e, i, j+1, total_falloff, cog, xa, xb, a, b) - SimpleGTO.get_Dij_rec(e, i+1, j, total_falloff, cog, xa, xb, a, b) - \
                    e*SimpleGTO.get_Dij_rec(e-1, i, j, total_falloff, cog, xa, xb, a, b)) / (xa - xb)

        if e == 0:
            return SimpleGTO.get_Sij(i, j, total_falloff, cog, xa, xb)
            
        if i <= j:
            return (cog - xb) * SimpleGTO.get_Dij_rec(e, i, j-1, total_falloff, cog, xa, xb, a, b) + \
                   (i*SimpleGTO.get_Dij_rec(e, i-1, j-1, total_falloff, cog, xa, xb, a, b) + \
                    (j - 1)*SimpleGTO.get_Dij_rec(e, i, j-2, total_falloff, cog, xa, xb, a, b) + \
                    2*a*e*SimpleGTO.get_Dij_rec(e-1, i, j-1, total_falloff, cog, xa, xb, a, b)) / (2*total_falloff)
                    
        if j < i:
            return (cog - xa) * SimpleGTO.get_Dij_rec(e, i-1, j, total_falloff, cog, xa, xb, a, b) + \
                   ((i - 1)*SimpleGTO.get_Dij_rec(e, i-2, j, total_falloff, cog, xa, xb, a, b) + \
                    j*SimpleGTO.get_Dij_rec(e, i-1, j-1, total_falloff, cog, xa, xb, a, b) - \
                    2*b*e*SimpleGTO.get_Dij_rec(e-1, i-1, j, total_falloff, cog, xa, xb, a, b)) / (2*total_falloff)
                   
            
    def inner_product(this, g2):
        total_falloff = this.falloff + g2.falloff
        mu = this.falloff*g2.falloff / (this.falloff + g2.falloff)
        cog = (this.falloff*this.center + g2.falloff*g2.center) / total_falloff
        Xab = this.falloff*this.center - g2.falloff*g2.center
        
        integral = (np.pi/total_falloff)**(3./2.)*np.exp(-mu*np.dot(Xab, Xab))
        
        for dim in range(0, 3):
            integral *= SimpleGTO.get_Sij(this.momentum[dim], g2.momentum[dim], total_falloff, cog[dim], this.center[dim], g2.center[dim])
            
        return integral
        
    def kinetic(this, g2):
        total_falloff = this.falloff + g2.falloff
        mu = this.falloff*g2.falloff / (this.falloff + g2.falloff)
        cog = (this.falloff*this.center + g2.falloff*g2.center) / total_falloff
        Xab = this.falloff*this.center - g2.falloff*g2.center

        T = np.zeros(3)
        S = np.zeros(3)
        for dim in range(0, 3):
            T[dim] = (np.pi/total_falloff)**(3./2.)*np.exp(-mu*np.dot(Xab, Xab)) * (\
                -2*g2.falloff**2.0*SimpleGTO.get_Sij(this.momentum[dim], g2.momentum[dim]+2, total_falloff, cog[dim], this.center[dim], g2.center[dim]) + \
                g2.falloff*(2*g2.momentum[dim] + 1)*SimpleGTO.get_Sij(this.momentum[dim], g2.momentum[dim], total_falloff, cog[dim], this.center[dim], g2.center[dim]) + \
                -1/2*g2.momentum[dim]*(g2.momentum[dim] - 1)*SimpleGTO.get_Sij(this.momentum[dim], g2.momentum[dim]-2, total_falloff, cog[dim], this.center[dim], g2.center[dim]))

            S[dim] = SimpleGTO.get_Sij(this.momentum[dim], g2.momentum[dim], total_falloff, cog[dim], this.center[dim], g2.center[dim])


        return (T[0]*S[1]*S[2] + S[0]*T[1]*S[2] + S[0]*S[1]*S[2])*(np.pi/total_falloff)**(3./2.)*np.exp(-mu*np.dot(Xab, Xab))
        
    def nuclear(this):
        return []
        
        
    # Calculate amplitude at point pos
    def value_at(this, pos):
        delta = np.array(pos) - this.center
        (2*this.falloff/np.pi)**(3./4.)*(delta[0]**this.momentum[0])*(delta[1]**this.momentum[1])*(delta[2]**this.momentum[2])*np.exp(-this.falloff*np.dot(delta, delta))
        
    def set_center(this, center):
        this.center = np.array(center)


class STO_nG:
    def __init__(this, center, momentum, falloffs, weights):
        this.gaussians = []
        
        for orb in zip(falloffs, weights):
            this.gaussians.append((SimpleGTO(center, momentum, orb[0]), orb[1]))

            
    def inner_product(this, g2):
        combinations = itertools.product(this.gaussians, g2.gaussians)
        
        result = 0
        
        for c in combinations:
            result += c[0][0].inner_product(c[1][0]) * c[0][1] * c[1][1]
            
        return result
        
    def kinetic(this, g2):
        combinations = itertools.product(this.gaussians, g2.gaussians)
        
        result = 0
        
        for c in combinations:
            result += c[0][0].kinetic(c[1][0]) * c[0][1] * c[1][1]

        return result


    def nuclear(this, nucleus_pos):
        result = 0
        
        for gs in this.gaussians:
            result += gs[0].nuclear(nucleus_pos) * gs[1]
            
        return result

    def coulomb(this, g2, g3, g4):
        combinations = itertools.product(this.gaussians, g2.gaussians, g3.gaussians, g4.gaussians)
        
        result = 0
        
        for c in combinations:
            result += c[0][0].coulomb(c[1][0], c[2][0], c[3][0])*c[0][1]*c[1][1]*c[2][1]*c[3][1]
            
        return result

    def value_at(this, pos):
        val = 0
        
        for orb in zip(this.gaussians, this.weights):
            val += orb[0].value_at(pos) * orb[1]
            
        return val

    def set_center(this, center):
        for g in this.gaussians:
            g[0].center = center
            

minimal_basis = [STO_nG([0., 0., 0.], [0., 0., 0.], [0.1516230, 0.851819], [0.678914, 0.4301290]), # 1s
                 STO_nG([0., 0., 0.], [0., 0., 0.], [0.0974545, 0.384244], [0.963782, 0.0494718]), # 2s
                 STO_nG([0., 0., 0.], [1., 0., 0.], [0.0974545, 0.384244], [0.612820, 0.5115410]), # 2p
                 STO_nG([0., 0., 0.], [0., 1., 0.], [0.0974545, 0.384244], [0.612820, 0.5115410]), # 2p
                 STO_nG([0., 0., 0.], [0., 0., 1.], [0.0974545, 0.384244], [0.612820, 0.5115410])] # 2p

# Calculate <g1|g2>
def inner_product(g1, g2):
    return g1.inner_product(g2)

# Calculate kinetic energy integral <g1|laplace|g2>
def kinetic(g1, g2):
    return g1.kinetic(g2)
    
# Calculate electron-nucleus interaction integral <g|1/r|g>
def nuclear(g1, g2, nucleus_pos):
    return g1.nuclear(nucleus_pos, g2)
    
#calculate electron-electron Coulomb integral <g1,g2|1/r|g1,g2>
def coulomb(g1, g2, g3, g4):
    return g1.coulomb(g2, g3, g4)
    
#calculate electron-electron exchange integral <g1,g2|1/r|g2,g1>
def exchange(g1, g2, g3, g4):
    return g1.coulomb(g3, g2, g4)
    

a = SimpleGTO([0., 1., 0.], [3., 0., 0.], 1.)
b = SimpleGTO([0., 1., 0.], [1., 0., 0.], 2.)
print a.inner_product(b)
print a.kinetic(a)
