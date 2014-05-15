# -*- coding: utf-8 -*-
# Released under BSD license

import molecule
import gto
import numpy

class operator_cache:
    
    def __init__(this, molecule):
        this.molecule = molecule
        
        this.inner_product_cache = numpy.zeros(molecule.basis_size, molecule.basis_size) * numpy.nan
        this.kinetic_cache = numpy.zeros(molecule.basis_size, molecule.basis_size) * numpy.nan
        this.nuclear_cache = numpy.zeros(molecule.basis_size, molecule.basis_size, molecule.nuclei_count) * numpy.nan
        this.coulomb_cache = numpy.zeros(molecule.basis_size, molecule.basis_size, molecule.basis_size, molecule.basis_size) * numpy.nan
        
    def precalculate(this):
        
        cnt = this.molecule.shell_count
        
        for i in range(0, cnt):
            for j in range(0, cnt):

                this.get_inner_product(i, j)
                this.get_kinetic(i, j)
                
                for k in range(0, this.molecule.nuclei_count)
                    this.get_nuclear(i, j, k)

                for k in range(0, cnt):
                    for l in range(0, cnt):
                        this.get_coulomb(i, j, k, l)
        
    def get_inner_product(this, i, j):
        
        if i > j:
            i, j = j, i
            
        if numpy.isnan(this.inner_product_cache[i, j]):
            this.inner_product_cache[i, j] = gto.inner_product(this.molecule.basis[i], this.molecule.basis[j])
            
        return this.inner_product_cache[i, j]
            
    def get_kinetic(this, i, j):
        
        if i > j:
            i, j = j, i
            
        if numpy.isnan(this.kinetic_cache[i, j]):
            this.kinetic_cache[i, j] = gto.kinetic(this.molecule.basis[i], this.molecule.basis[j])
            
        return this.kinetic_cache[i, j]
            
    def get_nuclear(this, i, j, nuclei):
        
        if i > j:
            i, j = j, i
            
        if numpy.isnan(this.nuclear_cache[i, j, nuclei]):
            this.nuclear_cache[i, j] = gto.nuclear(this.molecule.basis[i], this.molecule.basis[j], this.molecule.nuclei_positions[nuclei])
            
        return this.nuclear_cache[i, j]
            
    def get_coulomb(this, i, j, k, l):
        
        if i > j:
            i, j = j, i
            
        
        if(numpy.isnan(this.coulomb_cache[i, j, k, l])):
            this.coulomb_cache[i, j, k, l] = gto.coulomb(this.molecule.basis[i], this.molecule.basis[j], this.molecule.basis[k], this.molecule.basis[l])
            
        return this.coulomb_cache[i, j, k, l]

class FockMatrixBuilder:
    
    def __init__(this, molecule):
        
        this.molecule = molecule
        this.opcache = operator_cache(molecule)
        this.prepare_inv_part()
        
    def prepare_inv_part(this):
        this.fock_cached = numpy.zeros(this.molecule.basis_size, this.molecule.basis_size)
        this.inner_product_matrix = numpy.zeros(this.molecule.basis_size, this.molecule.basis_size)
        
        for i in range(0, this.molecule.basis_size):
            for j in range(0, this.molecule.basis_size):
                this.inner_product_matrix = this.opcache.get_inner_product(i, j)
                this.fock_cached[i, j] = -this.opcache.get_kinetic(i, j)
                
                for k in range(0, this.molecule.nuclei_count):
                    this.fock_cached[i, j] += -this.opcache.get_nuclear(i, j, k)
        
    def build_matrix(this, state):
        ret = fock_cached
        
        for i in range(0, this.molecule.basis_size):
            for j in range(0, this.molecule.basis_size):        
                for k in range()
        