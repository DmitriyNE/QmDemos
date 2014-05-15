# -*- coding: utf-8 -*-
# Released under BSD license

import gto
import copy

class Molecule:

    def __init__(this, nuclei_positions, charges):
        
        this.nuclei_positions = nuclei_positions
        this.nuclei_charges = charges
        this.nuclei_count = len(nuclei_positions)
        
        this.shell_count = sum(charges)
        this.basis = []
        
        for pos in this.nuclei_positions:
            for shell in gto.minimal_basis:
                s = copy.deepcopy(shell)
                s.set_center(pos)
                this.basis.append(s)
                
        this.basis_size = len(this.basis)

        
    def plot_density(this):
        []

    def plot_nmr_shielding(this):
        []

    def plot_orbital(this):
        []
