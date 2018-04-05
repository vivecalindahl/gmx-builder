#! /usr/bin/env python2.7

import os, sys
import numpy as np

import gmx_builder as gmxb
from gmx_builder import run_in_shell as xsh

# We should be able to set the gmx binary, e.g. from an environment variable.
# It should be set somewhere else than here though. Now we just assume gmx is in the path.
gmx='gmx'

def make_box_for_periodic_dna(gro='conf.gro'):
    # Make a box for a DNA molecule that is periodically connected in the z-direction.

    # The box angles (naming as in gmx manual).
    # Here vectors a and b lie in the xy-plane, and are 60 deg apart. 
    # vector c is parallel to the z axis
    bc, ac, ab = 90, 90, 60

    # Very simple measure of the diameter along each dimension:
    coords = gmxb.read_gro_coordinates(gro)
    diameter = coords.max(0) - coords.min(0)

    # Dna is roughly circular in x and y. Use the x diameter and add 2 nm
    # to ensure that periodic images are ~ 2 nm apart.

    a = diameter[0] + 2.0
    b = a

    # The z box length should equal the molecule length since we are making a
    # periodic molecule.
    c = diameter[-1]

    # Make the box with editconf
    angles = ' '.join([str(bc), str(ac), str(ab)])
    lengths =  ' '.join([str(a), str(b), str(c)])
    out='conf.gro'
    args = ' '.join(['-angles', angles, '-box', lengths, '-o', out])
    stdout = xsh(' '.join([gmx, 'editconf', args]))

def pdb2gmx_periodic(pdb, watermodel, forcefield):

    # Assume a pdb2gmx wrapper script for periodic molecules is in the same directory as this file.
    periodic_pdb2gmx = os.path.dirname(os.path.realpath(__file__)) + '/gmx-pdb2gmx-wrapper-periodic-dna.sh'

    args=' '.join(['-water', watermodel, '-ff', forcefield, '-f', pdb])
    xsh(' '.join([periodic_pdb2gmx, args]))

def solvate_box(gro='conf.gro', top='topol.top', tpr='topol.tpr'):
    args = ' '.join(['-cp', gro, '-p', top, '-cs','-o', 'conf.gro'])
    stdout = xsh(' '.join([gmx, 'solvate', args]))

def neutralize_box(gro='conf.gro', top='topol.top', tpr='topol.tpr'):
    args = ' '.join(['-neutral', '-pname', 'NA', '-p', 'topol.top', '-o', 'conf.gro'])
    interactive_args = 'SOL'
    cmd = ' '.join(['echo -e', interactive_args, '|', gmx, 'genion', args])

    stdout = xsh(cmd)

def build_periodic_dna(pdb_path, watermodel='tip3p', forcefield='charmm27', make_clean=True):

    # Generate gmx topology and config file.
    pdb2gmx_periodic(pdb_path, watermodel, forcefield)

    # Test if a gmx run/tpr file can be generated.
    # Also, the gmx tools often require a tpr as input.
    gmxb.make_tpr('.', nomdp=True)

    # Set the box size and shape.
    make_box_for_periodic_dna()

    # Add water to the box.
    gmxb.make_tpr('.', nomdp=True)
    solvate_box()

    # Neutralize system by adding ions.
    gmxb.make_tpr('.', nomdp=True)
    neutralize_box()

    # Make a final tpr.
    gmxb.make_tpr('.', nomdp=True)

    # Clean up
    if make_clean:
        gmxb.remove_temporary_files()

# Parameter file (mdp) settings
em_mdp = {
    'integrator': 'steep',
    'emtol': '1000.0',
    'emstep': '0.01',
    'nsteps': '50000',
    'coulombtype': 'pme',
    'vdwtype': 'cut-off',
    'cut-off-scheme': 'verlet',
}
        
npt_mdp = {
    'integrator': 'md',
    'dt': '0.002',
    'nsteps': '25000000',
    'nstlog': '5000000',
    'nstenergy': '50000',
    'nstxtcout': '500000',
    # pressure
    'pcoupl': 'parrinello-rahman',
    'tau-p': '5.0',
    'pcoupl-type': 'isotropic',
    'ref-p': '1.0',
    'compressibility': '4.5e-5',
    # temperature
    'tcoupl': 'v-rescale',
    'tau-t': '0.5',
    'ref-t': '300',
    'tc-grps': 'system',
    'gen-vel': 'yes',
    'gen-temp': '300',
    # electrostatics and vdw
    'coulombtype': 'pme',
    'vdwtype':'cut-off',
    'cut-off-scheme': 'verlet',
    # other
    'constraints': 'h-bonds'
}

def pull_basepair_distance_mdp(base_name1, base_name2):
    mdp = {
        'pull                     ': 'yes',
        'pull-print-ref-value     ': 'yes',
        'pull-nstxout             ': '5000',
        'pull-nstfout             ': '0',
        'pull-ngroups             ': '2',
        'pull-ncoords             ': '1',
        'pull-group1-name         ': base_name1,
        'pull-group2-name         ': base_name2,
        'pull-coord1-groups       ': '1 2',
        'pull-coord1-geometry    ' : 'distance'
    }
    return mdp

def awh_basepair_dist_mdp(name1, name2):

    pull_mdp = pull_basepair_distance_mdp(name1, name2)

    awh_mdp ={
        'pull-coord1-potential-provider':'awh',
        'pull-coord1-type        ': 'external-potential',
        'awh                      ': 'yes',
        'awh-nbias                ': '1',
        'awh-nstout               ': '50000',
        'awh-share-multisim       ': 'yes',
        'awh1-share-group         ': '1',
        'awh1-ndim                ': '1',
        'awh1-error-init          ': '5',
        'awh1-dim1-coord-index    ': '1',
        'awh1-dim1-diffusion      ': '5e-5',
        'awh1-dim1-start          ': '0.25',
        'awh1-dim1-end            ': '0.60',
        'awh1-dim1-force-constant ': '128000',
        'awh1-equilibrate-histogram' : 'yes'
    }
    
    return gmxb.merge_mdps([pull_mdp, awh_mdp])

def electrostatics_vdw_mdp(ff_name):
    # These are just examples, but should at least be reasonable.
    if ff_name == 'charmm':
        mdp = {
            'rcoulomb':'1.2',
            'fourierspacing': '0.14',
            'vdw-modifier': 'force-switch',
            'rvdw-switch': '0.8',
            'rvdw': '1.2'
        }
    elif ff_name == 'amber':
         mdp = {
       'rcoulomb':'1.0',
        'fourierspacing': '0.121',
        'vdw-modifier': 'potential-shift',
        'rvdw': '1.0',
        'dispcorr': 'ener-pres'
         }
    else:
        valid = {'charmm', 'amber'}
        raise ValueError("results: status must be one of %r." % valid)

    return mdp
            
def mdp_periodic_dna(ff_name, run_type):

    if run_type == 'npt':
        base = npt_mdp
    elif run_type == 'em':
        base = em_mdp
    else:    
        valid = {'npt', 'em'}
        raise ValueError("results: status must be one of %r." % valid)

    # Force field specific settings
    electrostatics_vdw = electrostatics_vdw_mdp(ff_name)

    # Extra definitions to add because of the molecule being periodi
    periodic_mol = {
        'periodic-molecules': 'yes',
        'pcoupl-type': 'semiisotropic',
        'ref-p': '1.0 1.0',
        'compressibility': '4.5e-5 4.5e-5',
    }

    # The merge order matters. Later mdps have precedence.
    mdp =  gmxb.merge_mdps([base, electrostatics_vdw, periodic_mol])

    return mdp

# Generate all the distance selections for a base pair in a DNA double helix, 
# i.e. two paired DNA chains of equal lengths N.
def basepair_distance_selections(resid1, resid2, name1='base_N1orN3', name2='partner_N1orN3'):

    # Select the two Watson-Crick hydrogen bonding N1 and N3 atom in the base pair.
    sel = ['base = resid '+ str(resid1),
           'partner = resid ' + str(resid2),
           name1 + ' = base and ((resname DA DG and name N1) or (resname DT DC and name N3))',
           name2 + ' = partner and ((resname DA DG and name N1) or (resname DT DC and name N3))',
           name1, name2]
           
    return sel

def get_basepair_resids(gro):
    # Figure out the length of the DNA chain (could alternatively give the sequence/length as input).
    # Maybe functionalities from e.g. MDAnalysis should be applied here instead...
    table = gmxb.read_gro_table(gro)
    col_resname = gmxb.gro_table_column('residue_name')
    col_atomname = gmxb.gro_table_column('atom_name')
    col_resid = gmxb.gro_table_column('resid')

    dna_table = filter(lambda row: row[col_resname] in {'DC', 'DT', 'DA', 'DG'}, table)

    # The number of base pairs equals the number of DNA residues/2
    nbp = len(set([row[col_resid] for row in dna_table]))/2

    # Assume that the bases are indexed and paired as:
    # 1    --   2N
    # 2    --   2N-1
    #     [..]  
    # n    --   2N+1-n     
    #     [..]       
    # N-1  --   N+2
    # N    --   N+1
    assert dna_table[0][col_resid] == '1' and dna_table[-1][col_resid] == str(nbp*2)

    # Resid indexing starts at 1
    resid_pairs = []
    for n in range(1,nbp+1):
        resid_pairs.append((n, 2*nbp + 1 - n))

    return resid_pairs

# Select the target base pairs to calculate opening free energy for.
def get_target_basepair_resids(gro):
    basepair_resids = get_basepair_resids(gro)

    # Take 3 in the middle, e.g.
    ntarget = 3
    nbp = len(basepair_resids)
    return basepair_resids[nbp/2:nbp/2 + ntarget]

# An example of how how one could build a simulation experiment for the periodic DNA system.
def example_build(make_clean=False):

    # External parameters are assumed to be present
    scripts_dir = os.path.dirname(os.path.realpath(__file__))
    external_params_dir =  scripts_dir + '/../external-parameters'

    # The pdb files could be generated from a sequence from within here calling some modeling tool like x3dna.
    pdb_dir = external_params_dir + '/pdbs'
    pdbs = ['/'.join([pdb_dir,f]) for f in os.listdir(pdb_dir) if f.endswith('.pdb') ]

    # Build specifications (model parameters)
    build_list = [
        {'name':'charmm', 'ff': 'charmm27', 'water':'tip3p', 'ffdir':None},
        {'name':'amber', 'ff': 'amber99bsc1', 'water':'spce', 'ffdir': external_params_dir + '/forcefields/amber99bsc1.ff'}
    ]

    def sysname(pdb):
        return pdb.split('.pdb')[0].split('/')[-1]

    startdir=os.getcwd()

    # Here, each pdb is built with each model.
    for pdb, specs in zip(len(build_list)*pdbs, len(pdbs)*build_list,):
        print 'Building:', sysname(pdb),  specs['name']

        # Define the directory hierarchy. Infer system name from pdb file.
        build_dir = '/'.join([startdir, specs['name'], sysname(pdb), 'build'])

        if make_clean:
            xsh('rm -rf ' + build_dir)
        xsh('mkdir -p ' + build_dir)

        # Add external parameters if needed
        if specs['ffdir']:
            xsh('cp -r ' + specs['ffdir'] + ' ' + build_dir)

        # Build the gmx system
        os.chdir(build_dir)
        build_periodic_dna(pdb, forcefield=specs['ff'], watermodel=specs['water'])
        os.chdir(startdir)

        # Add runs.  Here, each run is added to each build.
        run_list = [
            {'name':'em', 'mdp': mdp_periodic_dna(specs['name'], 'em'), 'selections':[]},          
            {'name':'npt', 'mdp': mdp_periodic_dna(specs['name'], 'npt'), 'selections':[]},
        ]
        
        # Add AWH runs, calculate PMF for some base pairs (bp).  The reaction coordinate requires
        # bp-specific atom selections. Here, generate all possible selections, then specify which
        # bp is targeted in the mdp-file (could also generate only the selections needed for the
        # target bp and use one mdp file).
        target_basepairs = get_target_basepair_resids(build_dir + '/conf.gro')

        for resid1, resid2 in target_basepairs: 
           name1, name2 = ['resid' + str(resid) for resid in [resid1, resid2]]
           awh_sel = basepair_distance_selections(resid1, resid2, name1=name1, name2=name2)
           awh_mdp = gmxb.merge_mdps([
               mdp_periodic_dna(specs['name'], 'npt'),
               awh_basepair_dist_mdp(name1, name2),
               {'nsteps':'100000000'}
           ])
           
           # Add the run specs (mdp and selections) for this base pair
           run_list.append({'name':'-'.join(['awh', name1, name2]), 'mdp': awh_mdp, 'selections': awh_sel})

        # Put the run directory on the same level as the build director
        print 'Adding runs:'
        for run in run_list:
            print run['name']
            run_dir = '/'.join([startdir, specs['name'], sysname(pdb), run['name']])
            if make_clean:
                xsh('rm -rf ' + run_dir)
            xsh('mkdir -p ' + run_dir)
            template_dir = run_dir + '/template'

            gmxb.make_run_template(build_dir, run['mdp'], template_dir, selections=run['selections'])
