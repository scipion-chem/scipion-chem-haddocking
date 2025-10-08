# coding: latin-1
# **************************************************************************
# *
# * Authors:  Daniel Del Hoyo Gomez (ddelhoyo@cnb.csic.es)
# *
# * Biocomputing Unit, CNB-CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

# Common constants
DEFAULT_VERSION = '1.0'

# Package dictionaries
HADDOCK_DIC = {'name': 'HADDOCK',    'version': '3.0', 'url': 'https://github.com/haddocking/haddock3.git',
               'home': 'HADDOCK_HOME', 'activation': 'HADDOCK_ACTIVATION_CMD'}

CNS_DIC = {'name': 'CNS',    'version': '1.3',
           'home': 'CNS_HOME', 'tar': 'CNS_TAR', 'pattern': 'cns'}

CPORT_DIC = {'name': 'CPORT',    'version': '0.2.0', 'url': 'https://github.com/haddocking/cport.git',
             'home': 'CPORT_HOME', 'activation': 'CPORT_ACTIVATION_CMD'}

GDOCK_DIC = {'name': 'GDock',    'version': '1.0', 'url': 'https://github.com/rvhonorato/gdock.git',
             'home': 'GDOCK_HOME'}

LIGHTDOCK_DIC = {'name': 'LightDock',    'version': '1.0', 'url': 'https://github.com/lightdock/lightdock.git',
                 'home': 'LIGHTDOCK_HOME'}


DEF_PROTOCOL = '''# ====================================================================
# Antibody-antigen docking example with restraints from the antibody
# paratope to the NMR-identified epitope on the antigen 
# ====================================================================

# Directory in which the scoring will be done
run_dir = "{}"

# Compute mode
mode = "local"
ncores = {}

# Self contained rundir (to avoid problems with long filename paths)
self_contained = true

# Post-processing to generate statistics and plots
postprocess = true
clean = false

molecules =  [
    "{}",
    "{}"
    ]

# ====================================================================
# Parameters for each stage are defined below, prefer full paths
# ====================================================================
[topoaa]

[rigidbody]
# CDR to NMR epitope ambig restraints
ambig_fname = "{}"
# Restraints to keep the antibody chains together
unambig_fname = "{}"
# Reduced sampling (100 instead of the default of 1000)
sampling = {}

[caprieval]

[seletop]
# Selection of the top 50 best scoring complexes (instead of the default of 200)
select = {}

[flexref]
tolerance = {}
# CDR to NMR epitope ambig restraints
ambig_fname = "{}"
# Restraints to keep the antibody chains together
unambig_fname = "{}"

[caprieval]

[emref]
tolerance = {}
# CDR to NMR epitope ambig restraints
ambig_fname = "{}"
# Restraints to keep the antibody chains together
unambig_fname = "{}"

[caprieval]

[clustfcc]
plot_matrix = true

[seletopclusts]
# Selection of the top 4 best scoring complexes from each cluster
top_models = {}

[caprieval]

[contactmap]

# ===================================================================='''
