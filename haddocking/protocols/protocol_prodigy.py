# **************************************************************************
# *
# * Authors:     Blanca Pueche (blanca.pueche@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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

import json
import os
import re
import shlex

from pyworkflow.protocol import params
from pwchem.protocols import ProtDefineStructROIs

from pwchem.utils import pdbFromASFile, getBaseName
from haddocking.constants import PRODIGY_DIC

from haddocking import Plugin

from Bio.PDB import PDBParser, PDBIO, Select
from pwchem.objects import SetOfStructROIs, StructROI
from pwchem.utils import cifFromASFile


class ProtPRODIGYPrediction(ProtDefineStructROIs):
  """Run protein-protein interface binding affinity predictions using PRODIGY software from Bonvin Lab"""
  _label = 'PRODIGY protein interface binding affinity prediction'

  def _defineParams(self, form):
    form.addSection(label='Input')
    iGroup = form.addGroup('Input')
    iGroup.addParam('inputAtomStruct', params.PointerParam, pointerClass="AtomStruct", label='Input atom struct: ',
                    help="Input AtomStruct to predict the interface on")

    iGroup.addParam('interaction_mode', params.EnumParam, choices=['All intermolecular interactions', 'Specific molecular groups'],
                    default=0, label='Interaction analysis: ', help='Choose whether to evaluate all intermolecular interactions '
                        'or specify the molecular groups to compare.')

    iGroup.addParam('chain_name1', params.StringParam, allowsNull=True, label='Chains of interest Group 1: ', condition='interaction_mode==1',
                    help='Specify the chains of the structure of interest to form group 1.')
    iGroup.addParam('chain_name2', params.StringParam, allowsNull=True, label='Chains of interest Group 2: ', condition='interaction_mode==1',
                    help='Specify the chains of the structure of interest to form group 2.')

    iGroup.addParam('dist_cutoff', params.FloatParam, label='Distance cutoff: ', default=0,
                    expertLevel=params.LEVEL_ADVANCED,
                    help='Distance cutoff to calculate ICs. Note: you can leave this as is to use the default.')
    iGroup.addParam('acc_threshold', params.FloatParam, label='Accessibility threshold: ', default=0,
                    expertLevel=params.LEVEL_ADVANCED,
                    help='Accessibility threshold for BSA analysis. This parameter defines the minimum surface area accessibility required to consider a region as part of the binding interface during analysis. It is used to refine the calculation of buried surface area (BSA), which is critical for evaluating interaction interfaces. Note: you can leave this as is to use the default.')

    iGroup.addParam('temperature', params.FloatParam, label='Temperature: ', default=25,
                    help='Temperature (C) for Kd prediction.')

    form.addParallelSection(threads=4, mpi=1)

  def _insertAllSteps(self):
      self._insertFunctionStep('predictStep')
      self._insertFunctionStep('createOutputStep')

  def predictStep(self):
    # Convert the input structure to PDB format
    inpFile = self.getInputPath()

    pdbFile = pdbFromASFile(
      inpFile,
      self._getPdbFile(),
      atomStruct=self.inputAtomStruct.get()
    )

    args = shlex.quote(pdbFile)

    # Define molecular groups if specific interaction analysis is selected
    if self.interaction_mode.get() == 1:
      group1 = self.parseChainSelection(self.chain_name1.get())
      group2 = self.parseChainSelection(self.chain_name2.get())

      if not group1 or not group2:
        raise ValueError(
          'Both molecular groups must contain at least one chain.'
        )

      if set(group1).intersection(set(group2)):
        raise ValueError(
          'The same chain cannot belong to both molecular groups.'
        )

      selection = '{} {}'.format(
        shlex.quote(','.join(group1)),
        shlex.quote(','.join(group2))
      )

      args += ' --selection {}'.format(selection)

    # Prediction temperature
    args += ' --temperature {}'.format(self.temperature.get())

    # Optional parameters: omit flags to use PRODIGY defaults
    if self.dist_cutoff.get() > 0:
      args += ' --distance-cutoff {}'.format(self.dist_cutoff.get())

    if self.acc_threshold.get() > 0:
      args += ' --acc-threshold {}'.format(self.acc_threshold.get())

    # Always generate the contact list and PyMOL selection
    args += ' --contact_list --pymol_selection'

    args += ' -np {}'.format(str(self.numberOfThreads.get()))

    Plugin.runCondaCommand(
      self,
      args=args,
      condaDic=PRODIGY_DIC,
      program='prodigy',
      cwd=os.path.abspath(Plugin.getVar(PRODIGY_DIC['home']))
    )

  def createOutputStep(self):

    # Find PRODIGY contact file
    contactFiles = [
        os.path.join(self._getExtraPath(), f)
        for f in os.listdir(self._getExtraPath())
        if f.endswith('.ic')
    ]

    if not contactFiles:
        raise FileNotFoundError(
            'No PRODIGY contact (.ic) file found in the extra directory.'
        )

    if len(contactFiles) > 1:
        raise ValueError(
            'Multiple PRODIGY contact files found: {}'.format(contactFiles)
        )

    contactFile = contactFiles[0]

    # Parse contacts and collect interacting residues by chain
    chainResidues = {}

    with open(contactFile) as f:
        for line in f:
            fields = line.split()

            if len(fields) < 6:
                continue

            try:
                resName1, resId1, chain1 = fields[:3]
                resName2, resId2, chain2 = fields[3:6]

                chainResidues.setdefault(chain1, set()).add(
                    (resName1, int(resId1))
                )

                chainResidues.setdefault(chain2, set()).add(
                    (resName2, int(resId2))
                )

            except ValueError:
                continue

    if not chainResidues:
        raise ValueError(
            'No valid residue contacts found in PRODIGY contact file.'
        )

    # Convert original structure to CIF
    pdbFile = self._getPdbFile()
    cifFile = self._getExtraPath('complex.cif')

    cifFromASFile(pdbFile, cifFile)

    # Parse original PDB structure
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('complex', pdbFile)

    # Create output SetOfStructROIs
    outROIs = SetOfStructROIs(
        filename=self._getPath('StructROIs.sqlite')
    )

    # Generate one ROI for each interacting chain
    for i, (chainId, residues) in enumerate(chainResidues.items(), 1):

        roiPdbFile = self._getExtraPath(
            'interface_chain_{}.pdb'.format(chainId)
        )

        roiCifFile = self._getExtraPath(
            'interface_chain_{}.cif'.format(chainId)
        )

        class ContactResidueSelect(Select):

            def __init__(self, chainId, residues):
                self.chainId = chainId
                self.residues = residues

            def accept_chain(self, chain):
                return chain.id == self.chainId

            def accept_residue(self, residue):
                resName = residue.resname.strip()
                resId = residue.id[1]

                return (resName, resId) in self.residues

        # Extract contacting residues from the original structure
        io = PDBIO()
        io.set_structure(structure)

        io.save(
            roiPdbFile,
            ContactResidueSelect(chainId, residues)
        )

        # Convert ROI to CIF
        cifFromASFile(roiPdbFile, roiCifFile)

        # Create StructROI
        pocket = StructROI(roiCifFile, cifFile)
        pocket.calculateContacts()

        outROIs.append(pocket)

    # Define outputs
    if len(outROIs) > 0:
        outROIs.buildPDBhetatmFile()
        self._defineOutputs(outputStructROIs=outROIs)

  def parseChainSelection(self, chainInfo):
    """Extract chain identifiers from the Scipion chain-selection wizard."""

    if not chainInfo:
      return []

    try:
      data = json.loads(chainInfo)
    except (ValueError, TypeError):
      raise ValueError(
        'Invalid chain selection. Please select the chains using '
        'the chain-selection wizard.'
      )

    chains = []

    # Multiple chain selection: {"model-chain": "0-A, 0-B"}
    if 'model-chain' in data:
      for item in data['model-chain'].split(','):
        item = item.strip()

        if '-' not in item:
          raise ValueError(
            'Invalid chain selection format: {}'.format(item)
          )

        model, chain = item.split('-', 1)
        chains.append(chain.strip())

    # Single chain selection: {"model": 0, "chain": "A"}
    elif 'chain' in data:
      chains.append(str(data['chain']).strip())

    else:
      raise ValueError(
        'Could not identify the selected chains. '
        'Please use the chain-selection wizard.'
      )

    # Remove duplicated chain identifiers while preserving order
    return list(dict.fromkeys(chains))

  def getInputPath(self):
    return self.inputAtomStruct.get().getFileName()

  def _getPdbFile(self):
    return os.path.abspath(
      self._getExtraPath(self._getInputName() + '.pdb')
    )

  def _getInputName(self):
    return getBaseName(self.getInputPath())

  def _validate(self):
    errors = []

    if self.interaction_mode.get() == 1:
      if not self.chain_name1.get():
        errors.append(
          'Please select the chains for molecular group 1.'
        )

      if not self.chain_name2.get():
        errors.append(
          'Please select the chains for molecular group 2.'
        )

      if self.chain_name1.get() and self.chain_name2.get():
        try:
          group1 = self.parseChainSelection(
            self.chain_name1.get()
          )
          group2 = self.parseChainSelection(
            self.chain_name2.get()
          )

          if set(group1).intersection(set(group2)):
            errors.append(
              'The same chain cannot belong to both '
              'molecular groups.'
            )

        except ValueError as e:
          errors.append(str(e))

    if self.temperature.get() < -273.15:
      errors.append(
        'Temperature cannot be below absolute zero (-273.15 °C).'
      )

    if self.dist_cutoff.get() < 0:
      errors.append(
        'Distance cutoff must be greater than or equal to 0.'
      )

    if self.acc_threshold.get() < 0:
      errors.append(
        'Accessibility threshold must be greater than or equal to 0.'
      )

    return errors

  def _summary(self):
    summary = []
    logFile = self._getLogsPath('run.stdout')
    if not os.path.exists(logFile):
        return ['PRODIGY output log not found.']
    with open(logFile) as f:
        log = f.read()
    # Number of intermolecular contacts
    contacts = re.search(
        r'No\. of intermolecular contacts:\s*(\d+)',
        log
    )
    # Predicted binding affinity
    affinity = re.search(
        r'Predicted binding affinity \(kcal\.mol-1\):\s*([-+]?\d*\.?\d+)',
        log
    )
    # Predicted dissociation constant
    kd = re.search(
        r'Predicted dissociation constant \(M\).*?:\s*([0-9.eE+-]+)',
        log
    )
    if contacts:
        summary.append(
            'Number of intermolecular contacts: {}'.format(contacts.group(1))
        )
    if affinity:
        summary.append(
            'Predicted binding affinity (ΔG): {} kcal/mol'.format(
                affinity.group(1)
            )
        )
    if kd:
        summary.append(
            'Predicted dissociation constant (Kd): {} M'.format(
                kd.group(1)
            )
        )
    if not summary:
        summary.append('No PRODIGY prediction results found in the log.')

    return summary