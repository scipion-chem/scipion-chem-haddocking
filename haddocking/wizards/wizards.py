# -*- coding: utf-8 -*-
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

"""
"""

from pwchem.wizards import SelectChainWizardQT, SelectResidueWizardQT, SelectElementWizard, AddElementWizard, \
  SelectAtomWizardQT
from ..protocols import ProtCPORTPrediction, ProtHaddockPrediction

SelectChainWizardQT().addTarget(protocol=ProtCPORTPrediction,
                                targets=['chain_name'],
                                inputs=['inputAtomStruct'],
                                outputs=['chain_name'])

SelectElementWizard().addTarget(protocol=ProtHaddockPrediction,
                                targets=['restrROI_1'],
                                inputs=['restrROIs_1'],
                                outputs=['restrROI_1'])
SelectElementWizard().addTarget(protocol=ProtHaddockPrediction,
                                targets=['restrROI_2'],
                                inputs=['restrROIs_2'],
                                outputs=['restrROI_2'])

# Ambiguous restraints
SelectChainWizardQT().addTarget(protocol=ProtHaddockPrediction,
                                targets=['chain_name_1'],
                                inputs=['inputAtomStruct_1'],
                                outputs=['chain_name_1'])
SelectChainWizardQT().addTarget(protocol=ProtHaddockPrediction,
                                targets=['chain_name_2'],
                                inputs=['inputAtomStruct_2'],
                                outputs=['chain_name_2'])

SelectResidueWizardQT().addTarget(protocol=ProtHaddockPrediction,
                                  targets=['restrRes_1'],
                                  inputs=['inputAtomStruct_1', 'chain_name_1'],
                                  outputs=['restrRes_1'])
SelectResidueWizardQT().addTarget(protocol=ProtHaddockPrediction,
                                  targets=['restrRes_2'],
                                  inputs=['inputAtomStruct_2', 'chain_name_2'],
                                  outputs=['restrRes_2'])

# Specific restraints
SelectChainWizardQT().addTarget(protocol=ProtHaddockPrediction,
                                targets=['chain_name_spec_1'],
                                inputs=['inputAtomStruct_1'],
                                outputs=['chain_name_spec_1'])
SelectChainWizardQT().addTarget(protocol=ProtHaddockPrediction,
                                targets=['chain_name_spec_2'],
                                inputs=['inputAtomStruct_2'],
                                outputs=['chain_name_spec_2'])
SelectResidueWizardQT().addTarget(protocol=ProtHaddockPrediction,
                                  targets=['restrRes_spec_1'],
                                  inputs=['inputAtomStruct_1', 'chain_name_spec_1'],
                                  outputs=['restrRes_spec_1'])
SelectResidueWizardQT().addTarget(protocol=ProtHaddockPrediction,
                                  targets=['restrRes_spec_2'],
                                  inputs=['inputAtomStruct_2', 'chain_name_spec_2'],
                                  outputs=['restrRes_spec_2'])
SelectAtomWizardQT().addTarget(protocol=ProtHaddockPrediction,
                               targets=['restrAtom_spec_1'],
                               inputs=['inputAtomStruct_1', 'chain_name_spec_1', 'restrRes_spec_1'],
                               outputs=['restrAtom_spec_1'])
SelectAtomWizardQT().addTarget(protocol=ProtHaddockPrediction,
                               targets=['restrAtom_spec_2'],
                               inputs=['inputAtomStruct_2', 'chain_name_spec_2', 'restrRes_spec_2'],
                               outputs=['restrAtom_spec_2'])


class AddRestraintWizard(AddElementWizard):
  """Add the content of a parameter to another"""
  _targets, _inputs, _outputs = [], {}, {}

  def show(self, form, *params):
    inputParam, outputParam = self.getInputOutput(form)
    protocol = form.protocol
    outIdx = getattr(protocol, inputParam[0]).get()  # 0 == active, 1 == passive

    prevList = self.curePrevList(getattr(protocol, outputParam[outIdx]).get())
    num = outputParam[outIdx].split('_')[-1]
    newList = protocol.createElementLine(num)

    prevList = prevList.strip().split(',')
    prevList = [i for i in prevList if i != '']
    form.setVar(outputParam[outIdx], ','.join(prevList + newList))

AddRestraintWizard().addTarget(protocol=ProtHaddockPrediction,
                               targets=['addRestraint_1'],
                               inputs=['restrType_1'],
                               outputs=['activeRestr_1', 'passiveRestr_1'])
AddRestraintWizard().addTarget(protocol=ProtHaddockPrediction,
                               targets=['addRestraint_2'],
                               inputs=['restrType_2'],
                               outputs=['activeRestr_2', 'passiveRestr_2'])


