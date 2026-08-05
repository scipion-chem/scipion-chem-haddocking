# **************************************************************************
# *
# * Authors:     Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
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

import csv, json, os
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.ResidueDepth import get_surface

from pyworkflow.protocol import params

from pwchem.protocols import ProtDefineStructROIs
from pwchem import Plugin as pwchemPlugin
from pwchem.constants import MGL_DIC
from pwchem.utils import clusterSurfaceCoords, createPocketFile, getBaseName, pdbFromASFile, cifFromASFile
from pwchem.objects import StructROI, SetOfStructROIs

from .. import Plugin as haddockPlugin


class ProtCPORTPrediction(ProtDefineStructROIs):
  """Run protein-protein interface predictions using CPORT software from Bonvin Lab"""
  _label = 'CPORT protein interface prediction'

  def _defineParams(self, form):
    form.addSection(label='Input')
    iGroup = form.addGroup('Input')
    iGroup.addParam('inputAtomStruct', params.PointerParam, pointerClass="AtomStruct", label='Input atom struct: ',
                    help="Input AtomStruct to predict the interface on")
    iGroup.addParam('chain_name', params.StringParam, allowsNull=False, label='Chain of interest',
                    help='Specify the chain of the structure of interest')

    iGroup.addParam('threshold', params.FloatParam, label='Score threshold: ', default=0.6,
                    help='Threshold to consider a residue positive in interface')

    group = form.addGroup('Pocket definition')
    group = self._defineClusterParams(group)


  def _insertAllSteps(self):
    self._insertFunctionStep(self.predictStep)
    if self.surfaceCoords:
      self._insertFunctionStep(self.getSurfaceStep)

    self._insertFunctionStep(self.definePocketsStep)
    self._insertFunctionStep(self.defineOutputStep)

  def predictStep(self):
    inpFile = self.getInputPath()
    pdbFile = pdbFromASFile(inpFile, self._getPdbFile(), atomStruct=self.inputAtomStruct.get())
    chainStr = json.loads(self.chain_name.get())['chain']

    args = f'{pdbFile} {chainStr} --pred sppider ispred4 scannet'
    haddockPlugin.runCPORT(self, args)

  def getSurfaceStep(self):
      parser = PDBParser()
      pdbFile = self._getPdbFile()

      structure = parser.get_structure(self._getInputName(), pdbFile)
      self.structModel = structure[0]

      self.structSurface = get_surface(self.structModel,
                                       MSMS=pwchemPlugin.getProgramHome(MGL_DIC, 'MGLToolsPckgs/binaries/msms'))

  def definePocketsStep(self):
    self.coordsClusters = []
    outRes = self.parseMLOutput(self.threshold.get())
    pocketCoords = self.getOutputCoords(outRes)
    if self.surfaceCoords:
      pocketCoords = self.mapSurfaceCoords(pocketCoords)

    self.coordsClusters += clusterSurfaceCoords(pocketCoords, self.maxIntraDistance.get())

  def defineOutputStep(self):
    inpStruct = self.inputAtomStruct.get()
    cifFile = cifFromASFile(self.getInputPath(), self._getCifFile(), inpStruct)

    outPockets = SetOfStructROIs(filename=self._getPath('StructROIs.sqlite'))
    for i, clust in enumerate(self.coordsClusters):
      pocketFile = self._getExtraPath(f'pocketFile_{i}.cif')
      createPocketFile(clust, i, pocketFile)
      pocket = StructROI(pocketFile, cifFile)
      if str(type(inpStruct).__name__) == 'SchrodingerAtomStruct':
        pocket._maeFile = params.String(os.path.relpath(inpStruct.getFileName()))
      pocket.calculateContacts()
      outPockets.append(pocket)

    if len(outPockets) > 0:
      outPockets.buildPDBhetatmFile()
      self._defineOutputs(outputStructROIs=outPockets)


  ################# UTIL FUNCTIONS ##############################

  def getOutputCoords(self, outRes):
    chainId = json.loads(self.chain_name.get())['chain']
    structModel = PDBParser().get_structure('receptor', self._getPdbFile())[0]

    oCoords = []
    for i, res in enumerate(structModel[chainId].get_residues()):
      if str(i+1) in outRes:
        for a in res.get_atoms():
          oCoords.append(list(a.get_coord()))
    return oCoords


  def getInputPath(self):
    return self.inputAtomStruct.get().getFileName()

  def getInputFileName(self):
    return self.getInputPath().split('/')[-1]

  def _getPdbFile(self):
    return os.path.abspath(self._getExtraPath(self._getInputName() + '.pdb'))

  def _getCifFile(self):
    return os.path.abspath(self._getExtraPath(self._getInputName() + '.cif'))

  def _getInputName(self):
    return getBaseName(self.getInputPath())

  def getOutputDir(self):
    return self._getPath('cportOutput')

  def getOutputRawFile(self):
    outDir = self.getOutputDir()
    for file in os.listdir(outDir):
      if file.startswith('predictors_'):
        return os.path.join(outDir, file)

  def parseMLOutput(self, threshold=0.6):
    label_scores = {'A': 1.0, 'AP': 1.0, 'P': 0.5, '-': None}
    with open(self.getOutputRawFile()) as f:
      reader = csv.reader(f)
      resIDs = next(reader)[1:]
      resScores = [[] for _ in resIDs]
      for row in reader:
        for resScore, value in zip(resScores, row[1:]):
          score = label_scores[value] if value in label_scores else float(value)
          if score is not None:
            resScore.append(score)

    outRes = [resID for resID, scores in zip(resIDs, resScores)
              if scores and sum(scores) / len(scores) >= threshold]
    return outRes

  def _validate(self):
    errors = []
    return errors

