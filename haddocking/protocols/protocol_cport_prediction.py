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

import json, os, shutil
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.ResidueDepth import get_surface

from pyworkflow.protocol import params

from pwchem.protocols import ProtDefineStructROIs
from pwchem import Plugin as pwchemPlugin
from pwchem.constants import MGL_DIC
from pwchem.utils import clusterSurfaceCoords, createPocketFile
from pwchem.objects import StructROI, SetOfStructROIs

from .. import Plugin as haddockPlugin


class ProtCPORTPrediction(ProtDefineStructROIs):
  """CPORT protein–protein interface prediction

  AI Generated:

      ProtCPORTPrediction - User Manual

      Overview
      --------
      The ProtCPORTPrediction protocol integrates the CPORT (Consensus Prediction
      Of interface Residues in Transient complexes) method into Scipion-Chem for
      predicting protein–protein interaction interfaces. It identifies residues
      likely involved in binding by combining multiple prediction approaches and
      scoring their consensus.

      This protocol is designed to detect interface regions on protein surfaces,
      which can then be used to define docking restraints, interaction hotspots,
      or structural regions of interest.

      Input Requirements
      ------------------
      1. **Protein Structure**:
         - AtomStruct object containing the receptor structure.
         - Must include valid atomic coordinates.

      2. **Chain Selection**:
         - A specific chain must be selected for interface prediction.
         - Only residues from this chain will be evaluated.

      Workflow
      --------
      1. **Interface Prediction**:
         - Executes CPORT to predict interface residues.
         - Outputs residue-level scores indicating likelihood of participation
           in protein–protein interactions.

      2. **Threshold Filtering**:
         - Residues with scores above the defined threshold are selected.
         - These residues are considered part of the predicted interface.

      3. **Coordinate Extraction**:
         - Atomic coordinates of selected residues are retrieved.
         - Optionally mapped onto the molecular surface.

      4. **Clustering**:
         - Coordinates are clustered into spatially coherent regions.
         - Each cluster represents a candidate interaction pocket.

      5. **Pocket Generation**:
         - Clusters are converted into StructROI objects.
         - Contacts and structural properties are computed.

      Outputs
      -------
      - **SetOfStructROIs**:
        - Each ROI corresponds to a predicted interface region.
        - Includes coordinates, contact information, and structural metadata.

      - **Auxiliary Files**:
        - CPORT prediction output files.
        - Intermediate clustering and pocket definition data.

      Advanced Options
      ----------------
      - Score threshold to control prediction strictness.
      - Clustering parameters (distance thresholds) to define pocket granularity.
      - Optional mapping to surface coordinates for more realistic interface regions.

      Validation & Warnings
      ---------------------
      - Input structure must contain the specified chain.
      - Prediction quality depends on structural completeness and correctness.
      - Very high thresholds may result in no predicted interface regions.

      Practical Recommendations
      -------------------------
      - Use moderate thresholds (e.g., ~0.6) for balanced sensitivity/specificity.
      - Combine with docking protocols (e.g., HADDOCK) for guided protein–protein docking.
      - Inspect resulting clusters to ensure biological relevance.

      Final Perspective
      -----------------
      ProtCPORTPrediction provides an efficient way to identify protein–protein
      interaction interfaces using consensus-based machine learning predictions.
      It facilitates downstream docking, interface analysis, and structural
      modeling workflows within Scipion-Chem.
  """
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
    inFile = self.getInputPDBFile()
    chainStr = json.loads(self.chain_name.get())['chain']

    args = f'{inFile} {chainStr}'
    haddockPlugin.runCPORT(self, args)

  def getSurfaceStep(self):
      parser = PDBParser()
      pdbFile = self.getProteinFileName()

      structure = parser.get_structure(self.getProteinName(), pdbFile)
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
    pdbFile = self.getProteinFileName()
    newPDBFile = self._getPath(os.path.basename(pdbFile))
    shutil.copy(pdbFile, newPDBFile)

    outPockets = SetOfStructROIs(filename=self._getPath('StructROIs.sqlite'))
    for i, clust in enumerate(self.coordsClusters):
      pocketFile = createPocketFile(clust, i, outDir=self._getExtraPath())
      pocket = StructROI(pocketFile, newPDBFile)
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
    structModel = PDBParser().get_structure('receptor', self.getInputPDBFile())[0]

    oCoords = []
    for i, res in enumerate(structModel[chainId].get_residues()):
      if str(i+1) in outRes:
        for a in res.get_atoms():
          oCoords.append(list(a.get_coord()))
    return oCoords


  def getInputPDBFile(self):
    return os.path.abspath(self.inputAtomStruct.get().getFileName())

  def getOutputDir(self):
    return self._getPath('cportOutput')

  def getOutputMLFile(self):
    outDir = self.getOutputDir()
    for file in os.listdir(outDir):
      if file.startswith('cport_ML'):
        return os.path.join(outDir, file)

  def parseMLOutput(self, threshold=0.6):
    '''Parses the output of CPORT and returns a list with the residues over the threshold
    '''
    outRes, outFile = [], self.getOutputMLFile()
    with open(outFile) as f:
      f.readline()
      for line in f:
        sline = line.split(',')
        resID, score = sline[1:3]
        if float(score) >= threshold:
          outRes.append(resID)
    return outRes

  def _validate(self):
    errors = []
    if not self.getProteinFileName(inProtocol=False):
      errors.append('You must specify an input structure')

    return errors

