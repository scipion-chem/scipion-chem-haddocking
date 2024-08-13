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

from pyworkflow.protocol import params
from pwem.protocols import EMProtocol

from pwchem.utils import clusterSurfaceCoords, createPocketFile
from pwchem.objects import StructROI, SetOfStructROIs

from ..constants import HADDOCK_DIC
from .. import Plugin as haddockPlugin


class ProtHaddockPrediction(EMProtocol):
  """Run protein-protein interaction predictions using HADDOCK3 software from Bonvin Lab"""
  _label = 'HADDOCK protein interaction prediction'

  def _defineParams(self, form):
    form.addSection(label='Input')
    for i in range(1, 3):
      form = self._defineInputParams(form, i)

    form.addSection(label='Restraints')
    form.addParam(f'filterSurface', params.FloatParam, default=0.15,
                  label=f'Minimum surface accessibility: ',
                  help=f'Filter the defined restraint residues for a minimum surface accessibility')

    for i in range(1, 3):
      form = self._defineRestraints(form, i)


  def _defineInputParams(self, form, num=1):
    iGroup = form.addGroup(f'Input_{num}')
    iGroup.addParam(f'inputAtomStruct_{num}', params.PointerParam, pointerClass="AtomStruct",
                    label=f'Input atom struct {num}: ',
                    help="Input AtomStruct to predict the interaction on")
    return form

  def _defineRestraints(self, form, num=1):
    rGroup = form.addGroup(f'Restraints_{num}')
    rGroup.addParam(f'restrType_{num}', params.EnumParam, label=f'Restraints type {num}: ',
                    choices=['Active', 'Passive'], default=0, display=params.EnumParam.DISPLAY_HLIST,
                    help=f"Type of defining restraints for input {num}")
    rGroup.addParam(f'restrSource{num}', params.EnumParam, label=f'Define restraints {num} from: ',
                    choices=['Structural ROI', 'Manual'], default=0, display=params.EnumParam.DISPLAY_HLIST,
                    help=f"Source of restraints for input {num}")
    rGroup.addParam(f'restrROIs_{num}', params.PointerParam, label=f'ROIs set for active restraints {num}: ',
                    condition=f'restrSource{num}==0', pointerClass='SetOfStructROIs', allowsNull=True,
                    help=f"Input ROIs active restraints for input {num}")
    rGroup.addParam(f'restrROI_{num}', params.StringParam, label=f'Chosen ROI for active restraints {num}: ',
                    condition=f'restrSource{num}==0',
                    help=f"Input ROIs active restraints for input {num}")

    rGroup.addParam(f'chain_name_{num}', params.StringParam, label=f'Chain of interest {num}',
                    condition=f'restrSource{num}==1',
                    help='Specify the chain of the structure of interest')
    rGroup.addParam(f'restrMan_{num}', params.StringParam, label=f'Active manual restraints {num}: ',
                    condition=f'restrSource{num}==1',
                    help=f"Manual active restraints for input {num}")

    rGroup.addParam(f'addRestraint_{num}', params.LabelParam, label=f'Add defined restraints {num}',
                    help=f'Add defined active restraints for input {num}')
    rGroup.addParam(f'activeRestr_{num}', params.TextParam, width=70, height=2, default='',
                    label=f'Active restraints summary {num}: ',
                    help=f'Summary of the active restraints defined for input {num}')
    rGroup.addParam(f'passiveRestr_{num}', params.TextParam, width=70, height=3, default='',
                    label=f'Passive restraints summary {num}: ',
                    help=f'Summary of the passive restraints defined for input {num}')
    return form

  def getInputROI(self, num=1):
    roiStr = getattr(self, f'restrROI_{num}').get()
    for roi in getattr(self, f'restrROIs_{num}').get():
      if roi.__str__() == roiStr:
        return roi

  def createElementLine(self, num=1):
    if getattr(self, f'restrSource{num}').get() == 0:
      inROI = self.getInputROI(num)
      rList = inROI.getDecodedCResidues()

    else:
      inChain, inRes = getattr(self, f'chain_name_{num}').get(), getattr(self, f'restrMan_{num}').get()
      chName = json.loads(inChain)['chain']
      resIdxs = list(map(int, json.loads(inRes)['index'].split('-')))
      rList = [f'{chName}_{i}' for i in range(resIdxs[0], resIdxs[1] + 1)]
    return rList

  def _insertAllSteps(self):
    self._insertFunctionStep(self.defineRestrictionsStep)
    self._insertFunctionStep(self.predictStep)
    self._insertFunctionStep(self.definePocketsStep)
    self._insertFunctionStep(self.defineOutputStep)

  def getInputFile(self, num):
    return os.path.abspath(getattr(self, f'inputAtomStruct_{num}').get().getFileName())

  def getMergedFile(self, num):
    return os.path.abspath(self._getExtraPath(f'inputFile_{num}.pdb'))

  def createMergedPDB(self, num):
    # Creating merged chains file needed for HADDOCK
    segId = 'A' if num == 1 else 'B'
    inFile, mergedFile = self.getInputFile(num), self.getMergedFile(num)

    mergeArgs = f'pdb_chain -{segId} {inFile} | pdb_reres -1 > {mergedFile}'
    haddockPlugin.runScript(self, '', mergeArgs, HADDOCK_DIC, scriptDir='', pyStr='')
    return mergedFile

  def getResNumPerChain(self, num):
    '''Returns a dictionary with the number of residues per chain in the model
    '''
    resNumDic = {}
    structModel = PDBParser().get_structure('receptor', self.getInputFile(num))[0]
    for chain in structModel.get_chains():
      resNumDic[chain.get_id()] = len(list(chain.get_residues()))
    return resNumDic

  def parseDefinedRes(self, num):
    '''Returns a dictionary {'active': [], 'passive': []} with the defined active and passive residues.
    Residues present in both lists are deleted from passive
    '''
    act, pas = getattr(self, f'activeRestr_{num}').get().strip().split(','), \
               getattr(self, f'passiveRestr_{num}').get().strip().split(',')
    pas = list(set(pas).difference(set(act)))
    return {'active': act, 'passive': pas}

  def parseSASAFile(self, sasaFile):
    '''Return a dictionary {chain_resnum: [backboneSASA, sideSASA]}'''
    sasaDic = {}
    with open(sasaFile) as f:
      for line in f:
        sline = line.split()
        sasaDic[sline[0]] = [sline[2], sline[3]]
    return sasaDic

  def filterSASA(self, defRes, num):
    sasaFile = self.calculateFreeSASA(num)
    sasaDic = self.parseSASAFile(sasaFile)
    minSASA = self.filterSurface.get()

    nDefRes = {}
    for k in defRes:
      nDefRes[k] = []
      for res in defRes[k]:
        if float(sasaDic[res][0]) >= minSASA or (sasaDic[res][1] != 'nan' and float(sasaDic[res][1]) >= minSASA):
          nDefRes[k].append(res)
    return nDefRes

  def getMappedResidues(self, num):
    '''Returns a dictionary {'active': [], 'passive': []} with the mapped numbers of the defined residues
    '''
    mapDic = {'active': [], 'passive': []}
    resNumDic = self.getResNumPerChain(num)
    defRes = self.parseDefinedRes(num)
    defRes = self.filterSASA(defRes, num)

    for k in defRes:
      for res in defRes[k]:
        chain, resnum = res.split('_')
        c = 0
        for mChain in resNumDic:
          if mChain == chain:
            c += int(resnum)
            break
          else:
            c += resNumDic[mChain]

        mapDic[k].append(c)
      mapDic[k].sort()
    return mapDic

  def createRestraintsFile(self):
    mapDic, segIds = {}, ['A', 'B']
    for i in range(2):
      mapDic[segIds[i]] = self.getMappedResidues(i+1)

    rFile = self._getExtraPath(f'restraints.tbl')
    with open(rFile, 'w') as f:
      for i in range(2):
        segId1, segId2 = segIds[i], segIds[i-1]
        f.write(f'! HADDOCK AIR restraints for partner {segId1}\n!\n')
        for actNum in mapDic[segId1]['active']:
          f.write(f'assign ( resid {actNum} and segid {segId1})\n\t(\n')
          actStrs = [f'\t ( resid {actNum} and segid {segId2})' for actNum in mapDic[segId2]['active']]
          pasStrs = [f'\t ( resid {pasNum} and segid {segId2})' for pasNum in mapDic[segId2]['passive']]
          f.write(f'\n\tor\n'.join(actStrs + pasStrs))
          f.write(f'\n\t)  2.0 2.0 0.0\n!\n')


  def calculateFreeSASA(self, num):
    inFile = self.getInputFile(num)
    oFile = self._getExtraPath(f'input_{num}_SASA.txt')
    args = f'{os.path.abspath(inFile)} {os.path.abspath(oFile)}'
    haddockPlugin.runScript(self, 'calculateFreeSASA.py', args, HADDOCK_DIC, cwd=self._getExtraPath())
    return oFile

  def defineRestrictionsStep(self):
    for num in range(1, 3):
      self.createMergedPDB(num)

    self.createRestraintsFile()


  def predictStep(self):
    inFile = self.getInputPDBFile()
    chainStr = json.loads(self.chain_name.get())['chain']

    args = f'{inFile} {chainStr}'
    haddockPlugin.runCPORT(self, args)

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
    return errors

