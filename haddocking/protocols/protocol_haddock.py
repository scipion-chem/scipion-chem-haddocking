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

import json, os
from Bio.PDB import PDBParser, MMCIFParser

from pyworkflow.protocol import params
from pwem.protocols import EMProtocol
from pwem.objects import SetOfAtomStructs, AtomStruct

from pwchem.utils import clusterSurfaceCoords, createPocketFile, getBaseName, pdbFromAS, getBaseFileName, \
  natural_sort, gunzipFile

from ..constants import HADDOCK_DIC, DEF_PROTOCOL, HADDOCK_HEADER
from .. import Plugin as haddockPlugin

FROM_ACTIVE = 'From active'
RIGDOC, FLEXREF, EMMIN, CLUST, SELMOD, SELCLUST, CAPRI, CONTACT = list(range(8))


class ProtHaddockPrediction(EMProtocol):
  """Run protein-protein interaction predictions using HADDOCK3 software from Bonvin Lab"""
  _label = 'HADDOCK protein interaction prediction'
  _haddockSteps = ['Rigid docking', 'Flexible refinement', 'Energy minimization', 'Clustering FCC',
   'Select top models', 'Select top clusters', 'CAPRI evaluation', 'Contact analysis']

  def _defineParams(self, form):
    form.addSection(label='Input')
    for i in range(1, 3):
      form = self._defineInputParams(form, i)

    # todo: optional: create unambiguous restraints
    # form.addSection(label='Specific restraints')
    # self._defineSpecificRestraints(form)

    form.addSection(label='Ambiguous restraints')
    form.addParam(f'filterSurface', params.FloatParam, default=0.15,
                  label=f'Minimum surface accessibility: ',
                  help=f'Filter the defined restraint residues for a minimum surface accessibility')

    for i in range(1, 3):
      form = self._defineAmbRestraints(form, i)

    form.addSection(label='Execution workflow')
    form.addParam(f'workflowType', params.EnumParam, label='Define workflow: ',
                    choices=['Default', 'Manually'], default=0, display=params.EnumParam.DISPLAY_HLIST,
                    help="Whether to use the default workflow (rigidDock+selection+flexDock+energyMin+analysis) or "
                         "one manually defined")

    form.addParam(f'stepType', params.EnumParam, label='Define step: ',
                  choices=self._haddockSteps, default=0, condition='workflowType==1',
                  help="Which workflow step to add. In any case, the topologies are calculated as first step.")

    iGroup = form.addGroup(f'Rigid docking', condition=self.getManualCondition([RIGDOC]))
    iGroup.addParam('sampling', params.IntParam, default=20, label='Number of generated models: ',
                    condition=self.getManualCondition([RIGDOC]),
                    help='Number of generated models with rigid body fitting.' 
                         'Rigid body energy minimization with CNS')

    sGroup = form.addGroup(f'Model selection', condition=self.getManualCondition([SELMOD]))
    sGroup.addParam('selectModel', params.IntParam, default=10, label='Number of selected models: ',
                    help='Number of selected models with best scores from the previous step')

    sGroup = form.addGroup(f'Flexible refinement/Energy minimization',
                           condition=self.getManualCondition([FLEXREF, EMMIN]))
    sGroup.addParam('tolerance', params.IntParam, default=5, label='Failure tolerance percentage: ',
                    condition=self.getManualCondition([FLEXREF, EMMIN]),
                    help='Percentage of allowed failures for a module to successfully complete')

    scGroup = form.addGroup(f'Cluster selection', condition=self.getManualCondition([SELCLUST]))
    scGroup.addParam('selectCluster', params.IntParam, default=5, label='Number of selected clusters: ',
                     help='Number of selected clusters with best scores from the previous step')

    form.addParam(f'addStep', params.LabelParam, label='Add defined step: ',
                  help="Add defined step to the workflow")

    sumGroup = form.addGroup(f'Workflow', condition='workflowType==1')
    sumGroup.addParam('workflowSum', params.TextParam, default='', label='Workflow steps: ',
                      condition='workflowType==1', width=70, height=10,
                      help='Workflow steps manually defined. In any case, the topologies are calculated as first step.'
                           '\nIf you want to add any HADDOCK parameter, you can do so following the format in '
                           'comma-separated:\n\nsteName:: paramName1: paramValue1, paramName2: paramValue2')

    form.addParallelSection(threads=4, mpi=1)


  def getManualCondition(self, stepNums):
    return f'workflowType==0 or (workflowType==1 and stepType in {stepNums})'

  def _defineInputParams(self, form, num=1):
    iGroup = form.addGroup(f'Input_{num}')
    iGroup.addParam(f'inputAtomStruct_{num}', params.PointerParam, pointerClass="AtomStruct",
                    label=f'Input atom struct {num}: ',
                    help="Input AtomStruct to predict the interaction on")
    return form

  def _defineSpecificRestraints(self, form):
    rGroup = form.addGroup(f'Specific distance restraints')
    for num in range(1, 3):
      rGroup = self._defineSelectionParams(rGroup, num, atomCond='True', specName='_spec')

  def _defineSelectionParams(self, group, num, cond='True', atomCond='False', specName=''):
    group.addParam(f'chain_name{specName}_{num}', params.StringParam, label=f'Chain of interest {num}', condition=cond,
                   help='Specify the chain of the structure of interest')
    group.addParam(f'restrRes{specName}_{num}', params.StringParam, label=f'Residue(s) of interest {num}: ',
                   condition=cond,
                   help=f"Select a list of residues of ambiguous restraints in {num} or a single residue to "
                        f"define specific atom restraints")

    group.addParam(f'restrAtom{specName}_{num}', params.StringParam, label=f'Atom of interest {num}: ',
                   condition=f'{cond} and {atomCond}',
                   help=f"Select an atom to define specific restraints for input {num}")
    return group


  def _defineAmbRestraints(self, form, num=1):
    rGroup = form.addGroup(f'Ambiguous restraints_{num}')
    rGroup.addParam(f'restrType_{num}', params.EnumParam, label=f'Restraints type {num}: ',
                    choices=['Active', 'Passive'], default=0, display=params.EnumParam.DISPLAY_HLIST,
                    help=f"Type of restraints to define for input {num}")
    rGroup.addParam(f'restrSource{num}', params.EnumParam, label=f'Define restraints {num} from: ',
                    choices=['Structural ROI', 'Manual'], default=0, display=params.EnumParam.DISPLAY_HLIST,
                    help=f"Source of restraints for input {num}")
    rGroup.addParam(f'restrROIs_{num}', params.PointerParam, label=f'ROIs set for active restraints {num}: ',
                    condition=f'restrSource{num}==0', pointerClass='SetOfStructROIs', allowsNull=True,
                    help=f"Input ROIs active restraints for input {num}")
    rGroup.addParam(f'restrROI_{num}', params.StringParam, label=f'Chosen ROI for active restraints {num}: ',
                    condition=f'restrSource{num}==0',
                    help=f"Input ROIs active restraints for input {num}")

    rGroup = self._defineSelectionParams(rGroup, num, cond=f'restrSource{num}==1')

    rGroup.addParam(f'addRestraint_{num}', params.LabelParam, label=f'Add defined restraints {num}',
                    help=f'Add defined active restraints for input {num}')
    rGroup.addParam(f'activeRestr_{num}', params.TextParam, width=70, height=2, default='',
                    label=f'Active restraints summary {num}: ',
                    help=f'Summary of the active restraints defined for input {num}')
    rGroup.addParam(f'passiveRestr_{num}', params.TextParam, width=70, height=3, default=f'{FROM_ACTIVE}\n',
                    label=f'Passive restraints summary {num}: ',
                    help=f'Summary of the passive restraints defined for input {num}.'
                         f'{FROM_ACTIVE} will automatically define those residues close to the active ones defined.'
                         f'You can be delete and define residues manually or add further passive residues')
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
      inChain, inRes = getattr(self, f'chain_name_{num}').get(), getattr(self, f'restrRes_{num}').get()
      chName = json.loads(inChain)['chain']
      resIdxs = list(map(int, json.loads(inRes)['index'].split('-')))
      rList = [f'{chName}_{i}' for i in range(resIdxs[0], resIdxs[1] + 1)]
    return rList

  def createStepLine(self):
      stepNum = self.stepType.get()
      stepName = self._haddockSteps[stepNum]

      s = f'{stepName}:: '
      if stepNum == RIGDOC:
          s += f'sampling: {self.sampling.get()}'
      elif stepNum in [FLEXREF, EMMIN]:
          s += f'tolerance: {self.tolerance.get()}'
      elif stepNum in [SELMOD]:
        s += f'select: {self.selectModel.get()}'
      elif stepNum in [SELCLUST]:
        s += f'top_models: {self.selectCluster.get()}'
      return s


  def _insertAllSteps(self):
    self._insertFunctionStep(self.convertStep)
    self._insertFunctionStep(self.defineRestrictionsStep)
    self._insertFunctionStep(self.predictStep)
    self._insertFunctionStep(self.defineOutputStep)
  
  ############### STEP FUNCTIONS #################### 
  
  def convertStep(self):
    for num in range(1, 3):
      pdbFromAS(self.getInputStruct(num), self.getPDBFile(num))
  
  def defineRestrictionsStep(self):
    for num in range(1, 3):
      self.createMergedPDB(num)

    # Create restraint files to keep the chains of each input together
    self.restraintChains()
    self.createRestraintsFile()
  
  def predictStep(self):
    confFile = self.buildConfFile()
    haddockPlugin.runHaddock(self, 'haddock3', confFile, cwd=self._getExtraPath())

  def defineOutputStep(self):

    modelFiles = self.getOutputModelFiles()
    outSet = SetOfAtomStructs.create(self._getPath())
    for mFile in modelFiles:
        scoresDic = self.parseOutputFile(mFile)

        outAS = AtomStruct(filename=mFile)
        for scoreName, score in scoresDic.items():
            setattr(outAS, scoreName, params.Float(score))

        outSet.append(outAS.clone())

    self._defineOutputs(**{'outputStructures': outSet})


  ############### MAIN FUNCTIONS #################### 
  
  def createMergedPDB(self, num):
    # Creating merged chains file needed for HADDOCK
    segId = self.getSegId(num)
    inFile, mergedFile = self.getPDBFile(num), self.getMergedFile(num)

    mergeArgs = f'pdb_chain -{segId} {inFile} | pdb_reres -1 > {mergedFile}'
    haddockPlugin.runScript(self, '', mergeArgs, HADDOCK_DIC, scriptDir='', pyStr='')
    return mergedFile
  
  def restraintChains(self):
    oStr = ''
    for num in range(1,3):
      mFile = self.getMergedFile(num)
      oStr += haddockPlugin.runHaddock(self, 'haddock3-restraints restrain_bodies ', mFile,
                                      cwd=self._getExtraPath(), popen=True).strip() + '\n'

    if oStr:
      with open(self.getChainRestrainstsFile(), 'w') as f:
        f.write(oStr)
  
  def createRestraintsFile(self):
    mapRes, segIds = {}, ['A', 'B']
    for i in range(2):
      mapRes[segIds[i]] = self.getMappedResidues(i+1)

    rFile = self.getAmbRestrFile()
    with open(rFile, 'w') as f:
      for i in range(2):
        segId1, segId2 = segIds[i], segIds[i-1]
        f.write(f'! HADDOCK AIR restraints for partner {segId1}\n!\n')
        for actNum in mapRes[segId1]['active']:
          f.write(f'assign ( resid {actNum} and segid {segId1})\n\t(\n')
          actStrs = [f'\t ( resid {actNum} and segid {segId2})' for actNum in mapRes[segId2]['active']]
          pasStrs = [f'\t ( resid {pasNum} and segid {segId2})' for pasNum in mapRes[segId2]['passive']]
          f.write(f'\n\tor\n'.join(actStrs + pasStrs))
          f.write(f'\n\t)  2.0 2.0 0.0\n!\n')
  
  def getMappedResidues(self, num):
    '''Returns a dictionary {'active': [], 'passive': []} with the mapped numbers of the defined residues
    '''
    restrDic = self.parseDefinedRes(num)
    restrDic = self.filterSASA(restrDic, num)
    return restrDic

  def buildConfFile(self):
    confFile = self.getConfFile()
    cwd = self.getRunDir()
    nCores = self.numberOfThreads.get() - 1 if self.numberOfThreads.get() > 1 else 1
    inFile1, inFile2 = self.getMergedFile(1), self.getMergedFile(2)
    ambRestrFile, unambRestrFile = self.getAmbRestrFile(), self.getChainRestrainstsFile()

    confStr = HADDOCK_HEADER.format(cwd, nCores, inFile1, inFile2)
    if self.workflowType.get() == 0:
        tol = self.tolerance.get()
        sampling, select = self.sampling.get(), self.selectModel.get()
        nPerCluster = self.selectCluster.get()
        confStr += DEF_PROTOCOL.format(ambRestrFile, unambRestrFile, sampling, select, tol, ambRestrFile,
                                       unambRestrFile, tol, ambRestrFile, unambRestrFile, nPerCluster)
    else:
        stepsList = self.parseManualWorkflow()
        confStr += self.buildConfStr(stepsList, ambRestrFile, unambRestrFile)

    with open(confFile, 'w') as f:
        f.write(confStr)
    return confFile

  def parseManualWorkflow(self):
      stepsList = []
      for line in self.workflowSum.get().split('\n'):
        line = line.split(')')[-1]
        if line.strip():
          stepName, argsStr = line.split('::')

          curStep = {}
          for argsPair in argsStr.split(','):
            if ':' in argsPair:
                argName, argValue = argsPair.split(':')
                curStep[argName] = argValue

          stepsList.append((stepName, curStep))

      return stepsList

  def buildConfStr(self, stepsList, ambRestrFile, unambRestrFile):
      s = '[topoaa]\n\n'
      h = self._haddockSteps
      for stepName, stepDic in stepsList:
          stepName = stepName.strip()
          if stepName == h[RIGDOC]:
              s += f'[rigidbody]\nambig_fname = "{ambRestrFile}"\nunambig_fname = "{unambRestrFile}"\n'
          elif stepName == h[FLEXREF]:
              s += f'[flexref]\nambig_fname = "{ambRestrFile}"\nunambig_fname = "{unambRestrFile}"\n'
          elif stepName == h[EMMIN]:
              s += f'[emref]\nambig_fname = "{ambRestrFile}"\nunambig_fname = "{unambRestrFile}"\n'
          elif stepName == h[CLUST]:
              s += f'[clustfcc]\n'
          elif stepName == h[SELMOD]:
              s += f'[seletop]\n'
          elif stepName == h[SELCLUST]:
              s += f'[seletopclusts]\n'
          elif stepName == h[CAPRI]:
              s += f'[caprieval]\n'
          elif stepName == h[CONTACT]:
              s += f'[contactmap]\n'

          for argName, argValue in stepDic.items():
              s += f'{argName} = {argValue}\n'

          s += '\n'
      return s

  ############### PARSER FUNCTIONS ####################
  
  def parseDefinedRes(self, num):
    '''Returns a dictionary {'active': [], 'passive': []} with the defined active and passive residues.
    Residues present in both lists are deleted from passive
    '''
    mapDic = self.buildMapDic(num)
    actStr, pasStr = getattr(self, f'activeRestr_{num}').get().strip(), \
                     getattr(self, f'passiveRestr_{num}').get().strip()

    act = [mapDic[s] for line in actStr.split('\n') for s in line.split(',')]
    pas = [mapDic[s] for line in pasStr.split('\n') for s in line.split(',') if s != FROM_ACTIVE]

    if FROM_ACTIVE in pasStr:
        segId = self.getSegId(num)
        mergedFile = self.getMergedFile(num)
        args = f'{mergedFile} -c {segId} {",".join(act)}'
        nPas = haddockPlugin.runHaddock(self, 'haddock3-restraints passive_from_active ', args,
                                        cwd=self._getExtraPath(), popen=True).split()
        pas = list(set(pas).union(set(nPas))) if pas else nPas

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

  def parseOutputFile(self, mFile):
    sDic = {}
    with open(mFile) as f:
        for line in f:
            if line.startswith('REMARK HADDOCK score:'):
                sDic['HaddockScore'] = float(line.split()[-1])
            elif line.startswith('REMARK Binding energy:'):
                sDic['BindingEnergy'] = float(line.split()[-1])
            elif line.startswith('REMARK buried surface area:'):
                sDic['BuriedSurfaceArea'] = float(line.split()[-1])
    return sDic


  ############### UTIL FUNCTIONS ####################

  def getSegId(self, num):
    return 'A' if num == 1 else 'B'

  def calculateFreeSASA(self, num):
    inFile = self.getMergedFile(num)
    oFile = self._getExtraPath(f'input_{num}_SASA.txt')
    args = f'{os.path.abspath(inFile)} {os.path.abspath(oFile)}'
    haddockPlugin.runScript(self, 'calculateFreeSASA.py', args, HADDOCK_DIC, cwd=self._getExtraPath())
    return oFile

  def filterSASA(self, defRes, num):
    segId = self.getSegId(num)
    sasaDic = self.parseSASAFile(self.calculateFreeSASA(num))
    minSASA = self.filterSurface.get()

    nDefRes = {}
    for k in defRes:
      nDefRes[k] = []
      for res in defRes[k]:
        if res.strip():
          resId = f'{segId}_{res}'
          if float(sasaDic[resId][0]) >= minSASA or \
                  (sasaDic[resId][1] != 'nan' and float(sasaDic[resId][1]) >= minSASA):
            nDefRes[k].append(res)

    for k in nDefRes:
      nDefRes[k].sort()
    return nDefRes

  def buildMapDic(self, num):
    '''Returns a dictionary with the mapped residues from the original to the merged form
    {"ori_Id": "mergeId"}
    '''
    resNumDic = {}
    inASFile = self.getInputFile(num)
    parser = PDBParser if inASFile.endswith('.pdb') else MMCIFParser
    structModel = parser().get_structure(getBaseName(inASFile), inASFile)[0]

    c = 1
    for chain in structModel.get_chains():
      for res in chain.get_residues():
        resNumDic[f'{chain.get_id()}_{res.get_id()[1]}'] = str(c)
        c += 1
    return resNumDic


  ############### FILES FUNCTIONS ####################
  
  def getInputStruct(self, num):
    return getattr(self, f'inputAtomStruct_{num}').get()

  def getInputFile(self, num):
    return os.path.abspath(self.getInputStruct(num).getFileName())

  def getPDBFile(self, num):
    pdbFile = getBaseFileName(self.getInputFile(num)).replace('.cif', '.pdb')
    return os.path.abspath(self._getExtraPath(pdbFile))

  def getMergedFile(self, num):
    return os.path.abspath(self._getExtraPath(f'inputFile_{num}.pdb'))

  def getChainRestrainstsFile(self):
    return os.path.abspath(self._getExtraPath(f'chainRestraints.tbl'))

  def getAmbRestrFile(self):
    return os.path.abspath(self._getExtraPath(f'restraints.tbl'))

  def getRunDir(self):
      return os.path.abspath(self._getExtraPath('execution'))

  def getConfFile(self):
    return os.path.abspath(self._getExtraPath('confFile.txt'))

  def getOutputModelFiles(self):
    runDir = self.getRunDir()
    oDirs = [os.path.join(runDir, oDir) for oDir in os.listdir(runDir) if len(oDir.split('_')) == 2]
    oDirs = natural_sort(oDirs, rev=True)

    mFiles, i = [], 0
    while i < len(oDirs):
        oDir = oDirs[i]
        for file in os.listdir(oDir):
            if file.endswith('.pdb.gz'):
                oFile = self._getPath(file.replace('.pdb.gz', '.pdb'))
                gunzipFile(os.path.join(oDir, file), oFile)
                mFiles.append(oFile)
                i = len(oDirs)
        i += 1
    return mFiles

  ############### MANAGER FUNCTIONS ####################,

  def _validate(self):
    errors = []
    return errors

