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

import json, os, re
from Bio.PDB import PDBParser, NeighborSearch

from pyworkflow.protocol import params
from pwem.protocols import EMProtocol
from pwem.objects import SetOfAtomStructs, AtomStruct

from pwchem.utils import clusterSurfaceCoords, createPocketFile, getBaseName, pdbFromASFile, getBaseFileName, \
  natural_sort, gunzipFile

from pwchem import Plugin as pwchemPlugin

from ..constants import HADDOCK_DIC, DEF_PROTOCOL, HADDOCK_HEADER, LIGAND_RESNAMES, WATER_RESNAMES, TOPOLOGY_TERMS
from .. import Plugin as haddockPlugin

FROM_ACTIVE = 'From active'
AUTO_TOP, USER_TOP = 0, 1
COVALENT_CUTOFF = 1.9
RIGDOC, FLEXREF, EMMIN, CLUST, SELMOD, SELCLUST, CAPRI, CONTACT = list(range(8))


class ProtHaddockPrediction(EMProtocol):
  """Run interaction predictions using HADDOCK3 software from Bonvin Lab between any pair of molecules
  (proteins, nucleic acids, glycans, lipids or small ligands)"""
  _label = 'HADDOCK interaction prediction'
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
                    help="Input AtomStruct to predict the interaction on. It can be a protein, nucleic acid, "
                         "glycan, lipid or small molecule, or any combination of them.")
    iGroup.addParam(f'ligTopSource_{num}', params.EnumParam, label=f'Topology of non-standard residues {num}: ',
                    choices=['Automatic (ACPYPE)', 'User files'], default=AUTO_TOP,
                    display=params.EnumParam.DISPLAY_HLIST,
                    help="HADDOCK only knows the residues of its CNS library (amino acids, nucleotides, common glycans, "
                         "cofactors and ions). Any other residue (e.g. lipid A, KDO, fatty acids, drugs) needs a CNS "
                         "topology and parameter file. The HADDOCK web server generates them internally, the local "
                         "HADDOCK3 does not. This parameter is ignored if all the residues of the input are supported."
                         "\n\nAutomatic: the non-library hetero residues (plus the library hetero residues covalently "
                         "bonded to them, like the phosphates of a lipid A) are merged into a single ligand residue, "
                         "protonated with OpenBabel and parametrized with ACPYPE (GAFF2). The topology is then converted "
                         "to united atom (non-polar hydrogens merged into their carbons, as in the HADDOCK library), "
                         "which keeps big ligands like polysaccharides within the CNS topology size limits. All "
                         "restraints defined on any of those residues are applied to the merged ligand."
                         "\n\nUser files: CNS topology and parameter files for the non-library residues, e.g. "
                         "ligand.top and ligand.param from the toppar folder of a HADDOCK web server run archive. "
                         "Residue and atom names must match those of the input structure.")
    iGroup.addParam(f'ligPH_{num}', params.FloatParam, default=7.4, label=f'Ligand protonation pH {num}: ',
                    condition=f'ligTopSource_{num}=={AUTO_TOP}', expertLevel=params.LEVEL_ADVANCED,
                    help="pH used by OpenBabel to add the hydrogens of the ligand")
    iGroup.addParam(f'ligCharge_{num}', params.StringParam, default='', label=f'Ligand net charge {num}: ',
                    condition=f'ligTopSource_{num}=={AUTO_TOP}', expertLevel=params.LEVEL_ADVANCED,
                    help="Net charge of the ligand. If empty, it is estimated from the OpenBabel protonation")
    iGroup.addParam(f'ligChargeMethod_{num}', params.EnumParam, default=0, label=f'Ligand charge method {num}: ',
                    choices=['gas', 'bcc'], condition=f'ligTopSource_{num}=={AUTO_TOP}',
                    display=params.EnumParam.DISPLAY_HLIST, expertLevel=params.LEVEL_ADVANCED,
                    help="ACPYPE partial charges method: gas (Gasteiger, fast) or bcc (AM1-BCC, more accurate but "
                         "very slow for big ligands like LPS)")
    iGroup.addParam(f'removeIons_{num}', params.BooleanParam, default=True, label=f'Remove isolated ions {num}: ',
                    condition=f'ligTopSource_{num}=={AUTO_TOP}', expertLevel=params.LEVEL_ADVANCED,
                    help="Remove the single atom ions (e.g. Mg2+ bridging LPS phosphates) of the input when a ligand "
                         "is parametrized. HADDOCK flexible refinement (torsion angle dynamics) fails with isolated, "
                         "not covalently bonded atoms in the same molecule as the ligand")
    iGroup.addParam(f'ligTopFile_{num}', params.FileParam, label=f'CNS topology file {num}: ',
                    condition=f'ligTopSource_{num}=={USER_TOP}', allowsNull=True,
                    help="CNS topology file (.top) of the non-standard residues of the input")
    iGroup.addParam(f'ligParFile_{num}', params.FileParam, label=f'CNS parameter file {num}: ',
                    condition=f'ligTopSource_{num}=={USER_TOP}', allowsNull=True,
                    help="CNS parameter file (.param/.par) of the non-standard residues of the input")
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
    self._insertFunctionStep(self.ligandTopologyStep)
    self._insertFunctionStep(self.defineRestrictionsStep)
    self._insertFunctionStep(self.predictStep)
    self._insertFunctionStep(self.defineOutputStep)
  
  ############### STEP FUNCTIONS #################### 
  
  def convertStep(self):
    for num in range(1, 3):
      pdbFromASFile(self.getInputFile(num), self.getPDBFile(num), self.getInputStruct(num))

  def ligandTopologyStep(self):
    supported = haddockPlugin.getSupportedResidues(self, cwd=self._getExtraPath())
    ligInfo = {}
    for num in range(1, 3):
      unsupported = self.getUnsupportedResnames(num, supported)
      ligKeys, removedKeys = [], []
      if unsupported and self.getTopSource(num) == AUTO_TOP:
        ligKeys = self.getLigandResidueKeys(num, unsupported)
        if getattr(self, f'removeIons_{num}').get():
          removedKeys = self.getIonResidueKeys(num)
        self.parametrizeLigand(num, ligKeys, removedKeys)
      ligInfo[str(num)] = {'unsupported': unsupported, 'ligandKeys': ligKeys, 'removedKeys': removedKeys}

    with open(self.getLigandInfoFile(), 'w') as f:
      json.dump(ligInfo, f, indent=2)
    self.joinLigandTopologies()

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
    segId = self.getSegId(num)
    inFile, mergedFile, ligFile = self.getPDBFile(num), self.getMergedFile(num), self.getLigandPDBFile(num)
    if not os.path.exists(ligFile):
      self.runChainReres(inFile, mergedFile, segId)
      return mergedFile

    keptFile, atomLines, nKept = self.getNoLigandFile(num), [], 0
    if self.getAtomLines(keptFile):
      self.runChainReres(keptFile, mergedFile, segId)
      atomLines = self.getAtomLines(mergedFile)
      nKept = len({line[22:27] for line in atomLines})

    resName, ligResid = LIGAND_RESNAMES[num - 1], nKept + 1
    for line in self.getAtomLines(ligFile):
      atomLines.append(f'HETATM{line[6:17]}{resName:>3} {segId}{ligResid:>4} {line[27:]}')

    with open(mergedFile, 'w') as f:
      for i, line in enumerate(atomLines):
        f.write(f'{line[:6]}{i + 1:>5}{line[11:]}')
      f.write('END\n')
    return mergedFile

  def runChainReres(self, inFile, outFile, segId):
    mergeArgs = f'pdb_chain -{segId} {inFile} | pdb_reres -1 > {outFile}'
    haddockPlugin.runScript(self, '', mergeArgs, HADDOCK_DIC, scriptDir='', pyStr='')

  def getUnsupportedResnames(self, num, supported):
    resNames = {line[17:20].strip() for line in self.getAtomLines(self.getPDBFile(num))}
    return sorted(resNames.difference(supported).difference(WATER_RESNAMES))

  def getLigandResidueKeys(self, num, unsupported):
    structModel = PDBParser(QUIET=True).get_structure('input', self.getPDBFile(num))[0]
    residues = list(structModel.get_residues())
    candidates = [res for res in residues if res.get_id()[0] != ' ' and
                  res.get_resname().strip() not in WATER_RESNAMES and len(res) > 1]
    nSearch = NeighborSearch([atom for res in candidates for atom in res if atom.element != 'H'])

    ligRes = {res for res in residues if res.get_resname().strip() in unsupported}
    queue = list(ligRes)
    while queue:
      res = queue.pop()
      for atom in res:
        for nRes in nSearch.search(atom.coord, COVALENT_CUTOFF, level='R'):
          if nRes not in ligRes:
            ligRes.add(nRes)
            queue.append(nRes)

    return [self.getResidueKey(res) for res in residues if res in ligRes]

  def getIonResidueKeys(self, num):
    structModel = PDBParser(QUIET=True).get_structure('input', self.getPDBFile(num))[0]
    return [self.getResidueKey(res) for res in structModel.get_residues()
            if res.get_id()[0] != ' ' and len(res) == 1 and res.get_resname().strip() not in WATER_RESNAMES]

  def parametrizeLigand(self, num, ligKeys, removedKeys):
    ligDir, resName = self.getLigandDir(num), LIGAND_RESNAMES[num - 1]
    os.makedirs(ligDir, exist_ok=True)
    heavyFile, mol2File = os.path.join(ligDir, 'ligand_heavy.pdb'), os.path.join(ligDir, 'ligand_H.mol2')

    ligKeys, removedKeys, ligLines, keptLines = set(ligKeys), set(removedKeys), [], []
    for line in self.getAtomLines(self.getPDBFile(num)):
      key = self.getLineKey(line)
      if key in removedKeys:
        continue
      elif key not in ligKeys:
        keptLines.append(line)
      elif self.getLineElement(line) != 'H':
        ligLines.append(line)

    if removedKeys:
      self.info(f'Input {num}: removed isolated ions {", ".join(sorted(removedKeys))}')

    for oFile, lines in [(heavyFile, ligLines), (self.getNoLigandFile(num), keptLines)]:
      with open(oFile, 'w') as f:
        f.write(''.join(lines) + 'END\n')

    pH = getattr(self, f'ligPH_{num}').get()
    pwchemPlugin.runOPENBABEL(self, args=f'-ipdb {heavyFile} -omol2 -O {mol2File} -p {pH}', cwd=ligDir)

    charge = self.getLigandCharge(num, mol2File)
    chargeMethod = ['gas', 'bcc'][getattr(self, f'ligChargeMethod_{num}').get()]
    acArgs = f'-i {mol2File} -b {resName} -c {chargeMethod} -n {charge} -a gaff2 -o cns'
    pwchemPlugin.runACPYPE(self, args=acArgs, cwd=ligDir)

    acBase = os.path.join(ligDir, f'{resName}.acpype', resName)
    with open(f'{acBase}_CNS.top') as f:
      topStr = re.sub(r'(RESIdue\s+)[^\s}]+', rf'\g<1>{resName}', f.read())
    topStr, droppedH = self.toUnitedAtomTopology(topStr)
    self.info(f'Input {num}: ligand {resName} converted to united atom ({len(droppedH)} non-polar hydrogens merged)')
    with open(self.getLigandTopFile(num), 'w') as f:
      f.write(topStr)
    with open(f'{acBase}_CNS.par') as fIn, open(self.getLigandParFile(num), 'w') as fOut:
      fOut.write(fIn.read())
    ligLines = [line for line in self.getAtomLines(f'{acBase}_NEW.pdb') if line[12:16].strip() not in droppedH]
    with open(self.getLigandPDBFile(num), 'w') as f:
      f.write(''.join(ligLines) + 'END\n')

  def toUnitedAtomTopology(self, topStr):
    lines = topStr.splitlines(keepends=True)
    types, charges, bonds = self.parseTopologyAtoms(lines)
    droppedH = self.getNonPolarHydrogens(types, bonds)
    for hAtom, cAtom in droppedH.items():
      charges[cAtom] += charges[hAtom]

    outLines = [self.setTopologyCharge(line, charges) for line in lines
                if not self.hasTopologyAtom(line, droppedH)]
    return ''.join(outLines), set(droppedH)

  def getTopologyKeyword(self, line):
    tokens = line.split()
    return (tokens[0][:4].upper(), tokens[1:]) if tokens else ('', [])

  def parseTopologyAtomLine(self, line):
    tokens = line.replace('=', ' = ').split()
    upperTokens = [token.upper() for token in tokens]
    atomType = tokens[upperTokens.index('TYPE') + 2].lower()
    charge = float(tokens[upperTokens.index('CHARGE') + 2])
    return tokens[1], atomType, charge

  def parseTopologyAtoms(self, lines):
    types, charges, bonds = {}, {}, []
    for line in lines:
      keyword, tokens = self.getTopologyKeyword(line)
      if keyword == 'ATOM':
        atomName, types[atomName], charges[atomName] = self.parseTopologyAtomLine(line)
      elif keyword == 'BOND':
        bonds.append(tokens[:2])
    return types, charges, bonds

  def getNonPolarHydrogens(self, types, bonds):
    droppedH = {}
    for atoms in bonds:
      for hAtom, cAtom in (atoms, atoms[::-1]):
        if types[hAtom].startswith('h') and types[cAtom].startswith('c'):
          droppedH[hAtom] = cAtom
    return droppedH

  def hasTopologyAtom(self, line, atomNames):
    keyword, tokens = self.getTopologyKeyword(line)
    if keyword == 'ATOM':
      return tokens[0] in atomNames
    return keyword in TOPOLOGY_TERMS and any(atom in atomNames for atom in tokens)

  def setTopologyCharge(self, line, charges):
    keyword, tokens = self.getTopologyKeyword(line)
    if keyword != 'ATOM':
      return line
    return re.sub(r'CHARGE=\s*\S+', f'CHARGE= {charges[tokens[0]]:9.4f}', line)

  def getLigandCharge(self, num, mol2File):
    userCharge = getattr(self, f'ligCharge_{num}').get()
    if userCharge and userCharge.strip():
      return int(userCharge)

    charge, inAtoms = 0.0, False
    with open(mol2File) as f:
      for line in f:
        if line.startswith('@<TRIPOS>'):
          inAtoms = line.strip() == '@<TRIPOS>ATOM'
        elif inAtoms and line.strip():
          charge += float(line.split()[8])
    return round(charge)

  def joinLigandTopologies(self):
    topFiles, parFiles = [], []
    for num in range(1, 3):
      if self.getTopSource(num) == USER_TOP:
        topFile, parFile = getattr(self, f'ligTopFile_{num}').get(), getattr(self, f'ligParFile_{num}').get()
      else:
        topFile, parFile = self.getLigandTopFile(num), self.getLigandParFile(num)

      if topFile and parFile and os.path.exists(topFile) and os.path.exists(parFile):
        topFiles.append(topFile)
        parFiles.append(parFile)

    for files, oFile in [(topFiles, self.getJoinedTopFile()), (parFiles, self.getJoinedParFile())]:
      if files:
        with open(oFile, 'w') as fOut:
          for iFile in files:
            with open(iFile) as fIn:
              fOut.write(fIn.read().rstrip() + '\n\n')
  
  def restraintChains(self):
    oStr = ''
    for num in range(1,3):
      mFile = self.getMergedFile(num)
      if self.hasBackboneAtoms(mFile):
        oStr += haddockPlugin.runHaddock(self, 'haddock3-restraints restrain_bodies ', mFile,
                                        cwd=self._getExtraPath(), popen=True).strip() + '\n'

    with open(self.getChainRestrainstsFile(), 'w') as f:
      f.write(oStr)

  def hasBackboneAtoms(self, pdbFile):
    '''haddock3-restraints restrain_bodies only accepts structures with protein (CA) or
    nucleic acid (P) backbone atoms. Ligand-like inputs (e.g. glycans such as LPS, made
    entirely of HETATM records) have none, so that step must be skipped for them, as done
    in the HADDOCK3 protein-glycan docking example, which never sets unambig_fname.'''
    with open(pdbFile) as f:
      for line in f:
        if line.startswith('ATOM') and line[12:16].strip() in ('CA', 'P'):
          return True
    return False
  
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
    ligLines = self.getLigandConfLines()
    moduleLines = self.getRestraintConfLines(ambRestrFile, unambRestrFile) + ligLines
    if self.workflowType.get() == 0:
        confStr += DEF_PROTOCOL.format(ligand=ligLines, restraints=moduleLines, sampling=self.sampling.get(),
                                       select=self.selectModel.get(), tolerance=self.tolerance.get(),
                                       topModels=self.selectCluster.get())
    else:
        stepsList = self.parseManualWorkflow()
        confStr += self.buildConfStr(stepsList, ligLines, moduleLines)

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

  def getRestraintConfLines(self, ambRestrFile, unambRestrFile):
      s = f'ambig_fname = "{ambRestrFile}"\n'
      if os.path.exists(unambRestrFile):
          with open(unambRestrFile) as f:
              if f.read().strip():
                  s += f'unambig_fname = "{unambRestrFile}"\n'
      return s

  def getLigandConfLines(self):
      topFile, parFile = self.getJoinedTopFile(), self.getJoinedParFile()
      if os.path.exists(topFile) and os.path.exists(parFile):
          return f'ligand_top_fname = "{topFile}"\nligand_param_fname = "{parFile}"\n'
      return ''

  def buildConfStr(self, stepsList, ligLines, moduleLines):
      s = f'[topoaa]\n{ligLines}\n'
      h = self._haddockSteps
      for stepName, stepDic in stepsList:
          stepName = stepName.strip()
          if stepName == h[RIGDOC]:
              s += f'[rigidbody]\n{moduleLines}'
          elif stepName == h[FLEXREF]:
              s += f'[flexref]\n{moduleLines}'
          elif stepName == h[EMMIN]:
              s += f'[emref]\n{moduleLines}'
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
    mapDic, ligIds = self.buildMapDic(num)
    pasStr = getattr(self, f'passiveRestr_{num}').get() or ''

    act = self.mapResidueIds(self.getDefinedResidueIds(num, 'active'), mapDic, num)
    pas = self.mapResidueIds(self.getDefinedResidueIds(num, 'passive'), mapDic, num)

    mergedFile = self.getMergedFile(num)
    polymerAct = [resId for resId in act if resId not in ligIds]
    if FROM_ACTIVE in pasStr and polymerAct and self.hasBackboneAtoms(mergedFile):
        segId = self.getSegId(num)
        args = f'{mergedFile} -c {segId} {",".join(polymerAct)}'
        nPas = haddockPlugin.runHaddock(self, 'haddock3-restraints passive_from_active ', args,
                                        cwd=self._getExtraPath(), popen=True).split()
        pas = list(set(pas).union(set(nPas))) if pas else nPas

    pas = list(set(pas).difference(set(act)))
    return {'active': act, 'passive': pas}

  def mapResidueIds(self, resIds, mapDic, num):
    missing = [s for s in resIds if s not in mapDic]
    if missing:
      raise ValueError(f'Restraint residues for input {num} not found in its structure: {", ".join(missing)}')
    return list(dict.fromkeys(mapDic[s] for s in resIds))

  def getDefinedResidueIds(self, num, kind):
    resStr = getattr(self, f'{kind}Restr_{num}').get() or ''
    resIds = [s.strip() for line in resStr.split('\n') for s in line.split(',')]
    return [s for s in resIds if s and s != FROM_ACTIVE]

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

  def getTopSource(self, num):
    return getattr(self, f'ligTopSource_{num}').get()

  def getAtomLines(self, pdbFile):
    with open(pdbFile) as f:
      return [line for line in f if line.startswith(('ATOM', 'HETATM'))]

  def getLineKey(self, line):
    return f'{line[21]}|{line[22:26].strip()}|{line[26].strip()}|{line[17:20].strip()}'

  def getResidueKey(self, res):
    _, resNum, iCode = res.get_id()
    return f'{res.get_parent().get_id()}|{resNum}|{iCode.strip()}|{res.get_resname().strip()}'

  def getLineElement(self, line):
    element = line[76:78].strip() if len(line) > 77 else ''
    return element.upper() if element else line[12:16].strip().lstrip('0123456789')[:1].upper()

  def getLigandInfo(self, num):
    with open(self.getLigandInfoFile()) as f:
      return json.load(f)[str(num)]

  def calculateFreeSASA(self, num):
    inFile = self.getMergedFile(num)
    oFile = self._getExtraPath(f'input_{num}_SASA.txt')
    args = f'{os.path.abspath(inFile)} {os.path.abspath(oFile)}'
    haddockPlugin.runScript(self, 'calculateFreeSASA.py', args, HADDOCK_DIC, cwd=self._getExtraPath())
    return oFile

  def filterSASA(self, defRes, num):
    if not self.hasBackboneAtoms(self.getMergedFile(num)):
      return {k: sorted(v) for k, v in defRes.items()}

    segId, ligIds = self.getSegId(num), self.buildMapDic(num)[1]
    sasaDic = self.parseSASAFile(self.calculateFreeSASA(num))
    minSASA = self.filterSurface.get()

    nDefRes = {}
    for k in defRes:
      nDefRes[k] = []
      for res in defRes[k]:
        if res.strip():
          resId = f'{segId}_{res}'
          if res in ligIds or resId not in sasaDic or self.isAccessible(sasaDic[resId], minSASA):
            nDefRes[k].append(res)

    for k in nDefRes:
      nDefRes[k].sort()
    return nDefRes

  def isAccessible(self, sasaValues, minSASA):
    return any(v != 'nan' and float(v) >= minSASA for v in sasaValues)

  def buildMapDic(self, num):
    '''Returns a dictionary with the mapped residues from the original to the merged form
    {"ori_Id": "mergeId"} and the set of merged ids corresponding to non-standard residues
    '''
    ligInfo = self.getLigandInfo(num)
    ligKeys, unsupported = set(ligInfo['ligandKeys']), set(ligInfo['unsupported'])
    removedKeys = set(ligInfo.get('removedKeys', []))

    resNumDic, ligOriIds, ligIds, prevKey, c = {}, [], set(), None, 0
    for line in self.getAtomLines(self.getPDBFile(num)):
      key, oriId = self.getLineKey(line), f'{line[21]}_{line[22:26].strip()}'
      if key in removedKeys:
        continue
      if key in ligKeys:
        ligOriIds.append(oriId)
        continue
      if key != prevKey:
        c, prevKey = c + 1, key
      resNumDic.setdefault(oriId, str(c))
      if line[17:20].strip() in unsupported:
        ligIds.add(resNumDic[oriId])

    if ligOriIds:
      ligIds.add(str(c + 1))
      for oriId in ligOriIds:
        resNumDic[oriId] = str(c + 1)
    return resNumDic, ligIds


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

  def getLigandDir(self, num):
    return os.path.abspath(self._getExtraPath(f'ligand_{num}'))

  def getNoLigandFile(self, num):
    return os.path.join(self.getLigandDir(num), 'nonLigand.pdb')

  def getLigandPDBFile(self, num):
    return os.path.join(self.getLigandDir(num), f'{LIGAND_RESNAMES[num - 1]}.pdb')

  def getLigandTopFile(self, num):
    return os.path.join(self.getLigandDir(num), f'{LIGAND_RESNAMES[num - 1]}.top')

  def getLigandParFile(self, num):
    return os.path.join(self.getLigandDir(num), f'{LIGAND_RESNAMES[num - 1]}.param')

  def getLigandInfoFile(self):
    return os.path.abspath(self._getExtraPath('ligandInfo.json'))

  def getJoinedTopFile(self):
    return os.path.abspath(self._getExtraPath('ligands.top'))

  def getJoinedParFile(self):
    return os.path.abspath(self._getExtraPath('ligands.param'))

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
    oDirs = [oDir for oDir in oDirs if os.path.isdir(oDir)]
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
    actDefined = [bool(self.getDefinedResidueIds(num, 'active')) for num in range(1, 3)]
    if not any(actDefined):
      errors.append('No active restraint residues defined for any input. Define the residues of interest '
                    '(Structural ROI or Manual chain + residues) in the "Ambiguous restraints" tab and press the '
                    '"Add defined restraints" wizard so they appear in the "Active restraints summary".')

    for num in range(1, 3):
      if not actDefined[num - 1] and not self.getDefinedResidueIds(num, 'passive'):
        errors.append(f'Input {num} has neither active nor passive restraint residues. HADDOCK needs at least '
                      f'passive residues on the partner of an active one. Note that "{FROM_ACTIVE}" only works '
                      f'if active residues are defined for that same input.')

      if self.getTopSource(num) == USER_TOP:
        topFile, parFile = getattr(self, f'ligTopFile_{num}').get(), getattr(self, f'ligParFile_{num}').get()
        if not (topFile and parFile and os.path.exists(topFile) and os.path.exists(parFile)):
          errors.append(f'Input {num}: both the CNS topology and parameter files must be provided and exist when '
                        f'using "User files" for the non-standard residues topology.')
    return errors

