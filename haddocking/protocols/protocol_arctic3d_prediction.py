# **************************************************************************
# *
# * Authors:     Judith Maestro (scipion@cnb.csic.es)
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
from glob import glob
from collections import Counter
from Bio.PDB.PDBParser import PDBParser

from pyworkflow.protocol import params

from pwchem.protocols import ProtDefineStructROIs
from pwchem.utils import createPocketFile, getBaseName, pdbFromASFile, cifFromASFile, cleanPDB, getChainIds, \
  parseResidueIndexes, groupResidueIndexes
from pwchem.objects import StructROI, SetOfStructROIs

from .. import Plugin as haddockPlugin
from ..constants import ARCTIC3D_OUTPUT_FILES


class ProtArctic3DPrediction(ProtDefineStructROIs):
  """Run protein-protein interface predictions using ARCTIC-3D software from Bonvin Lab. It retrieves the
  interfaces of a protein (UniProt ID) known in the PDBe-KB and groups them into binding surfaces, which are defined
  as structural regions of interest on the input structure.

  More information:
  [[https://www.bonvinlab.org/arctic3d/][ARCTIC-3D documentation]]
  [[https://www.bonvinlab.org/arctic3d/webserver.html][ARCTIC-3D web server manual]]"""
  _label = 'ARCTIC-3D protein interface prediction'
  _references = ['Giulini2024']

  def _defineParams(self, form):
    form.addSection(label='Input')
    iGroup = form.addGroup('Input')
    iGroup.addParam('inputAtomStruct', params.PointerParam, pointerClass="AtomStruct", label='Input atom struct: ', allowsNull=True,
                    help="Input AtomStruct where the binding surfaces will be defined. If it is left empty, "
                         "the reference structure selected by ARCTIC-3D (UniProt numbering) is used.")
    iGroup.addParam('chain_name', params.StringParam, allowsNull=True, label='Chain of interest',
                    help='Specify the chain of the structure of interest. Required only if an input atom '
                         'struct is given.')
    iGroup.addParam('restrictRegion', params.BooleanParam, default=True, condition='inputAtomStruct is not None',
                    label='Restrict the prediction to the chain of interest: ',
                    help='If set, the interfaces retrieved from PDBe-KB are trimmed to the region covered by the chain '
                         'of interest and the clustering is performed on that chain alone, so the binding surfaces are '
                         'specific to it (for example, the RBD of a Spike protein instead of the whole protein).\n'
                         'If not set, ARCTIC-3D clusters the interfaces of the whole UniProt entry on the reference '
                         'structure it selects, and only the resulting surfaces that fall in the chain are kept.')
    iGroup.addParam('resRange', params.StringParam, default='',
                    condition='inputAtomStruct is not None and restrictRegion',
                    label='Residues of interest: ',
                    help='Region of the chain the prediction is restricted to, as residue numbers of the input atom '
                         'struct (for example 333-526 or 333-400,450-526). If empty, all the residues present in the '
                         'chain of interest are used.')
    iGroup.addParam('uniprotId', params.StringParam, allowsNull=False, label='UniProt ID: ',
                    help='UniProt accession of the protein (e.g. P68871). If you do not know it, search UniProt for the '
                         'name of your protein: for Hemoglobin, the first result is Hemoglobin subunit beta, with the '
                         'accession (Entry) P68871. Entering the UniProt ID is enough to run the protocol; the '
                         'other parameters are optional. ARCTIC-3D retrieves from PDBe-KB all the known interfaces '
                         'of this protein.')

    rGroup = form.addGroup('Interface retrieval')
    rGroup.addParam('ligand', params.EnumParam, default=0, choices=['no', 'yes', 'both'],
                    label='Retrieve ligand binding residues: ',
                    help='Whether to retrieve the ligand binding residues instead of (yes) or together with (both) '
                         'the protein-protein interfaces.')
    rGroup.addParam('full', params.BooleanParam, default=False, label='Use full UniProt-PDB-chain information: ',
                    help='Consider the full UniProt-PDB-chain information in the retrieval.')
    rGroup.addParam('outPartner', params.StringParam, default='', label='Partners to exclude: ',
                    help='Comma-separated UniProt IDs of the partners to exclude from the search.')
    rGroup.addParam('outPdb', params.StringParam, default='', label='PDBs to exclude: ',
                    help='PDB IDs to be excluded: PDB entries to leave out of the search, comma-separated (for example '
                         '7pcq,5jdo): because you do not trust an interface, or to leave it out for benchmarking.')

    cGroup = form.addGroup('Interface clustering')
    cGroup.addParam('clustThreshold', params.FloatParam, default=0.866, label='Clustering threshold: ',
                    help='Interfaces that merge below the threshold (default: 0.866, corresponding to a 60° angle) form a '
                         'single binding surface. The threshold can be adjusted with --threshold to obtain finer or '
                         'coarser clustering.')
    cGroup.addParam('linkage', params.EnumParam, default=0,
                    choices=['average', 'single', 'complete', 'weighted', 'centroid', 'median', 'ward'],
                    label='Linkage strategy: ', help='Linkage strategy (default average): how the distance between two clusters is computed. With '
                         'average, it is the average distance between their interfaces.')
    cGroup.addParam('minClustSize', params.IntParam, default=0, label='Minimum cluster size: ',
                    help='Minimum number of residues in a cluster to be considered.')
    cGroup.addParam('intCovCutoff', params.FloatParam, default=0.7, label='Interface coverage cutoff: ',
                    help='Interface coverage cutoff (default 0.7): the fraction of the residues of an interface that must '
                         'be present in the structure for the interface to be kept. With 0.5, at least half of the '
                         'residues of every interface must be in the structure.')

    oGroup = form.addGroup('Regions definition')
    oGroup.addParam('threshold', params.FloatParam, default=0.0, label='Residue probability threshold: ',
                    help='Minimum contact probability (0-1) of a residue in its cluster to be included in the '
                         'region. Residues with probability 1 are the hotspots of the binding surface.')


  def _insertAllSteps(self):
    self._insertFunctionStep(self.predictStep)
    if self.restrictsToRegion():
      self._insertFunctionStep(self.restrictedPredictStep)
    self._insertFunctionStep(self.definePocketsStep)
    self._insertFunctionStep(self.defineOutputStep)

  def predictStep(self):
    if self.hasInputStruct():
      pdbFromASFile(self.getInputPath(), self._getPdbFile(), atomStruct=self.inputAtomStruct.get())

    runDir = self.getRetrievalDir() if self.restrictsToRegion() else self.getOutputDir()
    args = f'{self.uniprotId.get().strip()} --run_dir {runDir} ' \
           f'--ligand {self.getEnumText("ligand")} --threshold {self.clustThreshold.get()} ' \
           f'--linkage_strategy {self.getEnumText("linkage")} --min_clust_size {self.minClustSize.get()} ' \
           f'--int_cov_cutoff {self.intCovCutoff.get()} --no_report'
    if self.full.get():
      args += ' --full'
    if self.outPartner.get().strip():
      args += f' --out_partner {self.outPartner.get().strip().replace(" ", "")}'
    if self.outPdb.get().strip():
      args += f' --out_pdb {self.outPdb.get().strip().replace(" ", "")}'

    haddockPlugin.runARCTIC3D(self, args, cwd=os.path.abspath(self._getExtraPath()))

  def restrictedPredictStep(self):
    """Cluster, on the chain of interest alone, the interfaces retrieved for the whole UniProt entry, trimmed to
    the region of interest and renumbered to the input structure numbering"""
    chainRes = self.getChainResidues()
    self.parseClusters(outDir=self.getRetrievalDir())
    offset = self.getNumberingOffset(chainRes)
    self.info(f'Offset between the UniProt numbering and the input structure numbering: {offset}')

    regionRes = self.getRegionResidues(chainRes)
    if regionRes and regionRes[0] < 1:
      raise Exception(f'The residues of the chain of interest are numbered from {regionRes[0]}, and ARCTIC-3D only '
                      f'admits positive residue numbers. Renumber the input structure or deactivate the restriction '
                      f'to the chain.')
    with open(self.getRegionFile(), 'w') as f:
      f.write(groupResidueIndexes(regionRes))
    self.info(f'Restricting the prediction to {self.getRegionDescription()}')

    chainFile = self.getChainFile()
    cleanPDB(self._getPdbFile(), chainFile, waters=True, hetatm=True, chainIds=self.getChainIds(), remAltloc=True)
    intFile = self.writeRegionInterfaces(regionRes, offset)

    args = f'{chainFile} --interface_file {intFile} --run_dir {self.getOutputDir()} ' \
           f'--threshold {self.clustThreshold.get()} --linkage_strategy {self.getEnumText("linkage")} ' \
           f'--min_clust_size {self.minClustSize.get()} --int_cov_cutoff {self.intCovCutoff.get()} --no_report'
    haddockPlugin.runARCTIC3D(self, args, cwd=os.path.abspath(self._getExtraPath()))

  def definePocketsStep(self):
    self.coordsClusters, self.resClusters = [], []
    clusters = self.parseClusters(self.threshold.get())
    chainRes = self.getChainResidues()
    offset = self.getNumberingOffset(chainRes) if self.hasInputStruct() and not self.restrictsToRegion() else 0
    self.info(f'Offset between the ARCTIC-3D numbering and the input structure numbering: {offset}')

    for clustId in sorted(clusters):
      coords, resIds = [], []
      for resId in clusters[clustId]:
        residue = chainRes.get(resId - offset)
        if residue is not None:
          coords += [list(a.get_coord()) for a in residue.get_atoms()]
          resIds.append(f'{residue.get_parent().id}_{residue.get_id()[1]}')
      if coords:
        self.coordsClusters.append(coords)
        self.resClusters.append(sorted(resIds, key=lambda x: int(x.split('_')[-1])))
      else:
        self.info(f'Cluster {clustId} could not be mapped to the input structure')

  def defineOutputStep(self):
    inpStruct = self.inputAtomStruct.get()
    cifFile = cifFromASFile(self.getStructFile(), self._getCifFile(), inpStruct)

    outPockets = SetOfStructROIs(filename=self._getPath('StructROIs.sqlite'))
    for i, (clust, resIds) in enumerate(zip(self.coordsClusters, self.resClusters)):
      pocketFile = self._getExtraPath(f'pocketFile_{i}.cif')
      createPocketFile(clust, i, pocketFile)
      pocket = StructROI(pocketFile, cifFile)
      if str(type(inpStruct).__name__) == 'SchrodingerAtomStruct':
        pocket._maeFile = params.String(os.path.relpath(inpStruct.getFileName()))
      self.setClusterContacts(pocket, resIds)
      outPockets.append(pocket)

    if len(outPockets) > 0:
      outPockets.buildPDBhetatmFile()
      self._defineOutputs(outputStructROIs=outPockets)


  ################# UTIL FUNCTIONS ##############################

  def setClusterContacts(self, pocket, resIds):
    """Set the contact residues of the structural region to the residues of the ARCTIC-3D cluster, and the contact
    atoms to their atoms, instead of every residue at less than 5A from them"""
    cAtoms = [atom for atom in pocket.buildContactAtoms(calculate=True) if pocket.buildResidueId(atom) in resIds]
    pocket.setContactAtoms(pocket.encodeIds(pocket.getAtomsIds(cAtoms)))
    pocket.setContactResidues(pocket.encodeIds(resIds))

  def parseClusters(self, threshold=0.0, outDir=None):
    """Return {clusterId: [residue ids]} from the ARCTIC-3D residues probabilities file, keeping the
    residues with probability >= threshold. The residue names are stored in self.resNames"""
    self.resNames = {}
    clusters, clustId = {}, None
    with open(self.getOutputFile('clustered_residues_probs.out', outDir)) as f:
      for line in f:
        words = line.split()
        if line.startswith('Cluster'):
          clustId = int(words[1])
          clusters[clustId] = []
        elif clustId is not None and len(words) == 4 and words[0].isdigit():
          resId, resName, prob = int(words[1]), words[2], float(words[3])
          self.resNames[resId] = resName
          if prob >= threshold:
            clusters[clustId].append(resId)
    return {cId: resIds for cId, resIds in clusters.items() if resIds}

  def getChainResidues(self):
    structModel = PDBParser(QUIET=True).get_structure('receptor', self.getStructFile())[0]
    if self.hasInputStruct():
      chain = structModel[self.getChainIds()[0]]
    else:
      chain = next(structModel.get_chains())
    return {res.get_id()[1]: res for res in chain.get_residues() if res.get_id()[0] == ' '}

  def getChainIds(self):
    return getChainIds(self.chain_name.get())

  def getRegionDescription(self):
    """Chain and residues the prediction was restricted to, as written in the log and the summary"""
    descr = f'chain {",".join(self.getChainIds())}'
    if os.path.exists(self.getRegionFile()):
      with open(self.getRegionFile()) as f:
        descr += f' and residues {f.read().strip()}'
    return descr

  def getRegionResidues(self, chainRes):
    """Residue ids (input structure numbering) the prediction is restricted to: those of the chain of interest,
    narrowed down to the residue range defined by the user, if any"""
    regionRes = set(chainRes)
    if self.resRange.get().strip():
      regionRes = regionRes & set(parseResidueIndexes(self.resRange.get().strip()))
    return sorted(regionRes)

  def writeRegionInterfaces(self, regionRes, offset):
    """Write the ARCTIC-3D interface file with the retrieved interfaces trimmed to the region of interest and
    renumbered from the UniProt numbering to the input structure numbering"""
    regionRes, intFile = set(regionRes), os.path.abspath(self._getExtraPath('regionInterfaces.txt'))
    nInts = 0
    with open(self.getOutputFile('retrieved_interfaces.out', self.getRetrievalDir())) as fIn, \
            open(intFile, 'w') as fOut:
      for i, line in enumerate(fIn):
        name, resIds = line.split('->')
        resIds = [resId for resId in map(lambda x: int(x) - offset, resIds.split()) if resId in regionRes]
        if resIds:
          nInts += 1
          fOut.write(f'{name.split()[1]}_{i} {",".join(map(str, resIds))}\n')

    self.info(f'{nInts} of the retrieved interfaces have residues in the region of interest')
    if nInts < 2:
      raise Exception('Less than two of the interfaces retrieved by ARCTIC-3D reach the region of interest, so they '
                      'cannot be clustered. Try a larger region or deactivate the restriction to the chain.')
    return intFile

  def hasInputStruct(self):
    return self.inputAtomStruct.get() is not None

  def restrictsToRegion(self):
    return self.hasInputStruct() and self.restrictRegion.get()

  def getReferenceFile(self):
    uniprotId = self.uniprotId.get().strip()
    for file in sorted(glob(self.getOutputFile(f'{uniprotId}-*-*.pdb'))):
      if '_cl' not in os.path.basename(file):
        return file

  def getStructFile(self):
    return self._getPdbFile() if self.hasInputStruct() else self.getReferenceFile()

  def getNumberingOffset(self, chainRes):
    """Offset (UniProt id - structure id) with more residue names in agreement"""
    counts = Counter()
    for resId, resName in self.resNames.items():
      for strId, res in chainRes.items():
        if res.get_resname() == resName:
          counts[resId - strId] += 1
    return counts.most_common(1)[0][0] if counts else 0

  def getInputPath(self):
    return self.inputAtomStruct.get().getFileName()

  def _getPdbFile(self):
    return os.path.abspath(self._getExtraPath(self._getInputName() + '.pdb'))

  def _getCifFile(self):
    return os.path.abspath(self._getExtraPath(self._getInputName() + '.cif'))

  def _getInputName(self):
    if self.hasInputStruct():
      return getBaseName(self.getInputPath())
    return getBaseName(self.getReferenceFile())

  def getRegionFile(self):
    return self._getExtraPath('regionResidues.txt')

  def getChainFile(self):
    return os.path.abspath(self._getExtraPath(f'{self._getInputName()}_{"".join(self.getChainIds())}.pdb'))

  def getOutputDir(self):
    return os.path.abspath(self._getPath('arctic3dOutput'))

  def getRetrievalDir(self):
    return os.path.abspath(self._getExtraPath('arctic3dRetrieval'))

  def getOutputFile(self, fileName, outDir=None):
    return os.path.join(outDir if outDir else self.getOutputDir(), fileName)

  def _summary(self):
    summary = []
    outDir = self.getOutputDir()
    if os.path.exists(outDir) and self.restrictsToRegion():
      summary.append(f'The prediction was restricted to {self.getRegionDescription()} of the input structure: '
                     f'the interfaces retrieved for the whole UniProt entry (in {self.getRetrievalDir()}) were '
                     f'trimmed to that region and clustered on the chain alone.')
    if os.path.exists(outDir):
      summary.append(f'In: {outDir}\nyou can find the ARCTIC-3D results '
                     f'([[https://www.bonvinlab.org/arctic3d/][documentation]]):')
      for fileName, descr in ARCTIC3D_OUTPUT_FILES:
        if glob(os.path.join(outDir, fileName)):
          summary.append(f'- {fileName}: {descr}')
      summary.append('Residue probability: fraction of the interfaces of a cluster where the residue is observed '
                     '(1.0 = hotspot residue, present in every interface of the binding surface).')
      if hasattr(self, 'outputStructROIs'):
        summary.append(f'Binding surfaces defined as structural regions: {len(self.outputStructROIs)}')
    return summary

  def _validate(self):
    errors = []
    if not self.uniprotId.get() or not self.uniprotId.get().strip():
      errors.append('A UniProt ID is required')
    if self.hasInputStruct() and not self.chain_name.get():
      errors.append('The chain of interest is required when an input atom struct is given')
    if self.restrictsToRegion() and self.resRange.get().strip():
      try:
        parseResidueIndexes(self.resRange.get().strip())
      except ValueError:
        errors.append('The residues of interest must be numbers, ranges or a comma-separated combination of them, '
                      'such as 333-526 or 333-400,450-526')
    return errors
