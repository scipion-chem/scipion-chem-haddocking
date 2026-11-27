# **************************************************************************

from pyworkflow.tests import BaseTest, setupTestProject, DataSet
from pwem.protocols import ProtImportPdb
from ..protocols import ProtArctic3DPrediction

class TestArctic3D(BaseTest):
    @classmethod
    def setUpClass(cls):
        cls.ds = DataSet.getDataSet('model_building_tutorial')

        setupTestProject(cls)
        cls._runImportPDB()

    @classmethod
    def _runImportPDB(cls):
        protImportPDB = cls.newProtocol(
            ProtImportPdb,
            inputPdbData=1,
            pdbFile=cls.ds.getFile('PDBx_mmCIF/5ni1.pdb'))
        cls.launchProtocol(protImportPDB)
        cls.protImportPDB = protImportPDB

    def _runArctic3D(self, **kwargs):
        protArctic = self.newProtocol(
            ProtArctic3DPrediction,
            inputAtomStruct=self.protImportPDB.outputPdb,
            chain_name='{"model": 0, "chain": "C", "residues": 141}',
            uniprotId='P69905', **kwargs)

        self.launchProtocol(protArctic)
        return protArctic

    def testArctic3D(self):
        protArctic = self._runArctic3D(restrictRegion=False)
        self.assertIsNotNone(protArctic.outputStructROIs)

    def testArctic3DRestricted(self):
        protArctic = self._runArctic3D()
        self.assertIsNotNone(protArctic.outputStructROIs)

    def testArctic3DResidueRange(self):
        protArctic = self._runArctic3D(resRange='20-120')
        self.assertIsNotNone(protArctic.outputStructROIs)

        with open(protArctic.getOutputFile('clustered_residues.out')) as f:
            resIdxs = [int(idx) for line in f for idx in line.split('->')[1].split()]
        self.assertTrue(resIdxs and all([20 <= idx <= 120 for idx in resIdxs]))

    def testArctic3DOnlyUniprot(self):
        protArctic = self.newProtocol(ProtArctic3DPrediction, uniprotId='P69905')
        self.launchProtocol(protArctic)
        self.assertIsNotNone(protArctic.outputStructROIs)
