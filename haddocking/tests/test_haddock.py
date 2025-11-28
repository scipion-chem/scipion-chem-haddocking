# **************************************************************************
# *
# * Authors:     Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
# *
# * Unidad de Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 3 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307 USA
# *
# * All comments concerning this program package may be sent to the
# * e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

from pyworkflow.tests import BaseTest, setupTestProject, DataSet
from pwem.protocols import ProtImportPdb
from ..protocols import ProtHaddockPrediction

class TestHADDOCK(BaseTest):
    @classmethod
    def setUpClass(cls):
        cls.ds = DataSet.getDataSet('model_building_tutorial')
        setupTestProject(cls)

    @classmethod
    def _runImportPDB(cls):
        protImportPDB = cls.newProtocol(
            ProtImportPdb,
            inputPdbData=1,
            pdbFile=cls.ds.getFile('PDBx_mmCIF/1ake_start.pdb'))
        cls.launchProtocol(protImportPDB)
        return protImportPDB

    def _runHADDOCK(self, protImportPdb):
        protHAD = self.newProtocol(
            ProtHaddockPrediction, numberOfThreads=6,
            inputAtomStruct_1=protImportPdb.outputPdb, inputAtomStruct_2=protImportPdb.outputPdb,
            activeRestr_1='A_35,A_36,A_37,A_38,A_39,A_40', activeRestr_2='A_113,A_114,A_115,A_116,A_117',
            workflowType=0, sampling=5, selectModel=3, selectCluster=2
        )

        self.launchProtocol(protHAD)
        return protHAD

    def testHADDOCK(self):
        protImportPdb = self._runImportPDB()
        protHADDOC = self._runHADDOCK(protImportPdb)
        self.assertIsNotNone(protHADDOC.outputStructures)

