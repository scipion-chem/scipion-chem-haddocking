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

from pyworkflow.tests import setupTestProject, DataSet, BaseTest

from pwem.protocols import ProtImportSequence

from pwchem.utils import assertHandle

from ..protocols import ProtIIITDEpitopePrediction, ProtIIITDEvaluations
from ..constants import SELSUM, EVALSUM

class TestIIITDPrediction(BaseTest):
	NAME = 'USER_SEQ'
	DESCRIPTION = 'User description'
	AMINOACIDSSEQ1 = 'MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHG'

	@classmethod
	def setUpClass(cls):
		super().setUpClass()
		cls.ds = DataSet.getDataSet('model_building_tutorial')
		setupTestProject(cls)

		cls._runImportSeq()
		cls._waitOutput(cls.protImportSeq, 'outputSequences', sleepTime=5)

	@classmethod
	def _runImportSeq(cls):
		kwargs = {'inputSequenceName': cls.NAME,
							'inputSequenceDescription': cls.DESCRIPTION,
							'inputRawSequence': cls.AMINOACIDSSEQ1
							}

		cls.protImportSeq = cls.newProtocol(
			ProtImportSequence, **kwargs)
		cls.proj.launchProtocol(cls.protImportSeq, wait=False)

	def _runIIITDSelection(self):
		protSel = self.newProtocol(ProtIIITDEpitopePrediction,
																		inSels=SELSUM)

		protSel.inputSequence.set(self.protImportSeq)
		protSel.inputSequence.setExtended('outputSequence')

		self.proj.launchProtocol(protSel, wait=False)
		return protSel

	def test(self):
		protSel = self._runIIITDSelection()
		self._waitOutput(protSel, 'outputROIs', sleepTime=10)
		assertHandle(self.assertIsNotNone, getattr(protSel, 'outputROIs', None))


class TestIIITDEvaluation(TestIIITDPrediction):
	NAME = 'USER_SEQ'
	DESCRIPTION = 'User description'
	AMINOACIDSSEQ1 = 'MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHG'

	def _runIIITDEvaluation(self, protSel):
		protEval = self.newProtocol(ProtIIITDEvaluations,
															 inEvals=EVALSUM)

		protEval.inputROIs.set(protSel)
		protEval.inputROIs.setExtended('outputROIs')

		self.proj.launchProtocol(protEval, wait=False)
		return protEval

	def test(self):
		protSel = self._runIIITDSelection()
		self._waitOutput(protSel, 'outputROIs', sleepTime=10)
		protEval = self._runIIITDEvaluation(protSel)
		self._waitOutput(protEval, 'outputROIs', sleepTime=10)
		assertHandle(self.assertIsNotNone, getattr(protEval, 'outputROIs', None))