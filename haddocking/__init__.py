# **************************************************************************
# *
# * Authors:	Daniel Del Hoyo Gomez (ddelhoyo@cnb.csic.es)
# *
# * Unidad de Bioinformatica of Centro Nacional de Biotecnologia, CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
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
# * 02111-1307  USA
# *
# * All comments concerning this program package may be sent to the
# * e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************
"""
This package contains protocols for creating and using IIITD Raghava software
"""

import os
import subprocess as sp

from scipion.install.funcs import InstallHelper

from pwchem import Plugin as pwchemPlugin

from .constants import *

# Pluging variables
_logo = 'haddocking_logo.jpeg'


class Plugin(pwchemPlugin):
	"""
	"""

	@classmethod
	def _defineVariables(cls):
		cls._defineEmVar(CPORT_DIC['home'], cls.getEnvName(CPORT_DIC))

		cls._defineVar('CPORT_ENV_ACTIVATION', cls.getEnvActivationCommand(CPORT_DIC))


	@classmethod
	def defineBinaries(cls, env, default=True):
		"""This function defines the binaries for each package."""
		cls.addCPORTPackage(env)

	@classmethod
	def addCPORTPackage(cls, env, default=True):
		installer = InstallHelper(CPORT_DIC['name'], packageHome=cls.getVar(CPORT_DIC['home']),
															packageVersion=CPORT_DIC['version'])
		# Installing CPORT package
		installer.getCondaEnvCommand(pythonVersion='3.10', requirementsFile=False) \
			.addCommand(f'git clone {CPORT_URL} && cd cport && {cls.getEnvActivationCommand(CPORT_DIC)} && pip install .') \
			.addPackage(env, ['git', 'conda', 'pip'], default=default)



	# ---------------------------------- Protocol functions-----------------------
	@classmethod
	def runCPORT(cls, protocol, args, popen=False):
		""" Run CPORRT command from a given protocol. """
		cportHome, cportAct = os.path.join(cls.getVar(CPORT_DIC["home"]), 'cport'), cls.getEnvActivationCommand(CPORT_DIC)

		outPath = os.path.join(cportHome, 'output')
		if not os.path.exists(outPath):
			os.mkdir(outPath)

		fullProgram = f'{cportAct} && cport'
		if not popen:
			protocol.runJob(fullProgram, args, env=cls.getEnviron(), cwd=cportHome)
		else:
			sp.check_call(f'{fullProgram} {args}', cwd=cportHome, shell=True)

		os.rename(outPath, protocol.getOutputDir())
