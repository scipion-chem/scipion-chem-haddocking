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
import shutil
import subprocess as sp

from scipion.install.funcs import InstallHelper
from pwem import Config as emConfig

from pwchem import Plugin as pwchemPlugin

from .constants import *

# Pluging variables
_logo = 'haddocking_logo.jpeg'


class Plugin(pwchemPlugin):
	"""
	"""

	@classmethod
	def _defineVariables(cls):
		cls._defineEmVar(HADDOCK_DIC['home'], cls.getEnvName(HADDOCK_DIC))
		cls._defineEmVar(CPORT_DIC['home'], cls.getEnvName(CPORT_DIC))
		cls._defineEmVar(GDOCK_DIC['home'], cls.getEnvName(GDOCK_DIC))
		cls._defineEmVar(LIGHTDOCK_DIC['home'], cls.getEnvName(LIGHTDOCK_DIC))
		cls._defineEmVar(PRODIGY_DIC['home'], cls.getEnvName(PRODIGY_DIC))

		cls._defineVar(CNS_DIC['home'], None)
		cls._defineVar(CNS_DIC['tar'], None)

		cls._defineVar('HADDOCK_ENV_ACTIVATION', cls.getEnvActivationCommand(HADDOCK_DIC))
		cls._defineVar('PRODIGY_ACTIVATION_CMD', cls.getEnvActivationCommand(PRODIGY_DIC))
		cls._defineVar('CPORT_ENV_ACTIVATION', cls.getEnvActivationCommand(CPORT_DIC))


	@classmethod
	def defineBinaries(cls, env, default=True):
		"""This function defines the binaries for each package."""
		cls.addHADDOCKPackage(env, cnsHome=cls.getVar(CNS_DIC['home']))
		cls.addCPORTPackage(env)
		cls.addPRODIGYPackage(env)
		# cls.addGDockPackage(env)
		# cls.addLightDockPackage(env)

	@classmethod
	def addHADDOCKPackage(cls, env, cnsHome=None, tarPath=None, default=True):
		installer = InstallHelper(HADDOCK_DIC['name'], packageHome=cls.getVar(HADDOCK_DIC['home']),
															packageVersion=HADDOCK_DIC['version'])

		# Installing HADDOCK + CNS package
		installer.getCondaEnvCommand(pythonVersion='3.9', requirementsFile=False)\
			.addCommand(f'{cls.getEnvActivationCommand(HADDOCK_DIC)} && pip install haddock3', 'HADDOCK_INSTALLED') \
			.addPackage(env, ['git', 'pip'], default=default)

	@classmethod
	def addPRODIGYPackage(cls, env, cnsHome=None, tarPath=None, default=True):
		installer = InstallHelper(PRODIGY_DIC['name'], packageHome=cls.getVar(PRODIGY_DIC['home']),
		                          packageVersion=PRODIGY_DIC['version'])

		installer.getCondaEnvCommand(pythonVersion='3.11', requirementsFile=False) \
			.addCommand(f'{cls.getEnvActivationCommand(PRODIGY_DIC)} && pip install prodigy-prot', 'PRODIGY_INSTALLED') \
			.addPackage(env, ['git', 'pip'], default=default)

	@classmethod
	def addGDockPackage(cls, env, default=True):
		installer = InstallHelper(GDOCK_DIC['name'], packageHome=cls.getVar(GDOCK_DIC['home']),
															packageVersion=GDOCK_DIC['version'])

		installer.addCommand(f'git clone {GDOCK_DIC["url"]} && {cls.getEnvActivationCommand(HADDOCK_DIC)} && '
									f'pip install deap scipy mgzip biopython && {cls.addCondaVarCmd("GDOCK_HOME", "$(pwd)/gdock")}',
									'GDOCK_INSTALLED') \
			.addPackage(env, ['git', 'pip'], default=default)

	@classmethod
	def addLightDockPackage(cls, env, default=True):
		installer = InstallHelper(LIGHTDOCK_DIC['name'], packageHome=cls.getVar(LIGHTDOCK_DIC['home']),
															packageVersion=LIGHTDOCK_DIC['version'])

		installer.addCommand(f'git clone {LIGHTDOCK_DIC["url"]} && {cls.getEnvActivationCommand(HADDOCK_DIC)} && '
									f'cd lightdock && pip install -e . && {cls.addCondaVarCmd("LIGHTDOCK_HOME", "$(pwd)/lightdock")} && '
									f'{cls.addCondaVarCmd("PATH", "$PATH:$LIGHTDOCK_HOME/bin")} && '
									f'{cls.addCondaVarCmd("PYTHONPATH", "$PYTHONPATH:$LIGHTDOCK_HOME")}',
									'LIGHTDOCK_INSTALLED') \
			.addPackage(env, ['git', 'pip'], default=default)

	@classmethod
	def addCPORTPackage(cls, env, default=True):
		installer = InstallHelper(CPORT_DIC['name'], packageHome=cls.getVar(CPORT_DIC['home']),
															packageVersion=CPORT_DIC['version'])
		cportPatcher = os.path.join(os.path.dirname(__file__), 'scripts', 'patch_cport.py')

		# Installing CPORT package
		installer.getCondaEnvCommand(pythonVersion='3.10', requirementsFile=False) \
			.addCommand(f'git clone {CPORT_DIC["url"]} && python3 {cportPatcher} cport',
									'CPORT_PATCHED') \
			.addCommand(f'cd cport && {cls.getEnvActivationCommand(CPORT_DIC)} && pip install .',
									'CPORT_INSTALLED') \
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

		shutil.move(outPath, protocol.getOutputDir())

	@classmethod
	def runHaddock(cls, protocol, program, args, cwd, popen=False):
		""" Run HADDOCK commands from a given protocol. """
		fullProgram = f'{cls.getEnvActivationCommand(HADDOCK_DIC)} && {program} '
		if not popen:
			protocol.runJob(fullProgram, args, env=cls.getEnviron(), cwd=cwd)
		else:
			res = sp.check_output(f'{fullProgram} {args}', cwd=cwd, shell=True, text=True)
			return res


	@classmethod
	def checkVarPath(cls, softDic, var='home'):
		'''Check if a plugin variable exists and so do its path'''
		exists = False
		varValue = cls.getVar(softDic[var])
		if varValue and os.path.exists(varValue):
			exists = True
		return exists

	@classmethod
	def getDefaultDir(cls, softDic, fn=""):
		emDir = emConfig.EM_ROOT
		for file in os.listdir(emDir):
			if softDic['pattern'] in file.lower():
				foundDir = os.path.join(emDir, file, fn)
				return foundDir.rstrip('/')
		return os.path.join(emConfig.EM_ROOT, cls.getEnvName(softDic))

	@classmethod
	def addCondaVarCmd(cls, varName, varValue):
		return f"conda env config vars set {varName}={varValue}"

	@classmethod
	def getPluginHome(cls, path=""):
		import haddocking
		fnDir = os.path.split(haddocking.__file__)[0]
		return os.path.join(fnDir, path)
