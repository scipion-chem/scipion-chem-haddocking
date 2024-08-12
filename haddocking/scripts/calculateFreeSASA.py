# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors: Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
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
# *  e-mail address 'you@yourinstitution.email'
# *
# **************************************************************************

'''Calculates de SASA of a protein per residue using FreeSASA'''


import sys, freesasa

if __name__ == "__main__":
  inFile = sys.argv[1]
  sysName = inFile.split("/")[-1].split(".")[0]
  oFile = f'{sysName}_FreeSASA.txt' if len(sys.argv) < 3 else sys.argv[2]

  struct = freesasa.Structure(inFile)
  result = freesasa.calc(struct)

  resRes = result.residueAreas()
  with open(oFile, 'w') as f:
    for chain, chDic in resRes.items():
      for resNum, sasaObj in chDic.items():
        f.write(f'{chain}_{resNum} : {sasaObj.relativeMainChain}\t{sasaObj.relativeSideChain}\n')


