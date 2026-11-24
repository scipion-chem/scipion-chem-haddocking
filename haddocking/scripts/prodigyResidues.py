# **************************************************************************
# *
# * Authors:     Judith Maestro (scipionchem@cnb.csic.es)
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

'''Classifies PRODIGY intermolecular contacts and computes per-residue relative SASA
of the complex (NIS analysis) using the PRODIGY package itself.
Usage: python prodigyResidues.py <structure> <contacts.ic> <accThreshold> <output.json>'''

import json
import sys

from prodigy_prot.modules.parsers import parse_structure
from prodigy_prot.modules.freesasa_tools import execute_freesasa_api
from prodigy_prot.modules.prodigy import analyse_nis
from prodigy_prot.modules.aa_properties import aa_character_ic, aa_character_protorp, aa_character_hydro


def readContacts(icFile):
  contacts = []
  with open(icFile) as f:
    for line in f:
      fields = line.split()
      if len(fields) == 6:
        contacts.append(((fields[2], fields[0], fields[1]), (fields[5], fields[3], fields[4])))
  return contacts


def pairType(res1, res2, table):
  c1, c2 = table.get(res1[1]), table.get(res2[1])
  if c1 is None or c2 is None:
    return None
  return ''.join(sorted((c1, c2)))


if __name__ == '__main__':
  structFile, icFile, accThreshold, outFile = sys.argv[1], sys.argv[2], float(sys.argv[3]), sys.argv[4]

  contacts = readContacts(icFile)
  bins = dict.fromkeys(
    ['CC', 'CP', 'AC', 'PP', 'AP', 'AA', 'LL', 'BL', 'BB'], 0
  )
  contactRows, nContacts = [], {}
  for res1, res2 in contacts:
    icType, hydroType = pairType(res1, res2, aa_character_ic), pairType(res1, res2, aa_character_hydro)
    for t in (icType, hydroType):
      if t in bins:
        bins[t] += 1
    for res in (res1, res2):
      nContacts[res] = nContacts.get(res, 0) + 1
    contactRows.append({'res1': list(res1), 'char1': aa_character_ic.get(res1[1]),
                        'res2': list(res2), 'char2': aa_character_ic.get(res2[1]),
                        'icType': icType, 'hydroType': hydroType})

  models, _, _ = parse_structure(structFile)
  _, rsaDic = execute_freesasa_api(models[0])
  nisA, nisC, nisP = analyse_nis(rsaDic, acc_threshold=accThreshold)

  residueRows = []
  for (chain, resName, resNum), rsa in rsaDic.items():
    key = (chain, resName, str(resNum).strip())
    residueRows.append({'chain': chain, 'resName': resName, 'resNum': key[2],
                        'icCharacter': aa_character_ic.get(resName),
                        'nisCharacter': aa_character_protorp.get(resName), 'rsa': rsa,
                        'surface': rsa >= accThreshold, 'nContacts': nContacts.get(key, 0)})

  with open(outFile, 'w') as f:
    json.dump({'accThreshold': accThreshold, 'bins': bins,
               'nis': {'A': nisA, 'C': nisC, 'P': nisP},
               'contacts': contactRows, 'residues': residueRows}, f, indent=2)
