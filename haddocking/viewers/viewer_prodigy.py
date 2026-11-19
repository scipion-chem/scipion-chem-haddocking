import os
import csv
import glob
import json

from pyworkflow.viewer import ProtocolViewer, DESKTOP_TKINTER
from pyworkflow.protocol import params
from pyworkflow.gui.tree import TreeProvider
from pyworkflow.gui.dialog import ListDialog

from pwchem.viewers.viewers_data import PyMolViewer

from haddocking.constants import PRODIGY_IC_CHARACTER, PRODIGY_CONTACT_TYPES
from haddocking.protocols.protocol_prodigy import ProtPRODIGYPrediction

RSA_HELP = 'RSA (relative solvent accessibility) is the solvent accessible surface area (SASA) of a residue ' \
           'in the complex, computed by PRODIGY with FreeSASA and NACCESS atomic radii, divided by the maximum ' \
           'SASA of that amino acid type in an extended Ala-X-Ala tripeptide (NACCESS reference values). ' \
           '0 means fully buried and 1 means as exposed as in the reference; very exposed residues ' \
           '(e.g. termini) can have values above 1.'


class TableRow:
    def __init__(self, idx, values):
        self.idx = idx
        self.values = [str(v) for v in values]


class TableProvider(TreeProvider):
    def __init__(self, headers, rows, widths):
        TreeProvider.__init__(self)
        self._sortEnabled = True
        self.headers = headers
        self.widths = widths
        self.rows = [TableRow(i, row) for i, row in enumerate(rows)]

    def getColumns(self):
        return list(zip(self.headers, self.widths))

    def getObjects(self):
        colName = self.getSortingColumnName()
        if colName not in self.headers:
            return self.rows
        col = self.headers.index(colName)
        return sorted(self.rows, key=lambda r: self._sortKey(r.values[col]), reverse=not self.isSortingAscending())

    def getObjectInfo(self, obj):
        return {'key': obj.idx, 'text': obj.values[0], 'values': tuple(obj.values[1:])}

    @staticmethod
    def _sortKey(value):
        try:
            return 0, float(value), ''
        except ValueError:
            return 1, 0, value


class ProtPRODIGYPredictionViewer(ProtocolViewer):

    _label = 'PRODIGY viewer'
    _targets = [ProtPRODIGYPrediction]
    _environments = [DESKTOP_TKINTER]

    def _defineParams(self, form):
        form.addSection(label='PRODIGY visualization')

        form.addParam(
            'showInterface',
            params.LabelParam,
            label='Open PRODIGY interface visualization in PyMOL: '
        )

        form.addSection(label='Interface residues')
        form.addParam('showContacts', params.LabelParam,
                      label='Intermolecular contacts (ICs): ',
                      help='Residue pairs in contact between the two molecular groups, with the '
                           'physicochemical character of each residue and the contact class used by PRODIGY.\n'
                           'Use the Filter box to search (e.g. a chain, residue name or number) and click on a '
                           'column header to sort.')
        form.addParam('showInterfaceResidues', params.LabelParam,
                      label='Interface residues: ',
                      help='Residues involved in at least one intermolecular contact, with their character, '
                           'number of contacts and relative solvent accessibility (RSA) in the complex.\n\n'
                           + RSA_HELP)
        form.addParam('showSurface', params.LabelParam,
                      label='Surface residues (RSA): ',
                      help='Character and RSA of every solvent exposed residue of the complex (both molecular '
                           'groups), i.e. residues with RSA greater or equal than the accessibility threshold '
                           '(PRODIGY default 0.05). PRODIGY computes the NIS (non-interacting surface) '
                           'percentages of apolar, charged and polar residues over this set, which also includes '
                           'the interface residues (column "In interface").\n\n' + RSA_HELP)
        form.addParam('showComposition', params.LabelParam,
                      label='Contacts and NIS composition: ',
                      help='Number of contacts of each class and NIS percentages of apolar, charged and polar residues.')

    def _getVisualizeDict(self):
        return {
            'showInterface': self._showInterface,
            'showContacts': self._showContacts,
            'showInterfaceResidues': self._showInterfaceResidues,
            'showSurface': self._showSurface,
            'showComposition': self._showComposition
        }

    def _getResiduesData(self):
        resFile = self.protocol.getResiduesFile()
        if not os.path.exists(resFile):
            raise FileNotFoundError(
                'PRODIGY residue analysis file not found: {}. '
                'Please, run the protocol again to generate it.'.format(resFile))
        with open(resFile) as f:
            return json.load(f)

    def _showTable(self, title, headers, rows, widths, message, csvName):
        csvFile = os.path.abspath(self.protocol._getExtraPath('prodigy_{}.csv'.format(csvName)))
        with open(csvFile, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(headers)
            writer.writerows(rows)
        message += '\nTable saved in: {}'.format(csvFile)
        ListDialog(self.getTkRoot(), title, TableProvider(headers, rows, widths),
                   message=message, allowSelect=False, cancelButton=True)

    def _showContacts(self, paramName=None):
        data = self._getResiduesData()
        rows = []
        for c in data['contacts']:
            rows.append(tuple(c['res1']) + (PRODIGY_IC_CHARACTER.get(c['char1'], '-'),) +
                        tuple(c['res2']) + (PRODIGY_IC_CHARACTER.get(c['char2'], '-'),
                                            PRODIGY_CONTACT_TYPES.get(c['icType'], '-'),
                                            PRODIGY_CONTACT_TYPES.get(c['hydroType'], '-')))
        headers = ['Chain 1', 'Residue 1', 'Number 1', 'Character 1',
                   'Chain 2', 'Residue 2', 'Number 2', 'Character 2',
                   'Contact class', 'Hydrophobicity class']
        self._showTable('PRODIGY intermolecular contacts', headers, rows,
                        [70, 80, 80, 90, 70, 80, 80, 90, 150, 220],
                        '{} intermolecular contacts (ICs)'.format(len(rows)), 'contacts')

    def _showInterfaceResidues(self, paramName=None):
        data = self._getResiduesData()
        rows = [(r['chain'], r['resName'], r['resNum'], PRODIGY_IC_CHARACTER.get(r['icCharacter'], '-'),
                 r['nContacts'], round(r['rsa'], 3)) for r in data['residues'] if r['nContacts'] > 0]
        self._showTable('PRODIGY interface residues', ['Chain', 'Residue', 'Number', 'Character', 'Contacts', 'RSA'],
                        rows, [70, 80, 80, 90, 80, 80],
                        '{} residues in the interface. RSA: relative solvent accessibility in the complex '
                        '(see help).'.format(len(rows)), 'interface_residues')

    def _showSurface(self, paramName=None):
        data = self._getResiduesData()
        rows = [(r['chain'], r['resName'], r['resNum'], PRODIGY_IC_CHARACTER.get(r['nisCharacter'], '-'),
                 round(r['rsa'], 3), 'yes' if r['nContacts'] > 0 else 'no')
                for r in data['residues'] if r['surface']]
        self._showTable('PRODIGY surface residues', ['Chain', 'Residue', 'Number', 'Character', 'RSA', 'In interface'],
                        rows, [70, 80, 80, 90, 80, 100],
                        '{} surface residues with RSA >= {} used for the NIS percentages\n'
                        '(apolar {:.2f}%, charged {:.2f}%, polar {:.2f}%)'.format(
                            len(rows), data['accThreshold'], data['nis']['A'], data['nis']['C'], data['nis']['P']),
                        'surface_residues')

    def _showComposition(self, paramName=None):
        data = self._getResiduesData()
        rows = [('No. of {} contacts'.format(name), data['bins'][key])
                for key, name in PRODIGY_CONTACT_TYPES.items()]
        rows += [('Percentage of {} NIS residues'.format(PRODIGY_IC_CHARACTER[key]), round(data['nis'][key], 2))
                 for key in ['A', 'C', 'P']]
        self._showTable('PRODIGY contacts and NIS composition', ['Descriptor', 'Value'], rows, [330, 100],
                        '{} intermolecular contacts'.format(len(data['contacts'])), 'composition')

    def _showInterface(self, paramName=None):
        pmlFiles = glob.glob(
            self.protocol._getExtraPath('*.pml')
        )
        pmlFiles = [
            f for f in pmlFiles
            if os.path.basename(f) != 'prodigy_viewer.pml'
        ]

        if not pmlFiles:
            raise FileNotFoundError(
                'No PRODIGY PyMOL (.pml) file found in the extra directory.'
            )

        if len(pmlFiles) > 1:
            raise ValueError(
                'Multiple PyMOL scripts found: {}'.format(pmlFiles)
            )

        pmlFile = os.path.abspath(pmlFiles[0])

        # Original input structure (PDB or mmCIF)
        structureFile = os.path.abspath(
            self.protocol.inputAtomStruct.get().getFileName()
        )

        if not os.path.exists(structureFile):
            raise FileNotFoundError(
                'Input structure not found: {}'.format(structureFile)
            )

        # Create a wrapper PML to load the structure before
        # executing the PRODIGY visualization script.
        wrapperFile = os.path.abspath(self.protocol._getExtraPath('prodigy_viewer.pml'))

        with open(wrapperFile, 'w') as wrapper:
            wrapper.write('load "{}"\n'.format(structureFile))
            wrapper.write('@{}\n'.format(pmlFile))

        pymolV = PyMolViewer(project=self.getProject())

        return pymolV._visualize(
            wrapperFile,
            cwd=os.path.dirname(pmlFile)
        )