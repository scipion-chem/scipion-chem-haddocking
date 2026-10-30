import os
import glob

from pyworkflow.viewer import ProtocolViewer, DESKTOP_TKINTER
from pyworkflow.protocol import params

from pwchem.viewers.viewers_data import PyMolViewer

from haddocking.protocols.protocol_prodigy import ProtPRODIGYPrediction


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

    def _getVisualizeDict(self):
        return {
            'showInterface': self._showInterface
        }

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