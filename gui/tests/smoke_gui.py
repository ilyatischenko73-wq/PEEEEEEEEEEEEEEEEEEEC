"""Offscreen Qt smoke test: widgets, node picking, VTK fields and subprocess lifecycle."""
import os
os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')
import sys
import tempfile
from pathlib import Path
import numpy as np
from gui.app import Window, W, C, pv
from gui.tests.test_model import MESH

app = W.QApplication([])
errors = []
W.QMessageBox.critical = lambda parent, title, text: errors.append(str(text))
with tempfile.TemporaryDirectory() as tmp:
    root = Path(tmp)
    mesh = root/'mesh.msh'
    mesh.write_text(MESH)
    window = Window()
    window.paths['mesh'].setText(str(mesh))
    window.paths['closed-mesh'].setText(str(mesh))
    window.paths['open-mesh'].setText(str(mesh))
    window.root.setText(str(root))
    window.output.setText(str(root/'output'))
    window.plot_pulse()
    assert window.chart_data[0] == 'line'
    window.task.setCurrentIndex(4)
    window.shunt.setChecked(True)
    window.role.setCurrentIndex(1)
    window.picked(np.array([1., 0, 0]))
    assert window.fields['strike-node'].text() == '1', errors
    window.role.setCurrentIndex(3)
    window.picked(np.array([0., 1, 0]))
    assert window.fields['shunt-a-node'].text() == '3'
    assert window.mesh_key == 'open-mesh'
    xyz = np.array([[0., 0, 0], [1., 0, 0], [1., 1, 0], [0., 1, 0]])
    for i in range(2):
        poly = pv.PolyData(xyz, [4, 0, 1, 2, 3])
        poly['phi_V'] = np.arange(4.)+i
        poly.save(root/f'run_surface_{i:06d}.vtk')
    window.set_result(root/'run_surface_000000.vtk')
    assert len(window.frames) == 2
    window.frame.setValue(1)
    np.testing.assert_array_equal(window.dataset['phi_V'], np.arange(4.)+1)
    window.fixed_scale.setChecked(True)
    window.cmin.setText('0')
    window.cmax.setText('5')
    window.render_result()
    # QProcess with Python verifies output and exit delivery, separate from actual solver integration.
    loop = C.QEventLoop()
    window.process.finished.connect(loop.quit)
    window.process.setProgram(sys.executable)
    window.process.setArguments(['-c', 'print("PEEC_PROCESS_OK")'])
    window.process.start()
    C.QTimer.singleShot(10000, loop.quit)
    loop.exec()
    assert window.process.state() == C.QProcess.ProcessState.NotRunning
    assert 'PEEC_PROCESS_OK' in window.log.toPlainText()
    # Failed launch must restore the Run button.
    window.process.setProgram(str(root/'does-not-exist'))
    window.process.setArguments([])
    window.run_button.setEnabled(False)
    window.process.start()
    window.process.waitForStarted(2000)
    app.processEvents()
    assert window.run_button.isEnabled()
    window.close()
assert not errors, errors
print('PASS: Qt construction, pulse plot, node assignment, VTK frames, scale, process output and failure recovery')
