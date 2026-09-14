"""Run every GUI task against the real Fortran executable (CI, CPU)."""
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
from gui import model as m

root = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location('compare', root/'fortran/test/compare.py')
compare = importlib.util.module_from_spec(spec)
spec.loader.exec_module(compare)
exe = Path(sys.argv[1]).resolve()
with tempfile.TemporaryDirectory() as tmp:
    work = Path(tmp)
    mesh = work/'mesh.msh'
    edges = compare.fixture(mesh)
    aperture = work/'aperture.json'
    aperture.write_text(json.dumps({
        'node_mapping': [{'closed_node': v, 'open_node': v} for v in (1, 3, 5, 7)],
        'cover_edges': [i for i, ab in enumerate(edges) if 4 in ab]}))
    slots = work/'slots.dat'
    slots.write_text('0 8 0.5\n2 6 0.5\n')
    for task, excitation in [(i, 'lightning-current') for i in range(6)]+[(3, 'incident-pulse')]:
        state = m.defaults()
        state.update(task=task, shunt=True, excitation=excitation)
        state['paths'].update({k: str(mesh) for k in ('mesh', 'closed-mesh', 'open-mesh')})
        state['paths'].update({'aperture-map': str(aperture), 'slot-cells': str(slots)})
        state['values'].update({'t-end': '3e-7', 'return-node': '8', 'lightning-K': '1',
                                 'shunt-a-node': '1', 'shunt-b-node': '7', 'rcs-phi-step': '45'})
        output = work/f'run{task}_{excitation}'
        cfg = work/'run.cfg'
        cfg.write_text(m.config_text(state, root, output), encoding='utf-8')
        process = subprocess.run([str(exe), '--config', str(cfg)], cwd=root,
                                 capture_output=True, text=True, timeout=120)
        assert process.returncode == 0, process.stdout+process.stderr
        if task:
            assert list(output.glob('*.vtk')), output
        if task == 2:
            assert len(m.read_table(output/'rcs.csv', ['phi_deg', 'sigma_m2'])) == 9
        if task in (4, 5):
            csv = next(output.glob('*stage2.csv')) if task == 4 else next(output.glob('*.csv'))
            assert len(m.shunt_series(csv)) == 4
        print('PASS: GUI configuration → Fortran:', m.TASKS[task][0], excitation)
