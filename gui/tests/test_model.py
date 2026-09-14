import tempfile
import unittest
from pathlib import Path
import numpy as np
from matplotlib.figure import Figure
from gui import model as m, charts

MESH = '''$MeshFormat
4.1 0 8
$EndMeshFormat
$Nodes
1 4 10 90
2 1 0 4
90
10
50
20
0 0 0
1 0 0
1 1 0
0 1 0
$EndNodes
$Elements
1 1 1 1
2 1 3 1
1 90 10 50 20
$EndElements
'''


class ModelTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.mesh = self.root/'mesh with spaces.msh'
        self.mesh.write_text(MESH)
        self.state = m.defaults()
        self.state['paths'] = {k: str(self.mesh) for k in m.PATHS}

    def test_node_order_is_file_order_not_sorted_tag(self):
        xyz, q, tags = m.read_mesh(self.mesh)
        np.testing.assert_array_equal(tags, [90, 10, 50, 20])
        np.testing.assert_array_equal(q, [[0, 1, 2, 3]])
        np.testing.assert_array_equal(xyz[1], [1, 0, 0])

    def test_mesh_rejects_duplicate_and_parametric(self):
        for content in (MESH.replace('90\n10\n50\n20', '90\n10\n50\n90'),
                        MESH.replace('2 1 0 4', '2 1 1 4')):
            self.mesh.write_text(content)
            with self.assertRaises(ValueError):
                m.read_mesh(self.mesh)

    def test_pulse_matches_fortran_and_delay(self):
        t = np.array([-1., 0., 1., 2.])
        np.testing.assert_allclose(m.pulse(t, 10, 2, 3, 1),
                                   [0, 0, 0, 10*(np.exp(-2)-np.exp(-3))])

    def test_config_task_modes_and_shunt_scope(self):
        for task in range(6):
            self.state['task'] = task
            self.state['shunt'] = True
            text = m.config_text(self.state, self.root, self.root/'results')
            self.assertIn('--task '+m.TASKS[task][1], text)
            self.assertEqual('--shunt-R ' in text, task in (4, 5))
            self.assertIn('"'+str(self.mesh)+'"', text)

    def test_invalid_numbers_steps_and_nodes(self):
        self.state['task'] = 4
        for key, value in [('dt', '0'), ('t-end', '1.23e-7'), ('parallel', '1.5'),
                           ('strike-node', '4'), ('strike-node', '1'), ('lightning-K', 'nan')]:
            old = self.state['values'][key]
            self.state['values'][key] = value
            with self.assertRaises(ValueError, msg=key):
                m.validate(self.state, self.root)
            self.state['values'][key] = old

    def test_inactive_parameters_are_ignored(self):
        self.state['values']['dt'] = 'not a number'
        m.validate(self.state, self.root)

    def test_rcs_allowed_orders_match_solver(self):
        self.state['task'] = 2
        for order in (2, 3, 4, 5, 6, 8, 10, 12, 16):
            self.state['values']['rcs-order'] = str(order)
            m.validate(self.state, self.root)
        self.state['values']['rcs-order'] = '1'
        with self.assertRaises(ValueError):
            m.validate(self.state, self.root)

    def test_project_roundtrip(self):
        path = self.root/'project.json'
        m.save_project(path, self.state)
        self.assertEqual(self.state, m.load_project(path))

    def test_config_roundtrip(self):
        for i in range(6):
            self.state['task'] = i
            p = self.root/'run.cfg'
            p.write_text(m.config_text(self.state, self.root, self.root/'output'), encoding='utf-8')
            other = m.import_config(p, self.root)
            self.assertEqual(other['task'], i)
            for key in m.active_keys(other):
                self.assertEqual(m.number(other['values'][key]), m.number(self.state['values'][key]))

    def test_coordinate_import_uses_solver_indices(self):
        p = self.root/'run.cfg'
        p.write_text('--task lightning --closed-mesh "'+str(self.mesh)+'" --strike 1 0 0 --return 0 1 0')
        state = m.import_config(p, self.root)
        self.assertEqual(state['values']['strike-node'], '1')
        self.assertEqual(state['values']['return-node'], '3')

    def test_windows_backslashes_in_import(self):
        p = self.root/'run.cfg'
        p.write_text('--task scattering --mesh "meshes\\body.msh"')
        state = m.import_config(p, self.root)
        self.assertIn('meshes\\body.msh', state['paths']['mesh']) if __import__('os').name != 'nt' else None

    def test_shunt_dat_series_does_not_mix_stages(self):
        for stem, t, current in [('stage2_shunt_000002', 2, 3), ('stage2_shunt_000001', 1, 4), ('stage1_shunt_000001', 1, 99)]:
            (self.root/(stem+'.dat')).write_text(f'time {t}\nshunt_current_A {current}\n')
        data = m.shunt_series(self.root/'stage2_shunt_000002.dat')
        np.testing.assert_array_equal(data, [[1, 4], [2, 3]])

    def test_csv_and_rcs_export_with_zero(self):
        p = self.root/'rcs.csv'
        p.write_text('phi_deg,sigma_m2\n0,1\n90,0\n180,2\n')
        data = m.read_table(p, ['phi_deg', 'sigma_m2'])
        for mode in ('db', 'linear', 'polar'):
            fig = Figure(figsize=(7, 4), layout='constrained')
            ax = charts.rcs_figure(fig, data, mode)
            if mode == 'db':
                self.assertTrue(np.isnan(ax.lines[0].get_ydata()[1]))
            for suffix in ('pdf', 'svg', 'png'):
                out = self.root/f'{mode}.{suffix}'
                charts.export(fig, out)
                self.assertGreater(out.stat().st_size, 1000)


if __name__ == '__main__':
    unittest.main()
