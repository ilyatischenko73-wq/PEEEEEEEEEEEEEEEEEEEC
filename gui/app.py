"""PEEC Studio: Qt desktop front end for the standalone Fortran solver."""
from __future__ import annotations
import os
os.environ.setdefault('QT_API', 'pyside6')
import sys
import functools
import locale
from pathlib import Path
from datetime import datetime
import numpy as np
from PySide6 import QtCore as C, QtWidgets as W
import pyvista as pv
from pyvistaqt import QtInteractor
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT
from matplotlib.figure import Figure
from matplotlib.font_manager import findfont
from . import model as m
from . import charts

ROOT = Path(__file__).resolve().parents[1]


def safe(method):
    @functools.wraps(method)
    def wrapper(self, *args, **kwargs):
        try:
            return method(self, *args, **kwargs)
        except Exception as exc:
            W.QMessageBox.critical(self, 'PEEC — ошибка', str(exc))
    return wrapper


def button(text, callback):
    b = W.QPushButton(text)
    b.clicked.connect(lambda checked=False: callback())
    return b


class Window(W.QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle('PEEC Studio — расчёт и визуализация')
        self.resize(1450, 960)
        self.fields, self.paths = {}, {}
        self.mesh_data = None
        self.mesh_key = None
        self.mesh_path = None
        self.frames, self.dataset, self.chart_data = [], None, None
        self.current_run = None
        self.process = C.QProcess(self)
        self.process.setProcessChannelMode(C.QProcess.ProcessChannelMode.MergedChannels)
        self.process.readyReadStandardOutput.connect(self.process_output)
        self.process.finished.connect(self.finished)
        self.process.errorOccurred.connect(self.process_error)
        self.timer = C.QTimer(self)
        self.timer.setInterval(150)
        self.timer.timeout.connect(self.next_frame)
        self.tabs = W.QTabWidget()
        self.setCentralWidget(self.tabs)
        self.build_setup()
        self.build_mesh()
        self.build_results()
        self.build_charts()
        self.apply(m.defaults())
        plate = ROOT/'meshes/plate.msh'
        if plate.exists():
            self.paths['mesh'].setText(str(plate))
        self.statusBar().showMessage('Готово. Настройте задачу и выберите сетку.')
        self.setStyleSheet('QGroupBox {font-weight:600; margin-top:12px; padding-top:12px;} '
                           'QPushButton {padding:6px 10px;} QLineEdit,QComboBox {padding:4px;} '
                           'QTabBar::tab {padding:10px 16px;}')

    def panel(self, title):
        page = W.QWidget()
        layout = W.QVBoxLayout(page)
        self.tabs.addTab(page, title)
        return page, layout

    def file_row(self, edit, kind='file', filt='Все файлы (*)'):
        row = W.QWidget()
        layout = W.QHBoxLayout(row)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(edit)
        def choose():
            if kind == 'dir':
                p = W.QFileDialog.getExistingDirectory(self, 'Выберите каталог', edit.text())
            else:
                p, _ = W.QFileDialog.getOpenFileName(self, 'Выберите файл', edit.text(), filt)
            if p:
                edit.setText(p)
        layout.addWidget(button('Обзор…', choose))
        return row

    def build_setup(self):
        page, layout = self.panel('1 · Расчёт')
        split = W.QSplitter()
        layout.addWidget(split, 1)
        settings = W.QWidget()
        form = W.QFormLayout(settings)
        self.root = W.QLineEdit(str(ROOT))
        exe = ROOT/'build/windows/PeecSolverFortran.exe'
        if os.name != 'nt':
            exe = ROOT/'build/fortran/PeecSolverFortran'
        self.exe = W.QLineEdit(str(exe))
        self.output = W.QLineEdit(str(ROOT/'results/GUI'))
        form.addRow('Корень проекта', self.file_row(self.root, 'dir'))
        form.addRow('Решатель Fortran', self.file_row(self.exe))
        form.addRow('Каталог новых расчётов', self.file_row(self.output, 'dir'))
        self.task = W.QComboBox()
        self.task.addItems([t[0] for t in m.TASKS])
        self.task.currentIndexChanged.connect(self.refresh_fields)
        form.addRow('Задача', self.task)
        for key, label in m.PATHS.items():
            edit = W.QLineEdit()
            self.paths[key] = edit
            form.addRow(label, self.file_row(edit))
        form.addRow(button('Показать сетку / выбрать узлы', self.show_setup_mesh))
        self.excitation = W.QComboBox()
        self.excitation.addItem('Ток молнии', 'lightning-current')
        self.excitation.addItem('Падающий импульс поля', 'incident-pulse')
        self.excitation.currentIndexChanged.connect(self.refresh_fields)
        form.addRow('Возбуждение переходного процесса', self.excitation)
        self.polarization = W.QComboBox()
        self.polarization.addItem('Горизонтальная', 'horizontal')
        self.polarization.addItem('Вертикальная', 'vertical')
        form.addRow('Поляризация', self.polarization)
        self.shunt = W.QCheckBox('Включить внутренний шунт')
        self.shunt.toggled.connect(self.refresh_fields)
        self.dual = W.QCheckBox('ЭПР: интегрирование двойственных областей')
        self.vtk = W.QCheckBox('Сохранять поля VTK')
        self.matrices = W.QCheckBox('Сохранять матрицы')
        for w in (self.shunt, self.dual, self.vtk, self.matrices):
            form.addRow(w)
        tools = W.QHBoxLayout()
        tools.addWidget(button('Сохранить проект…', self.save_project))
        tools.addWidget(button('Открыть проект…', self.open_project))
        tools.addWidget(button('Импорт .cfg…', self.import_config))
        tools.addWidget(button('Экспорт .cfg…', self.export_config))
        form.addRow(tools)
        note = W.QLabel('Каждый запуск создаёт отдельную папку с конфигурацией, журналом и результатами.\n'
                        'Нумерация узлов — с нуля. Для closed-режима решатель не использует шунт.')
        note.setWordWrap(True)
        form.addRow(note)
        scroll = W.QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setWidget(settings)
        split.addWidget(scroll)
        paramtabs = W.QTabWidget()
        split.addWidget(paramtabs)
        groups = [
            ('Волна и ЭПР', ['frequency', 'field-amplitude', 'theta', 'phi', 'phase', 'ka',
                            'characteristic-length', 'rcs-observation-theta', 'rcs-phi-start',
                            'rcs-phi-end', 'rcs-phi-step', 'rcs-order']),
            ('Молния и шунт', ['lightning-K', 'lightning-alpha', 'lightning-beta', 'lightning-delay',
                              'strike-node', 'return-node', 'shunt-a-node', 'shunt-b-node', 'shunt-R',
                              'pulse-K', 'pulse-alpha', 'pulse-beta']),
            ('Время и щель', ['dt', 't-end', 'time-order', 'parallel', 'vtk-every', 'slot-width',
                             'slot-wall-span', 'slot-wall-thickness', 'slot-epsilon-r', 'slot-mu-r'])]
        for title, keys in groups:
            box = W.QWidget()
            f = W.QFormLayout(box)
            for key in keys:
                edit = W.QLineEdit(m.FIELDS[key][1])
                edit.setToolTip('--'+key)
                self.fields[key] = edit
                f.addRow(m.FIELDS[key][0], edit)
            area = W.QScrollArea()
            area.setWidgetResizable(True)
            area.setWidget(box)
            paramtabs.addTab(area, title)
        actions = W.QHBoxLayout()
        self.run_button = button('▶ Запустить расчёт', self.run)
        self.stop_button = button('■ Остановить', self.stop)
        self.stop_button.setEnabled(False)
        actions.addWidget(self.run_button)
        actions.addWidget(self.stop_button)
        actions.addWidget(button('Предпросмотр импульса', self.plot_pulse))
        self.run_status = W.QLabel('Расчёт не запущен')
        self.run_status.setWordWrap(True)
        actions.addWidget(self.run_status, 1)
        layout.addLayout(actions)
        self.log = W.QPlainTextEdit()
        self.log.setReadOnly(True)
        self.log.setFixedHeight(140)
        layout.addWidget(self.log)

    def state(self):
        return dict(values={k: v.text().strip() for k, v in self.fields.items()},
                    paths={k: v.text().strip() for k, v in self.paths.items()},
                    task=self.task.currentIndex(), polarization=self.polarization.currentData(),
                    excitation=self.excitation.currentData(), shunt=self.shunt.isChecked(),
                    dual=self.dual.isChecked(), vtk=self.vtk.isChecked(), matrices=self.matrices.isChecked())

    def apply(self, state):
        for key, val in state['values'].items():
            self.fields[key].setText(val)
        for key, val in state['paths'].items():
            self.paths[key].setText(val)
        self.task.setCurrentIndex(state['task'])
        self.polarization.setCurrentIndex(self.polarization.findData(state['polarization']))
        self.excitation.setCurrentIndex(self.excitation.findData(state['excitation']))
        for key in ('shunt', 'dual', 'vtk', 'matrices'):
            getattr(self, key).setChecked(state[key])
        self.refresh_fields()

    def refresh_fields(self, *args):
        if not self.fields or not hasattr(self, 'shunt'):
            return
        state = self.state()
        active = m.active_keys(state)
        for key, widget in self.fields.items():
            widget.setEnabled(key in active)
        task, case = m.TASKS[state['task']][1:]
        self.shunt.setEnabled(task == 'lightning' and case != 'closed')
        self.excitation.setEnabled(task == 'lightning')
        self.dual.setEnabled(task == 'rcs')
        self.polarization.setEnabled(task in ('scattering', 'rcs') or
                                     (task == 'lightning' and state['excitation'] == 'incident-pulse'))
        for key, edit in self.paths.items():
            edit.parentWidget().setEnabled(key in m.required_paths(state))

    @safe
    def save_project(self):
        path, _ = W.QFileDialog.getSaveFileName(self, 'Сохранить параметры', str(ROOT/'project.peec.json'), 'PEEC (*.json)')
        if path:
            state = self.state()
            # Absolute mesh paths survive moving the project JSON.
            state['paths'] = {k: str(m.resolve(self.root.text(), v)) if v else '' for k, v in state['paths'].items()}
            m.save_project(path, state)

    @safe
    def open_project(self):
        path, _ = W.QFileDialog.getOpenFileName(self, 'Открыть параметры', str(ROOT), 'PEEC (*.json)')
        if path:
            self.apply(m.load_project(path))

    @safe
    def import_config(self):
        path, _ = W.QFileDialog.getOpenFileName(self, 'Импорт конфигурации', self.root.text(), 'Конфигурация (*.cfg)')
        if path:
            self.apply(m.import_config(path, self.root.text()))
            self.statusBar().showMessage('Параметры импортированы. Результаты будут записаны в новую папку запуска.')

    @safe
    def export_config(self):
        output = m.resolve(self.root.text(), self.output.text())
        text = m.config_text(self.state(), self.root.text(), output)
        path, _ = W.QFileDialog.getSaveFileName(self, 'Экспорт конфигурации', str(ROOT/'gui.cfg'), 'Конфигурация (*.cfg)')
        if path:
            Path(path).write_text(text, encoding='utf-8')

    @safe
    def run(self):
        if self.process.state() != C.QProcess.ProcessState.NotRunning:
            return
        root = Path(self.root.text()).resolve()
        exe = m.resolve(root, self.exe.text())
        if not exe.is_file():
            raise ValueError('Не найден исполняемый файл решателя. Соберите Fortran и выберите EXE.')
        output = m.resolve(root, self.output.text()) / datetime.now().strftime('%Y%m%d_%H%M%S_%f')
        if os.name == 'nt' and any(c in str(output) for c in '"%!&|<>^'):
            raise ValueError('Выходной путь содержит неподдерживаемые символы Windows.')
        text = m.config_text(self.state(), root, output)
        output.mkdir(parents=True, exist_ok=False)
        cfg = output/'run.cfg'
        cfg.write_text(text, encoding='utf-8')
        snapshot = self.state()
        snapshot['paths'] = {k: str(m.resolve(root, v)) if v else '' for k, v in snapshot['paths'].items()}
        m.save_project(output/'project.json', snapshot)
        self.current_run = output
        self.log.clear()
        self.process.setWorkingDirectory(str(root))
        self.process.setProgram(str(exe))
        self.process.setArguments(['--config', str(cfg)])
        self.run_button.setEnabled(False)
        self.stop_button.setEnabled(True)
        self.run_status.setText('Расчёт выполняется…')
        self.process.start()

    def process_output(self):
        data = bytes(self.process.readAllStandardOutput())
        text = data.decode(locale.getpreferredencoding(False), errors='replace')
        self.log.moveCursor(self.log.textCursor().MoveOperation.End)
        self.log.insertPlainText(text)
        if self.current_run:
            with (self.current_run/'solver.log').open('ab') as f:
                f.write(data)

    def finished(self, code, status):
        self.process_output()
        self.run_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        ok = code == 0 and status == C.QProcess.ExitStatus.NormalExit
        self.run_status.setText(('Готово: ' if ok else 'Прервано / ошибка: ')+str(self.current_run))
        self.statusBar().showMessage(self.run_status.text())

    def process_error(self, error):
        self.log.appendPlainText(self.process.errorString())
        if error == C.QProcess.ProcessError.FailedToStart:
            self.run_button.setEnabled(True)
            self.stop_button.setEnabled(False)
            self.run_status.setText('Не удалось запустить решатель')

    def stop(self):
        # QProcess is the directly launched solver, no shell or child process tree.
        self.process.kill()

    def build_mesh(self):
        page, layout = self.panel('2 · Сетка и узлы')
        row = W.QHBoxLayout()
        self.mesh_choice = W.QComboBox()
        for key in ('mesh', 'closed-mesh', 'open-mesh'):
            self.mesh_choice.addItem(m.PATHS[key], key)
        row.addWidget(self.mesh_choice)
        row.addWidget(button('Загрузить сетку', self.load_mesh))
        self.role = W.QComboBox()
        for label, key in [('Осмотр', ''), ('Удар', 'strike-node'), ('Возврат', 'return-node'),
                           ('Шунт A', 'shunt-a-node'), ('Шунт B', 'shunt-b-node')]:
            self.role.addItem(label, key)
        row.addWidget(W.QLabel('Назначить выбранный узел:'))
        row.addWidget(self.role)
        self.role.currentIndexChanged.connect(self.role_changed)
        row.addWidget(button('Изометрия', lambda: self.mesh_view.view_isometric()))
        layout.addLayout(row)
        self.mesh_info = W.QLabel('Укажите файл на вкладке «Расчёт». Выбор точки — правой кнопкой мыши.')
        self.mesh_info.setWordWrap(True)
        layout.addWidget(self.mesh_info)
        self.mesh_view = QtInteractor(page)
        self.mesh_view.set_background('white')
        layout.addWidget(self.mesh_view.interactor)

    @safe
    def show_setup_mesh(self):
        key = m.required_paths(self.state())[0]
        self.mesh_choice.setCurrentIndex(self.mesh_choice.findData(key))
        self.tabs.setCurrentIndex(1)
        self.load_mesh()

    def role_mesh(self):
        role = self.role.currentData()
        case = m.TASKS[self.task.currentIndex()][2]
        if role in ('shunt-a-node', 'shunt-b-node'):
            return 'open-mesh'
        return 'open-mesh' if case == 'slot-cells' else 'closed-mesh'

    @safe
    def role_changed(self, *args):
        if self.role.currentData():
            key = self.role_mesh()
            self.mesh_choice.setCurrentIndex(self.mesh_choice.findData(key))
            self.load_mesh()

    @safe
    def load_mesh(self):
        key = self.mesh_choice.currentData()
        path = m.resolve(self.root.text(), self.paths[key].text())
        data = m.read_mesh(path)
        xyz, quads, tags = data
        poly = pv.PolyData(xyz, np.column_stack((np.full(len(quads), 4), quads)).ravel())
        self.mesh_view.disable_picking()
        self.mesh_view.clear()
        self.mesh_view.add_mesh(poly, color='#cbd5e1', show_edges=True)
        self.mesh_view.add_axes()
        self.mesh_data, self.mesh_key, self.mesh_path = data, key, path
        self.mesh_view.enable_point_picking(callback=self.picked, show_message=False,
                                            show_point=False, picker='point', left_clicking=False)
        self.mesh_view.view_isometric()
        self.mesh_info.setText(f'{path.name}: {len(xyz)} узлов, {len(quads)} четырёхугольников. '
                               'Правый щелчок выбирает узел; его индекс совпадает с Fortran.')

    @safe
    def picked(self, point):
        if self.mesh_data is None:
            return
        xyz, _, tags = self.mesh_data
        idx = int(np.argmin(np.linalg.norm(xyz-np.asarray(point), axis=1)))
        role = self.role.currentData()
        if role:
            if role not in m.active_keys(self.state()):
                raise ValueError('Выбранный вывод не используется текущей задачей. Проверьте задачу и включение шунта.')
            expected = self.role_mesh()
            if self.mesh_key != expected or self.mesh_path != m.resolve(self.root.text(), self.paths[expected].text()):
                raise ValueError('Сетка изменилась. Загрузите сетку заново для выбранного вывода.')
            self.fields[role].setText(str(idx))
        self.mesh_view.add_point_labels(xyz[[idx]], [f'{idx}'], name='picked_node',
                                        point_size=14, point_color='#ef4444', always_visible=True)
        self.mesh_info.setText(f'Индекс Fortran: {idx}; тег Gmsh: {tags[idx]}; '
                               f'x={xyz[idx,0]:.9g}, y={xyz[idx,1]:.9g}, z={xyz[idx,2]:.9g} м. '
                               + (f'Назначен: {self.role.currentText()}.' if role else ''))

    def build_results(self):
        page, layout = self.panel('3 · Поля результатов')
        row = W.QHBoxLayout()
        row.addWidget(button('Открыть VTK…', self.open_result))
        self.field = W.QComboBox()
        self.field.currentIndexChanged.connect(self.render_result)
        row.addWidget(self.field, 1)
        self.component = W.QComboBox()
        self.component.addItems(['Модуль', 'X', 'Y', 'Z'])
        self.component.currentIndexChanged.connect(self.render_result)
        row.addWidget(self.component)
        self.colormap = W.QComboBox()
        self.colormap.addItems(['viridis', 'coolwarm', 'plasma', 'cividis'])
        self.colormap.currentIndexChanged.connect(self.render_result)
        row.addWidget(self.colormap)
        self.edges = W.QCheckBox('Рёбра сетки')
        self.edges.toggled.connect(self.render_result)
        row.addWidget(self.edges)
        self.arrows = W.QCheckBox('Стрелки векторов')
        self.arrows.toggled.connect(self.render_result)
        row.addWidget(self.arrows)
        row.addWidget(button('Снимок PNG…', self.screenshot))
        layout.addLayout(row)
        scale = W.QHBoxLayout()
        self.fixed_scale = W.QCheckBox('Фиксированная цветовая шкала')
        self.cmin, self.cmax = W.QLineEdit('0'), W.QLineEdit('1')
        for edit in (self.cmin, self.cmax):
            edit.setMaximumWidth(120)
        scale.addWidget(self.fixed_scale)
        scale.addWidget(W.QLabel('Мин.'))
        scale.addWidget(self.cmin)
        scale.addWidget(W.QLabel('Макс.'))
        scale.addWidget(self.cmax)
        scale.addWidget(button('Применить', self.render_result))
        scale.addWidget(button('Изометрия', lambda: self.result_view.view_isometric()))
        scale.addStretch()
        layout.addLayout(scale)
        self.result_info = W.QLabel('Выберите surface.vtk для заряда/потенциала или edges.vtk для токов.')
        self.result_info.setWordWrap(True)
        layout.addWidget(self.result_info)
        self.result_view = QtInteractor(page)
        self.result_view.set_background('white')
        layout.addWidget(self.result_view.interactor)
        animation = W.QHBoxLayout()
        self.play = button('▶ / ❚❚', self.toggle_animation)
        animation.addWidget(self.play)
        self.frame = W.QSlider(C.Qt.Orientation.Horizontal)
        self.frame.valueChanged.connect(self.load_frame)
        animation.addWidget(self.frame, 1)
        self.frame_info = W.QLabel('Нет кадров')
        animation.addWidget(self.frame_info)
        layout.addLayout(animation)

    @safe
    def open_result(self):
        path, _ = W.QFileDialog.getOpenFileName(self, 'Результат расчёта', str(self.current_run or ROOT), 'VTK (*.vtk)')
        if path:
            self.set_result(Path(path))

    def set_result(self, path):
        self.timer.stop()
        self.frames = m.result_series(path)
        self.frame.blockSignals(True)
        self.frame.setRange(0, len(self.frames)-1)
        self.frame.setValue(self.frames.index(path))
        self.frame.blockSignals(False)
        self.load_frame(self.frame.value(), reset=True)

    @safe
    def load_frame(self, index, reset=False):
        self.timer.stop() if not self.frames else None
        if not self.frames:
            return
        try:
            data = pv.read(self.frames[index])
            previous = self.field.currentData()
            self.dataset = data
            self.field.blockSignals(True)
            self.field.clear()
            for assoc, attrs in [('point', data.point_data), ('cell', data.cell_data)]:
                for key in attrs.keys():
                    array = np.asarray(attrs[key])
                    if np.issubdtype(array.dtype, np.number) and (array.ndim == 1 or array.shape[1] == 3):
                        self.field.addItem(self.field_label(key)+' · '+('узлы' if assoc == 'point' else 'рёбра/ячейки'), (assoc, key))
            found = self.field.findData(previous)
            if found >= 0:
                self.field.setCurrentIndex(found)
            self.field.blockSignals(False)
            self.frame_info.setText(f'{index+1} / {len(self.frames)}')
            self.render_result()
            if reset:
                self.result_view.view_isometric()
        except Exception:
            self.field.blockSignals(False)
            self.timer.stop()
            raise

    @staticmethod
    def field_label(key):
        for prefix, name in [('phi', 'Потенциал'), ('charge', 'Заряд'), ('sigma', 'Плотность заряда'),
                             ('current', 'Ток'), ('J', 'Поверхностный ток')]:
            if key.startswith(prefix+'_'):
                part = ' · Re' if '_re_' in key else ' · Im' if '_im_' in key else ' · |·|' if '_abs_' in key else ''
                unit = {'phi': 'В', 'charge': 'Кл', 'sigma': 'Кл/м²', 'current': 'А', 'J': 'А/м'}[prefix]
                return name+part+', '+unit+' ['+key+']'
        return 'Время, с' if key == 'time' else key

    @safe
    def render_result(self, *args):
        if self.dataset is None or self.field.currentData() is None:
            return
        assoc, key = self.field.currentData()
        attrs = self.dataset.point_data if assoc == 'point' else self.dataset.cell_data
        values = np.asarray(attrs[key])
        vector = values.ndim == 2
        self.component.setEnabled(vector)
        self.arrows.setEnabled(vector)
        if vector:
            c = self.component.currentIndex()
            scalars = np.linalg.norm(values, axis=1) if c == 0 else values[:, c-1]
        else:
            scalars = values
        if not np.isfinite(scalars).all():
            self.timer.stop()
            raise ValueError('Поле содержит NaN или бесконечность.')
        clim = None
        if self.fixed_scale.isChecked():
            clim = (m.number(self.cmin.text()), m.number(self.cmax.text()))
            if clim[0] >= clim[1]:
                self.timer.stop()
                raise ValueError('Минимум цветовой шкалы должен быть меньше максимума.')
        camera = self.result_view.camera_position
        self.result_view.clear()
        # Companion surface provides geometry context for edge currents.
        path = self.frames[self.frame.value()]
        companion = path.with_name(path.name.replace('_edges', '_surface'))
        if companion != path and companion.exists():
            self.result_view.add_mesh(pv.read(companion), color='lightgray', opacity=.15, pickable=False)
        self.result_view.add_mesh(self.dataset, scalars=scalars, preference=assoc,
                                  cmap=self.colormap.currentText(), clim=clim,
                                  show_edges=self.edges.isChecked(), line_width=3,
                                  scalar_bar_args={'title': self.field_label(key), 'fmt': '%.2e'})
        # VTK built-in fonts may omit Cyrillic; use Matplotlib's bundled font.
        bar = self.result_view.scalar_bar
        for prop in (bar.GetTitleTextProperty(), bar.GetLabelTextProperty()):
            prop.SetFontFamily(4)  # VTK_FONT_FILE
            prop.SetFontFile(findfont('DejaVu Sans'))
        if vector and self.arrows.isChecked():
            centers = self.dataset.points if assoc == 'point' else self.dataset.cell_centers().points
            magnitude = np.linalg.norm(values, axis=1)
            maximum = float(np.max(magnitude))
            if maximum > 0:
                mask = magnitude >= .1*maximum
                cloud = pv.PolyData(centers[mask])
                cloud['direction'] = values[mask]/maximum
                glyph = cloud.glyph(orient='direction', scale='direction', factor=.08*self.dataset.length)
                self.result_view.add_mesh(glyph, color='#be123c', show_scalar_bar=False)
        self.result_view.add_axes()
        self.result_view.camera_position = camera
        t = ''
        if 'time' in self.dataset.point_data:
            t = f'; t={self.dataset.point_data["time"][0]:.6g} с'
        self.result_info.setText(f'{path.name}{t} · минимум {np.min(scalars):.6g}, максимум {np.max(scalars):.6g}')

    def toggle_animation(self):
        if self.timer.isActive():
            self.timer.stop()
        elif len(self.frames) > 1:
            self.timer.start()

    def next_frame(self):
        if self.frames:
            self.frame.setValue((self.frame.value()+1) % len(self.frames))

    @safe
    def screenshot(self):
        path, _ = W.QFileDialog.getSaveFileName(self, 'Сохранить вид', str(ROOT/'field.png'), 'PNG (*.png)')
        if path:
            self.result_view.screenshot(path, scale=2)

    def build_charts(self):
        page, layout = self.panel('4 · Графики')
        row = W.QHBoxLayout()
        row.addWidget(button('ЭПР из CSV…', self.open_rcs))
        row.addWidget(button('Заданный импульс', self.plot_pulse))
        row.addWidget(button('Ток шунта из CSV / DAT…', self.open_shunt))
        row.addWidget(button('Ток источника из CSV…', self.open_source))
        self.rcs_mode = W.QComboBox()
        self.rcs_mode.addItem('ЭПР: дБм²', 'db')
        self.rcs_mode.addItem('ЭПР: м²', 'linear')
        self.rcs_mode.addItem('ЭПР: полярная, м²', 'polar')
        self.rcs_mode.currentIndexChanged.connect(self.redraw_chart)
        row.addWidget(self.rcs_mode)
        layout.addLayout(row)
        row = W.QHBoxLayout()
        self.chart_title = W.QLineEdit()
        self.chart_title.setPlaceholderText('Свой заголовок графика (необязательно)')
        row.addWidget(self.chart_title, 1)
        row.addWidget(button('Применить заголовок', self.redraw_chart))
        row.addWidget(button('Экспорт PDF / SVG / PNG…', self.export_chart))
        layout.addLayout(row)
        self.figure = Figure(figsize=(8, 5), layout='constrained')
        self.canvas = FigureCanvasQTAgg(self.figure)
        layout.addWidget(NavigationToolbar2QT(self.canvas, self))
        layout.addWidget(self.canvas, 1)
        self.chart_info = W.QLabel('Экспорт: векторные PDF/SVG или PNG 600 dpi; шрифт с поддержкой кириллицы.')
        self.chart_info.setWordWrap(True)
        layout.addWidget(self.chart_info)

    @safe
    def plot_pulse(self):
        prefix = 'pulse' if self.excitation.currentData() == 'incident-pulse' else 'lightning'
        k, a, b = (m.number(self.fields[prefix+'-'+s].text()) for s in ('K', 'alpha', 'beta'))
        delay = m.number(self.fields['lightning-delay'].text()) if prefix == 'lightning' else 0.
        end = m.number(self.fields['t-end'].text())
        if min(a, b, end) <= 0 or delay < 0:
            raise ValueError('α, β и время должны быть положительными, задержка — неотрицательной.')
        t = np.linspace(0, end, 2001)
        data = np.column_stack((t*1e6, m.pulse(t, k, a, b, delay)))
        electric = prefix == 'pulse'
        self.chart_data = ('line', data, 'Заданный импульс поля' if electric else 'Заданный ток молнии',
                           'Время t, мкс', 'Напряжённость E, В/м' if electric else 'Ток I, А')
        self.chart_info.setText('Предпросмотр источника: K [exp(−ατ) − exp(−βτ)], τ=t−задержка, ноль при τ<0. '
                               'K — коэффициент, не пиковое значение. Для поля показан импульс в начале координат.')
        self.chart_title.clear()
        self.redraw_chart()
        self.tabs.setCurrentIndex(3)

    @safe
    def open_rcs(self):
        path, _ = W.QFileDialog.getOpenFileName(self, 'ЭПР', str(self.current_run or ROOT), 'CSV (*.csv)')
        if path:
            data = m.read_table(path, ['phi_deg', 'sigma_m2'])
            self.chart_data = ('rcs', data)
            self.chart_title.clear()
            self.chart_info.setText(path+' · нулевые значения в дБ показаны разрывами; сглаживание не применяется.')
            self.redraw_chart()

    @safe
    def open_shunt(self):
        path, _ = W.QFileDialog.getOpenFileName(self, 'Ток шунта', str(self.current_run or ROOT), 'Данные (*.csv *.dat)')
        if path:
            data = m.shunt_series(path)
            data[:, 0] *= 1e6
            self.chart_data = ('line', data, 'Ток во внутреннем шунте', 'Время t, мкс', 'Ток I, А')
            self.chart_title.clear()
            self.chart_info.setText(path+' · ток из результатов расчёта; направление A → B. '
                                   'При выборе DAT читаются кадры того же префикса. CSV содержит все шаги.')
            self.redraw_chart()

    @safe
    def open_source(self):
        path, _ = W.QFileDialog.getOpenFileName(self, 'Ток источника', str(self.current_run or ROOT), 'CSV (*.csv)')
        if path:
            data = m.read_table(path, ['time_s', 'source_A'])
            if not np.isfinite(data).all():
                raise ValueError('Некорректные значения источника.')
            data[:, 0] *= 1e6
            self.chart_data = ('line', data, 'Ток источника из расчёта', 'Время t, мкс', 'Ток I, А')
            self.chart_title.clear()
            self.chart_info.setText(path+' · для stage2 source_A — половина суммы модулей токов связи, не ток молнии.')
            self.redraw_chart()

    @safe
    def redraw_chart(self, *args):
        if self.chart_data is None:
            return
        kind, data, *labels = self.chart_data
        custom = self.chart_title.text().strip()
        if kind == 'rcs':
            charts.rcs_figure(self.figure, data, self.rcs_mode.currentData(),
                              custom or 'Диаграмма эффективной площади рассеяния')
        else:
            title, xlabel, ylabel = labels
            charts.draw(self.figure, data[:, 0], data[:, 1], custom or title, xlabel, ylabel)
        self.canvas.draw_idle()

    @safe
    def export_chart(self):
        if self.chart_data is None:
            raise ValueError('Сначала постройте график.')
        path, _ = W.QFileDialog.getSaveFileName(self, 'Экспорт графика', str(ROOT/'figure.pdf'),
                                                'PDF (*.pdf);;SVG (*.svg);;PNG (*.png)')
        if path:
            charts.export(self.figure, path)

    def closeEvent(self, event):
        if self.process.state() != C.QProcess.ProcessState.NotRunning:
            answer = W.QMessageBox.question(self, 'Расчёт выполняется', 'Остановить расчёт и закрыть окно?')
            if answer != W.QMessageBox.StandardButton.Yes:
                event.ignore()
                return
            self.process.kill()
            self.process.waitForFinished(3000)
        self.timer.stop()
        self.mesh_view.close()
        self.result_view.close()
        event.accept()


def main():
    app = W.QApplication(sys.argv)
    app.setStyle('Fusion')
    window = Window()
    window.show()
    sys.exit(app.exec())


if __name__ == '__main__':
    main()
