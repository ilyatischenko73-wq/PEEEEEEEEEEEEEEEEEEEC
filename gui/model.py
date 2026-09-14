"""PEEC data model; no Qt or VTK imports, solver-compatible node ordering."""
from pathlib import Path
import csv
import json
import math
import re
import shlex
import numpy as np

# key: (Russian label, default). Values are strings to support scientific notation.
FIELDS = {
    'frequency': ('Частота, Гц', '1e6'), 'field-amplitude': ('Амплитуда E, В/м', '1'),
    'theta': ('Угол θ, град', '90'), 'phi': ('Угол φ, град', '0'),
    'phase': ('Начальная фаза, рад', '0'), 'ka': ('Электрический размер ka', '0.02'),
    'characteristic-length': ('Характерный размер a, м', '1'),
    'rcs-observation-theta': ('Угол наблюдения θ, град', '90'),
    'rcs-phi-start': ('Начало φ, град', '0'), 'rcs-phi-end': ('Конец φ, град', '360'),
    'rcs-phi-step': ('Шаг φ, град', '1'), 'rcs-order': ('Порядок квадратуры ЭПР', '4'),
    'lightning-K': ('Коэффициент K, А', '102900'),
    'lightning-alpha': ('Коэффициент α, 1/с', '1500'),
    'lightning-beta': ('Коэффициент β, 1/с', '1e6'),
    'lightning-delay': ('Задержка, с', '0'),
    'pulse-K': ('Амплитуда импульса E, В/м', '1'),
    'pulse-alpha': ('α импульса E, 1/с', '1500'), 'pulse-beta': ('β импульса E, 1/с', '1e6'),
    'dt': ('Шаг времени, с', '1e-7'), 't-end': ('Время окончания, с', '1e-5'),
    'time-order': ('Схема: 1 — Эйлер, 2 — трапеции', '2'),
    'strike-node': ('Узел удара (с нуля)', '0'), 'return-node': ('Узел возврата (с нуля)', '1'),
    'shunt-a-node': ('Узел шунта A (с нуля)', '0'), 'shunt-b-node': ('Узел шунта B (с нуля)', '1'),
    'shunt-R': ('Сопротивление шунта, Ом', '1000'),
    'slot-width': ('Ширина щели, м', '0.1'), 'slot-wall-span': ('Ширина стенки, м', '2'),
    'slot-wall-thickness': ('Толщина стенки, м', '0.002'),
    'slot-epsilon-r': ('Относительная ε', '1'), 'slot-mu-r': ('Относительная μ', '1'),
    'parallel': ('Потоки OpenMP', '2'), 'vtk-every': ('Сохранять каждый N-й шаг', '5'),
}
PATHS = {'mesh': 'Основная сетка', 'closed-mesh': 'Закрытая сетка',
         'open-mesh': 'Открытая сетка', 'aperture-map': 'Связь этапов (JSON)',
         'slot-cells': 'Ячейки щели (DAT)'}
TASKS = [('Информация о сетке', 'mesh-info', ''), ('Рассеяние', 'scattering', ''),
         ('ЭПР', 'rcs', ''), ('Молния: замкнутое тело', 'lightning', 'closed'),
         ('Молния: два этапа', 'lightning', 'two-stage'),
         ('Молния: ячейки щели', 'lightning', 'slot-cells')]


def defaults():
    return dict(values={k: v[1] for k, v in FIELDS.items()},
                paths={k: '' for k in PATHS}, task=1, polarization='horizontal',
                excitation='lightning-current', shunt=False, dual=True, vtk=True, matrices=False)


def read_mesh(path):
    """ASCII Gmsh 4.1; keep block/file ordering exactly as peec_mesh.f90."""
    lines = Path(path).read_text(encoding='utf-8-sig').splitlines()
    def section(name):
        a, b = lines.index('$'+name), lines.index('$End'+name)
        return lines[a+1:b]
    fmt = section('MeshFormat')[0].split()
    if not (4.1 <= float(fmt[0]) < 5) or fmt[1] != '0':
        raise ValueError('Поддерживается только ASCII Gmsh 4.1.')
    tokens = iter(' '.join(section('Nodes')).split())
    blocks, total, _, _ = (int(next(tokens)) for _ in range(4))
    tags, points = [], []
    for _ in range(blocks):
        dim, entity, parametric, count = (int(next(tokens)) for _ in range(4))
        if parametric:
            raise ValueError('Параметрические узлы не поддерживаются решателем.')
        tags.extend(int(next(tokens)) for _ in range(count))
        points.extend([float(next(tokens)) for _ in range(3)] for _ in range(count))
    if len(tags) != total or len(set(tags)) != total or total == 0 or min(tags) <= 0:
        raise ValueError('Некорректное число или теги узлов.')
    xyz = np.array(points)
    if not np.isfinite(xyz).all():
        raise ValueError('Координаты должны быть конечными.')
    index = {tag: i for i, tag in enumerate(tags)}
    elements = iter(section('Elements'))
    blocks, total, _, _ = map(int, next(elements).split())
    quads, count_all = [], 0
    for _ in range(blocks):
        dim, entity, kind, count = map(int, next(elements).split())
        for _ in range(count):
            row = list(map(int, next(elements).split()))
            count_all += 1
            if kind == 3:
                if len(row) != 5 or len(set(row[1:])) != 4:
                    raise ValueError('Некорректный четырёхугольник.')
                quads.append([index[n] for n in row[1:]])
    if count_all != total or not quads or len(set(np.ravel(quads))) != len(xyz):
        raise ValueError('Нужны четырёхугольники без изолированных узлов.')
    return xyz, np.array(quads, dtype=int), np.array(tags)


def pulse(time, k, alpha, beta, delay=0.):
    t = np.maximum(np.asarray(time) - delay, 0.)
    return k * (np.exp(-alpha*t) - np.exp(-beta*t))


def number(s):
    value = float(str(s).replace('D', 'e').replace('d', 'e'))
    if not math.isfinite(value):
        raise ValueError('Требуется конечное число.')
    return value


def active_keys(state):
    task, case = TASKS[state['task']][1:]
    keys = ['parallel']
    if task in ('scattering', 'rcs'):
        keys += ['field-amplitude', 'theta', 'phi', 'phase']
        keys += ['frequency'] if task == 'scattering' else [
            'ka', 'characteristic-length', 'rcs-observation-theta',
            'rcs-phi-start', 'rcs-phi-end', 'rcs-phi-step', 'rcs-order']
    if task == 'lightning':
        keys += ['dt', 't-end', 'time-order', 'vtk-every']
        if state['excitation'] == 'lightning-current':
            keys += ['lightning-K', 'lightning-alpha', 'lightning-beta',
                     'lightning-delay', 'strike-node', 'return-node']
        else:
            keys += ['pulse-K', 'pulse-alpha', 'pulse-beta', 'theta', 'phi']
        if case == 'slot-cells':
            keys += ['slot-width', 'slot-wall-span', 'slot-wall-thickness', 'slot-epsilon-r', 'slot-mu-r']
        if state['shunt'] and case != 'closed':
            keys += ['shunt-R', 'shunt-a-node', 'shunt-b-node']
    return keys


def required_paths(state):
    task, case = TASKS[state['task']][1:]
    if task != 'lightning':
        return ['mesh']
    if case == 'closed':
        return ['closed-mesh']
    return ['open-mesh', 'slot-cells'] if case == 'slot-cells' else [
        'closed-mesh', 'open-mesh', 'aperture-map']


def resolve(root, path):
    if not path.strip():
        raise ValueError('Не указан путь к файлу.')
    p = Path(path).expanduser()
    return (p if p.is_absolute() else Path(root)/p).resolve()


def validate(state, root):
    vals = {}
    for key in active_keys(state):
        try:
            vals[key] = number(state['values'][key])
        except (ValueError, KeyError) as e:
            raise ValueError(f'{FIELDS[key][0]}: неверное значение.') from e
    integers = {'parallel', 'vtk-every', 'time-order', 'rcs-order',
                'strike-node', 'return-node', 'shunt-a-node', 'shunt-b-node'}
    for key, val in vals.items():
        if key in integers and (not val.is_integer() or val < 0 or val > 2147483646):
            raise ValueError(f'{FIELDS[key][0]}: требуется неотрицательное целое в диапазоне INTEGER(4).')
    positive = {'parallel', 'vtk-every', 'frequency', 'ka', 'characteristic-length', 'rcs-phi-step',
                'dt', 't-end', 'lightning-alpha', 'lightning-beta', 'pulse-alpha', 'pulse-beta',
                'slot-width', 'slot-wall-span', 'slot-epsilon-r', 'slot-mu-r'}
    nonnegative = {'field-amplitude', 'shunt-R', 'lightning-delay', 'slot-wall-thickness'}
    for key in positive & vals.keys():
        if vals[key] <= 0:
            raise ValueError(f'{FIELDS[key][0]} должно быть больше нуля.')
    for key in nonnegative & vals.keys():
        if vals[key] < 0:
            raise ValueError(f'{FIELDS[key][0]} не должно быть отрицательным.')
    if 'time-order' in vals and vals['time-order'] not in (1, 2):
        raise ValueError('Выберите схему времени 1 или 2.')
    if 'dt' in vals:
        n = vals['t-end']/vals['dt']
        if n < 1 or n > 2147483646 or abs(n-round(n)) > 1e-10*max(1, n):
            raise ValueError('Время окончания должно быть целым числом шагов dt, не менее одного.')
    if 'rcs-order' in vals:
        if vals['rcs-order'] not in (2, 3, 4, 5, 6, 8, 10, 12, 16):
            raise ValueError('Неподдерживаемый порядок квадратуры ЭПР.')
        n = (vals['rcs-phi-end']-vals['rcs-phi-start'])/vals['rcs-phi-step']
        if not 0 <= n <= 1e6 or vals['field-amplitude'] <= 0:
            raise ValueError('Проверьте диапазон ЭПР и ненулевую амплитуду поля.')
    meshes = {}
    for key in required_paths(state):
        path = resolve(root, state['paths'][key])
        if not path.is_file():
            raise ValueError(f'{PATHS[key]}: файл не найден: {path}')
        if 'mesh' in key:
            meshes[key] = read_mesh(path)
    case = TASKS[state['task']][2]
    for a, b, meshkey in [('strike-node', 'return-node', 'open-mesh' if case == 'slot-cells' else 'closed-mesh'),
                           ('shunt-a-node', 'shunt-b-node', 'open-mesh')]:
        if a in vals:
            if vals[a] == vals[b] or max(vals[a], vals[b]) >= len(meshes[meshkey][0]):
                raise ValueError(f'{FIELDS[a][0]} и {FIELDS[b][0]}: нужны разные узлы выбранной сетки.')
    return vals


def quote(value):
    value = str(value)
    if any(c in value for c in '\n\r"'):
        raise ValueError('Путь содержит кавычку или перевод строки.')
    return '"'+value+'"'


def config_text(state, root, output):
    vals = validate(state, root)
    task, case = TASKS[state['task']][1:]
    rows = ['# Generated by PEEC Studio', '--task '+task, '--physics quasistatic']
    rows += ['--'+key+' '+quote(resolve(root, state['paths'][key])) for key in required_paths(state)]
    rows += ['--'+key+' '+format(value, '.17g') for key, value in vals.items()]
    if task in ('scattering', 'rcs') or state['excitation'] == 'incident-pulse':
        rows += ['--polarization '+state['polarization']]
    if task == 'lightning':
        rows += ['--lightning-case '+case, '--excitation '+state['excitation']]
    if task == 'rcs':
        rows += ['--rcs-use-dual' if state['dual'] else '--rcs-no-dual',
                 '--output '+quote(Path(output)/'rcs.csv')]
    if state['vtk']:
        rows += ['--vtk']
    rows += ['--vtk-dir '+quote(output)]
    if state['matrices']:
        rows += ['--save-matrices', '--matrix-dir '+quote(Path(output)/'matrices')]
    return '\n'.join(rows)+'\n'


def read_table(path, columns):
    with Path(path).open(encoding='utf-8-sig', newline='') as f:
        reader = csv.DictReader(f)
        if not set(columns) <= set(reader.fieldnames or []):
            raise ValueError('Ожидаются столбцы: '+', '.join(columns))
        data = np.array([[float(row[c]) for c in columns] for row in reader], dtype=float)
    if not data.size:
        raise ValueError('Файл не содержит данных.')
    if not np.isfinite(data[:, 0]).all() or np.isnan(data).any():
        raise ValueError('Файл содержит некорректные значения.')
    return data


def shunt_series(path):
    path = Path(path)
    if path.suffix.lower() == '.csv':
        result = read_table(path, ['time_s', 'shunt_current_A'])
    else:
        match = re.match(r'(.*_shunt)_\d+\.dat$', path.name)
        files = sorted(path.parent.glob(match[1]+'_*.dat')) if match else [path]
        values = []
        for f in files:
            data = dict(line.split(maxsplit=1) for line in f.read_text().splitlines() if line.strip())
            values.append([float(data['time']), float(data['shunt_current_A'])])
        result = np.array(values)
    if not np.isfinite(result).all():
        raise ValueError('Некорректные значения тока шунта.')
    return result[np.argsort(result[:, 0])]


def result_series(path):
    path = Path(path)
    match = re.match(r'(.+)_(\d+)\.vtk$', path.name)
    if not match:
        return [path]
    return sorted(p for p in path.parent.glob(match[1]+'_*.vtk')
                  if re.fullmatch(re.escape(match[1])+r'_\d+\.vtk', p.name))


def save_project(path, state):
    Path(path).write_text(json.dumps({'version': 1, 'state': state}, ensure_ascii=False, indent=2), encoding='utf-8')


def load_project(path):
    data = json.loads(Path(path).read_text(encoding='utf-8'))
    if data.get('version') != 1:
        raise ValueError('Неизвестная версия проекта.')
    state = defaults()
    incoming = data['state']
    for key in ('values', 'paths'):
        if not isinstance(incoming.get(key), dict):
            raise ValueError('Неверная структура проекта.')
        state[key].update({k: str(v) for k, v in incoming[key].items() if k in state[key]})
    for key in ('task', 'polarization', 'excitation', 'shunt', 'dual', 'vtk', 'matrices'):
        if key in incoming:
            state[key] = incoming[key]
    if type(state['task']) is not int or not 0 <= state['task'] < len(TASKS):
        raise ValueError('Неизвестная задача.')
    if state['polarization'] not in ('horizontal', 'vertical') or state['excitation'] not in ('lightning-current', 'incident-pulse'):
        raise ValueError('Неизвестное возбуждение.')
    for key in ('shunt', 'dual', 'vtk', 'matrices'):
        if type(state[key]) is not bool:
            raise ValueError('Неверный тип переключателя.')
    return state


def import_config(path, root):
    """Import solver configs without silently dropping unsupported options."""
    lexer = shlex.shlex(Path(path).read_text(encoding='utf-8-sig'), posix=True)
    lexer.whitespace_split = True
    lexer.escape = ''  # Windows backslashes are literal, as in the Fortran tokenizer.
    args = list(lexer)
    state = defaults()
    state['vtk'] = False
    aliases = {'E': 'field-amplitude', 'a': 'characteristic-length', 'obs-theta': 'rcs-observation-theta',
               'phi-start': 'rcs-phi-start', 'phi-end': 'rcs-phi-end', 'phi-step': 'rcs-phi-step'}
    flags = {'vtk': ('vtk', True), 'save-matrices': ('matrices', True),
             'rcs-dual': ('dual', True), 'rcs-use-dual': ('dual', True), 'rcs-no-dual': ('dual', False)}
    task, case, coordinates = 'mesh-info', 'closed', {}
    i = 0
    while i < len(args):
        key = args[i].removeprefix('--')
        key = aliases.get(key, key)
        if key in flags:
            k, v = flags[key]
            state[k] = v
            i += 1
            continue
        count = 3 if key in ('strike', 'return', 'shunt-a', 'shunt-b') else 1
        values = args[i+1:i+1+count]
        if len(values) != count:
            raise ValueError('Нет значения параметра --'+key)
        value = values[0]
        if key in state['values']:
            state['values'][key] = value
            if key == 'shunt-R':
                state['shunt'] = True
        elif key in state['paths']:
            state['paths'][key] = str(resolve(root, value))
        elif key == 'task':
            task = value
        elif key == 'lightning-case':
            case = value
        elif key in ('polarization', 'excitation'):
            state[key] = value
        elif key == 'physics':
            if value != 'quasistatic':
                raise ValueError('Поддерживается только quasistatic.')
        elif key in ('output', 'matrix-dir', 'vtk-dir'):
            pass  # GUI deliberately creates a fresh per-run output directory.
        elif key in ('strike', 'return', 'shunt-a', 'shunt-b'):
            coordinates[key+'-node'] = np.array([number(x) for x in values])
        else:
            raise ValueError('Неизвестный параметр --'+key)
        i += count+1
    matches = [i for i, t in enumerate(TASKS) if t[1] == task and (task != 'lightning' or t[2] == case)]
    if not matches:
        raise ValueError('Неизвестная задача / вариант молнии.')
    state['task'] = matches[0]
    if state['polarization'] not in ('horizontal', 'vertical') or state['excitation'] not in ('lightning-current', 'incident-pulse'):
        raise ValueError('Неизвестное возбуждение.')
    fallback = 'open-mesh' if case == 'slot-cells' else 'closed-mesh'
    if task == 'lightning' and not state['paths'][fallback]:
        state['paths'][fallback] = state['paths']['mesh']
    for key, coord in coordinates.items():
        meshkey = 'open-mesh' if key.startswith('shunt') else fallback
        xyz = read_mesh(state['paths'][meshkey])[0]
        state['values'][key] = str(int(np.argmin(np.linalg.norm(xyz-coord, axis=1))))
    return state
