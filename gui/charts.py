"""Publication-style Matplotlib figures with Cyrillic and explicit physical units."""
import numpy as np
from matplotlib.figure import Figure
from matplotlib import rc_context

STYLE = {'font.family': 'DejaVu Serif', 'font.size': 11, 'axes.labelsize': 12,
         'axes.titlesize': 13, 'axes.linewidth': .8, 'lines.linewidth': 1.8,
         'xtick.direction': 'in', 'ytick.direction': 'in', 'savefig.dpi': 600,
         'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none'}


def draw(figure, x, y, title, xlabel, ylabel, polar=False):
    with rc_context(STYLE):
        figure.clear()
        ax = figure.add_subplot(111, projection='polar' if polar else None)
        ax.plot(x, y, color='#155e75')
        ax.set_title(title, pad=18)
        ax.set_xlabel(xlabel, labelpad=10)
        ax.set_ylabel(ylabel, labelpad=12)
        ax.grid(True, color='#b8c1cc', alpha=.55, linewidth=.6)
        if not polar:
            ax.minorticks_on()
            ax.tick_params(which='both', top=True, right=True)
            ax.ticklabel_format(axis='y', style='sci', scilimits=(-3, 4), useMathText=True)
        else:
            ax.set_theta_zero_location('E')
        return ax


def rcs_figure(figure, data, mode, title='Диаграмма эффективной площади рассеяния'):
    # Input columns: phi_deg, sigma_m2. Never replace zero by an arbitrary dB floor.
    if not np.isfinite(data).all() or np.any(data[:, 1] < 0):
        raise ValueError('ЭПР должна быть конечной и неотрицательной.')
    x, y = data[:, 0], data[:, 1]
    if mode == 'db':
        mask = y > 0
        if not mask.any():
            raise ValueError('Все значения ЭПР равны нулю: выберите линейный масштаб.')
        y = np.full_like(y, np.nan)
        y[mask] = 10*np.log10(data[mask, 1])
    return draw(figure, np.deg2rad(x) if mode == 'polar' else x, y, title,
                '' if mode == 'polar' else 'Угол наблюдения φ, град',
                'ЭПР σ, дБм²' if mode == 'db' else 'ЭПР σ, м²', mode == 'polar')


def export(figure, path):
    with rc_context(STYLE):
        figure.savefig(path, dpi=600, bbox_inches='tight', facecolor='white')
