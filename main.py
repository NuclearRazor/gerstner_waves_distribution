from pyqtgraph.Qt import QtCore, QtGui
import pyqtgraph.opengl as gl
import pyqtgraph as pg
import numpy as np
import sys

from config import ArgProcessor

#import math core constants and functions
from algorithms import CARTESIAN_POINTS_DIMENSION
from algorithms import CartesianDistibutionpoints
from algorithms import b_n_series
from algorithms import VXCartesianGerstner
from algorithms import VYCartesianGerstner
from algorithms import CartesianOXtoCylindricalU
from algorithms import CartesianOYtoCylindricalV
from algorithms import Grid
from algorithms import DFT
from algorithms import PhillipsPotential
from algorithms import Gaussian
from algorithms import THETA
from algorithms import SURFACE_DIM

#miscellaneous options
CAMERA_DISTANCE = 850
OCEAN_COLOR = 66

#example
#python main.py -w 1.5 -d 100 -a 0.60 -l 0.006

class Visualizer(object):
    def __init__(self):

        super(Visualizer, self).__init__()

        self.traces = dict()
        self.app = QtGui.QApplication(sys.argv)
        self.w = gl.GLViewWidget()

        _pxl_density = gl.GLViewWidget().devicePixelRatio()
        _width = 2*self.w.width() + _pxl_density
        _height = 2*self.w.height() + _pxl_density

        print('Screen width: {}'.format(_width))
        print('Screen height: {}'.format(_height))

        self.w.opts['distance'] = CAMERA_DISTANCE
        self.w.opts['azimuth'] = np.tan(_width/_height)
        self.w.opts['elevation'] = np.tan(_width/_height)
        self.w.opts['bgcolor'] = pg.glColor((249, 97)) #background color

        self.w.setWindowTitle('Statistical distribution of big & small wave amplitudes')
        self.w.setGeometry(30, 40, _width, _height)

        self.w.show()

        self.n = CARTESIAN_POINTS_DIMENSION
        self.m = CARTESIAN_POINTS_DIMENSION

        grid_step = 0.05

        print('Grid step: {}'.format(grid_step))

        self.y = np.array([i*SURFACE_DIM for i in np.sort(CartesianDistibutionpoints(grid_step)[1], axis=None, kind='mergesort')])
        self.x = np.array([j*SURFACE_DIM for j in np.sort(CartesianDistibutionpoints(grid_step)[0], axis=None, kind='mergesort')])

        self.y = np.linspace(np.min(self.y), np.max(self.y), self.n)
        self.x = np.linspace(np.min(self.x), np.max(self.x), self.m)

        self.bound_x = np.mean(self.x)

        self.grad_x = VXCartesianGerstner(self.x)
        self.grad_y = VYCartesianGerstner(self.y)

        self.dft_x_phase = [x_p.real for x_p in DFT(VXCartesianGerstner(self.x))]
        self.dft_y_phase = [y_p.real for y_p in DFT(VYCartesianGerstner(self.y))]

        self.xoy_grid_width = ((np.max(self.x) - np.min(self.x))/2.0 + (np.max(self.y) - np.min(self.y))/2.0)/2.0

        self.phase = 0.0

        self.z = np.zeros(self.n)

        for i in range(self.n):

            yi = np.array([self.y[i]] * self.m)

            u_current = CartesianOXtoCylindricalU(self.grad_x[i], yi)[i]
            y_plot_c = CartesianOYtoCylindricalV(yi, self.x[i])

            alpha = np.mean([x_p.real for x_p in DFT(PhillipsPotential(self.x[i] + self.phase, yi[0]))])

            grad_u = u_current*np.cos(THETA)*u_current + y_plot_c[i]*np.sin(THETA)
            grad_v = u_current*np.sin(THETA)*u_current - y_plot_c[i]*np.cos(THETA)

            #scalar distibution function with current slizes of vectors as scalar values
            #mini wave generates by phases shifting (as animation "step") of sin and cos product
            z = Grid(u_current, y_plot_c) * alpha * Gaussian() * 1.0 / np.sqrt(2.0) + alpha * Grid(u_current, y_plot_c) + np.cos(np.pi*(2.0*grad_v - 0.5) * b_n_series(grad_u, self.phase)) + \
                np.cos(self.x[i] - self.phase) * np.sin(self.y[i] + self.phase)

            self.z[i] = z[i]

            points = np.vstack([self.x, y_plot_c, self.z]).transpose()
            self.traces[i] = gl.GLLinePlotItem(pos = points, color=pg.glColor((OCEAN_COLOR, 100)), width=self.xoy_grid_width*2, antialias=False)
            self.w.addItem(self.traces[i])

    def start(self):
        if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
            QtGui.QApplication.instance().exec_()

    def set_plotdata(self, name, points, color, width):
        self.traces[name].setData(pos=points, color=color, width=width)

    def update(self):

        for i in range(self.n):

            yi = np.array([self.y[i]] * self.m)

            u_current = CartesianOXtoCylindricalU(self.grad_x[i], yi)[i]
            y_plot_c = CartesianOYtoCylindricalV(yi, self.x[i])

            alpha = np.mean([x_p.real for x_p in DFT(PhillipsPotential(self.x[i] + self.phase, yi[0]))])

            grad_u = u_current*np.cos(THETA)*u_current + y_plot_c[i]*np.sin(THETA)
            grad_v = u_current*np.sin(THETA)*u_current - y_plot_c[i]*np.cos(THETA)

            z = Grid(u_current, y_plot_c) * alpha * Gaussian() * 1.0 / np.sqrt(2.0) + alpha * Grid(u_current, y_plot_c) + np.cos(np.pi*(2.0*grad_v - 0.5) * b_n_series(grad_u, self.phase)) + \
                np.cos(self.x[i] - self.phase) * np.sin(self.y[i] + self.phase)

            self.z[i] = z[i]

            #print('X = {}\tY = {}\tZ = {}'.format(self.x[i], y_plot_c[i], self.z[i]))

            points = np.vstack([self.x, y_plot_c, self.z]).transpose()
            self.set_plotdata(name=i, points=points, color=pg.glColor((OCEAN_COLOR, 100)), width=self.xoy_grid_width*2)

            #diffusion ~= -E_system velocity
            self.phase -= .004

    def animation(self):
        timer = QtCore.QTimer()
        timer.timeout.connect(self.update)
        timer.start(100)
        self.start()


if __name__ == '__main__':
    v = Visualizer()
    v.animation()
