"""
平面Stewart平台正向运动学求解与可视化
给定三个支杆的长度，找到平台的位置
输入：三角形平台的三条边长$L_1, L_2, L_3$，边$L_2$与边$L_3$的夹角$\gamma$，
      三个支杆的长度$p_1, p_2, p_3$，支杆$p_2$的横坐标$x_1$，支杆$p_3$的坐标$(x_2, y_2)$
输出：三角形平台的坐标$(x, y)$和边$L_3$与水平方向的夹角$\theta$
"""

import math

import sys
from PyQt5.QtCore import Qt
from PyQt5.QtCore import QPointF
from PyQt5.QtGui import QPolygonF
from PyQt5.QtGui import QPen, QBrush, QColor

from PyQt5.QtWidgets import QWidget
from PyQt5.QtWidgets import QLabel
from PyQt5.QtWidgets import QSlider
from PyQt5.QtWidgets import QMainWindow
from PyQt5.QtWidgets import QGraphicsView
from PyQt5.QtWidgets import QGraphicsScene
from PyQt5.QtWidgets import QHBoxLayout, QVBoxLayout
from PyQt5.QtWidgets import QGraphicsPolygonItem, QGraphicsLineItem, QGraphicsEllipseItem
from PyQt5.QtWidgets import QApplication


class StewartPlatform():
    def __init__(self, L1, L2, L3, gamma, p1, p2, p3, x1, x2, y2):
        self.L1 = L1
        self.L2 = L2
        self.L3 = L3
        self.gamma = gamma
        self.p1 = p1
        self.p2 = p2
        self.p3 = p3
        self.x1 = x1
        self.x2 = x2
        self.y2 = y2

    def set_p(self, p1, p2, p3):
        self.p1 = p1
        self.p2 = p2
        self.p3 = p3

    def func(self, theta):
        A2 = self.L3 * math.cos(theta) - self.x1
        B2 = self.L3 * math.sin(theta)
        A3 = self.L2 * math.cos(theta + self.gamma) - self.x2
        B3 = self.L2 * math.sin(theta + self.gamma) - self.y2

        D = 2 * (A2 * B3 - B2 * A3)
        E = self.p2 ** 2 - self.p1 ** 2 - A2 ** 2 - B2 ** 2
        F = self.p3 ** 2 - self.p1 ** 2 - A3 ** 2 - B3 ** 2

        N1 = B3 * E - B2 * F
        N2 = -A3 * E + A2 * F

        f = N1 ** 2 + N2 ** 2 - self.p1 ** 2 * D ** 2

        return f

    def func_prime(self, theta):
        A2 = self.L3 * math.cos(theta) - self.x1
        B2 = self.L3 * math.sin(theta)
        A3 = self.L2 * math.cos(theta + self.gamma) - self.x2
        B3 = self.L2 * math.sin(theta + self.gamma) - self.y2

        D = 2 * (A2 * B3 - B2 * A3)
        E = self.p2 ** 2 - self.p1 ** 2 - A2 ** 2 - B2 ** 2
        F = self.p3 ** 2 - self.p1 ** 2 - A3 ** 2 - B3 ** 2

        N1 = B3 * E - B2 * F
        N2 = -A3 * E + A2 * F

        A2_prime = -self.L3 * math.sin(theta)
        B2_prime = self.L3 * math.cos(theta)
        A3_prime = -self.L2 * math.sin(theta + self.gamma)
        B3_prime = self.L2 * math.cos(theta + self.gamma)

        E_prime = -2 * A2 * A2_prime - 2 * B2 * B2_prime
        F_prime = -2 * A3 * A3_prime - 2 * B3 * B3_prime

        D_prime = 2 * (A2_prime * B3 + A2 * B3_prime - B2_prime * A3 - B2 * A3_prime)

        N1_prime = B3_prime * E + B3 * E_prime - B2_prime * F - B2 * F_prime
        N2_prime = -(A3_prime * E + A3 * E_prime) + A2_prime * F + A2 * F_prime

        f_prime = 2 * N1 * N1_prime + 2 * N2 * N2_prime - 2 * self.p1 ** 2 * D * D_prime

        return f_prime

    def get_config(self):
        config = {"L1": self.L1,
                  "L2": self.L2,
                  "L3": self.L3,
                  "gamma": self.gamma,
                  "x1": self.x1,
                  "x2": self.x2,
                  "y2": self.y2}
        return config

    def newton_raphson(self, f, f_prime, x0, k, tol=1e-8):
        xc = [x0]

        while k and f(xc[-1]) > tol:
            xc.append(xc[-1] - f(xc[-1]) / f_prime(xc[-1]))
            k = k - 1

        return xc[-1]

    def stewart_forward(self, p1, p2, p3):
        self.set_p(p1, p2, p3)

        theta = self.newton_raphson(self.func, self.func_prime, math.pi / 2, 25)

        A2 = self.L3 * math.cos(theta) - self.x1
        B2 = self.L3 * math.sin(theta)
        A3 = self.L2 * math.cos(theta + self.gamma) - self.x2
        B3 = self.L2 * math.sin(theta + self.gamma) - self.y2

        D = 2 * (A2 * B3 - B2 * A3)
        E = self.p2 ** 2 - self.p1 ** 2 - A2 ** 2 - B2 ** 2
        F = self.p3 ** 2 - self.p1 ** 2 - A3 ** 2 - B3 ** 2

        N1 = B3 * E - B2 * F
        N2 = -A3 * E + A2 * F

        x = N1 / D
        y = N2 / D

        return x, y, theta


class GraphicsScene(QGraphicsScene):
    def __init__(self, parent=None, stewart=None, radius=6, zoomIn=100):
        super(GraphicsScene, self).__init__(parent=parent)

        self.stewart = stewart
        self.zoomIn = zoomIn
        self.radius = radius
        self.diameter = radius * 2

        self.triangle = QGraphicsPolygonItem()
        self.triangle.setPen(QPen(QColor(0, 0, 255)))
        self.triangle.setBrush(QBrush(QColor(0, 0, 255, 128)))
        self.addItem(self.triangle)

        self.p1 = QGraphicsLineItem()
        self.p1.setPen(QPen(QColor(0, 0, 0)))
        self.addItem(self.p1)
        self.p2 = QGraphicsLineItem()
        self.p2.setPen(QPen(QColor(0, 0, 0)))
        self.addItem(self.p2)
        self.p3 = QGraphicsLineItem()
        self.p3.setPen(QPen(QColor(0, 0, 0)))
        self.addItem(self.p3)

        self.vertex1 = QGraphicsEllipseItem()
        self.vertex1.setPen(QPen(QColor(0, 0, 255)))
        self.vertex1.setBrush(QBrush(QColor(0, 0, 255, 128)))
        self.addItem(self.vertex1)
        self.vertex2 = QGraphicsEllipseItem()
        self.vertex2.setPen(QPen(QColor(0, 0, 255)))
        self.vertex2.setBrush(QBrush(QColor(0, 0, 255, 128)))
        self.addItem(self.vertex2)
        self.vertex3 = QGraphicsEllipseItem()
        self.vertex3.setPen(QPen(QColor(0, 0, 255)))
        self.vertex3.setBrush(QBrush(QColor(0, 0, 255, 128)))
        self.addItem(self.vertex3)

        self.anchor1 = QGraphicsEllipseItem()
        self.anchor1.setPen(QPen(QColor(0, 0, 0)))
        self.anchor1.setBrush(QBrush(QColor(0, 0, 0, 128)))
        self.addItem(self.anchor1)
        self.anchor2 = QGraphicsEllipseItem()
        self.anchor2.setPen(QPen(QColor(0, 0, 0)))
        self.anchor2.setBrush(QBrush(QColor(0, 0, 0, 128)))
        self.addItem(self.anchor2)
        self.anchor3 = QGraphicsEllipseItem()
        self.anchor3.setPen(QPen(QColor(0, 0, 0)))
        self.anchor3.setBrush(QBrush(QColor(0, 0, 0, 128)))
        self.addItem(self.anchor3)

        self.update(math.sqrt(5), math.sqrt(5), math.sqrt(5))

    def update(self, p1, p2, p3):
        x, y, theta = stewart.stewart_forward(p1, p2, p3)
        x = x * self.zoomIn
        y = y * self.zoomIn
        config = stewart.get_config()
        L2 = config["L2"]
        L3 = config["L3"]
        gamma = config["gamma"]
        x1 = config["x1"]
        x2 = config["x2"]
        y2 = config["y2"]

        x1 = x1 * self.zoomIn
        x2 = x2 * self.zoomIn
        y2 = y2 * self.zoomIn
        L2 = L2 * self.zoomIn
        L3 = L3 * self.zoomIn

        px1, py1 = x, -y
        px2, py2 = x + L2 * math.cos(theta + gamma), -(y + L2 * math.sin(theta + gamma))
        px3, py3 = x + L3 * math.cos(theta), -(y + L3 * math.sin(theta))

        self.vertex1.setRect(px1 - self.radius, py1 - self.radius, self.diameter, self.diameter)
        self.vertex2.setRect(px2 - self.radius, py2 - self.radius, self.diameter, self.diameter)
        self.vertex3.setRect(px3 - self.radius, py3 - self.radius, self.diameter, self.diameter)

        self.anchor1.setRect(-self.radius, -self.radius, self.diameter, self.diameter)
        self.anchor2.setRect(x1 - self.radius, -self.radius, self.diameter, self.diameter)
        self.anchor3.setRect(x2 - self.radius, -y2 - self.radius, self.diameter, self.diameter)

        self.triangle.setPolygon(QPolygonF([QPointF(px1, py1), QPointF(px2, py2), QPointF(px3, py3)]))

        self.p1.setLine(0, 0, x, -y)
        self.p3.setLine(x2, -y2, px2, py2)
        self.p2.setLine(x1, 0, px3, py3)


class CentralWidget(QWidget):
    def __init__(self, stewart=None, factor=1.0):
        super(CentralWidget, self).__init__()

        self.factor = factor
        self.setLayout(QHBoxLayout())

        self.p1_slider = QSlider()
        self.p1_slider.setMinimum(0)
        self.p1_slider.setMaximum(100)
        self.p1_slider.setValue(math.sqrt(5) * self.factor)
        self.p1_slider.setOrientation(Qt.Horizontal)
        self.p1_slider.valueChanged.connect(self.update_p1)
        self.p1_label = QLabel("p1=%.3f" % math.sqrt(5))
        self.p1_layout = QHBoxLayout(self)
        self.p1_layout.addWidget(self.p1_slider)
        self.p1_layout.addWidget(self.p1_label)

        self.p2_slider = QSlider()
        self.p2_slider.setMinimum(0)
        self.p2_slider.setMaximum(100)
        self.p2_slider.setValue(math.sqrt(5) * self.factor)
        self.p2_slider.setOrientation(Qt.Horizontal)
        self.p2_slider.valueChanged.connect(self.update_p2)
        self.p2_label = QLabel("p2=%.3f" % math.sqrt(5))
        self.p2_layout = QHBoxLayout(self)
        self.p2_layout.addWidget(self.p2_slider)
        self.p2_layout.addWidget(self.p2_label)

        self.p3_slider = QSlider()
        self.p3_slider.setMinimum(0)
        self.p3_slider.setMaximum(100)
        self.p3_slider.setValue(math.sqrt(5) * self.factor)
        self.p3_slider.setOrientation(Qt.Horizontal)
        self.p3_slider.valueChanged.connect(self.update_p3)
        self.p3_label = QLabel("p3=%.3f" % math.sqrt(5))
        self.p3_layout = QHBoxLayout(self)
        self.p3_layout.addWidget(self.p3_slider)
        self.p3_layout.addWidget(self.p3_label)

        self.sliderLayout = QVBoxLayout(self)
        self.sliderLayout.addLayout(self.p1_layout)
        self.sliderLayout.addLayout(self.p2_layout)
        self.sliderLayout.addLayout(self.p3_layout)

        self.scene = GraphicsScene(stewart=stewart)
        self.view = QGraphicsView(self)
        self.view.setScene(self.scene)

        self.layout().addWidget(self.view)
        self.layout().addLayout(self.sliderLayout)

    def update_p1(self, value):
        self.p1_label.setText("p1=%.3f" % (value / self.factor))
        self.update()

    def update_p2(self, value):
        self.p2_label.setText("p2=%.3f" % (value / self.factor))
        self.update()

    def update_p3(self, value):
        self.p3_label.setText("p3=%.3f" % (value / self.factor))
        self.update()

    def update(self):
        p1 = self.p1_slider.value() / self.factor
        p2 = self.p2_slider.value() / self.factor
        p3 = self.p3_slider.value() / self.factor
        self.scene.update(p1, p2, p3)


class MainWindow(QMainWindow):
    def __init__(self, stewart=None):
        super(MainWindow, self).__init__()
        self.view = CentralWidget(stewart, factor=10)
        self.setCentralWidget(self.view)
        self.statusBar().showMessage("Ready!")


if __name__ == '__main__':
    app = QApplication(sys.argv)

    stewart = StewartPlatform(2, math.sqrt(2), math.sqrt(2), math.pi / 2, math.sqrt(5),
                              math.sqrt(5), math.sqrt(5), 4, 0, 4)

    mainwindow = MainWindow(stewart)
    mainwindow.setWindowTitle("2D Stewart Platform")
    mainwindow.show()
    app.exec_()
