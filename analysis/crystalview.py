from scipy.stats import maxwell
import pygame
import pygame.color
from pygame.locals import *
import crystalgen
import basisgen
import particles
import helpers
import numpy as np
import math
from PyQt5.QtWidgets import *

key_to_function = {
    pygame.K_LEFT:   (lambda x: x.translateAll([-10,0,0])),
    pygame.K_RIGHT:  (lambda x: x.translateAll([ 10,0,0])),
    pygame.K_DOWN:   (lambda x: x.translateAll([0, 10,0])),
    pygame.K_UP:     (lambda x: x.translateAll([0,-10,0])),
    pygame.K_EQUALS: (lambda x: x.scaleAll(1.25)),
    pygame.K_MINUS:  (lambda x: x.scaleAll( 0.8)),
    pygame.K_q:      (lambda x: x.rotateAll([ 0.0005,0,0])),
    pygame.K_w:      (lambda x: x.rotateAll([-0.0005,0,0])),
    pygame.K_a:      (lambda x: x.rotateAll([0, 0.0005,0])),
    pygame.K_s:      (lambda x: x.rotateAll([0,-0.0005,0])),
    pygame.K_z:      (lambda x: x.rotateAll([0,0, 0.0005])),
    pygame.K_x:      (lambda x: x.rotateAll([0,0,-0.0005]))}

class Points:
    def __init__(self):
        # These are the points that render on the screen
        self.points = []
        self.custom_points = []

        self.printed = False
        self.points_render = np.zeros((0,4))

        self.colour = (255,255,125)
        self.colour_custom = (255,55,0)
        
        self.transform = np.array([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]])
    
    def update(self, particles):
        r = particles.positions
        r = np.vstack((r, self.custom_points))

        i = np.ones((len(r),1))
        
        self.points_render = np.hstack((r, i))

        trans = np.copy(self.transform)
        self.transform = np.array([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]])
        self.applyTransform(trans)
        return

    def addPoint(self, point):
        self.custom_points.append(point)

    def draw(self, screen, nodeRadius):
        num = len(self.points_render)
        cp = num - len(self.custom_points)
        for i in range(num):
            point = self.points_render[i]
            colour = self.colour;
            if i > cp:
                colour = self.colour_custom
            pygame.draw.circle(screen, colour, (int(point[0]), int(point[1])), nodeRadius, 0)
    
    def applyTransform(self, matrix):
        self.transform = np.dot(self.transform, matrix)
        self.points_render = np.dot(self.points_render,matrix)
    
    def translate(self, dx=0, dy=0, dz=0):
        self.applyTransform(self.translationMatrix(dx,dy,dz))
        
    def rotate(self, rx, ry, rz):
        self.applyTransform(self.rotateMatrix(rx, ry, rz))
        
    def findCentre(self):
        num = len(self.points_render)
        meanX = sum([point[0] for point in self.points_render]) / num
        meanY = sum([point[1] for point in self.points_render]) / num
        meanZ = sum([point[2] for point in self.points_render]) / num
        var = True
        if var:
            return (0,0,0)
        return (meanX, meanY, meanZ)
        
    def scaleMatrix(self,sx=0, sy=0, sz=0):
        return np.array([[sx, 0,  0,  0],
                         [0,  sy, 0,  0],
                         [0,  0,  sz, 0],
                         [0,  0,  0,  1]])
    
    def translationMatrix(self,dx=0, dy=0, dz=0):
        return np.array([[1,0,0,0],
                         [0,1,0,0],
                         [0,0,1,0],
                         [dx,dy,dz,1]])
    
    def rotateXMatrix(self,radians):
        c = np.cos(radians)
        s = np.sin(radians)
        return np.array([[1, 0, 0, 0],
                         [0, c,-s, 0],
                         [0, s, c, 0],
                         [0, 0, 0, 1]])
    def rotateYMatrix(self,radians):
        c = np.cos(radians)
        s = np.sin(radians)
        return np.array([[ c, 0, s, 0],
                         [ 0, 1, 0, 0],
                         [-s, 0, c, 0],
                         [ 0, 0, 0, 1]])
    def rotateZMatrix(self,radians):
        c = np.cos(radians)
        s = np.sin(radians)
        return np.array([[c,-s, 0, 0],
                         [s, c, 0, 0],
                         [0, 0, 1, 0],
                         [0, 0, 0, 1]])
        
    def rotateMatrix(self,rx,ry,rz):
        return self.rotateXMatrix(rx).dot(self.rotateYMatrix(ry).dot(self.rotateZMatrix(rz)))

class PointViewer:
    """ Displays 3D objects on a Pygame screen """

    def __init__(self, width, height):
        self.width = width
        self.height = height
        self.background = (10,10,10)
        self.points = Points()
        self.lines = Points()
        self.particles = particles.Particles()
        self.nodeColour = (255,255,255)
        self.nodeRadius = 4
        self.doTick = True
        self.tick_step = 0.01
        self.outputfile = None
        self.load()
        self.dir = [0,0,1]
        self.axis = [0,0,1]
        pygame.font.init() 
        self.myfont = pygame.font.SysFont('Arial', 30)

    def addRef(self, direction):

        #Normalize the direction.
        ss = math.sqrt(direction[0]**2 + direction[1]**2 + direction[2]**2);
        direction[0]/=ss;
        direction[1]/=ss;
        direction[2]/=ss;

        print(direction)

        im = np.array(self.axis)
        #im = np.array(direction)
        ix = np.array(self.dir)
        #ix = np.array(direction)
        R = helpers.rotate(ix, im)
        #R_inv = R#np.linalg.inv(R)
        R_inv = np.linalg.inv(R)

        ex = np.asarray(np.matmul(R_inv, np.array([1,0,0])))[0]
        ey = np.asarray(np.matmul(R_inv, np.array([0,1,0])))[0]
        ez = np.asarray(np.matmul(R_inv, np.array([0,0,1])))[0]

        #direction = self.axis
        #Convert the direction to the local coordinate ststem
        tmp = [0,0,0]
        tmp[0] = ex[0] * direction[0] + ex[1] * direction[1] + ex[2] * direction[2]
        tmp[1] = ey[0] * direction[0] + ey[1] * direction[1] + ey[2] * direction[2]
        tmp[2] = ez[0] * direction[0] + ez[1] * direction[1] + ez[2] * direction[2]

        for i in range(10):
            p = [0,0,0]
            p[0] = p[0] + tmp[0]*i/1.0
            p[1] = p[1] + tmp[1]*i/1.0
            p[2] = p[2] + tmp[2]*i/1.0
            self.points.addPoint(p)

    def tick(self):
        if self.doTick:
            self.particles.step(self.tick_step)
            self.points.update(self.particles)
            self.outputfile.write(str(self.particles.T())+'\n')
        return

    def save(self):
        if self.particles.steps:
            self.particles.save('crystal.input')
        self.outputfile.close()

    def load(self):
        self.outputfile = open('T.output', 'w')
        size = 4.0786
        self.dir = [0,0,1]
        self.axis = [7,8,8]
        atom = basisgen.Atom(196.967,79)
        #crystal = crystalgen.gen(size, self.dir, self.axis, basisgen.fccBasis(atom), 10, 0.1, -1.75*size)
        #n = 5
        #crystalgen.clearOutOfBounds(crystal, -size * n, size * n, -size * n, size * n)

        self.addRef([-1,-1,1])
        self.addRef([2,-1,1])
        self.addRef([0,1,1])



        #self.addRef([-7,-8,8])

        '''
        self.addRef([1,0,0])
        self.addRef([0,0,1])
        self.addRef([0,1,0])#'''
        
        self.particles.coupling = False
        self.particles.steps = False
        self.particles.load('crystal.input')
        self.points.update(self.particles)

        self.translateAll([self.width/2,self.height/2,0])
        self.scaleAll(15)

    def onEvent(self, event):
        if event.type == pygame.KEYDOWN:
            if event.key in key_to_function:
                key_to_function[event.key](self)
            if event.key == pygame.K_KP8:
                self.particles.couplingMult *= 2
            if event.key == pygame.K_KP2:
                self.particles.couplingMult /= 2
            if event.key == pygame.K_KP6:
                self.particles.latticeMult *= 2
            if event.key == pygame.K_KP4:
                self.particles.latticeMult /= 2
        if event.type == pygame.MOUSEBUTTONDOWN:
            #scroll in
            if event.button == 4:
                self.scaleAll(1.05)
            # scroll out
            elif event.button == 5:
                self.scaleAll(1/1.05)
        if event.type == pygame.MOUSEMOTION:
            if pygame.mouse.get_pressed()[0]:
                rel = event.rel
                if rel[0] != 0:
                    self.rotateAll([0,-0.001*rel[0],0])
                if rel[1] != 0:
                    self.rotateAll([0.001 * rel[1],0,0])
            if pygame.mouse.get_pressed()[2]:
                rel = event.rel
                if rel[0] != 0:
                    self.translateAll([1*rel[0],0,0])
                if rel[1] != 0:
                    self.translateAll([0, 1*rel[1],0])
            if pygame.mouse.get_pressed()[1]:
                rel = event.rel
                if rel[0] != 0:
                    self.rotateAll([0,0,0.001 * rel[0]])

    def display(self):
        self.screen.fill(self.background)

        self.points.draw(self.screen, self.nodeRadius)

        T = math.trunc(self.particles.T() * 10000)/10000
        textsurface = self.myfont.render(str(T), False, (255, 0, 0))
        self.screen.blit(textsurface,(5,0))
        pygame.display.flip()

    def translateAll(self, vector):
        self.temp = Points()
        matrix = self.temp.translationMatrix(*vector)
        self.points.applyTransform(matrix)
        return

    def scaleAll(self, scale):
        self.translateAll([-self.width/2,-self.height/2,0])
        self.temp = Points()
        matrix = self.temp.scaleMatrix(scale,scale,scale)
        shift = self.points.findCentre()
        offset = self.points.translationMatrix(-shift[0],-shift[1],-shift[2])
        self.points.applyTransform(offset)
        self.points.applyTransform(matrix)
        offset = self.points.translationMatrix(shift[0],shift[1],shift[2])
        self.points.applyTransform(offset)
        self.translateAll([self.width/2,self.height/2,0])
        return
    
    def rotateAll(self, vector):
        self.translateAll([-self.width/2, -self.height/2, 0])
        self.temp = Points()
        matrix = self.temp.rotateMatrix(*vector)
        shift = self.points.findCentre()
        offset = self.points.translationMatrix(-shift[0],-shift[1],-shift[2])
        self.points.applyTransform(offset)
        self.points.applyTransform(matrix)
        offset = self.points.translationMatrix(shift[0],shift[1],shift[2])
        self.points.applyTransform(offset)
        self.translateAll([self.width/2, self.height/2, 0])
        return

class App:
    def __init__(self):
        self._running = True
        self._display_surf = None
        self.size = self.weight, self.height = 800, 800
        self.view = PointViewer(800,800)

    def on_init(self):
        pygame.init()
        self._display_surf = pygame.display.set_mode(self.size, pygame.HWSURFACE | pygame.DOUBLEBUF)
        self.view.screen = self._display_surf
        self._running = True
 
    def on_event(self, event):
        self.view.onEvent(event)
        if event.type == pygame.QUIT:
            self._running = False

    def on_loop(self):
        self.view.tick()

    def on_render(self):
        self.view.display()

    def on_cleanup(self):
        self.view.save()
        pygame.quit()
 
    def on_execute(self):
        if self.on_init() == False:
            self._running = False
 
        while( self._running ):
            for event in pygame.event.get():
                self.on_event(event)
            self.on_loop()
            self.on_render()

        self.on_cleanup()

def makeInputBoxes(theApp):
    
    app = QApplication([])
    window = QWidget()
    layout = QGridLayout()

    #Make some buttons
    run = QPushButton('Button')
    def push():
        #Replace with whatever for running safari later.
        print(theApp.view.particles.temp())
    run.clicked.connect(push)
    
    box = QVBoxLayout()
    box.addWidget(run)

    x = 0
    y = 0
    layout.addLayout(box, x, y + 00)
    
    window.setLayout(layout)
    window.show()
    app.exec_()
 
if __name__ == "__main__" :
    theApp = App()
    theApp.on_execute()
