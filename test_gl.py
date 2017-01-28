# from pymol.opengl.gl import *
from pymol.callback import Callback
from pymol import cmd

# class myCallback(Callback):
#    def __call__(self):
#       glBegin(GL_LINES)
#       glColor3f(1.0,1.0,1.0)
#       glVertex3f(0.0,0.0,0.0)
#       glVertex3f(1.0,0.0,0.0)
#       glVertex3f(0.0,0.0,0.0)
#       glVertex3f(0.0,2.0,0.0)
#       glVertex3f(0.0,0.0,0.0)
#       glVertex3f(0.0,0.0,3.0)
#       glEnd()
#    def get_extent(self):
#       return [[0.0,0.0,0.0],[1.0,2.0,3.0]]

# cmd.load_callback(myCallback(),'gl01')


def random_walk(N):
    x = 10.0/sqrt(N) * (random.rand(N+1).astype('f')-0.5).cumsum()
    y = 10.0/sqrt(N) * (random.rand(N+1).astype('f')-0.5).cumsum()
    z = 10.0/sqrt(N) * (random.rand(N+1).astype('f')-0.5).cumsum()
    x -= mean(x)
    y -= mean(y)
    z -= mean(z)
    walk = vstack((x,y,z)).transpose()
    return walk / amax(walk)



import time
from OpenGLContext import testingcontext
BaseContext = testingcontext.getInteractive()
from OpenGL.GL import *
from OpenGL.arrays import vbo
from OpenGLContext.arrays import *
from OpenGL.GL import shaders


N = 10000000

class myCallback( Callback ):
    """Creates a simple vertex shader..."""
    def __init__(self):
        mygeom = random_walk(N) * 5
        self.vbo = vbo.VBO(mygeom)
    def __call__(self):
        try:
            self.vbo.bind()
            try:
                glEnableClientState(GL_VERTEX_ARRAY);
                glVertexPointerf( self.vbo )
                glDrawArrays(GL_LINE_STRIP, 0, N+1)
            finally:
                self.vbo.unbind()
                glDisableClientState(GL_VERTEX_ARRAY);
        finally:
            pass

cmd.load_callback(myCallback(), 'gl01')
