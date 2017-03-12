# pip install PyOpenGL OpenGLContext PyVRML97 pydispatcher


from pymol import cmd
from pymol.callback import Callback

import time
from OpenGLContext import testingcontext
BaseContext = testingcontext.getInteractive()
from OpenGL.GL import *
from OpenGL.arrays import vbo
from OpenGLContext.arrays import *
from OpenGL.GL import shaders

import numpy as np


def random_walk(N):
    x = np.random.randn(N + 1, 3).astype('f').cumsum(axis=0)
    x -= x.mean(axis=0)
    return 0.5 * x / x.std()


class myCallback(Callback):
    """Creates a simple vertex shader..."""

    def __init__(self):
        self.rebuild(10 * 1000 * 1000)

    def rebuild(self, nlines):
        self.nlines = nlines
        mygeom = random_walk(nlines) * 10
        self.vbo = vbo.VBO(mygeom)

    def __call__(self):
        try:
            self.vbo.bind()
            try:
                glEnableClientState(GL_VERTEX_ARRAY)
                glVertexPointerf(self.vbo)
                glDrawArrays(GL_LINE_STRIP, 0, self.nlines + 1)
            finally:
                self.vbo.unbind()
                glDisableClientState(GL_VERTEX_ARRAY)
        finally:
            pass


walker = myCallback()
cmd.load_callback(walker, 'gl01')
