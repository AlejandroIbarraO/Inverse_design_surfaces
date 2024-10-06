import polyscope as ps
import polyscope.imgui as psim
import numpy as np
import tkinter as tk
import matplotlib.pyplot as plt
from tkinter import filedialog
from scipy import spatial
from igl import igl
from local_global import local_global
from phasor import phasor
from shapely.geometry.polygon import Polygon
from shapely.geometry import Point

class pattern:
    def __init__(self,phasor_position,phasor_direction,period):
        self.pp = phasor_position
        self.pd = phasor_direction
        self.T = period
        self.phasor = phasor()
        self.phasor.phasor_compute(self.pp,self.pd,1/self.T,self.T)
    
    def 