import genomewindow as gw
import numpy as np
from matplotlib.patches import Polygon

class GenePatch():
    def __init__(self, gene_name = None, gene_track = None, x_extents = None, y_extents = None, direction = None, head_length = None, gene_patch_kwargs = {}):
        self.gene_name = gene_name
        self.gene_track = gene_track
        self.gene_track_extents = [self.gene_track.genome_window.x_array[0], self.gene_track.genome_window.x_array[-1]]
        self.dinchperx, self.dinchpery = gw.ddatatodinchaxis(self.gene_track.ax)
        self.x_extents = x_extents
        self.y_extents = y_extents
        self.patch_height = self.dinchpery * abs(self.y_extents[0] - self.y_extents[1])
        self.direction = direction
        self.head_length = head_length
        self.gene_patch_kwargs = gene_patch_kwargs
        self.patch = self.draw()
            
class ArrowGenePatch(GenePatch):
    def draw(self):
        head_length_x = self.head_length / self.dinchperx
        if self.direction is 'left':
            coordinates = [[min(self.x_extents[1],self.gene_track_extents[1]+head_length_x), self.y_extents[0]], # bottom right
                           [min(self.x_extents[1],self.gene_track_extents[1]+head_length_x), self.y_extents[1]], # top right
                           [max(min(self.x_extents[1],self.x_extents[0]+head_length_x), self.gene_track_extents[0]), self.y_extents[1]], # top left
                           [max(self.x_extents[0], self.gene_track_extents[0] - head_length_x), np.mean(self.y_extents)], # triangle point
                           [max(min(self.x_extents[1],self.x_extents[0]+head_length_x), self.gene_track_extents[0]), self.y_extents[0]],] # bottom left
        if self.direction is 'right':
            coordinates = [[max(self.x_extents[0],self.gene_track_extents[0]-head_length_x), self.y_extents[0]], # bottom left
                           [min(max(self.x_extents[0], self.x_extents[1] - head_length_x), self.gene_track_extents[1]), self.y_extents[0]], # bottom right
                           [min(self.x_extents[1], self.gene_track_extents[1] + head_length_x), np.mean(self.y_extents)], # triangle point
                           [min(max(self.x_extents[0], self.x_extents[1] - head_length_x), self.gene_track_extents[1]), self.y_extents[1]], # top right 
                           [max(self.x_extents[0],self.gene_track_extents[0]-head_length_x), self.y_extents[1]],] # top left
        self.coordinates = coordinates
        return self.gene_track.ax.add_patch(Polygon(self.coordinates, **self.gene_patch_kwargs ))