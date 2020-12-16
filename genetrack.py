import genomewindow as gw
import numpy as np

class DoubleStrandGeneTrack():
    def drawgenes(self):
        overlap_mask = np.all([self.genome_window.gene_regions[:,1] <= self.genome_window.window_right, 
                               self.genome_window.gene_regions[:,2] >= self.genome_window.window_left],axis=0)
        self.overlapping_dataframe = self.genome_window.gene_table.loc[overlap_mask]
        self.overlapping_names = self.genome_window.gene_names[overlap_mask]
        self.overlapping_regions = self.genome_window.gene_regions[overlap_mask]
        self.gene_patches = []
        for name, region in zip(self.overlapping_names, self.overlapping_regions):
            strand, start, end = region
            name = str(name)
            if self.genome_window.top_positive is True:
                gene_direction = 'right'
                y_extents = [.05, .95]
                if strand == 1:
                    gene_direction = 'left'
                    y_extents = [-.95, -.05]
                x_extents = [start - self.genome_window.window_zero, end - self.genome_window.window_zero]
                self.gene_patches.append(self.patch_type(gene_name = name, x_extents = x_extents, y_extents = y_extents,
                                                         gene_track = self, direction = gene_direction, head_length = self.head_length, gene_patch_kwargs = self.gene_patch_kwargs))
            if self.genome_window.top_positive is False:
                gene_direction = 'right'
                y_extents = [.05, .95]
                if strand == 0:
                    gene_direction = 'left'
                    y_extents = [-.95, -.05]
                x_extents = [(self.genome_window.window_right - end)   - (self.genome_window.window_right - self.genome_window.window_zero),
                             (self.genome_window.window_right - start) - (self.genome_window.window_right - self.genome_window.window_zero)]
                self.gene_patches.append(self.patch_type(gene_name = name, x_extents = x_extents, y_extents = y_extents,
                                                         gene_track = self, direction = gene_direction, head_length = self.head_length, gene_patch_kwargs = self.gene_patch_kwargs))

    def drawlabels(self):
        for patch in self.gene_patches:
            # try to draw text in the patch
            horizontal_position = np.mean([max(patch.x_extents[0], self.genome_window.x_array[0]), min(patch.x_extents[1], self.genome_window.x_array[-1])])
            text = self.ax.text(horizontal_position, np.mean(patch.y_extents) - abs(patch.y_extents[0]-patch.y_extents[1])*.05, patch.gene_name,
                                va='center', ha='center', fontsize=patch.patch_height*72*.8/(patch.gene_name.count('\n')+1))
            # now check text width against size of gene, if gene box is too small, redraw the gene information below
            bbox = text.get_window_extent(renderer = self.ax.get_figure().canvas.get_renderer()).inverse_transformed(self.ax.transData)
            text_width = bbox.x1 - bbox.x0
            patch_width = min(patch.x_extents[1], self.genome_window.x_array[-1]) - max(patch.x_extents[0], self.genome_window.x_array[0])
            if text_width > patch_width * .85:
                text.remove()
                text = self.ax.text(horizontal_position, -1.05, patch.gene_name,
                                    va='top', ha='center', fontsize=patch.patch_height*72*.8/(patch.gene_name.count('\n')+1))
                self.ax.plot([horizontal_position, horizontal_position], [patch.y_extents[0],-1], c='k', lw=.5)
                
            
            continue
            # does not fit in
            if not(patch.x_extents[0] < self.genome_window.x_array[0] or patch.x_extents[1] > self.genome_window.x_array[-1]):
                self.ax.text(np.mean(patch.x_extents), np.mean(patch.y_extents) - abs(patch.y_extents[0]-patch.y_extents[1])*.05, patch.gene_name,
                             va='center', ha='center', fontsize=patch.patch_height*72*.8/(patch.gene_name.count('\n')+1))
            

    def __init__(self, genome_window = None, patch_type = gw.ArrowGenePatch, head_length = .11, gene_patch_kwargs = {'linewidth':1,'facecolor':'w','edgecolor':'k'}):
        self.genome_window = genome_window
        self.ax = self.genome_window.ax_genes
        self.ax.set_ylim(-1,1)
        self.patch_type = patch_type
        self.head_length = head_length
        self.gene_patch_kwargs = gene_patch_kwargs