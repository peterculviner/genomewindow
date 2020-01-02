import numpy as np
import genomewindow as gw

class GenomeWindow():
    def setpositionbygene(self, name = None, nt_spacer = 0, nt_5 = 0, nt_3 = 0, swap_strand = False, verbose = True):
        # attempt to find the name of the gene provided
        found_genes = [gn for gn in self.gene_names if name in str(gn)]
        if len(found_genes) == 0:
            raise(ValueError('Gene not found.'))
        if len(found_genes) > 1:
            if verbose:
                print ('Multiple genes found matching or partially matching %s. Plotting the first one in the list, complete list is: %s.'%
                        (name, ', '.join(found_genes)))
        strand, start, end = self.gene_regions[np.isin(self.gene_names,found_genes)][0]
        # plot directly from given coordinates
        nt_spacer = int((end - start)*nt_spacer) # convert nt spacer to actual nt rather than gene fraction
        if (strand == 0 and not swap_strand) or (strand == 1 and swap_strand):
            self.top_positive = True
            self.window_zero = start
            self.window_left = start - nt_5 - nt_spacer
            self.window_right = end + nt_3 + nt_spacer
        else:
            self.top_positive = False
            self.window_zero = end
            self.window_right = end + nt_5 + nt_spacer
            self.window_left = start - nt_3 - nt_spacer
        self._setx()
    
    def setpositionbycoordinates(self, coordinates = None, strand = None):
        pass
    
    def plotdatastreams(self, **kwargs):
        """ Description.
        
        Parameters:
        ----------
        axis_n : required
        
        stream_type : required
        
        input_data : required
        
        Returns:
        ----------
        adds plots
        """
        # accepts multiple streams of data or a single stream (all most be the same stream type)
        # kwargs may be assigned to all streams (single value) or to each stream (list of values)
        try:
            iter(kwargs['axis_n']) # check if axes is iterable
            # if no exception, split kwargs and input each into function separately
            for data_set in gw.splitkwargs(kwargs):
                self.plotdatastreams(**data_set)
        except: # treat as a single instance
            stream_type = kwargs.pop('stream_type')
            stream_type(genome_window = self, **kwargs)
        
    
    def _setx(self):
        # if user asked left or right coordinates to beyond genome size, set to edge of genome
        self.window_left = max(self.window_left,0)
        self.window_right = min(self.window_right,len(self.genome)-1)
        # now set up a mask and the shown x-coordinates
        window_mask = np.zeros((2,len(self.genome))).astype(bool)
        if self.top_positive:
            self.x_array = np.arange(self.window_left - self.window_zero, self.window_right - self.window_zero + 1)
            window_mask[np.zeros(self.x_array.shape).astype(int), np.arange(self.window_left, self.window_right+1)] = True
        else:
            self.x_array = np.arange(self.window_zero - self.window_right, self.window_zero - self.window_left + 1)
            window_mask[np.zeros(self.x_array.shape).astype(int)+1, np.arange(self.window_left, self.window_right+1)] = True
        self.window_mask = window_mask
        # draw and label genes
        self.ax_genes.set_xlim(self.x_array[0],self.x_array[-1])
        self.gene_track.drawgenes()
        self.gene_track.drawlabels()
        # draw sequence in sequence track if present
        if self.ax_sequence is not None:
            self.ax_sequence.set_xlim(self.x_array[0],self.x_array[-1])
            self.ax_sequence.set_ylim(-1.5,1.5)
            self.sequence_track.drawnucleotides()
    
    def __init__(self, gene_table = None, genome = None, axes_heights = None, fig_width = 4, show_sequence = True,
                       gene_track = gw.DoubleStrandGeneTrack, gene_track_kwargs = {},
                       sequence_track = gw.SequenceTrack, sequence_track_kwargs = {},
                       coordinate_columns = ['strand', 'start', 'end'], gene_name_column = 'plt_tag', vpad=.15):
        # store important variables for object function
        self.gene_table = gene_table
        self.gene_names = gene_table.loc[:,gene_name_column].values
        self.gene_regions = gene_table.loc[:,coordinate_columns].values.astype(int)
        self.top_positive = None # True if target strand is positive strand, define later based on position of plot
        self.genome = genome
        # generate figure and axes based on above parameters
        fig, all_axes = gw.verticalplots(axes_heights, fig_width, vpad = vpad)
        if show_sequence:
            self.axesdata = all_axes[:-2]
            self.ax_sequence = all_axes[-2]
            self.ax_genes = all_axes[-1]
            self.sequence_track = sequence_track(genome_window = self, **sequence_track_kwargs)
        if show_sequence is False:
            self.axesdata = all_axes[:-1]
            self.ax_sequence = None
            self.sequence_track = None
            self.ax_genes = all_axes[-1]
        self.gene_track = gene_track(genome_window = self, **gene_track_kwargs)
        self.figure = fig
        # remove ticks for all but last data axis, remove top and right spines for all data axes
        [(ax.spines['top'].set_visible(False), ax.spines['right'].set_visible(False), plt.setp(ax.get_xticklabels(),visible=False)) for ax in self.axesdata[:-1]]
        (self.axesdata[-1].spines['top'].set_visible(False), self.axesdata[-1].spines['right'].set_visible(False))
        # for gene axis and sequence axis remove everything
        (self.ax_genes.set_frame_on(False), self.ax_genes.axes.get_xaxis().set_visible(False), self.ax_genes.axes.get_yaxis().set_visible(False))
        if self.ax_sequence is not None:
            (self.ax_sequence.set_frame_on(False), self.ax_sequence.axes.get_xaxis().set_visible(False), self.ax_sequence.axes.get_yaxis().set_visible(False))
        # pad labels
        gw.padxlabels(fig, all_axes)