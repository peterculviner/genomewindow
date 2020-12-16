import numpy as np

class DataStream():
    # basic class through which all data is displayed
    def plotdata(self, **kwargs):
        # requried things to do after every derived class's data is plotted
        # move x-limits to extent of data
        self.axis.set_xlim(self.genome_window.x_array[0], self.genome_window.x_array[-1])
    
    def __init__(self, input_data = None, axis_n = None, genome_window = None, on_strand=True, rasterize_data=True, **kwargs):
        # store key variables for all types of DataStream
        self.genome_window = genome_window # the plot to which the data stream belongs, this allows pulling all position information etc
        self.axis = self.genome_window.axesdata[axis_n] # pull axis from zero-genome window
        self.on_strand = on_strand
        self.rasterize_data = rasterize_data
        # plot the data, all derived classes have their own version of this function
        self.plotdata(input_data, **kwargs)


class ArrayAsLine(DataStream):
    # plot a numpy array (2 x len(genome)) or len(genome) as a line
    def plotdata(self, input_data, **kwargs):
        # narrow the window_mask (which is stranded) to which strand of data to include (if data is stranded)
        data_mask = self.genome_window.window_mask[int(not self.genome_window.top_positive)]
        if len(input_data.shape) == 1:
            plot_data = input_data[data_mask]
        elif len(input_data.shape) == 2:
            plot_data = input_data[int(np.logical_xor(self.genome_window.top_positive, self.on_strand)),data_mask]
        elif len(input_data.shape) != 1:
            raise(ValueError('Input data shape of %s not supported.'%str(input_data.shape)))
        if not self.genome_window.top_positive: # data needs to be reversed, as genome view has been reversed 
            plot_data = np.flip(plot_data)
        # now plot the line
        self.axis.plot(self.genome_window.x_array, plot_data, rasterized=self.rasterize_data, **kwargs)
        # call DataStream's plotdata function to do final plot clean up, remaining kwargs are passed
        DataStream.plotdata(self, **kwargs)

class MotifArrows(DataStream):
    # plot a numpy array (2 x len(genome)) or len(genome) as a line
    def plotdata(self, input_data, autoscale = False, **kwargs):
        # determine left and right based on lengths: [[strand ->], [positions (left) ->], [lengths ->], [scores (y-values) ->]]
        strands, left_positions, lengths, scores = input_data
        right_positions = left_positions + lengths - 1
        # determine which inputs fall into current window
        plot_mask = np.all([self.genome_window.window_right >= left_positions, self.genome_window.window_left <= right_positions], axis=0)
        # pare down data set to only those we're plotting, convert into plot coordinates
        strands = strands[plot_mask]
        left_positions = self.genome_window._genometoplotposition(left_positions[plot_mask])
        right_positions = self.genome_window._genometoplotposition(right_positions[plot_mask])
        scores = scores[plot_mask]
        if autoscale == True:
            # rescale ylimits based on scores
            ydelta = np.max(scores) - np.min(scores)
            yrange = [np.min(scores)-ydelta*.25,np.max(scores)+ydelta*.25]
            self.axis.set_ylim(yrange)
        # check if kwargs has arrowprops
        try:
            arrowprops = kwargs.pop('arrowprops')
        except KeyError:
            arrowprops = {'width':3,'headwidth':8,'headlength':8,'facecolor':'k','edgecolor':None}
        # now plot motifs as arrows
        for st,l,r,sc in zip(strands, left_positions, right_positions, scores):
            # if self.genome_window.top_positive:
            if st == 0: # right is arrow end
                self.axis.annotate('', xy=[r,sc], xytext=[l,sc], arrowprops=arrowprops, **kwargs)
            else: # left is arrow end
                self.axis.annotate('', xy=[l,sc], xytext=[r,sc], arrowprops=arrowprops, **kwargs)
            # else:
            #     if st == 0: # left is arrow end
            #         self.axis.annotate('', xy=[l,sc], xytext=[r,sc], arrowprops=arrowprops, **kwargs)
            #     else: # right is arrow end
            #         self.axis.annotate('', xy=[r,sc], xytext=[l,sc], arrowprops=arrowprops, **kwargs)
        # call DataStream's plotdata function to do final plot clean up, remaining kwargs are passed
        DataStream.plotdata(self, **kwargs)

class LocalMotifHatch(DataStream):
    # plot a numpy array (2 x len(genome)) or len(genome) as a line
    def plotdata(self, input_data, autoscale = False, **kwargs):
        # determine left and right based on lengths: [[strand ->], [positions (left) ->], [lengths ->], [scores (y-values) ->]]
        strands, left_positions, lengths, scores = input_data
        # call DataStream's plotdata function to do final plot clean up, remaining kwargs are passed
        DataStream.plotdata(self, **kwargs)

class ShadeRegions(DataStream):
    # plot a numpy array (2 x len(genome)) or len(genome) as a line
    def plotdata(self, input_data, intype = 'ss', color='k', alpha=.25,  **kwargs):
        if intype is 'ss':
            # identify overlapping regions from the input
            overlap_i = np.where(np.all([self.genome_window.window_left <= input_data[:,1],
                                         self.genome_window.window_right >= input_data[:,0]],axis=0))
            overlapping_regions = input_data[overlap_i]
            # update coordinates from genomic to plot
            left_positions = self.genome_window._genometoplotposition(overlapping_regions[:,0])
            right_positions = self.genome_window._genometoplotposition(overlapping_regions[:,1])
            for l,r in zip(left_positions, right_positions):
                # fill between defined by where function
                self.axis.fill_between(self.genome_window.x_array, self.axis.get_ylim()[0],self.axis.get_ylim()[1],
                                                     where=np.all([self.genome_window.x_array >= np.min([l,r]), 
                                                                   self.genome_window.x_array <= np.max([l,r])], axis=0),
                                                     linewidth=0, color=color, alpha=alpha, **kwargs)
        else:
            raise(ValueError('only ss regions supported right now (shape: (n_regions, 2))'))
        # call DataStream's plotdata function to do final plot clean up, remaining kwargs are passed
        DataStream.plotdata(self, **kwargs)