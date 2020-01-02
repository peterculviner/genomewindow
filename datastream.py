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
        self.axis.plot(self.genome_window.x_array, plot_data, rasterized=True, **kwargs)
        # call DataStream's plotdata function to do final plot clean up, remaining kwargs are passed
        DataStream.plotdata(self, **kwargs)
