import genomewindow as gw
import numpy as np
from Bio.Seq import complement

class SequenceTrack():
    def drawnucleotides(self):
        self.calculatefontsize()
        if self.display_sequence:
            sequence = self.genome_window.genome.seq[self.genome_window.window_left:self.genome_window.window_right+1]
            if not self.genome_window.top_positive:
                sequence = sequence.reverse_complement()
            for x,nt,nt_comp in zip(self.genome_window.x_array, sequence, complement(sequence)):
                text = self.ax.text(x,0,nt, ha='center', va = 'bottom', fontsize = self.fontsize - self.font_spacer, color = self.color_dictionary[nt], **self.font_kwargs)
                # gw.drawpolygon(self.ax, gw.plottedobjectbbox(text,self.ax), color=self.color_dictionary[nt])
                # text = self.ax.text(x,0,nt_comp, ha='center', va = 'top', fontsize = self.fontsize - self.font_spacer, **self.font_kwargs)
                # gw.drawpolygon(self.ax, gw.plottedobjectbbox(text,self.ax), color=self.color_dictionary[nt_comp])
                
    def calculatefontsize(self):
        # calculate maximum fontsize required for single nucleotide spacing
        try:
            self.font_kwargs.pop('fontsize')
        except KeyError: pass
        calculated_fontsize = np.max(gw.calculatenucleotidefontsize(self.ax, **self.font_kwargs))
        self.fontsize = min(calculated_fontsize, self.max_fontsize)
        # decide if sequence will be dispalyed or not
        self.display_sequence = True
        if self.fontsize < self.min_fontsize:
            self.display_sequence = False

    def __init__(self, genome_window = None, max_fontsize = 12, min_fontsize = 4, font_spacer = 1,
                 color_dictionary = {'A':'g','T':'r','G':'k','C':'b'}, **kwargs):
        self.genome_window = genome_window
        self.ax = self.genome_window.ax_sequence
        self.max_fontsize = max_fontsize
        self.min_fontsize = min_fontsize
        self.font_spacer = font_spacer
        self.font_kwargs = kwargs
        self.color_dictionary = color_dictionary