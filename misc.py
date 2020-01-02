import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Polygon

def calculatenucleotidefontsize(axis, **kwargs):
    nt_widths = []
    for nt in 'ATGC':
        # for given font, figure out the maximum size allowing single character per position plotting with no overlap
        bbox10 = stringbboxextents(nt, axis, fontsize = 10)
        bbox9 = stringbboxextents(nt, axis, fontsize = 9)
        width10 = bbox10.x1 - bbox10.x0
        width9 = bbox9.x1 - bbox9.x0
        slope = width10 - width9
        # now determine what fontsize woudld achieve single nucleotide width
        nt_widths.append(10 - (width10 - 1)/slope)
    return nt_widths

def stringbboxextents(string, axis, fontname = plt.rcParams['font.family'][0], fontsize = 12):
    # return data space bounding box of a string in the given font family and size before removing generated object
    text = axis.text(0, 0, string, fontname = fontname, fontsize = fontsize)
    bbox = text.get_window_extent(renderer = axis.get_figure().canvas.get_renderer()).inverse_transformed(axis.transData)
    text.remove()
    return bbox

def plottedobjectbbox(plotted_object, axis):
    # return data space bounding box of a string in the given font family and size before removing generated object
    bbox = plotted_object.get_window_extent(renderer = axis.get_figure().canvas.get_renderer()).inverse_transformed(axis.transData)
    return bbox

def pointtoinchesaxis(ax, x_ax, y_ax):
    xdis, ydis = ax.transData.transform_point((x_ax, y_ax))
    xinch, yinch = ax.get_figure().dpi_scale_trans.inverted().transform_point((xdis, ydis))
    return xinch, yinch

def drawpolygon(axis, bbox, **kwargs):
    if isinstance(bbox,list):
        polygon_coordinates = [[min([b.x0 for b in bbox]),min([b.y0 for b in bbox])],
                               [max([b.x1 for b in bbox]),min([b.y0 for b in bbox])],
                               [max([b.x1 for b in bbox]),max([b.y1 for b in bbox])],
                               [min([b.x0 for b in bbox]),max([b.y1 for b in bbox])]]
    else:
        polygon_coordinates = [[bbox.x0,bbox.y0],
                               [bbox.x1,bbox.y0],
                               [bbox.x1,bbox.y1],
                               [bbox.x0,bbox.y1]]
    axis.add_patch(Polygon(polygon_coordinates, **kwargs))

def ddatatodinchaxis(ax):
    # starting x and y
    x = ax.get_xlim()[0]
    y = ax.get_ylim()[0]
    dinchperx = pointtoinchesaxis(ax, x + 1, y)[0] - pointtoinchesaxis(ax, x, y)[0]
    dinchpery = pointtoinchesaxis(ax, x, y + 1)[1] - pointtoinchesaxis(ax, x, y)[1]
    return dinchperx, dinchpery

def verticalplots(heights, width, vpad=.15):
    fig, axes = plt.subplots(len(heights),1,figsize=(width,1))
    prev_y1 = 0 # store previous y1 for determining location
    for h,ax in zip(heights[::-1],axes): # step through axes vertically and reassign veritcal coordinates
        coor = ax.get_position()
        ax.set_position([coor.x0, prev_y1, coor.x1-coor.x0, h])
        prev_y1 = ax.get_position().y1 + vpad
    return fig, axes[::-1]

def moveaxis(ax, dx=0, dy=0):
    # where x and y are in inches on final figure
    fig = ax.get_figure()
    ax_coor = ax.get_position()
    dx = dx / fig.get_figwidth() # convert to figure coordinates
    dy = dy / fig.get_figheight() #  .5 = 2 1 # .25  = 2 .5
    ax.set_position([ax_coor.x0+dx, ax_coor.y0+dy, ax_coor.x1-ax_coor.x0, ax_coor.y1-ax_coor.y0])

def padxlabels(figure, axes):
    for i in range(len(axes)-1): # assumes axes are in vertical order (top -> bottom)
        label = axes[i].get_xmajorticklabels()[0]
        if label.get_visible(): # is visible, move down next axis to accomdate full label
            label_bottom = label.get_window_extent().inverse_transformed(figure.transFigure).y0
            ax_bottom = axes[i].get_position().y0
            # now move all axes below to accomdate x-pad
            delta_y = ax_bottom - label_bottom
            for j in range(i+1,len(axes)):
                moveaxis(axes[j], dy=-delta_y)

def splitkwargs(kwarg_dict):
    # split into static and dynamic kwargs based on if each is a list or not
    static_kwargs = {}
    dynamic_kwargs = {}
    for key, term in kwarg_dict.items():
        if isinstance(term, list): # check if the term is a list
            dynamic_kwargs[key] = term
        else:
            static_kwargs[key] = term
    # now generate a kwarg dictionary for each set of dynamic kwargs
    split_kwargs = []
    for i in range(len(list(dynamic_kwargs.items())[0][1])):
        single_set = static_kwargs.copy()
        for key, term in dynamic_kwargs.items():
            single_set[key] = term[i]
        split_kwargs.append(single_set)
    return split_kwargs