import sys
import numpy as np
import matplotlib.pyplot as plt
import aplpy
from astropy.io import fits
from astropy import wcs


def cutout(inputfile, xcenter, ycenter, width, height, outfile, cdelt = ''):

    ### This function creates a new fits file based on rectangular region information.
    ### xcenter, ycenter are the centers of the region in galactic longitude and latitude in degrees
    ### width and height are similarly longitude and latitude measurements in degrees
    
    w = wcs.WCS(inputfile).celestial
    print(w)
    data, header = fits.getdata(inputfile, header = True)
    print(data.shape)
    data = data[0,:,:]
    print(data.shape)
    #if (header['CRVAL1'] != 0) or (header['CRVAL2'] != 0):
    #    raise ValueError('\nThe input CRVAL is not 0!\n')
    
    center_x, center_y = w.wcs_world2pix(xcenter, ycenter, 1)
    print(center_x,center_y)
    center_x = int(center_x)
    center_y = int(center_y)
    
    if cdelt == '':
        cdelt = np.mean( [np.absolute(header['CDELT1']), np.absolute(header['CDELT2'])] )
        
    size_pixels_x = int(np.rint((width / 2.) / cdelt))
    size_pixels_y = int(np.rint((height / 2.) / cdelt))
    
    delta_y = np.array([center_x - size_pixels_x, center_x + size_pixels_x], dtype = int)
    delta_x = np.array([center_y - size_pixels_y, center_y + size_pixels_y], dtype = int)
    
    center_x_new, center_y_new  = w.wcs_pix2world( delta_y[0], delta_x[0], 1)
    
    data_cutout = data[delta_x[0]:delta_x[1], delta_y[0]:delta_y[1]]

    center_x = ((delta_x[1] - delta_x[0]) / 2.) + delta_x[0]
    center_y = ((delta_y[1] - delta_y[0]) / 2.) + delta_y[0]
    
    del header[33:] ### strictly for FORCAST files, needs to be 28 for cmzoom mosaic I think...
    del header['PC*']
        
    header['CRPIX1'] = size_pixels_x
    header['CRPIX2'] = size_pixels_y
    header['CRVAL1'] = xcenter
    header['CRVAL2'] = ycenter
    
    fits.writeto(outfile, data = data_cutout, header = header, overwrite = True)
    
def make_single(inputfits, outputfile, xcenter, ycenter, label, vmin=1.0e06, vmax =3.0e+08, contourfile=None,
                sizex=4,sizey=4):
    
    fg_color='white'
    bg_color='black'
    fig = plt.figure(figsize=(sizex, sizey))#,facecolor=bg_color, edgecolor=fg_color)
    subplot = aplpy.FITSFigure(inputfits, figure = fig, convention='calabretta')
    subplot.show_colorscale(vmin=vmin, vmax=vmax, stretch='log',cmap='inferno')
    #subplot.frame.set_color(fg_color)
 
    subplot.set_nan_color(fg_color)

    #subplot.ticks.set_xspacing(0.02)
    #subplot.ticks.set_yspacing(0.02)
    subplot.ticks.set_color('black')
    subplot.tick_labels.set_xformat('d.dd')
    subplot.tick_labels.set_yformat('d.dd')
    #subplot.recenter(xcenter, ycenter, width = 50. / 3600., height = 50. / 3600.)
    subplot.add_label(0.35, 0.95, str(label), relative = True, weight = 'bold', size = 10, color = 'black')
    if contourfile!=None:
        subplot.show_contour(contourfile,
                             colors = 'white', levels = [0.0,1.0], linewidths = 0.7,
                             convention = 'calabretta', zorder = 10, linestyle = 'solid')
    subplot.axis_labels.set_xtext('GLON')
    subplot.axis_labels.set_ytext('GLAT')
    subplot.ticks.show()
    #subplot.add_scalebar(length=24./3600.)
    #subplot.scalebar.set_label('1 pc')
    #subplot.scalebar.set_color('k')
    fig.tight_layout()
    fig.savefig(outputfile,bbox_inches='tight',dpi=300)
    plt.show()

        
