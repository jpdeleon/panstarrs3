#!/usr/bin/env python

"""
Reference:
https://ps1images.stsci.edu/ps1image.html
"""
import requests
import numpy as np
import matplotlib.pyplot as pl
from astropy import units as u
from astropy.table import Table
from PIL import Image
from io import BytesIO
from astropy.io import fits


__all__ = ['Panstarrs']

PIXEL_SCALE = 0.25*u.arcsec/u.pix

class Panstarrs:
    def __init__(self, ra, dec, fov=30*u.arcsec, filters="grizy", color=False,
                 format="fits", output_fov=None):
        """
        Attributes
        ----------
        ra, dec = position in degrees
        fov = field of view
        output_size = output (display) image size in pixels (default = size).
                      output_size has no effect for fits format images.
        filters = string with filters to include
        format = data format (options are "jpg", "png" or "fits")
        """
        self.ra = ra
        self.dec = dec
        self.fov = fov
        self.size = int((self.fov/PIXEL_SCALE).value)
        self.color = color
        self.filters= filters
        self.format = format
        self.output_fov = output_fov
        if self.output_fov is not None:
            self.output_size = self.output_fov/PIXEL_SCALE
        else:
            self.output_size = self.size
        msg = 'The PanSTARRS 3PI survey images cover the sky north of dec-30 deg.'
        if self.dec<-30:
            raise ValueError(msg)

    def get_images(self):
        """Query ps1filenames.py service to get a list of images
        Returns a table with the results
        """

        service = "https://ps1images.stsci.edu/cgi-bin/ps1filenames.py"
        url = (f"{service}?ra={self.ra}&dec={self.dec}&size={self.size}"
               f"&format={self.format}&filters={self.filters}")
        table = Table.read(url, format='ascii')
        return table

    def get_url(self):
        """Get URL for images in the table
        Returns a string with the URL
        """
        if self.color and self.format == "fits":
            raise ValueError("color images are available only for jpg or png formats")
        if self.format not in ("jpg","png","fits"):
            raise ValueError("format must be one of jpg, png, fits")
        table = self.get_images()
        url = ("https://ps1images.stsci.edu/cgi-bin/fitscut.cgi?"
               f"ra={self.ra}&dec={self.dec}&size={self.size}&format={self.format}")
        if self.output_size:
            url = url + f"&output_size={self.output_size}"
        # sort filters from red to blue
        flist = ["yzirg".find(x) for x in table['filter']]
        table = table[np.argsort(flist)]
        if self.color:
            if len(table) > 3:
                # pick 3 filters
                table = table[[0,len(table)//2,len(table)-1]]
            for i, param in enumerate(["red","green","blue"]):
                url = url + f"&{param}={table['filename'][i]}"
        else:
            urlbase = url + "&red="
            url = []
            for filename in table['filename']:
                url.append(urlbase+filename)
        assert len(url)>0, 'empty url'
        return url


    def get_color_img(self):
        """Get color image at a sky position
        Returns the image
        """
        if self.format not in ("jpg","png"):
            raise ValueError("format must be jpg or png")
        url = self.get_url()
        if not self.color and isinstance(url,list):
            r = requests.get(url[0])
        else:
            r = requests.get(url)
        im = Image.open(BytesIO(r.content))
        return im


    def get_gray_img(self, filter="g"):
        """Get grayscale image at a sky position
        Returns the image
        """
        self.color = False
        if self.format not in ("jpg","png"):
            raise ValueError("format must be jpg or png")
        if filter not in list("grizy"):
            raise ValueError("filter must be one of grizy")
        url = self.get_url()
        r = requests.get(url[0])
        im = Image.open(BytesIO(r.content))
        return im

    def get_fits(self, filter=None, verbose=True):
        """Get fits image at a sky position
        Returns the image and header
        """
        self.format = "fits"
        self.filters = filter if filter is not None else "g"
        fitsurl = self.get_url()
        hl = fits.open(fitsurl[0])
        if verbose:
            print(hl.info())
        img, hdr = hl[0].data, hl[0].header
        return img, hdr

if __name__=='__main__':
    from astropy.coordinates import SkyCoord
    from astropy.visualization import ZScaleInterval
    from astropy.wcs import WCS

    tcoord = SkyCoord.from_name('TRAPPIST-1')
    ps = Panstarrs(ra=tcoord.ra.deg, dec=tcoord.dec.deg,
                   filters="g", format="fits", color=False)
    interval = ZScaleInterval(contrast=0.5)
    img, hdr = ps.get_fits(verbose=False)
    zmin, zmax = interval.get_limits(img)

    fig = pl.figure(figsize=(5,5))
    ax = fig.add_subplot(111, projection=WCS(hdr))
    ax.imshow(img, vmin=zmin, vmax=zmax, cmap='viridis', origin='lower')
    ax.grid(color='white', ls='solid')
    pl.show()
