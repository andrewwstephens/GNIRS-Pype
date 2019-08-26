#!/usr/bin/env python

################################################################################
# CHANGE LOG                                                                   #
#                                                                              #
# 23 Jul 2014 - * Moved the nc term to the end of the SNR equation to increase #
#               readability.                                                   #
#		* Removed the 1/2 factor multiplying the nc term, because the  #
#               number of images in NSCHL0* cards is already the amount of     #
#               source exposures.                                              #
#               * Removed the import statement for spectools. Since the        #
#               the functions in this script no longer write ASCII files,      #
#               there is no need to build an array of wavelengths.             #
#                                                                              #
# 16 Apr 2015 - * Added the joinorders_noresampling function to merge the      #
#               spectra without resampling to a linear dispersion across       #
#               the different orders of the cross dispersed exposure.          #
#                                                                              #
# KNOWN ISSUES                                                                 #
#                                                                              #
# - Limit to 99 combined images: The nc term is the result of a search in the  #
# header for NSCHL0* cards, which end in a two digit counter, e.g. NSCHL001.   # 
# Therefore, there is an artificial limitation to 99 combined source images    #
# until a better way of determining the number of combined exposures is        #
# implemented. Nevertheless, this should no be a problem since the number of   #
# combined exposures is usually around 10.                                     #
#                                                                              #
################################################################################



from numpy import *
from astropy.io import fits as pf
from copy import deepcopy
from pyraf import iraf
from scipy.interpolate import interp1d

def get_wl(image, dimension=0, hdrext=0, dataext=0, dwlkey='CD1_1', wl0key='CRVAL1', pix0key='CRPIX1'):
  
  """
  Obtains the wavelenght coordinates from the header keywords of the
  FITS file image. The default keywords are CD1_1 for the delta lambda,
  CRVAL for the value of the first pixel and CRPIX1 for the number of
  the first pixel. These keywords are the standard for GEMINI images.
  
  The function is prepared to work with Multi-Extesion FITS (MEF) files.

  Parameters:
  -----------
  image : string
    Name of the FITS file containing the spectrum
  dimension : integer
    Dimension of the dispersion direction
  hdrext : number
    Extension that contains the header
  dataext : number
    Extension that contains the actual spectral data
  dwlkey : string
    Header keyword for the interval between two data points
  wl0key : string
    Header keyword for the first pixel value
  pix0key : string
    Header keyword for the first pixel coordinate
  
  Returns:
  --------
  wl : numpy.array
    Wavelength coordinates for each data point
  """
 
  h = pf.getheader(image,ext=hdrext)
  dwl,wl1,pix0 = [float(h[i]) for i in [dwlkey,wl0key,pix0key]]
  npoints = shape(pf.getdata(image,dataext))[dimension]
  wl = wl1 + (arange(1,npoints+1)-pix0)*dwl 
    
  return wl

def joinorders_noresampling(inlist,merged_spec,outspec,writefits,snrlist=None):
  """
  Join the different spectral orders from the cross-dispersed exposure
  into a single spectrum, having the wavelength coordinates as a
  lookup table in the header.

  Parameters:
  -----------
  inlist : list
    List of the names of FITS files containing the 1D spectra of each order.
    In the standard nomenclature of the XDGNIRS pipeline these would be the
    final?.fits files.
  merged_spec : string
    Name of FITS file containing the final merged spectrum, combined with 
    ODCOMBINE. This is required to get the wavelength limits of each order.
  outspec : string
    Name of the final output spectrum. The code also writes an ASCII file with
    the same root name, but with '.fits' replaced with '.txt'. Therefore
    outspec needs to have '.fits' at the end.
  writefits : boolean
    Writes the resulting stitched spectrum in a FITS file.
  snrlist : list or string
    If snrlist is a list, it is expected to have six elements, each one
    being the name of a FITS 1D signal-to-noise ratio spectrum. The order of
    this list has equal to the inlist, from longer to shorter wavelengths.
    snrlist can also be the name of a MEF FITS file, where the signal-to-noise
    ratio is written in the 3n+2 extension, with n from 1 to 6.


  Returns:
  --------
  jspec : numpy.ndarray  
    jspec[:,0] - Wavelength coordinates
    jspec[:,1] - Flux coordinates

  Description:
  ------------
  This function was developed to produce a stitched spectrum spanning the full
  wavelength range available in the cross dispersed exposure of GNIRS without
  resampling the different orders with the same wavelength step. Instead each
  order carries the dispersion it had when it was flux calibrated, and the
  final stitched spectrum has the wavelength coordinates written in a lookup
  table at the image header, with 6 differente dispersions.

    -- Motivation --
  The ODCOMBINE routine that is employed by the pipeline can only produce a
  spectrum with a linear dispersion function, defined by a fixed step and
  starting coordinate. This of course leads to oversampling of the regions
  with lower dispersion, or longer wavelengths. Originally the oversampling
  was done with the highest dispersion between the input spectra, which led
  to a sampling frequency of roughly 3 at the 6th order, and still smaller
  frequencies for every other order. Such a low frequency caused the flux
  values to be 'smoothed', specially when using high order polynomials for the
  interpolation.

    -- Methodology --
  The main problem with stitching the different difraction orders was the
  appropriate treatment of the overlaping regions. The method employed here
  consist in storing each spectrum as a scipy.interpolate.interp1d function
  and averaging both spectra at the wavelength coordinates of the one with the
  higher dispersion.
  """
  
  if type(snrlist) == list:
    snrmef = False
  elif type(snrlist) == str:
    if len(pf.open(snrlist)) == 19:
      snrmef = True
    else:
      print 'ERROR! Expecting 19 extension MEF as snrlist.'
      return 

  #
  #  Reading general info from the input files
  #
  for i in inlist:
    print i
  cdelts = array([pf.getheader(i)['CD1_1'] for i in inlist])
  # cdelts is not actually needed, but who knows?
  wls = array([get_wl(i) for i in inlist])

  #
  # Storing each spectrum as a linear interpolation function
  #
  f = array([interp1d(wls[i],pf.getdata(inlist[i])) for i in range(6)])
  if snrlist != None:
    if snrmef:
      n = array([interp1d(wls[i],pf.getdata(snrlist,ext=(3*i+2))) for i in range(6)])
    else:
      n = array([interp1d(wls[i],pf.getdata(snrlist[i])) for i in range(6)])
 
  #
  # Getting the limits of each order
  #
  h = pf.getheader(merged_spec)['IMC*']

  # gb stands for stringed good bits
  # If it takes longer than one line it is wrong.

  gb = array([h[i][h[i].find('[')+1:h[i].find(']')].split(':') for i in range(6)],dtype='int')
  wllims = array([[wls[i][gb[i,0]],wls[i][gb[i,1]]] for i in range(6)],dtype='int')

  # js stands for joined spectrum
  js = array([])
  jwl = array([])

  if snrlist != None:
    jsnr = array([])

  for i in arange(6)[::-1]:

    # shortest wavelength
    if i == 5:
      singlecondition = (wls[i] > wllims[i,0])&(wls[i] < wllims[i-1,0])
      overlapcondition = (wls[i] > wllims[i-1,0])&(wls[i] < wllims[i,1])

    # longest wavelength
    elif i == 0:
      # only the single spectrum part
      singlecondition = (wls[i] > wllims[i+1,1])&(wls[i] < wllims[i,1])
      js = append(js,f[i](wls[i][singlecondition]))
      jwl = append(jwl,wls[i][singlecondition])
      if snrlist != None:
        jsnr = append(jsnr,n[i](wls[i][singlecondition]))
      break

    # everything else
    else:
      singlecondition = (wls[i] > wllims[i+1,1])&(wls[i] < wllims[i-1,0])
      overlapcondition = (wls[i] > wllims[i-1,0])&(wls[i] < wllims[i,1])

    js = append(js,f[i](wls[i][singlecondition]))
    js = append(js,average([f[i](wls[i][overlapcondition]),f[i-1](wls[i][overlapcondition])],0))

    jwl = append(jwl,wls[i][singlecondition])
    jwl = append(jwl,wls[i][overlapcondition])
    
    if snrlist != None:
      jsnr = append(jsnr,n[i](wls[i][singlecondition]))
      # this should be revised as soon as the code starts working
      jsnr = append(jsnr,average([n[i](wls[i][overlapcondition])+n[i-1](wls[i][overlapcondition])],0)/sqrt(2.))
   
  #
  #  Writing the 1D spectrum to a FITS file with the aid of noao.onedspec.rspectext
  #
  outtext = outspec.replace('.fits','.txt')

  if snrlist == None:
    jspec = column_stack([jwl,js])
    if writefits:
      savetxt(outtext,jspec,fmt='%.6e\t%.6e')
  else:
    jspec = column_stack([jwl,js,jsnr])
    if writefits:
      savetxt(outtext,jspec,fmt='%.6e\t%.6e\t%.6e')

  if writefits:
    iraf.noao()
    iraf.onedspec()
    iraf.rspectext(outtext,outspec,dtype='nonlinear',flux='yes',title=pf.getheader(merged_spec)['OBJECT'])

  return jspec 

def build_unc(inimage,outimage=None,clobber=False,writefits=False,aperture_width=12,skyglow=True,readnoise=None,dark=0.15,verbose=False,skyspec='/media/storage/Dropbox/palomar/xdgnirs/skyext_persec_perpix.fits'):
  """
  Builds a Poisson noise estimate for the all the
  extensions in the image ''inimage''.

  Parameters:
  -----------
  inimage : string
    Name of the input image.
  outimage : string
    Name of the output image.
  clobber : boolean
    Overwrite existing images.
  writefits : boolean
    Writes the results in a FITS image.
  aperture_width : float
    Default aperture width to use when there is no
    APNUM1 card in the header.
  skyglow : boolean
    Includes the emission lines from the atmosphere
    in the noise estimate.
  readnoise : float
    Reading noise in electrons per pixel. If None
    its value will be read from the 'rdnoise' card in
    the image header.
  dark : float
    Detector's dark current in electrons per pixel
    per second. The standard value comes from the manual
    of GNIRS, and is 0.15.
  skyspec : string
    Atmospheric emission spectrum in units of ADU per second
    per pixel.
  
  Returns:
  --------
  snrlst : list
    A list containing all the SNR spectra, beggining with the
    one with longer wavelength (K).


  Description:
  ------------
  This function takes a FITS file with 18 extensions and
  writes a copy of it with every 2+3n extension replaced
  by a S/N estimate.

  The estimate relies on the 1+3n extensions being the
  spectrum in electrons s^-1 units. According to Steve B. Howell,
  in Handbook of CCD Astronomy (2nd ed. page 73), and also
  Mortara & Fowler 1981 (SPIE Conference Proceedings, Vol. 290.,
  page 28), the signal to noise ratio is given by

  S/N = (Nx)/sqrt(Nx + npix*(Nd*t + Ns*t + Nr**2))

  where Nx is the number of electrons from the astronomical
  source, t is the exposure time, npix is the number of pixels
  in the aperture, Nd is the dark current in electrons per pixel
  per second, Ns is the number of electrons from the
  atmospheric emission per pixel and per second and Nr is the
  read noise times the square root of the number of exposures.

  By comparison with data from the GNIRS Integration Time Calculator
  it has been found that the data at the 'vtarget_comb.fits' stage
  is already in units of electrons, therefore no gain multiplication
  is being performed by this function.
  """
  
  im = pf.open(inimage)
  h =  im[0].header
  nc = len(h['NSCHL0*'])  # number of combined frames
  sky = pf.getdata(skyspec)
  sky[sky < 1.e-3] = 1.e-3  # to avoid divisions by zero
  snrlst = []


  for i in arange(1,17,3):
    hext = im[i].header
    t = float(hext['exptime'])
    s = deepcopy(im[i].data)  
    if readnoise == None:
      readnoise = 7.  # assuming the best possible read noise from GNIRS' manual
#    rn = readnoise*sqrt(nc)  # the read noise is not linear with the number of exposures
    rn = readnoise
    s[s < 1.e-3] = 1.e-3  # to avoid divisions by zero

    # getting the aperture width
    try:
      npix = diff([float(k) for k in hext['apnum1'].split()[2:]])
    except KeyError:
      if verbose:
        print 'WARNING! APNUM1 card not found in {:s}.'.format(inimage+'['+str(i)+']')
        print 'Proceeding with default value of {:.1f} pixels.'.format(aperture_width)
      npix = aperture_width

    if skyglow:
      snr = s/sqrt(s+npix*(sky[(i-1)/3]*t+dark*t+rn**2)) * sqrt(nc)
    else:
      snr = s/sqrt(s+npix*(dark*t+rn**2)) * sqrt(nc)

    im[i+1].data = snr
    snrlst.append(snr)
    if verbose:
      print 'npix = {:.2f}; nc = {:d}; Combined read noise = {:.2f}; Texp = {:.2f}'.format(float(npix),nc,rn,t)

  if writefits:
    if outimage == None:
      outimage = inimage[:inimage.find('.fits')]+'_unc.fits'
    im.writeto(outimage,output_verify='silentfix',clobber=clobber)

  return snrlst

def join_extensions(infile,reference=None,outjoined=None,outname=None):
  """
  Joins in to a single spectrum the different orders in a GNIRS XD 
  spectrum, with one order from each FITS extension. This functions
  is specifically aimed at bringing together the SNR spectrum with
  the flux calibrated one.

  Parameters:
  -----------
  infile : string
    Input FITS file with SNR spectra in the 2+3n extensions.
  reference : string
    Name of the FITS file containing the combined, flux calibrated
    and redshift corrected spectrum.
  outjoined : string
    Name of the "iraf.noao.onedspec.odcombine" output.
  output : string
    Name of the final output file, which consists of an ASCII with
    three columns: wavelength, flux density, error in flux density.

  Returns:
  --------
  Nothing.


  """

  s = infile[:infile.find('_unc.fits')]

  if reference == None:
    try:
      srest = s+'_rest.fits'
    except IOError:
      srest = srest.replace('NGC','ngc')
  else:
    srest = reference

  h = pf.getheader(srest)['IMC*']
  hfull = pf.getheader(srest)

  if pf.getheader(infile)['object'] == pf.getheader(srest)['object']:
    pass
  else:
    print 'ERROR! Objects in input and reference files do not match.'
    print 'Input object: {:s}'.format(pf.getheader(infile)['object'])
    print 'Reference object: {:s}'.format(pf.getheader(srest)['object'])
    print 'Exiting function now.'
    return

  gb = [h[i][h[i].find('[')+1:h[i].find(']')] for i in range(6)]
  
  l = open('listod.lst','w')

  for i in arange(6):
    l.write(infile+'['+str(2+3*i)+']['+gb[i]+']\n')
  l.close()

  iraf.noao()
  iraf.onedspec()
  if outjoined == None:
    outjoined = s+'_unc_joined.fits'
  iraf.odcombine(input='@listod.lst',output=outjoined,headers='',bpmasks='',rejmasks='',nrejmasks='',expmasks='',sigmas='',logfile='STDOUT',apertures='',group='all',first='no',w1='INDEF',w2='INDEF',dw='INDEF',nw='INDEF',log='no',combine='average',reject='none',outtype='real',outlimits='',smaskformat='bpmspectrum',smasktype='none',smaskvalue=0.0,blank=0.0,scale='none',zero='none',weight='none',statsec='',expname='',lthreshold='INDEF',hthreshold='INDEF',nlow=1,nhigh=1,nkeep=1,mclip='yes',lsigma=3.0,hsigma=3.0,rdnoise='0.0',gain='1.0',snoise='0.0',sigscale=0.1,pclip=-0.5,grow=0.0,offsets='physical',masktype='none',maskvalue=0.0,mode='al')

  srestdata = pf.getdata(srest)
  suncdata = pf.getdata(outjoined)
#  suncdata[suncdata < 1.e-3] = 1.e-3

  if len(srestdata) != len(suncdata):
    print '*****************'
    print 'ERROR! Length of flux calibrated spectrum does not match the'
    print 'length of combined uncertainty spectrum.'
    print '{:s} : {:d}'.format(srest,len(srestdata))
    print '{:s} : {:d}'.format(outjoined,len(suncdata))
    print '*****************'
    return

  if outname == None:
    outname = s+'_flam_unc.fits'
#  savetxt(output,column_stack([wl,srestdata,srestdata/suncdata]),fmt='%.6e')
  outh = deepcopy(hfull)
  outh['CTYPE2'] = 'LINEAR'
  outh['CRVAL2'] = 1
  outh['CRPIX2'] = 1
  outh['CDELT2'] = 1
  outh['CD2_2'] = 1
  outh['WCSDIM'] = 2
  outh['LTM2_2'] = 1
  outh['APNUM2'] = '2 2'
  outh['BANDID1'] = 'Spectrum'
  outh['BANDID2'] = 'Flux uncertainty'

  pf.writeto(outname,data=row_stack([srestdata,srestdata/suncdata]),header=outh)

  return
