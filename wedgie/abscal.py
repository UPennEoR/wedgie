# From CASA Memo Calibrating HERA-19 with a Galactic Centre Observation
import os, shutil
import clplot.clquants as clquants

"""
# Add observatory "HERA" at same position as PAPER_SA.  This changes the
# CASA installation -- do this only once per installation.
obstablename = os.getenv("CASAPATH").split()[0] + "/data/geodetic/Observatories/"
tb.open(obstablename, nomodify=False)
paperi = (tb.getcol("Name") == "PAPER_SA").nonzero()[0]
tb.copyrows(obstablename, startrowin=paperi, startrowout=-1, nrow=1)
tb.putcell("Name", tb.nrows()-1, "HERA")
tb.close()
"""

indata = [ os.path.join("/fast/temp/heratest/", x) for x in ["zen.2457545.48011.xx.HH.uvcU.uvfits", "zen.2457545.48707.xx.HH.uvcU.uvfits"]]
procdir = "/fast/temp/heratest/casareduction"

def calname(m, c):
    return os.path.basename(m) + c + ".cal"

def cleanspace():
    shutil.rmtree(procdir, ignore_errors=True)
    os.mkdir(procdir)

def impdata(din):
    msname = os.path.splitext(os.path.basename(din))[0] + ".ms"
    importuvfits(din, msname)
    clquants.reorrder(msname)
    return msname

def fringerot(din):
    """Fringe Rotate to Galactic Centre"""
    fixvis(din, din, phasecenter="J2000 17h45m40s -29d00m28s")
    return din

def gcflagdata(msin):
    flagdata(msin, flagbackup=T, mode="manual", antenna="82")

    flagdata(msin,flagbackup=T,mode=’manual’,spw="0:0~65")
    flagdata(msin,flagbackup=T,mode=’manual’,spw="0:377~387")
    flagdata(msin,flagbackup=T,mode=’manual’,spw="0:850~854")
    flagdata(msin,flagbackup=T,mode=’manual’,spw="0:930~1024")
    flagdata(msin,flagbackup=T,mode=’manual’,spw="0:831")
    flagdata(msin,flagbackup=T,mode=’manual’,spw="0:769")
    flagdata(msin,flagbackup=T,mode=’manual’,spw="0:511")
    flagdata(msin,flagbackup=T,mode=’manual’,spw="0:913")
    flagdata(msin,autocorr=True)
    return msin

def mkinitmodel():
    """Initial model: just a point source in GC direction"""
    cl.addcomponent(flux=1.0,
    fluxunit="Jy",
    shape="point",
    dir="J2000 17h45m40.0409s -29d0m28.118s")
    cl.rename("GC.cl")
    cl.close()

def calinitial(msin):
    """Initial calibration of the dataset"""
    #Fill the model column
    ft(msin, complist="GC.cl", usescratch=True)
    kc = calname(msin, "K")
    gc = calname(msin, "G")
    gaincal(msin,
        caltable = kc,
        gaintype="K",
        solint="inf",
        refant="10",
        minsnr=1,
        spw="0:100~130,0:400~600")
    gaincal(msin,
        caltable=gc,
        gaintype="G",
        solint="inf",
        refant="10",
        minsnr=2,
        calmode="ap",
        gaintable=kc)
    applycal(msin, gaintable=[kc,gc])
    return[kc, gc]

def dosplit(msin, inf, datacolumn="corrected", spw=""):
    """Split the initial ia calibrated data"""
    newms = os.path.basename(msin) + inf + ".ms"
    split(msin, newms, datacolumn=datacolumn, spw=spw)
    return newms

def cleaninit(msin):
    imgname = os.path.basename(msin) + ".init.img"
    clean(vis=msin,
    imagename=imgname,
    niter=500,
    weighting="briggs",
    robust=0,
    imsize=[512,512],
    cell=["500arcsec"],
    mode="mfs",
    nterms=1,
    spw="0:150~900",
    mask="circle[[17h45m00.0s,-29d00m00.00s],32000arcsec]")

def cleanfinal(msl):
    imgname = "GC.combined.img"
    clean(vis=msl,
    imagename=imgname,
    spw="0:60~745",
    niter=5000,
    weighting="briggs",
    robust=-1,
    imsize=[512,512],
    cell=["250arcsec"],
    mode="mfs",
    nterms=1,
    mask="circle[[17h45m00.0s,-29d00m00.00s],32000arcsec]")
    imggal="GC.combined.galcord"
    imregrid(imagename=imgname + ".image", output=imggal, template="GALACTIC")

def dobandpass(msin):
    bc = calname(msin, "B")
    bandpass(vis=msin, spw="", minsnr=1, solnorm=F, bandtype="B", caltable=bc)
    applycal(msin, gaintable=[bc])
    return bc

def setAntennaD(ms, D):
    """Set Antenna Dish Diameter in Antenna Table"""
    tb.open(os.path.join(ms, "ANTENNA"), nomodify=False)
    d = tb.getcol("DISH_DIAMETER")
    d[:] = D
    tb.putcol("DISH_DIAMETER", d)
    tb.flush()
    tb.close()

def main():
    cleanspace()
    os.chdir(procdir)
    ms = [impdata(x) for x in indata]
    ms=[gcflagdata(x) for x in ms]
    if 0:
        #Don’t fringe rotate for mosaic later
        ms = [fringerot(x) for x in ms]
    mkinitmodel()
    cals1 = [calinitial(x) for x in ms]
    msical = [dosplit(x, "ical") for x in ms]
    imgs1 = [cleaninit(x) for x in msical]
    cals2 = [dobandpass(x) for x in msical]
    ms2 = [dosplit(x, "c2", spw="0:100~880") for x in msical]
    imgs2 = [cleaninit(x) for x in ms2]
    cals3 = [dobandpass(x) for x in ms2]
    cleanfinal(ms2)

if 1:
    main()