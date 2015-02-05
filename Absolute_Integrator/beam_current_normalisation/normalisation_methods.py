The normalization of ADF data requires two values to be known, the D.C. offset and the effective detector sensitivity. The effective detector sensitivity is the combined effect of the intrinsic sensitivity of the detector and the amplifier gain.

user options should be defined in analyse_CCD_flux.py


currentNormalization = faradayCurrentImaging / faradayCurrentDetector  # Increase in current used for imaging.




def file_opener():
	#a function with a bunch of if statements on what should be done to the detector image

#to look up... change the following below to a detector class...
#------------------------------------------------------------------------


def detectorThresholded(detectorImage):
	#function that measures detector threshold

def detectorCentreFinder(detectorImage):
	#function that finds the centre of the detector

def detectorSensitivity():
	#procedures for mean, radial and azimuthal sensitivity

def fluenceCone():
	#calculate the fluence cone


#------------------------------------------------------------------------
def noiseTargetSNR():

def filesaver():
	#saves the important/selected datafiles of the images from above (like matlab workspace)

def plotter():
	#plots the desired figures that the user wants


#main body that pulls the functions from above depending on user choices goes here   

	#includes useroptions linked to current normalisation, spectrum integrator
	#analyse CCD flux, or neither
	#detector threshold

















#add file opener / message prompt for loading the detectorImage 
