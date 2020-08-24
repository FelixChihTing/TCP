# Two-color Pyrometry (TCP)
## What is?
TCP is a temperature imaging, or planar meausrement, method that applies the principles of Planck's law.
By monitoring the radiation intensity from an object surface at two different wavelengths, the temperature can be derived.
## Only use a camera and some simple optics
The photo sensor employed in this experiment is a monochromatic CCD camera (PCO pixelfly USB) and the image field was divided by a image doubler aligned on the light-of-sight.
The two image fields are then filtered by two different filters centered at 690 nm and 750 nm, respectively, before entering the camera lens, and ultimately, the photodetector in the camera. The calibration of the current TCP was implemeneted by the use of a blackbody stove with a circular cavity.
The temperature of the blackbody stove was set between 840 and 1000 Celsius degrees with an interval of 10 degrees. The radiation intensity at different temperatures were recorded by the camera and converted to gray scales.
With this calibration data and considering the emissivity of the target of interest, say a metal flame, the true temperature of the flame can be estimated.
## Goal
The main goal of this repository is to convert an tested MATLAB code to a python version.

