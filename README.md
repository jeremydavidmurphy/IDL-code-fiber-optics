# IDL-code-fiber-optics

This directory includes several IDL routines I wrote to characterize various properties of the fiber optics used on the Mitchell Spectrograph (formally VIRUS-p) in use on the Harlan J. Smith telescope at McDonald Observatory.

There are a number of routines here, many of which are simply one-off pieces of code to perform specific tasks, often some type of visual comparison of data.

The two primary properties of the fiber optics I was testing transmission and what's called 'focal ratio degradation' or FRD, both as a function of wavelength. Therefore, you'll see a lot of mention of 'trans' and 'FRD' in the names of the routines and the comments.

FRD is the affect where light sent into a fiber at a given focal ratio has a tendency to have a higher focal ratio when it leaves the fiber. This 'degradation' is important for astronomy as fiber optics are commonly used to feed instruments (often spectrographs) and characterizing how your fibers increase the focal ratio is critical for designing the instrument's optics.

The two primary pieces of code are named trans_FRD.pro and fratioF.pro, so if you are poking around, I suggest you start there. A good followup is lifetime2.pro which reduced data from a series of flexure tests we ran on the fibers than simulated 10 years of life on the telescope. Cool stuff!

As this code was never intended to be public, much of the parameters are placed directly into the code rather than being passed. The routines accept in a list of .fits files (a standard format for CCD readouts in astronomy), then iterate through those lists to perform a variety of tasks.

