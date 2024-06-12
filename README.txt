The material in this folder is taken from Marius folder in the share. I (Martino) was looking for the original GSA code and found something and started working on it. At the moment (2024.06.10) the one that i have used most is SmallPixel, for which I prepare the image using 2gauss-Npx.

Possible GSAs

THE MOST PROMISING
\\share.univie.ac.at\a511-quantum\Work\microscopy\Personal Share\Marius\Dissertation data\Gerchberg_Saxton_Algorithm

THE LATEST VERSION 23.06.2022
\\share.univie.ac.at\a511-quantum\Work\microscopy\Personal Share\Marius\Dissertation data\GSA ligh\Tiger...

OLDER VERSIONS summer.2021
\\share.univie.ac.at\a511-quantum\Work\microscopy\Personal Share\Marius\GSA_small\*.m

,-------------------------------,
|     NOTES on the codes        |
'-------------------------------'

I copied marius codes and give them new names:
TigerInteratively -> GSA1
SmallPixel -> GSA2
smiley -> GSA3

The aim is to merge them and then translate them into python.

Marius codes are in a separate folder and cannot run without the functions that are now in the lib folder.

,-------------------------------,
| Thomas suggestions 2024.06.12 |
'-------------------------------'

Implement our GSA in python (starting from marius code) with following features (* = priority):
Possibility to define:
    - *target image (intensity)
    - *initial phase on SLM
    - *laser beam on SLM (both phase [flat] and intensity [gauss])
    - noise region in the target plane
Possibility to:
    - add aberrations (can be done directly on the SLM -> Tilman))
    - shift the target image in X and Y (idem)
Implement:
    - *GSA with FFT
    - GSA with propagation and BB


