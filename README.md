# RAPID-ZTF

# FastFinder

To begin using RAPID-ZTF:

1. Follow the [ZTF ForcedPhot handbook](http://web.ipac.caltech.edu/staff/fmasci/ztf/forcedphot.pdf) to create an account on the ZTF Forced Photometry Server and request forced photometry on your object.
2. Create a .txt file of the forced photometry as is and add it to to the candidates directory.
3. In a Terminal window, cd into the RAPID-ZTF directory and run startup.py.

Assuming you already have all the Python dependencies installed, it may be ran like the following example:

    python startup.py -n ZTF21aaqytjr -z 0.003319 -ra 186.423654 -dec 7.228395 -mwebv 0.2851 -v
    python startup.py -h          # to get syntax help


The .pdf and .mp4 classifiaction result files for your object will be saved in the outputs directory.


List of Python dependencies RAPID-ZTF needs to run:
```
astrorapid
matplotlib
numpy 
pandas
```
