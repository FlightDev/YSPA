import os
import glob
from astropy.io import fits

#/home/student/Desktop/Images/iTelescope/20180716-California-T24-GOOD

#   Yo Neal. When you use this program, you have to change a few things between iTelescope and LFOP
#   FIRST, remember to change the file path or you'll be a dummy. Also for LFOP -13 and -12 while
#   for iTelescope it should be -9 and -8. Hopefully you know what to do with those numbers...

dir = '20180714-LFOP-OOFOKAY'
path = '/home/student/Desktop/Images/LFOP/' + dir + '/'
dict = {}
date = ""
for filename in os.listdir(path):
    if filename.endswith(".fit"):
        file = path + str(filename)
        image = fits.open(file)
        s = image[0].header.get("DATE-OBS")
        date = s[:len(s) - 13]
        dict.update({s[len(s) - 12:]: filename})
for key, value in sorted(dict.items()):
    print value + "\t\t" + str(key)
print date
