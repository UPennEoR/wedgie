import imageio
import glob

images = sorted(glob.glob("./*png"))
gif = []

for filename in images:
    gif.append(imageio.imread(filename))
 
kargs = { 'duration': 1.5 } # in seconds
imageio.mimsave("./all_day.gif", gif,**kargs)