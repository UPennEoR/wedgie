import imageio
import glob
import argparse
import pprint

parser = argparse.ArgumentParser()
parser.add_argument('-F', '--files', help='Files to be giffed. Must use "-F=" notation.', nargs='*', required=True)
parser.add_argument('-s', '--save', help='Input path for save location.', default='./')
parser.add_argument('-d', '--duration', help='Set duration of each from of gif in seconds.', default=2, type=float)
args = parser.parse_args()


images = sorted(glob.glob("{}".format(" ".join(args.files))))

gif =[]
for image in images:
    gif.append(imageio.imread(image))

day = images[0].split('/')[-1].split('.')[1]
start = images[0].split('/')[-1].split('.')[2]
end = images[-1].split('/')[-1].split('.')[2]
length = len(images)
 
kargs = {'duration': args.duration} # in seconds
imageio.mimsave("{}day{}__{}files__start{}__end{}.gif".format(args.save, day, length, start, end), gif,**kargs)
