import os

days = ['2457548', '2457549', '2457550', '2457551', '2457552', '2457553', '2457554', '2457555']
rng = '580_680'

for day in days:
    os.system("python2.7 getWedge.py -F /data4/paper/plaplant/HERA2015/upenn_projects/{day}/*RO -X=81,22,43 -R={rng} -t".format(
        day=day, rng=rng))
    folder = "{day}_{rng}".format(day=day, rng=rng)
    os.mkdir(folder)
    os.system("mv *{day}* {folder}".format(day=day, rng=rng))
