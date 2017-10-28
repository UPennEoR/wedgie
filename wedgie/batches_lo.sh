#!/bin/bash
#$ -N HolyDaysLow
#$ -l h_vmem=45G
#$ -cwd
#$ -m bea
#$ -V

source activate hera
python2.7 getWedge.py -F /data4/paper/plaplant/HERA2015/upenn_projects/2457548/*.[3-6]*RO -X=81,22,43 -R=200_300 -t
python2.7 getWedge.py -F /data4/paper/plaplant/HERA2015/upenn_projects/2457549/*.[3-6]*RO -X=81,22,43 -R=200_300 -t
python2.7 getWedge.py -F /data4/paper/plaplant/HERA2015/upenn_projects/2457550/*.[3-6]*RO -X=81,22,43 -R=200_300 -t
python2.7 getWedge.py -F /data4/paper/plaplant/HERA2015/upenn_projects/2457551/*.[3-6]*RO -X=81,22,43 -R=200_300 -t
python2.7 getWedge.py -F /data4/paper/plaplant/HERA2015/upenn_projects/2457552/*.[3-6]*RO -X=81,22,43 -R=200_300 -t
python2.7 getWedge.py -F /data4/paper/plaplant/HERA2015/upenn_projects/2457553/*.[3-6]*RO -X=81,22,43 -R=200_300 -t
python2.7 getWedge.py -F /data4/paper/plaplant/HERA2015/upenn_projects/2457554/*.[3-6]*RO -X=81,22,43 -R=200_300 -t
python2.7 getWedge.py -F /data4/paper/plaplant/HERA2015/upenn_projects/2457555/*.[3-6]*RO -X=81,22,43 -R=200_300 -t