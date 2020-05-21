[![Build Status](https://travis-ci.com/I2PC/xmipp.svg?branch=devel)](https://travis-ci.com/I2PC/xmipp)
<!---  [![Quality Gate](https://sonarcloud.io/api/project_badges/measure?project=Xmipp&metric=alert_status)](https://sonarcloud.io/dashboard?id=Xmipp)
[![Technical debt](https://sonarcloud.io/api/project_badges/measure?project=Xmipp&metric=sqale_index)](https://sonarcloud.io/component_measures?id=Xmipp&metric=sqale_index)
[![Bugs](https://sonarcloud.io/api/project_badges/measure?project=Xmipp&metric=bugs)](https://sonarcloud.io/project/issues?id=Xmipp&resolved=false&types=BUG)
--->
# Xmipp-Lite

Xmipp is a suite of image processing programs, primarily aimed at single-particle 3D electron microscopy. 

This is a light version of Xmipp, which only contains the basic tools addres to estimate resolution measurements.

## Getting started

Start by cloning the repository from GitHub and go there.
```
git clone https://github.com/Vilax/xmipp-lite.git
cd xmipp-lite
```
Run `xmipp` script in the root folder (it might be necessary to add execute permission via `chmod +x xmipp`)
This script will checkout additional repositories and build Xmipp for you.

## Troubleshooting

There are some dependencies that are required. If in the Xmipp compilation some dependences are needed run
```
sudo apt-get install libsqlite3-dev libtiff5-dev libhdf5-dev
```
