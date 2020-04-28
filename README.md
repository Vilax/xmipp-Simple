[![Build Status](https://travis-ci.com/I2PC/xmipp.svg?branch=devel)](https://travis-ci.com/I2PC/xmipp)
<!---  [![Quality Gate](https://sonarcloud.io/api/project_badges/measure?project=Xmipp&metric=alert_status)](https://sonarcloud.io/dashboard?id=Xmipp)
[![Technical debt](https://sonarcloud.io/api/project_badges/measure?project=Xmipp&metric=sqale_index)](https://sonarcloud.io/component_measures?id=Xmipp&metric=sqale_index)
[![Bugs](https://sonarcloud.io/api/project_badges/measure?project=Xmipp&metric=bugs)](https://sonarcloud.io/project/issues?id=Xmipp&resolved=false&types=BUG)
--->
# Xmipp-Simple

Welcome to Xmipp. Xmipp is a suite of image processing programs, primarily aimed at single-particle 3D electron microscopy. 

This is a light version of Xmipp, which only contains the basic tools addres to estimate resolution measuremetns

## Getting started

Start by cloning the repository from GitHub and go there.
```
git clone https://github.com/I2PC/xmipp xmipp-bundle
cd xmipp-bundle
```

Run `xmipp` script in the root folder via Scipion (it might be necessary to add execute permission via `chmod +x xmipp`)
This script will checkout additional repositories and build Xmipp for you.

You can see the whole usage of the script with `./xmipp --help`. The most useful options are `br=branch_name` to select a specific branch to be checkout-ed, and `N=#processors` to use for the build.
