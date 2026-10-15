# This directory contains NonLinLoc configuration files that can be used for generating the travel-time grids for ISNet.
# PRESTo only needs and uses the final binary files in the time/ directory to locate all the earthquakes.

# How to modify and regenerate the grids if needed:
# Install NonLinLoc [http://alomax.free.fr/nlloc/ - Cygwin can be used under Windows]
# Make sure the model/ and time/ dirs exist and are empty
mkdir model time
rm model/* time/*

# - Modify run/ISNet.run: TRANS  LAMBERT  Clarke-1880  44.5 10.2 43.1993 44.9961  0.0
#   with:
#    projection type, RefEllipsoid, [e.g. Lambert Conformal Conic, Clarke 1880]
#    LatOrig, LongOrig, [e.g. center of the grid]
#    FirstStdParal, SecondStdParal, [best set as minLat, maxLat for minimun distortion]
#    RotCW [grid rotation angle, usually 0]
#
# - Modify Vel2Grid_P.in: VGGRID  161 111 43 -110.0 -60.0 -2.0  1.0 1.0 1.0  SLOW_LEN
#   with:
#    numPointsX(=numKmLon/dX+1), numPointsY(=numKmLat/dY+1), numPointsZ(=numKmVertical/dZ+1),
#    minX_km(minLon-origLon in km), minY_km, minZ_km, [origin is the at the center of the grid, depth 0]
#    dX(in km), dY, dZ [grid step in km]
#    SLOW_LEN [<- i.e. time grid]
#   and change velocity model

cp Vel2Grid_P.in Vel2Grid_S.in
# Modify Vel2Grid_S.in: VGTYPE P -> VGTYPE S

# Change stations in Gri2Time_P.in
cp Grid2Time_P.in Grid2Time_S.in
# Modify Grid2Time_S.in: GTFILES  ./model/layer  ./time/layer P ->  S

# Uncomment P includes in run/ISNet.run
# Comment S includes in run/ISNet.run
Vel2Grid run/ISNet.run
Grid2Time run/ISNet.run

# Uncomment S includes in run/ISNet.run
# Comment P includes in run/ISNet.run
Vel2Grid run/ISNet.run
Grid2Time run/ISNet.run
