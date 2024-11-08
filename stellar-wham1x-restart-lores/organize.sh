if [ ! -d "PositivityFourMoments" ]; then
  mkdir PositivityFourMoments
fi
mv *positivity_shift_FourMoments_* PositivityFourMoments/

if [ ! -d "BiMaxwellianMoments" ]; then
  mkdir BiMaxwellianMoments
fi
mv *_BiMaxwellianMoments_* BiMaxwellianMoments/

if [ ! -d "Field" ]; then
  mkdir Field
fi
mv *-field_* Field/

if [ ! -d "Geometry" ]; then
  mkdir Geometry
fi

mv *-b_i* Geometry/
mv *-bcart* Geometry/
mv *-bmag* Geometry/
mv *-cmag* Geometry/
mv *-dxdz* Geometry/
mv *-dzdx* Geometry/
mv *-eps2* Geometry/
mv *-g_ij* Geometry/
mv *-gij* Geometry/
mv *-gxxj* Geometry/
mv *-gxyj* Geometry/
mv *-gxzj* Geometry/
mv *-gyyj* Geometry/
mv *-jacobgeo* Geometry/
mv *-jacobtot* Geometry/
mv *-mapc2p* Geometry/
mv *-nodes* Geometry/
mv *-normals* Geometry/

if [ ! -d "Slurmscripts" ]; then
  mkdir Slurmscripts
fi
mv *.out Slurmscripts/