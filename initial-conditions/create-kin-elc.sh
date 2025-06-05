simFolder="../stellar-wham1x-init-kinetic-from-boltzmann"
targetFolder="kinet-elc-288z-nu2000"
frame=144

if [ ! -d "$simFolder" ]; then
  echo "Simulation folder $simFolder does not exist."
  exit 1
fi
if [ ! -d "$targetFolder" ]; then
  echo "Creating target folder $targetFolder."
  mkdir "$targetFolder"
else
  echo "Target folder $targetFolder already exists."
fi

cp "$simFolder/Field/gk_wham-field_${frame}.gkyl" "$targetFolder"
cp "$simFolder/Distributions/gk_wham-ion_${frame}.gkyl" "$targetFolder"
cp "$simFolder/Distributions/gk_wham-elc_${frame}.gkyl" "$targetFolder"
cp "$simFolder/BiMaxwellianMoments/gk_wham-ion_BiMaxwellianMoments_${frame}.gkyl" "$targetFolder"
cp "$simFolder/BiMaxwellianMoments/gk_wham-elc_BiMaxwellianMoments_${frame}.gkyl" "$targetFolder"
cp "$simFolder/M/gk_wham-ion_M0_${frame}.gkyl" "$targetFolder"
cp "$simFolder/M/gk_wham-elc_M0_${frame}.gkyl" "$targetFolder"
cp "$simFolder/misc/gk_wham-ion_jacobvel.gkyl" "$targetFolder"
cp "$simFolder/misc/gk_wham-elc_jacobvel.gkyl" "$targetFolder"
cp "$simFolder/Geometry/gk_wham-bmag.gkyl" "$targetFolder"
cp "$simFolder/Geometry/gk_wham-jacobgeo.gkyl" "$targetFolder"
cp "$simFolder/Geometry/gk_wham-jacobtot.gkyl" "$targetFolder"
cp "$simFolder/Geometry/gk_wham-mapc2p.gkyl" "$targetFolder"
cp "$simFolder/Geometry/gk_wham-mc2nu_pos.gkyl" "$targetFolder"
cp "$simFolder/Geometry/gk_wham-nodes.gkyl" "$targetFolder"