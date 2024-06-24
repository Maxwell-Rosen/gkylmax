
# name1="outputs/gk_mirror_adiabatic_elc_1x2v_p1_nosource_uniform"
# name2="outputs/gk_mirror_adiabatic_elc_1x2v_p1_nosource_nonuniform"
name="outputs/gk_wham"
species="ion"

pgkyl "$name"-"$species"_0.gkyl -t f0 \
  "$name"-"$species"_1.gkyl -t f1 \
  interp -b gkhyb -p1 \
  ev -t dfdt "f1 f0 - 6.26254e-12 /"\
  activate -t dfdt select --z0 0.001 --z1 0.0 \
  pl --title "dfdt at psi = 0.001 z = 0.0"&

pgkyl "$name"-"$species"_0.gkyl -t f0 \
  "$name"-"$species"_1.gkyl -t f1 \
  interp -b gkhyb -p1 \
  ev -t dfdt "f1 f0 - 6.26254e-12 /"\
  activate -t dfdt select --z0 0.001 integrate 3 \
  pl --title "dfdt at psi = 0.001 integrate 3"&

  pgkyl "$name"-"$species"_0.gkyl -t f0 \
  "$name"-"$species"_1.gkyl -t f1 \
  interp -b gkhyb -p1 \
  ev -t dfdt "f1 f0 - 6.26254e-12 /"\
  activate -t dfdt select --z0 0.001 integrate 2 \
  pl --title "dfdt at psi = 0.001 integrate 2"&
# name="outputs/gk_mirror_adiabatic_elc_1x2v_p1_nosource_nonuniform"
# species="ion"
# Make a loop to set species to "ion" and "elc"

# Make a loop over frame values between 0 and 10
# for frame in {0..42}
# do
#   echo "frame $frame"

# for species in "elc" "ion"
# do
#   echo $species
# saveLoc="python-plots-frame-order/gk_wham-$frame"



# if [ "$species" = "elc" ]; then mass=9.11e-31
# elif [ "$species" = "ion" ]; then mass=3.34e-27
# else echo "species must be ion or elc"
# fi
# psival=0.001


# if [ "$species" = "elc" ]; then
#   pgkyl "$name-elc_M0_$frame.gkyl" -t M0 "$name-elc_M1_$frame.gkyl" -t M1 \
#     "$name-elc_M2par_$frame.gkyl" -t M2par "$name-elc_M2perp_$frame.gkyl" -t M2perp \
#     "$name-field_$frame.gkyl" -t phi \
#     activate -t M0,M1,M2par,M2perp interp -b ms -p1 select --z0 $psival ev -a -t res "M2par M1 M1 * M0 / - $mass * M0 / 1.6e-19 /" \
#     activate -t M0,M2perp ev -a -t res2 "$mass M2perp M0 / * 1.6e-19 /"\
#     activate -t res,res2 ev -a -t res3 "res2 2 * res + 3 /" \
#     activate -t phi interp -b ms -p1 select -t phi_interp --z0 $psival \
#     activate -t res3 select -t T0 --z0 0.0 --z1 0.0 \
#     activate -t T0,phi_interp ev -a -t ephiTe "phi T0 /" \
#     activate -t ephiTe select -t plotfit --z0 $psival \
#     activate -t plotfit pl --title "Potential at psi = $psival" -y "e phi / T" -x "Field line length (m)"\
#      --saveas "$saveLoc-1d-ephiTe-$frame.png" --no-show&
    
#   pgkyl "$name-elc_M0_$frame.gkyl" -t M0 "$name-elc_M1_$frame.gkyl" -t M1 \
#     "$name-elc_M2par_$frame.gkyl" -t M2par "$name-elc_M2perp_$frame.gkyl" -t M2perp \
#     "$name-field_$frame.gkyl" -t phi \
#     activate -t M0,M1,M2par,M2perp interp -b ms -p1 ev -a -t res "M2par M1 M1 * M0 / - $mass * M0 / 1.6e-19 /" \
#     activate -t M0,M2perp ev -a -t res2 "$mass M2perp M0 / * 1.6e-19 /" \
#     activate -t res,res2 ev -a -t res3 "res2 2 * res + 3 /" \
#     activate -t phi interp -t phi_interp -b ms -p1  \
#     activate -t res3 select -t T0 --z1 0.0 \
#     activate -t T0,phi_interp ev -a -t ephiTe "phi_interp T0 /" \
#     activate -t ephiTe pl --title "Potential" -y "Field line length (m)" -x "psi" --clabel "e phi / T" \
#     --saveas "$saveLoc-2d-ephiTe-$frame.png" --no-show &

#   pgkyl "$name-field_$frame.gkyl" interp -b ms -p1 select --z0 $psival pl --title "phi" \
#     -x "Field line length (m)" -y "Electric potential (V)" --saveas "$saveLoc-1d-field-$frame.png" --no-show&

#   pgkyl "$name-field_$frame.gkyl" interp -b ms -p1 pl --title "phi" \
#     -y "Field line length (m)" -x "Psi" --clabel "Electric potential (V)" --saveas "$saveLoc-2d-field-$frame.png" --no-show&

  # pgkyl "test_phi_pol.gkyl" interp -b mt -p2 select --z0 $psival pl --title "phi" \
  #   -x "Field line length (m)" -y "Electric potential (V)" --saveas "python-plots/test_phi_pol-1d-field.png" --no-show&

  # pgkyl "test_phi_pol.gkyl" interp -b mt -p2 pl --title "phi" \
  #   -y "Field line length (m)" -x "Psi" --clabel "Electric potential (V)" --saveas "python-plots/test_phi_pol-2d-field.png" --no-show&

  # pgkyl "$name-field_$frame.gkyl" -t Field test_phi_pol.gkyl -t phi activate -t Field interp -b ms -p1 select -t Field_proj \
  # --z1 0.0 activate -t phi interp -b mt -p2 select -t phi_proj --z1 0.0 activate -t phi_proj,Field_proj pl -f0 --title \
  # "phi (blue) vs phi_pol (orange)" --saveas "$saveLoc-$frame-1d-center-field.png" --no-show &
# fi
# # 1D plots of distribution function at a certain psival
# pgkyl "$name-"$species"_$frame.gkyl" interp -b gkhyb -p1 select --z0 $psival integrate 2 pl --title "$species distribution function integrating over vpar"  \
#   --xscale 0.79233226837 -x "Field line length (m)" -y "Magnetic moment" --saveas "$saveLoc-"$species"_$frame-1d-vpar.png" --no-show &
# pgkyl "$name-"$species"_$frame.gkyl" interp -b gkhyb -p1 select --z0 $psival integrate 3 pl --title "$species distribution function integrating over magnetic moment" \
#   --xscale 0.79233226837 -x "Field line length (m)" -y "Parallel velocity (m/s)" --saveas "$saveLoc-"$species"_$frame-1d-mu.png" --no-show &

# pgkyl "$name-"$species"_$frame.gkyl" interp -b gkhyb -p1 select --z0 $psival --z1 0.0001 pl --title "$species distribution function at z=0.0001" \
#   --xscale 0.79233226837 -x "Parallel velocity (m/s)" -y "Magnetic moment" --logz --saveas "$saveLoc-"$species"_$frame-1d-log-z=0.0001.png" --no-show &
# pgkyl "$name-"$species"_$frame.gkyl" interp -b gkhyb -p1 select --z0 $psival --z1 0.743 pl --title "$species distribution function at z=0.59 m" \
#   --xscale 0.79233226837 -x "Parallel velocity (m/s)" -y "Magnetic moment" --logz --saveas "$saveLoc-"$species"_$frame-1d-log-z=0.59.png" --no-show &
# pgkyl "$name-"$species"_$frame.gkyl" interp -b gkhyb -p1 select --z0 $psival --z1 1.224 pl --title "$species distribution function at z=0.97 m" \
#   --xscale 0.79233226837 -x "Parallel velocity (m/s)" -y "Magnetic moment" --logz --saveas "$saveLoc-"$species"_$frame-1d-log-z=0.97.png" --no-show &

# pgkyl "$name-"$species"_$frame.gkyl" interp -b gkhyb -p1 select --z0 $psival --z1 0.0001 pl --title "$species distribution function at z=0.0001" \
#   --xscale 0.79233226837 -x "Parallel velocity (m/s)" -y "Magnetic moment" --saveas "$saveLoc-"$species"_$frame-1d-z=0.0001.png" --no-show &
# pgkyl "$name-"$species"_$frame.gkyl" interp -b gkhyb -p1 select --z0 $psival --z1 0.0743 pl --title "$species distribution function at z=0.59 m" \
#   --xscale 0.79233226837 -x "Parallel velocity (m/s)" -y "Magnetic moment" --saveas "$saveLoc-"$species"_$frame-1d-z=0.59.png" --no-show &
# pgkyl "$name-"$species"_$frame.gkyl" interp -b gkhyb -p1 select --z0 $psival --z1 1.224 pl --title "$species distribution function at z=0.97 m" \
#   --xscale 0.79233226837 -x "Parallel velocity (m/s)" -y "Magnetic moment" --saveas "$saveLoc-"$species"_$frame-1d-z=0.97.png" --no-show &

# 1D plots of moments at a certain psival
# pgkyl "$name-$species"_prim_moms"_$frame.gkyl" interp -b ms -p1 select --z0 $psival pl --title "prim_moms" \
#   --xscale 0.79233226837 -x "Field line length (m)" \
#   --saveas "$saveLoc-"$species"-1d-prim_moms-$frame.png" --no-show&
# pgkyl "$name-field_$frame.gkyl" interp -b ms -p1 select --z0 $psival pl --title "phi" \
#   --xscale 0.79233226837 -x "Field line length (m)" -y "Electric potential (V)" --saveas "$saveLoc-1d-field-$frame.png" --no-show &
# pgkyl "$name-"$species"_M0_$frame.gkyl" -t M0 "$name-"$species"_M1_$frame.gkyl" -t M1 \
#   "$name-"$species"_M2par_$frame.gkyl" -t M2par activate -t M0,M1,M2par\
#   interp -b ms -p1 select --z0 $psival ev -a -t res "M2par M1 M1 * M0 / - $mass * M0 / 1.6e-19 /" \
#   activate -t res pl --title "$species Tpar (eV)" \
#   --xscale 0.79233226837 -x "Field line length (m)" -y "Temperature (eV)" --no-show \
#    --saveas "$saveLoc-"$species"-1d-Tpar-$frame.png" &
# pgkyl "$name-"$species"_M0_$frame.gkyl" -t M0 "$name-"$species"_M2perp_$frame.gkyl" -t M2perp \
#   activate -t M0,M2perp interp -b ms -p1 select --z0 $psival ev -a -t res "$mass M2perp M0 / * 1.6e-19 /" \
#   activate -t res pl --title "$species Tperp (eV)" \
#   --xscale 0.79233226837 -x "Field line length (m)" -y "Temperature (eV)" \
#   --saveas "$saveLoc-"$species"-1d-Tperp-$frame.png" --no-show&
# pgkyl "$name-"$species"_M0_$frame.gkyl" -t M0 "$name-"$species"_M1_$frame.gkyl" -t M1 \
#   "$name-"$species"_M2par_$frame.gkyl" -t M2par "$name-"$species"_M2perp_$frame.gkyl" -t M2perp \
#   activate -t M0,M1,M2par,M2perp interp -b ms -p1 select --z0 $psival ev -a -t res "M2par M1 M1 * M0 / - $mass * M0 / 1.6e-19 /" \
#   activate -t M0,M2perp ev -a -t res2 "$mass M2perp M0 / * 1.6e-19 /"\
#   activate -t res,res2 ev -a -t res3 "res2 2 * res + 3 /" \
#   activate -t res3 pl --title "$species T (eV)" \
#   --xscale 0.79233226837 -x "Field line length (m)" -y "Temperature (eV)"\
#    --saveas "$saveLoc-"$species"-1d-T-$frame.png" --no-show &
# pgkyl "$name-"$species"_M0_$frame.gkyl" interp -b ms -p1 select --z0 $psival pl --title "$species Density" \
#   --xscale 0.79233226837 -x "Field line length (m)" -y "Density (m^-3)" --saveas "$saveLoc-"$species"-1d-density-$frame.png" --no-show &
# pgkyl "$name-"$species"_M0_$frame.gkyl" "$name-"$species"_M1_$frame.gkyl" \
#   interp -b ms -p1 select --z0 $psival ev "f[1] f[0] /" pl --title "$species Upar" \
#   --xscale 0.79233226837 -x "Field line length (m)" -y "Parallel velocity (m/s)" --saveas "$saveLoc-"$species"-1d-upar-$frame.png" --no-show &

# # 2D plots of moments
# pgkyl "$name-field_$frame.gkyl" interp -b ms -p1 pl --title "phi" \
#   --yscale 0.79233226837 -y "Field line length (m)" -x "Psi" --clabel "Electric potential (V)" --saveas "$saveLoc-2d-field-$frame.png" --no-show &
# pgkyl "$name-"$species"_M0_$frame.gkyl" -t M0 "$name-"$species"_M1_$frame.gkyl" -t M1 \
#   "$name-"$species"_M2par_$frame.gkyl" -t M2par activate -t M0,M1,M2par\
#   interp -b ms -p1 ev -a -t res "M2par M1 M1 * M0 / - $mass * M0 / 1.6e-19 /" \
#   activate -t res pl --title "$species Tpar (eV)" \
#   --yscale 0.79233226837 -y "Field line length (m)" -x "Psi" --clabel "Temperature (eV)" --saveas "$saveLoc-"$species"-2d-Tpar-$frame.png" --no-show &
# pgkyl "$name-"$species"_M0_$frame.gkyl" -t M0 "$name-"$species"_M2perp_$frame.gkyl" -t M2perp \
#   activate -t M0,M2perp interp -b ms -p1 ev -a -t res "$mass M2perp M0 / * 1.6e-19 /" \
#   activate -t res pl --title "$species Tperp (eV)" \
#   --yscale 0.79233226837 -y "Field line length (m)" -x "Psi" --clabel "Temperature (eV)" --saveas "$saveLoc-"$species"-2d-Tperp-$frame.png" --no-show &
# pgkyl "$name-"$species"_M0_$frame.gkyl" -t M0 "$name-"$species"_M1_$frame.gkyl" -t M1 \
#   "$name-"$species"_M2par_$frame.gkyl" -t M2par "$name-"$species"_M2perp_$frame.gkyl" -t M2perp \
#   activate -t M0,M1,M2par,M2perp interp -b ms -p1 ev -a -t res "M2par M1 M1 * M0 / - $mass * M0 / 1.6e-19 /" \
#   activate -t M0,M2perp ev -a -t res2 "$mass M2perp M0 / * 1.6e-19 /" \
#   activate -t res,res2 ev -a -t res3 "res2 2 * res + 3 /" \
#   activate -t res3 pl --title "$species T (eV)" \
#   --yscale 0.79233226837 -y "Field line length (m)" -x "Psi" --clabel "Temperature (eV)"\
#   --saveas "$saveLoc-"$species"-2d-T-$frame.png" --no-show &
# pgkyl "$name-"$species"_M0_$frame.gkyl" interp -b ms -p1 pl --title "$species Density" \
#   --yscale 0.79233226837 -y "Field line length (m)" -x "Psi" --clabel "Density (m^-3)" --logz --saveas "$saveLoc-"$species"-2d-log-density-$frame.png" --no-show &
# pgkyl "$name-"$species"_M0_$frame.gkyl" interp -b ms -p1 pl --title "$species Density" \
#   --yscale 0.79233226837 -y "Field line length (m)" -x "Psi" --clabel "Density (m^-3)" --saveas "$saveLoc-"$species"-2d-density-$frame.png" --no-show &
# pgkyl "$name-"$species"_M0_$frame.gkyl" "$name-"$species"_M1_$frame.gkyl" \
#   interp -b ms -p1  ev "f[1] f[0] /" pl --title "$species Upar" \
#   --yscale 0.79233226837 -y "Field line length (m)" -x "Psi" --clabel "Parallel velocity (m/s)" --saveas "$saveLoc-"$species"-2d-upar-$frame.png" --no-show &
# sleep 2s
# done
# done

# # Plot tperp at all frames
# pgkyl "$name-"$species"_M0_[0-9]*.gkyl" -t M0 "$name-"$species"_M2perp_[0-9]*.gkyl" -t M2perp\
#   activate -t M0,M2perp interp -b ms -p1 ev -a -t res \
#   "$mass M2perp M0 / * 1.6e-19 /" activate -t res collect pl  --title 'Tperp (ev)'&

# # Plot tpar at all frames
# pgkyl "$name-"$species"_M0_[0-9]*.gkyl" -t M0 "$name-"$species"_M1_[0-9]*.gkyl"\
#  -t M1 "$name-"$species"_M2par_[0-9]*.gkyl" -t M2par\
#   activate -t M0,M1,M2par interp -b ms -p1 ev -a -t res \
#   "M2par M1 M1 * M0 / - $mass * M0 / 1.6e-19 /" activate -t res\
#   collect pl  --title 'Tpar (ev)' --zmin 0.0 &



# pgkyl "$name-"$species"_$frame.gkyl" interp -b gkhyb -p1 select --z0 0.005 --z1 0.00001 pl --logz --title "z=0 distribution function" &
# pgkyl "$name-"$species"_$frame.gkyl" interp -b gkhyb -p1 select --z0 0.005 --z1 0.1 pl --logz --title "z=0.1 distribution function" &
# pgkyl "$name-"$species"_$frame.gkyl" interp -b gkhyb -p1 select --z0 0.005 --z1 0.2 pl --logz --title "z=0.2 distribution function" &
# pgkyl "$name-"$species"_$frame.gkyl" interp -b gkhyb -p1 select --z0 0.005 --z1 0.3 pl --logz --title "z=0.3 distribution function" &
# pgkyl "$name-"$species"_$frame.gkyl" interp -b gkhyb -p1 select --z0 0.005 --z1 0.4 pl --title "z=0.4 distribution function" &
# pgkyl "$name-"$species"_$frame.gkyl" interp -b gkhyb -p1 select --z0 0.005 --z1 0.5 pl --title "z=0.5 distribution function" &
# pgkyl "$name-"$species"_$frame.gkyl" interp -b gkhyb -p1 select --z0 0.005 --z1 0.6 pl --title "z=0.6 distribution function" &
# pgkyl "$name-"$species"_$frame.gkyl" interp -b gkhyb -p1 select --z0 0.005 --z1 0.7 pl --title "z=0.7 distribution function" &
# pgkyl "$name-"$species"_$frame.gkyl" interp -b gkhyb -p1 select --z0 0.005 --z1 0.8 pl --title "z=0.8 distribution function" &
# pgkyl "$name-"$species"_$frame.gkyl" interp -b gkhyb -p1 select --z0 0.005 --z1 0.9 pl --title "z=0.9 distribution function" &

# pgkyl outputs/gk_mirror_adiabatic_elc_1x2v_p1_true_maxwellian-ion_0.gkyl --c2p outputs/mapc2p.gkyl \
#  interp -b gkhyb -p1 integrate 1 pl &
# pgkyl outputs/gk_mirror_adiabatic_elc_1x2v_p1_true_maxwellian-ion_0.gkyl \
#  interp -b gkhyb -p1 integrate 1 pl &


# Animations of distribution functions with vpar on z and mu on z
# pgkyl "$name-ion_[0-9]*.gkyl"\
#   interp -b gkhyb -p1 integrate 2 ev 'f[:] abs' animate --logz --zmin 1e-20 --fps 4 \
#   --saveas "vpar.mp4" &
# pgkyl "$name-ion_[0-9]*.gkyl"\
#   interp -b gkhyb -p1 integrate 1 ev 'f[:] abs' animate --logz --zmin 1e-4 --fps 4\
#   --saveas "mu.mp4" &
# pgkyl "$name-ion_[0-9]*.gkyl"\
#   interp -b gkhyb -p1 sel --z0 0.0 ev 'f[:] abs' animate --logz --zmin 1e-10 --fps 4\
#   --saveas "z=0,0.mp4" &
# pgkyl "$name-ion_[0-9]*.gkyl"\
#   interp -b gkhyb -p1 sel --z0 0.70 ev 'f[:] abs' animate --logz --zmin 1e-10 --fps 4\
#   --saveas "z=0,7.mp4" &

# pgkyl "$name-ion_0.gkyl" "$name-ion_1.gkyl" "$name-ion_2.gkyl" "$name-ion_3.gkyl" \
# "$name-ion_4.gkyl" "$name-ion_5.gkyl" "$name-ion_6.gkyl" "$name-ion_7.gkyl" \
# "$name-ion_8.gkyl" "$name-ion_9.gkyl" "$name-ion_10.gkyl"\
#   interp -b gkhyb -p1 sel --z0 0.7 ev 'f[:] abs' \
#   animate --logz --zmin 1e-10 --fps 4 --title "z=7" \
#   --saveas "$name z=0,7.mp4"&


# Plot single frames of the distribution function
# pgkyl "$name-ion_71.gkyl" interp -b gkhyb -p1 integrate 2 ev 'f[:] abs' sel --z0 0.0 pl --logy &
frame=0


# pgkyl "$name-"$species"_$frame.gkyl" interp -b gkhyb -p1 sel --z0 0.0 ev 'f[:] abs' pl --logz --zmin 1e-10 --title "frame $frame z=0"&
# pgkyl "$name-"$species"_$frame.gkyl" interp -b gkhyb -p1 sel --z0 0.90 ev 'f[:] abs' pl --logz --zmin 1e-10 --title "frame $frame z=.9"&
# pgkyl "$name-"$species"_$frame.gkyl" interp -b gkhyb -p1 sel --z0 0.7 ev 'f[:] abs' pl --logz --zmin 1e-10 --title "frame $frame z=1.0"&
#Taking absolute value
# pgkyl "$name-"$species"_$frame.gkyl" interp -b gkhyb -p1 integrate 2 ev 'f[:] abs' pl --logz --zmin 1e-20 --title "frame $frame vpar"&
# pgkyl "$name-"$species"_$frame.gkyl" interp -b gkhyb -p1 integrate 1 ev 'f[:] abs' pl --logz --zmin 1e-4 --title "frame $frame mu"&
# pgkyl "$name-"$species"_$frame.gkyl" interp -b gkhyb -p1 sel --z0 0.0 ev 'f[:] abs' pl --logz --zmin 1e-10 --title "frame $frame z=0"&
#Without absolute value
# pgkyl "$name-"$species"_$frame.gkyl" interp -b gkhyb -p1 integrate 2 pl --logz --zmin 1e-20 --title "frame $frame vpar"&
# pgkyl "$name-"$species"_$frame.gkyl" interp -b gkhyb -p1 integrate 1 pl --logz --zmin 1e-4 --title "frame $frame mu"&
# pgkyl "$name-"$species"_$frame.gkyl" interp -b gkhyb -p1 sel --z0 0.0 pl --logz --zmin 1e-10 --title "frame $frame z=0"&

# Plot geometry quantities
# echo "geometry"
# pgkyl "$name-"jacobgeo.gkyl interp -b ms -p1 pl --title "jacobgeo" --saveas "$saveLoc-geo-jacobgeo.png" --no-show &
# pgkyl "$name-"jacobtot.gkyl interp -b ms -p1 pl --title "jacobtot" --saveas "$saveLoc-geo-jacobtot.png" --no-show &
# pgkyl "$name-"jacobtot_inv.gkyl interp -b ms -p1 pl --title "jacobtot_inv" --saveas "$saveLoc-geo-jacobtot_inv.png" --no-show &
# pgkyl "$name-"jacobgeo_inv.gkyl interp -b ms -p1 pl --title "jacobgeo_inv" --saveas "$saveLoc-geo-jacobgeo_inv.png" --no-show &
# pgkyl "$name-"b_i.gkyl interp -b ms -p1 pl --title "b_i" --saveas "$saveLoc-geo-b_i.png" --no-show &
# pgkyl "$name-"mapc2p.gkyl interp -b ms -p1 pl --title "mapc2p" --saveas "$saveLoc-geo-mapc2p.png" --no-show &
# pgkyl "$name-"bmag.gkyl interp -b ms -p1 pl --title "bmag" --saveas "$saveLoc-geo-bmag.png" --no-show &
# pgkyl "$name-"bmag_inv.gkyl interp -b ms -p1 pl --title "bmag_inv" --saveas "$saveLoc-geo-bmag_inv.png" --no-show &
# pgkyl "$name-"bmag_inv_sq.gkyl interp -b ms -p1 pl --title "bmag_inv_sq" --saveas "$saveLoc-geo-bmag_inv_sq.png" --no-show &
# pgkyl "$name-"cmag.gkyl interp -b ms -p1 pl --title "cmag" --saveas "$saveLoc-geo-cmag.png" --no-show &

# Distribution function at a single velocity space point at all z
# "$name-"$species"_[0-9]*.gkyl"\
#   interp -b gkhyb -p1 sel --z1 5e6 --z2 3.5e-15 ev 'f[:] abs' collect pl --logz --zmin 1e-4 \

#Plot phi
frame=117
# pgkyl "$name-field_[0-9]*.gkyl" interp -b ms -p1 collect sel --z1 0 pl --title 'phi' &
# pgkyl "$name-field_[0-9]*.gkyl" interp -b ms -p1 collect pl --title 'phi' &
# pgkyl "$name-field_$frame.gkyl" interp -b ms -p1 ev 'f[:] 940 /' pl --title "phi at frame $frame in unite e phi/Te" &



# Moments of the distribution function at a single frame
frame=1
# pgkyl $name-"$species"_M0_$frame.gkyl interp -b ms -p1 pl --title 'Density' &
# pgkyl $name-"$species"_M0_$frame.gkyl $name-"$species"_M1_$frame.gkyl interp -b ms -p1 ev "f[1] f[0] /" pl --title 'Upar' &
# pgkyl $name-"$species"_prim_moms_$frame.gkyl interp -b ms -p1 pl --title 'prim--moms' &
# pgkyl $name-"$species"_M0_$frame.gkyl $name-"$species"_M2perp_$frame.gkyl interp -b ms -p1 ev "$mass f[1] f[0] / * 1.6e-19 /" pl --title 'Tperp (eV)' &
# pgkyl $name-"$species"_M0_$frame.gkyl $name-"$species"_M1_$frame.gkyl $name-"$species"_M2par_$frame.gkyl interp -b ms -p1 ev "f[2] f[1] f[1] * f[0] / - $mass * f[0] / 1.6e-19 /" pl --title 'Tpar (eV)' &

# # Plot tperp at all frames
# pgkyl "$name-"$species"_M0_[0-9]*.gkyl" -t M0 "$name-"$species"_M2perp_[0-9]*.gkyl" -t M2perp\
#   activate -t M0,M2perp interp -b ms -p1 ev -a -t res \
#   "$mass M2perp M0 / * 1.6e-19 /" activate -t res collect pl  --title 'Tperp (ev)'&

# # Plot tpar at all frames
# pgkyl "$name-"$species"_M0_[0-9]*.gkyl" -t M0 "$name-"$species"_M1_[0-9]*.gkyl"\
#  -t M1 "$name-"$species"_M2par_[0-9]*.gkyl" -t M2par\
#   activate -t M0,M1,M2par interp -b ms -p1 ev -a -t res \
#   "M2par M1 M1 * M0 / - $mass * M0 / 1.6e-19 /" activate -t res\
#   collect pl  --title 'Tpar (ev)' --zmin 0.0 &

# # Plot density at all frames
# pgkyl "$name-"$species"_M0_[0-9]*.gkyl" interp -b ms -p1 collect pl --title 'Density' &

# # Plot Upar at all frames
# pgkyl "$name-"$species"_M0_[0-9]*.gkyl" -t M0 "$name-"$species"_M1_[0-9]*.gkyl"\
#  -t M1 activate -t M0,M1 interp -b ms -p1 ev -a -t res "M1 M0 /" activate -t res\
#   collect pl  --title 'Upar'&

# Compare moments at two frames
# frame1=0
# frame2=40
# #Density
# pgkyl "$name-"$species"_M0_$frame1.gkyl" -t M00 "$name-"$species"_M0_$frame2.gkyl" \
#   -t M01 interp -b ms -p1 pl --title 'Density' -f0 &
# #Upar
# pgkyl "$name-"$species"_M0_$frame1.gkyl" -t M00 "$name-"$species"_M1_$frame1.gkyl" -t M10\
#   "$name-"$species"_M0_$frame2.gkyl" -t M01 "$name-"$species"_M1_$frame2.gkyl" -t M11\
#   activate -t M00,M10 interp -b ms -p1 ev -a -t res0 -l "Frame $frame1" "M10 M00 /"\
#   activate -t M01,M11 interp -b ms -p1 ev -a -t res1 -l "Frame $frame2" "M11 M01 /"\
#   activate -t res0,res1 pl  --title 'Upar' -f0 &
# #Tperp
# pgkyl "$name-"$species"_M0_$frame1.gkyl" -t M00 "$name-"$species"_M2perp_$frame1.gkyl" -t M20\
#   "$name-"$species"_M0_$frame2.gkyl" -t M01 "$name-"$species"_M2perp_$frame2.gkyl" -t M21\
#   activate -t M00,M20 interp -b ms -p1 ev -a -t res0 -l "Frame $frame1" "$mass M20 M00 / * 1.6e-19 /"\
#   activate -t M01,M21 interp -b ms -p1 ev -a -t res1 -l "Frame $frame2" "$mass M21 M01 / * 1.6e-19 /"\
#   activate -t res0,res1 pl  --title 'Tperp (eV)' -f0 &
# #Tpar
# pgkyl "$name-"$species"_M0_$frame1.gkyl" -t M00 "$name-"$species"_M1_$frame1.gkyl" -t M10\
#   "$name-"$species"_M2par_$frame1.gkyl" -t M20 "$name-"$species"_M0_$frame2.gkyl" -t M01\
#   "$name-"$species"_M1_$frame2.gkyl" -t M11 "$name-"$species"_M2par_$frame2.gkyl" -t M21\
#   activate -t M00,M10,M20 interp -b ms -p1 ev -a -t res0 -l \
#   "Frame $frame1" "M20 M10 M10 * M00 / - $mass * M00 / 1.6e-19 /"\
#   activate -t M01,M11,M21 interp -b ms -p1 ev -a -t res1 -l \
#   "Frame $frame2" "M21 M11 M11 * M01 / - $mass * M01 / 1.6e-19 /"\
#   activate -t res0,res1 pl  --title 'Tpar (eV)' -f0 &

# Compare two files at the same frame
frame=32
# #Density
# pgkyl "$name1-"$species"_M0_$frame.gkyl" -t M00 "$name2-"$species"_M0_$frame.gkyl" \
#   -t M01 interp -b ms -p1 pl --title 'Density' -f0 &
# #Upar
# pgkyl "$name1-"$species"_M0_$frame.gkyl" -t M00 "$name1-"$species"_M1_$frame.gkyl" -t M10\
#   "$name2-"$species"_M0_$frame.gkyl" -t M01 "$name2-"$species"_M1_$frame.gkyl" -t M11\
#   activate -t M00,M10 interp -b ms -p1 ev -a -t res0 -l "$name1" "M10 M00 /"\
#   activate -t M01,M11 interp -b ms -p1 ev -a -t res1 -l "$name2" "M11 M01 /"\
#   activate -t res0,res1 pl  --title 'Upar' -f0 &
# #Tperp
# pgkyl "$name1-"$species"_M0_$frame.gkyl" -t M00 "$name1-"$species"_M2perp_$frame.gkyl" -t M20\
#   "$name2-"$species"_M0_$frame.gkyl" -t M01 "$name2-"$species"_M2perp_$frame.gkyl" -t M21\
#   activate -t M00,M20 interp -b ms -p1 ev -a -t res0 -l "$name1" "$mass M20 M00 / * 1.6e-19 /"\
#   activate -t M01,M21 interp -b ms -p1 ev -a -t res1 -l "$name2" "$mass M21 M01 / * 1.6e-19 /"\
#   activate -t res0,res1 pl  --title 'Tperp (eV)' -f0 &
# #Tpar
# pgkyl "$name1-"$species"_M0_$frame.gkyl" -t M00 "$name1-"$species"_M1_$frame.gkyl" -t M10\
#   "$name1-"$species"_M2par_$frame.gkyl" -t M20 "$name2-"$species"_M0_$frame.gkyl" -t M01\
#   "$name2-"$species"_M1_$frame.gkyl" -t M11 "$name2-"$species"_M2par_$frame.gkyl" -t M21\
#   activate -t M00,M10,M20 interp -b ms -p1 ev -a -t res0 -l \
#   "$name1" "M20 M10 M10 * M00 / - $mass * M00 / 1.6e-19 /"\
#   activate -t M01,M11,M21 interp -b ms -p1 ev -a -t res1 -l \
#   "$name2" "M21 M11 M11 * M01 / - $mass * M01 / 1.6e-19 /"\
#   activate -t res0,res1 pl  --title 'Tpar (eV)' -f0 &

