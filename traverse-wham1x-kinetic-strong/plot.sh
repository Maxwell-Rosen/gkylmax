
# name1="outputs/gk_mirror_adiabatic_elc_1x2v_p1_nosource_uniform"
# name2="outputs/gk_mirror_adiabatic_elc_1x2v_p1_nosource_nonuniform"
name="gk_wham"
# name="outputs/gk_mirror_adiabatic_elc_1x2v_p1_nosource_nonuniform"
species="ion"

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
frame=32


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
pgkyl "$name-"jacobgeo.gkyl interp -b ms -p1 pl --title "jacobgeo"&
pgkyl "$name-"jacobtot.gkyl interp -b ms -p1 pl --title "jacobtot"&
pgkyl "$name-"jacobtot_inv.gkyl interp -b ms -p1 pl --title "jacobtot_inv"&
pgkyl "$name-"jacogeo_inv.gkyl interp -b ms -p1 pl --title "jacobgeo_inv"&
pgkyl "$name-"b_i.gkyl interp -b ms -p1 pl --title "b_i"&
pgkyl "$name-"mapc2p.gkyl interp -b ms -p1 pl --title "mapc2p"&
pgkyl "$name-"bmag.gkyl interp -b ms -p1 pl --title "bmag"&
pgkyl "$name-"bmag_inv.gkyl interp -b ms -p1 pl --title "bmag_inv"&
pgkyl "$name-"bmag_inv_sq.gkyl interp -b ms -p1 pl --title "bmag_inv_sq"&
pgkyl "$name-"cmag.gkyl interp -b ms -p1 pl --title "cmag"&

# Distribution function at a single velocity space point at all z
# "$name-"$species"_[0-9]*.gkyl"\
#   interp -b gkhyb -p1 sel --z1 5e6 --z2 3.5e-15 ev 'f[:] abs' collect pl --logz --zmin 1e-4 \

#Plot phi
frame=117
# pgkyl "$name-field_[0-9]*.gkyl" interp -b ms -p1 collect sel --z1 0 pl --title 'phi' &
# pgkyl "$name-field_[0-9]*.gkyl" interp -b ms -p1 collect pl --title 'phi' &
# pgkyl "$name-field_$frame.gkyl" interp -b ms -p1 ev 'f[:] 940 /' pl --title "phi at frame $frame in unite e phi/Te" &


if [ "$species" = "elc" ]; then mass=9.11e-31
elif [ "$species" = "ion" ]; then mass=3.34e-27
else echo "species must be ion or elc"
fi
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

