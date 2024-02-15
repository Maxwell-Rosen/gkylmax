
# name="outputs/gk_wham_1x2v_p1_adiabatic"
name="outputs/gk_mirror_adiabatic_elc_1x2v_p1_nonuniform_nosource"
# name="outputs/gk_mirror_adiabatic_elc_1x2v_p1_nonuniform_nosource_mirrorIC"
# name="outputs/gk_mirror_adiabatic_elc_1x2v_p1_nonuniform_nosource_bimaxIC"
species="ion"

# Animations of distribution functions with vpar on z and mu on z
# pgkyl "$name-ion_[0-9]*.gkyl"\
#   interp -b gkhyb -p1 integrate 2 ev 'f[:] abs' animate --logz --zmin 1e-20 --fps 4 \
#   --saveas "$name vpar.mp4" &
# pgkyl "$name-ion_[0-9]*.gkyl"\
#   interp -b gkhyb -p1 integrate 1 ev 'f[:] abs' animate --logz --zmin 1e-4 --fps 4\
#   --saveas "$name mu.mp4" &
# pgkyl "$name-ion_[0-9]*.gkyl"\
#   interp -b gkhyb -p1 sel --z0 0.0 ev 'f[:] abs' animate --logz --zmin 1e-10 --fps 4\
#   --saveas "$name z=0,0.mp4" &
# pgkyl "$name-ion_[0-9]*.gkyl"\
#   interp -b gkhyb -p1 sel --z0 0.70 ev 'f[:] abs' animate --logz --zmin 1e-10 --fps 4\
#   --saveas "$name z=0,7.mp4" &

# pgkyl "$name-ion_0.gkyl" "$name-ion_1.gkyl" "$name-ion_2.gkyl" "$name-ion_3.gkyl" \
# "$name-ion_4.gkyl" "$name-ion_5.gkyl" "$name-ion_6.gkyl" "$name-ion_7.gkyl" \
# "$name-ion_8.gkyl" "$name-ion_9.gkyl" "$name-ion_10.gkyl"\
#   interp -b gkhyb -p1 sel --z0 0.0 ev 'f[:] abs' \
#   animate --logz --zmin 1e-10 --fps 4 --title "z=7" \
#   --saveas "$name z=0,0.mp4"&


# Plot single frames of the distribution function
# pgkyl "$name-ion_71.gkyl" interp -b gkhyb -p1 integrate 2 ev 'f[:] abs' sel --z0 0.0 pl --logy &
frame=1


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
# pgkyl jacobgeo.gkyl interp -b ms -p1 pl --title "jacobgeo"&
# pgkyl jacobtot.gkyl interp -b ms -p1 pl --title "jacobtot"&
# pgkyl jacobtot_inv.gkyl interp -b ms -p1 pl --title "jacobtot_inv"&
# pgkyl jacogeo_inv.gkyl interp -b ms -p1 pl --title "jacobgeo_inv"&
# pgkyl b_i.gkyl interp -b ms -p1 pl --title "b_i"&
# pgkyl mapc2p.gkyl interp -b ms -p1 select --z0 0.5 --z1 0 pl --title "mapc2p"&
# pgkyl bmag.gkyl interp -b ms -p1 pl --title "bmag"&
# pgkyl bmag_inv.gkyl interp -b ms -p1 pl --title "bmag_inv"&
# pgkyl bmag_inv_sq.gkyl interp -b ms -p1 pl --title "bmag_inv_sq"&
# pgkyl cmag.gkyl interp -b ms -p1 pl --title "cmag"&

# Distribution function at a single velocity space point at all z
# "$name-"$species"_[0-9]*.gkyl"\
#   interp -b gkhyb -p1 sel --z1 5e6 --z2 3.5e-15 ev 'f[:] abs' collect pl --logz --zmin 1e-4 \

#Plot phi
frame=117
# pgkyl "$name-field_[0-9]*.gkyl" interp -b ms -p1 collect sel --z1 0 pl --title 'phi' &
# pgkyl "$name-field_[0-9]*.gkyl" interp -b ms -p1 collect pl --title 'phi' &
# pgkyl "$name-field_$frame.gkyl" interp -b ms -p1 ev 'f[:] 940 /' pl --title "phi at frame $frame in unite e phi/Te" &

# Moments of the distribution function
frame=3
if [ "$species" = "elc" ]; then mass=9.11e-31
elif [ "$species" = "ion" ]; then mass=3.34e-27
else echo "species must be ion or elc"
fi
# pgkyl $name-"$species"_M0_$frame.gkyl $name-"$species"_M1_$frame.gkyl interp -b ms -p1 ev "f[1] f[0] /" pl --title 'Upar' &
# pgkyl $name-"$species"_prim_moms_$frame.gkyl interp -b ms -p1 pl --title 'prim--moms' &
# pgkyl $name-"$species"_M0_$frame.gkyl $name-"$species"_M2perp_$frame.gkyl interp -b ms -p1 ev "$mass f[1] f[0] / * 1.6e-19 /" pl --title 'Tperp (eV)' &
# pgkyl $name-"$species"_M0_$frame.gkyl $name-"$species"_M1_$frame.gkyl $name-"$species"_M2par_$frame.gkyl interp -b ms -p1 ev "f[2] f[1] f[1] * f[0] / - $mass * f[0] / 1.6e-19 /" pl --title 'Tpar (eV)' &


## compute moments of distribution function
# for frame in $(seq 0 200); do
#   pgkyl $name-"$species"_M0_$frame.gkyl interp -b ms -p1 writ -f "$name-"$species"-dens_$frame.gkyl" &
#   pgkyl $name-"$species"_M0_$frame.gkyl $name-"$species"_M1_$frame.gkyl interp -b ms -p1 ev "f[1] f[0] /"  writ -f "$name-"$species"-upar_$frame.gkyl" &
#   pgkyl $name-"$species"_M0_$frame.gkyl $name-"$species"_M2perp_$frame.gkyl interp -b ms -p1 ev "$mass f[1] f[0] / * 1.6e-19 /" writ -f "$name-"$species"-tperp_$frame.gkyl" &
#   pgkyl $name-"$species"_M0_$frame.gkyl $name-"$species"_M1_$frame.gkyl $name-"$species"_M2par_$frame.gkyl interp -b ms -p1 ev "f[2] f[1] f[1] * f[0] / - $mass * f[0] / 1.6e-19 /" writ -f "$name-"$species"-tpar_$frame.gkyl" &
# done

# Plot moments of two distribution functions
frame1=0
frame2=200
# pgkyl $name-"$species"_dens_$frame1.gkyl $name-"$species"_dens_$frame2.gkyl pl --title 'Density' -f0 &
# pgkyl $name-"$species"_upar_$frame1.gkyl $name-"$species"_upar_$frame2.gkyl pl --title 'Upar' -f0 &
# pgkyl $name-"$species"_tperp_$frame1.gkyl $name-"$species"_tperp_$frame2.gkyl pl --title 'Tperp (eV)' -f0 &
# pgkyl $name-"$species"_tpar_$frame1.gkyl $name-"$species"_tpar_$frame2.gkyl pl --title 'Tpar (eV)' -f0 &

# Plot moments at all frames
# pgkyl "$name-$species-tperp_[0-9]*.gkyl" collect &
# pgkyl $name-"$species"_tpar_[0-9]*.gkyl collect pl --title 'Tpar (eV)' &

frame=[0-9]*
pgkyl $name-"$species"_M0_$frame.gkyl $name-"$species"_M2perp_$frame.gkyl\
 interp -b ms -p1 ev "$mass f[1] f[0] / * 1.6e-19 /" collect pl &

