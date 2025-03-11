
# name1="outputs/gk_mirror_adiabatic_elc_1x2v_p1_nosource_uniform"
# name2="outputs/gk_mirror_adiabatic_elc_1x2v_p1_nosource_nonuniform"
name="gk_wham"
# name="outputs/gk_mirror_adiabatic_elc_1x2v_p1_nosource_nonuniform"
species="elc"
saveLoc="python-plots/Distribution-functions-figures-for-movies/$name"
saveMovieLoc="python-plots/Distribution-functions-movies/$name"
hiresLoc="../stellar-wham1x-kinetic-compare-collision-frequencies"

pgkyl gk_wham-ion_0.gkyl interp sel --z0 288 pl --title "initial" &
pgkyl initial-condition/gk_wham-ion_26.gkyl interp sel --z0 288 pl --title "donor" &
pgkyl gk_wham-ion_M0_0.gkyl interp pl --title "M0 initial" &
pgkyl initial-condition/gk_wham-ion_M0_26.gkyl interp pl --title "M0 donor" &
# pgkyl gk_wham-ion_M0_0.gkyl interp pl --title "M0" --saveas "M0.png" --no-show --title "128 cells"&
# pgkyl "64-cells/gk_wham-ion_M0_300.gkyl" interp pl --title "M0" --saveas "M0-300.png" --no-show --title "64 cells"&
# # Negativity of the T_par
# pgkyl "BiMaxwellianMoments/"$name"-"$species"_BiMaxwellianMoments_[0-9]*.gkyl" interp anim --float &
# pgkyl "Field/$name-field_[0-9]*.gkyl" interp anim --float &

# pgkyl "BiMaxwellianMoments/"$name"-"$species"_BiMaxwellianMoments_[0-9]*.gkyl" interp sel --z0 0.0 -c0 collect pl --logy --xlabel "Time, s" --ylabel "$species density, m^{-3}" --title "Midplane $species density 32 cells" &
# pgkyl "$hiresLoc/BiMaxwellianMoments/gk_wham_modified-"$species"_BiMaxwellianMoments_[0-9]*.gkyl" interp sel --z0 0.0 -c0 collect pl --logy --xlabel "Time, s" --ylabel ""$species" density, m^{-3}" --title "Midplane "$species" density 288 cells" &

# pgkyl "Field/"$name"-field_[0-9]*.gkyl" interp sel --z0 0.0 -c0 collect pl --xlabel "Time, s" --ylabel "Potential, V" --title "Midplane potential 32 cells" &
# pgkyl "$hiresLoc/Field/gk_wham_modified-field_[0-9]*.gkyl" interp sel --z0 0.0 -c0 collect pl --xlabel "Time, s" --ylabel "Potential, V" --title "Midplane potential 288 cells" &

# pgkyl --c2p-vel "$name"-"$species"_mapc2p_vel.gkyl "Distributions/"$name"-"$species"_[0-9]*.gkyl" interp sel --z0 7 anim --logz --zmin 1e-20 &
# pgkyl --c2p-vel "$name"-"$species"_mapc2p_vel.gkyl "Distributions/"$name"-"$species"_[0-9]*.gkyl" interp sel --z0 7 anim --float &

# pgkyl --c2p-vel gk_wham-elc_mapc2p_vel.gkyl Distributions/gk_wham-elc_16.gkyl interp sel --z0 7 --z2 0.2e-15 pl

# pgkyl --c2p-vel "$hiresLoc/gk_wham_modified-"$species"_mapc2p_vel.gkyl" "$hiresLoc/Distributions/gk_wham_modified-"$species"_[0-9]*.gkyl" interp sel --z0 2 anim --logz --zmin 1e-20 &
# pgkyl --c2p-vel "$hiresLoc/gk_wham_modified-"$species"_mapc2p_vel.gkyl" "$hiresLoc/Distributions/gk_wham_modified-"$species"_[0-9]*.gkyl" interp sel --z0 2 anim --float &

# frame=16
# pgkyl --c2p-vel "$hiresLoc/gk_wham_modified-"$species"_mapc2p_vel.gkyl" "$hiresLoc/Distributions/gk_wham_modified-"$species"_"$frame".gkyl" interp sel --z2 0. --z0 -2. pl --title "288 cells $species" --logz --zmin 1e-20 &
# pgkyl --c2p-vel "$name"-"$species"_mapc2p_vel.gkyl Distributions/"$name"-"$species"_"$frame".gkyl interp sel --z2 0. --z0 -2. pl --title "32 cells $species" --logz --zmin 1e-20 &

# pgkyl --c2p-vel "$hiresLoc/gk_wham_modified-"$species"_mapc2p_vel.gkyl" "$hiresLoc/Distributions/gk_wham_modified-"$species"_"$frame".gkyl" interp sel --z0 0.0 pl --logz --zmin 1e-20 --title "288 cells midplane $species"&
# pgkyl --c2p-vel "$name"-ion_mapc2p_vel.gkyl Distributions/"$name"-"$species"_"$frame".gkyl interp sel --z0 0.0 pl --title "32 cells midplane $species" --logz --zmin 1e-20 &

# pgkyl BiMaxwellianMoments/"$name"-ion_BiMaxwellianMoments_62.gkyl interp pl --title "Initial" &
# pgkyl "$name"-ion_BiMaxwellianMoments_62.gkyl interp pl --title "Restart" &
# pgkyl "$name"-ion_BiMaxwellianMoments_63.gkyl interp pl --title "63" &
# pgkyl "$name"-ion_BiMaxwellianMoments_64.gkyl interp pl --title "64" &



# pgkyl --c2p-vel "$name"-ion_mapc2p_vel.gkyl Distributions/"$name"-ion_62.gkyl interp sel --z0 2 pl --title "Initial" &
# pgkyl --c2p-vel "$name"-ion_mapc2p_vel.gkyl Distributions/"$name"-ion_63.gkyl interp sel --z0 2 pl --title "63" &
# pgkyl --c2p-vel "$name"-ion_mapc2p_vel.gkyl Distributions/"$name"-ion_64.gkyl interp sel --z0 2 pl --title "64" &
# pgkyl --c2p-vel "$name"-ion_mapc2p_vel.gkyl Distributions/"$name"-ion_68.gkyl interp sel --z0 2 pl --title "68" --logz --zmin 1e-20 &

# pgkyl --c2p-vel "$name"-"$species"_mapc2p_vel.gkyl Distributions/"$name"-"$species"_150.gkyl interp sel --z0 0.0 pl --title "150 microseconds 32 cells" --logz --zmin 1e-20 &
# pgkyl --c2p-vel $hiresLoc/gk_wham_modified-"$species"_mapc2p_vel.gkyl $hiresLoc/Distributions/gk_wham_modified-"$species"_150.gkyl interp sel --z0 0.0 pl --title "150 microseconds 288 cells" --logz --zmin 1e-20 &

# pgkyl Distributions/"$name"-ion_62.gkyl interp sel --z0 0.0 pl --title "Initial" &
# pgkyl "$name"-ion_62.gkyl interp sel --z0 0.0 pl --title "Restart" &



# # Early frame when T_par < 0
# pgkyl BiMaxwellianMoments/gk_wham_64_npos-ion_BiMaxwellianMoments_10.gkyl interp pl

# # f(v)  
# pgkyl gk_wham_64_npos-ion_10.gkyl interp sel --z0 1.86 pl --xlabel "\$v_\parallel\$" --ylabel "\$\mu\$" -d
# pgkyl gk_wham_64_npos-ion_10.gkyl interp sel --z0 1.86 --z1 0.133 pl

# # n*upar
# pgkyl BiMaxwellianMoments/gk_wham_64_npos-ion_BiMaxwellianMoments_10.gkyl interp ev "f[0][0] f[0][1] *" pl --ylabel "\$n u_\parallel\$"


# # Late frame when T_par < 0
# pgkyl BiMaxwellianMoments/gk_wham_64_npos-ion_BiMaxwellianMoments_100.gkyl interp pl   

# # f(v)
# pgkyl gk_wham_64_npos-ion_100.gkyl interp sel --z0 1.86 pl --xlabel "\$v_\parallel\$" --ylabel "\$\mu\$" -d

# pgkyl outputs/gk_mirror_adiabatic_elc_1x2v_p1_true_maxwellian-ion_0.gkyl --c2p outputs/mapc2p.gkyl \
#  interp -b gkhyb -p1 integrate 1 pl &
# pgkyl outputs/gk_mirror_adiabatic_elc_1x2v_p1_true_maxwellian-ion_0.gkyl \
#  interp -b gkhyb -p1 integrate 1 pl &


# pgkyl "Geometry/$name-"b_i.gkyl interp -b ms -p1 pl --title "b_i" --saveas "$saveLoc-geo-b_i.png" --no-show &
# pgkyl "Geometry/$name-"bmag_inv_sq.gkyl interp -b ms -p1 pl --title "bmag_inv_sq" --saveas "$saveLoc-geo-bmag_inv_sq.png" --no-show &
# pgkyl "Geometry/$name-"bmag_inv.gkyl interp -b ms -p1 pl --title "bmag_inv" --saveas "$saveLoc-geo-bmag_inv.png" --no-show &
# pgkyl "Geometry/$name-"bmag.gkyl interp -b ms -p1 pl --title "bmag" --saveas "$saveLoc-geo-bmag.png" --no-show &
# pgkyl "Geometry/$name-"cmag.gkyl interp -b ms -p1 pl --title "cmag" --saveas "$saveLoc-geo-cmag.png" --no-show &
# pgkyl "Geometry/$name-"jacobtot_inv.gkyl interp -b ms -p1 pl --title "jacobtot_inv" --saveas "$saveLoc-geo-jacobtot_inv.png" --no-show &
# pgkyl "Geometry/$name-"jacobgeo.gkyl interp -b ms -p1 pl --title "jacobgeo" --saveas "$saveLoc-geo-jacobgeo.png" --no-show &
# pgkyl "Geometry/$name-"jacobtot.gkyl interp -b ms -p1 pl --title "jacobtot" --saveas "$saveLoc-geo-jacobtot.png" --no-show &
# pgkyl "Geometry/$name-"mapc2p.gkyl interp -b ms -p1 pl --title "mapc2p" --saveas "$saveLoc-geo-mapc2p.png" --no-show &

frame=98
# pgkyl "BiMaxwellianMoments/$name-"$species"_BiMaxwellianMoments_[0-9]*.gkyl" interp -b ms -p1 anim &
# for frame in {0..100}
# do
#   pgkyl "BiMaxwellianMoments/$name-"$species"_BiMaxwellianMoments_0.gkyl" interp -b ms -p1 pl --saveas "$saveLoc-$species-BiMaxwellianMoments-$frame.png" --no-show &
# done

# for frame in {0..98}
# do
#   framepad=$(printf "%03d" $frame)
#   echo "frame $framepad"
#   pgkyl --c2p-vel $name-"$species"_mapc2p_vel.gkyl $name-"$species"_$frame.gkyl -t jf "$name-"$species"_jacobvel.gkyl" -t jac ev -t df 'jf jac /' activ -t df interp -b gkhyb -p1 integrate 2 pl --title "Frame $frame vpar" --saveas "$saveLoc-$species-$framepad-vpar.png" --no-show --xlabel "z" --ylabel "$ v_{||}$" &
#   pgkyl --c2p-vel $name-"$species"_mapc2p_vel.gkyl $name-"$species"_$frame.gkyl -t jf "$name-"$species"_jacobvel.gkyl" -t jac ev -t df 'jf jac /' activ -t df interp -b gkhyb -p1 integrate 1 pl --title "Frame $frame mu" --saveas "$saveLoc-$species-$framepad-mu.png" --no-show --xlabel "z" --ylabel "$\mu$" &
#   pgkyl --c2p-vel $name-"$species"_mapc2p_vel.gkyl $name-"$species"_$frame.gkyl -t jf "$name-"$species"_jacobvel.gkyl" -t jac ev -t df 'jf jac /' activ -t df interp -b gkhyb -p1 sel --z0 0.0 pl --title "Frame $frame z=0" --saveas "$saveLoc-$species-$framepad-z0.png" --no-show --xlabel "$ v_{||}$" --ylabel "$\mu$" &
#   pgkyl --c2p-vel $name-"$species"_mapc2p_vel.gkyl $name-"$species"_$frame.gkyl -t jf "$name-"$species"_jacobvel.gkyl" -t jac ev -t df 'jf jac /' activ -t df interp -b gkhyb -p1 sel --z0 0.70 pl --title "Frame $frame z=0.7" --saveas "$saveLoc-$species-$framepad-z0.7.png" --no-show --xlabel "$ v_{||}$" --ylabel "$\mu$" &
#   pgkyl --c2p-vel $name-"$species"_mapc2p_vel.gkyl $name-"$species"_$frame.gkyl -t jf "$name-"$species"_jacobvel.gkyl" -t jac ev -t df 'jf jac /' activ -t df interp -b gkhyb -p1 sel --z0 0.98 pl --title "Frame $frame z=0.98" --saveas "$saveLoc-$species-$framepad-z0.98.png" --no-show --xlabel "$ v_{||}$" --ylabel "$\mu$" &
#   pgkyl --c2p-vel $name-"$species"_mapc2p_vel.gkyl $name-"$species"_$frame.gkyl -t jf "$name-"$species"_jacobvel.gkyl" -t jac ev -t df 'jf jac /' activ -t df interp -b gkhyb -p1 sel --z0 2.00 pl --title "Frame $frame z=2.0" --saveas "$saveLoc-$species-$framepad-z2.0.png" --no-show --xlabel "$ v_{||}$" --ylabel "$\mu$" &

#   pgkyl --c2p-vel $name-"$species"_mapc2p_vel.gkyl $name-"$species"_$frame.gkyl -t jf "$name-"$species"_jacobvel.gkyl" -t jac ev -t df 'jf jac /' activ -t df interp -b gkhyb -p1 integrate 2 pl --logz --zmin 1e-30 --title "Frame $frame vpar" --saveas "$saveLoc-$species-$framepad-vpar-logz.png" --no-show --xlabel "z" --ylabel "$ v_{||}$" &
#   pgkyl --c2p-vel $name-"$species"_mapc2p_vel.gkyl $name-"$species"_$frame.gkyl -t jf "$name-"$species"_jacobvel.gkyl" -t jac ev -t df 'jf jac /' activ -t df interp -b gkhyb -p1 integrate 1 pl --logz --zmin 1e-30 --title "Frame $frame mu" --saveas "$saveLoc-$species-$framepad-mu-logz.png" --no-show --xlabel "z" --ylabel "$\mu$" &
#   pgkyl --c2p-vel $name-"$species"_mapc2p_vel.gkyl $name-"$species"_$frame.gkyl -t jf "$name-"$species"_jacobvel.gkyl" -t jac ev -t df 'jf jac /' activ -t df interp -b gkhyb -p1 sel --z0 0.0 pl --logz --zmin 1e-30 --title "Frame $frame z=0" --saveas "$saveLoc-$species-$framepad-z0-logz.png" --no-show --xlabel "$ v_{||}$" --ylabel "$\mu$" &
#   pgkyl --c2p-vel $name-"$species"_mapc2p_vel.gkyl $name-"$species"_$frame.gkyl -t jf "$name-"$species"_jacobvel.gkyl" -t jac ev -t df 'jf jac /' activ -t df interp -b gkhyb -p1 sel --z0 0.70 pl --logz --zmin 1e-30 --title "Frame $frame z=0.7" --saveas "$saveLoc-$species-$framepad-z0.7-logz.png" --no-show --xlabel "$ v_{||}$" --ylabel "$\mu$" &
#   pgkyl --c2p-vel $name-"$species"_mapc2p_vel.gkyl $name-"$species"_$frame.gkyl -t jf "$name-"$species"_jacobvel.gkyl" -t jac ev -t df 'jf jac /' activ -t df interp -b gkhyb -p1 sel --z0 0.98 pl --logz --zmin 1e-30 --title "Frame $frame z=0.98" --saveas "$saveLoc-$species-$framepad-z0.98-logz.png" --no-show --xlabel "$ v_{||}$" --ylabel "$\mu$" &
#   pgkyl --c2p-vel $name-"$species"_mapc2p_vel.gkyl $name-"$species"_$frame.gkyl -t jf "$name-"$species"_jacobvel.gkyl" -t jac ev -t df 'jf jac /' activ -t df interp -b gkhyb -p1 sel --z0 2.00 pl --logz --zmin 1e-30 --title "Frame $frame z=2.0" --saveas "$saveLoc-$species-$framepad-z2.0-logz.png" --no-show --xlabel "$ v_{||}$" --ylabel "$\mu$" &
#   sleep 3
# done
# bash make-movie.sh -r 10 -o "$saveMovieLoc-$species-vpar.mp4" -i "$saveLoc-$species"-%03d-vpar.png
# bash make-movie.sh -r 10 -o "$saveMovieLoc-$species-vpar-logz.mp4" -i "$saveLoc-$species"-%03d-vpar-logz.png
# bash make-movie.sh -r 10 -o "$saveMovieLoc-$species-mu.mp4" -i "$saveLoc-$species"-%03d-mu.png
# bash make-movie.sh -r 10 -o "$saveMovieLoc-$species-mu-logz.mp4" -i "$saveLoc-$species"-%03d-mu-logz.png
# bash make-movie.sh -r 10 -o "$saveMovieLoc-$species-z0.mp4" -i "$saveLoc-$species"-%03d-z0.png
# bash make-movie.sh -r 10 -o "$saveMovieLoc-$species-z0-logz.mp4" -i "$saveLoc-$species"-%03d-z0-logz.png
# bash make-movie.sh -r 10 -o "$saveMovieLoc-$species-z0.7.mp4" -i "$saveLoc-$species"-%03d-z0.7.png
# bash make-movie.sh -r 10 -o "$saveMovieLoc-$species-z0.7-logz.mp4" -i "$saveLoc-$species"-%03d-z0.7-logz.png
# bash make-movie.sh -r 10 -o "$saveMovieLoc-$species-z0.98.mp4" -i "$saveLoc-$species"-%03d-z0.98.png
# bash make-movie.sh -r 10 -o "$saveMovieLoc-$species-z0.98-logz.mp4" -i "$saveLoc-$species"-%03d-z0.98-logz.png
# bash make-movie.sh -r 10 -o "$saveMovieLoc-$species-z2.0.mp4" -i "$saveLoc-$species"-%03d-z2.0.png
# bash make-movie.sh -r 10 -o "$saveMovieLoc-$species-z2.0-logz.mp4" -i "$saveLoc-$species"-%03d-z2.0-logz.png

# pgkyl --c2p-vel $name-"$species"_mapc2p_vel.gkyl $name-"$species"_[0-9]*.gkyl interp -b gkhyb -p1 integrate 2 anim --logz --zmin 1e-30 --saveas "$saveLoc-$species-vpar-logz.mp4" --no-show & 
# pgkyl --c2p-vel $name-"$species"_mapc2p_vel.gkyl $name-"$species"_[0-9]*.gkyl interp -b gkhyb -p1 integrate 1 anim --saveas "$saveLoc-$species-mu.mp4" --no-show &
# pgkyl --c2p-vel $name-"$species"_mapc2p_vel.gkyl $name-"$species"_[0-9]*.gkyl interp -b gkhyb -p1 integrate 1 anim --logz --zmin 1e-30 --saveas "$saveLoc-$species-mu-logz.mp4" --no-show &
# pgkyl --c2p-vel $name-"$species"_mapc2p_vel.gkyl $name-"$species"_[0-9]*.gkyl interp -b gkhyb -p1 sel --z0 0.0 anim --saveas "$saveLoc-$species-z0.mp4" --no-show &
# pgkyl --c2p-vel $name-"$species"_mapc2p_vel.gkyl $name-"$species"_[0-9]*.gkyl interp -b gkhyb -p1 sel --z0 0.0 anim --logz --zmin 1e-30 --saveas "$saveLoc-$species-z0-logz.mp4" --no-show &
# pgkyl --c2p-vel $name-"$species"_mapc2p_vel.gkyl $name-"$species"_[0-9]*.gkyl interp -b gkhyb -p1 sel --z0 0.70 anim --saveas "$saveLoc-$species-z0.7.mp4" --no-show &
# pgkyl --c2p-vel $name-"$species"_mapc2p_vel.gkyl $name-"$species"_[0-9]*.gkyl interp -b gkhyb -p1 sel --z0 0.70 anim --logz --zmin 1e-30 --saveas "$saveLoc-$species-z0.7-logz.mp4" --no-show &
# pgkyl --c2p-vel $name-"$species"_mapc2p_vel.gkyl $name-"$species"_[0-9]*.gkyl interp -b gkhyb -p1 sel --z0 0.98 anim --saveas "$saveLoc-$species-z0.98.mp4" --no-show &
# pgkyl --c2p-vel $name-"$species"_mapc2p_vel.gkyl $name-"$species"_[0-9]*.gkyl interp -b gkhyb -p1 sel --z0 0.98 anim --logz --zmin 1e-30 --saveas "$saveLoc-$species-z0.98-logz.mp4" --no-show &
# pgkyl --c2p-vel $name-"$species"_mapc2p_vel.gkyl $name-"$species"_[0-9]*.gkyl interp -b gkhyb -p1 sel --z0 2.00 anim --saveas "$saveLoc-$species-z2.0.mp4" --no-show &  
# pgkyl --c2p-vel $name-"$species"_mapc2p_vel.gkyl $name-"$species"_[0-9]*.gkyl interp -b gkhyb -p1 sel --z0 2.00 anim --logz --zmin 1e-30 --saveas "$saveLoc-$species-z2.0-logz.mp4" --no-show &


# pgkyl "$name-ion_[0-9]*.gkyl" interp -b gkhyb -p1 sel --z0 0.0 anim --title "ion distribution function" --xlabel "vpar" --ylabel "mu" --logz --zmin 1e-20 --fps 4 &
# pgkyl "$name-ion_[0-9]*.gkyl" interp -b gkhyb -p1 sel --z0 0.98 anim --title "ion distribution function" --xlabel "vpar" --ylabel "mu" --logz --zmin 1e-20 --fps 4 &
# pgkyl "$name-ion_[0-9]*.gkyl" interp -b gkhyb -p1 sel --z0 0.98 anim --title "ion distribution function" --xlabel "vpar" --ylabel "mu" --zmax 1e-8 --fps 4 &
# pgkyl --c2p-vel "$name-ion_mapc2p_vel.gkyl" "$name-ion_[0-9]*.gkyl" interp -b gkhyb -p1 sel --z0 2.00 anim --xlabel "vpar" --ylabel "mu, z=2.0, $name" --logz --zmin 1e-20 --fps 4 --saveas "$saveLoc-ion-z2.0_logz.mp4" --no-show &
# pgkyl --c2p-vel "$name-ion_mapc2p_vel.gkyl" "$name-ion_[0-9]*.gkyl" interp -b gkhyb -p1 sel --z0 0.98 anim --xlabel "vpar" --ylabel "mu, z=0.98, $name" --logz --zmin 1e-20 --fps 4 --saveas "$saveLoc-ion-z0.98_logz.mp4" --no-show &
# pgkyl --c2p-vel "$name-ion_mapc2p_vel.gkyl" "$name-ion_[0-9]*.gkyl" interp -b gkhyb -p1 sel --z0 0.00 anim --xlabel "vpar" --ylabel "mu, z=0, $name" --logz --zmin 1e-20 --fps 4 --saveas "$saveLoc-ion-z0.0_logz.mp4" --no-show &


# pgkyl --c2p-vel "$name-ion_mapc2p_vel.gkyl" "$name-ion_[0-9]*.gkyl" interp -b gkhyb -p1 sel --z0 2.00 anim --xlabel "vpar" --ylabel "mu, z=2.0, $name" --fps 4 --saveas "$saveLoc-ion-z2.0.mp4" --no-show &
# pgkyl --c2p-vel "$name-ion_mapc2p_vel.gkyl" "$name-ion_[0-9]*.gkyl" interp -b gkhyb -p1 sel --z0 0.98 anim --xlabel "vpar" --ylabel "mu, z=0.98, $name" --fps 4 --saveas "$saveLoc-ion-z0.98.mp4" --no-show &
# pgkyl --c2p-vel "$name-ion_mapc2p_vel.gkyl" "$name-ion_[0-9]*.gkyl" interp -b gkhyb -p1 sel --z0 0.00 anim --xlabel "vpar" --ylabel "mu, z=0, $name" --fps 4 --saveas "$saveLoc-ion-z0.0.mp4" --no-show &


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
# frame=32


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
# frame=117
# pgkyl "$name-field_[0-9]*.gkyl" interp -b ms -p1 collect sel --z1 0 pl --title 'phi' &
# pgkyl "$name-field_[0-9]*.gkyl" interp -b ms -p1 collect pl --title 'phi' &
# pgkyl "$name-field_$frame.gkyl" interp -b ms -p1 ev 'f[:] 940 /' pl --title "phi at frame $frame in unite e phi/Te" &


# if [ "$species" = "elc" ]; then mass=9.11e-31
# elif [ "$species" = "ion" ]; then mass=3.34e-27
# else echo "species must be ion or elc"
# fi
# Moments of the distribution function at a single frame
# frame=32
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
# frame=32
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

