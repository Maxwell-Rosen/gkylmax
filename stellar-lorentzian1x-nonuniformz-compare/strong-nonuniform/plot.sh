mkdir -p plots

# frame=$(ls gk_lorentzian_mirror-ion_[0-9]*.gkyl | sed -E 's/.*_([0-9]+)\.gkyl/\1/' | sort -n | tail -1)
frame=0

pgkyl --c2p-vel gk_lorentzian_mirror-ion_mapc2p_vel.gkyl gk_lorentzian_mirror-ion_${frame}.gkyl -t JJf gk_lorentzian_mirror-jacobgeo_inv.gkyl -t jacobgeo_inv gk_lorentzian_mirror-ion_jacobvel.gkyl -t jacobvel \
	ev -t Jf "JJf jacobvel /" \
	activ -t Jf interp -b gkhyb -p1 -t Jfint \
	activ -t jacobgeo_inv interp -t jacobgeoint \
	activ -t jacobgeoint,Jfint ev "Jfint jacobgeoint *" integ 2 pl --title "Ion f integrated over mu" --xlabel "z computational" --ylabel "$ v_\parallel $" --saveas plots/ion_f_25.pdf

# frame=$(ls gk_lorentzian_mirror-ion_BiMaxwellianMoments_*.gkyl | sed -E 's/.*_([0-9]+)\.gkyl/\1/' | sort -n | tail -1)

# pgkyl --c2p gk_lorentzian_mirror-mc2nu_pos_deflated.gkyl gk_lorentzian_mirror-ion_BiMaxwellianMoments_${frame}.gkyl interp pl -f0 --title "Ion BiMaxwellian Moments" --xlabel "z cylindrical" --saveas plots/ion_BiMaxMoments.pdf --no-show &

# pgkyl --c2p gk_lorentzian_mirror-mc2nu_pos_deflated.gkyl gk_lorentzian_mirror-bmag.gkyl interp pl -f0 --title "Magnetic Field Strength" --xlabel "z cylindrical" --saveas plots/bmag.pdf --no-show &
# pgkyl gk_lorentzian_mirror-bmag.gkyl interp pl -f0 --title "Magnetic Field Strength" --xlabel "z computational" --saveas plots/bmag_comp.pdf --no-show &

# #Jacobian
# pgkyl --c2p gk_lorentzian_mirror-mc2nu_pos_deflated.gkyl gk_lorentzian_mirror-jacobgeo.gkyl interp pl -f0 --title "Jacobian" --xlabel "z cylindrical" --saveas plots/jacobian.pdf --no-show &
# pgkyl gk_lorentzian_mirror-jacobgeo.gkyl interp pl -f0 --title "Jacobian" --xlabel "z computational" --saveas plots/jacobian_comp.pdf --no-show &

# # Jacobian inverse
# pgkyl --c2p gk_lorentzian_mirror-mc2nu_pos_deflated.gkyl gk_lorentzian_mirror-jacobgeo_inv.gkyl interp pl -f0 --title "Jacobian Inverse" --xlabel "z cylindrical" --saveas plots/jacobian_inv.pdf --no-show &
# pgkyl gk_lorentzian_mirror-jacobgeo_inv.gkyl interp pl -f0 --title "Jacobian Inverse" --xlabel "z computational" --saveas plots/jacobian_inv_comp.pdf --no-show &

# pgkyl gk_lorentzian_mirror-mc2nu_pos_deflated.gkyl interp ev 'f grad' pl --title "Gradient of non-uniform map Mapping" --xlabel "z computational" --ylabel "z physical" --saveas plots/grad_map.pdf --no-show &
# pgkyl gk_lorentzian_mirror-mc2nu_pos_deflated.gkyl interp ev '1 f grad /' pl --title "Inverse gradient of non-uniform map Mapping" --xlabel "z computational" --ylabel "z physical" --saveas plots/grad_map_inv.pdf --no-show &