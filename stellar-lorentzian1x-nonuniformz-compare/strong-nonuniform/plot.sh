mkdir -p plots

frame=$(ls gk_lorentzian_mirror-ion_*.gkyl | sed -E 's/.*_([0-9]+)\.gkyl/\1/' | sort -n | tail -1)

pgkyl --c2p-vel gk_lorentzian_mirror-ion_mapc2p_vel.gkyl gk_lorentzian_mirror-ion_${frame}.gkyl gk_lorentzian_mirror-jacobgeo_inv.gkyl interp ev "f[0] f[1] *" integ 2 pl --title "Ion f integrated over mu" --xlabel "z cylindrical" --ylabel "$v_\parallel$" --saveas plots/ion_f_25.pdf &

frame=$(ls gk_lorentzian_mirror-ion_BiMaxwellianMoments_*.gkyl | sed -E 's/.*_([0-9]+)\.gkyl/\1/' | sort -n | tail -1)


pgkyl --c2p gk_lorentzian_mirror-mc2nu_pos_deflated.gkyl gk_lorentzian_mirror-ion_BiMaxwellianMoments_${frame}.gkyl interp pl -f0 --title "Ion BiMaxwellian Moments" --xlabel "z cylindrical" --saveas plots/ion_BiMaxMoments.pdf &