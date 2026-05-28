# pgkyl with-nu-ie/gk_lorentzian_mirror-ion_integrated_moms.gkyl  without-nu-ie/gk_lorentzian_mirror-ion_integrated_moms.gkyl pl -f0 --legend "with nuie, without nuie" -x "t" --subplot-ylabels "Integrated M0,Integrated M1,Integrated M2par,Integrated M2perp" -y "" --saveas "integ_compare.png" --title "Integrated moments comparison" &

# pgkyl with-nu-ie/gk_lorentzian_mirror-elc_integrated_moms.gkyl pl -f0 -x "t" --subplot-ylabels "Integrated M0,Integrated M1,Integrated M2" -y "" --saveas "elc_integ.png" --title "Electron integrated moments" &

# pgkyl --c2p with-nu-ie/gk_lorentzian_mirror-mapc2p_deflated.gkyl with-nu-ie/gk_lorentzian_mirror-ion_BiMaxwellianMoments_0.gkyl with-nu-ie/gk_lorentzian_mirror-ion_BiMaxwellianMoments_18.gkyl without-nu-ie/gk_lorentzian_mirror-ion_BiMaxwellianMoments_29.gkyl interp sel --z0 0 pl -f0 --legend "initial,with nuie,without nuie" -x "t" --subplot-ylabels "Density,Upar,Tpar,Tperp" -y "" --saveas "bimax_compare.png" --title "Bimaxwellian moments comparison" &

# pgkyl --c2p-vel with-nu-ie/gk_lorentzian_mirror-ion_mapc2p_vel.gkyl with-nu-ie/gk_lorentzian_mirror-ion_18.gkyl interp sel --z0 0 --z1 0.0 pl -x "vpar" -y "vperp" --saveas "distf-with-nu-ie.png" --title "Final distf with nuie" &

# pgkyl --c2p-vel without-nu-ie/gk_lorentzian_mirror-ion_mapc2p_vel.gkyl without-nu-ie/gk_lorentzian_mirror-ion_29.gkyl interp sel --z0 0 --z1 0.0 pl -x "vpar" -y "vperp" --saveas "distf-without-nu-ie.png" --title "Final distf without nuie" &

# pgkyl --c2p with-nu-ie/gk_lorentzian_mirror-mapc2p_deflated.gkyl with-nu-ie/gk_lorentzian_mirror-elc_MaxwellianMoments_18.gkyl interp sel --z0 0 pl -x "Z" -y "" --subplot-ylabels "Density,Upar,T" --saveas "elc_max_compare.png" --title "Electron Maxwellian moments comparison" &

# pgkyl --c2p with-nu-ie/gk_lorentzian_mirror-mapc2p_deflated.gkyl with-nu-ie/gk_lorentzian_mirror-elc_lbo_nu_sum_18.gkyl interp sel --z0 0 pl --saveas "elc_lbo_sum.png" --title "Electron LBO sum" &

# pgkyl --c2p with-nu-ie/gk_lorentzian_mirror-mapc2p_deflated.gkyl with-nu-ie/gk_lorentzian_mirror-ion_lbo_nu_sum_18.gkyl interp sel --z0 0 pl --saveas "ion_lbo_sum.png" --title "Ion LBO sum" &

# pgkyl --c2p without-nu-ie/gk_lorentzian_mirror-mapc2p_deflated.gkyl without-nu-ie/gk_lorentzian_mirror-ion_lbo_nu_sum_18.gkyl interp sel --z0 0 pl --saveas "ion_lbo_sum_without_nuie.png" --title "Ion LBO sum without nuie" &

pgkyl --c2p with-nu-ie/gk_lorentzian_mirror-mapc2p_deflated.gkyl with-nu-ie/gk_lorentzian_mirror-ion_lbo_nu_prim_moms_1.gkyl interp sel --z0 0 pl -x "Z" --saveas "ion_lbo_prim_compare.png" --title "Ion LBO primitive moments with nuie" &

pgkyl --c2p without-nu-ie/gk_lorentzian_mirror-mapc2p_deflated.gkyl without-nu-ie/gk_lorentzian_mirror-ion_lbo_nu_prim_moms_1.gkyl interp sel --z0 0 pl -x "Z" --saveas "ion_lbo_prim_compare_without_nuie.png" --title "Ion LBO primitive moments without nuie" &