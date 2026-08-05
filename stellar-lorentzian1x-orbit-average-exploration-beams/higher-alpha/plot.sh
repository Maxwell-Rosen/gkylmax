pgkyl "gk_lorentzian_mirror-ion_BiMaxwellianMoments_65.gkyl" \
      "../../stellar-lorentzian1x-orbit-average-beams/gk_lorentzian_mirror-ion_BiMaxwellianMoments_65.gkyl"\
      interp ev 'f 2,3 1.67e-27 2.014 * 1.6e-19 / scale_comp' \
      pl --xlabel 'z [m]' -f3 --saveas plots/high-vs-low-alpha.png --figsize "12,10" --no-show \
      --subplot-ylabels 'Density $ m^3 $ , $ U_\parallel m/s $, $ T_\parallel eV $, $T_\perp eV $' --title "Resolution scan in the direction of $\mu$" --legend "alpha = 2e-4, alpha = 2e-5"