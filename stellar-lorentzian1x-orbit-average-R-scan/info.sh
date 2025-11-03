echo "Generating data for confinement time plot..."

cd R-scan/R-3
echo "Processing R=3 data..."
echo "Particle M0"
pgkyl --c2p gk_lorentzian_mirror-mc2nu_pos_deflated.gkyl gk_lorentzian_mirror-ion_M0_100.gkyl interp sel --z0 -0.98:0.98 integ 0 info
echo "Source M0"
pgkyl --c2p gk_lorentzian_mirror-mc2nu_pos_deflated.gkyl gk_lorentzian_mirror-ion_source_M0_0.gkyl interp sel --z0 -0.98:0.98 integ 0 info
echo "collision frequency"
pgkyl --c2p gk_lorentzian_mirror-mc2nu_pos_deflated.gkyl gk_lorentzian_mirror-ion_lbo_nu_sum_100.gkyl interp sel --z0 0.0 info
cd ../..
cd R-scan/R-5
echo "Processing R=5 data..."
echo "Particle M0"
pgkyl --c2p gk_lorentzian_mirror-mc2nu_pos_deflated.gkyl gk_lorentzian_mirror-ion_M0_100.gkyl interp sel --z0 -0.98:0.98 integ 0 info
echo "Source M0"
pgkyl --c2p gk_lorentzian_mirror-mc2nu_pos_deflated.gkyl gk_lorentzian_mirror-ion_source_M0_0.gkyl interp sel --z0 -0.98:0.98 integ 0 info
echo "collision frequency"
pgkyl --c2p gk_lorentzian_mirror-mc2nu_pos_deflated.gkyl gk_lorentzian_mirror-ion_lbo_nu_sum_100.gkyl interp sel --z0 0.0 info
cd ../..
cd R-scan/R-10
echo "Processing R=10 data..."
echo "Particle M0"
pgkyl --c2p gk_lorentzian_mirror-mc2nu_pos_deflated.gkyl gk_lorentzian_mirror-ion_M0_100.gkyl interp sel --z0 -0.98:0.98 integ 0 info
echo "Source M0"
pgkyl --c2p gk_lorentzian_mirror-mc2nu_pos_deflated.gkyl gk_lorentzian_mirror-ion_source_M0_0.gkyl interp sel --z0 -0.98:0.98 integ 0 info
echo "collision frequency"
pgkyl --c2p gk_lorentzian_mirror-mc2nu_pos_deflated.gkyl gk_lorentzian_mirror-ion_lbo_nu_sum_100.gkyl interp sel --z0 0.0 info
cd ../..
cd R-scan/R-15
echo "Processing R=15 data..."
echo "Particle M0"
pgkyl --c2p gk_lorentzian_mirror-mc2nu_pos_deflated.gkyl gk_lorentzian_mirror-ion_M0_100.gkyl interp sel --z0 -0.98:0.98 integ 0 info
echo "Source M0"
pgkyl --c2p gk_lorentzian_mirror-mc2nu_pos_deflated.gkyl gk_lorentzian_mirror-ion_source_M0_0.gkyl interp sel --z0 -0.98:0.98 integ 0 info
echo "collision frequency"
pgkyl --c2p gk_lorentzian_mirror-mc2nu_pos_deflated.gkyl gk_lorentzian_mirror-ion_lbo_nu_sum_100.gkyl interp sel --z0 0.0 info
cd ../..
cd R-scan/R-20
echo "Processing R=20 data..."
echo "Particle M0"
pgkyl --c2p gk_lorentzian_mirror-mc2nu_pos_deflated.gkyl gk_lorentzian_mirror-ion_M0_100.gkyl interp sel --z0 -0.98:0.98 integ 0 info
echo "Source M0"
pgkyl --c2p gk_lorentzian_mirror-mc2nu_pos_deflated.gkyl gk_lorentzian_mirror-ion_source_M0_0.gkyl interp sel --z0 -0.98:0.98 integ 0 info
echo "collision frequency"
pgkyl --c2p gk_lorentzian_mirror-mc2nu_pos_deflated.gkyl gk_lorentzian_mirror-ion_lbo_nu_sum_100.gkyl interp sel --z0 0.0 info
