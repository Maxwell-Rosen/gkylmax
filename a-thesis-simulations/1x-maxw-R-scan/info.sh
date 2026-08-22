echo "Generating data for confinement time plot..."
echo ""

# Function to extract a value from pgkyl info output (used after sel'ing a single point)
extract_max() {
    grep "Maximum:" | awk '{print $3}'
}

# R=3
cd R-scan/R-3
M0_R3=$(pgkyl gk_lorentzian_mirror-ion_integrated_moms.gkyl sel --z0 -1 -c 0 info 2>&1 | extract_max)
M0S_R3=$(pgkyl gk_lorentzian_mirror-ion_source_integrated_moms.gkyl sel --z0 -1 -c 0 info 2>&1 | extract_max)
nu_R3=$(pgkyl gk_lorentzian_mirror-ion_lbo_nu_sum_65.gkyl interp sel --z0 0.0 info 2>&1 | extract_max)
cd ../..

# R=5
cd R-scan/R-5
M0_R5=$(pgkyl gk_lorentzian_mirror-ion_integrated_moms.gkyl sel --z0 -1 -c 0 info 2>&1 | extract_max)
M0S_R5=$(pgkyl gk_lorentzian_mirror-ion_source_integrated_moms.gkyl sel --z0 -1 -c 0 info 2>&1 | extract_max)
nu_R5=$(pgkyl gk_lorentzian_mirror-ion_lbo_nu_sum_65.gkyl interp sel --z0 0.0 info 2>&1 | extract_max)
cd ../..

# R=10
cd R-scan/R-10
M0_R10=$(pgkyl gk_lorentzian_mirror-ion_integrated_moms.gkyl sel --z0 -1 -c 0 info 2>&1 | extract_max)
M0S_R10=$(pgkyl gk_lorentzian_mirror-ion_source_integrated_moms.gkyl sel --z0 -1 -c 0 info 2>&1 | extract_max)
nu_R10=$(pgkyl gk_lorentzian_mirror-ion_lbo_nu_sum_65.gkyl interp sel --z0 0.0 info 2>&1 | extract_max)
cd ../..

# R=15
cd R-scan/R-15
M0_R15=$(pgkyl gk_lorentzian_mirror-ion_integrated_moms.gkyl sel --z0 -1 -c 0 info 2>&1 | extract_max)
M0S_R15=$(pgkyl gk_lorentzian_mirror-ion_source_integrated_moms.gkyl sel --z0 -1 -c 0 info 2>&1 | extract_max)
nu_R15=$(pgkyl gk_lorentzian_mirror-ion_lbo_nu_sum_65.gkyl interp sel --z0 0.0 info 2>&1 | extract_max)
cd ../..

# R=22
cd R-scan/R-22
M0_R22=$(pgkyl gk_lorentzian_mirror-ion_integrated_moms.gkyl sel --z0 -1 -c 0 info 2>&1 | extract_max)
M0S_R22=$(pgkyl gk_lorentzian_mirror-ion_source_integrated_moms.gkyl sel --z0 -1 -c 0 info 2>&1 | extract_max)
nu_R22=$(pgkyl gk_lorentzian_mirror-ion_lbo_nu_sum_65.gkyl interp sel --z0 0.0 info 2>&1 | extract_max)
cd ../..

cd ../stellar-lorentzian1x-orbit-average
M0_R32=$(pgkyl gk_lorentzian_mirror-ion_integrated_moms.gkyl sel --z0 -1 -c 0 info 2>&1 | extract_max)
M0S_R32=$(pgkyl gk_lorentzian_mirror-ion_source_integrated_moms.gkyl sel --z0 -1 -c 0 info 2>&1 | extract_max)
nu_R32=$(pgkyl gk_lorentzian_mirror-ion_lbo_nu_sum_65.gkyl interp sel --z0 0.0 info 2>&1 | extract_max)
cd ../stellar-lorentzian1x-orbit-average-R-scan

# Print Python-ready output
echo "# Copy the lines below to Python:"
echo ""
echo "intM0dx_R3 = $M0_R3"
echo "intM0dx_R5 = $M0_R5"
echo "intM0dx_R10 = $M0_R10"
echo "intM0dx_R15 = $M0_R15"
echo "intM0dx_R22 = $M0_R22"
echo "intM0dx_R32 = $M0_R32"
echo ""
echo "intM0Sdx_R3 = $M0S_R3"
echo "intM0Sdx_R5 = $M0S_R5"
echo "intM0Sdx_R10 = $M0S_R10"
echo "intM0Sdx_R15 = $M0S_R15"
echo "intM0Sdx_R22 = $M0S_R22"
echo "intM0Sdx_R32 = $M0S_R32"
echo ""
echo "nu_ii_R3 = $nu_R3"
echo "nu_ii_R5 = $nu_R5"
echo "nu_ii_R10 = $nu_R10"
echo "nu_ii_R15 = $nu_R15"
echo "nu_ii_R22 = $nu_R22"
echo "nu_ii_R32 = $nu_R32"
