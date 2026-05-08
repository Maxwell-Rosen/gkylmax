echo "Generating data for confinement time plot..."
echo ""

# Function to extract maximum value from pgkyl info output
extract_max() {
    grep "Maximum:" | awk '{print $3}'
}

M0=$(pgkyl --c2p gk_lorentzian_mirror-mc2nu_pos_deflated.gkyl gk_lorentzian_mirror-ion_M0_130.gkyl interp sel --z0 -4.0:4.0 integ 0 info 2>&1 | extract_max)
M0S=$(pgkyl --c2p gk_lorentzian_mirror-mc2nu_pos_deflated.gkyl gk_lorentzian_mirror-ion_source_M0_0.gkyl interp sel --z0 -4.0:4.0 integ 0 info 2>&1 | extract_max)
nu=$(pgkyl --c2p gk_lorentzian_mirror-mc2nu_pos_deflated.gkyl gk_lorentzian_mirror-ion_lbo_nu_sum_130.gkyl interp sel --z0 3.0 info 2>&1 | extract_max)


# Print Python-ready output
echo "# Copy the lines below to Python:"
echo ""
echo "intM0dx_tandem = $M0"
echo ""
echo "intM0Sdx_tandem = $M0S"
echo ""
echo "nu_ii_tandem = $nu"