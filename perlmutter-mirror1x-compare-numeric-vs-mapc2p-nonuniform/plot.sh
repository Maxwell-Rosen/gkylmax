name1="mc2p/gk_mirror"
name2="numeric/gk_mirror"
saveLoc="./python-plots/"

# pgkyl "$name1"-b_i.gkyl interp -b ms -p1 pl --title "b_i" --saveas "$saveLoc-geo-b_i.png" --no-show &
# pgkyl "$name1"-bmag_inv_sq.gkyl interp -b ms -p1 pl --title "bmag_inv_sq" --saveas "$saveLoc-geo-bmag_inv_sq.png" --no-show &
# pgkyl "$name1"-bmag_inv.gkyl interp -b ms -p1 pl --title "bmag_inv" --saveas "$saveLoc-geo-bmag_inv.png" --no-show &
# pgkyl "$name1"-bmag.gkyl interp -b ms -p1 pl --title "bmag" --saveas "$saveLoc-geo-bmag.png" --no-show &
# pgkyl "$name1"-cmag.gkyl interp -b ms -p1 pl --title "cmag" --saveas "$saveLoc-geo-cmag.png" --no-show &
# pgkyl "$name1"-jacobtot_inv.gkyl interp -b ms -p1 pl --title "jacobtot_inv" --saveas "$saveLoc-geo-jacobtot_inv.png" --no-show &
# pgkyl "$name1"-jacobgeo.gkyl interp -b ms -p1 pl --title "jacobgeo" --saveas "$saveLoc-geo-jacobgeo.png" --no-show &
# pgkyl "$name1"-jacobtot.gkyl interp -b ms -p1 pl --title "jacobtot" --saveas "$saveLoc-geo-jacobtot.png" --no-show &
# pgkyl "$name1"-mapc2p.gkyl interp -b ms -p1 pl --title "mapc2p" --saveas "$saveLoc-geo-mapc2p.png" --no-show &
# pgkyl "$name1"-ion_BiMaxwellianMoments_0.gkyl interp pl --title "ion_BiMaxwellianMoments_0" --saveas "$saveLoc-ion_BiMaxwellianMoments_0.png" --no-show &

# pgkyl "$name1"-b_i.gkyl "$name2"-b_i.gkyl interp -b ms -p1 pl --title "b_i" --saveas "$saveLoc-geo-b_i.png" --no-show -f0 &
# pgkyl "$name1"-bmag_inv_sq.gkyl "$name2"-bmag_inv_sq.gkyl interp -b ms -p1 pl --title "bmag_inv_sq" --saveas "$saveLoc-geo-bmag_inv_sq.png" --no-show -f0 &
# pgkyl "$name1"-bmag_inv.gkyl "$name2"-bmag_inv.gkyl interp -b ms -p1 pl --title "bmag_inv" --saveas "$saveLoc-geo-bmag_inv.png" --no-show -f0 &
# # pgkyl "$name1"-bmag.gkyl "$name2"-bmag.gkyl interp -b ms -p1 pl --title "bmag" --saveas "$saveLoc-geo-bmag.png" --no-show -f0 &
# pgkyl "$name1"-cmag.gkyl "$name2"-cmag.gkyl interp -b ms -p1 pl --title "cmag" --saveas "$saveLoc-geo-cmag.png" --no-show -f0 &
# pgkyl "$name1"-jacobtot_inv.gkyl "$name2"-jacobtot_inv.gkyl interp -b ms -p1 pl --title "jacobtot_inv" --saveas "$saveLoc-geo-jacobtot_inv.png" --no-show -f0 &
# pgkyl "$name1"-jacobgeo.gkyl "$name2"-jacobgeo.gkyl interp -b ms -p1 pl --title "jacobgeo" --saveas "$saveLoc-geo-jacobgeo.png" --no-show -f0 &
# pgkyl "$name1"-jacobtot.gkyl "$name2"-jacobtot.gkyl interp -b ms -p1 pl --title "jacobtot" --saveas "$saveLoc-geo-jacobtot.png" --no-show -f0 &
# pgkyl "$name1"-mapc2p.gkyl "$name2"-mapc2p.gkyl interp -b ms -p1 pl --title "mapc2p" --saveas "$saveLoc-geo-mapc2p.png" --no-show -f0 &
# pgkyl "$name1"-ion_BiMaxwellianMoments_0.gkyl "$name2"-ion_BiMaxwellianMoments_0.gkyl interp pl --title "ion_BiMaxwellianMoments_0" --saveas "$saveLoc-ion_BiMaxwellianMoments_0.png" --no-show -f0 &

# pgkyl "$name1"-b_i.gkyl -t mc2p "$name2"-b_i.gkyl -t numeric interp -b ms -p1 activ -t mc2p ev "f 1.25069595987 /" pl --title "b_i" --xscale 1.25069595987 --no-show -f0 activ -t numeric pl -f0 --no-show --saveas "$saveLoc-geo-b_i.png" &
# pgkyl "$name1"-bmag.gkyl -t mc2p "$name2"-bmag.gkyl -t numeric interp -b ms -p1 activ -t mc2p pl --title "bmag" --xscale 1.25069595987 --no-show -f0 activ -t numeric pl -f0 --no-show --saveas "$saveLoc-geo-bmag.png" &
# pgkyl "$name1"-bmag_inv.gkyl -t mc2p "$name2"-bmag_inv.gkyl -t numeric interp -b ms -p1 activ -t mc2p pl --title "bmag_inv" --xscale 1.25069595987 --no-show -f0 activ -t numeric pl -f0 --no-show --saveas "$saveLoc-geo-bmag_inv.png" &
# pgkyl "$name1"-bmag_inv_sq.gkyl -t mc2p "$name2"-bmag_inv_sq.gkyl -t numeric interp -b ms -p1 activ -t mc2p pl --title "bmag_inv_sq" --xscale 1.25069595987 --no-show -f0 activ -t numeric pl -f0 --no-show --saveas "$saveLoc-geo-bmag_inv_sq.png" &
# pgkyl "$name1"-cmag.gkyl -t mc2p "$name2"-cmag.gkyl -t numeric interp -b ms -p1 activ -t mc2p pl --title "cmag" --xscale 1.25069595987 --no-show -f0 activ -t numeric pl -f0 --no-show --saveas "$saveLoc-geo-cmag.png" &

# pgkyl "$name1"-jacobgeo.gkyl -t mc2p "$name2"-jacobgeo.gkyl -t numeric interp -b ms -p1 activ -t mc2p ev "f 1.25069595987 /" pl --title "jacobgeo" --xscale 1.25069595987 --no-show -f0 activ -t numeric pl -f0 --no-show --saveas "$saveLoc-geo-jacobgeo.png" &

# pgkyl "$name1"-jacobtot.gkyl -t mc2p "$name2"-jacobtot.gkyl -t numeric interp -b ms -p1 activ -t mc2p ev "f 1.25069595987 /" pl --title "jacobtot" --xscale 1.25069595987 --no-show -f0 activ -t numeric pl -f0 --no-show --saveas "$saveLoc-geo-jacobtot.png" &
# pgkyl "$name1"-jacobtot_inv.gkyl -t mc2p "$name2"-jacobtot_inv.gkyl -t numeric interp -b ms -p1 activ -t mc2p ev "f 1.25069595987 *" pl --title "jacobtot_inv" --xscale 1.25069595987 --no-show -f0 activ -t numeric pl -f0 --no-show --saveas "$saveLoc-geo-jacobtot_inv.png" &
# pgkyl "$name1"-mapc2p.gkyl -t mc2p "$name2"-mapc2p.gkyl -t numeric interp -b ms -p1 activ -t mc2p pl --title "mapc2p" --xscale 1.25069595987 --no-show -f0 activ -t numeric pl -f0 --no-show --saveas "$saveLoc-geo-mapc2p.png" &
# pgkyl "$name1"-ion_BiMaxwellianMoments_0.gkyl -t mc2p "$name2"-ion_BiMaxwellianMoments_0.gkyl -t numeric interp -b ms -p1 activ -t mc2p pl --title "ion_BiMaxwellianMoments_0" --xscale 1.25069595987 --no-show -f0 activ -t numeric pl -f0 --no-show --saveas "$saveLoc-ion_BiMaxwellianMoments_0.png" --linestyle dashed --legend "mc2p" "numeric" &

frame=100
for frame in {0..100}
do
pgkyl mc2p/BiMaxwellianMoments/gk_mirror-ion_BiMaxwellianMoments_$frame.gkyl -t "mc2p" numeric/BiMaxwellianMoments/gk_mirror-ion_BiMaxwellianMoments_$frame.gkyl -t numeric interp -b ms -p1 activ -t mc2p pl --title "frame $frame" --xscale 1.25069595987 --no-show -f0 activ -t numeric pl -f0 --no-show --saveas "$saveLoc-ion_BiMaxwellianMoments_"$frame".png" --linestyle dashed --legend "mc2p" "numeric" &
done