# Ask user for confirmation before cleaning up
echo "This script will remove all files and directories related to the WHAM2x simulation."
read -p "Are you sure you want to proceed? (y/n): " confirmation
if [[ "$confirmation" != "y" ]]; then
  echo "Cleanup aborted."
  exit 1
fi

rm -rf BiMaxwellianMoments
rm -rf Distributions
rm -rf Field
rm -rf Geometry
rm -rf M
rm -rf PositivityFourMoments
rm -rf misc
rm -rf Coll
rm -rf PrimMoms
rm -rf Source
rm *.gkyl
rm *.json
rm *.out
rm *.d