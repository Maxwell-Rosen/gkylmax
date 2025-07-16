# Ask for user confirmation before deleting files
echo "This script will delete categorized and .gkyl files in directory."
read -p "Are you sure you want to continue? (y/n): " confirm
if [[ "$confirm" != "y" ]]; then
  echo "Aborting deletion."
  exit 1
fi

rm -rf BiMaxwellianMoments
rm -rf Distributions
rm -rf Field
rm -rf Geometry
rm -rf M
rm -rf PositivityFourMoments
rm -rf PrimMoms
rm -rf Source
rm -rf Coll
rm -rf misc
rm *.gkyl
rm *.json
rm *.out
rm *.d