# Ask if the user is sure they want to delete all files
read -p "Are you sure you want to delete all files? (y/n): "
if [[ ! $REPLY =~ ^[Yy]$ ]]; then
    echo "Aborting deletion."
    exit 1
fiin

rm -rf BiMaxwellianMoments
rm -rf Distributions
rm -rf Field
rm -rf Geometry
rm -rf M
rm -rf PositivityFourMoments
rm -rf misc
rm *.gkyl
rm *.json
rm *.out
rm *.d