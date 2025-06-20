# Copy the file clean.sh from stellar-wham1x-init-kinetic-from-boltzmann to all subdirectories in this directory
for dir in */ ; do
  if [ -d "$dir" ]; then
    cp stellar-wham1x-init-kinetic-from-boltzmann/clean.sh "$dir"
  fi
done
# This command finds all directories in the current directory (excluding hidden directories)
# and copies the file clean.sh from the stellar-wham1x-init-kinetic-from-boltzmann directory into each of those directories.
# The -not -path '*/\.*' part ensures that hidden directories (those starting with a dot) are not included in the search.
# The {} is a placeholder for each directory found, and \; indicates the end of the command to be executed for each directory.
# The cp command copies the file clean.sh into each of those directories.
# Make sure to run this script from the directory where you want to copy the clean.sh file
# and that the stellar-wham1x-init-kinetic-from-boltzmann directory exists and contains the clean.sh file.
# Also, ensure that you have the necessary permissions to write to the target directories. 

