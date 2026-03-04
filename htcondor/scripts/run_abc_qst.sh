#!/bin/bash
# Run ABC QST estimation (trait or neutral mode)
# Arguments: mode input ext_sd_or_file output_file num_sim summary_stats
#
# Supports two R environment transfer modes:
# 1. Traditional: r_env.tar.gz transferred and unpacked here
# 2. Auto-unpack: osdf:// with ?pack=auto unpacks during transfer (lib/, bin/ in cwd)

mode=$1
input=$2
ext_sd_or_file=$3
output_file=$4
num_sim=${5:-100000}
summary_stats=${6:-"QST,F_within_pop"}

wd=$(pwd)
echo "Running ABC QST estimation"
echo "Mode: $mode"
echo "Input: $input"
echo "Ext SD or file: $ext_sd_or_file"
echo "Num sim: $num_sim"
echo "Working directory: $wd"

# Determine R environment path based on transfer mode
if [ -f "r_env.tar.gz" ]; then
    # Traditional mode: tarball exists, unpack it
    echo "Found r_env.tar.gz - unpacking..."
    mkdir -p r_env
    tar -xzf r_env.tar.gz -C r_env
    tar_exit=$?
    echo "Tar extraction exit code: $tar_exit"
    if [ $tar_exit -ne 0 ]; then
        echo "ERROR: Failed to extract tarball"
        exit 1
    fi
    R_ENV_PATH="$wd/r_env"
    
    # List contents to verify extraction
    echo "Extracted contents:"
    ls -la "$R_ENV_PATH/" | head -10
    ls -la "$R_ENV_PATH/bin/" | head -10
    
    # Run conda-unpack to fix absolute paths (required for conda-pack'd environments)
    if [ -f "$R_ENV_PATH/bin/conda-unpack" ]; then
        echo "Running conda-unpack to fix paths..."
        source "$R_ENV_PATH/bin/activate"
        "$R_ENV_PATH/bin/conda-unpack" 2>&1 || echo "conda-unpack completed"
    fi
elif [ -d "lib/R" ]; then
    # Auto-unpack mode: contents already extracted to current directory
    echo "Auto-unpack detected - using current directory as R environment"
    R_ENV_PATH="$wd"
    
    if [ -f "$R_ENV_PATH/bin/conda-unpack" ]; then
        $R_ENV_PATH/bin/python $R_ENV_PATH/bin/conda-unpack 2>&1 || true
    fi
else
    echo "ERROR: No R environment found!"
    echo "Expected either r_env.tar.gz or lib/R directory"
    ls -la
    exit 1
fi

echo "R environment path: $R_ENV_PATH"

# Set environment variables
export PATH=$R_ENV_PATH/bin:$PATH
export R_HOME=$R_ENV_PATH/lib/R
export R_LIBS=$R_ENV_PATH/lib/R/library
export R_LIBS_USER=$R_ENV_PATH/lib/R/library
export LD_LIBRARY_PATH=$R_ENV_PATH/lib:$R_ENV_PATH/lib/R/lib:$LD_LIBRARY_PATH

echo "Environment:"
echo "  PATH: $PATH"
echo "  R_HOME: $R_HOME"
echo "  R_LIBS: $R_LIBS"
echo "  LD_LIBRARY_PATH: $LD_LIBRARY_PATH"
echo "  R_HOME exists: $(test -d $R_HOME && echo yes || echo no)"
echo "  R_LIBS exists: $(test -d $R_LIBS && echo yes || echo no)"

# Verify Rscript is available
if [ ! -f "$R_ENV_PATH/bin/Rscript" ]; then
    echo "ERROR: Rscript not found at $R_ENV_PATH/bin/Rscript"
    ls -la "$R_ENV_PATH/bin/" | head -20
    exit 1
fi

# Check library dependencies
echo "Checking library dependencies..."
ldd "$R_ENV_PATH/bin/Rscript" 2>&1 | grep -i "not found" && echo "WARNING: Some libraries missing"

# Verify R script exists
if [ ! -f "$wd/qst_abc_sim.R" ]; then
    echo "ERROR: R script not found at $wd/qst_abc_sim.R"
    ls -la "$wd/"
    exit 1
fi

# Show working directory contents
echo "Working directory contents:"
ls -la "$wd/"

echo "Running Rscript..."
echo "Command: $R_ENV_PATH/bin/Rscript $wd/qst_abc_sim.R $mode $input $ext_sd_or_file $output_file $num_sim $summary_stats"

# Verify Rscript is executable
echo "Rscript file info:"
file "$R_ENV_PATH/bin/Rscript"
ls -la "$R_ENV_PATH/bin/Rscript"

# Test if we can run Rscript at all
echo "Testing Rscript --version:"
"$R_ENV_PATH/bin/Rscript" --version 2>&1 || echo "Rscript --version failed with exit code $?"

# Check library dependencies
echo "Checking R library shared objects..."
ldd "$R_ENV_PATH/lib/R/library/abc/libs/abc.so" 2>&1 | grep -i "not found" && echo "WARNING: abc has missing libraries"
ldd "$R_ENV_PATH/lib/R/library/e1071/libs/e1071.so" 2>&1 | grep -i "not found" && echo "WARNING: e1071 has missing libraries"

# Try loading R packages
echo "Testing R package loading..."
"$R_ENV_PATH/bin/Rscript" -e "library(abc); cat('abc loaded\n')" 2>&1 || echo "Package loading failed"

# Debug: List library contents
echo "Checking library contents..."
ls "$R_ENV_PATH/lib/R/library/abc/libs/" 2>&1 || echo "No abc libs dir"

# Debug: Check if R can start at all with verbose output
echo "Testing R startup..."
"$R_ENV_PATH/bin/R" --vanilla --slave -e "cat('R startup test OK\n')" 2>&1 || echo "R startup failed"

# Run R script with explicit error capture - combine stdout and stderr  
echo "Starting R script execution..."
set +e  # Don't exit on error
# Use R CMD BATCH style for better error reporting
"$R_ENV_PATH/bin/Rscript" --vanilla "$wd/qst_abc_sim.R" \
    "$mode" \
    "$input" \
    "$ext_sd_or_file" \
    "$output_file" \
    "$num_sim" \
    "$summary_stats" 2>&1
exit_code=$?
set -e

echo "R script finished."
echo "Exit code: $exit_code"

if [ -f "$output_file" ]; then
    echo "Output file created: $output_file"
    ls -la "$output_file"
else
    echo "ERROR: Output file not created"
    exit 1
fi
