# Auto-IchorCNA

# Usage 
chmod +x Auto-IchorCNA.sh

./Auto-IchorCNA.sh [-i|n|s|r|c|o|t|h]

# Example 
./Auto-IchorCNA.sh -n WGS-001 -s WGS -r hg38 -t True -c num -o ~/Cancer/path/WGS-001/

# Optional argruments
    -i              Input file path (BAM file) *Default [all BAM in the current directory]
                    If you fill out the BAM path, you should have the .bam.bai file in the same
                    directory. If you choose the default value, this script will check for the
                    existence of .bam.bai and generate it if it does not exist.
    -n              Id name of patient *Default [same as BAM file name]
    -s              The type of next generation sequencing methods ('WGS' or 'WES') *Required
    -r              Reference genome type ('hg19' or 'hg38') *Required
    -c              chromosome name ('num' or 'chr') *Default [chr]
    -o              Output directory path *Default [current directory]
    -t              Temporary output (True or False) *Default [False]
    -h              Show this help message and exit

# Author's email
    songphon_sutthittha@cmu.ac.th
    
** Please contact us if you have any questions or problems with this script.
