# Copy to config.R and set DATA_DIR to your own MIMIC-IV v2.2 directory.
# MIMIC-IV requires PhysioNet credentialed access and is not redistributable.

DATA_DIR <- "/path/to/physionet.org/files/mimiciv/2.2"

# A fit is called divergent when the interaction coefficient runs past this.
# The cutoff is set against the scale of the target: full-pool estimates are
# around 0.05, so 2 is forty times the truth.
EXPLODE_THRESHOLD <- 2
