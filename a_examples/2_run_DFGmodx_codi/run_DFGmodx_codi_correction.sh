
# run the first-pass to generate the initial setup files and sequence
# comparison

../1_auto_DFGmodx.py \
  ulk4_and_kit.fasta \
  codi               \
  10 5               \
  5                  \
  1 > run_DFGmodx_codi_correction.result.txt

grep '> #.#' run_DFGmodx_codi_correction.result.txt > run_DFGmodx_codi_correction.warning.txt

../1_auto_DFGmodx.py \
  ulk4_and_kit.fasta \
  codi               \
  10 5               \
  5                  \
  2 > run_DFGmodx_codi_step_2.result.txt

grep '> #.#' run_DFGmodx_codi_step_2.result.txt > run_DFGmodx_codi_step_2.warning.txt

cat run_DFGmodx_codi_correction.warning.txt
cat run_DFGmodx_codi_step_2.warning.txt

