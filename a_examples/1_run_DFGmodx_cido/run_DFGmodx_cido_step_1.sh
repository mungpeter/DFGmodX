
# run the first-pass to generate the initial setup files and sequence
# comparison

../1_auto_DFGmodx.py \
  ulk4_and_kit.fasta \
  cido               \
  10 5               \
  5                  \
  0 > run_DFGmodx_cido_step_1.result.txt

grep '> #.#' run_DFGmodx_cido_step_1.result.txt > run_DFGmodx_cido_step_1.warning.txt

cat run_DFGmodx_cido_step_1.warning.txt

