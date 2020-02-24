load /home/pmung/Dropbox/9_scripts/3_program/structures/3_DFGmodx/x_dataset/1atp.pdb, template
sele ref_resid, template and resi 122-138+162-183

load None, None
super None, ref_resid
save None.mod.pdb, None
hide everything
show ribbon
color white, all
 color red, template
save _TEMP.struct_super.super.pse