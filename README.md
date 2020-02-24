# DFGModx

Generate Typical Protein Kinase Catalytic Domain in Inactive Conformations (CIDO/CODI/CODO)

- Step 1: compare the input kinase FASTA sequence to an internal library of protein kinase with atomic structures available
- Step 2: use the kinase structures that is most similar to the input kinase sequence as Base structure
- Step 3: align input kinase sequence to the selected kinase structures and a set of template kinase structures in alternative conformation (CIDI/CIDO/CODI/CODO/wCD)
- Step 4: generate homology models of input kinase based on the Base kinase structure and the template structures in alternative conformation
- Step 5: calculate the volume of the active site, rank and select ones with the largest volume

CIDO (aC-helix-in/DFG-out)
This is the typical DFG-out conformation described in the literature, where DFG-motif has a 180-degree flip and brings the DFG-Phe phenyl ring out from the hydrophobic deep pocket in the kinase, leaving the deep DFG-pocket vacant and available for small molecule binding.
This range of movement can be well-captured by **one population** (set) of kinase structures in DFG-out conformation, where DFG-motif has a clear flip, while the set captures the range of small-lobe movement observed. Because the TK family has been heavily studied and with the most number of knownn CIDO structures, there are **2** sets of CIDO templates, one of TK family and another one for all other kinases.

* When doing molecular docking to these CIDO models, I found that default settings for OpenEye FRED works very well, while Schrö
dinger GLIDE will need a 0.75x adjustment to the GRID van der Waals radius to get similar results. In both cases, should also set a hydrophobic sphere in the DFG-pocket to enforce docking to that site while having H-bond constraints to the hinge residue.

CODI (aC-helix-out/DFG-in)
This is the C-helix-out conformation described in the literature as the "CDK-like" or "Src-like" inactive conformation. This conformation encompasses a range of aC-helix movements - translational, rotational, angular - that is more heterogeneous than DFG-motif movement (mostly just -in and -out). 
Due to the large range of aC-helix movemnents, there are actually **several (4+) sub-populations** of CODI conformations. **4** of them are often seen in multiple kinases in different families, while some are family-specific. To generalize the modeling of CODI conformation, **4 sub-populations** (sets) are used in the modeling process as compared to CIDO model generation.
Notably, MEK has been seen to adopt a unique CODI conformation when bound to type-III inhibitor and ATP together, while cMET adopts another unique CODI conformation held up by an electrostatic interaction on its activation loop. These case-specific CODI conformations are not included in the general CODI sub-populations used for modeling but are also available for use.


_Reference 1_: [Ung, P.M.U.; Schlessinger, A. DFGmodel: predicting protein kinase structures in inactive states for structure-based discovery of type-II inhibitors. ACS Chem. Biol. 2015, 10(1):269-278.](https://doi.org/10.1021/cb500696t)

_Reference 2_: [\*Ung PMU, \*Rahman R, Schlessinger A. Redefining the Protein Kinase Conformational Space with Machine Learning. Cell Chemical Biology (2018) 25(7), 916-924.](https://doi.org/10.1016/j.chembiol.2018.05.002)

_Reference 3_: [\*Rahman R, \*Ung PMU, Schlessinger A. KinaMetrix: a web resource to investigate kinase conformations and inhibitor space. Nucleic Acids Research (2019) 47(D1), D361–D366.](https://doi.org/10.1093/nar/gky916)

_Website_: [KinaMetrix.org](http://kinametrix.org)

######################################################################################
# Primary script: Modeling multiple kinases with simple default settings (multi-steps)
```
> 1_auto_DFGmodx.py
    [ FASTA file of kinase(s) to be modeled ]
    [ Conformation: cido | codi | codo | cidi | wcd ]   ** capital or lowercase
    [ Number of models to generate ]
    [ Number of top-models to take ]
    [ CPU per Model generation run ]
    [ Step: 0 - Setup homology modeling (default)                   ]
    [       1 - Setup homology modeling (_TEMP.y.fasta corrected)   ]
    [           #( rename *.y1.fasta to corrected *.y1.corr.fasta ) ]
    [       2 - Run homology modeling                               ]
    [       3 - Run volume calculation and selection                ]
    
    [ Optional Conformation: III\t- CODI MEK1-like for type-III inh ]
    [                        MET\t- CODI cMET-like altered A-loop   ]
    [                        EGFR\t-CODI EGFR-like altered A-loop   ]
    
e.g.> ./1_auto_DFGmodx.py      \
          ulk4_and_kit.fasta   \
          cido  10 5 5         \
          0           # "0" will stop after generate initial setup
                      # "1" will restart those x.fasta ones with _TEMP.y.fasta correction
                      # "2" will *only* generate actual models + volumes with successful setup
                      # "3" will *only* use generated models to calculate volumes
```
For general purpose and simply automated setup/running of DFGModx, use **1_auto_DFGMod.py**. This is the wrapper script of the backbone script *1_run_single_DFGmodx.py* with all default settings for modeling kinase structures in different conformations. The entire modeling takes 2-3 steps:
- **Step 0**: generate the initial setup files and folders for all kinases being modeled with default setting
- **Step 2**: generate kinase models with default settings, results in the **/1_dfgmod** folder
- *step 1*: an optional step to fix problems appeared in "Step 0" by introducing a **fasta file** with **corrected** alignment, with nameing convention: \_TEMP.{0}.y1.fasta, where {0} = model output prefix.
- *step 3*: an optional step to run homology model binding site volume calculation. Generally won't need to use this unless Step 2 failed.

######################################################################################
# Optional script: Run DFGmodx on individual kinase with more options (multi-step)
```
> 1_run_single_DFGmodx.py
      [ Setup Script: formatted setup file | Pickle ]\n
      [ Mode: -full   Perform modified alignment, DFGmodel, and selection
              -pir    Run modified alignment only
              -mod    Run DFGmodel only
              -vol    Select top models only
              -set    Generate template setup script\n
    ** Always use -pir then -mod, and ALWAYS check .pir for additional chain
       breaks. Only use -full if you have checked the kinase and it has no
       additional chain break or missing residues\n
      [ Opt:   -force   Forced restart of /dfgworking directory ]
      [ Opt:   -restart Forced restart of /1_result   directory ]
      [ Opt:   -pass    Pass along; always use this unless otherwise ]\n

e.g.>  ./1_run_single_DFGmodx.py setup.file -pir -pass
```
This is the backbone script to manually setup and modeling of a single kinase. To run this, several steps are needed:
- **Step 0**: generate the initial template setup script when there is nothing
```> ./1_run_single_DFGmodx.py -set ```
- **Step 1**: after modifying the setup file, generate the .pir file used by Modeller modeling. Recommended to check .pir file for errors before proceeding to model generation
```> ./1_run_single_DFGmodx.py setup.txt -pir -pass```
- **Step 2**: generate models using the **.pir**  file. If no error occurs, volume calculation will start right after modeling.
```> ./1_run_single_DFGmodx.py setup.txt -mod -pass```
- _Step 3_: optional step, use if volume was not calculated right after the modeling step
```> ./1_run_single_DFGmodx.py setup.txt -vol -pass```

###################################################################################
# Optional script: Superimpose all kinase structures to a common reference frame
```
> 1_pymol_alignment.py
      [Template PDB]                      * usually 1ATP.pdb
      [Template Residues to superimpose]  * 'resi 122-138+162-181' PyMOL naming
      [List of Model PDBs] 
      [Output PyMOL session prefix]
      [Name of Superimposed PDB List]
      [Number of Model] (for Modeller PDB, B999000xx)
      
e.g.> ./1_pymol_alignment.py 1atp.pdb 'resi 122-138+162-181' model.list align_output _tmp.list 1
```
This script uses PyMOL as the main driver to superimpose all kinase PDB structures to a common reference frame, often 1ATP.pdb large-lobe, then save the superimposed stuctures with **.mod.pdb** ending. All structures used in DFGmodx, should use this to superimpose the structures to the common frame.
* *
- PyMOL has a pretty good structure superimpose algorithm that is independent of the atom number and sequence matching, which is great for superimposing structures from different kinases with different lengths and residues.
- The **'super'** function is sequence-independent while **'align'** is sequence-dependent, and **'super'** is the primary function to be used. Becuase some kinases have different C-lobe helices arrangement, it is safest to align to two upmost conserved region: hinge region + D-helix, and catalytic loop + beta-hairpin, '1atp and resid 122-138+162-181'.
- If the input structures are homology models from Modeller, they will have **.B999000** in the filename. For these model structures, the filename will be condensed to generate {0}.B*.mod.pdb filename. If there are 10+ files, the first 9 models will have zero, i.e. .B0*, to match up for easier sorting.
- If the first round of **'super'** superimposing fails due to 'Insufficient Number of Atom for Structure Superposition', that is a problem of older PyMOL that should be fixed by now but if that happens, **'align'** will be automatically used instead to superimpose the structures. Although **'align'** gives slightly different superimposed structure than **'super'**, the difference is minor.

###################################################################################
# Other scripts:

```
> 2_pir_result_info.py
  [ fasta file with Kinase name ] [ conformation ] [ output prefix ]
  [ conformation: cidi -- cidi_s_templ or cidi_y_templ ]
  [               cido -- cido_s_templ or cido_y_templ ]
  [               codi -- codi_templ                   ]
  [               codo -- codo_templ                   ]
  [               wcd  -- wcd_templ                    ]

```


```
> 3_extract_dope_vol.py

```


######################################################################################
# Example on running DFGmodx
```
/a_examples
     |-----/1_run_DFGmodx_cido
     |           |..... run_DFGmodx_cido_step_1.sh
     |           |..... run_DFGmodx_cido_step_2.sh
     |           |..... run_DFGmodx_cido_correction.sh
     |           |..... troubleshoot_common_issues.txt
     |           |..... ulk4_kit.fasta
     |           |...../correct_results/
     |                         |...../KIT
     |                         |...../ULK4
     |                                 |...../0_
     |                                 |...../1_
     |                                 |..... 
     |                                 |..... 
     |                                 |..... 
     |                                 |..... 
     |                                 |..... 
     |-----/2_run_DFGmodx_codi
                 |..... run_DFGmodx_codi_step_1.sh
                 |..... run_DFGmodx_codi_step_2.sh
                 |..... run_DFGmodx_codi_correction.sh
                 |..... troubleshoot_common_issues.txt
                 |..... ulk4_kit.fasta
                 |...../correct_results/
                               |...../KIT
                               |...../ULK4
                                        |...../0_
                                        |...../1_
                                        |..... ...

```
In this example, 2 kinases (ULK4 and KIT) are used. 
* KIT modeling should go without a hitch since it has a fairly typical sequence.
* Standard sequence alignment will fail for ULK4 because it has 1 fewer large-lobe helix than other typical kinases and will generate bad models (zDOPE > +0.0), while the pseudokinase nature makes the automated alignment of conserved region bad, hence will need to do a alignment correction before proceeding.

At the end of every step, the log file should be checked for **WARNING** and **ERROR** flags, especially when many kinases are being modeled at the same time. Most common errors are misalignment of the conserved region, mismatching between sequence and template structures.

###################################################################################
* Required packages:
```
  POVME          # 2.1     standalone program
  PyMOL          # 1.8+    standalone program
  Modeller       # 9.20+   standalone program or Conda
  Muscle         # 3.8.51
  T-Coffee       # 13

  Python         # 3.7.2+
    numpy          # 1.17.3
    pandas         # 0.25.3
    pathos         # 0.2.3
    biopython      # 1.74
    zlib           # 1.2.11
    bzip2          # 1.0.6
    plotnine       # 
```
